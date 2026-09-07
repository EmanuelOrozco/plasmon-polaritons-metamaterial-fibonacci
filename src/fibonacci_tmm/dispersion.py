"""Relación de dispersión cos(k Lm) = Rm y detección de bandas/gaps."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from fibonacci_tmm.constants import BAND_ABS_R_TOLERANCE, GHZ
from fibonacci_tmm.electromagnetics import Polarization
from fibonacci_tmm.fibonacci import cell_length
from fibonacci_tmm.params import SuperlatticeSpec
from fibonacci_tmm.transfer_matrix import semitrace_by_product, semitrace_by_recurrence


@dataclass(frozen=True)
class FrequencyInterval:
    nu_min_ghz: float
    nu_max_ghz: float
    n_points: int

    @property
    def bandwidth_ghz(self) -> float:
        return self.nu_max_ghz - self.nu_min_ghz

    @property
    def center_ghz(self) -> float:
        return 0.5 * (self.nu_min_ghz + self.nu_max_ghz)


@dataclass(frozen=True)
class DispersionScan:
    nu_ghz: NDArray[np.float64]
    omega: NDArray[np.float64]
    r: NDArray[np.complex128]
    abs_r: NDArray[np.float64]
    allowed: NDArray[np.bool_]
    k_lm_over_pi: NDArray[np.float64]
    method: str
    m: int
    theta: float
    polarization: str
    clipped_near_unit: int


def omega_from_nu_ghz(nu_ghz: ArrayLike) -> NDArray[np.float64]:
    return 2.0 * np.pi * np.asarray(nu_ghz, dtype=np.float64) * GHZ


def nu_ghz_from_omega(omega: ArrayLike) -> NDArray[np.float64]:
    return np.asarray(omega, dtype=np.float64) / (2.0 * np.pi * GHZ)


def real_semitrace(r: NDArray[np.complex128]) -> NDArray[np.float64]:
    """Parte real de Rm. Im(R) debe ser ~0 en el modelo sin pérdidas."""
    return np.real(r)


def allowed_mask(
    r: NDArray[np.complex128],
    *,
    atol: float = BAND_ABS_R_TOLERANCE,
) -> tuple[NDArray[np.bool_], int]:
    """|R| ≤ 1 + atol. El entero es cuántos puntos cayeron en (1, 1+atol]."""
    abs_r = np.abs(r)
    strictly_allowed = abs_r <= 1.0
    clipped = np.count_nonzero((abs_r > 1.0) & (abs_r <= 1.0 + atol))
    allowed = abs_r <= 1.0 + atol
    del strictly_allowed
    return allowed, int(clipped)


def reduced_wavevector(
    r: NDArray[np.complex128],
    allowed: NDArray[np.bool_],
) -> NDArray[np.float64]:
    """k Lm / π = arccos(R)/π ∈ [0, 1] donde hay banda; NaN en gaps."""
    r_real = np.real(r)
    # Clip solo donde already allowed, para arccos del dominio real.
    r_clipped = np.clip(r_real, -1.0, 1.0)
    k_over_pi = np.full(r_real.shape, np.nan, dtype=np.float64)
    k_over_pi[allowed] = np.arccos(r_clipped[allowed]) / np.pi
    return k_over_pi


def scan_dispersion(
    spec: SuperlatticeSpec,
    m: int,
    theta: float,
    nu_ghz: ArrayLike,
    polarization: Polarization = Polarization.TE,
    *,
    method: str = "recurrence",
    atol: float = BAND_ABS_R_TOLERANCE,
) -> DispersionScan:
    nu = np.asarray(nu_ghz, dtype=np.float64)
    omega = omega_from_nu_ghz(nu)
    # Evitar los polos exactos de Drude (μ=0 o ε=0), donde Ψ no está definido.
    pole_e = spec.medium_b.omega_e
    pole_m = spec.medium_b.omega_m
    near_pole = (np.abs(omega - pole_e) < 1.0e-6 * pole_e) | (np.abs(omega - pole_m) < 1.0e-6 * pole_m)
    if method == "recurrence":
        r = semitrace_by_recurrence(spec, m, omega, theta, polarization)
    elif method == "product":
        r = semitrace_by_product(spec, m, omega, theta, polarization)
    else:
        raise ValueError("method debe ser 'recurrence' o 'product'")
    allowed, n_clipped = allowed_mask(r, atol=atol)
    allowed = allowed & ~near_pole
    k_red = reduced_wavevector(r, allowed)
    return DispersionScan(
        nu_ghz=nu,
        omega=omega,
        r=r,
        abs_r=np.abs(r),
        allowed=allowed,
        k_lm_over_pi=k_red,
        method=method,
        m=m,
        theta=theta,
        polarization=polarization.value,
        clipped_near_unit=n_clipped,
    )


def allowed_intervals(scan: DispersionScan) -> list[FrequencyInterval]:
    """Componentes conexas de frecuencias con |R| ≤ 1."""
    allowed = scan.allowed
    nu = scan.nu_ghz
    intervals: list[FrequencyInterval] = []
    n = allowed.size
    i = 0
    while i < n:
        if not allowed[i]:
            i += 1
            continue
        j = i
        while j + 1 < n and allowed[j + 1]:
            j += 1
        intervals.append(
            FrequencyInterval(
                nu_min_ghz=float(nu[i]),
                nu_max_ghz=float(nu[j]),
                n_points=int(j - i + 1),
            )
        )
        i = j + 1
    return intervals


def average_epsilon_mu(
    spec: SuperlatticeSpec,
    m: int,
    omega: ArrayLike,
) -> tuple[NDArray[np.complex128], NDArray[np.complex128]]:
    """⟨ε⟩_m y ⟨μ⟩_m del paper."""
    from fibonacci_tmm.fibonacci import n_layers_a, n_layers_b

    omega_arr = np.asarray(omega, dtype=np.float64)
    eps_a, mu_a = spec.medium_a.epsilon_mu(omega_arr)
    eps_b, mu_b = spec.medium_b.epsilon_mu(omega_arr)
    length = cell_length(m, spec.thickness_a, spec.thickness_b)
    n_a = n_layers_a(m)
    n_b = n_layers_b(m)
    avg_eps = (n_a * eps_a * spec.thickness_a + n_b * eps_b * spec.thickness_b) / length
    avg_mu = (n_a * mu_a * spec.thickness_a + n_b * mu_b * spec.thickness_b) / length
    return avg_eps, avg_mu


def zero_average_index_frequency_ghz(spec: SuperlatticeSpec, m: int) -> float | None:
    """Frecuencia donde ⟨ε⟩=⟨μ⟩=0 cuando ω_e = ω_m y A es aire, a = b.

    ν = ν_p √(F(m-2)/F(m)). None si las hipótesis no se cumplen.
    """
    from fibonacci_tmm.fibonacci import n_layers_b, n_layers_total

    equal_plasma = np.isclose(spec.omega_e, spec.omega_m)
    equal_thickness = np.isclose(spec.thickness_a, spec.thickness_b)
    air = np.isclose(float(np.real(spec.medium_a.epsilon)), 1.0) and np.isclose(
        float(np.real(spec.medium_a.mu)), 1.0
    )
    if not (equal_plasma and equal_thickness and air):
        return None
    nu_p = spec.nu_e_ghz()
    return nu_p * float(np.sqrt(n_layers_b(m) / n_layers_total(m)))
