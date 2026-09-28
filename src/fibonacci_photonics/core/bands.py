"""Representación de un barrido de dispersión, independiente del método numérico.

Tanto la TMM como el PWE entregan la semitraza R(ν) = cos(k Lm) de la celda Sm.
A partir de ella se clasifican bandas (|R| ≤ 1) y gaps (|R| > 1), se obtiene
k Lm / π y se extraen los intervalos de frecuencia permitidos.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from fibonacci_photonics.core.constants import BAND_ABS_R_TOLERANCE, GHZ
from fibonacci_photonics.core.electromagnetics import Polarization
from fibonacci_photonics.core.params import SuperlatticeSpec


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
    clipped = np.count_nonzero((abs_r > 1.0) & (abs_r <= 1.0 + atol))
    allowed = abs_r <= 1.0 + atol
    return allowed, int(clipped)


def reduced_wavevector(
    r: NDArray[np.complex128],
    allowed: NDArray[np.bool_],
) -> NDArray[np.float64]:
    """k Lm / π = arccos(R)/π ∈ [0, 1] donde hay banda; NaN en gaps."""
    r_real = np.real(r)
    r_clipped = np.clip(r_real, -1.0, 1.0)
    k_over_pi = np.full(r_real.shape, np.nan, dtype=np.float64)
    k_over_pi[allowed] = np.arccos(r_clipped[allowed]) / np.pi
    return k_over_pi


def near_drude_pole(spec: SuperlatticeSpec, omega: NDArray[np.float64]) -> NDArray[np.bool_]:
    """Frecuencias pegadas a μ_B = 0 o ε_B = 0, donde Ψ no está definido."""
    pole_e = spec.medium_b.omega_e
    pole_m = spec.medium_b.omega_m
    return (np.abs(omega - pole_e) < 1.0e-6 * pole_e) | (np.abs(omega - pole_m) < 1.0e-6 * pole_m)


def scan_from_semitrace(
    spec: SuperlatticeSpec,
    m: int,
    theta: float,
    nu_ghz: ArrayLike,
    r: ArrayLike,
    polarization: Polarization,
    method: str,
    *,
    atol: float = BAND_ABS_R_TOLERANCE,
) -> DispersionScan:
    """Empaqueta R(ν) de cualquier método en un DispersionScan."""
    nu = np.asarray(nu_ghz, dtype=np.float64)
    omega = omega_from_nu_ghz(nu)
    r_arr = np.asarray(r, dtype=np.complex128)
    allowed, n_clipped = allowed_mask(r_arr, atol=atol)
    allowed = allowed & ~near_drude_pole(spec, omega)
    return DispersionScan(
        nu_ghz=nu,
        omega=omega,
        r=r_arr,
        abs_r=np.abs(r_arr),
        allowed=allowed,
        k_lm_over_pi=reduced_wavevector(r_arr, allowed),
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
