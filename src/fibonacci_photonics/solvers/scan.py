"""Barrido de dispersión independiente del método numérico.

Tanto la TMM como el PWE entregan la semitraza de la celda S_m,

    R_m(ν) = ½ Tr T_m(ν) = cos(k L_m),

donde k es el vector de Bloch de la superred. A partir de R se clasifican
bandas (|R| ≤ 1, k real) y gaps (|R| > 1, k complejo), se obtiene
k L_m/π = arccos(R)/π ∈ [0, 1] y se extraen los intervalos permitidos.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from fibonacci_photonics.physics.electromagnetics import Polarization
from fibonacci_photonics.physics.superlattice import SuperlatticeSpec
from fibonacci_photonics.physics.units import omega_from_nu_ghz

BAND_ABS_R_TOLERANCE = 1.0e-12
"""Tolerancia de |R| ≤ 1 + atol: absorbe el redondeo sin recortar en silencio. [IMPL]"""

DRUDE_POLE_RELATIVE_WIDTH = 1.0e-6
"""Semiancho relativo |ω − ω_p|/ω_p de la vecindad excluida alrededor de ε_B = 0 o μ_B = 0."""


@dataclass(frozen=True)
class FrequencyInterval:
    """Intervalo de frecuencias permitido [ν_min, ν_max] en GHz.

    ``n_points`` es el número de muestras de la malla dentro del intervalo
    (0 cuando los bordes vienen de un refinamiento sin malla asociada).
    """

    nu_min_ghz: float
    nu_max_ghz: float
    n_points: int

    @property
    def bandwidth_ghz(self) -> float:
        """Ancho Δν = ν_max − ν_min en GHz."""
        return self.nu_max_ghz - self.nu_min_ghz

    @property
    def center_ghz(self) -> float:
        """Centro (ν_min + ν_max)/2 en GHz."""
        return 0.5 * (self.nu_min_ghz + self.nu_max_ghz)


@dataclass(frozen=True)
class DispersionScan:
    """R(ν) en una malla de frecuencias y todo lo que se deriva de ella.

    Attributes
    ----------
    nu_ghz : ndarray
        Malla de frecuencias en GHz (estrictamente creciente).
    omega : ndarray
        ω = 2πν en rad/s.
    r : ndarray
        Semitraza R_m(ν), compleja (Im R ≈ 0 sin pérdidas).
    abs_r : ndarray
        |R|.
    allowed : ndarray of bool
        Banda: |R| ≤ 1 + atol y lejos de los polos de Drude.
    k_lm_over_pi : ndarray
        k L_m/π en bandas, NaN en gaps.
    method : str
        Método que produjo R (``"recurrence"``, ``"product"``, ``"pwe-inverse"``, ...).
    m : int
        Orden de Fibonacci de la celda.
    theta : float
        Ángulo de incidencia en rad.
    polarization : str
        ``"TE"`` o ``"TM"``.
    clipped_near_unit : int
        Muestras con 1 < |R| ≤ 1 + atol, tratadas como banda.
    """

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


def real_semitrace(r: NDArray[np.complex128]) -> NDArray[np.float64]:
    """Parte real de R_m. Im R debe ser ~0 en el modelo sin pérdidas."""
    return np.real(r)


def allowed_mask(
    r: NDArray[np.complex128],
    *,
    atol: float = BAND_ABS_R_TOLERANCE,
) -> tuple[NDArray[np.bool_], int]:
    """Máscara |R| ≤ 1 + atol y número de muestras en (1, 1 + atol]."""
    abs_r = np.abs(r)
    clipped = np.count_nonzero((abs_r > 1.0) & (abs_r <= 1.0 + atol))
    allowed = abs_r <= 1.0 + atol
    return allowed, int(clipped)


def reduced_wavevector(
    r: NDArray[np.complex128],
    allowed: NDArray[np.bool_],
) -> NDArray[np.float64]:
    """k L_m/π = arccos(R)/π ∈ [0, 1] en bandas; NaN en gaps."""
    r_real = np.real(r)
    r_clipped = np.clip(r_real, -1.0, 1.0)
    k_over_pi = np.full(r_real.shape, np.nan, dtype=np.float64)
    k_over_pi[allowed] = np.arccos(r_clipped[allowed]) / np.pi
    return k_over_pi


def near_drude_pole(spec: SuperlatticeSpec, omega: NDArray[np.float64]) -> NDArray[np.bool_]:
    """Frecuencias pegadas a ε_B = 0 o μ_B = 0, donde ψ o 1/χ no están definidos."""
    pole_e = spec.medium_b.omega_e
    pole_m = spec.medium_b.omega_m
    width = DRUDE_POLE_RELATIVE_WIDTH
    return (np.abs(omega - pole_e) < width * pole_e) | (np.abs(omega - pole_m) < width * pole_m)


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
    """Empaqueta R(ν) de cualquier método en un ``DispersionScan``."""
    nu = np.asarray(nu_ghz, dtype=np.float64)
    omega = omega_from_nu_ghz(nu)
    r_arr = np.asarray(r, dtype=np.complex128)
    if r_arr.shape != nu.shape:
        raise ValueError(f"R tiene forma {r_arr.shape} y la malla {nu.shape}")
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
    """Componentes conexas de muestras en banda, como intervalos de la malla."""
    padded = np.concatenate(([False], np.asarray(scan.allowed, dtype=bool), [False]))
    changes = np.flatnonzero(padded[1:] != padded[:-1])
    starts, stops = changes[0::2], changes[1::2]
    nu = scan.nu_ghz
    return [
        FrequencyInterval(nu_min_ghz=float(nu[i]), nu_max_ghz=float(nu[j - 1]), n_points=int(j - i))
        for i, j in zip(starts, stops, strict=True)
    ]
