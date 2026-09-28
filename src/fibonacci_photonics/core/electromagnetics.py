"""Cantidades electromagnéticas para incidencia oblicua TE/TM.

q se implementa como (ω/c) n_A sin θ  [DEDUC + REF16].
El paper escribe q = n_A sin θ omitiendo ω/c (notación incompleta).
"""

from __future__ import annotations

from enum import Enum

import numpy as np
from numpy.typing import ArrayLike, NDArray

from fibonacci_photonics.core.constants import SPEED_OF_LIGHT


class Polarization(str, Enum):
    TE = "TE"
    TM = "TM"


def refractive_index(epsilon: ArrayLike, mu: ArrayLike) -> NDArray[np.complex128]:
    """n = √μ √ε con las dos raíces principales (paper)."""
    eps = np.asarray(epsilon, dtype=np.complex128)
    mu_arr = np.asarray(mu, dtype=np.complex128)
    return np.sqrt(mu_arr) * np.sqrt(eps)


def n_squared(epsilon: ArrayLike, mu: ArrayLike) -> NDArray[np.complex128]:
    """n² = με. Es lo que entra en Q; evita ramas superfluas del índice."""
    return np.asarray(epsilon, dtype=np.complex128) * np.asarray(mu, dtype=np.complex128)


def in_plane_wavevector(
    omega: ArrayLike,
    theta: float,
    n_incident: ArrayLike = 1.0,
    speed_of_light: float = SPEED_OF_LIGHT,
) -> NDArray[np.complex128]:
    """q = (ω/c) n_A sin θ. [IMPL: factor ω/c restaurado.]"""
    omega_arr = np.asarray(omega, dtype=np.float64)
    n_a = np.asarray(n_incident, dtype=np.complex128)
    return (omega_arr / speed_of_light) * n_a * np.sin(theta)


def longitudinal_wavevector(
    omega: ArrayLike,
    epsilon: ArrayLike,
    mu: ArrayLike,
    q: ArrayLike,
    speed_of_light: float = SPEED_OF_LIGHT,
) -> NDArray[np.complex128]:
    """Q = sqrt((ω/c)² n² − q²) con rama principal compleja."""
    omega_arr = np.asarray(omega, dtype=np.float64)
    k0 = omega_arr / speed_of_light
    q_arr = np.asarray(q, dtype=np.complex128)
    return np.sqrt(k0**2 * n_squared(epsilon, mu) - q_arr**2)


def characteristic_parameter(polarization: Polarization, epsilon: ArrayLike, mu: ArrayLike) -> NDArray[np.complex128]:
    """χ = μ en TE y χ = ε en TM (el factor que aparece en Ψ y en R2)."""
    if polarization is Polarization.TE:
        return np.asarray(mu, dtype=np.complex128)
    if polarization is Polarization.TM:
        return np.asarray(epsilon, dtype=np.complex128)
    raise ValueError(f"polarización desconocida: {polarization}")
