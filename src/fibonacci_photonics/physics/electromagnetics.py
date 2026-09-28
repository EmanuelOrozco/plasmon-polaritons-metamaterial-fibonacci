"""Cantidades electromagnéticas para incidencia oblicua TE/TM.

Con la superred apilada en z y el plano de incidencia xz, la componente
tangencial del vector de onda se conserva en todas las capas:

    q = (ω/c) n_A sen θ,                         [DEDUC + REF16]

y en la capa j el vector de onda longitudinal es

    Q_j = √((ω/c)² ε_j μ_j − q²),

con la rama principal de la raíz compleja. El paper escribe q = n_A sen θ,
omitiendo el factor ω/c (notación incompleta); aquí se restaura.

El parámetro característico χ que entra en la matriz de capa es χ = μ en TE
(E ∥ y) y χ = ε en TM (H ∥ y); ψ es el otro: ψ = ε en TE y ψ = μ en TM.
"""

from __future__ import annotations

from enum import Enum

import numpy as np
from numpy.typing import ArrayLike, NDArray

from fibonacci_photonics.errors import InvalidParameterError
from fibonacci_photonics.physics.constants import SPEED_OF_LIGHT


class Polarization(str, Enum):
    """Polarización de la onda: TE (E ∥ y) o TM (H ∥ y)."""

    TE = "TE"
    TM = "TM"


def as_polarization(value: Polarization | str) -> Polarization:
    """Convierte ``"TE"``/``"te"``/``Polarization.TE`` en ``Polarization``.

    Raises
    ------
    InvalidParameterError
        Si el nombre no es TE ni TM.
    """
    if isinstance(value, Polarization):
        return value
    try:
        return Polarization[str(value).strip().upper()]
    except KeyError as exc:
        raise InvalidParameterError(f"polarización desconocida: {value!r} (use 'TE' o 'TM')") from exc


def refractive_index(epsilon: ArrayLike, mu: ArrayLike) -> NDArray[np.complex128]:
    """n = √μ √ε con las raíces principales (convención del paper)."""
    eps = np.asarray(epsilon, dtype=np.complex128)
    mu_arr = np.asarray(mu, dtype=np.complex128)
    return np.sqrt(mu_arr) * np.sqrt(eps)


def n_squared(epsilon: ArrayLike, mu: ArrayLike) -> NDArray[np.complex128]:
    """n² = ε μ. Es lo que entra en Q y evita elegir la rama del índice."""
    return np.asarray(epsilon, dtype=np.complex128) * np.asarray(mu, dtype=np.complex128)


def in_plane_wavevector(
    omega: ArrayLike,
    theta: float,
    n_incident: ArrayLike = 1.0,
    speed_of_light: float = SPEED_OF_LIGHT,
) -> NDArray[np.complex128]:
    """Componente tangencial q = (ω/c) n_A sen θ, en 1/m.

    Parameters
    ----------
    omega : array_like
        Frecuencia angular en rad/s.
    theta : float
        Ángulo de incidencia en el medio A, en rad.
    n_incident : array_like
        Índice del medio de incidencia n_A.
    speed_of_light : float
        c en m/s.
    """
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
    """Q = √((ω/c)² ε μ − q²) en 1/m, con la rama principal compleja.

    R_m es par en cada Q_j, así que la rama elegida no cambia la semitraza.
    """
    omega_arr = np.asarray(omega, dtype=np.float64)
    k0 = omega_arr / speed_of_light
    q_arr = np.asarray(q, dtype=np.complex128)
    return np.sqrt(k0**2 * n_squared(epsilon, mu) - q_arr**2)


def characteristic_parameter(
    polarization: Polarization, epsilon: ArrayLike, mu: ArrayLike
) -> NDArray[np.complex128]:
    """χ = μ en TE y χ = ε en TM (el factor de impedancia de la matriz de capa)."""
    if polarization is Polarization.TE:
        return np.asarray(mu, dtype=np.complex128)
    if polarization is Polarization.TM:
        return np.asarray(epsilon, dtype=np.complex128)
    raise InvalidParameterError(f"polarización desconocida: {polarization}")


def chi_psi(polarization: Polarization, epsilon: complex, mu: complex) -> tuple[complex, complex]:
    """(χ, ψ) de una capa: (μ, ε) en TE y (ε, μ) en TM.

    Son los coeficientes de la ecuación de onda escalar
    d/dz[(1/χ) du/dz] + ((ω/c)² ψ − q²/χ) u = 0, con u = E_y (TE) o H_y (TM).
    """
    if polarization is Polarization.TE:
        return complex(mu), complex(epsilon)
    if polarization is Polarization.TM:
        return complex(epsilon), complex(mu)
    raise InvalidParameterError(f"polarización desconocida: {polarization}")
