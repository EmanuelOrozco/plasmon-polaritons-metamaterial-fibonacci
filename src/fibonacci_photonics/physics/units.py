"""Conversiones de unidades: única fuente para GHz ↔ rad/s y mm ↔ m."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from fibonacci_photonics.physics.constants import GHZ, MM


def omega_from_nu_ghz(nu_ghz: ArrayLike) -> NDArray[np.float64]:
    """Frecuencia angular ω = 2π ν.

    Parameters
    ----------
    nu_ghz : array_like
        Frecuencia ν en GHz.

    Returns
    -------
    ndarray
        ω en rad/s, con la forma de ``nu_ghz``.
    """
    return 2.0 * np.pi * np.asarray(nu_ghz, dtype=np.float64) * GHZ


def nu_ghz_from_omega(omega: ArrayLike) -> NDArray[np.float64]:
    """Frecuencia ν = ω / 2π en GHz a partir de ω en rad/s."""
    return np.asarray(omega, dtype=np.float64) / (2.0 * np.pi * GHZ)


def mm_to_m(length_mm: float) -> float:
    """Longitud en metros a partir de milímetros."""
    return length_mm * MM
