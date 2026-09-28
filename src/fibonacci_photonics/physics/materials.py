"""Respuestas electromagnéticas de los medios A y B.

Medio A: homogéneo y no dispersivo; en el paper es aire, ε_A = μ_A = 1.

Medio B: metamaterial de Drude sin pérdidas (ec. 10 del paper),

.. math::

    \\varepsilon_B(\\omega) = \\varepsilon_0 - \\frac{\\omega_e^2}{\\omega^2}, \\qquad
    \\mu_B(\\omega) = \\mu_0 - \\frac{\\omega_m^2}{\\omega^2},

con ε_0 = μ_0 = 1 en todos los resultados del paper y sin amortiguamiento γ.
Las permitividades y permeabilidades son relativas (adimensionales).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from fibonacci_photonics.errors import InvalidParameterError


@dataclass(frozen=True)
class HomogeneousMedium:
    """Medio no dispersivo.

    Attributes
    ----------
    epsilon, mu : complex
        Permitividad y permeabilidad relativas; ninguna puede ser nula.
    label : str
        Etiqueta de la capa en la palabra de Fibonacci.
    """

    epsilon: complex
    mu: complex
    label: str = "A"

    def __post_init__(self) -> None:
        for name in ("epsilon", "mu"):
            value = complex(getattr(self, name))
            if not (math.isfinite(value.real) and math.isfinite(value.imag)) or value == 0:
                raise InvalidParameterError(f"{name} del medio {self.label} debe ser finito y no nulo: {value}")

    def epsilon_mu(self, omega: ArrayLike) -> tuple[NDArray[np.complex128], NDArray[np.complex128]]:
        """ε y μ constantes, con la forma de ``omega`` (rad/s)."""
        omega_arr = np.asarray(omega, dtype=np.float64)
        eps = np.full(omega_arr.shape, self.epsilon, dtype=np.complex128)
        mu = np.full(omega_arr.shape, self.mu, dtype=np.complex128)
        return eps, mu


@dataclass(frozen=True)
class DrudeMetamaterial:
    """Modelo de Drude sin pérdidas del paper, ec. (10).

    Attributes
    ----------
    omega_e, omega_m : float
        Frecuencias de plasma eléctrica y magnética en rad/s (≥ 0).
    epsilon_0, mu_0 : float
        Fondos ε_0 y μ_0 (> 0). Con ω_e = ω_m = 0 el medio es un dieléctrico
        no dispersivo, como la capa ε_2 del benchmark del libro.
    label : str
        Etiqueta de la capa en la palabra de Fibonacci.
    """

    omega_e: float
    omega_m: float
    epsilon_0: float = 1.0
    mu_0: float = 1.0
    label: str = "B"

    def __post_init__(self) -> None:
        for name in ("omega_e", "omega_m"):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value < 0.0:
                raise InvalidParameterError(f"{name} debe ser finita y ≥ 0 (rad/s): {value}")
        for name in ("epsilon_0", "mu_0"):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0:
                raise InvalidParameterError(f"{name} debe ser finito y > 0: {value}")

    def epsilon(self, omega: ArrayLike) -> NDArray[np.complex128]:
        """ε_B(ω) = ε_0 − ω_e²/ω², con ω en rad/s.

        Raises
        ------
        ZeroDivisionError
            Si algún ω es 0 con omega_e > 0: el modelo diverge.
        """
        omega_arr = np.asarray(omega, dtype=np.float64)
        if self.omega_e == 0.0:
            return np.full(omega_arr.shape, self.epsilon_0, dtype=np.complex128)
        if np.any(omega_arr == 0.0):
            raise ZeroDivisionError("el modelo de Drude diverge en ω = 0")
        return np.asarray(self.epsilon_0 - (self.omega_e**2) / omega_arr**2, dtype=np.complex128)

    def mu(self, omega: ArrayLike) -> NDArray[np.complex128]:
        """μ_B(ω) = μ_0 − ω_m²/ω², con ω en rad/s.

        Raises
        ------
        ZeroDivisionError
            Si algún ω es 0 con omega_m > 0: el modelo diverge.
        """
        omega_arr = np.asarray(omega, dtype=np.float64)
        if self.omega_m == 0.0:
            return np.full(omega_arr.shape, self.mu_0, dtype=np.complex128)
        if np.any(omega_arr == 0.0):
            raise ZeroDivisionError("el modelo de Drude diverge en ω = 0")
        return np.asarray(self.mu_0 - (self.omega_m**2) / omega_arr**2, dtype=np.complex128)

    def epsilon_mu(self, omega: ArrayLike) -> tuple[NDArray[np.complex128], NDArray[np.complex128]]:
        """(ε_B(ω), μ_B(ω)) con ω en rad/s."""
        return self.epsilon(omega), self.mu(omega)

    def electric_plasmon_frequency_hz(self) -> float:
        """ν_e = ω_e / (2π √ε_0) en Hz: ε_B(2π ν_e) = 0."""
        return self.omega_e / (2.0 * np.pi * np.sqrt(self.epsilon_0))

    def magnetic_plasmon_frequency_hz(self) -> float:
        """ν_m = ω_m / (2π √μ_0) en Hz: μ_B(2π ν_m) = 0 (plasmón magnético)."""
        return self.omega_m / (2.0 * np.pi * np.sqrt(self.mu_0))
