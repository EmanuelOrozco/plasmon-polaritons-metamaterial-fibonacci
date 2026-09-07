"""Respuestas electromagnéticas de los medios A y B.

Medio A: aire, ε = μ = 1 (paper).
Medio B: Drude sin pérdidas (eq. 10 del paper):

    ε_B(ω) = ε0 - ω_e² / ω²
    μ_B(ω) = μ0 - ω_m² / ω²

El paper usa ε0 = μ0 = 1 en todos los resultados numéricos.
No hay término de amortiguamiento γ.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray


@dataclass(frozen=True)
class HomogeneousMedium:
    epsilon: complex
    mu: complex
    label: str = "A"

    def epsilon_mu(self, omega: ArrayLike) -> tuple[NDArray[np.complex128], NDArray[np.complex128]]:
        omega_arr = np.asarray(omega, dtype=np.float64)
        eps = np.full(omega_arr.shape, self.epsilon, dtype=np.complex128)
        mu = np.full(omega_arr.shape, self.mu, dtype=np.complex128)
        return eps, mu


@dataclass(frozen=True)
class DrudeMetamaterial:
    """Modelo de Drude del paper, eq. (10). Frecuencias angulares en rad/s."""

    omega_e: float
    omega_m: float
    epsilon_0: float = 1.0
    mu_0: float = 1.0
    label: str = "B"

    def epsilon(self, omega: ArrayLike) -> NDArray[np.complex128]:
        omega_arr = np.asarray(omega, dtype=np.float64)
        if np.any(omega_arr == 0.0):
            raise ZeroDivisionError("el modelo de Drude diverge en ω = 0")
        return np.asarray(self.epsilon_0 - (self.omega_e**2) / omega_arr**2, dtype=np.complex128)

    def mu(self, omega: ArrayLike) -> NDArray[np.complex128]:
        omega_arr = np.asarray(omega, dtype=np.float64)
        if np.any(omega_arr == 0.0):
            raise ZeroDivisionError("el modelo de Drude diverge en ω = 0")
        return np.asarray(self.mu_0 - (self.omega_m**2) / omega_arr**2, dtype=np.complex128)

    def epsilon_mu(self, omega: ArrayLike) -> tuple[NDArray[np.complex128], NDArray[np.complex128]]:
        return self.epsilon(omega), self.mu(omega)

    def electric_plasmon_frequency_hz(self) -> float:
        """ν_e = ω_e / (2π √ε0)."""
        return self.omega_e / (2.0 * np.pi * np.sqrt(self.epsilon_0))

    def magnetic_plasmon_frequency_hz(self) -> float:
        """ν_m = ω_m / (2π √μ0)."""
        return self.omega_m / (2.0 * np.pi * np.sqrt(self.mu_0))
