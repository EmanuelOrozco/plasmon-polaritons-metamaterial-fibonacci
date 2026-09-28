"""Medio efectivo de la celda S_m: ⟨ε⟩_m, ⟨μ⟩_m y el gap ⟨n⟩ = 0.

Promedios ponderados por espesor sobre la celda,

    ⟨ε⟩_m = (N_A ε_A a + N_B ε_B b) / L_m,   ⟨μ⟩_m = (N_A μ_A a + N_B μ_B b) / L_m.

El gap de índice promedio nulo aparece donde ⟨ε⟩_m y ⟨μ⟩_m cambian de signo.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from fibonacci_photonics.physics.fibonacci import cell_length, n_layers_a, n_layers_b, n_layers_total
from fibonacci_photonics.physics.superlattice import SuperlatticeSpec


def average_epsilon_mu(
    spec: SuperlatticeSpec,
    m: int,
    omega: ArrayLike,
) -> tuple[NDArray[np.complex128], NDArray[np.complex128]]:
    """⟨ε⟩_m y ⟨μ⟩_m (adimensionales) a la frecuencia angular ``omega`` (rad/s)."""
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
    """Frecuencia (GHz) donde ⟨ε⟩_m = ⟨μ⟩_m = 0.

    Con ω_e = ω_m, A = aire y a = b, ambos promedios se anulan juntos en
    ν = ν_p √(F_{m−2}/F_m). Devuelve ``None`` si no se cumplen esas hipótesis.
    """
    equal_plasma = np.isclose(spec.omega_e, spec.omega_m)
    equal_thickness = np.isclose(spec.thickness_a, spec.thickness_b)
    air = np.isclose(float(np.real(spec.medium_a.epsilon)), 1.0) and np.isclose(
        float(np.real(spec.medium_a.mu)), 1.0
    )
    if not (equal_plasma and equal_thickness and air):
        return None
    nu_p = spec.nu_e_ghz()
    return nu_p * float(np.sqrt(n_layers_b(m) / n_layers_total(m)))
