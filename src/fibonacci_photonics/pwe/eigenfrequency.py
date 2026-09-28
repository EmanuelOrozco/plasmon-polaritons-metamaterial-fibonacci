"""PWE en la forma ω(k) del libro: frecuencias propias para k dado.

1. Medios no dispersivos (benchmark del capítulo 1D de Sukhoivanov y Guryev).
   Con q fijo, la ecuación de onda da el problema generalizado hermítico

       (K P K + q² [[1/χ]]) u = κ² [[ψ]] u,

   que con χ = ε, ψ = 1, q = 0 y regla de Laurent (P = [[1/ε]]) es exactamente
   la "master equation" (4.37)/(4.44) del libro: Σ χ(G−G')(k+G)(k+G') h = (ω/c)² h.

2. Metamaterial de Drude (TE). Como μ_B y ε_B dependen de λ = κ², el problema
   deja de ser lineal. Multiplicando por D = μ0 λ − λm se obtiene un problema
   cuadrático en λ,

       (λ² M2 + λ M1 + M0) u = 0,  M2 = [[r2]],  M1 = [[r1]] − K[[p1]]K,
                                    M0 = [[r0]] − K[[p0]]K,

   con p = D/μ y r = D(ε − s²/μ) λ desarrollados por potencias de λ. Esta forma
   sirve para mostrar por qué no se usa en la superred: los coeficientes
   escalonados de D/μ se desarrollan con la regla de Laurent y aparecen modos
   espurios y acumulación de autovalores en ν_m (contaminación espectral).

Longitudes adimensionalizadas con el período L: k̃ = kL, κ̃ = ωL/c.
"""

from __future__ import annotations

import numpy as np
import scipy.linalg
from numpy.typing import ArrayLike, NDArray

from fibonacci_photonics.core.electromagnetics import Polarization
from fibonacci_photonics.core.params import SuperlatticeSpec
from fibonacci_photonics.pwe.bloch import FactorizationRule
from fibonacci_photonics.pwe.fourier import CellFourier, cell_fourier


def _hermitian(matrix: NDArray[np.complex128]) -> NDArray[np.complex128]:
    return 0.5 * (matrix + matrix.conj().T)


def nondispersive_operator(
    cell: CellFourier,
    epsilon: tuple[complex, complex],
    mu: tuple[complex, complex],
    k_tilde: float,
    q_tilde: float = 0.0,
    polarization: Polarization = Polarization.TM,
    rule: FactorizationRule = "laurent",
) -> tuple[NDArray[np.complex128], NDArray[np.complex128]]:
    """Matrices (A, B) del problema A u = κ̃² B u para k̃ y q̃ = qL dados."""
    if polarization == Polarization.TE:
        chi, psi = mu, epsilon
    else:
        chi, psi = epsilon, mu
    g = 2.0 * np.pi * np.arange(-cell.n_max, cell.n_max + 1)
    kg = np.diag(k_tilde + g).astype(np.complex128)
    t_inv_chi = cell.toeplitz(1.0 / chi[0], 1.0 / chi[1])
    if rule == "laurent":
        p = t_inv_chi
    elif rule == "inverse":
        p = np.linalg.inv(cell.toeplitz(chi[0], chi[1]))
    else:
        raise ValueError(f"regla de factorización desconocida: {rule!r}")
    a = kg @ p @ kg + (q_tilde**2) * t_inv_chi
    b = cell.toeplitz(psi[0], psi[1])
    return _hermitian(a), _hermitian(b)


def nondispersive_bands(
    cell: CellFourier,
    epsilon: tuple[complex, complex],
    mu: tuple[complex, complex],
    k_tilde: ArrayLike,
    *,
    q_tilde: float = 0.0,
    polarization: Polarization = Polarization.TM,
    rule: FactorizationRule = "laurent",
    n_bands: int = 8,
    return_modes: bool = False,
):
    """Frecuencias normalizadas ωL/(2πc) de las primeras ``n_bands`` bandas.

    Con ``return_modes`` también devuelve los coeficientes u_n de cada modo.
    """
    ks = np.atleast_1d(np.asarray(k_tilde, dtype=np.float64))
    freqs = np.empty((ks.size, n_bands))
    modes = []
    for i, k in enumerate(ks):
        a, b = nondispersive_operator(cell, epsilon, mu, float(k), q_tilde, polarization, rule)
        values, vectors = scipy.linalg.eigh(a, b, subset_by_index=(0, n_bands - 1))
        freqs[i] = np.sqrt(np.clip(values, 0.0, None)) / (2.0 * np.pi)
        if return_modes:
            modes.append(vectors)
    if return_modes:
        return freqs, modes
    return freqs


def mode_profile(
    cell: CellFourier,
    coefficients: NDArray[np.complex128],
    k_tilde: float,
    z_over_period: NDArray[np.float64],
) -> NDArray[np.complex128]:
    """u(z) = Σ_n u_n e^{i(k̃ + 2πn) z/L}."""
    g = 2.0 * np.pi * np.arange(-cell.n_max, cell.n_max + 1)
    return np.exp(1j * np.outer(z_over_period, k_tilde + g)) @ coefficients


def drude_qep_matrices(
    spec: SuperlatticeSpec,
    cell: CellFourier,
    k_tilde: float,
    theta: float,
) -> tuple[NDArray[np.complex128], NDArray[np.complex128], NDArray[np.complex128]]:
    """(M2, M1, M0) del QEP en λ̃ = κ̃² para la superred aire/Drude en TE."""
    length = cell.period
    c = spec.speed_of_light
    eps_a = complex(spec.medium_a.epsilon)
    mu_a = complex(spec.medium_a.mu)
    eps0 = spec.medium_b.epsilon_0
    mu0 = spec.medium_b.mu_0
    lam_e = (spec.medium_b.omega_e * length / c) ** 2
    lam_m = (spec.medium_b.omega_m * length / c) ** 2
    s2 = eps_a * mu_a * np.sin(theta) ** 2

    def t(value_a: complex, value_b: complex) -> NDArray[np.complex128]:
        return cell.toeplitz(value_a, value_b)

    alpha_a = eps_a - s2 / mu_a
    p1 = t(mu0 / mu_a, 1.0)
    p0 = t(-lam_m / mu_a, 0.0)
    r2 = t(alpha_a * mu0, eps0 * mu0 - s2)
    r1 = t(-alpha_a * lam_m, -(eps0 * lam_m + mu0 * lam_e))
    r0 = t(0.0, lam_e * lam_m)

    g = 2.0 * np.pi * np.arange(-cell.n_max, cell.n_max + 1)
    kg = np.diag(k_tilde + g).astype(np.complex128)
    return r2, r1 - kg @ p1 @ kg, r0 - kg @ p0 @ kg


def drude_qep_frequencies_ghz(
    spec: SuperlatticeSpec,
    m: int,
    k_tilde: float,
    theta: float,
    n_max: int,
    *,
    imag_tolerance: float = 1.0e-6,
) -> NDArray[np.float64]:
    """Frecuencias reales ν (GHz) del QEP ω(k) de Drude, ordenadas.

    Se descartan autovalores con parte imaginaria relativa mayor que
    ``imag_tolerance`` o con λ̃ ≤ 0.
    """
    cell = cell_fourier(m, spec.thickness_a, spec.thickness_b, n_max)
    m2, m1, m0 = drude_qep_matrices(spec, cell, k_tilde, theta)
    size = m2.shape[0]
    eye = np.eye(size)
    zero = np.zeros((size, size))
    lhs = np.block([[zero, eye], [-m0, -m1]])
    rhs = np.block([[eye, zero], [zero, m2]])
    lam = scipy.linalg.eigvals(lhs, rhs)
    lam = lam[np.isfinite(lam)]
    real = lam[(np.abs(lam.imag) <= imag_tolerance * np.abs(lam)) & (lam.real > 0)].real
    nu_hz = np.sqrt(real) * spec.speed_of_light / (2.0 * np.pi * cell.period)
    return np.sort(nu_hz / 1.0e9)
