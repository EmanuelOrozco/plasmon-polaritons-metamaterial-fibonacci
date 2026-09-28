"""PWE en la forma k(ω): vectores de Bloch de la superred para ω dada.

La componente transversal u(z) (E_y en TE, H_y en TM) cumple

    d/dz[(1/χ) du/dz] + (κ² ψ − q²/χ) u = 0,   κ = ω/c,  q = κ n_A sen θ,

con (χ, ψ) = (μ, ε) en TE y (ε, μ) en TM. Con u(z) = Σ_n u_n e^{i(k+G_n)z} y
K = diag(k + G_n) queda

    K P K u = F u,   F = κ² [[ψ]] − q² [[1/χ]],

donde P representa 1/χ. La regla de Laurent usa P = [[1/χ]]; la regla inversa
de Li usa P = [[χ]]^{−1}, correcta porque el flujo (1/χ) du/dz es continuo en
las interfaces. Desarrollando K = k I + G sale un problema cuadrático en k,

    k² P + k (P G + G P) + (G P G − F) = 0,

que se resuelve linealizado en forma compañera de tamaño 2N. A diferencia de
la forma ω(k) del libro, ε(ω) y μ(ω) de Drude entran como números, sin
aumentar el grado del problema, y no aparecen modos espurios.

Todo se adimensionaliza con L = L_m: k̃ = kL, G̃_n = 2πn, κ̃ = κL.
Costo por frecuencia O((2N)³) en tiempo y O((2N)²) en memoria.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray

from fibonacci_photonics.errors import InvalidParameterError
from fibonacci_photonics.parallel import starmap_processes
from fibonacci_photonics.physics.electromagnetics import Polarization, chi_psi
from fibonacci_photonics.physics.fibonacci import n_layers_total
from fibonacci_photonics.physics.superlattice import SuperlatticeSpec
from fibonacci_photonics.solvers.base import validate_positive_int
from fibonacci_photonics.solvers.pwe.fourier import CellFourier, cell_fourier

FactorizationRule = Literal["inverse", "laurent"]
FACTORIZATION_RULES: tuple[str, ...] = ("inverse", "laurent")

BRILLOUIN_WINDOW = 1.25
"""Se buscan autovalores con |Re k̃| ≤ 1.25π (zona reducida más un margen)."""

CHUNKS_PER_WORKER = 4
"""Trozos del barrido por proceso: equilibra la carga sin multiplicar el costo de arranque."""


@dataclass(frozen=True)
class BlochSolution:
    """Resultado por frecuencia.

    Attributes
    ----------
    k_tilde : ndarray
        k̃ = k L_m del modo físico (complejo en gaps).
    r : ndarray
        R = Re cos k̃.
    decay : ndarray
        |Im k̃|, tasa de decaimiento por celda (0 en bandas).
    """

    k_tilde: NDArray[np.complex128]
    r: NDArray[np.float64]
    decay: NDArray[np.float64]


def validate_rule(rule: str) -> FactorizationRule:
    """``"inverse"`` (Li) o ``"laurent"``."""
    if rule not in FACTORIZATION_RULES:
        raise InvalidParameterError(f"regla de factorización desconocida: {rule!r}")
    return rule  # type: ignore[return-value]


def n_max_for(m: int, harmonics_per_layer: int) -> int:
    """Truncamiento n = −n_max … n_max con n_max = h · F_m (proporcional al número de capas)."""
    return validate_positive_int(harmonics_per_layer, "harmonics_per_layer") * n_layers_total(m)


def _layer_coefficients(
    spec: SuperlatticeSpec, omega: float, polarization: Polarization
) -> tuple[tuple[complex, complex], tuple[complex, complex], complex]:
    """((χ_A, χ_B), (ψ_A, ψ_B), ε_A μ_A) a la frecuencia ω (rad/s)."""
    eps_a, mu_a = (complex(v) for v in spec.medium_a.epsilon_mu(omega))
    eps_b, mu_b = (complex(v) for v in spec.medium_b.epsilon_mu(omega))
    chi_a, psi_a = chi_psi(polarization, eps_a, mu_a)
    chi_b, psi_b = chi_psi(polarization, eps_b, mu_b)
    return (chi_a, chi_b), (psi_a, psi_b), eps_a * mu_a


def qep_matrices(
    spec: SuperlatticeSpec,
    cell: CellFourier,
    omega: float,
    theta: float,
    polarization: Polarization = Polarization.TE,
    rule: FactorizationRule = "inverse",
) -> tuple[NDArray[np.complex128], NDArray[np.complex128], NDArray[np.complex128], NDArray[np.complex128]]:
    """(P, A1, A0, P^{−1}) del problema cuadrático k̃² P + k̃ A1 + A0 = 0."""
    rule = validate_rule(rule)
    (chi_a, chi_b), (psi_a, psi_b), n2_a = _layer_coefficients(spec, omega, polarization)
    kappa = omega / spec.speed_of_light * cell.period
    q2 = (kappa**2) * n2_a * np.sin(theta) ** 2
    g = np.diag(2.0 * np.pi * np.arange(-cell.n_max, cell.n_max + 1)).astype(np.complex128)

    t_chi = cell.toeplitz(chi_a, chi_b)
    t_inv_chi = cell.toeplitz(1.0 / chi_a, 1.0 / chi_b)
    if rule == "inverse":
        p_inv = t_chi
        p = np.linalg.inv(t_chi)
    else:
        p = t_inv_chi
        p_inv = np.linalg.inv(t_inv_chi)

    f = kappa**2 * cell.toeplitz(psi_a, psi_b) - q2 * t_inv_chi
    a1 = p @ g + g @ p
    a0 = g @ p @ g - f
    return p, a1, a0, p_inv


def bloch_eigenvalues(
    spec: SuperlatticeSpec,
    cell: CellFourier,
    omega: float,
    theta: float,
    polarization: Polarization = Polarization.TE,
    rule: FactorizationRule = "inverse",
) -> NDArray[np.complex128]:
    """Los 2N autovalores k̃ del problema linealizado [[0, I], [−P⁻¹A0, −P⁻¹A1]]."""
    _, a1, a0, p_inv = qep_matrices(spec, cell, omega, theta, polarization, rule)
    size = a1.shape[0]
    companion = np.zeros((2 * size, 2 * size), dtype=np.complex128)
    companion[:size, size:] = np.eye(size)
    companion[size:, :size] = -p_inv @ a0
    companion[size:, size:] = -p_inv @ a1
    return np.linalg.eigvals(companion)


def select_bloch_wavevector(eigenvalues: NDArray[np.complex128]) -> complex:
    """Modo físico: el de menor |Im k̃| con |Re k̃| ≤ 1.25π (si no hay, entre todos)."""
    inside = eigenvalues[np.abs(eigenvalues.real) <= BRILLOUIN_WINDOW * np.pi]
    if inside.size == 0:
        inside = eigenvalues
    return complex(inside[np.argmin(np.abs(inside.imag))])


def _solve_chunk(
    spec: SuperlatticeSpec,
    cell: CellFourier,
    omegas: NDArray[np.float64],
    theta: float,
    polarization: Polarization,
    rule: FactorizationRule,
) -> NDArray[np.complex128]:
    out = np.full(omegas.shape, np.nan + 1j * np.nan, dtype=np.complex128)
    with np.errstate(divide="ignore", invalid="ignore"):
        for i, omega in enumerate(omegas):
            try:
                eigenvalues = bloch_eigenvalues(spec, cell, float(omega), theta, polarization, rule)
            except (np.linalg.LinAlgError, ValueError):
                continue
            finite = eigenvalues[np.isfinite(eigenvalues)]
            if finite.size:
                out[i] = select_bloch_wavevector(finite)
    return out


def solve_bloch(
    spec: SuperlatticeSpec,
    m: int,
    omega: ArrayLike,
    theta: float,
    polarization: Polarization = Polarization.TE,
    *,
    n_max: int,
    rule: FactorizationRule = "inverse",
    workers: int = 1,
) -> BlochSolution:
    """k̃, R = Re cos k̃ y |Im k̃| para cada ω (rad/s), repartido en ``workers`` procesos.

    Las frecuencias donde el problema es singular (polos de Drude) quedan en NaN.
    """
    rule = validate_rule(rule)
    workers = validate_positive_int(workers, "workers")
    omegas = np.atleast_1d(np.asarray(omega, dtype=np.float64))
    cell = cell_fourier(m, spec.thickness_a, spec.thickness_b, n_max)
    if workers == 1 or omegas.size < 2 * workers:
        k_tilde = _solve_chunk(spec, cell, omegas, theta, polarization, rule)
    else:
        chunks = np.array_split(omegas, CHUNKS_PER_WORKER * workers)
        parts = starmap_processes(
            _solve_chunk,
            [(spec, cell, chunk, theta, polarization, rule) for chunk in chunks],
            workers,
        )
        k_tilde = np.concatenate(parts)
    return BlochSolution(k_tilde=k_tilde, r=np.cos(k_tilde).real, decay=np.abs(k_tilde.imag))
