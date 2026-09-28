"""PWE en la forma k(ω): vectores de Bloch de la superred para ω dada.

La ecuación de onda para la componente transversal u(z) (E_y en TE, H_y en TM),

    d/dz[(1/χ) du/dz] + (κ² ψ − q²/χ) u = 0,   κ = ω/c,  q = κ n_A sen θ,

con (χ, ψ) = (μ, ε) en TE y (ε, μ) en TM, se expande en ondas planas
u(z) = Σ_n u_n e^{i(k+G_n)z}. Con K = diag(k + G_n) queda

    K P K u = F u,   F = κ² [[ψ]] − q² [[1/χ]],

donde P es la representación de 1/χ. La regla de Laurent usa P = [[1/χ]];
la regla inversa de Li usa P = [[χ]]^{-1}, que es la correcta porque el flujo
(1/χ) du/dz es continuo en las interfaces. Desarrollando K = k I + G sale un
problema cuadrático de autovalores en k,

    k² P + k (P G + G P) + (G P G − F) = 0,

que se resuelve linealizado en forma compañera (tamaño 2N). A diferencia de
la forma ω(k) del libro, aquí ε(ω) y μ(ω) de Drude entran como números, sin
aumentar el grado del problema, y no aparecen modos espurios.

Todo se adimensionaliza con L = Lm: k̃ = kL, G̃_n = 2πn, κ̃ = κL.
"""

from __future__ import annotations

import os
from collections.abc import Sequence
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray

from fibonacci_photonics.core.electromagnetics import Polarization
from fibonacci_photonics.core.fibonacci import n_layers_total
from fibonacci_photonics.core.params import SuperlatticeSpec
from fibonacci_photonics.pwe.fourier import CellFourier, cell_fourier

FactorizationRule = Literal["inverse", "laurent"]

BRILLOUIN_WINDOW = 1.25
"""Se buscan autovalores con |Re k̃| ≤ 1.25 π (zona reducida más un margen)."""

_BLAS_ENV = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")


@dataclass(frozen=True)
class BlochSolution:
    """Resultado por frecuencia: k̃ elegido, R = Re cos k̃ y |Im k̃| (decaimiento)."""

    k_tilde: NDArray[np.complex128]
    r: NDArray[np.float64]
    decay: NDArray[np.float64]


def n_max_for(m: int, harmonics_per_layer: int) -> int:
    """Truncamiento n = −n_max … n_max con n_max proporcional al número de capas."""
    return int(harmonics_per_layer) * n_layers_total(m)


def _responses(spec: SuperlatticeSpec, omega: float, polarization: Polarization):
    eps_a, mu_a = (complex(v) for v in spec.medium_a.epsilon_mu(omega))
    eps_b, mu_b = (complex(v) for v in spec.medium_b.epsilon_mu(omega))
    if polarization == Polarization.TE:
        return (mu_a, mu_b), (eps_a, eps_b), (eps_a, mu_a)
    return (eps_a, eps_b), (mu_a, mu_b), (eps_a, mu_a)


def qep_matrices(
    spec: SuperlatticeSpec,
    cell: CellFourier,
    omega: float,
    theta: float,
    polarization: Polarization = Polarization.TE,
    rule: FactorizationRule = "inverse",
) -> tuple[NDArray[np.complex128], NDArray[np.complex128], NDArray[np.complex128], NDArray[np.complex128]]:
    """Devuelve (P, A1, A0, P^{-1}) del QEP k̃² P + k̃ A1 + A0 = 0."""
    (chi_a, chi_b), (psi_a, psi_b), (eps_a, mu_a) = _responses(spec, omega, polarization)
    kappa = omega / spec.speed_of_light * cell.period
    q2 = (kappa**2) * complex(eps_a * mu_a) * np.sin(theta) ** 2
    g = np.diag(2.0 * np.pi * np.arange(-cell.n_max, cell.n_max + 1)).astype(np.complex128)

    t_chi = cell.toeplitz(chi_a, chi_b)
    t_inv_chi = cell.toeplitz(1.0 / chi_a, 1.0 / chi_b)
    if rule == "inverse":
        p_inv = t_chi
        p = np.linalg.inv(t_chi)
    elif rule == "laurent":
        p = t_inv_chi
        p_inv = np.linalg.inv(t_inv_chi)
    else:
        raise ValueError(f"regla de factorización desconocida: {rule!r}")

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
    """Los 2N autovalores k̃ = k L del QEP linealizado."""
    _, a1, a0, p_inv = qep_matrices(spec, cell, omega, theta, polarization, rule)
    size = a1.shape[0]
    companion = np.zeros((2 * size, 2 * size), dtype=np.complex128)
    companion[:size, size:] = np.eye(size)
    companion[size:, :size] = -p_inv @ a0
    companion[size:, size:] = -p_inv @ a1
    return np.linalg.eigvals(companion)


def select_bloch_wavevector(eigenvalues: NDArray[np.complex128]) -> complex:
    """El modo físico: el de menor |Im k̃| dentro de la zona reducida ampliada."""
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


def _single_threaded_blas() -> None:
    # Los procesos hijos heredan el entorno: un hilo de BLAS por proceso evita
    # la sobre-suscripción cuando se reparte el barrido en varios núcleos.
    for name in _BLAS_ENV:
        os.environ.setdefault(name, "1")


def default_workers() -> int:
    return max(1, (os.cpu_count() or 1))


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
    """Resuelve el QEP para cada ω y devuelve k̃, R = Re cos k̃ y |Im k̃|."""
    omegas = np.atleast_1d(np.asarray(omega, dtype=np.float64))
    cell = cell_fourier(m, spec.thickness_a, spec.thickness_b, n_max)
    if workers <= 1 or omegas.size < 2 * workers:
        k_tilde = _solve_chunk(spec, cell, omegas, theta, polarization, rule)
    else:
        _single_threaded_blas()
        chunks: Sequence[NDArray[np.float64]] = np.array_split(omegas, 4 * workers)
        with ProcessPoolExecutor(max_workers=workers) as pool:
            parts = pool.map(
                _solve_chunk,
                [spec] * len(chunks),
                [cell] * len(chunks),
                chunks,
                [theta] * len(chunks),
                [polarization] * len(chunks),
                [rule] * len(chunks),
            )
            k_tilde = np.concatenate(list(parts))
    return BlochSolution(
        k_tilde=k_tilde,
        r=np.cos(k_tilde).real,
        decay=np.abs(k_tilde.imag),
    )
