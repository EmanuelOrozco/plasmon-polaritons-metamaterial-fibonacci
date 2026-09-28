"""Series de Fourier de perfiles escalonados en la celda S_m.

Generaliza la ec. (4.38) de Sukhoivanov y Guryev (bicapa) a una celda con
cualquier número de capas. Para f(z) constante por tramos, f = f_j en
[z_j, z_{j+1}], con período L y G_n = 2πn/L,

    f_n = (1/L) ∫_0^L f(z) e^{−i G_n z} dz,

    f_0 = Σ_j f_j d_j / L,
    f_n = Σ_j f_j · i/(G_n L) · (e^{−i G_n z_{j+1}} − e^{−i G_n z_j}),   n ≠ 0.

Los perfiles de la superred solo toman dos valores (capa A o capa B), así que
basta con los coeficientes de las funciones indicadoras 1_A(z) y 1_B(z):
f_n = f_A (1_A)_n + f_B (1_B)_n. Los coeficientes decaen como 1/|n|.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from fibonacci_photonics.errors import InvalidParameterError
from fibonacci_photonics.physics.fibonacci import LAYER_A, sequence


def _validate_n_max(n_max: int) -> int:
    if isinstance(n_max, bool) or int(n_max) != n_max or int(n_max) < 0:
        raise InvalidParameterError(f"n_max debe ser un entero ≥ 0: {n_max!r}")
    return int(n_max)


@dataclass(frozen=True)
class CellFourier:
    """Coeficientes de 1_A y 1_B de la celda para los órdenes n = −2 n_max … 2 n_max.

    Con ellos se arma cualquier matriz de Toeplitz [[f]] de tamaño N × N,
    N = 2 n_max + 1 (las diferencias n_i − n_j llegan a ±2 n_max).
    """

    n_max: int
    period: float
    orders: NDArray[np.int64]
    indicator_a: NDArray[np.complex128]
    indicator_b: NDArray[np.complex128]

    @property
    def size(self) -> int:
        """Número de ondas planas N = 2 n_max + 1."""
        return 2 * self.n_max + 1

    @property
    def reciprocal_vectors(self) -> NDArray[np.float64]:
        """G_n = 2πn/L para n = −n_max … n_max (la base de ondas planas), en 1/[L]."""
        return np.asarray(2.0 * np.pi * np.arange(-self.n_max, self.n_max + 1) / self.period, dtype=np.float64)

    def profile_coefficients(self, value_a: complex, value_b: complex) -> NDArray[np.complex128]:
        """f_n del perfil que vale ``value_a`` en A y ``value_b`` en B."""
        return value_a * self.indicator_a + value_b * self.indicator_b

    def toeplitz(self, value_a: complex, value_b: complex) -> NDArray[np.complex128]:
        """Matriz de Toeplitz [[f]]_{ij} = f_{n_i − n_j}, de tamaño N × N."""
        return toeplitz_matrix(self.profile_coefficients(value_a, value_b), self.n_max)


def layer_edges(m: int, thickness_a: float, thickness_b: float) -> tuple[str, NDArray[np.float64]]:
    """Palabra S_m y posiciones de las interfaces z_0 = 0 < z_1 < … < z_N = L_m."""
    word = sequence(m)
    widths = np.array([thickness_a if s == LAYER_A else thickness_b for s in word])
    return word, np.concatenate(([0.0], np.cumsum(widths)))


def step_fourier_coefficients(
    values: ArrayLike,
    edges: NDArray[np.float64],
    orders: NDArray[np.int64],
) -> NDArray[np.complex128]:
    """Coeficientes f_n de un perfil constante por tramos (un valor por capa)."""
    period = float(edges[-1])
    values = np.asarray(values, dtype=np.complex128)
    g = 2.0 * np.pi * np.asarray(orders, dtype=np.float64) / period
    nonzero = orders != 0
    g_safe = np.where(nonzero, g, 1.0)
    phase_hi = np.exp(-1j * np.outer(g, edges[1:]))
    phase_lo = np.exp(-1j * np.outer(g, edges[:-1]))
    integrals = (1j / (g_safe[:, None] * period)) * (phase_hi - phase_lo)
    widths = np.diff(edges) / period
    integrals[~nonzero, :] = widths
    return np.asarray(integrals @ values, dtype=np.complex128)


def cell_fourier(m: int, thickness_a: float, thickness_b: float, n_max: int) -> CellFourier:
    """Coeficientes analíticos de las indicadoras de la celda S_m."""
    n_max = _validate_n_max(n_max)
    word, edges = layer_edges(m, thickness_a, thickness_b)
    orders = np.arange(-2 * n_max, 2 * n_max + 1)
    is_a = np.array([s == LAYER_A for s in word], dtype=float)
    return CellFourier(
        n_max=n_max,
        period=float(edges[-1]),
        orders=orders,
        indicator_a=step_fourier_coefficients(is_a, edges, orders),
        indicator_b=step_fourier_coefficients(1.0 - is_a, edges, orders),
    )


def bilayer_fourier(thickness_1: float, thickness_2: float, n_max: int) -> CellFourier:
    """Celda bicapa del libro (capa 1 en [0, l1], capa 2 en [l1, l1 + l2]); ec. (4.38)."""
    n_max = _validate_n_max(n_max)
    edges = np.array([0.0, thickness_1, thickness_1 + thickness_2])
    orders = np.arange(-2 * n_max, 2 * n_max + 1)
    return CellFourier(
        n_max=n_max,
        period=float(edges[-1]),
        orders=orders,
        indicator_a=step_fourier_coefficients(np.array([1.0, 0.0]), edges, orders),
        indicator_b=step_fourier_coefficients(np.array([0.0, 1.0]), edges, orders),
    )


def toeplitz_matrix(coefficients: NDArray[np.complex128], n_max: int) -> NDArray[np.complex128]:
    """T_{ij} = c_{n_i − n_j}, con ``coefficients`` indexado por n = −2 n_max … 2 n_max."""
    index = np.arange(-n_max, n_max + 1)
    return coefficients[index[:, None] - index[None, :] + 2 * n_max]


def synthesize_profile(
    coefficients: NDArray[np.complex128],
    orders: NDArray[np.int64],
    period: float,
    z: NDArray[np.float64],
) -> NDArray[np.complex128]:
    """Serie truncada f(z) = Σ_n f_n e^{i G_n z} (validación de la fig. 4.6 del libro)."""
    g = 2.0 * np.pi * np.asarray(orders, dtype=np.float64) / period
    return np.asarray(np.exp(1j * np.outer(z, g)) @ coefficients, dtype=np.complex128)
