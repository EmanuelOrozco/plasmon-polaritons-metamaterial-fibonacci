"""Coeficientes de Fourier analíticos de la celda escalonada."""

import numpy as np
import pytest

from fibonacci_photonics.core.fibonacci import n_layers_a, n_layers_b
from fibonacci_photonics.pwe.fourier import (
    bilayer_fourier,
    cell_fourier,
    layer_edges,
    synthesize_profile,
    toeplitz_matrix,
)


def test_book_bilayer_matches_eq_4_38():
    """Mismos χ(G−G') que el notebook del libro (ε1 = 1, ε2 = 9, l1 = 0.2, l2 = 0.8)."""
    l1, l2, eps1, eps2, n_max = 0.2, 0.8, 1.0, 9.0, 6
    a = l1 + l2
    cell = bilayer_fourier(l1, l2, n_max)
    chi = cell.toeplitz(1 / eps1, 1 / eps2)
    g = np.arange(-n_max, n_max + 1) * 2 * np.pi / a
    for i, gi in enumerate(g):
        for j, gj in enumerate(g):
            if i == j:
                ref = (l1 / eps1 + l2 / eps2) / a
            else:
                d = gi - gj
                ref = 1j / a / d * (
                    (np.exp(-1j * d * l1) - 1) / eps1
                    + (np.exp(-1j * d * a) - np.exp(-1j * d * l1)) / eps2
                )
            assert chi[i, j] == pytest.approx(ref, abs=1e-14)


@pytest.mark.parametrize("m", [2, 3, 4, 5, 6])
def test_zero_order_is_filling_fraction(m):
    cell = cell_fourier(m, 1.0, 1.0, 4)
    zero = 2 * cell.n_max
    total = n_layers_a(m) + n_layers_b(m)
    assert cell.indicator_a[zero].real == pytest.approx(n_layers_a(m) / total)
    assert cell.indicator_b[zero].real == pytest.approx(n_layers_b(m) / total)


@pytest.mark.parametrize("m", [3, 5])
def test_indicators_are_complementary(m):
    cell = cell_fourier(m, 0.7, 1.3, 10)
    total = cell.indicator_a + cell.indicator_b
    expected = np.zeros_like(total)
    expected[2 * cell.n_max] = 1.0
    np.testing.assert_allclose(total, expected, atol=1e-14)


def test_toeplitz_is_hermitian_for_real_profile():
    cell = cell_fourier(4, 1.0, 2.0, 8)
    t = cell.toeplitz(1.0, -0.37)
    np.testing.assert_allclose(t, t.conj().T, atol=1e-14)
    assert t.shape == (17, 17)


def test_toeplitz_indexing():
    n_max = 2
    coeffs = np.arange(-2 * n_max, 2 * n_max + 1).astype(complex)
    t = toeplitz_matrix(coeffs, n_max)
    for i in range(2 * n_max + 1):
        for j in range(2 * n_max + 1):
            assert t[i, j] == (i - j)


def test_synthesis_converges_to_profile_away_from_interfaces():
    word, edges = layer_edges(4, 1.0, 1.0)
    cell = cell_fourier(4, 1.0, 1.0, 400)
    coeffs = cell.profile_coefficients(1.0, 9.0)
    centers = 0.5 * (edges[1:] + edges[:-1])
    values = synthesize_profile(coeffs, cell.orders, cell.period, centers)
    expected = np.array([1.0 if s == "A" else 9.0 for s in word])
    np.testing.assert_allclose(values.real, expected, atol=2e-2)
    np.testing.assert_allclose(values.imag, 0.0, atol=1e-10)
