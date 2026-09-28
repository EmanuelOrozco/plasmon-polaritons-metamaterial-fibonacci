"""Matriz de capa, producto de celda, semitraza analítica y recurrencia."""

import numpy as np
import pytest

from fibonacci_photonics.benchmark.compare import compare_scans
from fibonacci_photonics.physics.electromagnetics import Polarization
from fibonacci_photonics.physics.units import omega_from_nu_ghz
from fibonacci_photonics.solvers import TMMSolver
from fibonacci_photonics.solvers.tmm.transfer_matrix import (
    analytic_r0_r1_r2,
    cell_matrix_by_recurrence_blocks,
    cell_transfer_matrix,
    semitrace,
    semitrace_by_product,
    semitrace_by_recurrence,
)

OMEGAS = omega_from_nu_ghz(np.array([0.5, 1.2, 2.1, 2.7, 3.4, 4.5]))
POLARIZATIONS = [Polarization.TE, Polarization.TM]


@pytest.mark.parametrize("polarization", POLARIZATIONS)
def test_analytic_r0_r1_r2_match_product(paper_equal_plasma, polarization):
    theta = np.pi / 3
    r0, r1, r2 = analytic_r0_r1_r2(paper_equal_plasma, OMEGAS, theta, polarization)
    for m, r in enumerate((r0, r1, r2)):
        t = cell_transfer_matrix(paper_equal_plasma, m, OMEGAS, theta, polarization)
        np.testing.assert_allclose(r, semitrace(t), atol=1e-10)


@pytest.mark.parametrize("polarization", POLARIZATIONS)
@pytest.mark.parametrize("m", [3, 4, 5, 6, 7])
def test_recurrence_matches_direct_product(paper_equal_plasma, polarization, m):
    r_rec = semitrace_by_recurrence(paper_equal_plasma, m, OMEGAS, np.pi / 12, polarization)
    r_prod = semitrace_by_product(paper_equal_plasma, m, OMEGAS, np.pi / 12, polarization)
    np.testing.assert_allclose(r_rec, r_prod, atol=1e-8, rtol=1e-8)


@pytest.mark.parametrize("m", [2, 3, 4, 5])
def test_block_recurrence_matches_layer_product(paper_split_plasma, m):
    t_layers = cell_transfer_matrix(paper_split_plasma, m, OMEGAS, np.pi / 6, Polarization.TE)
    t_blocks = cell_matrix_by_recurrence_blocks(paper_split_plasma, m, OMEGAS, np.pi / 6, Polarization.TE)
    np.testing.assert_allclose(t_layers, t_blocks, atol=1e-8)


def test_solver_algorithms_agree_on_a_scan(paper_equal_plasma):
    nu = np.linspace(0.5, 5.0, 300)
    recurrence = TMMSolver(paper_equal_plasma).scan(4, np.pi / 6, nu)
    product = TMMSolver(paper_equal_plasma, algorithm="product").scan(4, np.pi / 6, nu)
    assert recurrence.method == "recurrence" and product.method == "product"
    assert compare_scans(recurrence, product).max_abs_delta_r_in_bands < 1e-10


def test_gaps_exist_at_oblique_incidence(paper_equal_plasma):
    scan = TMMSolver(paper_equal_plasma).scan(3, np.pi / 3, np.linspace(0.05, 5.0, 2000))
    assert np.any(scan.allowed) and np.any(~scan.allowed)
