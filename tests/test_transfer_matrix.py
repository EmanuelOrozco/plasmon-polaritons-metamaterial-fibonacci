"""Tests 4–5: matriz de capa, producto, semitraza analítica y recurrencia."""

import numpy as np
import pytest

from fibonacci_tmm.electromagnetics import Polarization
from fibonacci_tmm.params import ghz_to_omega
from fibonacci_tmm.transfer_matrix import (
    analytic_r0_r1_r2,
    cell_matrix_by_recurrence_blocks,
    cell_transfer_matrix,
    determinant,
    semitrace,
    semitrace_by_product,
    semitrace_by_recurrence,
)


def _sample_omegas():
    return ghz_to_omega(np.array([0.5, 1.2, 2.1, 2.7, 3.4, 4.5]))


@pytest.mark.parametrize("polarization", [Polarization.TE, Polarization.TM])
@pytest.mark.parametrize("m", [0, 1, 2, 3, 4, 5, 6])
def test_unimodular_transfer_matrix(paper_equal_plasma, polarization, m):
    omega = _sample_omegas()
    t = cell_transfer_matrix(paper_equal_plasma, m, omega, np.pi / 6, polarization)
    det = determinant(t)
    assert np.allclose(det, 1.0, atol=1e-8)


@pytest.mark.parametrize("polarization", [Polarization.TE, Polarization.TM])
def test_analytic_r0_r1_r2_match_product(paper_equal_plasma, polarization):
    omega = _sample_omegas()
    theta = np.pi / 3
    r0, r1, r2 = analytic_r0_r1_r2(paper_equal_plasma, omega, theta, polarization)
    t0 = cell_transfer_matrix(paper_equal_plasma, 0, omega, theta, polarization)
    t1 = cell_transfer_matrix(paper_equal_plasma, 1, omega, theta, polarization)
    t2 = cell_transfer_matrix(paper_equal_plasma, 2, omega, theta, polarization)
    assert np.allclose(r0, semitrace(t0), atol=1e-10)
    assert np.allclose(r1, semitrace(t1), atol=1e-10)
    assert np.allclose(r2, semitrace(t2), atol=1e-10)


@pytest.mark.parametrize("polarization", [Polarization.TE, Polarization.TM])
@pytest.mark.parametrize("m", [3, 4, 5, 6, 7])
def test_recurrence_matches_direct_product(paper_equal_plasma, polarization, m):
    omega = _sample_omegas()
    theta = np.pi / 12
    r_rec = semitrace_by_recurrence(paper_equal_plasma, m, omega, theta, polarization)
    r_prod = semitrace_by_product(paper_equal_plasma, m, omega, theta, polarization)
    assert np.allclose(r_rec, r_prod, atol=1e-8, rtol=1e-8)


@pytest.mark.parametrize("m", [2, 3, 4, 5])
def test_block_recurrence_matches_layer_product(paper_split_plasma, m):
    omega = _sample_omegas()
    theta = np.pi / 6
    t_layers = cell_transfer_matrix(paper_split_plasma, m, omega, theta, Polarization.TE)
    t_blocks = cell_matrix_by_recurrence_blocks(
        paper_split_plasma, m, omega, theta, Polarization.TE
    )
    assert np.allclose(t_layers, t_blocks, atol=1e-8)


def test_normal_incidence_te_equals_tm(paper_equal_plasma):
    omega = _sample_omegas()
    r_te = semitrace_by_recurrence(paper_equal_plasma, 4, omega, 0.0, Polarization.TE)
    r_tm = semitrace_by_recurrence(paper_equal_plasma, 4, omega, 0.0, Polarization.TM)
    assert np.allclose(r_te, r_tm, atol=1e-10)
