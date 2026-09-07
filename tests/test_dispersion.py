"""Tests 6–7: condición de banda e incidencia normal."""

import numpy as np
import pytest

from fibonacci_tmm.dispersion import (
    allowed_mask,
    scan_dispersion,
    zero_average_index_frequency_ghz,
)
from fibonacci_tmm.electromagnetics import Polarization


def test_zero_average_index_frequency_formula(paper_equal_plasma):
    nu3 = zero_average_index_frequency_ghz(paper_equal_plasma, 3)
    nu4 = zero_average_index_frequency_ghz(paper_equal_plasma, 4)
    assert nu3 == pytest.approx(3.0 * np.sqrt(1 / 3), rel=1e-12)
    assert nu4 == pytest.approx(3.0 * np.sqrt(2 / 5), rel=1e-12)


def test_allowed_mask_counts_near_unit_clipping():
    r = np.array([0.0, 0.999, 1.0, 1.0 + 0.5e-12, 1.1], dtype=np.complex128)
    allowed, n_clipped = allowed_mask(r, atol=1e-12)
    assert allowed.tolist() == [True, True, True, True, False]
    assert n_clipped == 1


def test_gaps_exist_at_oblique_incidence(paper_equal_plasma):
    nu = np.linspace(0.05, 5.0, 2000)
    scan = scan_dispersion(paper_equal_plasma, 3, np.pi / 3, nu, Polarization.TE)
    assert np.any(scan.allowed)
    assert np.any(~scan.allowed)
    assert np.all((scan.k_lm_over_pi[scan.allowed] >= 0.0) & (scan.k_lm_over_pi[scan.allowed] <= 1.0))
    assert np.all(np.isnan(scan.k_lm_over_pi[~scan.allowed]))


def test_normal_incidence_has_allowed_states_near_closed_n0_gap(paper_equal_plasma):
    """A θ = 0 el gap ⟨n⟩=0 está cerrado: debe haber banda en ν = 3/√3 GHz."""
    nu_gap = zero_average_index_frequency_ghz(paper_equal_plasma, 3)
    nu = np.linspace(nu_gap - 0.05, nu_gap + 0.05, 400)
    scan0 = scan_dispersion(paper_equal_plasma, 3, 0.0, nu, Polarization.TE)
    scan_oblique = scan_dispersion(paper_equal_plasma, 3, np.pi / 3, nu, Polarization.TE)
    assert np.any(scan0.allowed)
    # A incidencia oblicua el gap no Bragg se abre: menos puntos permitidos cerca del cierre.
    assert scan_oblique.allowed.mean() <= scan0.allowed.mean() + 1e-12
