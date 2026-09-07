"""Drude y medio A."""

import numpy as np
import pytest

from fibonacci_tmm.constants import GHZ
from fibonacci_tmm.materials import HomogeneousMedium


def test_air_is_vacuum_like(paper_equal_plasma):
    eps, mu = paper_equal_plasma.medium_a.epsilon_mu(2 * np.pi * 1.0e9)
    assert np.allclose(eps, 1.0)
    assert np.allclose(mu, 1.0)


def test_drude_vanishes_at_plasma_with_unit_background(paper_equal_plasma):
    omega_e = paper_equal_plasma.omega_e
    eps = paper_equal_plasma.medium_b.epsilon(omega_e)
    mu = paper_equal_plasma.medium_b.mu(omega_e)
    assert np.allclose(eps, 0.0, atol=1e-12)
    assert np.allclose(mu, 0.0, atol=1e-12)


def test_plasmon_frequencies_ghz(paper_equal_plasma, paper_split_plasma):
    assert paper_equal_plasma.nu_e_ghz() == pytest.approx(3.0)
    assert paper_equal_plasma.nu_m_ghz() == pytest.approx(3.0)
    assert paper_split_plasma.nu_e_ghz() == pytest.approx(3.0)
    assert paper_split_plasma.nu_m_ghz() == pytest.approx(1.0)


def test_negative_index_below_plasma(paper_equal_plasma):
    omega = 2 * np.pi * 1.5 * GHZ
    eps, mu = paper_equal_plasma.medium_b.epsilon_mu(omega)
    assert np.real(eps) < 0
    assert np.real(mu) < 0


def test_drude_diverges_at_zero_frequency(paper_equal_plasma):
    with pytest.raises(ZeroDivisionError):
        paper_equal_plasma.medium_b.epsilon(0.0)


def test_homogeneous_broadcast():
    medium = HomogeneousMedium(1.0, 1.0)
    omega = np.array([1.0, 2.0, 3.0])
    eps, mu = medium.epsilon_mu(omega)
    assert eps.shape == (3,)
    assert mu.shape == (3,)
