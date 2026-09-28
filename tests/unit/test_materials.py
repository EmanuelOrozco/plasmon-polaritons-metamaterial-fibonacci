"""Medio A homogéneo y metamaterial B de Drude."""

import numpy as np
import pytest

from fibonacci_photonics.errors import InvalidParameterError
from fibonacci_photonics.physics.constants import GHZ
from fibonacci_photonics.physics.materials import DrudeMetamaterial, HomogeneousMedium
from fibonacci_photonics.physics.superlattice import SuperlatticeSpec


def test_air_is_vacuum_like(paper_equal_plasma):
    eps, mu = paper_equal_plasma.medium_a.epsilon_mu(2 * np.pi * 1.0e9)
    assert np.allclose(eps, 1.0)
    assert np.allclose(mu, 1.0)


def test_drude_vanishes_at_plasma_with_unit_background(paper_equal_plasma):
    omega_e = paper_equal_plasma.omega_e
    assert np.allclose(paper_equal_plasma.medium_b.epsilon(omega_e), 0.0, atol=1e-12)
    assert np.allclose(paper_equal_plasma.medium_b.mu(omega_e), 0.0, atol=1e-12)


def test_plasmon_frequencies_ghz(paper_equal_plasma, paper_split_plasma):
    assert paper_equal_plasma.nu_e_ghz() == pytest.approx(3.0)
    assert paper_equal_plasma.nu_m_ghz() == pytest.approx(3.0)
    assert paper_split_plasma.nu_e_ghz() == pytest.approx(3.0)
    assert paper_split_plasma.nu_m_ghz() == pytest.approx(1.0)


def test_negative_index_below_plasma(paper_equal_plasma):
    eps, mu = paper_equal_plasma.medium_b.epsilon_mu(2 * np.pi * 1.5 * GHZ)
    assert np.real(eps) < 0
    assert np.real(mu) < 0


def test_drude_diverges_at_zero_frequency(paper_equal_plasma):
    with pytest.raises(ZeroDivisionError):
        paper_equal_plasma.medium_b.epsilon(0.0)


def test_nondispersive_drude_is_finite_at_zero_frequency():
    eps, mu = DrudeMetamaterial(omega_e=0.0, omega_m=0.0, epsilon_0=9.0).epsilon_mu(np.array([0.0, 1.0]))
    np.testing.assert_array_equal(eps, [9.0, 9.0])
    np.testing.assert_array_equal(mu, [1.0, 1.0])


def test_homogeneous_broadcast():
    eps, mu = HomogeneousMedium(1.0, 1.0).epsilon_mu(np.ones((2, 3)))
    assert eps.shape == mu.shape == (2, 3)


@pytest.mark.parametrize("kwargs", [{"epsilon": 0.0, "mu": 1.0}, {"epsilon": 1.0, "mu": np.inf}])
def test_homogeneous_rejects_degenerate_values(kwargs):
    with pytest.raises(InvalidParameterError):
        HomogeneousMedium(**kwargs)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"omega_e": -1.0, "omega_m": 1.0},
        {"omega_e": 1.0, "omega_m": np.nan},
        {"omega_e": 1.0, "omega_m": 1.0, "epsilon_0": 0.0},
        {"omega_e": 1.0, "omega_m": 1.0, "mu_0": -2.0},
    ],
)
def test_drude_rejects_invalid_parameters(kwargs):
    with pytest.raises(InvalidParameterError):
        DrudeMetamaterial(**kwargs)


@pytest.mark.parametrize("field", ["thickness_a", "thickness_b", "speed_of_light"])
@pytest.mark.parametrize("value", [0.0, -1e-3, np.inf])
def test_spec_rejects_non_positive_dimensions(paper_split_plasma, field, value):
    kwargs = {
        "thickness_a": paper_split_plasma.thickness_a,
        "thickness_b": paper_split_plasma.thickness_b,
        "medium_a": paper_split_plasma.medium_a,
        "medium_b": paper_split_plasma.medium_b,
        field: value,
    }
    with pytest.raises(InvalidParameterError):
        SuperlatticeSpec(**kwargs)
