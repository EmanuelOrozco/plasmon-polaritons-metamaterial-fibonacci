"""PWE frente a la TMM: benchmark del libro y superred de Drude."""

import numpy as np
import pytest

from fibonacci_photonics.core.bands import omega_from_nu_ghz
from fibonacci_photonics.core.electromagnetics import Polarization
from fibonacci_photonics.core.params import spec_from_mapping
from fibonacci_photonics.pwe.bloch import n_max_for, select_bloch_wavevector
from fibonacci_photonics.pwe.dispersion import pwe_semitrace, scan_dispersion as pwe_scan
from fibonacci_photonics.pwe.eigenfrequency import (
    drude_qep_frequencies_ghz,
    mode_profile,
    nondispersive_bands,
)
from fibonacci_photonics.pwe.fourier import bilayer_fourier
from fibonacci_photonics.tmm.dispersion import scan_dispersion as tmm_scan
from fibonacci_photonics.tmm.transfer_matrix import semitrace_by_recurrence

BOOK_PERIOD = 1.0e-6


@pytest.fixture
def book_bilayer():
    """Cristal del capítulo 1D: ε1 = 1 (0.2 µm), ε2 = 9 (0.8 µm), μ = 1."""
    return spec_from_mapping(
        {
            "layer_a_thickness_mm": 0.2e-3,
            "layer_b_thickness_mm": 0.8e-3,
            "omega_e_over_2pi_ghz": 0.0,
            "omega_m_over_2pi_ghz": 0.0,
            "epsilon_0": 9.0,
        }
    )


def _tmm_residual(spec, normalized_freq, k_tilde):
    nu_ghz = normalized_freq * spec.speed_of_light / BOOK_PERIOD / 1e9
    r = semitrace_by_recurrence(spec, 2, omega_from_nu_ghz(nu_ghz), 0.0, Polarization.TM).real
    return r - np.cos(k_tilde)


@pytest.mark.parametrize("rule,tol", [("laurent", 1e-1), ("inverse", 1e-4)])
def test_book_band_structure_satisfies_tmm_dispersion(book_bilayer, rule, tol):
    ks = np.array([0.25, 0.5, 0.75]) * np.pi
    freqs = nondispersive_bands(bilayer_fourier(0.2, 0.8, 50), (1, 9), (1, 1), ks, rule=rule, n_bands=4)
    for k, row in zip(ks, freqs, strict=True):
        assert np.max(np.abs(_tmm_residual(book_bilayer, row, k))) < tol


def test_inverse_rule_converges_faster_than_laurent():
    k = np.array([0.5 * np.pi])
    ref = nondispersive_bands(bilayer_fourier(0.2, 0.8, 200), (1, 9), (1, 1), k, rule="inverse", n_bands=4)
    errors = {}
    for rule in ("laurent", "inverse"):
        f = nondispersive_bands(bilayer_fourier(0.2, 0.8, 20), (1, 9), (1, 1), k, rule=rule, n_bands=4)
        errors[rule] = np.max(np.abs(f - ref) / ref)
    assert errors["inverse"] < errors["laurent"] / 20


def test_mode_profile_is_bloch_periodic():
    cell = bilayer_fourier(0.2, 0.8, 30)
    k = np.pi / 15
    _, modes = nondispersive_bands(cell, (1, 9), (1, 1), [k], n_bands=3, return_modes=True)
    u0 = mode_profile(cell, modes[0][:, 1], k, np.array([0.1]))
    u1 = mode_profile(cell, modes[0][:, 1], k, np.array([1.1]))
    np.testing.assert_allclose(u1, u0 * np.exp(1j * k), rtol=1e-10)


@pytest.mark.parametrize("polarization", [Polarization.TE, Polarization.TM])
def test_drude_semitrace_matches_tmm(paper_split_plasma, polarization):
    nu = np.linspace(0.3, 5.0, 24)
    tmm = tmm_scan(paper_split_plasma, 3, np.pi / 6, nu, polarization)
    pwe = pwe_scan(paper_split_plasma, 3, np.pi / 6, nu, polarization, harmonics_per_layer=16)
    assert np.array_equal(tmm.allowed, pwe.allowed)
    both = tmm.allowed & pwe.allowed
    assert np.max(np.abs(tmm.r.real - pwe.r.real)[both]) < 1e-5
    assert pwe.method == "pwe-inverse"


def test_plasmon_window_semitrace_matches_tmm(paper_split_plasma):
    nu = np.linspace(0.9955, 0.9994, 10)
    r_tmm = tmm_scan(paper_split_plasma, 3, np.pi / 12, nu).r.real
    r_pwe = pwe_semitrace(paper_split_plasma, 3, np.pi / 12, nu, harmonics_per_layer=16)
    assert np.max(np.abs(r_tmm - r_pwe)) < 1e-5


def test_laurent_rule_is_worse_for_drude(paper_split_plasma):
    nu = np.linspace(0.3, 5.0, 24)
    r_tmm = tmm_scan(paper_split_plasma, 3, np.pi / 6, nu, Polarization.TM).r.real
    kwargs = dict(polarization=Polarization.TM, harmonics_per_layer=8)
    inverse = pwe_semitrace(paper_split_plasma, 3, np.pi / 6, nu, rule="inverse", **kwargs)
    laurent = pwe_semitrace(paper_split_plasma, 3, np.pi / 6, nu, rule="laurent", **kwargs)
    in_band = np.abs(r_tmm) <= 1
    assert np.max(np.abs(inverse - r_tmm)[in_band]) < np.max(np.abs(laurent - r_tmm)[in_band]) / 100


def test_parallel_matches_serial(paper_split_plasma):
    nu = np.linspace(1.5, 4.5, 16)
    serial = pwe_semitrace(paper_split_plasma, 3, np.pi / 3, nu, harmonics_per_layer=4)
    parallel = pwe_semitrace(paper_split_plasma, 3, np.pi / 3, nu, harmonics_per_layer=4, workers=2)
    np.testing.assert_allclose(parallel, serial, rtol=0, atol=1e-12)


def test_selection_prefers_propagating_mode():
    eigenvalues = np.array([0.4 + 1e-12j, 5.0 + 0.0j, 0.1 + 0.3j, -9.0 + 0.0j])
    assert select_bloch_wavevector(eigenvalues) == pytest.approx(0.4)


def test_truncation_scales_with_layers():
    assert n_max_for(3, 16) == 16 * 3
    assert n_max_for(5, 4) == 4 * 8


def test_drude_omega_k_has_spectral_pollution(paper_split_plasma):
    """Con Laurent aparecen autovalores que se acumulan en ν_m = 1 GHz."""
    nu = drude_qep_frequencies_ghz(paper_split_plasma, 3, 0.5 * np.pi, np.pi / 12, 20)
    near = nu[(nu > 0.99) & (nu < 1.0)]
    assert near.size >= 5
