"""PWE frente a la TMM, cierre de bandas y forma ω(k) del libro."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from fibonacci_photonics.analysis.band_closure import (
    EDGE,
    EXTREMUM,
    close_band_edges,
    edge_brackets,
)
from fibonacci_photonics.config import spec_from_mapping
from fibonacci_photonics.physics.electromagnetics import Polarization
from fibonacci_photonics.physics.units import omega_from_nu_ghz
from fibonacci_photonics.solvers import PWESolver, TMMSolver
from fibonacci_photonics.solvers.pwe.bloch import n_max_for, select_bloch_wavevector
from fibonacci_photonics.solvers.pwe.eigenfrequency import (
    drude_qep_frequencies_ghz,
    mode_profile,
    nondispersive_bands,
)
from fibonacci_photonics.solvers.pwe.fourier import bilayer_fourier
from fibonacci_photonics.solvers.tmm.transfer_matrix import semitrace_by_recurrence
from fibonacci_photonics.viz.dispersion import plot_dispersion_branches

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


@pytest.mark.parametrize("polarization", ["TE", "TM"])
def test_drude_semitrace_matches_tmm(paper_split_plasma, polarization):
    nu = np.linspace(0.3, 5.0, 24)
    tmm = TMMSolver(paper_split_plasma).scan(3, np.pi / 6, nu, polarization)
    pwe = PWESolver(paper_split_plasma).scan(3, np.pi / 6, nu, polarization)
    assert np.array_equal(tmm.allowed, pwe.allowed)
    assert np.max(np.abs(tmm.r.real - pwe.r.real)[tmm.allowed]) < 1e-5
    assert pwe.method == "pwe-inverse"


def test_plasmon_window_semitrace_matches_tmm(paper_split_plasma):
    nu = np.linspace(0.9955, 0.9994, 10)
    r_tmm = TMMSolver(paper_split_plasma).semitrace(3, np.pi / 12, nu).real
    r_pwe = PWESolver(paper_split_plasma).semitrace(3, np.pi / 12, nu).real
    assert np.max(np.abs(r_tmm - r_pwe)) < 1e-5


def test_laurent_rule_is_worse_for_drude(paper_split_plasma):
    nu = np.linspace(0.3, 5.0, 24)
    r_tmm = TMMSolver(paper_split_plasma).semitrace(3, np.pi / 6, nu, "TM").real
    errors = {}
    for rule in ("inverse", "laurent"):
        r = PWESolver(paper_split_plasma, harmonics_per_layer=8, rule=rule).semitrace(3, np.pi / 6, nu, "TM").real
        errors[rule] = np.max(np.abs(r - r_tmm)[np.abs(r_tmm) <= 1])
    assert errors["inverse"] < errors["laurent"] / 100


def test_parallel_matches_serial(paper_split_plasma):
    nu = np.linspace(1.5, 4.5, 16)
    serial = PWESolver(paper_split_plasma, harmonics_per_layer=4).semitrace(3, np.pi / 3, nu)
    parallel = PWESolver(paper_split_plasma, harmonics_per_layer=4, workers=2).semitrace(3, np.pi / 3, nu)
    np.testing.assert_allclose(parallel, serial, rtol=0, atol=1e-12)


@pytest.fixture
def coarse_solver(paper_split_plasma):
    return PWESolver(paper_split_plasma, harmonics_per_layer=4)


@pytest.fixture
def coarse_scan(coarse_solver):
    return coarse_solver.scan(3, np.pi / 3, np.linspace(0.5, 5.0, 60))


def test_band_edges_lie_on_unit_semitrace(coarse_solver, coarse_scan):
    edges = close_band_edges(coarse_solver, coarse_scan)
    brackets = edge_brackets(coarse_scan)
    assert edges.n_edges >= len(brackets) > 0
    is_edge = edges.kind == EDGE
    evaluator = coarse_solver.evaluator(3, np.pi / 3)
    np.testing.assert_allclose(evaluator.many(edges.nu_ghz[is_edge]), edges.r[is_edge], rtol=0, atol=1e-5)
    np.testing.assert_array_equal(edges.k_lm_over_pi[is_edge], np.where(edges.r[is_edge] > 0, 0.0, 1.0))
    for inside, outside, _ in brackets:
        lo, hi = sorted((inside, outside))
        assert np.any(is_edge & (edges.nu_ghz >= lo) & (edges.nu_ghz <= hi))


def test_extrema_close_touching_bands(paper_equal_plasma):
    """Con impedancias iguales y θ = 0 no hay gaps: R = cos φ(ν) toca ±1 entre muestras."""
    solver = PWESolver(paper_equal_plasma, harmonics_per_layer=4)
    scan = solver.scan(3, 0.0, np.linspace(0.5, 5.0, 80))
    edges = close_band_edges(solver, scan)
    assert edges.n_extrema > 0
    k_tips = edges.k_lm_over_pi[edges.kind == EXTREMUM]
    assert np.all(np.minimum(k_tips, 1.0 - k_tips) < 1e-3)
    sampled = np.where(scan.allowed, scan.k_lm_over_pi, np.nan)
    assert np.nanmin(np.minimum(sampled, 1.0 - sampled)) > np.max(np.minimum(k_tips, 1.0 - k_tips))


def test_band_edges_parallel_matches_serial(coarse_solver, coarse_scan):
    serial = close_band_edges(coarse_solver, coarse_scan)
    parallel = close_band_edges(coarse_solver, coarse_scan, workers=2)
    np.testing.assert_allclose(parallel.nu_ghz, serial.nu_ghz, rtol=0, atol=1e-12)
    np.testing.assert_array_equal(parallel.r, serial.r)


def test_band_closure_works_with_the_tmm(paper_split_plasma):
    solver = TMMSolver(paper_split_plasma)
    edges = close_band_edges(solver, solver.scan(3, np.pi / 3, np.linspace(0.5, 5.0, 60)))
    is_edge = edges.kind == EDGE
    r = solver.semitrace(3, np.pi / 3, edges.nu_ghz[is_edge]).real
    np.testing.assert_allclose(np.abs(r), 1.0, atol=1e-9)


def test_branches_close_at_band_edges(coarse_solver, coarse_scan):
    edges = close_band_edges(coarse_solver, coarse_scan)
    fig, ax = plt.subplots()
    plot_dispersion_branches(
        ax, coarse_scan, color="k", linestyle="-", edge_nu_ghz=edges.nu_ghz, edge_k_lm_over_pi=edges.k_lm_over_pi
    )
    k_line, nu_line = ax.lines[1].get_xdata(), ax.lines[1].get_ydata()
    for nu_edge, k_edge in zip(edges.nu_ghz, edges.k_lm_over_pi, strict=True):
        if not np.isfinite(k_edge):
            continue
        i = int(np.flatnonzero(nu_line == nu_edge)[0])
        assert k_line[i] == k_edge
        assert np.count_nonzero(np.isfinite(k_line[max(i - 1, 0) : i + 2])) >= 2
    plt.close(fig)


def test_selection_prefers_propagating_mode():
    eigenvalues = np.array([0.4 + 1e-12j, 5.0 + 0.0j, 0.1 + 0.3j, -9.0 + 0.0j])
    assert select_bloch_wavevector(eigenvalues) == pytest.approx(0.4)


def test_truncation_scales_with_layers():
    assert n_max_for(3, 16) == 16 * 3
    assert n_max_for(5, 4) == 4 * 8


def test_drude_omega_k_has_spectral_pollution(paper_split_plasma):
    """Con Laurent aparecen autovalores que se acumulan en ν_m = 1 GHz."""
    nu = drude_qep_frequencies_ghz(paper_split_plasma, 3, 0.5 * np.pi, np.pi / 12, 20)
    assert np.count_nonzero((nu > 0.99) & (nu < 1.0)) >= 5
