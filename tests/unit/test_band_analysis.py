"""Raíces, bordes de banda, subbandas plasmónicas y comparación de barridos."""

import numpy as np
import pytest

from fibonacci_photonics.analysis.band_edges import (
    locate_bands,
    plasmon_grid,
    refine_band_edges,
    refine_grid_interval,
    validate_bands,
)
from fibonacci_photonics.analysis.plasmon_modes import merge_touching_intervals
from fibonacci_photonics.analysis.roots import refine_crossing
from fibonacci_photonics.benchmark.compare import compare_scans, match_intervals
from fibonacci_photonics.errors import InvalidParameterError
from fibonacci_photonics.solvers import FrequencyInterval, TMMSolver, allowed_intervals
from fibonacci_photonics.solvers.scan import allowed_mask


def test_allowed_mask_counts_near_unit_clipping():
    r = np.array([0.0, 0.999, 1.0, 1.0 + 0.5e-12, 1.1], dtype=np.complex128)
    allowed, n_clipped = allowed_mask(r, atol=1e-12)
    assert allowed.tolist() == [True, True, True, True, False]
    assert n_clipped == 1


def test_identical_scans_agree(paper_equal_plasma):
    scan = TMMSolver(paper_equal_plasma).scan(3, np.pi / 6, np.linspace(0.5, 5.0, 300))
    agreement = compare_scans(scan, scan)
    assert agreement.classification_agreement == 1.0
    assert agreement.max_abs_delta_r_in_bands == 0.0


def test_different_grids_are_rejected(paper_equal_plasma):
    solver = TMMSolver(paper_equal_plasma)
    with pytest.raises(ValueError):
        compare_scans(solver.scan(3, 0.1, np.linspace(1, 2, 10)), solver.scan(3, 0.1, np.linspace(1, 2, 11)))


def test_match_intervals_by_overlap():
    ref = [FrequencyInterval(1.0, 2.0, 5), FrequencyInterval(3.0, 4.0, 5)]
    cand = [FrequencyInterval(3.1, 4.1, 5), FrequencyInterval(0.9, 1.9, 5)]
    pairs = match_intervals(ref, cand)
    assert pairs[0][1] is cand[1]
    assert pairs[1][1] is cand[0]


def test_refine_crossing():
    assert refine_crossing(lambda x: x**2 - 2.0, 1.0, 2.0) == pytest.approx(np.sqrt(2.0), abs=1e-12)
    assert refine_crossing(lambda x: x**2 + 1.0, -1.0, 1.0) is None
    assert refine_crossing(lambda x: np.nan, -1.0, 1.0) is None


def test_locate_bands_finds_band_between_two_samples():
    """Una banda más estrecha que el paso de malla se detecta por el cambio de signo de R."""

    def r(nu):
        return 50.0 * (nu - 0.5)

    grid = np.array([0.0, 0.3, 0.7, 1.0])
    bands = locate_bands(r, grid, np.array([r(x) for x in grid]))
    assert len(bands) == 1
    assert bands[0].nu_min_ghz == pytest.approx(0.48, abs=1e-12)
    assert bands[0].nu_max_ghz == pytest.approx(0.52, abs=1e-12)


@pytest.mark.parametrize("grid", [[0.0, 0.47, 0.53, 1.0], [0.0, 0.45, 0.52, 1.0]])
def test_locate_bands_finds_gap_between_two_samples(grid):
    """Un gap más estrecho que la malla, con muestras en banda a ambos lados, separa dos bandas."""

    def r(nu):
        return 1.0 + 1e-6 - 400.0 * (nu - 0.5) ** 2

    grid = np.array(grid)
    bands = locate_bands(r, grid, np.array([r(x) for x in grid]))
    half_gap = np.sqrt(1e-6 / 400.0)
    assert len(bands) == 2
    assert bands[0].nu_max_ghz == pytest.approx(0.5 - half_gap, abs=1e-10)
    assert bands[1].nu_min_ghz == pytest.approx(0.5 + half_gap, abs=1e-10)


def test_validate_bands_rejects_unconverged():
    good, bad = FrequencyInterval(0.1, 0.2, 0), FrequencyInterval(0.8, 0.9, 0)
    kept, rejected = validate_bands([good, bad], lambda nu: 0.3, lambda nu: 0.3 if nu < 0.5 else 1e6)
    assert kept == [good] and rejected == [bad]


def test_plasmon_grid_accumulates_at_nu_m():
    grid = plasmon_grid(0.9, 1.0, 50, 1e-8)
    assert grid[0] == pytest.approx(0.9) and grid[-1] == pytest.approx(1.0 - 1e-8)
    assert np.all(np.diff(grid) > 0)
    with pytest.raises(InvalidParameterError):
        plasmon_grid(0.9, 1.0, 50, 0.5)


def test_refined_edges_are_where_abs_r_is_one(paper_equal_plasma):
    solver = TMMSolver(paper_equal_plasma)
    nu = np.linspace(0.5, 5.0, 2000)
    interval = allowed_intervals(solver.scan(3, np.pi / 6, nu))[1]
    r_of_nu = solver.evaluator(3, np.pi / 6)
    lo, hi = refine_band_edges(r_of_nu, interval, nu[1] - nu[0])
    assert abs(abs(r_of_nu(lo)) - 1.0) < 1e-9
    assert abs(abs(r_of_nu(hi)) - 1.0) < 1e-9
    assert lo <= interval.nu_min_ghz and hi >= interval.nu_max_ghz


def test_refine_grid_interval_on_nonuniform_grid():
    """Los bordes pasan de las muestras extremas al cruce |R| = 1 aunque la malla no sea uniforme."""

    def r(nu):
        return 10.0 * (nu - 0.5)

    grid = np.array([0.0, 0.2, 0.41, 0.45, 0.5, 0.58, 0.8, 1.0])
    inside = FrequencyInterval(0.41, 0.58, 4)
    refined = refine_grid_interval(r, grid, inside)
    assert refined.nu_min_ghz == pytest.approx(0.4, abs=1e-12)
    assert refined.nu_max_ghz == pytest.approx(0.6, abs=1e-12)
    at_edge = refine_grid_interval(r, grid[2:6], inside)
    assert (at_edge.nu_min_ghz, at_edge.nu_max_ghz) == (0.41, 0.58)


def test_merge_touching_intervals():
    parts = [FrequencyInterval(1.0, 1.1, 3), FrequencyInterval(1.1005, 1.2, 3), FrequencyInterval(1.5, 1.6, 3)]
    merged = merge_touching_intervals(parts, gap_ghz=1e-3)
    assert [(i.nu_min_ghz, i.nu_max_ghz) for i in merged] == [(1.0, 1.2), (1.5, 1.6)]
