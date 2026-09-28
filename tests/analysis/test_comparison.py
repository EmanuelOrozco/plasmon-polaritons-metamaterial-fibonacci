"""Herramientas de comparación entre métodos."""

import numpy as np
import pytest

from fibonacci_photonics.analysis.band_edges import locate_bands, plasmon_grid, refine_grid_interval, validate_bands
from fibonacci_photonics.analysis.comparison import (
    compare_scans,
    match_intervals,
    refine_band_edges,
    refine_crossing,
)
from fibonacci_photonics.core.bands import FrequencyInterval, allowed_intervals, omega_from_nu_ghz
from fibonacci_photonics.core.electromagnetics import Polarization
from fibonacci_photonics.tmm.dispersion import scan_dispersion
from fibonacci_photonics.tmm.transfer_matrix import semitrace_by_recurrence


def test_identical_scans_agree(paper_equal_plasma):
    nu = np.linspace(0.5, 5.0, 300)
    scan = scan_dispersion(paper_equal_plasma, 3, np.pi / 6, nu)
    agreement = compare_scans(scan, scan)
    assert agreement.classification_agreement == 1.0
    assert agreement.max_abs_delta_r_in_bands == 0.0


def test_product_and_recurrence_agree(paper_equal_plasma):
    nu = np.linspace(0.5, 5.0, 300)
    a = scan_dispersion(paper_equal_plasma, 4, np.pi / 6, nu, method="recurrence")
    b = scan_dispersion(paper_equal_plasma, 4, np.pi / 6, nu, method="product")
    assert compare_scans(a, b).max_abs_delta_r_in_bands < 1e-10


def test_different_grids_are_rejected(paper_equal_plasma):
    a = scan_dispersion(paper_equal_plasma, 3, 0.1, np.linspace(1, 2, 10))
    b = scan_dispersion(paper_equal_plasma, 3, 0.1, np.linspace(1, 2, 11))
    with pytest.raises(ValueError):
        compare_scans(a, b)


def test_match_intervals_by_overlap():
    ref = [FrequencyInterval(1.0, 2.0, 5), FrequencyInterval(3.0, 4.0, 5)]
    cand = [FrequencyInterval(3.1, 4.1, 5), FrequencyInterval(0.9, 1.9, 5)]
    pairs = match_intervals(ref, cand)
    assert pairs[0][1] is cand[1]
    assert pairs[1][1] is cand[0]


def test_refine_crossing():
    assert refine_crossing(lambda x: x**2 - 2.0, 1.0, 2.0) == pytest.approx(np.sqrt(2.0), abs=1e-12)
    assert refine_crossing(lambda x: x**2 + 1.0, -1.0, 1.0) is None


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

    def r1(nu):
        return 0.3

    def r2(nu):
        return 0.3 if nu < 0.5 else 1e6

    kept, rejected = validate_bands([good, bad], r1, r2)
    assert kept == [good] and rejected == [bad]


def test_plasmon_grid_accumulates_at_nu_m():
    grid = plasmon_grid(0.9, 1.0, 50, 1e-8)
    assert grid[0] == pytest.approx(0.9) and grid[-1] == pytest.approx(1.0 - 1e-8)
    assert np.all(np.diff(grid) > 0)


def test_refined_edges_are_where_abs_r_is_one(paper_equal_plasma):
    nu = np.linspace(0.5, 5.0, 2000)
    scan = scan_dispersion(paper_equal_plasma, 3, np.pi / 6, nu)
    interval = allowed_intervals(scan)[1]

    def r_of_nu(value):
        omega = omega_from_nu_ghz(np.array([value]))
        return float(semitrace_by_recurrence(paper_equal_plasma, 3, omega, np.pi / 6, Polarization.TE).real[0])

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
