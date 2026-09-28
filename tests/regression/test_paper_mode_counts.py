"""Resultado central del PRB: F_{m−2} subbandas plasmón-polaritón bajo ν_m."""

import numpy as np
import pytest

from fibonacci_photonics.analysis.plasmon_modes import detect_plasmon_modes, merge_touching_intervals
from fibonacci_photonics.physics.fibonacci import n_layers_b
from fibonacci_photonics.solvers import TMMSolver


@pytest.mark.parametrize("m, expected", [(3, 1), (4, 2), (5, 3), (6, 5)])
def test_mode_count_figure_2_window(paper_equal_plasma, m, expected):
    nu = np.linspace(2.0, 2.995, 8000)
    scan = TMMSolver(paper_equal_plasma).scan(m, np.pi / 3, nu)
    report = detect_plasmon_modes(scan, nu_min_ghz=2.0, nu_max_ghz=2.995, min_points=5)
    assert len(merge_touching_intervals(list(report.intervals), gap_ghz=2e-3)) == expected


@pytest.mark.parametrize("m, expected", [(3, 1), (4, 2)])
def test_mode_count_figure_4_window(paper_split_plasma, m, expected):
    nu = np.linspace(0.94, 0.9995, 8000)
    scan = TMMSolver(paper_split_plasma).scan(m, np.pi / 3, nu)
    report = detect_plasmon_modes(scan, nu_min_ghz=0.94, nu_max_ghz=0.9995, min_points=5)
    assert len(merge_touching_intervals(list(report.intervals), gap_ghz=5e-5)) == expected


@pytest.mark.parametrize("m", range(2, 9))
@pytest.mark.parametrize("theta, window", [(np.pi / 12, (0.994, 1.0)), (np.pi / 3, (0.95, 1.0))])
def test_mode_count_figure_5_windows(paper_split_plasma, m, theta, window):
    nu = np.linspace(window[0], window[1] - 1e-4, 20000)
    scan = TMMSolver(paper_split_plasma).scan(m, theta, nu)
    report = detect_plasmon_modes(scan, nu_min_ghz=window[0], nu_max_ghz=window[1])
    assert report.detected_modes == n_layers_b(m)


@pytest.mark.parametrize("m", [6, 7])
@pytest.mark.parametrize("theta_deg", [3.0, 13.5, 45.0])
def test_mode_count_figure_6_grid(paper_split_plasma, m, theta_deg):
    nu_m = paper_split_plasma.nu_m_ghz()
    nu = np.unique(
        np.concatenate(
            [
                nu_m - np.logspace(-10, np.log10(nu_m - 0.90), 60000),
                np.linspace(0.95, nu_m - 1e-10, 250000),
            ]
        )
    )
    scan = TMMSolver(paper_split_plasma).scan(m, np.radians(theta_deg), nu)
    report = detect_plasmon_modes(scan, nu_min_ghz=0.90, nu_max_ghz=1.02, min_points=2)
    assert report.detected_modes == n_layers_b(m)
