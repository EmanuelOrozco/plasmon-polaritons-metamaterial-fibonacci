"""Test 8: número de modos plasmon-polaritón = F(m-2)."""

import numpy as np
import pytest

from fibonacci_tmm.dispersion import scan_dispersion
from fibonacci_tmm.electromagnetics import Polarization
from fibonacci_tmm.plasmon_modes import detect_plasmon_modes, merge_touching_intervals


@pytest.mark.parametrize("m, expected", [(3, 1), (4, 2), (5, 3), (6, 5)])
def test_mode_count_figure_2_window(paper_equal_plasma, m, expected):
    nu = np.linspace(2.0, 2.995, 8000)
    scan = scan_dispersion(paper_equal_plasma, m, np.pi / 3, nu, Polarization.TE)
    report = detect_plasmon_modes(scan, nu_min_ghz=2.0, nu_max_ghz=2.995, min_points=5)
    intervals = merge_touching_intervals(list(report.intervals), gap_ghz=2e-3)
    assert len(intervals) == expected


@pytest.mark.parametrize("m, expected", [(3, 1), (4, 2)])
def test_mode_count_figure_4_window(paper_split_plasma, m, expected):
    nu = np.linspace(0.94, 0.9995, 8000)
    scan = scan_dispersion(paper_split_plasma, m, np.pi / 3, nu, Polarization.TE)
    report = detect_plasmon_modes(scan, nu_min_ghz=0.94, nu_max_ghz=0.9995, min_points=5)
    intervals = merge_touching_intervals(list(report.intervals), gap_ghz=5e-5)
    assert len(intervals) == expected
