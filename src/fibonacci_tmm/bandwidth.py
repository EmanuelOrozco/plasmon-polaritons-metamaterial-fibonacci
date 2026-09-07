"""Anchos de banda en función del ángulo de incidencia."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from fibonacci_tmm.dispersion import FrequencyInterval, scan_dispersion
from fibonacci_tmm.electromagnetics import Polarization
from fibonacci_tmm.params import SuperlatticeSpec
from fibonacci_tmm.plasmon_modes import (
    detect_plasmon_modes,
    intervals_in_window,
    merge_touching_intervals,
)
from fibonacci_tmm.dispersion import allowed_intervals


@dataclass(frozen=True)
class BandwidthPoint:
    theta: float
    m: int
    bandwidths_ghz: tuple[float, ...]
    centers_ghz: tuple[float, ...]
    n_modes: int


def bandwidth_versus_angle(
    spec: SuperlatticeSpec,
    m: int,
    thetas: np.ndarray,
    nu_ghz: np.ndarray,
    polarization: Polarization,
    *,
    window: tuple[float, float],
    min_points: int = 3,
    merge_gap_ghz: float = 1.0e-5,
) -> list[BandwidthPoint]:
    points: list[BandwidthPoint] = []
    nu_min, nu_max = window
    for theta in thetas:
        scan = scan_dispersion(spec, m, float(theta), nu_ghz, polarization)
        intervals = merge_touching_intervals(
            intervals_in_window(allowed_intervals(scan), nu_min, nu_max, min_points=min_points),
            gap_ghz=merge_gap_ghz,
        )
        points.append(
            BandwidthPoint(
                theta=float(theta),
                m=m,
                bandwidths_ghz=tuple(item.bandwidth_ghz for item in intervals),
                centers_ghz=tuple(item.center_ghz for item in intervals),
                n_modes=len(intervals),
            )
        )
    return points
