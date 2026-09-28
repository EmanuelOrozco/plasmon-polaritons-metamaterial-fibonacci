"""Detección automática de modos plasmon-polaritón y anchos de banda."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from fibonacci_photonics.core.bands import DispersionScan, FrequencyInterval, allowed_intervals
from fibonacci_photonics.core.fibonacci import n_layers_b


@dataclass(frozen=True)
class PlasmonModeReport:
    m: int
    expected_modes: int
    detected_modes: int
    intervals: tuple[FrequencyInterval, ...]
    window_ghz: tuple[float, float]


def intervals_in_window(
    intervals: list[FrequencyInterval],
    nu_min: float,
    nu_max: float,
    *,
    min_points: int = 3,
) -> list[FrequencyInterval]:
    selected = [
        interval
        for interval in intervals
        if interval.n_points >= min_points
        and interval.nu_max_ghz >= nu_min
        and interval.nu_min_ghz <= nu_max
    ]
    return selected


def detect_plasmon_modes(
    scan: DispersionScan,
    *,
    nu_min_ghz: float,
    nu_max_ghz: float,
    min_points: int = 3,
) -> PlasmonModeReport:
    """Cuenta subbandas conexas en una ventana de frecuencia."""
    intervals = intervals_in_window(
        allowed_intervals(scan), nu_min_ghz, nu_max_ghz, min_points=min_points
    )
    expected = n_layers_b(scan.m)
    return PlasmonModeReport(
        m=scan.m,
        expected_modes=expected,
        detected_modes=len(intervals),
        intervals=tuple(intervals),
        window_ghz=(nu_min_ghz, nu_max_ghz),
    )


def bandwidths_ghz(intervals: list[FrequencyInterval]) -> list[float]:
    return [interval.bandwidth_ghz for interval in intervals]


def merge_touching_intervals(
    intervals: list[FrequencyInterval],
    *,
    gap_ghz: float = 1.0e-6,
) -> list[FrequencyInterval]:
    """Une intervalos separados por menos de gap_ghz (artefacto de malla)."""
    if not intervals:
        return []
    ordered = sorted(intervals, key=lambda item: item.nu_min_ghz)
    merged = [ordered[0]]
    for current in ordered[1:]:
        previous = merged[-1]
        if current.nu_min_ghz - previous.nu_max_ghz <= gap_ghz:
            merged[-1] = FrequencyInterval(
                nu_min_ghz=previous.nu_min_ghz,
                nu_max_ghz=max(previous.nu_max_ghz, current.nu_max_ghz),
                n_points=previous.n_points + current.n_points,
            )
        else:
            merged.append(current)
    return merged
