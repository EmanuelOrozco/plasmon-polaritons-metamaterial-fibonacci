"""Conteo de subbandas plasmon-polaritón en una ventana de frecuencia.

El paper predice F_{m−2} subbandas (una por capa B de la celda) bajo ν_m.
"""

from __future__ import annotations

from dataclasses import dataclass

from fibonacci_photonics.physics.fibonacci import n_layers_b
from fibonacci_photonics.solvers.scan import DispersionScan, FrequencyInterval, allowed_intervals


@dataclass(frozen=True)
class PlasmonModeReport:
    """Subbandas detectadas frente a las F_{m−2} esperadas."""

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
    """Intervalos con al menos ``min_points`` muestras que tocan [nu_min, nu_max] (GHz)."""
    return [
        interval
        for interval in intervals
        if interval.n_points >= min_points and interval.nu_max_ghz >= nu_min and interval.nu_min_ghz <= nu_max
    ]


def detect_plasmon_modes(
    scan: DispersionScan,
    *,
    nu_min_ghz: float,
    nu_max_ghz: float,
    min_points: int = 3,
) -> PlasmonModeReport:
    """Cuenta las subbandas conexas del barrido en la ventana [nu_min_ghz, nu_max_ghz]."""
    intervals = intervals_in_window(allowed_intervals(scan), nu_min_ghz, nu_max_ghz, min_points=min_points)
    return PlasmonModeReport(
        m=scan.m,
        expected_modes=n_layers_b(scan.m),
        detected_modes=len(intervals),
        intervals=tuple(intervals),
        window_ghz=(nu_min_ghz, nu_max_ghz),
    )


def bandwidths_ghz(intervals: list[FrequencyInterval]) -> list[float]:
    """Anchos Δν (GHz) de cada intervalo."""
    return [interval.bandwidth_ghz for interval in intervals]


def merge_touching_intervals(
    intervals: list[FrequencyInterval],
    *,
    gap_ghz: float = 1.0e-6,
) -> list[FrequencyInterval]:
    """Une intervalos separados por menos de ``gap_ghz`` (artefacto de malla)."""
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
