"""Anchos de banda en función del ángulo de incidencia."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np

from fibonacci_photonics.analysis.band_edges import refine_grid_interval
from fibonacci_photonics.analysis.plasmon_modes import intervals_in_window, merge_touching_intervals
from fibonacci_photonics.core.bands import DispersionScan, allowed_intervals
from fibonacci_photonics.core.electromagnetics import Polarization
from fibonacci_photonics.core.params import SuperlatticeSpec
from fibonacci_photonics.tmm.dispersion import scan_dispersion

Scanner = Callable[[SuperlatticeSpec, int, float, np.ndarray, Polarization], DispersionScan]


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
    merge_gap_ghz: float | None = 1.0e-5,
    scanner: Scanner = scan_dispersion,
    refine: bool = False,
) -> list[BandwidthPoint]:
    """Subbandas por ángulo; ``merge_gap_ghz=None`` desactiva la fusión de intervalos.

    ``scanner`` es el método que calcula R(ν) (TMM por defecto, o PWE). Con
    ``refine`` cada borde pasa de la muestra de malla al cruce |R| = 1 (Brent);
    sin refinar, el ancho queda subestimado hasta en dos pasos de malla.
    """
    points: list[BandwidthPoint] = []
    nu_min, nu_max = window
    for theta in thetas:
        scan = scanner(spec, m, float(theta), nu_ghz, polarization)
        intervals = intervals_in_window(
            allowed_intervals(scan), nu_min, nu_max, min_points=min_points
        )
        if refine:
            def r_of_nu(nu: float, theta: float = float(theta)) -> float:
                return float(scanner(spec, m, theta, np.array([nu]), polarization).r.real[0])

            intervals = [refine_grid_interval(r_of_nu, scan.nu_ghz, item) for item in intervals]
        if merge_gap_ghz is not None:
            intervals = merge_touching_intervals(intervals, gap_ghz=merge_gap_ghz)
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
