"""Anchos de banda de las subbandas en función del ángulo de incidencia."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from fibonacci_photonics.analysis.band_edges import refine_grid_interval
from fibonacci_photonics.analysis.plasmon_modes import intervals_in_window, merge_touching_intervals
from fibonacci_photonics.physics.electromagnetics import Polarization
from fibonacci_photonics.solvers.base import DispersionSolver
from fibonacci_photonics.solvers.scan import FrequencyInterval, allowed_intervals


@dataclass(frozen=True)
class BandwidthPoint:
    """Subbandas de un ángulo θ (rad): anchos, centros y bordes en GHz."""

    theta: float
    m: int
    bandwidths_ghz: tuple[float, ...]
    centers_ghz: tuple[float, ...]
    n_modes: int
    intervals: tuple[FrequencyInterval, ...] = ()


def bandwidth_versus_angle(
    solver: DispersionSolver,
    m: int,
    thetas: ArrayLike,
    nu_ghz: ArrayLike,
    polarization: Polarization,
    *,
    window: tuple[float, float],
    min_points: int = 3,
    merge_gap_ghz: float | None = 1.0e-5,
    refine: bool = False,
) -> list[BandwidthPoint]:
    """Subbandas por ángulo en la malla ``nu_ghz``.

    Parameters
    ----------
    solver : DispersionSolver
        Método que calcula R(ν).
    m : int
        Orden de Fibonacci.
    thetas : array_like
        Ángulos en rad.
    nu_ghz : array_like
        Malla creciente de frecuencias (GHz), uniforme o no.
    polarization : Polarization
        TE o TM.
    window : (float, float)
        Ventana de frecuencias (GHz) donde se cuentan las subbandas.
    min_points : int
        Muestras mínimas de un intervalo para contarlo.
    merge_gap_ghz : float or None
        Une intervalos separados por menos que esto; ``None`` no fusiona.
    refine : bool
        Lleva cada borde de la muestra de malla al cruce |R| = 1 (Brent); sin
        refinar, el ancho queda subestimado hasta en dos pasos de malla.
    """
    points: list[BandwidthPoint] = []
    nu = np.asarray(nu_ghz, dtype=np.float64)
    for theta in np.asarray(thetas, dtype=np.float64):
        scan = solver.scan(m, float(theta), nu, polarization)
        intervals = intervals_in_window(allowed_intervals(scan), window[0], window[1], min_points=min_points)
        if refine:
            r_of_nu = solver.evaluator(m, float(theta), polarization)
            intervals = [refine_grid_interval(r_of_nu, nu, item) for item in intervals]
        if merge_gap_ghz is not None:
            intervals = merge_touching_intervals(intervals, gap_ghz=merge_gap_ghz)
        points.append(
            BandwidthPoint(
                theta=float(theta),
                m=m,
                bandwidths_ghz=tuple(item.bandwidth_ghz for item in intervals),
                centers_ghz=tuple(item.center_ghz for item in intervals),
                n_modes=len(intervals),
                intervals=tuple(intervals),
            )
        )
    return points
