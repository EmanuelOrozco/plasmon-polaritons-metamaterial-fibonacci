"""Estrategias para localizar subbandas plasmon-polaritón con bordes refinados.

Las subbandas bajo ν_m se estrechan geométricamente al acercarse al plasmón
magnético. Cada método usa la estrategia que su costo permite:

- ``dense_grid_bands`` (TMM): malla uniforme densa, vectorizada, y Brent solo
  en los bordes de cada intervalo; barato porque R cuesta microsegundos.
- ``adaptive_bands`` (PWE): malla logarítmica en ν_m − ν con pocos puntos,
  ``locate_bands`` para bandas y gaps más estrechos que la malla, y validación
  de cada banda con más ondas planas para descartar las espurias.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from fibonacci_photonics.analysis.band_edges import locate_bands, refine_grid_interval, validate_bands
from fibonacci_photonics.analysis.plasmon_modes import intervals_in_window
from fibonacci_photonics.physics.electromagnetics import Polarization
from fibonacci_photonics.solvers.base import DispersionSolver
from fibonacci_photonics.solvers.scan import FrequencyInterval, allowed_intervals


@dataclass(frozen=True)
class AdaptiveBands:
    """Subbandas aceptadas y descartadas por la validación, con su costo."""

    bands: list[FrequencyInterval]
    rejected: list[FrequencyInterval]
    evaluations: int


def dense_grid_bands(
    solver: DispersionSolver,
    m: int,
    theta: float,
    polarization: Polarization,
    grid: ArrayLike,
    window: tuple[float, float],
    *,
    min_points: int = 3,
) -> list[FrequencyInterval]:
    """Intervalos de la malla en ``window`` con sus bordes llevados a |R| = 1."""
    grid = np.asarray(grid, dtype=np.float64)
    scan = solver.scan(m, theta, grid, polarization)
    intervals = intervals_in_window(allowed_intervals(scan), window[0], window[1], min_points=min_points)
    r_of_nu = solver.evaluator(m, theta, polarization)
    return [refine_grid_interval(r_of_nu, grid, item) for item in intervals]


def adaptive_bands(
    solver: DispersionSolver,
    check: DispersionSolver,
    m: int,
    theta: float,
    polarization: Polarization,
    grid: ArrayLike,
) -> AdaptiveBands:
    """Bandas en una malla no uniforme, validadas con la discretización ``check``."""
    evaluator = solver.evaluator(m, theta, polarization)
    control = check.evaluator(m, theta, polarization)
    grid = np.asarray(grid, dtype=np.float64)
    r_grid = np.array([evaluator(float(nu)) for nu in grid], dtype=np.float64)
    candidates = locate_bands(evaluator, grid, r_grid)
    bands, rejected = validate_bands(candidates, evaluator, control)
    return AdaptiveBands(bands=bands, rejected=rejected, evaluations=evaluator.calls + control.calls)
