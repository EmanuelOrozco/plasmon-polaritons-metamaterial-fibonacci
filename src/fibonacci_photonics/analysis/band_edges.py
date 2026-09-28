"""Localización de bandas con bordes refinados, válida para cualquier método.

Las subbandas plasmon-polaritón pueden ser mucho más estrechas que el paso de
la malla. Con R(ν) muestreada en una malla se detectan:

- bordes: tramos donde |R| − 1 cambia de signo;
- bandas enteras entre dos muestras: |R| > 1 en ambos extremos pero R cambia
  de signo (R pasa por 0 y por lo tanto por la banda |R| ≤ 1);
- gaps enteros entre dos muestras: dentro de una banda R(ν) es monótona, así
  que un extremo local de R entre muestras con |R| ≤ 1 indica que R tocó ±1;
  se maximiza |R| en ese tramo y, si supera 1, hay un gap que separa dos
  subbandas.

Cada borde se refina con Brent sobre |R(ν)| − 1 y los intervalos se arman
evaluando el signo entre bordes consecutivos.
"""

from __future__ import annotations

from collections.abc import Callable
from itertools import pairwise

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import minimize_scalar

from fibonacci_photonics.analysis.roots import DEFAULT_XTOL_GHZ, refine_crossing
from fibonacci_photonics.errors import InvalidParameterError
from fibonacci_photonics.solvers.scan import FrequencyInterval

RealFunction = Callable[[float], float]


def _hidden_gap_edges(
    r_of_nu: RealFunction,
    grid: NDArray[np.float64],
    r_grid: NDArray[np.float64],
    excess: RealFunction,
    xtol: float,
) -> list[float]:
    """Bordes de gaps más estrechos que la malla, rodeados de muestras en banda."""
    edges: list[float] = []
    inside = np.abs(r_grid) <= 1.0
    for i in range(1, grid.size - 1):
        if not inside[i]:
            continue
        slope_lo, slope_hi = r_grid[i] - r_grid[i - 1], r_grid[i + 1] - r_grid[i]
        if slope_lo * slope_hi > 0.0 or slope_lo == slope_hi == 0.0 or not np.isfinite(slope_lo * slope_hi):
            continue
        sign = float(np.sign(slope_lo)) if slope_lo != 0.0 else -float(np.sign(slope_hi))
        found = minimize_scalar(
            lambda nu, s=sign: -s * r_of_nu(nu),
            bounds=(grid[i - 1], grid[i + 1]),
            method="bounded",
            options={"xatol": max(xtol, 1.0e-15 * abs(grid[i + 1]))},
        )
        peak = float(found.x)
        if not (np.isfinite(found.fun) and -found.fun > 1.0):
            continue
        # El gap debe quedar entre dos muestras en banda; si no, lo ven los cambios de signo.
        j = i + 1 if peak > grid[i] else i - 1
        if not inside[j]:
            continue
        a, b = sorted((grid[i], grid[j]))
        lo = refine_crossing(excess, a, peak, xtol=xtol)
        hi = refine_crossing(excess, peak, b, xtol=xtol)
        if lo is not None and hi is not None:
            edges.extend((lo, hi))
    return edges


def locate_bands(
    r_of_nu: RealFunction,
    grid: ArrayLike,
    r_grid: ArrayLike,
    *,
    xtol: float = DEFAULT_XTOL_GHZ,
) -> list[FrequencyInterval]:
    """Bandas |R| ≤ 1 en [grid[0], grid[−1]] con bordes refinados por Brent.

    Parameters
    ----------
    r_of_nu : callable
        Re R(ν) escalar, ν en GHz.
    grid : array_like
        Malla creciente de frecuencias (GHz); puede ser no uniforme.
    r_grid : array_like
        Re R en la malla.
    xtol : float
        Tolerancia absoluta de los bordes, en GHz.
    """
    grid = np.asarray(grid, dtype=np.float64)
    r_grid = np.asarray(r_grid, dtype=np.float64)
    if grid.shape != r_grid.shape or grid.ndim != 1:
        raise InvalidParameterError("la malla y R deben ser 1-D y de igual tamaño")

    def excess(nu: float) -> float:
        return abs(r_of_nu(nu)) - 1.0

    f = np.abs(r_grid) - 1.0
    edges: list[float] = []
    for i in range(grid.size - 1):
        a, b = grid[i], grid[i + 1]
        fa, fb, ra, rb = f[i], f[i + 1], r_grid[i], r_grid[i + 1]
        if not (np.isfinite(fa) and np.isfinite(fb)):
            continue
        if fa * fb < 0.0:
            root = refine_crossing(excess, a, b, xtol=xtol)
            if root is not None:
                edges.append(root)
        elif fa > 0.0 and fb > 0.0 and ra * rb < 0.0:
            zero = refine_crossing(r_of_nu, a, b, xtol=xtol)
            if zero is None:
                continue
            lo = refine_crossing(excess, a, zero, xtol=xtol)
            hi = refine_crossing(excess, zero, b, xtol=xtol)
            if lo is not None and hi is not None:
                edges.extend((lo, hi))
    edges.extend(_hidden_gap_edges(r_of_nu, grid, r_grid, excess, xtol))

    points = np.unique(np.concatenate(([grid[0]], np.sort(edges), [grid[-1]])))
    # Un mismo gap oculto puede detectarse desde dos tripletas vecinas.
    points = points[np.concatenate(([True], np.diff(points) > 10.0 * xtol))]
    intervals: list[FrequencyInterval] = []
    for lo, hi in pairwise(points):
        if hi <= lo:
            continue
        if excess(0.5 * (lo + hi)) <= 0.0:
            if intervals and intervals[-1].nu_max_ghz == lo:
                lo = intervals.pop().nu_min_ghz
            intervals.append(FrequencyInterval(float(lo), float(hi), 0))
    return intervals


def refine_grid_interval(
    r_of_nu: RealFunction,
    grid: ArrayLike,
    interval: FrequencyInterval,
    *,
    xtol: float = DEFAULT_XTOL_GHZ,
) -> FrequencyInterval:
    """Lleva los bordes de un intervalo de malla (uniforme o no) al cruce |R| = 1.

    Los extremos de ``interval`` son muestras con |R| ≤ 1; cada cruce está entre
    esa muestra y su vecina exterior. Sin vecina (borde de la malla) o sin cambio
    de signo, el extremo se conserva.
    """
    grid = np.asarray(grid, dtype=np.float64)

    def excess(nu: float) -> float:
        return abs(r_of_nu(nu)) - 1.0

    i = int(np.searchsorted(grid, interval.nu_min_ghz))
    j = int(np.searchsorted(grid, interval.nu_max_ghz))
    lo = refine_crossing(excess, grid[i - 1], grid[i], xtol=xtol) if i > 0 else None
    hi = refine_crossing(excess, grid[j], grid[j + 1], xtol=xtol) if j + 1 < grid.size else None
    return FrequencyInterval(
        interval.nu_min_ghz if lo is None else lo,
        interval.nu_max_ghz if hi is None else hi,
        interval.n_points,
    )


def refine_band_edges(
    r_of_nu: RealFunction,
    interval: FrequencyInterval,
    step_ghz: float,
    *,
    xtol: float = DEFAULT_XTOL_GHZ,
) -> tuple[float, float]:
    """Bordes exactos de una banda de malla uniforme de paso ``step_ghz``.

    ``interval`` son las muestras extremas con |R| ≤ 1; los cruces están a menos
    de un paso hacia fuera.
    """

    def excess(nu: float) -> float:
        return abs(r_of_nu(nu)) - 1.0

    lo = refine_crossing(excess, interval.nu_min_ghz - step_ghz, interval.nu_min_ghz, xtol=xtol)
    hi = refine_crossing(excess, interval.nu_max_ghz, interval.nu_max_ghz + step_ghz, xtol=xtol)
    return (
        interval.nu_min_ghz if lo is None else lo,
        interval.nu_max_ghz if hi is None else hi,
    )


def validate_bands(
    bands: list[FrequencyInterval],
    r_of_nu: RealFunction,
    r_check: RealFunction,
    *,
    tolerance: float = 1.0e-2,
) -> tuple[list[FrequencyInterval], list[FrequencyInterval]]:
    """Separa bandas convergidas de espurias comparando dos discretizaciones.

    Una banda se acepta si en su centro la evaluación de control (p. ej. con más
    ondas planas) también da |R| ≤ 1 y difiere de la original en menos de
    ``tolerance``. Las bandas espurias del PWE, donde |Im k L_m| es grande, no
    sobreviven al cambio de truncamiento.
    """
    kept: list[FrequencyInterval] = []
    rejected: list[FrequencyInterval] = []
    for band in bands:
        center = band.center_ghz
        r1, r2 = r_of_nu(center), r_check(center)
        ok = np.isfinite(r1) and np.isfinite(r2) and abs(r2) <= 1.0 and abs(r1 - r2) <= tolerance
        (kept if ok else rejected).append(band)
    return kept, rejected


def plasmon_grid(nu_min: float, nu_m: float, points: int, min_offset: float) -> NDArray[np.float64]:
    """Malla logarítmica en ν_m − ν ∈ [min_offset, ν_m − nu_min] (GHz).

    Resuelve las subbandas que colapsan hacia ν_m, cuyo ancho decrece
    geométricamente al acercarse al plasmón magnético.
    """
    if not 0.0 < min_offset < nu_m - nu_min:
        raise InvalidParameterError(
            f"se requiere 0 < min_offset < ν_m − ν_min: {min_offset} y {nu_m - nu_min}"
        )
    offsets = np.logspace(np.log10(min_offset), np.log10(nu_m - nu_min), points)
    return np.sort(nu_m - offsets)
