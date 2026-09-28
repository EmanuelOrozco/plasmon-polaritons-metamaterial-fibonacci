"""Cierre de las curvas ν(k) en los bordes de banda, para cualquier método.

Una malla uniforme en ν casi nunca cae sobre un borde |R| = 1, donde la curva
debe llegar a k = 0 (R = +1) o k = π/L_m (R = −1): cerca del borde
k ∝ √|ν − ν_borde| y la última muestra queda a k finito. Aquí se refinan esos
puntos con el mismo método que produjo el barrido (Brent sobre su propia
R(ν)), sin mezclar métodos.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Union

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import brentq, minimize_scalar

from fibonacci_photonics.parallel import starmap_processes
from fibonacci_photonics.physics.electromagnetics import Polarization, as_polarization
from fibonacci_photonics.solvers.base import DispersionSolver
from fibonacci_photonics.solvers.scan import BAND_ABS_R_TOLERANCE, DispersionScan

EDGE, EXTREMUM, HIDDEN_GAP = 0, 1, 2
"""Tipos de punto en ``BandEdges.kind``."""

DEFAULT_EDGE_XTOL_FRACTION = 1.0e-6
"""Tolerancia de Brent para los bordes, como fracción del menor paso de la malla."""

EdgeTask = tuple[str, float, float, float]
ExtremumTask = tuple[str, float, float, float, float]
Task = Union[EdgeTask, ExtremumTask]
Point = tuple[float, float, float, int]


@dataclass(frozen=True)
class BandEdges:
    """Puntos de cierre de las bandas que la malla uniforme en ν no muestrea.

    - ``EDGE``: borde de banda, |R| = 1; con R = +1, k L_m = 0 y con R = −1, k L_m = π.
    - ``EXTREMUM``: extremo interior de |R| ≤ 1 (punta de la curva); si |R| toca 1
      sin abrir gap, como en incidencia normal con impedancias iguales, da k = 0 o 1.
    - ``HIDDEN_GAP``: máximo de |R| > 1 entre dos muestras en banda; k = NaN marca el
      gap, que queda delimitado por dos puntos ``EDGE``.

    Attributes
    ----------
    nu_ghz, r, k_lm_over_pi : ndarray
        Frecuencia (GHz), semitraza y k L_m/π de cada punto, ordenados en ν.
    kind : ndarray of int8
        ``EDGE``, ``EXTREMUM`` o ``HIDDEN_GAP``.
    evaluations : int
        Evaluaciones de R gastadas en el refinamiento.
    """

    nu_ghz: NDArray[np.float64]
    r: NDArray[np.float64]
    k_lm_over_pi: NDArray[np.float64]
    kind: NDArray[np.int8]
    evaluations: int

    @property
    def size(self) -> int:
        return int(self.nu_ghz.size)

    @property
    def n_edges(self) -> int:
        return int(np.count_nonzero(self.kind == EDGE))

    @property
    def n_extrema(self) -> int:
        return int(np.count_nonzero(self.kind == EXTREMUM))

    @property
    def n_hidden_gaps(self) -> int:
        return int(np.count_nonzero(self.kind == HIDDEN_GAP))

    @classmethod
    def empty(cls) -> BandEdges:
        nothing = np.array([], dtype=np.float64)
        return cls(nu_ghz=nothing, r=nothing, k_lm_over_pi=nothing, kind=np.array([], dtype=np.int8), evaluations=0)


def edge_brackets(scan: DispersionScan) -> list[tuple[float, float, float]]:
    """(ν en banda, ν en gap, signo de R en banda) para cada cambio banda/gap de la malla.

    Solo se aceptan pares con |R| − 1 de signo opuesto y valores finitos: los
    puntos excluidos junto a los polos de Drude o sin solución no son bordes.
    """
    nu = scan.nu_ghz
    r = np.real(scan.r)
    excess = np.abs(r) - 1.0
    allowed = scan.allowed
    brackets: list[tuple[float, float, float]] = []
    for i in np.flatnonzero(allowed[:-1] != allowed[1:]):
        inside, outside = (i, i + 1) if allowed[i] else (i + 1, i)
        if not (np.isfinite(excess[inside]) and np.isfinite(excess[outside])):
            continue
        if excess[inside] > 0.0 or excess[outside] <= 0.0:
            continue
        brackets.append((float(nu[inside]), float(nu[outside]), float(np.sign(r[inside]) or 1.0)))
    return brackets


def extremum_brackets(scan: DispersionScan) -> list[tuple[float, float, float, float]]:
    """(ν_{i−1}, ν_i, ν_{i+1}, |R_i|) para cada máximo local de |R| dentro de una banda.

    Ahí la curva k(ν) tiene su punta (k mínimo o máximo); el extremo verdadero cae
    entre las muestras vecinas y puede tocar k = 0 o 1, o esconder un gap estrecho.
    """
    nu = scan.nu_ghz
    abs_r = np.abs(np.real(scan.r))
    allowed = scan.allowed & np.isfinite(abs_r)
    inner = allowed[:-2] & allowed[1:-1] & allowed[2:]
    peak = (abs_r[1:-1] > abs_r[:-2]) & (abs_r[1:-1] >= abs_r[2:])
    return [
        (float(nu[i - 1]), float(nu[i]), float(nu[i + 1]), float(abs_r[i]))
        for i in np.flatnonzero(inner & peak) + 1
    ]


def _refine_chunk(
    solver: DispersionSolver,
    m: int,
    theta: float,
    polarization: Polarization,
    tasks: Sequence[Task],
    xtol: float,
) -> tuple[list[Point], int]:
    evaluator = solver.evaluator(m, theta, polarization)

    def excess(nu: float) -> float:
        return abs(evaluator(nu)) - 1.0

    def edge(lo: float, hi: float, fallback_sign: float) -> Point | None:
        f_lo, f_hi = excess(lo), excess(hi)
        if not (np.isfinite(f_lo) and np.isfinite(f_hi)) or f_lo * f_hi > 0.0:
            return None
        if f_lo == 0.0 or f_hi == 0.0:
            root = lo if f_lo == 0.0 else hi
        else:
            root = float(brentq(excess, lo, hi, xtol=xtol, rtol=4.0 * np.finfo(float).eps))
        # En una banda más estrecha que el paso, R puede cruzar casi todo [−1, 1]
        # entre dos muestras: el signo del borde se toma en la raíz, no en la muestra.
        r_root = evaluator(root)
        sign = float(np.sign(r_root)) if np.isfinite(r_root) and r_root != 0.0 else fallback_sign
        return (root, sign, 0.0 if sign > 0.0 else 1.0, EDGE)

    points: list[Point] = []
    for task in tasks:
        if task[0] == "edge":
            _, inside, outside, sign = task  # type: ignore[misc]
            lo, hi = sorted((inside, outside))
            point = edge(lo, hi, sign)
            if point is not None:
                points.append(point)
            continue
        _, lo, _mid, hi, abs_r_mid = task  # type: ignore[misc]
        best = minimize_scalar(
            lambda nu: -abs(evaluator(nu)), bounds=(lo, hi), method="bounded", options={"xatol": xtol}
        )
        nu_star = float(best.x)
        r_star = evaluator(nu_star)
        if not np.isfinite(r_star) or abs(r_star) <= abs_r_mid:
            continue
        if abs(r_star) <= 1.0 + BAND_ABS_R_TOLERANCE:
            k_star = float(np.arccos(np.clip(r_star, -1.0, 1.0)) / np.pi)
            points.append((nu_star, float(r_star), k_star, EXTREMUM))
            continue
        sign = float(np.sign(r_star))
        for a, b in ((lo, nu_star), (nu_star, hi)):
            point = edge(a, b, sign)
            if point is not None:
                points.append(point)
        points.append((nu_star, float(r_star), float("nan"), HIDDEN_GAP))
    return points, evaluator.calls


def close_band_edges(
    solver: DispersionSolver,
    scan: DispersionScan,
    *,
    workers: int = 1,
    xtol_fraction: float = DEFAULT_EDGE_XTOL_FRACTION,
) -> BandEdges:
    """Bordes, puntas y gaps ocultos de un barrido, con el método que lo produjo.

    - Bordes: Brent sobre |R(ν)| − 1 entre la última muestra en banda y la primera en gap.
    - Extremos: Brent acotado sobre −|R(ν)| alrededor de cada máximo local de |R|
      muestreado en banda; si el máximo supera 1 hay un gap oculto y se buscan sus
      dos bordes.

    Parameters
    ----------
    solver : DispersionSolver
        El método (y truncamiento) con que se calculó ``scan``.
    scan : DispersionScan
        Barrido en una malla creciente.
    workers : int
        Procesos; las tareas se reparten en bloques contiguos.
    xtol_fraction : float
        Tolerancia de Brent como fracción del menor paso de la malla.
    """
    polarization = as_polarization(scan.polarization)
    tasks: list[Task] = [("edge", *b) for b in edge_brackets(scan)]
    tasks += [("extremum", *b) for b in extremum_brackets(scan)]
    if not tasks:
        return BandEdges.empty()
    steps = np.diff(scan.nu_ghz)
    xtol = xtol_fraction * float(np.min(np.abs(steps[steps != 0.0]))) if steps.size else 1.0e-12
    if workers <= 1 or len(tasks) == 1:
        points, calls = _refine_chunk(solver, scan.m, scan.theta, polarization, tasks, xtol)
    else:
        split = np.array_split(np.arange(len(tasks)), min(workers, len(tasks)))
        parts = starmap_processes(
            _refine_chunk,
            [(solver, scan.m, scan.theta, polarization, [tasks[i] for i in chunk], xtol) for chunk in split],
            workers,
        )
        points = [point for part, _ in parts for point in part]
        calls = sum(part_calls for _, part_calls in parts)
    points.sort(key=lambda p: p[0])
    return BandEdges(
        nu_ghz=np.array([p[0] for p in points], dtype=np.float64),
        r=np.array([p[1] for p in points], dtype=np.float64),
        k_lm_over_pi=np.array([p[2] for p in points], dtype=np.float64),
        kind=np.array([p[3] for p in points], dtype=np.int8),
        evaluations=int(calls),
    )
