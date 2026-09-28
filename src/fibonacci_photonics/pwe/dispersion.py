"""Relación de dispersión cos(k Lm) = Rm obtenida con el método de ondas planas."""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import brentq, minimize_scalar

from fibonacci_photonics.core.bands import DispersionScan, omega_from_nu_ghz, scan_from_semitrace
from fibonacci_photonics.core.constants import BAND_ABS_R_TOLERANCE
from fibonacci_photonics.core.electromagnetics import Polarization
from fibonacci_photonics.core.params import SuperlatticeSpec
from fibonacci_photonics.pwe.bloch import (
    FactorizationRule,
    _single_threaded_blas,
    bloch_eigenvalues,
    n_max_for,
    select_bloch_wavevector,
    solve_bloch,
)
from fibonacci_photonics.pwe.fourier import cell_fourier

DEFAULT_HARMONICS_PER_LAYER = 16
DEFAULT_EDGE_XTOL_FRACTION = 1.0e-6
"""Tolerancia de Brent para los bordes, como fracción del paso de la malla."""


class SemitraceEvaluator:
    """R_PWE(ν) escalar con la celda de Fourier precalculada (para Brent)."""

    def __init__(
        self,
        spec: SuperlatticeSpec,
        m: int,
        theta: float,
        polarization: Polarization = Polarization.TE,
        *,
        n_max: int,
        rule: FactorizationRule = "inverse",
    ) -> None:
        self.spec = spec
        self.theta = theta
        self.polarization = polarization
        self.rule = rule
        self.cell = cell_fourier(m, spec.thickness_a, spec.thickness_b, n_max)
        self.calls = 0

    def __call__(self, nu_ghz: float) -> float:
        self.calls += 1
        omega = float(omega_from_nu_ghz(nu_ghz))
        try:
            eigenvalues = bloch_eigenvalues(
                self.spec, self.cell, omega, self.theta, self.polarization, self.rule
            )
        except np.linalg.LinAlgError:
            return float("nan")
        eigenvalues = eigenvalues[np.isfinite(eigenvalues)]
        if eigenvalues.size == 0:
            return float("nan")
        return float(np.cos(select_bloch_wavevector(eigenvalues)).real)

    def many(self, nu_ghz: ArrayLike) -> np.ndarray:
        return np.array([self(float(nu)) for nu in np.atleast_1d(nu_ghz)])


def pwe_semitrace(
    spec: SuperlatticeSpec,
    m: int,
    theta: float,
    nu_ghz: ArrayLike,
    polarization: Polarization = Polarization.TE,
    *,
    harmonics_per_layer: int = DEFAULT_HARMONICS_PER_LAYER,
    n_max: int | None = None,
    rule: FactorizationRule = "inverse",
    workers: int = 1,
) -> np.ndarray:
    """R_PWE(ν) = Re cos(k̃), con k̃ el vector de Bloch físico del QEP."""
    omega = omega_from_nu_ghz(np.asarray(nu_ghz, dtype=np.float64))
    n_max = n_max_for(m, harmonics_per_layer) if n_max is None else int(n_max)
    return solve_bloch(
        spec, m, omega, theta, polarization, n_max=n_max, rule=rule, workers=workers
    ).r


def scan_dispersion(
    spec: SuperlatticeSpec,
    m: int,
    theta: float,
    nu_ghz: ArrayLike,
    polarization: Polarization = Polarization.TE,
    *,
    harmonics_per_layer: int = DEFAULT_HARMONICS_PER_LAYER,
    n_max: int | None = None,
    rule: FactorizationRule = "inverse",
    workers: int = 1,
    atol: float = BAND_ABS_R_TOLERANCE,
) -> DispersionScan:
    """Mismo contrato que ``fibonacci_photonics.tmm.dispersion.scan_dispersion``."""
    nu = np.asarray(nu_ghz, dtype=np.float64)
    r = pwe_semitrace(
        spec,
        m,
        theta,
        nu,
        polarization,
        harmonics_per_layer=harmonics_per_layer,
        n_max=n_max,
        rule=rule,
        workers=workers,
    )
    return scan_from_semitrace(spec, m, theta, nu, r, polarization, f"pwe-{rule}", atol=atol)


EDGE, EXTREMUM, HIDDEN_GAP = 0, 1, 2
"""Tipos de punto en ``BandEdges.kind``."""


@dataclass(frozen=True)
class BandEdges:
    """Puntos de cierre de las bandas del PWE que la malla uniforme en ν no muestrea.

    - ``EDGE``: borde de banda, |R| = 1; con R = +1, k Lm = 0 y con R = −1, k Lm = π.
      Cerca del borde k ∝ √|ν − ν_borde|, así que la última muestra queda a k finito.
    - ``EXTREMUM``: extremo interior de |R| ≤ 1 (punta de la curva); si |R| toca 1
      sin abrir gap, como en incidencia normal con impedancias iguales, da k = 0 o 1.
    - ``HIDDEN_GAP``: máximo de |R| > 1 entre dos muestras en banda; k = NaN marca el
      gap, que queda delimitado por dos puntos ``EDGE``.
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
    brackets: list[tuple[float, float, float, float]] = []
    for i in range(1, nu.size - 1):
        if not (allowed[i - 1] and allowed[i] and allowed[i + 1]):
            continue
        if abs_r[i] > abs_r[i - 1] and abs_r[i] >= abs_r[i + 1]:
            brackets.append((float(nu[i - 1]), float(nu[i]), float(nu[i + 1]), float(abs_r[i])))
    return brackets


Point = tuple[float, float, float, int]


def _refine_edges_chunk(
    spec: SuperlatticeSpec,
    m: int,
    theta: float,
    polarization: Polarization,
    n_max: int,
    rule: FactorizationRule,
    tasks: list[tuple],
    xtol: float,
) -> tuple[list[Point], int]:
    evaluator = SemitraceEvaluator(spec, m, theta, polarization, n_max=n_max, rule=rule)

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
            _, inside, outside, sign = task
            point = edge(*sorted((inside, outside)), sign)
            if point is not None:
                points.append(point)
            continue
        _, lo, _mid, hi, abs_r_mid = task
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


def pwe_band_edges(
    spec: SuperlatticeSpec,
    m: int,
    theta: float,
    scan: DispersionScan,
    polarization: Polarization = Polarization.TE,
    *,
    harmonics_per_layer: int = DEFAULT_HARMONICS_PER_LAYER,
    n_max: int | None = None,
    rule: FactorizationRule = "inverse",
    workers: int = 1,
    xtol_fraction: float = DEFAULT_EDGE_XTOL_FRACTION,
) -> BandEdges:
    """Puntos de cierre de las bandas de un barrido PWE, con el mismo truncamiento.

    - Bordes: Brent sobre |R_PWE(ν)| − 1 entre la última muestra en banda y la
      primera en gap.
    - Extremos: Brent acotado sobre −|R_PWE(ν)| alrededor de cada máximo local de
      |R| muestreado en banda; si el máximo supera 1 hay un gap oculto y se buscan
      sus dos bordes.

    Solo interviene R_PWE; no se usa ningún otro método.
    """
    n_max = n_max_for(m, harmonics_per_layer) if n_max is None else int(n_max)
    tasks: list[tuple] = [("edge", *b) for b in edge_brackets(scan)]
    tasks += [("extremum", *b) for b in extremum_brackets(scan)]
    if not tasks:
        empty = np.array([], dtype=np.float64)
        return BandEdges(nu_ghz=empty, r=empty, k_lm_over_pi=empty,
                         kind=np.array([], dtype=np.int8), evaluations=0)
    steps = np.diff(scan.nu_ghz)
    xtol = xtol_fraction * float(np.min(np.abs(steps[steps != 0.0]))) if steps.size else 1.0e-12
    if workers <= 1 or len(tasks) == 1:
        points, calls = _refine_edges_chunk(spec, m, theta, polarization, n_max, rule, tasks, xtol)
    else:
        _single_threaded_blas()
        split = np.array_split(np.arange(len(tasks)), min(workers, len(tasks)))
        chunks = [[tasks[i] for i in c] for c in split]
        n = len(chunks)
        with ProcessPoolExecutor(max_workers=workers) as pool:
            parts = list(pool.map(
                _refine_edges_chunk,
                [spec] * n, [m] * n, [theta] * n, [polarization] * n, [n_max] * n, [rule] * n,
                chunks, [xtol] * n,
            ))
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
