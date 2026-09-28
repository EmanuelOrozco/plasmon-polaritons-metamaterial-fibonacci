"""Medición de tiempos y comparación de métodos con una interfaz común."""

from __future__ import annotations

import time
from collections.abc import Callable, Iterator, Sequence
from contextlib import contextmanager
from dataclasses import dataclass
from functools import partial
from typing import TypeVar

import numpy as np
from numpy.typing import ArrayLike

from fibonacci_photonics.comparison.compare import ScanAgreement, compare_scans
from fibonacci_photonics.physics.electromagnetics import Polarization
from fibonacci_photonics.solvers.base import DispersionSolver
from fibonacci_photonics.solvers.scan import DispersionScan

T = TypeVar("T")


def timed(func: Callable[[], T], repeats: int = 1) -> tuple[T, float]:
    """Resultado de la última llamada y tiempo medio por llamada (s)."""
    if repeats < 1:
        raise ValueError(f"repeats debe ser ≥ 1: {repeats}")
    start = time.perf_counter()
    for _ in range(repeats):
        out = func()
    return out, (time.perf_counter() - start) / repeats


@contextmanager
def stopwatch(store: dict[str, float], key: str) -> Iterator[None]:
    """Acumula en ``store[key]`` el tiempo del bloque (s)."""
    start = time.perf_counter()
    try:
        yield
    finally:
        store[key] = store.get(key, 0.0) + time.perf_counter() - start


def power_fit(x: ArrayLike, y: ArrayLike) -> tuple[float, float]:
    """y ≈ c·x^α por mínimos cuadrados en log-log; devuelve (c, α)."""
    alpha, log_c = np.polyfit(np.log(np.asarray(x, dtype=float)), np.log(np.asarray(y, dtype=float)), 1)
    return float(np.exp(log_c)), float(alpha)


@dataclass(frozen=True)
class MethodRun:
    """Barrido de un método con su tiempo y su acuerdo con la referencia."""

    label: str
    scan: DispersionScan
    seconds: float
    agreement: ScanAgreement | None

    @property
    def seconds_per_frequency(self) -> float:
        return self.seconds / self.scan.nu_ghz.size


def benchmark_methods(
    solvers: Sequence[DispersionSolver],
    m: int,
    theta: float,
    nu_ghz: ArrayLike,
    polarization: Polarization = Polarization.TE,
    *,
    repeats: Sequence[int] | None = None,
) -> list[MethodRun]:
    """Barre la misma malla con cada método y compara todos con el primero.

    Parameters
    ----------
    solvers : sequence of DispersionSolver
        El primero es la referencia (normalmente ``TMMSolver``).
    m, theta, nu_ghz, polarization
        Caso a comparar: orden, ángulo (rad), malla creciente (GHz) y polarización.
    repeats : sequence of int, optional
        Repeticiones por método para el tiempo medio (1 por defecto).

    Examples
    --------
    >>> runs = benchmark_methods([TMMSolver(spec), PWESolver(spec, 8)], 3, 0.5, nu)  # doctest: +SKIP
    >>> runs[1].agreement.max_abs_delta_r_in_bands  # doctest: +SKIP
    """
    if not solvers:
        raise ValueError("se necesita al menos un método")
    counts = list(repeats) if repeats is not None else [1] * len(solvers)
    if len(counts) != len(solvers):
        raise ValueError("repeats debe tener un valor por método")
    runs: list[MethodRun] = []
    for solver, count in zip(solvers, counts, strict=True):
        scan, seconds = timed(partial(solver.scan, m, theta, nu_ghz, polarization), count)
        agreement = compare_scans(runs[0].scan, scan) if runs else None
        runs.append(MethodRun(label=getattr(solver, "label", solver.method), scan=scan, seconds=seconds,
                              agreement=agreement))
    return runs
