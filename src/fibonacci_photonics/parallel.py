"""Paralelismo por procesos para los barridos costosos del PWE.

Cada frecuencia (o cada caso m, θ) es independiente, así que el trabajo se
reparte en procesos con ``ProcessPoolExecutor``. Los procesos hijos heredan el
entorno: con un hilo de BLAS por proceso se evita la sobre-suscripción de
núcleos. El número de procesos se toma de ``PWE_WORKERS`` o, si no está
definido, de todos los núcleos.
"""

from __future__ import annotations

import os
from collections.abc import Callable, Iterable, Sequence
from concurrent.futures import ProcessPoolExecutor
from typing import TypeVar

T = TypeVar("T")
R = TypeVar("R")

BLAS_THREAD_VARIABLES = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")


def single_threaded_blas() -> None:
    """Fija un hilo de BLAS para este proceso y sus hijos (si no se fijó antes)."""
    for name in BLAS_THREAD_VARIABLES:
        os.environ.setdefault(name, "1")


def worker_count(requested: int | None = None) -> int:
    """Procesos a usar: ``requested``, o ``PWE_WORKERS``, o todos los núcleos."""
    if requested is not None:
        return max(1, int(requested))
    env = os.environ.get("PWE_WORKERS")
    if env:
        return max(1, int(env))
    return max(1, os.cpu_count() or 1)


def parallel_map(func: Callable[[T], R], tasks: Sequence[T], workers: int | None = None) -> list[R]:
    """``map`` en procesos, preservando el orden; en serie si hay un solo proceso o tarea."""
    single_threaded_blas()
    n = min(worker_count(workers), max(1, len(tasks)))
    if n == 1:
        return [func(task) for task in tasks]
    with ProcessPoolExecutor(max_workers=n) as pool:
        return list(pool.map(func, tasks))


def starmap_processes(func: Callable[..., R], argument_lists: Iterable[Sequence[object]], workers: int) -> list[R]:
    """``func(*args)`` para cada lista de argumentos, en ``workers`` procesos."""
    single_threaded_blas()
    columns = list(zip(*argument_lists, strict=True))
    with ProcessPoolExecutor(max_workers=workers) as pool:
        return list(pool.map(func, *columns))
