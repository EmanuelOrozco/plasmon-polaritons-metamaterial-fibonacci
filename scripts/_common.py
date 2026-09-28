"""Utilidades compartidas por los scripts de TMM, PWE y comparación.

Organización de salidas: ``results/<método>/<figura>/{output,data}`` con
``<método>`` ∈ {tmm, pwe, comparison}.
"""

from __future__ import annotations

import json
import os
import sys
import time
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import Any

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from fibonacci_photonics.core.params import SuperlatticeSpec, spec_from_mapping  # noqa: E402
from fibonacci_photonics.core.provenance import run_metadata, write_metadata  # noqa: E402

RESULTS = ROOT / "results"
METHODS = ("tmm", "pwe", "comparison")


def load_figure_config(name: str) -> tuple[SuperlatticeSpec, dict[str, Any], Path]:
    """Parámetros físicos (y malla TMM) de ``configs/<name>``."""
    path = ROOT / "configs" / name
    with path.open(encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    return spec_from_mapping(data), data, path


def load_pwe_config(figure_id: str) -> tuple[SuperlatticeSpec, dict[str, Any], dict[str, Any], Path]:
    """Física de ``configs/<figura>.yaml`` + parámetros numéricos de ``configs/pwe/<figura>.yaml``."""
    spec, physics, _ = load_figure_config(f"{figure_id}.yaml")
    path = ROOT / "configs" / "pwe" / f"{figure_id}.yaml"
    with path.open(encoding="utf-8") as handle:
        numerics = yaml.safe_load(handle)
    return spec, physics, numerics, path


def result_dirs(method: str, figure_id: str) -> dict[str, Path]:
    if method not in METHODS:
        raise ValueError(f"método desconocido: {method!r}")
    base = RESULTS / method / figure_id
    dirs = {"root": base, "data": base / "data", "output": base / "output"}
    for path in dirs.values():
        path.mkdir(parents=True, exist_ok=True)
    return dirs


def figure_dirs(figure_id: str) -> dict[str, Path]:
    """Directorios de la reproducción TMM (``results/tmm/<figura>``)."""
    return result_dirs("tmm", figure_id)


def save_scan_npz(path: Path, **arrays: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, **arrays)


def dump_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def relative(path: Path) -> str:
    """Ruta relativa a la raíz del repositorio (para JSON portables)."""
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path)


def write_run_sidecar(path: Path, spec: SuperlatticeSpec, extra: dict[str, Any]) -> None:
    write_metadata(path, run_metadata(spec, extra))


def write_figure_readme(path: Path, title: str, body: str) -> None:
    if path.exists():
        return
    path.write_text(f"# {title}\n\n{body.strip()}\n", encoding="utf-8")


def worker_count(requested: int | None = None) -> int:
    """Procesos para barridos PWE: ``PWE_WORKERS`` o todos los núcleos."""
    if requested is not None:
        return max(1, int(requested))
    env = os.environ.get("PWE_WORKERS")
    if env:
        return max(1, int(env))
    return max(1, os.cpu_count() or 1)


def single_threaded_blas() -> None:
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ.setdefault(name, "1")


def parallel_map(func, tasks: list, workers: int | None = None) -> list:
    """``map`` en procesos (un hilo de BLAS por proceso), preservando el orden."""
    from concurrent.futures import ProcessPoolExecutor

    single_threaded_blas()
    n = min(worker_count(workers), max(1, len(tasks)))
    if n == 1:
        return [func(task) for task in tasks]
    with ProcessPoolExecutor(max_workers=n) as pool:
        return list(pool.map(func, tasks))


@contextmanager
def stopwatch(store: dict[str, float], key: str) -> Iterator[None]:
    start = time.perf_counter()
    try:
        yield
    finally:
        store[key] = store.get(key, 0.0) + time.perf_counter() - start
