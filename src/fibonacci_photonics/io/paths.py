"""Rutas del proyecto: raíz, configuraciones y resultados.

Organización de salidas: ``results/<método>/<figura>/{output,data}`` con
``<método>`` ∈ {tmm, pwe, comparison}. La raíz se toma de
``FIBONACCI_PHOTONICS_ROOT`` o, por defecto, del árbol del código fuente.
"""

from __future__ import annotations

import os
from pathlib import Path

from fibonacci_photonics.errors import InvalidParameterError

RESULT_METHODS = ("tmm", "pwe", "comparison")


def project_root() -> Path:
    """Raíz del repositorio (contiene ``configs/`` y ``results/``)."""
    env = os.environ.get("FIBONACCI_PHOTONICS_ROOT")
    if env:
        return Path(env).resolve()
    return Path(__file__).resolve().parents[3]


def configs_dir() -> Path:
    return project_root() / "configs"


def results_dir() -> Path:
    return project_root() / "results"


def result_dirs(method: str, name: str) -> dict[str, Path]:
    """Carpetas ``root``, ``data`` y ``output`` de ``results/<method>/<name>``, ya creadas."""
    if method not in RESULT_METHODS:
        raise InvalidParameterError(f"método de resultados desconocido: {method!r}")
    base = results_dir() / method / name
    dirs = {"root": base, "data": base / "data", "output": base / "output"}
    for path in dirs.values():
        path.mkdir(parents=True, exist_ok=True)
    return dirs


def relative(path: Path) -> str:
    """Ruta relativa a la raíz del repositorio (JSON portables)."""
    try:
        return str(path.resolve().relative_to(project_root()))
    except ValueError:
        return str(path)
