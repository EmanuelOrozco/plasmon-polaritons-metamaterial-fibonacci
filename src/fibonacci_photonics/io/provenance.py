"""Metadatos de reproducibilidad para figuras y datos."""

from __future__ import annotations

import platform
import subprocess
import sys
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from typing import Any

from fibonacci_photonics.io.paths import project_root
from fibonacci_photonics.physics.superlattice import SuperlatticeSpec, spec_as_dict
from fibonacci_photonics.solvers.scan import BAND_ABS_R_TOLERANCE


def git_commit() -> str | None:
    """Commit actual del repositorio, o ``None`` fuera de git."""
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            check=False,
            capture_output=True,
            text=True,
            cwd=project_root(),
        )
    except OSError:
        return None
    if result.returncode != 0:
        return None
    return result.stdout.strip() or None


def _package_version(name: str) -> str:
    try:
        return version(name)
    except PackageNotFoundError:
        return "unknown"


def jsonable(value: Any) -> Any:
    """Convierte complejos y tuplas anidados en tipos JSON."""
    if isinstance(value, complex):
        return {"real": value.real, "imag": value.imag}
    if isinstance(value, dict):
        return {str(key): jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(item) for item in value]
    return value


def run_metadata(spec: SuperlatticeSpec, extra: dict[str, Any] | None = None) -> dict[str, Any]:
    """Entorno, versión del código y parámetros físicos de una corrida."""
    payload: dict[str, Any] = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "numpy_version": _package_version("numpy"),
        "scipy_version": _package_version("scipy"),
        "git_commit": git_commit(),
        "speed_of_light": spec.speed_of_light,
        "speed_of_light_note": "IMPL: el paper no especifica c; se usa SI salvo que el YAML diga otra cosa",
        "band_abs_r_tolerance": BAND_ABS_R_TOLERANCE,
        "q_convention": "q = (omega/c) n_A sin(theta)  [omega/c restaurado]",
        "spec": jsonable(spec_as_dict(spec)),
        "author": "Emanuel Orozco Gallego",
    }
    if extra:
        payload["run"] = extra
    return payload
