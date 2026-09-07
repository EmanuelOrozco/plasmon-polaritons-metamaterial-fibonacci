"""Metadatos de reproducibilidad para figuras y datos."""

from __future__ import annotations

import json
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from fibonacci_tmm.constants import BAND_ABS_R_TOLERANCE, SPEED_OF_LIGHT
from fibonacci_tmm.params import SuperlatticeSpec, spec_as_dict


def git_commit() -> str | None:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            check=False,
            capture_output=True,
            text=True,
            cwd=Path(__file__).resolve().parents[2],
        )
    except OSError:
        return None
    if result.returncode != 0:
        return None
    return result.stdout.strip() or None


def run_metadata(spec: SuperlatticeSpec, extra: dict[str, Any] | None = None) -> dict[str, Any]:
    payload: dict[str, Any] = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "numpy_version": _package_version("numpy"),
        "git_commit": git_commit(),
        "speed_of_light": spec.speed_of_light,
        "speed_of_light_note": "IMPL: el paper no especifica c; se usa SI salvo que el YAML diga otra cosa",
        "band_abs_r_tolerance": BAND_ABS_R_TOLERANCE,
        "q_convention": "q = (omega/c) n_A sin(theta)  [omega/c restaurado]",
        "spec": _jsonable(spec_as_dict(spec)),
        "author": "Emanuel Orozco Gallego",
    }
    if extra:
        payload["run"] = extra
    return payload


def write_metadata(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")


def _package_version(name: str) -> str:
    try:
        from importlib.metadata import version

        return version(name)
    except Exception:
        return "unknown"


def _jsonable(value: Any) -> Any:
    if isinstance(value, complex):
        return {"real": value.real, "imag": value.imag}
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    return value
