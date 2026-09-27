"""Utilidades compartidas por los scripts de reproducción de figuras."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from fibonacci_tmm.params import SuperlatticeSpec, spec_from_mapping
from fibonacci_tmm.provenance import run_metadata, write_metadata


def load_figure_config(name: str) -> tuple[SuperlatticeSpec, dict[str, Any], Path]:
    path = ROOT / "configs" / name
    with path.open(encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    return spec_from_mapping(data), data, path


def figure_dirs(figure_id: str) -> dict[str, Path]:
    base = ROOT / "figures" / figure_id
    dirs = {
        "root": base,
        "data": base / "data",
        "output": base / "output",
        "reproduced": ROOT / "figures" / "reproduced" / figure_id,
    }
    for path in dirs.values():
        path.mkdir(parents=True, exist_ok=True)
    return dirs


def save_scan_npz(path: Path, **arrays: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, **arrays)


def dump_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")


def write_run_sidecar(path: Path, spec: SuperlatticeSpec, extra: dict[str, Any]) -> None:
    write_metadata(path, run_metadata(spec, extra))


def write_figure_readme(path: Path, title: str, body: str) -> None:
    if path.exists():
        return
    path.write_text(f"# {title}\n\n{body.strip()}\n", encoding="utf-8")
