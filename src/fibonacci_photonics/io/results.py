"""Escritura y lectura de resultados: npz, JSON y README de cada figura."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from fibonacci_photonics.io.provenance import run_metadata
from fibonacci_photonics.physics.superlattice import SuperlatticeSpec


def save_npz(path: Path, **arrays: Any) -> None:
    """Arreglos comprimidos en ``path``."""
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, **arrays)


def dump_json(path: Path, payload: dict[str, Any]) -> None:
    """JSON indentado; lo no serializable se escribe con ``str``."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")


def load_json(path: Path) -> dict[str, Any]:
    data: dict[str, Any] = json.loads(path.read_text(encoding="utf-8"))
    return data


def write_run_metadata(path: Path, spec: SuperlatticeSpec, extra: dict[str, Any]) -> None:
    """``metadata.json`` con entorno, commit y parámetros de la corrida."""
    dump_json(path, run_metadata(spec, extra))


def write_readme(path: Path, title: str, body: str) -> None:
    """README de la carpeta de una figura (no sobrescribe uno existente)."""
    if path.exists():
        return
    path.write_text(f"# {title}\n\n{body.strip()}\n", encoding="utf-8")
