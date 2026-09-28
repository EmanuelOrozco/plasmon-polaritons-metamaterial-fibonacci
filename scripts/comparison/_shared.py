"""Rutas comunes de los scripts de comparación TMM vs PWE."""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
for path in (ROOT / "src", ROOT / "scripts"):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

PWE_RESULTS = ROOT / "results" / "pwe"
TMM_RESULTS = ROOT / "results" / "tmm"
