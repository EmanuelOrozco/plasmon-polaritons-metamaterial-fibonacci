#!/usr/bin/env python3
"""Compatibilidad: ejecuta ``scripts/tmm/run_all.py`` (figuras 1–6 con la TMM)."""

import runpy
from pathlib import Path

runpy.run_path(str(Path(__file__).resolve().parent / "tmm" / "run_all.py"), run_name="__main__")
