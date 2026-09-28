#!/usr/bin/env python3
"""Compatibilidad: ejecuta ``scripts/tmm/figure_02.py`` (Figura 2 con la TMM)."""

import runpy
from pathlib import Path

runpy.run_path(str(Path(__file__).resolve().parent / "tmm" / "figure_02.py"), run_name="__main__")
