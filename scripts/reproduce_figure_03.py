#!/usr/bin/env python3
"""Compatibilidad: ejecuta ``scripts/tmm/figure_03.py`` (Figura 3 con la TMM)."""

import runpy
from pathlib import Path

runpy.run_path(str(Path(__file__).resolve().parent / "tmm" / "figure_03.py"), run_name="__main__")
