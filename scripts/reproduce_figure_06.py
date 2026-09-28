#!/usr/bin/env python3
"""Compatibilidad: ejecuta ``scripts/tmm/figure_06.py`` (Figura 6 con la TMM)."""

import runpy
from pathlib import Path

runpy.run_path(str(Path(__file__).resolve().parent / "tmm" / "figure_06.py"), run_name="__main__")
