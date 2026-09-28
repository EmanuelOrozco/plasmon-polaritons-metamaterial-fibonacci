#!/usr/bin/env python3
"""Compatibilidad: ejecuta ``scripts/tmm/figure_05.py`` (Figura 5 con la TMM)."""

import runpy
from pathlib import Path

runpy.run_path(str(Path(__file__).resolve().parent / "tmm" / "figure_05.py"), run_name="__main__")
