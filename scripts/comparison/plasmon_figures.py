#!/usr/bin/env python3
"""Comparación TMM vs PWE de las Figuras 5 y 6 (bordes de subbandas refinados)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from fibonacci_photonics.comparison import plasmon_report

if __name__ == "__main__":
    plasmon_report.main()
