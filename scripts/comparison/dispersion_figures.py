#!/usr/bin/env python3
"""Comparación TMM vs PWE de las Figuras 1–4 (requiere los barridos del PWE)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from fibonacci_photonics.comparison import dispersion_report

if __name__ == "__main__":
    dispersion_report.main()
