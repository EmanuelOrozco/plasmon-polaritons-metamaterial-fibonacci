#!/usr/bin/env python3
"""Figura 2 con la TMM (equivale a ``scripts/run_figure.py figure_02 --method tmm``)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from fibonacci_photonics.reproduction import run_figure

if __name__ == "__main__":
    run_figure("figure_02", "tmm")
