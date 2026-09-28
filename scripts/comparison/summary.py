#!/usr/bin/env python3
"""Resumen global TMM vs PWE: JSON, tablas LaTeX y figura de tiempos."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from fibonacci_photonics.benchmark import summary

if __name__ == "__main__":
    summary.main()
