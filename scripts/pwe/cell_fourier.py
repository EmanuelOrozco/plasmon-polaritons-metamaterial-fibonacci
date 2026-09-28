#!/usr/bin/env python3
"""Perfiles ε(z), μ(z) de la celda S4 y su serie de Fourier truncada."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from fibonacci_photonics.studies.pwe import cell_fourier

if __name__ == "__main__":
    cell_fourier.main()
