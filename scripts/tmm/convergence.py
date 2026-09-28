#!/usr/bin/env python3
"""Convergencia del ancho de banda de la Fig. 4 (TMM) respecto a la malla en frecuencia."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from fibonacci_photonics.studies.tmm import convergence as tmm_convergence

if __name__ == "__main__":
    tmm_convergence.main()
