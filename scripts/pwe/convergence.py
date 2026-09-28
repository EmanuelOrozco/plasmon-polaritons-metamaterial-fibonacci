#!/usr/bin/env python3
"""Convergencia, costo y límite de validez del PWE (``--plot-only`` redibuja sin recalcular)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from fibonacci_photonics.studies.pwe import convergence as pwe_convergence

if __name__ == "__main__":
    pwe_convergence.main(plot_only="--plot-only" in sys.argv)
