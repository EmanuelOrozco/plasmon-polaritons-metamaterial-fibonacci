#!/usr/bin/env python3
"""Eficiencia TMM vs PWE en un núcleo (``--plot-only`` redibuja desde el JSON)."""

import os
import sys
from pathlib import Path

# Un hilo de BLAS para que los tiempos sean comparables; debe fijarse antes de importar NumPy.
for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_name] = "1"

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from fibonacci_photonics.benchmark import efficiency  # noqa: E402

if __name__ == "__main__":
    efficiency.main(plot_only="--plot-only" in sys.argv)
