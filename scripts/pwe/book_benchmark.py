#!/usr/bin/env python3
"""Benchmark del capítulo 1D de Sukhoivanov y Guryev con PWE y TMM."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from fibonacci_photonics.studies.pwe import book_benchmark

if __name__ == "__main__":
    book_benchmark.main()
