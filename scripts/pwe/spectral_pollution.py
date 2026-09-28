#!/usr/bin/env python3
"""Forma ω(k) del libro extendida a Drude: contaminación espectral."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from fibonacci_photonics.studies import spectral_pollution

if __name__ == "__main__":
    spectral_pollution.main()
