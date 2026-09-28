#!/usr/bin/env python3
"""Reproduce figuras de PRB 81, 153101 con la TMM, el PWE o ambos.

Ejemplos::

    python scripts/run_figure.py figure_01 --method tmm
    python scripts/run_figure.py 2 4 --method pwe --workers 8
    python scripts/run_figure.py all --method both
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from fibonacci_photonics.reproduction.cli import main

if __name__ == "__main__":
    main()
