#!/usr/bin/env python3
"""Todos los resultados del PWE: benchmark del libro, estudios numéricos y Figs. 1–6.

Usa todos los núcleos por defecto; ``PWE_WORKERS=n`` limita el número de procesos.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCRIPTS = [
    "book_benchmark.py",
    "cell_fourier.py",
    "spectral_pollution.py",
    "convergence.py",
    *[f"figure_0{n}.py" for n in range(1, 7)],
]


def main() -> None:
    for name in SCRIPTS:
        print("=" * 60)
        print("PWE:", name)
        subprocess.run([sys.executable, "-u", str(HERE / name)], check=True)


if __name__ == "__main__":
    main()
