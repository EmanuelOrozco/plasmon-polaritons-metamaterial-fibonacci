#!/usr/bin/env python3
"""Comparación completa TMM vs PWE (requiere haber corrido scripts/tmm y scripts/pwe)."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCRIPTS = ["dispersion_figures.py", "plasmon_figures.py", "efficiency.py", "summary.py"]


def main() -> None:
    for name in SCRIPTS:
        print("=" * 60)
        print("Comparación:", name)
        subprocess.run([sys.executable, "-u", str(HERE / name)], check=True)


if __name__ == "__main__":
    main()
