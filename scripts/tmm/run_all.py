#!/usr/bin/env python3
"""Reproduce las figuras 1–6 con la TMM y el estudio de convergencia de la Fig. 4."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCRIPTS = [f"figure_0{n}.py" for n in range(1, 7)] + ["convergence.py"]


def main() -> None:
    for name in SCRIPTS:
        print("=" * 60)
        print("TMM:", name)
        subprocess.run([sys.executable, str(HERE / name)], check=True)


if __name__ == "__main__":
    main()
