#!/usr/bin/env python3
"""Genera todas las figuras 1–6."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = [
    "reproduce_figure_01.py",
    "reproduce_figure_02.py",
    "reproduce_figure_03.py",
    "reproduce_figure_04.py",
    "reproduce_figure_05.py",
    "reproduce_figure_06.py",
]


def main() -> None:
    python = sys.executable
    for name in SCRIPTS:
        print("=" * 60)
        print("Ejecutando", name)
        subprocess.run([python, str(ROOT / "scripts" / name)], check=True)


if __name__ == "__main__":
    main()
