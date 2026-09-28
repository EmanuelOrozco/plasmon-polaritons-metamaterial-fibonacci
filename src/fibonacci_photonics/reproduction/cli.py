"""Reproduce figuras de PRB 81, 153101 con la TMM, el PWE o ambos.

Ejemplos::

    fibonacci-figure figure_01 --method tmm
    fibonacci-figure 2 4 --method pwe --workers 8
    fibonacci-figure all --method both
"""

from __future__ import annotations

import argparse

from fibonacci_photonics.reproduction import FIGURES, METHODS, normalize_figure_id, run_figure


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    parser.add_argument("figures", nargs="+", help="figure_01 … figure_06, 1 … 6, o 'all'")
    parser.add_argument("--method", choices=(*METHODS, "both"), default="tmm", help="método (defecto: tmm)")
    parser.add_argument("--workers", type=int, default=None, help="procesos del PWE (defecto: PWE_WORKERS o todos)")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    figures = list(FIGURES) if "all" in args.figures else [normalize_figure_id(f) for f in args.figures]
    methods = METHODS if args.method == "both" else (args.method,)
    for method in methods:
        for figure_id in figures:
            print("=" * 60)
            print(f"{figure_id} con {method.upper()}")
            run_figure(figure_id, method, workers=args.workers)
