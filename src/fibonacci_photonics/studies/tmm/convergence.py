"""Convergencia del ancho de banda de la Fig. 4 (TMM) respecto a la malla en frecuencia."""

from __future__ import annotations

import json
from typing import Any

import numpy as np

from fibonacci_photonics.analysis.plasmon_modes import detect_plasmon_modes, merge_touching_intervals
from fibonacci_photonics.config import load_spec
from fibonacci_photonics.io.paths import results_dir
from fibonacci_photonics.physics.electromagnetics import Polarization
from fibonacci_photonics.solvers.tmm import TMMSolver

GRID_POINTS = (500, 1000, 2000, 4000, 8000, 16000)
THETA = np.pi / 3
WINDOW_GHZ = (0.94, 0.999)


def main() -> list[dict[str, Any]]:
    """Ancho de la subbanda de m = 3, θ = π/3, frente al número de puntos de malla."""
    solver = TMMSolver(load_spec("figure_04.yaml"))
    rows = []
    for npts in GRID_POINTS:
        nu = np.linspace(*WINDOW_GHZ, npts)
        scan = solver.scan(3, THETA, nu, Polarization.TE)
        report = detect_plasmon_modes(scan, nu_min_ghz=WINDOW_GHZ[0], nu_max_ghz=WINDOW_GHZ[1], min_points=3)
        intervals = merge_touching_intervals(list(report.intervals), gap_ghz=5e-5)
        width = intervals[0].bandwidth_ghz if intervals else float("nan")
        rows.append({"n_points": npts, "n_modes": len(intervals), "width_ghz": width})
        print(rows[-1])
    path = results_dir() / "tmm" / "convergence" / "convergence_figure04.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(rows, indent=2), encoding="utf-8")
    return rows
