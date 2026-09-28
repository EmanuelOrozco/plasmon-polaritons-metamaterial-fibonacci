#!/usr/bin/env python3
"""Convergencia del ancho de banda de la Fig. 4 respecto a la malla en frecuencia."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from _common import load_figure_config
from fibonacci_photonics.tmm.dispersion import scan_dispersion
from fibonacci_photonics.core.electromagnetics import Polarization
from fibonacci_photonics.analysis.plasmon_modes import detect_plasmon_modes, merge_touching_intervals


def main() -> None:
    spec, cfg, _ = load_figure_config("figure_04.yaml")
    out = ROOT / "results" / "tmm" / "convergence"
    out.mkdir(parents=True, exist_ok=True)
    theta = np.pi / 3
    window = (0.94, 0.999)
    rows = []
    for npts in (500, 1000, 2000, 4000, 8000, 16000):
        nu = np.linspace(*window, npts)
        scan = scan_dispersion(spec, 3, theta, nu, Polarization.TE)
        intervals = merge_touching_intervals(
            list(detect_plasmon_modes(scan, nu_min_ghz=window[0], nu_max_ghz=window[1], min_points=3).intervals),
            gap_ghz=5e-5,
        )
        width = intervals[0].bandwidth_ghz if intervals else float("nan")
        rows.append({"n_points": npts, "n_modes": len(intervals), "width_ghz": width})
        print(rows[-1])
    (out / "convergence_figure04.json").write_text(json.dumps(rows, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
