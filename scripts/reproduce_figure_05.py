#!/usr/bin/env python3
"""Figura 5: evolución de las subbandas plasmon-polaritón con m."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from _common import dump_json, figure_dirs, load_figure_config, write_figure_readme, write_run_sidecar
from fibonacci_tmm.dispersion import allowed_intervals, scan_dispersion
from fibonacci_tmm.electromagnetics import Polarization
from fibonacci_tmm.plasmon_modes import intervals_in_window, merge_touching_intervals
from fibonacci_tmm.plotting import BRANCH_COLORS, apply_prb_style, latex_theta, save_figure


def main() -> None:
    spec, cfg, cfg_path = load_figure_config("figure_05.yaml")
    dirs = figure_dirs("figure_05")
    polarization = Polarization[cfg["polarization"]]
    orders = list(cfg["fibonacci_orders"])
    thetas = list(cfg["thetas_rad"])
    labels = list(cfg["theta_labels"])
    windows = cfg["frequency_windows_ghz"]

    apply_prb_style()
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 4.4))
    all_intervals = {}
    panel_ids = ["a", "b"]
    for ax, theta, label, panel in zip(axes, thetas, labels, panel_ids, strict=True):
        window = windows[label]
        nu = np.linspace(window[0], min(window[1], spec.nu_m_ghz() - 1e-4), cfg["frequency_points"])
        stored = []
        for m in orders:
            scan = scan_dispersion(spec, m, theta, nu, polarization)
            intervals = merge_touching_intervals(
                intervals_in_window(allowed_intervals(scan), window[0], window[1], min_points=3),
                gap_ghz=2e-5,
            )
            stored.append({"m": m, "n": len(intervals),
                           "bands": [{"min": i.nu_min_ghz, "max": i.nu_max_ghz} for i in intervals]})
            for idx, interval in enumerate(intervals):
                ax.plot(
                    [m, m],
                    [interval.nu_min_ghz, interval.nu_max_ghz],
                    color=BRANCH_COLORS[idx % len(BRANCH_COLORS)],
                    linewidth=4.0,
                    solid_capstyle="butt",
                )
        all_intervals[label] = stored
        ax.set_xlim(1.5, 8.5)
        ax.set_ylim(window[0], window[1])
        ax.set_xlabel("Fibonacci order")
        ax.set_ylabel("bandwidth (GHz)")
        ax.set_xticks(orders)
        ax.text(0.04, 0.95, rf"({panel})", transform=ax.transAxes, va="top")
        ax.text(
            0.96,
            0.95,
            rf"$\theta = {latex_theta(label)}$",
            transform=ax.transAxes,
            va="top",
            ha="right",
        )
    fig.tight_layout()
    paths = save_figure(fig, dirs["output"], "figure_05")
    for dest in paths.values():
        (dirs["reproduced"] / dest.name).write_bytes(dest.read_bytes())
    extra = {"config": str(cfg_path), "intervals": all_intervals, "outputs": {k: str(p) for k, p in paths.items()}}
    write_run_sidecar(dirs["data"] / "metadata.json", spec, extra)
    dump_json(dirs["data"] / "summary.json", extra)
    write_figure_readme(
        dirs["root"] / "README.md",
        "Figura 5",
        "Estructura de bandas plasmon-polaritón frente al orden de Fibonacci.\n\n"
        "Script: `python scripts/reproduce_figure_05.py`",
    )
    print("Figura 5 escrita")


if __name__ == "__main__":
    main()
