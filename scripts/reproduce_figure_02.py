#!/usr/bin/env python3
"""Figura 2: polaritones cerca de ν_m = 3 GHz para m = 3..6, θ = π/3."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from _common import dump_json, figure_dirs, load_figure_config, save_scan_npz, write_figure_readme, write_run_sidecar
from fibonacci_tmm.dispersion import scan_dispersion
from fibonacci_tmm.electromagnetics import Polarization
from fibonacci_tmm.fibonacci import n_layers_b
from fibonacci_tmm.plasmon_modes import detect_plasmon_modes, merge_touching_intervals
from fibonacci_tmm.plotting import apply_prb_style, format_dispersion_axes, plot_dispersion_branches, save_figure


def main() -> None:
    spec, cfg, cfg_path = load_figure_config("figure_02.yaml")
    dirs = figure_dirs("figure_02")
    polarization = Polarization[cfg["polarization"]]
    nu = np.linspace(cfg["frequency_min_ghz"], cfg["frequency_max_ghz"], cfg["frequency_points"])
    theta = float(cfg["thetas_rad"][0])
    window = tuple(cfg["plasmon_window_ghz"])
    orders = list(cfg["fibonacci_orders"])

    apply_prb_style()
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 6.0), sharex=True, sharey=True)
    reports = []
    panels = ["a", "b", "c", "d"]
    for ax, m, panel in zip(axes.ravel(), orders, panels, strict=True):
        scan = scan_dispersion(spec, m, theta, nu, polarization)
        save_scan_npz(
            dirs["data"] / f"m{m}.npz",
            nu_ghz=scan.nu_ghz,
            abs_r=scan.abs_r,
            allowed=scan.allowed.astype(np.uint8),
            k_lm_over_pi=scan.k_lm_over_pi,
        )
        plot_dispersion_branches(ax, scan, color="#2c4d8c", linestyle="-", linewidth=0.9)
        ax.axhline(spec.nu_m_ghz(), color="#c23b22", linestyle="--", linewidth=0.8)
        format_dispersion_axes(ax, nu_min=2.0, nu_max=4.0, panel=panel, m_label=str(m))
        report = detect_plasmon_modes(scan, nu_min_ghz=window[0], nu_max_ghz=window[1], min_points=4)
        intervals = merge_touching_intervals(list(report.intervals), gap_ghz=3e-3)
        reports.append(
            {
                "m": m,
                "expected": n_layers_b(m),
                "detected": len(intervals),
                "intervals": [
                    {"nu_min": i.nu_min_ghz, "nu_max": i.nu_max_ghz, "width": i.bandwidth_ghz}
                    for i in intervals
                ],
            }
        )

    fig.tight_layout()
    paths = save_figure(fig, dirs["output"], "figure_02")
    for dest in paths.values():
        (dirs["reproduced"] / dest.name).write_bytes(dest.read_bytes())
    extra = {"config": str(cfg_path), "mode_counts": reports, "outputs": {k: str(p) for k, p in paths.items()}}
    write_run_sidecar(dirs["data"] / "metadata.json", spec, extra)
    dump_json(dirs["data"] / "summary.json", extra)
    write_figure_readme(
        dirs["root"] / "README.md",
        "Figura 2",
        "Dispersión TE cerca de ν_m = 3 GHz. El número de subbandas debe ser F(m-2).\n\n"
        "Script: `python scripts/reproduce_figure_02.py`",
    )
    print("Figura 2:", reports)


if __name__ == "__main__":
    main()
