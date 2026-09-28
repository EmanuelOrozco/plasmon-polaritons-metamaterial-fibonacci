#!/usr/bin/env python3
"""Figura 4: zoom de las subbandas plasmon-polaritón, m = 3 y 4."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from _common import dump_json, relative, figure_dirs, load_figure_config, save_scan_npz, write_figure_readme, write_run_sidecar
from fibonacci_photonics.tmm.dispersion import scan_dispersion
from fibonacci_photonics.core.electromagnetics import Polarization
from fibonacci_photonics.analysis.plasmon_modes import detect_plasmon_modes, merge_touching_intervals
from fibonacci_photonics.plotting import apply_prb_style, format_dispersion_axes, plot_dispersion_branches, save_figure


def main() -> None:
    spec, cfg, cfg_path = load_figure_config("figure_04.yaml")
    dirs = figure_dirs("figure_04")
    polarization = Polarization[cfg["polarization"]]
    windows = cfg["frequency_windows_ghz"]
    panels = [
        (3, "pi/12", 0.2617993877991494, windows["pi/12"], "#c23b22", "a"),
        (3, "pi/3", 1.0471975511965976, windows["pi/3"], "#c23b22", "b"),
        (4, "pi/12", 0.2617993877991494, windows["pi/12"], "#2c4d8c", "c"),
        (4, "pi/3", 1.0471975511965976, windows["pi/3"], "#2c4d8c", "d"),
    ]
    apply_prb_style()
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 6.0))
    summaries = []
    theta_tex = {"pi/12": r"\pi/12", "pi/3": r"\pi/3"}
    for ax, (m, label, theta, window, color, panel) in zip(axes.ravel(), panels, strict=True):
        nu = np.linspace(window[0], window[1], cfg["frequency_points"])
        scan = scan_dispersion(spec, m, theta, nu, polarization)
        save_scan_npz(
            dirs["data"] / f"m{m}_{label.replace('/', '_')}.npz",
            nu_ghz=scan.nu_ghz,
            abs_r=scan.abs_r,
            allowed=scan.allowed.astype(np.uint8),
            k_lm_over_pi=scan.k_lm_over_pi,
        )
        plot_dispersion_branches(ax, scan, color=color, linestyle="-", linewidth=1.0)
        format_dispersion_axes(
            ax,
            nu_min=window[0],
            nu_max=window[1],
            panel=panel,
            theta_label=theta_tex[label],
            m_label=str(m),
        )
        report = detect_plasmon_modes(scan, nu_min_ghz=window[0], nu_max_ghz=window[1], min_points=4)
        intervals = merge_touching_intervals(list(report.intervals), gap_ghz=2e-5)
        summaries.append(
            {
                "m": m,
                "theta": label,
                "detected": len(intervals),
                "intervals": [
                    {"nu_min": i.nu_min_ghz, "nu_max": i.nu_max_ghz, "width": i.bandwidth_ghz}
                    for i in intervals
                ],
            }
        )
    fig.tight_layout()
    paths = save_figure(fig, dirs["output"], "figure_04")
    extra = {"config": relative(cfg_path), "modes": summaries, "outputs": {k: relative(p) for k, p in paths.items()}}
    write_run_sidecar(dirs["data"] / "metadata.json", spec, extra)
    dump_json(dirs["data"] / "summary.json", extra)
    write_figure_readme(
        dirs["root"] / "README.md",
        "Figura 4",
        "Zoom de modos plasmon-polaritón. m=3: 1 subbanda; m=4: 2.\n\n"
        "Script: `python scripts/tmm/figure_04.py`",
    )
    print("Figura 4:", summaries)


if __name__ == "__main__":
    main()
