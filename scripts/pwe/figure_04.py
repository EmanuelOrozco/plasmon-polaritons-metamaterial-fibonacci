#!/usr/bin/env python3
"""Figura 4 con el PWE: zoom de las subbandas plasmon-polaritón, m = 3 y 4."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from _shared import scan_arrays, timed_pwe_scan
from _common import (
    dump_json,
    load_pwe_config,
    relative,
    result_dirs,
    save_scan_npz,
    write_figure_readme,
    write_run_sidecar,
)
from fibonacci_photonics.analysis.plasmon_modes import detect_plasmon_modes, merge_touching_intervals
from fibonacci_photonics.core.electromagnetics import Polarization
from fibonacci_photonics.plotting import apply_prb_style, format_dispersion_axes, plot_dispersion_branches, save_figure

PANELS = [
    (3, "pi/12", np.pi / 12, "#c23b22", "a"),
    (3, "pi/3", np.pi / 3, "#c23b22", "b"),
    (4, "pi/12", np.pi / 12, "#2c4d8c", "c"),
    (4, "pi/3", np.pi / 3, "#2c4d8c", "d"),
]
THETA_TEX = {"pi/12": r"\pi/12", "pi/3": r"\pi/3"}


def main() -> None:
    spec, physics, numerics, cfg_path = load_pwe_config("figure_04")
    dirs = result_dirs("pwe", "figure_04")
    polarization = Polarization[physics["polarization"]]
    windows = physics["frequency_windows_ghz"]

    apply_prb_style()
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 6.0))
    summaries = []
    for ax, (m, label, theta, color, panel) in zip(axes.ravel(), PANELS, strict=True):
        window = windows[label]
        nu = np.linspace(window[0], window[1], numerics["frequency_points"])
        scan, seconds = timed_pwe_scan(spec, m, theta, nu, polarization, numerics)
        save_scan_npz(dirs["data"] / f"m{m}_{label.replace('/', '_')}.npz", **scan_arrays(scan, seconds))
        plot_dispersion_branches(ax, scan, color=color, linestyle="-", linewidth=1.0)
        format_dispersion_axes(ax, nu_min=window[0], nu_max=window[1], panel=panel,
                               theta_label=THETA_TEX[label], m_label=str(m))
        report = detect_plasmon_modes(scan, nu_min_ghz=window[0], nu_max_ghz=window[1], min_points=4)
        intervals = merge_touching_intervals(list(report.intervals), gap_ghz=2e-5)
        summaries.append(
            {
                "m": m,
                "theta": label,
                "detected": len(intervals),
                "seconds": seconds,
                "intervals": [
                    {"nu_min": i.nu_min_ghz, "nu_max": i.nu_max_ghz, "width": i.bandwidth_ghz}
                    for i in intervals
                ],
            }
        )
        print(f"  m={m} θ={label}: {len(intervals)} subbandas, {seconds:.1f} s")

    fig.tight_layout()
    paths = save_figure(fig, dirs["output"], "figure_04", formats=("pdf", "png"))
    extra = {"config": relative(cfg_path), "numerics": numerics, "modes": summaries,
             "outputs": {k: relative(p) for k, p in paths.items()}}
    write_run_sidecar(dirs["data"] / "metadata.json", spec, extra)
    dump_json(dirs["data"] / "summary.json", extra)
    write_figure_readme(
        dirs["root"] / "README.md",
        "Figura 4 (PWE)",
        "Zoom de modos plasmon-polaritón con ondas planas. m=3: 1 subbanda; m=4: 2.\n\n"
        "Script: `python scripts/pwe/figure_04.py`",
    )


if __name__ == "__main__":
    main()
