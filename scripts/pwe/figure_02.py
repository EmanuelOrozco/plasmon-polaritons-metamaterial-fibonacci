#!/usr/bin/env python3
"""Figura 2 con el PWE: polaritones cerca de ν_m = 3 GHz para m = 3..6, θ = π/3."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from _shared import edge_record, edge_text, scan_arrays, timed_pwe_solution
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
from fibonacci_photonics.core.fibonacci import n_layers_b
from fibonacci_photonics.plotting import apply_prb_style, format_dispersion_axes, plot_dispersion_branches, save_figure


def main() -> None:
    spec, physics, numerics, cfg_path = load_pwe_config("figure_02")
    dirs = result_dirs("pwe", "figure_02")
    polarization = Polarization[physics["polarization"]]
    nu = np.linspace(physics["frequency_min_ghz"], physics["frequency_max_ghz"], numerics["frequency_points"])
    theta = float(physics["thetas_rad"][0])
    window = tuple(physics["plasmon_window_ghz"])

    apply_prb_style()
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 6.0), sharex=True, sharey=True)
    reports = []
    for ax, m, panel in zip(axes.ravel(), physics["fibonacci_orders"], "abcd", strict=True):
        scan, edges, seconds, edge_seconds = timed_pwe_solution(spec, m, theta, nu, polarization, numerics)
        save_scan_npz(dirs["data"] / f"m{m}.npz", **scan_arrays(scan, seconds, edges, edge_seconds))
        plot_dispersion_branches(ax, scan, color="#2c4d8c", linestyle="-", linewidth=0.9,
                                 edge_nu_ghz=edges.nu_ghz, edge_k_lm_over_pi=edges.k_lm_over_pi)
        ax.axhline(spec.nu_m_ghz(), color="#c23b22", linestyle="--", linewidth=0.8)
        format_dispersion_axes(ax, nu_min=2.0, nu_max=4.0, panel=panel, m_label=str(m))
        report = detect_plasmon_modes(scan, nu_min_ghz=window[0], nu_max_ghz=window[1], min_points=4)
        intervals = merge_touching_intervals(list(report.intervals), gap_ghz=3e-3)
        reports.append(
            {
                "m": m,
                "expected": n_layers_b(m),
                "detected": len(intervals),
                "seconds": seconds,
                **edge_record(edges, edge_seconds),
                "intervals": [
                    {"nu_min": i.nu_min_ghz, "nu_max": i.nu_max_ghz, "width": i.bandwidth_ghz}
                    for i in intervals
                ],
            }
        )
        print(f"  m={m}: {len(intervals)} subbandas (esperadas {n_layers_b(m)}), {seconds:.1f} s, "
              f"{edge_text(edges, edge_seconds)}")

    fig.tight_layout()
    paths = save_figure(fig, dirs["output"], "figure_02", formats=("pdf", "png"))
    extra = {"config": relative(cfg_path), "numerics": numerics, "mode_counts": reports,
             "outputs": {k: relative(p) for k, p in paths.items()}}
    write_run_sidecar(dirs["data"] / "metadata.json", spec, extra)
    dump_json(dirs["data"] / "summary.json", extra)
    write_figure_readme(
        dirs["root"] / "README.md",
        "Figura 2 (PWE)",
        "Dispersión TE cerca de ν_m = 3 GHz calculada con ondas planas; F(m-2) subbandas. "
        "Los bordes de banda (k = 0 y k = ±1) se refinan con Brent sobre R_PWE.\n\n"
        "Script: `python scripts/pwe/figure_02.py`",
    )


if __name__ == "__main__":
    main()
