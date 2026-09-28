#!/usr/bin/env python3
"""Figura 5 con el PWE: subbandas plasmon-polaritón frente al orden m (m ≤ 6)."""

from __future__ import annotations

import matplotlib.pyplot as plt

from _shared import harmonics_for, plasmon_bands_task
from _common import (
    dump_json,
    load_pwe_config,
    parallel_map,
    relative,
    result_dirs,
    write_figure_readme,
    write_run_sidecar,
)
from fibonacci_photonics.core.electromagnetics import Polarization
from fibonacci_photonics.core.fibonacci import n_layers_b
from fibonacci_photonics.plotting import BRANCH_COLORS, apply_prb_style, latex_theta, save_figure


def main() -> None:
    spec, physics, numerics, cfg_path = load_pwe_config("figure_05")
    dirs = result_dirs("pwe", "figure_05")
    polarization = Polarization[physics["polarization"]]
    windows = physics["frequency_windows_ghz"]
    tasks = [
        {
            "spec": spec,
            "m": m,
            "theta": theta,
            "label": label,
            "nu_min": windows[label][0],
            "grid_points": numerics["grid_points"],
            "min_offset": numerics["nu_m_min_offset_ghz"],
            "harmonics": harmonics_for(numerics, m),
            "rule": numerics["rule"],
            "polarization": polarization,
        }
        for label, theta in zip(physics["theta_labels"], physics["thetas_rad"], strict=True)
        for m in numerics["fibonacci_orders"]
    ]
    results = parallel_map(plasmon_bands_task, tasks)

    apply_prb_style()
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 4.4))
    intervals: dict[str, list] = {}
    for ax, label, panel in zip(axes, physics["theta_labels"], "ab", strict=True):
        window = windows[label]
        stored = []
        for task, res in zip(tasks, results, strict=True):
            if task["label"] != label:
                continue
            m = res["m"]
            bands = [b for b in res["bands"] if b["max"] >= window[0]]
            stored.append({"m": m, "n": len(bands), "expected": n_layers_b(m), "bands": bands,
                           "harmonics_per_layer": res["harmonics_per_layer"], "n_max": res["n_max"],
                           "rejected_bands": res["rejected_bands"],
                           "evaluations": res["evaluations"], "seconds": res["seconds"]})
            for idx, band in enumerate(bands):
                ax.plot([m, m], [band["min"], band["max"]], color=BRANCH_COLORS[idx % len(BRANCH_COLORS)],
                        linewidth=4.0, solid_capstyle="butt")
            print(f"  θ={label} m={m}: {len(bands)} subbandas (esperadas {n_layers_b(m)}), "
                  f"{len(res['rejected_bands'])} espurias descartadas, "
                  f"{res['evaluations']} evaluaciones, {res['seconds']:.0f} s")
        intervals[label] = stored
        ax.set_xlim(1.5, 8.5)
        ax.set_ylim(window[0], window[1])
        ax.set_xticks(range(2, 9))
        ax.set_xlabel("Fibonacci order")
        ax.set_ylabel("bandwidth (GHz)")
        ax.text(0.04, 0.95, rf"({panel})", transform=ax.transAxes, va="top")
        ax.text(0.96, 0.95, rf"$\theta = {latex_theta(label)}$", transform=ax.transAxes, va="top", ha="right")
        ax.text(0.96, 0.05, "PWE: m ≤ 6", transform=ax.transAxes, ha="right", fontsize=7, color="#555555")
    fig.tight_layout()
    paths = save_figure(fig, dirs["output"], "figure_05", formats=("pdf", "png"))
    extra = {"config": relative(cfg_path), "numerics": numerics, "intervals": intervals,
             "outputs": {k: relative(p) for k, p in paths.items()}}
    write_run_sidecar(dirs["data"] / "metadata.json", spec, extra)
    dump_json(dirs["data"] / "summary.json", extra)
    write_figure_readme(
        dirs["root"] / "README.md",
        "Figura 5 (PWE)",
        "Subbandas plasmon-polaritón frente al orden de Fibonacci, m = 2..6, con ondas planas.\n"
        "Bordes refinados con Brent sobre |R_PWE| − 1; cada banda se valida con 1.5× ondas planas.\n"
        "m = 7, 8 se omiten por costo O(N³).\n\n"
        "Script: `python scripts/pwe/figure_05.py`",
    )


if __name__ == "__main__":
    main()
