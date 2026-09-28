#!/usr/bin/env python3
"""Figura 6 con el PWE: subbandas plasmon-polaritón frente al ángulo (m ≤ 6)."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

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
from fibonacci_photonics.plotting import BRANCH_COLORS, apply_prb_style, save_figure


def main() -> None:
    spec, physics, numerics, cfg_path = load_pwe_config("figure_06")
    dirs = result_dirs("pwe", "figure_06")
    polarization = Polarization[physics["polarization"]]
    thetas = np.linspace(numerics["theta_min_rad"], numerics["theta_max_rad"], numerics["theta_points"])
    orders = list(numerics["fibonacci_orders"])
    tasks = [
        {
            "spec": spec,
            "m": m,
            "theta": float(theta),
            "nu_min": physics["plasmon_window_ghz"][0],
            "grid_points": numerics["grid_points"],
            "min_offset": numerics["nu_m_min_offset_ghz"],
            "harmonics": harmonics_for(numerics, m),
            "rule": numerics["rule"],
            "polarization": polarization,
        }
        # Primero las tareas caras para repartir mejor la carga.
        for m in sorted(orders, reverse=True)
        for theta in thetas
    ]
    results = parallel_map(plasmon_bands_task, tasks)

    apply_prb_style()
    fig, axes = plt.subplots(3, 2, figsize=(7.0, 7.6), sharex=True, sharey=True)
    payload: dict[str, list] = {}
    width = 0.6 * (thetas[1] - thetas[0])
    for ax, m, panel in zip(axes.ravel(), orders, "abcdef", strict=False):
        rows = sorted((r for r in results if r["m"] == m), key=lambda r: r["theta"])
        payload[str(m)] = rows
        expected = n_layers_b(m)
        incomplete = [round(float(np.degrees(r["theta"])), 1) for r in rows if len(r["bands"]) != expected]
        if incomplete:
            print(f"Aviso m={m}: sin F_(m-2)={expected} subbandas en θ(°) = {incomplete}")
        for row in rows:
            for idx, band in enumerate(sorted(row["bands"], key=lambda b: b["min"])):
                ax.bar(row["theta"], band["max"] - band["min"], bottom=band["min"], width=width,
                       color=BRANCH_COLORS[idx % len(BRANCH_COLORS)], linewidth=0)
        ax.set_xlim(0.0, float(physics["theta_max_rad"]) + width)
        ax.set_ylim(0.95, 1.002)
        ax.set_xticks([0.0, np.pi / 12, np.pi / 6, np.pi / 3])
        ax.set_xticklabels([r"$0$", r"$\pi/12$", r"$\pi/6$", r"$\pi/3$"])
        ax.text(0.04, 0.92, rf"({panel})", transform=ax.transAxes, va="top")
        ax.text(0.96, 0.92, rf"$m={m}$", transform=ax.transAxes, va="top", ha="right")
    for ax in axes.ravel()[len(orders):]:
        ax.axis("off")
        ax.text(0.5, 0.5, "m = 7: omitido en PWE\n(costo O(N³), N ∝ F$_m$)", transform=ax.transAxes,
                ha="center", va="center", fontsize=8, color="#555555")
    for ax in axes[:, 0]:
        ax.set_ylabel("frequency (GHz)")
    for ax in axes[-1, :]:
        ax.set_xlabel(r"$\theta$")
    fig.tight_layout()
    paths = save_figure(fig, dirs["output"], "figure_06", formats=("pdf", "png"))
    extra = {"config": relative(cfg_path), "numerics": numerics, "bands": payload,
             "outputs": {k: relative(p) for k, p in paths.items()}}
    write_run_sidecar(dirs["data"] / "metadata.json", spec, extra)
    dump_json(dirs["data"] / "summary.json", extra)
    write_figure_readme(
        dirs["root"] / "README.md",
        "Figura 6 (PWE)",
        "Subbandas plasmon-polaritón frente al ángulo con ondas planas (m = 2..6, θ = 5°..60°).\n"
        "Cada barra va de ν_min a ν_max de una subbanda, con bordes refinados por Brent.\n\n"
        "Script: `python scripts/pwe/figure_06.py`",
    )


if __name__ == "__main__":
    main()
