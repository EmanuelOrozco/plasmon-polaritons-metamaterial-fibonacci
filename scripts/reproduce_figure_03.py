#!/usr/bin/env python3
"""Figura 3: dispersión TE con ωe/2π = 3 GHz, ωm/2π = 1 GHz."""

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
from fibonacci_tmm.plotting import (
    S3_STYLE,
    S4_STYLE,
    apply_prb_style,
    format_dispersion_axes,
    plot_dispersion_branches,
    save_figure,
)


def main() -> None:
    spec, cfg, cfg_path = load_figure_config("figure_03.yaml")
    dirs = figure_dirs("figure_03")
    polarization = Polarization[cfg["polarization"]]
    nu = np.linspace(cfg["frequency_min_ghz"], cfg["frequency_max_ghz"], cfg["frequency_points"])
    orders = list(cfg["fibonacci_orders"])
    thetas = list(cfg["thetas_rad"])
    labels = list(cfg["theta_labels"])
    styles = {3: S3_STYLE, 4: S4_STYLE}
    panels = ["a", "b", "c", "d"]

    apply_prb_style()
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 6.2), sharex=True, sharey=True)
    for ax, theta, label, panel in zip(axes.ravel(), thetas, labels, panels, strict=True):
        for m in orders:
            scan = scan_dispersion(spec, m, theta, nu, polarization)
            save_scan_npz(
                dirs["data"] / f"m{m}_theta_{label.replace('/', '_')}.npz",
                nu_ghz=scan.nu_ghz,
                abs_r=scan.abs_r,
                allowed=scan.allowed.astype(np.uint8),
                k_lm_over_pi=scan.k_lm_over_pi,
            )
            plot_dispersion_branches(ax, scan, **styles[m])
        format_dispersion_axes(ax, nu_min=0.0, nu_max=5.0, panel=panel, theta_label=label)

    fig.tight_layout()
    paths = save_figure(fig, dirs["output"], "figure_03")
    for dest in paths.values():
        (dirs["reproduced"] / dest.name).write_bytes(dest.read_bytes())
    extra = {"config": str(cfg_path), "outputs": {k: str(p) for k, p in paths.items()}}
    write_run_sidecar(dirs["data"] / "metadata.json", spec, extra)
    dump_json(dirs["data"] / "summary.json", extra)
    write_figure_readme(
        dirs["root"] / "README.md",
        "Figura 3",
        "Como Fig. 1 con ν_m = 1 GHz dentro del gap ⟨n⟩=0.\n\nScript: `python scripts/reproduce_figure_03.py`",
    )
    print("Figura 3 escrita en", dirs["output"])


if __name__ == "__main__":
    main()
