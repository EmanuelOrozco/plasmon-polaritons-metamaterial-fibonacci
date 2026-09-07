#!/usr/bin/env python3
"""Reproduce la Figura 1: dispersión TE para S3 y S4, varios ángulos."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from _common import (  # noqa: E402
    dump_json,
    figure_dirs,
    load_figure_config,
    save_scan_npz,
    write_figure_readme,
    write_run_sidecar,
)
from fibonacci_tmm.dispersion import scan_dispersion, zero_average_index_frequency_ghz
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
    spec, cfg, cfg_path = load_figure_config("figure_01.yaml")
    dirs = figure_dirs("figure_01")
    polarization = Polarization[cfg["polarization"]]
    nu = np.linspace(cfg["frequency_min_ghz"], cfg["frequency_max_ghz"], cfg["frequency_points"])
    orders = list(cfg["fibonacci_orders"])
    thetas = list(cfg["thetas_rad"])
    labels = list(cfg["theta_labels"])
    styles = {3: S3_STYLE, 4: S4_STYLE}
    panels = ["a", "b", "c", "d"]

    apply_prb_style()
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 6.2), sharex=True, sharey=True)
    summaries = []

    for ax, theta, label, panel in zip(axes.ravel(), thetas, labels, panels, strict=True):
        for m in orders:
            scan = scan_dispersion(spec, m, theta, nu, polarization, method="recurrence")
            save_scan_npz(
                dirs["data"] / f"m{m}_theta_{label.replace('/', '_')}.npz",
                nu_ghz=scan.nu_ghz,
                abs_r=scan.abs_r,
                allowed=scan.allowed.astype(np.uint8),
                k_lm_over_pi=scan.k_lm_over_pi,
                r_real=np.real(scan.r),
                r_imag=np.imag(scan.r),
            )
            plot_dispersion_branches(ax, scan, **styles[m])
            summaries.append(
                {
                    "m": m,
                    "theta": theta,
                    "theta_label": label,
                    "n_allowed": int(scan.allowed.sum()),
                    "clipped_near_unit": scan.clipped_near_unit,
                    "max_abs_imag_r": float(np.max(np.abs(np.imag(scan.r)))),
                }
            )
        format_dispersion_axes(ax, nu_min=0.0, nu_max=5.0, panel=panel, theta_label=label)

    fig.tight_layout()
    paths = save_figure(fig, dirs["output"], "figure_01")
    for dest in paths.values():
        target = dirs["reproduced"] / dest.name
        target.write_bytes(dest.read_bytes())

    extra = {
        "config": str(cfg_path),
        "n0_gap_closed_ghz": {
            str(m): zero_average_index_frequency_ghz(spec, m) for m in orders
        },
        "panels": summaries,
        "outputs": {key: str(path) for key, path in paths.items()},
    }
    write_run_sidecar(dirs["data"] / "metadata.json", spec, extra)
    dump_json(dirs["data"] / "summary.json", extra)
    write_figure_readme(
        dirs["root"] / "README.md",
        "Figura 1",
        """
Relación de dispersión TE ν(k) para celdas S3 (discontinua) y S4 (continua).

Parámetros del paper: a = b = 12 mm, εA = μA = 1, ωe/2π = ωm/2π = 3 GHz.
Ángulos: 0, π/12, π/6, π/3.

Script: `python scripts/reproduce_figure_01.py`
""",
    )
    print("Figura 1 escrita en", dirs["output"])


if __name__ == "__main__":
    main()
