"""Figuras 1 y 3 con el PWE: ν(k) de S3 y S4 en cuatro ángulos."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from _shared import edge_record, edge_text, scan_arrays, timed_pwe_solution  # noqa: F401  (ajusta sys.path)
from _common import (
    dump_json,
    load_pwe_config,
    relative,
    result_dirs,
    save_scan_npz,
    write_figure_readme,
    write_run_sidecar,
)
from fibonacci_photonics.core.electromagnetics import Polarization
from fibonacci_photonics.plotting import (
    S3_STYLE,
    S4_STYLE,
    apply_prb_style,
    format_dispersion_axes,
    plot_dispersion_branches,
    save_figure,
)


def run(figure_id: str, title: str, description: str) -> None:
    spec, physics, numerics, cfg_path = load_pwe_config(figure_id)
    dirs = result_dirs("pwe", figure_id)
    polarization = Polarization[physics["polarization"]]
    nu = np.linspace(physics["frequency_min_ghz"], physics["frequency_max_ghz"], numerics["frequency_points"])
    styles = {3: S3_STYLE, 4: S4_STYLE}

    apply_prb_style()
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 6.2), sharex=True, sharey=True)
    runs = []
    for ax, theta, label, panel in zip(
        axes.ravel(), physics["thetas_rad"], physics["theta_labels"], "abcd", strict=True
    ):
        for m in physics["fibonacci_orders"]:
            scan, edges, seconds, edge_seconds = timed_pwe_solution(spec, m, theta, nu, polarization, numerics)
            save_scan_npz(dirs["data"] / f"m{m}_theta_{label.replace('/', '_')}.npz",
                          **scan_arrays(scan, seconds, edges, edge_seconds))
            plot_dispersion_branches(ax, scan, edge_nu_ghz=edges.nu_ghz,
                                     edge_k_lm_over_pi=edges.k_lm_over_pi, **styles[m])
            runs.append({"m": m, "theta": theta, "theta_label": label, "seconds": seconds,
                         "n_allowed": int(scan.allowed.sum()), **edge_record(edges, edge_seconds)})
            print(f"  m={m} θ={label}: {seconds:.1f} s, {edge_text(edges, edge_seconds)}")
        format_dispersion_axes(ax, nu_min=0.0, nu_max=5.0, panel=panel, theta_label=label)
    fig.tight_layout()
    paths = save_figure(fig, dirs["output"], figure_id, formats=("pdf", "png"))
    extra = {"config": relative(cfg_path), "numerics": numerics, "runs": runs,
             "outputs": {k: relative(p) for k, p in paths.items()}}
    write_run_sidecar(dirs["data"] / "metadata.json", spec, extra)
    dump_json(dirs["data"] / "summary.json", extra)
    write_figure_readme(
        dirs["root"] / "README.md",
        f"{title} (PWE)",
        f"{description}\n\nCalculada con el método de ondas planas, forma k(ω) y regla inversa; "
        "los bordes de banda (k = 0 y k = ±1) se refinan con Brent sobre R_PWE.\n\n"
        f"Script: `python scripts/pwe/{figure_id}.py`",
    )
    print(title, "(PWE) escrita en", dirs["output"])
