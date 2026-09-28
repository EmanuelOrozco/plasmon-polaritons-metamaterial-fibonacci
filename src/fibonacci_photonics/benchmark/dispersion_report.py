"""Comparación TMM vs PWE de las Figuras 1–4 (barridos de dispersión).

Para cada barrido guardado por el PWE se recalcula la TMM en la misma malla y
se comparan R(ν), la clasificación banda/gap y k L_m/π.
"""

from __future__ import annotations

from functools import partial
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

from fibonacci_photonics.benchmark.compare import compare_scans
from fibonacci_photonics.benchmark.timing import timed
from fibonacci_photonics.config import load_figure
from fibonacci_photonics.io.latex import theta_tex
from fibonacci_photonics.io.paths import relative, result_dirs, results_dir
from fibonacci_photonics.io.results import dump_json, load_json
from fibonacci_photonics.solvers.scan import DispersionScan, scan_from_semitrace
from fibonacci_photonics.solvers.tmm import TMMSolver
from fibonacci_photonics.viz.dispersion import format_dispersion_axes, plot_dispersion_branches
from fibonacci_photonics.viz.style import M_COLORS, PWE_FORMATS, PWE_STYLE, TMM_STYLE, apply_prb_style, save_figure

TMM_TIMING_SAMPLES = 2e5
"""Frecuencias totales por medición de tiempo de la TMM (se repite el barrido)."""

FIGURES = {
    "figure_01": "theta",
    "figure_02": "m",
    "figure_03": "theta",
    "figure_04": "m_theta",
}
"""Cómo se agrupan los barridos en paneles en cada figura."""

FIXED_YLIM = {"figure_01": (0.0, 5.0), "figure_02": (2.0, 4.0), "figure_03": (0.0, 5.0)}


def pwe_markers(ax: plt.Axes, scan: DispersionScan, every: int) -> None:
    idx = np.flatnonzero(scan.allowed)[::every]
    ax.plot(scan.k_lm_over_pi[idx], scan.nu_ghz[idx], **PWE_STYLE)
    ax.plot(-scan.k_lm_over_pi[idx], scan.nu_ghz[idx], **PWE_STYLE)


def compare_figure(figure_id: str) -> dict[str, Any]:
    """Superposición, error relativo y resumen de una figura de dispersión."""
    setup = load_figure(figure_id)
    spec, polarization = setup.spec, setup.config.polarization
    tmm_solver = TMMSolver(spec)
    dirs = result_dirs("comparison", figure_id)
    files = sorted((results_dir() / "pwe" / figure_id / "data").glob("*.npz"))
    if not files:
        raise SystemExit(f"Faltan los datos PWE de {figure_id}: ejecute scripts/run_figure.py {figure_id} --method pwe")

    pairs: list[dict[str, Any]] = []
    for path in files:
        with np.load(path) as data:
            m, theta = int(data["m"]), float(data["theta"])
            nu = data["nu_ghz"]
            pwe = scan_from_semitrace(spec, m, theta, nu, data["r_real"], polarization, "pwe-inverse")
            pwe_seconds = float(data["seconds"])
        repeats = max(1, int(TMM_TIMING_SAMPLES // nu.size))
        tmm, tmm_seconds = timed(partial(tmm_solver.scan, m, theta, nu, polarization), repeats)
        pairs.append({"m": m, "theta": theta, "tmm": tmm, "pwe": pwe, "agreement": compare_scans(tmm, pwe),
                      "tmm_seconds": tmm_seconds, "pwe_seconds": pwe_seconds})

    key = FIGURES[figure_id]
    if key == "theta":
        groups = [[p for p in pairs if p["theta"] == t] for t in sorted({p["theta"] for p in pairs})]
    elif key == "m":
        groups = [[p for p in pairs if p["m"] == m] for m in sorted({p["m"] for p in pairs})]
    else:
        keys = sorted({(p["m"], p["theta"]) for p in pairs})
        groups = [[p for p in pairs if (p["m"], p["theta"]) == k] for k in keys]

    apply_prb_style()
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 6.2))
    for ax, group, letter in zip(axes.ravel(), groups, "abcd", strict=False):
        for p in sorted(group, key=lambda q: q["m"]):
            plot_dispersion_branches(ax, p["tmm"], **{**TMM_STYLE, "color": M_COLORS.get(p["m"], "black")})
            pwe_markers(ax, p["pwe"], every=max(1, p["pwe"].nu_ghz.size // 250))
        first = group[0]
        nu_first = first["tmm"].nu_ghz
        ylim = FIXED_YLIM.get(figure_id, (float(nu_first[0]), float(nu_first[-1])))
        m_label = str(first["m"]) if key != "theta" else None
        format_dispersion_axes(ax, nu_min=ylim[0], nu_max=ylim[1], panel=letter,
                               theta_label=theta_tex(first["theta"]), m_label=m_label)
    handles = [plt.Line2D([], [], **TMM_STYLE, label="TMM"), plt.Line2D([], [], **PWE_STYLE, label="PWE")]
    fig.legend(handles=handles, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.01))
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    overlay = save_figure(fig, dirs["output"], f"{figure_id}_overlay", formats=PWE_FORMATS)

    fig, axes = plt.subplots(2, 2, figsize=(7.0, 5.4), sharey=True)
    for ax, group, letter in zip(axes.ravel(), groups, "abcd", strict=False):
        for p in sorted(group, key=lambda q: q["m"]):
            r_t, r_p = p["tmm"].r.real, p["pwe"].r.real
            rel = np.abs(r_p - r_t) / np.maximum(1.0, np.abs(r_t))
            ax.semilogy(p["tmm"].nu_ghz, np.maximum(rel, 1e-17), color=M_COLORS.get(p["m"], "black"),
                        linewidth=0.6, label=f"m={p['m']}")
            ax.fill_between(p["tmm"].nu_ghz, 1e-17, 1, where=~p["tmm"].allowed, color="#eeeeee", step="mid",
                            linewidth=0)
        ax.set_ylim(1e-14, 1)
        ax.set_xlabel(r"$\nu$ (GHz)")
        ax.text(0.03, 0.93, f"({letter})", transform=ax.transAxes, va="top")
        ax.text(0.97, 0.93, rf"$\theta={theta_tex(group[0]['theta'])}$", transform=ax.transAxes, va="top",
                ha="right")
        ax.legend(fontsize=6, loc="lower right")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$|R_{\rm PWE}-R_{\rm TMM}|/\max(1,|R_{\rm TMM}|)$")
    fig.tight_layout()
    errors = save_figure(fig, dirs["output"], f"{figure_id}_error", formats=PWE_FORMATS)

    tmm_summary = load_json(results_dir() / "tmm" / figure_id / "data" / "summary.json")
    pwe_summary = load_json(results_dir() / "pwe" / figure_id / "data" / "summary.json")
    rows = [{"m": p["m"], "theta": p["theta"], "n_points": p["agreement"].n_points,
             "tmm_seconds": p["tmm_seconds"], "pwe_seconds": p["pwe_seconds"],
             **p["agreement"].as_dict()} for p in pairs]
    summary = {
        "figure": figure_id,
        "scans": rows,
        "worst_max_abs_delta_r_in_bands": max(r["max_abs_delta_r_in_bands"] for r in rows),
        "worst_classification_agreement": min(r["classification_agreement"] for r in rows),
        "mode_counts": {"tmm": tmm_summary.get("mode_counts") or tmm_summary.get("modes"),
                        "pwe": pwe_summary.get("mode_counts") or pwe_summary.get("modes")},
        "outputs": {"overlay": {k: relative(v) for k, v in overlay.items()},
                    "error": {k: relative(v) for k, v in errors.items()}},
    }
    dump_json(dirs["data"] / "summary.json", summary)
    print(f"{figure_id}: max|ΔR| en bandas = {summary['worst_max_abs_delta_r_in_bands']:.2e}, "
          f"acuerdo de clasificación mínimo = {summary['worst_classification_agreement']:.5f}")
    return summary


def main() -> dict[str, dict[str, Any]]:
    return {figure_id: compare_figure(figure_id) for figure_id in FIGURES}
