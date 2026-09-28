"""Comparación TMM vs PWE de las Figuras 1–4 (barridos de dispersión).

La superposición dibuja las mismas curvas que ``results/tmm`` y ``results/pwe``
(``reproduction.dispersion``), cerradas en sus bordes exactos. Para medir el
acuerdo se recalcula la TMM en la malla de cada barrido PWE y se comparan R(ν),
la clasificación banda/gap, k L_m/π y los bordes |R| = 1 refinados.
"""

from __future__ import annotations

from functools import partial
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

from fibonacci_photonics.analysis.band_closure import EDGE, BandEdges, close_band_edges
from fibonacci_photonics.comparison.compare import compare_scans
from fibonacci_photonics.comparison.timing import timed
from fibonacci_photonics.config import load_figure
from fibonacci_photonics.io.latex import theta_tex
from fibonacci_photonics.io.paths import relative, result_dirs, results_dir
from fibonacci_photonics.io.results import dump_json, load_json
from fibonacci_photonics.reproduction.dispersion import branch_style, load_solved, plot_solved
from fibonacci_photonics.solvers.tmm import TMMSolver
from fibonacci_photonics.viz.dispersion import format_dispersion_axes
from fibonacci_photonics.viz.style import (
    M_COLORS,
    PWE_FORMATS,
    TMM_UNDERLAY,
    apply_prb_style,
    lighten,
    save_figure,
)

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


def edge_deltas(pwe: BandEdges | None, tmm: BandEdges) -> dict[str, Any]:
    """Distancia de cada borde |R| = 1 del PWE al borde TMM más cercano (misma malla)."""
    pwe_nu = pwe.nu_ghz[pwe.kind == EDGE] if pwe is not None else np.array([])
    tmm_nu = tmm.nu_ghz[tmm.kind == EDGE]
    if pwe_nu.size == 0 or tmm_nu.size == 0:
        return {"n_edges_pwe": int(pwe_nu.size), "n_edges_tmm": int(tmm_nu.size), "max_abs_edge_delta_ghz": None}
    deltas = np.min(np.abs(pwe_nu[:, None] - tmm_nu[None, :]), axis=1)
    return {"n_edges_pwe": int(pwe_nu.size), "n_edges_tmm": int(tmm_nu.size),
            "max_abs_edge_delta_ghz": float(deltas.max()), "median_abs_edge_delta_ghz": float(np.median(deltas))}


def figure_scans(method: str, figure_id: str) -> list[Path]:
    files = sorted((results_dir() / method / figure_id / "data").glob("*.npz"))
    if not files:
        raise SystemExit(
            f"Faltan los barridos {method.upper()} de {figure_id}: "
            f"ejecute scripts/run_figure.py {figure_id} --method {method}"
        )
    return files


def compare_figure(figure_id: str) -> dict[str, Any]:
    """Superposición, error relativo y resumen de una figura de dispersión."""
    setup = load_figure(figure_id)
    spec, polarization = setup.spec, setup.config.polarization
    tmm_solver = TMMSolver(spec)
    dirs = result_dirs("comparison", figure_id)
    tmm_files = {path.name: path for path in figure_scans("tmm", figure_id)}

    pairs: list[dict[str, Any]] = []
    for path in figure_scans("pwe", figure_id):
        if path.name not in tmm_files:
            raise SystemExit(f"{relative(path)} no tiene barrido TMM equivalente: vuelva a correr la TMM")
        pwe_solved = load_solved(path, spec, polarization, "pwe-inverse")
        tmm_figure = load_solved(tmm_files[path.name], spec, polarization, "recurrence")
        pwe = pwe_solved.scan
        m, theta, nu = pwe.m, pwe.theta, pwe.nu_ghz
        repeats = max(1, int(TMM_TIMING_SAMPLES // nu.size))
        tmm, tmm_seconds = timed(partial(tmm_solver.scan, m, theta, nu, polarization), repeats)
        tmm_edges = close_band_edges(tmm_solver, tmm)
        pairs.append({"m": m, "theta": theta, "tmm": tmm, "pwe": pwe, "agreement": compare_scans(tmm, pwe),
                      "edges": edge_deltas(pwe_solved.edges, tmm_edges),
                      "tmm_figure": tmm_figure, "pwe_figure": pwe_solved,
                      "tmm_seconds": tmm_seconds, "pwe_seconds": pwe_solved.seconds})

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
            style = branch_style(figure_id, p["m"])
            plot_solved(ax, p["tmm_figure"], **{**style, **TMM_UNDERLAY, "color": lighten(style["color"])})
            plot_solved(ax, p["pwe_figure"], **style)
        first = group[0]
        nu_first = first["tmm"].nu_ghz
        ylim = FIXED_YLIM.get(figure_id, (float(nu_first[0]), float(nu_first[-1])))
        m_label = str(first["m"]) if key != "theta" else None
        format_dispersion_axes(ax, nu_min=ylim[0], nu_max=ylim[1], panel=letter,
                               theta_label=theta_tex(first["theta"]), m_label=m_label)
    handles = [plt.Line2D([], [], color=lighten("#1a1a1a"), **TMM_UNDERLAY, label="TMM (trazo ancho)"),
               plt.Line2D([], [], color="#1a1a1a", linewidth=0.8, label="PWE (trazo fino)")]
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
             **p["agreement"].as_dict(), **p["edges"]} for p in pairs]
    edge_errors = [r["max_abs_edge_delta_ghz"] for r in rows if r["max_abs_edge_delta_ghz"] is not None]
    summary = {
        "figure": figure_id,
        "scans": rows,
        "worst_max_abs_delta_r_in_bands": max(r["max_abs_delta_r_in_bands"] for r in rows),
        "worst_classification_agreement": min(r["classification_agreement"] for r in rows),
        "worst_max_abs_edge_delta_ghz": max(edge_errors, default=None),
        "mode_counts": {"tmm": tmm_summary.get("mode_counts") or tmm_summary.get("modes"),
                        "pwe": pwe_summary.get("mode_counts") or pwe_summary.get("modes")},
        "outputs": {"overlay": {k: relative(v) for k, v in overlay.items()},
                    "error": {k: relative(v) for k, v in errors.items()}},
    }
    dump_json(dirs["data"] / "summary.json", summary)
    edge_text = "sin bordes" if not edge_errors else f"{max(edge_errors):.2e} GHz"
    print(f"{figure_id}: max|ΔR| en bandas = {summary['worst_max_abs_delta_r_in_bands']:.2e}, "
          f"acuerdo de clasificación mínimo = {summary['worst_classification_agreement']:.5f}, "
          f"max|Δν_borde| = {edge_text}")
    return summary


def main() -> dict[str, dict[str, Any]]:
    return {figure_id: compare_figure(figure_id) for figure_id in FIGURES}
