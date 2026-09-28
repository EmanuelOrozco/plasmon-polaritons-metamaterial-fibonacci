#!/usr/bin/env python3
"""Comparación TMM vs PWE de las Figuras 5 y 6 (bordes de subbandas refinados).

Se usa el mismo localizador de bandas (``analysis.band_edges``) y la misma malla
para ambos métodos; solo cambia quién evalúa R(ν). Así la diferencia de bordes
mide únicamente el error del PWE frente a la TMM.
"""

from __future__ import annotations

import time

import matplotlib.pyplot as plt
import numpy as np

from _shared import PWE_RESULTS, TMM_RESULTS
from _common import dump_json, load_figure_config, load_json, load_pwe_config, relative, result_dirs
from fibonacci_photonics.analysis.band_edges import locate_bands, plasmon_grid
from fibonacci_photonics.core.bands import omega_from_nu_ghz
from fibonacci_photonics.core.electromagnetics import Polarization
from fibonacci_photonics.core.fibonacci import n_layers_b
from fibonacci_photonics.plotting import BRANCH_COLORS, apply_prb_style, latex_theta, save_figure
from fibonacci_photonics.tmm.transfer_matrix import semitrace_by_recurrence


def tmm_bands(spec, m, theta, grid, polarization):
    def r_of_nu(nu: float) -> float:
        return float(semitrace_by_recurrence(spec, m, omega_from_nu_ghz(np.array([nu])), theta, polarization).real[0])

    start = time.perf_counter()
    r_grid = semitrace_by_recurrence(spec, m, omega_from_nu_ghz(grid), theta, polarization).real
    bands = locate_bands(r_of_nu, grid, r_grid)
    return [{"min": b.nu_min_ghz, "max": b.nu_max_ghz} for b in bands], time.perf_counter() - start


def edge_errors(tmm: list[dict], pwe: list[dict]) -> list[dict]:
    rows = []
    for t, p in zip(sorted(tmm, key=lambda b: b["min"]), sorted(pwe, key=lambda b: b["min"]), strict=False):
        width = t["max"] - t["min"]
        rows.append({
            "tmm": t, "pwe": p,
            "delta_min_ghz": p["min"] - t["min"],
            "delta_max_ghz": p["max"] - t["max"],
            "relative_width_error": ((p["max"] - p["min"]) - width) / width if width > 0 else float("nan"),
        })
    return rows


def compare_case(spec, m, theta, nu_min, numerics, polarization, pwe_bands, pwe_seconds, pwe_evals):
    grid = plasmon_grid(nu_min, spec.nu_m_ghz(), numerics["grid_points"], numerics["nu_m_min_offset_ghz"])
    tmm, tmm_seconds = tmm_bands(spec, m, theta, grid, polarization)
    tmm = [b for b in tmm if b["max"] >= nu_min]
    errors = edge_errors(tmm, pwe_bands)
    max_edge = max((max(abs(e["delta_min_ghz"]), abs(e["delta_max_ghz"])) for e in errors), default=0.0)
    max_width = max((abs(e["relative_width_error"]) for e in errors), default=0.0)
    return {"m": m, "theta": theta, "expected": n_layers_b(m), "n_tmm": len(tmm), "n_pwe": len(pwe_bands),
            "tmm_bands": tmm, "pwe_bands": pwe_bands, "edges": errors,
            "max_abs_edge_error_ghz": max_edge, "max_relative_width_error": max_width,
            "tmm_seconds": tmm_seconds, "pwe_seconds": pwe_seconds, "pwe_evaluations": pwe_evals}


def figure_05() -> dict:
    spec, physics, numerics, _ = load_pwe_config("figure_05")
    polarization = Polarization[physics["polarization"]]
    pwe = load_json(PWE_RESULTS / "figure_05" / "data" / "summary.json")["intervals"]
    tmm_fig = load_json(TMM_RESULTS / "figure_05" / "data" / "summary.json")["intervals"]
    dirs = result_dirs("comparison", "figure_05")
    windows = physics["frequency_windows_ghz"]
    cases = {}
    for label, theta in zip(physics["theta_labels"], physics["thetas_rad"], strict=True):
        cases[label] = [compare_case(spec, row["m"], theta, windows[label][0], numerics, polarization,
                                     row["bands"], row["seconds"], row["evaluations"]) for row in pwe[label]]

    apply_prb_style()
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.4), gridspec_kw={"height_ratios": [1.6, 1]})
    for col, label in enumerate(physics["theta_labels"]):
        ax = axes[0, col]
        window = windows[label]
        for row in tmm_fig[label]:
            for idx, band in enumerate(row["bands"]):
                ax.plot([row["m"] - 0.14] * 2, [band["min"], band["max"]], color=BRANCH_COLORS[idx % 8],
                        linewidth=3.2, solid_capstyle="butt")
        for case in cases[label]:
            for idx, band in enumerate(sorted(case["pwe_bands"], key=lambda b: b["min"])):
                ax.plot([case["m"] + 0.14] * 2, [band["min"], band["max"]], color=BRANCH_COLORS[idx % 8],
                        linewidth=3.2, solid_capstyle="butt", alpha=0.55)
        ax.set_xlim(1.5, 8.5)
        ax.set_ylim(*window)
        ax.set_xticks(range(2, 9))
        ax.set_xlabel("Fibonacci order m")
        ax.set_ylabel("frequency (GHz)")
        ax.text(0.04, 0.95, "izq.: TMM   der.: PWE", transform=ax.transAxes, va="top", fontsize=7)
        ax.text(0.96, 0.95, rf"$\theta={latex_theta(label)}$", transform=ax.transAxes, va="top", ha="right")
        ax = axes[1, col]
        for case in cases[label]:
            errs = [abs(e["relative_width_error"]) for e in case["edges"]]
            ax.semilogy(np.full(len(errs), case["m"]) + np.linspace(-0.15, 0.15, max(1, len(errs)))[: len(errs)],
                        np.maximum(errs, 1e-14), "o", markersize=3, color="#c23b22")
        ax.set_xlim(1.5, 6.5)
        ax.set_xticks(range(2, 7))
        ax.set_xlabel("Fibonacci order m")
        ax.set_ylabel(r"$|\Delta w|/w$ por subbanda")
    fig.tight_layout()
    paths = save_figure(fig, dirs["output"], "figure_05_comparison", formats=("pdf", "png"))
    summary = {"figure": "figure_05", "cases": cases, "outputs": {k: relative(v) for k, v in paths.items()}}
    dump_json(dirs["data"] / "summary.json", summary)
    worst = max(c["max_relative_width_error"] for rows in cases.values() for c in rows)
    counts_ok = all(c["n_tmm"] == c["n_pwe"] == c["expected"] for rows in cases.values() for c in rows)
    print(f"figure_05: peor error relativo de ancho = {worst:.2e}; conteos coinciden = {counts_ok}")
    return summary


def figure_06() -> dict:
    spec, physics, numerics, _ = load_pwe_config("figure_06")
    polarization = Polarization[physics["polarization"]]
    pwe = load_json(PWE_RESULTS / "figure_06" / "data" / "summary.json")["bands"]
    tmm_fig = load_json(TMM_RESULTS / "figure_06" / "data" / "summary.json")["bandwidths"]
    dirs = result_dirs("comparison", "figure_06")
    nu_min = physics["plasmon_window_ghz"][0]
    cases = {m: [compare_case(spec, int(m), row["theta"], nu_min, numerics, polarization, row["bands"],
                              row["seconds"], row["evaluations"]) for row in rows] for m, rows in pwe.items()}

    apply_prb_style()
    orders = sorted(int(m) for m in pwe)
    fig, axes = plt.subplots(3, 2, figsize=(7.0, 7.6), sharex=True, sharey=True)
    for ax, m, letter in zip(axes.ravel(), orders, "abcdef", strict=False):
        rows = tmm_fig[str(m)]
        thetas = np.array([r["theta"] for r in rows])
        n = n_layers_b(m)
        lo = np.full((n, thetas.size), np.nan)
        hi = np.full((n, thetas.size), np.nan)
        for j, r in enumerate(rows):
            if r["n_modes"] != n:
                continue
            pairs = sorted((c - w / 2, c + w / 2) for c, w in zip(r["centers_ghz"], r["bandwidths_ghz"], strict=True))
            for i, (a, b) in enumerate(pairs):
                lo[i, j], hi[i, j] = a, b
        for i in range(n):
            ok = np.isfinite(lo[i])
            ax.fill_between(thetas[ok], lo[i, ok], hi[i, ok], color=BRANCH_COLORS[i % 8], alpha=0.35, linewidth=0)
        for case in cases[str(m)]:
            for i, band in enumerate(sorted(case["pwe_bands"], key=lambda b: b["min"])):
                ax.plot([case["theta"]] * 2, [band["min"], band["max"]], color=BRANCH_COLORS[i % 8], linewidth=1.6)
                ax.plot([case["theta"]] * 2, [band["min"], band["max"]], "_", color="black", markersize=4)
        ax.set_ylim(0.95, 1.002)
        ax.set_xticks([0.0, np.pi / 12, np.pi / 6, np.pi / 3])
        ax.set_xticklabels([r"$0$", r"$\pi/12$", r"$\pi/6$", r"$\pi/3$"])
        ax.text(0.04, 0.92, f"({letter})", transform=ax.transAxes, va="top")
        ax.text(0.96, 0.92, rf"$m={m}$", transform=ax.transAxes, va="top", ha="right")
    last = axes.ravel()[len(orders)]
    last.axis("off")
    last.text(0.5, 0.5, "Relleno: TMM (41 ángulos)\nBarras: PWE (bordes refinados)", transform=last.transAxes,
              ha="center", va="center", fontsize=8)
    for ax in axes[:, 0]:
        ax.set_ylabel("frequency (GHz)")
    for ax in axes[-1, :]:
        ax.set_xlabel(r"$\theta$")
    fig.tight_layout()
    overlay = save_figure(fig, dirs["output"], "figure_06_overlay", formats=("pdf", "png"))

    fig, ax = plt.subplots(figsize=(4.6, 3.4))
    for idx, m in enumerate(orders):
        rows = sorted(cases[str(m)], key=lambda c: c["theta"])
        ax.semilogy([np.degrees(c["theta"]) for c in rows], [max(c["max_abs_edge_error_ghz"], 1e-16) for c in rows],
                    marker="o", markersize=3, linewidth=0.8, color=BRANCH_COLORS[idx % 8], label=f"m={m}")
    ax.set_xlabel(r"$\theta$ (grados)")
    ax.set_ylabel(r"max $|\nu_{\rm borde}^{\rm PWE}-\nu_{\rm borde}^{\rm TMM}|$ (GHz)")
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    errors = save_figure(fig, dirs["output"], "figure_06_edge_error", formats=("pdf", "png"))
    summary = {"figure": "figure_06", "cases": cases,
               "outputs": {"overlay": {k: relative(v) for k, v in overlay.items()},
                           "edge_error": {k: relative(v) for k, v in errors.items()}}}
    dump_json(dirs["data"] / "summary.json", summary)
    worst = max(c["max_abs_edge_error_ghz"] for rows in cases.values() for c in rows)
    mismatches = [(c["m"], round(np.degrees(c["theta"]), 1), c["n_tmm"], c["n_pwe"])
                  for rows in cases.values() for c in rows if c["n_tmm"] != c["n_pwe"]]
    print(f"figure_06: peor error de borde = {worst:.2e} GHz; conteos distintos: {mismatches or 'ninguno'}")
    return summary


def main() -> None:
    figure_05()
    figure_06()


if __name__ == "__main__":
    main()
