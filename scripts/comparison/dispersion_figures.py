#!/usr/bin/env python3
"""Comparación TMM vs PWE de las Figuras 1–4 (barridos de dispersión).

Para cada barrido guardado por ``scripts/pwe/figure_0N.py`` se recalcula la TMM
en la misma malla y se comparan R(ν), la clasificación banda/gap y k Lm/π.
"""

from __future__ import annotations

import time

import matplotlib.pyplot as plt
import numpy as np

from _shared import PWE_RESULTS, TMM_RESULTS
from _common import dump_json, load_figure_config, load_json, relative, result_dirs
from fibonacci_photonics.analysis.comparison import compare_scans
from fibonacci_photonics.core.bands import scan_from_semitrace
from fibonacci_photonics.core.electromagnetics import Polarization
from fibonacci_photonics.plotting import PWE_STYLE, TMM_STYLE, apply_prb_style, format_dispersion_axes, plot_dispersion_branches, save_figure
from fibonacci_photonics.tmm.dispersion import scan_dispersion

THETA_TEX = {0.0: "0", round(np.pi / 12, 6): r"\pi/12", round(np.pi / 6, 6): r"\pi/6", round(np.pi / 3, 6): r"\pi/3"}
M_COLORS = {3: "#c23b22", 4: "#2c4d8c", 5: "#2e8b4a", 6: "#c44ec0"}

FIGURES = {
    "figure_01": {"panel_key": "theta", "ylim": lambda cfg, scan: (0.0, 5.0)},
    "figure_02": {"panel_key": "m", "ylim": lambda cfg, scan: (2.0, 4.0)},
    "figure_03": {"panel_key": "theta", "ylim": lambda cfg, scan: (0.0, 5.0)},
    "figure_04": {"panel_key": "m_theta", "ylim": lambda cfg, scan: (float(scan.nu_ghz[0]), float(scan.nu_ghz[-1]))},
}


def theta_tex(theta: float) -> str:
    return THETA_TEX.get(round(theta, 6), f"{theta:.3f}")


def tmm_timed(spec, m, theta, nu, polarization):
    repeats = max(1, int(2e5 // nu.size))
    start = time.perf_counter()
    for _ in range(repeats):
        scan = scan_dispersion(spec, m, theta, nu, polarization)
    return scan, (time.perf_counter() - start) / repeats


def pwe_markers(ax, scan, every: int) -> None:
    idx = np.flatnonzero(scan.allowed)[::every]
    ax.plot(scan.k_lm_over_pi[idx], scan.nu_ghz[idx], **PWE_STYLE)
    ax.plot(-scan.k_lm_over_pi[idx], scan.nu_ghz[idx], **PWE_STYLE)


def compare_figure(figure_id: str) -> dict:
    spec, physics, _ = load_figure_config(f"{figure_id}.yaml")
    polarization = Polarization[physics["polarization"]]
    dirs = result_dirs("comparison", figure_id)
    options = FIGURES[figure_id]
    files = sorted((PWE_RESULTS / figure_id / "data").glob("*.npz"))
    if not files:
        raise SystemExit(f"Faltan los datos PWE de {figure_id}: ejecute scripts/pwe/{figure_id}.py")

    pairs = []
    for path in files:
        data = np.load(path)
        m, theta = int(data["m"]), float(data["theta"])
        nu = data["nu_ghz"]
        pwe = scan_from_semitrace(spec, m, theta, nu, data["r_real"], polarization, "pwe-inverse")
        tmm, tmm_seconds = tmm_timed(spec, m, theta, nu, polarization)
        agreement = compare_scans(tmm, pwe)
        pairs.append({"m": m, "theta": theta, "tmm": tmm, "pwe": pwe, "agreement": agreement,
                      "tmm_seconds": tmm_seconds, "pwe_seconds": float(data["seconds"])})

    key = options["panel_key"]
    if key == "theta":
        panels = sorted({p["theta"] for p in pairs})
        members = [[p for p in pairs if p["theta"] == t] for t in panels]
    elif key == "m":
        panels = sorted({p["m"] for p in pairs})
        members = [[p for p in pairs if p["m"] == m] for m in panels]
    else:
        panels = sorted({(p["m"], p["theta"]) for p in pairs})
        members = [[p for p in pairs if (p["m"], p["theta"]) == k] for k in panels]

    apply_prb_style()
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 6.2))
    for ax, group, letter in zip(axes.ravel(), members, "abcd", strict=False):
        for p in sorted(group, key=lambda q: q["m"]):
            color = M_COLORS.get(p["m"], "black")
            plot_dispersion_branches(ax, p["tmm"], **{**TMM_STYLE, "color": color})
            pwe_markers(ax, p["pwe"], every=max(1, p["pwe"].nu_ghz.size // 250))
        first = group[0]
        ylim = options["ylim"](physics, first["tmm"])
        m_label = str(first["m"]) if key != "theta" else None
        format_dispersion_axes(ax, nu_min=ylim[0], nu_max=ylim[1], panel=letter,
                               theta_label=theta_tex(first["theta"]), m_label=m_label)
    handles = [plt.Line2D([], [], **TMM_STYLE, label="TMM"), plt.Line2D([], [], **PWE_STYLE, label="PWE")]
    fig.legend(handles=handles, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.01))
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    overlay = save_figure(fig, dirs["output"], f"{figure_id}_overlay", formats=("pdf", "png"))

    fig, axes = plt.subplots(2, 2, figsize=(7.0, 5.4), sharey=True)
    for ax, group, letter in zip(axes.ravel(), members, "abcd", strict=False):
        for p in sorted(group, key=lambda q: q["m"]):
            r_t, r_p = p["tmm"].r.real, p["pwe"].r.real
            rel = np.abs(r_p - r_t) / np.maximum(1.0, np.abs(r_t))
            ax.semilogy(p["tmm"].nu_ghz, np.maximum(rel, 1e-17), color=M_COLORS.get(p["m"], "black"),
                        linewidth=0.6, label=f"m={p['m']}")
            gaps = ~p["tmm"].allowed
            ax.fill_between(p["tmm"].nu_ghz, 1e-17, 1, where=gaps, color="#eeeeee", step="mid", linewidth=0)
        ax.set_ylim(1e-14, 1)
        ax.set_xlabel(r"$\nu$ (GHz)")
        ax.text(0.03, 0.93, f"({letter})", transform=ax.transAxes, va="top")
        ax.text(0.97, 0.93, rf"$\theta={theta_tex(group[0]['theta'])}$", transform=ax.transAxes, va="top", ha="right")
        ax.legend(fontsize=6, loc="lower right")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$|R_{\rm PWE}-R_{\rm TMM}|/\max(1,|R_{\rm TMM}|)$")
    fig.tight_layout()
    errors = save_figure(fig, dirs["output"], f"{figure_id}_error", formats=("pdf", "png"))

    tmm_summary = load_json(TMM_RESULTS / figure_id / "data" / "summary.json")
    pwe_summary = load_json(PWE_RESULTS / figure_id / "data" / "summary.json")
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


def main() -> None:
    for figure_id in FIGURES:
        compare_figure(figure_id)


if __name__ == "__main__":
    main()
