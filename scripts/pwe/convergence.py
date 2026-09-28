#!/usr/bin/env python3
"""Convergencia, costo y límite de validez del PWE, con la TMM como referencia.

1. Bicapa del libro: error relativo de las frecuencias propias frente a N.
2. Superred de Drude: max |R_PWE − R_TMM| en bandas frente a ondas planas por capa,
   reglas de Laurent e inversa, TE y TM, m = 3, 4, 5.
3. Costo: segundos por frecuencia frente a N (m = 3..8), con la TMM como referencia.
4. Validez cerca de ν_m: |R| de ambos métodos frente a ν_m − ν.

Con ``--plot-only`` redibuja las figuras desde data/summary.json sin recalcular.
"""

from __future__ import annotations

import sys
import time

import matplotlib.pyplot as plt
import numpy as np
import yaml
from scipy.optimize import brentq

from _shared import ROOT
from _common import dump_json, load_figure_config, load_json, parallel_map, relative, result_dirs, single_threaded_blas
from fibonacci_photonics.core.bands import omega_from_nu_ghz
from fibonacci_photonics.core.electromagnetics import Polarization
from fibonacci_photonics.core.fibonacci import n_layers_total
from fibonacci_photonics.plotting import apply_prb_style, save_figure
from fibonacci_photonics.pwe.bloch import n_max_for
from fibonacci_photonics.pwe.dispersion import SemitraceEvaluator
from fibonacci_photonics.pwe.eigenfrequency import nondispersive_bands
from fibonacci_photonics.pwe.fourier import bilayer_fourier
from fibonacci_photonics.tmm.transfer_matrix import semitrace_by_recurrence

RULE_STYLE = {"inverse": {"color": "#c23b22", "marker": "o"}, "laurent": {"color": "#2e8b4a", "marker": "s"}}


def _bilayer_r(f, fractions, eps):
    kappa = 2 * np.pi * np.asarray(f, dtype=float)
    q1, q2 = kappa * np.sqrt(eps[0]), kappa * np.sqrt(eps[1])
    ratio = 0.5 * (q1 * eps[1] / (q2 * eps[0]) + q2 * eps[0] / (q1 * eps[1]))
    return np.cos(q1 * fractions[0]) * np.cos(q2 * fractions[1]) - ratio * np.sin(q1 * fractions[0]) * np.sin(
        q2 * fractions[1]
    )


def book_study(cfg):
    fractions = (cfg["thickness_1_um"], cfg["thickness_2_um"])
    total = sum(fractions)
    fractions = (fractions[0] / total, fractions[1] / total)
    eps = (cfg["epsilon_1"], cfg["epsilon_2"])
    ks = np.array(cfg["k_over_pi"]) * np.pi
    n_bands = cfg["n_bands"]
    # Frecuencias exactas: raíces de R(f) = cos k (TMM analítica) refinadas con Brent.
    grid = np.linspace(1e-4, 1.0, 20001)
    exact = []
    for k in ks:
        d = _bilayer_r(grid, fractions, eps) - np.cos(k)
        idx = np.flatnonzero(np.sign(d[1:]) != np.sign(d[:-1]))[:n_bands]
        exact.append([brentq(lambda f: _bilayer_r(f, fractions, eps) - np.cos(k), grid[i], grid[i + 1], xtol=1e-15)
                      for i in idx])
    exact = np.array(exact)
    rows = []
    for rule in ("laurent", "inverse"):
        for n_max in cfg["n_max_values"]:
            cell = bilayer_fourier(*fractions, n_max)
            f = nondispersive_bands(cell, eps, (1, 1), ks, rule=rule, n_bands=n_bands)
            rows.append({"rule": rule, "n_max": n_max, "n_plane_waves": 2 * n_max + 1,
                         "max_rel_error": float(np.max(np.abs(f - exact) / exact))})
    return rows


def drude_case(task):
    spec, m, theta, pol, h, rule, nu = task
    r_tmm = semitrace_by_recurrence(spec, m, omega_from_nu_ghz(nu), theta, pol).real
    evaluator = SemitraceEvaluator(spec, m, theta, pol, n_max=n_max_for(m, h), rule=rule)
    r_pwe = evaluator.many(nu)
    band = np.abs(r_tmm) <= 1
    return {"m": m, "polarization": pol.value, "harmonics_per_layer": h, "rule": rule,
            "n_plane_waves": 2 * n_max_for(m, h) + 1,
            "max_abs_delta_r_in_bands": float(np.nanmax(np.abs(r_pwe - r_tmm)[band])),
            "classification_mismatch": int(np.count_nonzero((np.abs(r_pwe) <= 1) != band))}


def cost_study(cfg):
    spec, _, _ = load_figure_config(cfg["physics"])
    nu = np.linspace(*cfg["window_ghz"], cfg["frequency_points"])
    theta = cfg["theta_rad"]
    rows = []
    for m, h in cfg["cases"]:
        start = time.perf_counter()
        repeats = 200
        for _ in range(repeats):
            r_tmm = semitrace_by_recurrence(spec, m, omega_from_nu_ghz(nu), theta, Polarization.TE).real
        tmm_seconds = (time.perf_counter() - start) / (repeats * nu.size)
        evaluator = SemitraceEvaluator(spec, m, theta, n_max=n_max_for(m, h))
        start = time.perf_counter()
        r_pwe = evaluator.many(nu)
        pwe_seconds = (time.perf_counter() - start) / nu.size
        rows.append({"m": m, "layers": n_layers_total(m), "harmonics_per_layer": h,
                     "n_plane_waves": 2 * n_max_for(m, h) + 1, "pwe_seconds_per_frequency": pwe_seconds,
                     "tmm_seconds_per_frequency": tmm_seconds,
                     "max_abs_delta_r": float(np.max(np.abs(r_pwe - r_tmm))),
                     "max_abs_r_tmm": float(np.max(np.abs(r_tmm)))})
        print(f"  costo m={m} h={h}: N={rows[-1]['n_plane_waves']} PWE {pwe_seconds:.3f} s/frec, "
              f"TMM {tmm_seconds * 1e6:.1f} µs/frec, max|ΔR|={rows[-1]['max_abs_delta_r']:.1e}")
    return rows


def near_pole_study(cfg):
    spec, _, _ = load_figure_config(cfg["physics"])
    m, h, theta = cfg["fibonacci_order"], cfg["harmonics_per_layer"], cfg["theta_rad"]
    offsets = np.logspace(np.log10(cfg["offsets_ghz"][0]), np.log10(cfg["offsets_ghz"][1]), cfg["points"])
    nu = spec.nu_m_ghz() - offsets
    r_tmm = semitrace_by_recurrence(spec, m, omega_from_nu_ghz(nu), theta, Polarization.TE).real
    r_pwe = SemitraceEvaluator(spec, m, theta, n_max=n_max_for(m, h)).many(nu)
    return {"offsets_ghz": offsets.tolist(), "r_tmm": r_tmm.tolist(), "r_pwe": r_pwe.tolist(),
            "m": m, "harmonics_per_layer": h, "theta": theta}


RULE_LABEL = {"laurent": "Laurent (libro)", "inverse": "inversa (Li)"}


def main() -> None:
    single_threaded_blas()
    cfg = yaml.safe_load((ROOT / "configs" / "pwe" / "convergence.yaml").read_text(encoding="utf-8"))
    dirs = result_dirs("pwe", "convergence")
    if "--plot-only" in sys.argv:
        data = load_json(dirs["data"] / "summary.json")
        plot_and_save(data["book"], data["drude"], data["cost"], data["near_pole"], cfg["drude"], dirs)
        return

    print("1) bicapa del libro")
    book = book_study(cfg["book"])

    print("2) superred de Drude")
    dcfg = cfg["drude"]
    spec, _, _ = load_figure_config(dcfg["physics"])
    nu = np.linspace(dcfg["frequency_min_ghz"], dcfg["frequency_max_ghz"], dcfg["frequency_points"])
    tasks = [(spec, m, dcfg["theta_rad"], Polarization[p], h, rule, nu)
             for m in dcfg["fibonacci_orders"] for p in dcfg["polarizations"]
             for h in dcfg["harmonics_per_layer"] for rule in ("inverse", "laurent")]
    tasks.sort(key=lambda t: -n_max_for(t[1], t[4]))
    drude = sorted(parallel_map(drude_case, tasks),
                   key=lambda r: (r["m"], r["polarization"], r["rule"], r["harmonics_per_layer"]))

    print("3) costo (un proceso, un hilo)")
    cost = cost_study(cfg["cost"])

    print("4) validez cerca de ν_m")
    pole = near_pole_study(cfg["near_pole"])
    plot_and_save(book, drude, cost, pole, dcfg, dirs)


def plot_and_save(book, drude, cost, pole, dcfg, dirs) -> None:
    apply_prb_style()
    outputs = {}
    fig, ax = plt.subplots(figsize=(4.2, 3.4))
    for rule in ("laurent", "inverse"):
        rows = [r for r in book if r["rule"] == rule]
        ax.loglog([r["n_plane_waves"] for r in rows], [r["max_rel_error"] for r in rows],
                  **RULE_STYLE[rule], linewidth=0.8, markersize=3.5, label=RULE_LABEL[rule])
    n = np.array([11, 321])
    ax.loglog(n, 0.3 * (n / 11.0) ** -1, color="gray", linestyle=":", linewidth=0.8, label=r"$\propto N^{-1}$")
    ax.loglog(n, 3e-3 * (n / 11.0) ** -3, color="gray", linestyle="--", linewidth=0.8, label=r"$\propto N^{-3}$")
    ax.set_xlabel("número de ondas planas N")
    ax.set_ylabel(r"error relativo máx. de $\omega$")
    ax.set_title("bicapa del libro (4 bandas, 3 valores de k)")
    ax.legend(fontsize=7)
    fig.tight_layout()
    outputs["book_convergence"] = save_figure(fig, dirs["output"], "book_convergence", formats=("pdf", "png"))

    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.8), sharey=True)
    colors = {3: "#1a1a1a", 4: "#2c4d8c", 5: "#c44ec0"}
    for ax, pol in zip(axes, dcfg["polarizations"], strict=True):
        for m in dcfg["fibonacci_orders"]:
            for rule, ls in (("inverse", "-"), ("laurent", "--")):
                rows = [r for r in drude if r["m"] == m and r["polarization"] == pol and r["rule"] == rule]
                ax.loglog([r["n_plane_waves"] for r in rows], [max(r["max_abs_delta_r_in_bands"], 1e-16) for r in rows],
                          color=colors[m], linestyle=ls, marker=RULE_STYLE[rule]["marker"], markersize=3,
                          linewidth=0.8, label=f"m={m}, {RULE_LABEL[rule]}")
        ax.set_xlabel("número de ondas planas N")
        ax.set_title(f"Drude, {pol}, " + r"$\theta=\pi/6$")
    axes[0].set_ylabel(r"max $|R_{\rm PWE}-R_{\rm TMM}|$ en bandas")
    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels, fontsize=6.5, ncol=3, loc="lower center", frameon=False)
    fig.tight_layout(rect=(0, 0.13, 1, 1))
    outputs["drude_convergence"] = save_figure(fig, dirs["output"], "drude_convergence", formats=("pdf", "png"))

    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.3))
    for h, marker in ((8, "s"), (16, "o")):
        rows = [r for r in cost if r["harmonics_per_layer"] == h]
        axes[0].loglog([r["n_plane_waves"] for r in rows], [r["pwe_seconds_per_frequency"] for r in rows],
                       marker=marker, linewidth=0.8, label=f"PWE, {h} armónicos/capa")
        for r in rows:
            axes[0].annotate(f"m={r['m']}", (r["n_plane_waves"], r["pwe_seconds_per_frequency"]), fontsize=6,
                             xytext=(3, -8), textcoords="offset points")
        axes[1].semilogy([r["m"] for r in rows], [r["max_abs_delta_r"] for r in rows], marker=marker,
                         linewidth=0.8, label=f"{h} armónicos/capa")
    tmm = sorted({(r["m"], r["tmm_seconds_per_frequency"]) for r in cost})
    axes[0].axhline(np.mean([t for _, t in tmm]), color="#2c4d8c", linestyle="--", linewidth=0.8, label="TMM")
    nn = np.array([100, 1100])
    axes[0].loglog(nn, 0.05 * (nn / 100.0) ** 3, color="gray", linestyle=":", linewidth=0.8, label=r"$\propto N^3$")
    axes[0].set_xlabel("número de ondas planas N")
    axes[0].set_ylabel("segundos por frecuencia")
    axes[0].legend(fontsize=6)
    axes[1].set_xlabel("orden de Fibonacci m")
    axes[1].set_ylabel(r"max $|R_{\rm PWE}-R_{\rm TMM}|$ (ventana Fig. 5)")
    axes[1].legend(fontsize=6)
    fig.tight_layout()
    outputs["cost"] = save_figure(fig, dirs["output"], "cost_and_accuracy", formats=("pdf", "png"))

    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.2))
    off = np.array(pole["offsets_ghz"])
    rt, rp = np.array(pole["r_tmm"]), np.array(pole["r_pwe"])
    axes[0].loglog(off, np.abs(rt), color="#2c4d8c", linewidth=0.9, label="TMM")
    axes[0].loglog(off, np.abs(rp), color="#c23b22", marker="o", markersize=2, linestyle="none", label="PWE")
    axes[0].axhline(1.0, color="gray", linewidth=0.6)
    axes[0].set_xlabel(r"$\nu_m-\nu$ (GHz)")
    axes[0].set_ylabel(r"$|R|$")
    axes[0].legend(fontsize=7)
    axes[1].loglog(np.abs(rt), np.abs(rp - rt) / np.maximum(np.abs(rt), 1), color="#c23b22", marker="o",
                   markersize=2, linestyle="none")
    axes[1].set_xlabel(r"$|R_{\rm TMM}|$")
    axes[1].set_ylabel(r"$|R_{\rm PWE}-R_{\rm TMM}| / \max(1,|R_{\rm TMM}|)$")
    fig.suptitle(rf"m={pole['m']}, $\theta=\pi/12$, {pole['harmonics_per_layer']} armónicos por capa", fontsize=8)
    fig.tight_layout()
    outputs["near_pole"] = save_figure(fig, dirs["output"], "validity_near_nu_m", formats=("pdf", "png"))

    dump_json(dirs["data"] / "summary.json", {
        "config": "configs/pwe/convergence.yaml",
        "book": book, "drude": drude, "cost": cost, "near_pole": pole,
        "outputs": {name: {k: relative(p) for k, p in paths.items()} for name, paths in outputs.items()},
    })
    print("Convergencia PWE escrita en", dirs["output"])


if __name__ == "__main__":
    main()
