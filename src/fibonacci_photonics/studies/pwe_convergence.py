"""Convergencia, costo y límite de validez del PWE, con la TMM como referencia.

1. Bicapa del libro: error relativo de las frecuencias propias frente a N.
2. Superred de Drude: max |R_PWE − R_TMM| en bandas frente a ondas planas por
   capa, reglas de Laurent e inversa, TE y TM, m = 3, 4, 5.
3. Costo: segundos por frecuencia frente a N (m = 3..8), con la TMM como referencia.
4. Validez cerca de ν_m: |R| de ambos métodos frente a ν_m − ν.

Con ``plot_only=True`` se redibuja desde ``data/summary.json`` sin recalcular.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import brentq

from fibonacci_photonics.config import load_model, load_spec
from fibonacci_photonics.config.schema import (
    BookConvergence,
    ConvergenceConfig,
    CostStudy,
    DrudeConvergence,
    NearPoleStudy,
)
from fibonacci_photonics.io.paths import relative, result_dirs
from fibonacci_photonics.io.results import dump_json, load_json
from fibonacci_photonics.parallel import parallel_map, single_threaded_blas
from fibonacci_photonics.physics.electromagnetics import Polarization
from fibonacci_photonics.physics.fibonacci import n_layers_total
from fibonacci_photonics.physics.superlattice import SuperlatticeSpec
from fibonacci_photonics.physics.units import omega_from_nu_ghz
from fibonacci_photonics.solvers.pwe import PWESolver, n_max_for
from fibonacci_photonics.solvers.pwe.eigenfrequency import nondispersive_bands
from fibonacci_photonics.solvers.pwe.fourier import bilayer_fourier
from fibonacci_photonics.solvers.tmm import TMMSolver
from fibonacci_photonics.solvers.tmm.transfer_matrix import semitrace_by_recurrence
from fibonacci_photonics.viz.style import PWE_FORMATS, apply_prb_style, save_figure

CONFIG = "pwe/convergence.yaml"
RULES = ("inverse", "laurent")
RULE_STYLE = {"inverse": {"color": "#c23b22", "marker": "o"}, "laurent": {"color": "#2e8b4a", "marker": "s"}}
RULE_LABEL = {"laurent": "Laurent (libro)", "inverse": "inversa (Li)"}
COST_REPEATS = 200


def _bilayer_r(f: ArrayLike, fractions: tuple[float, float], eps: tuple[float, float]) -> NDArray[np.float64]:
    """R(f) analítica de la bicapa en incidencia normal, f = ωa/2πc."""
    kappa = 2 * np.pi * np.asarray(f, dtype=float)
    q1, q2 = kappa * np.sqrt(eps[0]), kappa * np.sqrt(eps[1])
    ratio = 0.5 * (q1 * eps[1] / (q2 * eps[0]) + q2 * eps[0] / (q1 * eps[1]))
    return np.cos(q1 * fractions[0]) * np.cos(q2 * fractions[1]) - ratio * np.sin(q1 * fractions[0]) * np.sin(
        q2 * fractions[1]
    )


def book_study(cfg: BookConvergence) -> list[dict[str, Any]]:
    total = cfg.thickness_1_um + cfg.thickness_2_um
    fractions = (cfg.thickness_1_um / total, cfg.thickness_2_um / total)
    eps = (cfg.epsilon_1, cfg.epsilon_2)
    ks = np.array(cfg.k_over_pi) * np.pi
    # Frecuencias exactas: raíces de R(f) = cos k (TMM analítica) refinadas con Brent.
    grid = np.linspace(1e-4, 1.0, 20001)
    exact = []
    for k in ks:
        d = _bilayer_r(grid, fractions, eps) - np.cos(k)
        idx = np.flatnonzero(np.sign(d[1:]) != np.sign(d[:-1]))[: cfg.n_bands]
        exact.append([brentq(lambda f, k=k: _bilayer_r(f, fractions, eps) - np.cos(k), grid[i], grid[i + 1],
                             xtol=1e-15) for i in idx])
    exact_array = np.array(exact)
    rows = []
    for rule in ("laurent", "inverse"):
        for n_max in cfg.n_max_values:
            cell = bilayer_fourier(*fractions, n_max)
            f = nondispersive_bands(cell, eps, (1, 1), ks, rule=rule, n_bands=cfg.n_bands)
            rows.append({"rule": rule, "n_max": n_max, "n_plane_waves": 2 * n_max + 1,
                         "max_rel_error": float(np.max(np.abs(f - exact_array) / exact_array))})
    return rows


DrudeTask = tuple[SuperlatticeSpec, int, float, Polarization, int, str, NDArray[np.float64]]


def drude_case(task: DrudeTask) -> dict[str, Any]:
    """max|R_PWE − R_TMM| en bandas para un (m, polarización, h, regla)."""
    spec, m, theta, pol, h, rule, nu = task
    r_tmm = TMMSolver(spec).semitrace(m, theta, nu, pol).real
    solver = PWESolver(spec, harmonics_per_layer=h, rule=rule)  # type: ignore[arg-type]
    r_pwe = solver.evaluator(m, theta, pol).many(nu)
    band = np.abs(r_tmm) <= 1
    return {"m": m, "polarization": pol.value, "harmonics_per_layer": h, "rule": rule,
            "n_plane_waves": solver.n_plane_waves(m),
            "max_abs_delta_r_in_bands": float(np.nanmax(np.abs(r_pwe - r_tmm)[band])),
            "classification_mismatch": int(np.count_nonzero((np.abs(r_pwe) <= 1) != band))}


def drude_study(cfg: DrudeConvergence) -> list[dict[str, Any]]:
    spec = load_spec(cfg.physics)
    nu = np.linspace(cfg.frequency_min_ghz, cfg.frequency_max_ghz, cfg.frequency_points)
    tasks = [(spec, m, cfg.theta_rad, Polarization(p), h, rule, nu)
             for m in cfg.fibonacci_orders for p in cfg.polarizations
             for h in cfg.harmonics_per_layer for rule in RULES]
    tasks.sort(key=lambda t: -n_max_for(t[1], t[4]))
    return sorted(parallel_map(drude_case, tasks),
                  key=lambda r: (r["m"], r["polarization"], r["rule"], r["harmonics_per_layer"]))


def cost_study(cfg: CostStudy) -> list[dict[str, Any]]:
    spec = load_spec(cfg.physics)
    nu = np.linspace(*cfg.window_ghz, cfg.frequency_points)
    omega = omega_from_nu_ghz(nu)
    theta = cfg.theta_rad
    rows = []
    for m, h in cfg.cases:
        start = time.perf_counter()
        for _ in range(COST_REPEATS):
            r_tmm = semitrace_by_recurrence(spec, m, omega, theta, Polarization.TE).real
        tmm_seconds = (time.perf_counter() - start) / (COST_REPEATS * nu.size)
        pwe = PWESolver(spec, harmonics_per_layer=h)
        evaluator = pwe.evaluator(m, theta, Polarization.TE)
        start = time.perf_counter()
        r_pwe = evaluator.many(nu)
        pwe_seconds = (time.perf_counter() - start) / nu.size
        rows.append({"m": m, "layers": n_layers_total(m), "harmonics_per_layer": h,
                     "n_plane_waves": pwe.n_plane_waves(m), "pwe_seconds_per_frequency": pwe_seconds,
                     "tmm_seconds_per_frequency": tmm_seconds,
                     "max_abs_delta_r": float(np.max(np.abs(r_pwe - r_tmm))),
                     "max_abs_r_tmm": float(np.max(np.abs(r_tmm)))})
        print(f"  costo m={m} h={h}: N={rows[-1]['n_plane_waves']} PWE {pwe_seconds:.3f} s/frec, "
              f"TMM {tmm_seconds * 1e6:.1f} µs/frec, max|ΔR|={rows[-1]['max_abs_delta_r']:.1e}")
    return rows


def near_pole_study(cfg: NearPoleStudy) -> dict[str, Any]:
    spec = load_spec(cfg.physics)
    m, h, theta = cfg.fibonacci_order, cfg.harmonics_per_layer, cfg.theta_rad
    offsets = np.logspace(np.log10(cfg.offsets_ghz[0]), np.log10(cfg.offsets_ghz[1]), cfg.points)
    nu = spec.nu_m_ghz() - offsets
    r_tmm = TMMSolver(spec).semitrace(m, theta, nu, Polarization.TE).real
    r_pwe = PWESolver(spec, harmonics_per_layer=h).evaluator(m, theta, Polarization.TE).many(nu)
    return {"offsets_ghz": offsets.tolist(), "r_tmm": r_tmm.tolist(), "r_pwe": r_pwe.tolist(),
            "m": m, "harmonics_per_layer": h, "theta": theta}


def main(*, plot_only: bool = False) -> dict[str, Any]:
    single_threaded_blas()
    loaded = load_model(CONFIG, ConvergenceConfig)
    cfg = loaded.model
    dirs = result_dirs("pwe", "convergence")
    if plot_only:
        data = load_json(dirs["data"] / "summary.json")
        return plot_and_save(data["book"], data["drude"], data["cost"], data["near_pole"], cfg.drude, dirs,
                             loaded.path)
    print("1) bicapa del libro")
    book = book_study(cfg.book)
    print("2) superred de Drude")
    drude = drude_study(cfg.drude)
    print("3) costo (un proceso, un hilo)")
    cost = cost_study(cfg.cost)
    print("4) validez cerca de ν_m")
    pole = near_pole_study(cfg.near_pole)
    return plot_and_save(book, drude, cost, pole, cfg.drude, dirs, loaded.path)


def plot_and_save(
    book: list[dict[str, Any]],
    drude: list[dict[str, Any]],
    cost: list[dict[str, Any]],
    pole: dict[str, Any],
    dcfg: DrudeConvergence,
    dirs: dict[str, Path],
    config_path: Path,
) -> dict[str, Any]:
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
    outputs["book_convergence"] = save_figure(fig, dirs["output"], "book_convergence", formats=PWE_FORMATS)

    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.8), sharey=True)
    colors = {3: "#1a1a1a", 4: "#2c4d8c", 5: "#c44ec0"}
    for ax, pol in zip(axes, dcfg.polarizations, strict=False):
        for m in dcfg.fibonacci_orders:
            for rule, ls in (("inverse", "-"), ("laurent", "--")):
                rows = [r for r in drude if r["m"] == m and r["polarization"] == pol.value and r["rule"] == rule]
                ax.loglog([r["n_plane_waves"] for r in rows], [max(r["max_abs_delta_r_in_bands"], 1e-16) for r in rows],
                          color=colors.get(m, "gray"), linestyle=ls, marker=RULE_STYLE[rule]["marker"], markersize=3,
                          linewidth=0.8, label=f"m={m}, {RULE_LABEL[rule]}")
        ax.set_xlabel("número de ondas planas N")
        ax.set_title(f"Drude, {pol.value}, " + r"$\theta=\pi/6$")
    axes[0].set_ylabel(r"max $|R_{\rm PWE}-R_{\rm TMM}|$ en bandas")
    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels, fontsize=6.5, ncol=3, loc="lower center", frameon=False)
    fig.tight_layout(rect=(0, 0.13, 1, 1))
    outputs["drude_convergence"] = save_figure(fig, dirs["output"], "drude_convergence", formats=PWE_FORMATS)

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
    outputs["cost"] = save_figure(fig, dirs["output"], "cost_and_accuracy", formats=PWE_FORMATS)

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
    outputs["near_pole"] = save_figure(fig, dirs["output"], "validity_near_nu_m", formats=PWE_FORMATS)

    summary = {
        "config": relative(config_path),
        "book": book, "drude": drude, "cost": cost, "near_pole": pole,
        "outputs": {name: {k: relative(p) for k, p in paths.items()} for name, paths in outputs.items()},
    }
    dump_json(dirs["data"] / "summary.json", summary)
    print("Convergencia PWE escrita en", relative(dirs["output"]))
    return summary
