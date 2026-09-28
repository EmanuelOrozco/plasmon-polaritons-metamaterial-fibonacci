#!/usr/bin/env python3
"""Eficiencia TMM vs PWE: tiempo, escalamiento, trabajo-precisión y memoria.

Todas las mediciones usan un núcleo y un hilo de BLAS, para que los tiempos de
los dos métodos sean comparables. Estudios:

(a) tiempo por frecuencia frente al orden de Fibonacci m (TMM por recurrencia
    de semitrazas, TMM por producto de capas y PWE con h armónicos por capa),
    con ajuste de potencias y extrapolación;
(b) diagrama trabajo-precisión: error máximo en R dentro de las bandas frente
    al tiempo por frecuencia, con una referencia TMM en precisión extendida;
(c) costo de un núcleo de cada figura del paper;
(d) un caso favorable al PWE: bandas ω(k) del cristal no dispersivo del libro;
(e) memoria del problema de autovalores del PWE.

Con ``--plot-only`` se redibuja desde ``data/summary.json`` sin medir de nuevo.
"""

from __future__ import annotations

import os

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_name] = "1"

import time  # noqa: E402

import sys  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.ticker import MaxNLocator  # noqa: E402
import yaml  # noqa: E402
from scipy.optimize import brentq  # noqa: E402

from _shared import ROOT  # noqa: E402
from _common import dump_json, load_figure_config, load_json, load_pwe_config, relative, result_dirs  # noqa: E402
from fibonacci_photonics.core.bands import omega_from_nu_ghz  # noqa: E402
from fibonacci_photonics.core.electromagnetics import Polarization  # noqa: E402
from fibonacci_photonics.core.fibonacci import n_layers_total  # noqa: E402
from fibonacci_photonics.core.params import spec_from_mapping  # noqa: E402
from fibonacci_photonics.plotting import apply_prb_style, save_figure  # noqa: E402
from fibonacci_photonics.pwe.bloch import n_max_for, solve_bloch  # noqa: E402
from fibonacci_photonics.pwe.eigenfrequency import nondispersive_bands  # noqa: E402
from fibonacci_photonics.pwe.fourier import bilayer_fourier  # noqa: E402
from fibonacci_photonics.tmm.transfer_matrix import semitrace_by_product, semitrace_by_recurrence  # noqa: E402

COMPARISON = ROOT / "results" / "comparison"
TABLES = COMPARISON / "tables"
COMPANION_BYTES_PER_N2 = 64
"""La matriz compañera es 2N×2N complex128: 4N² entradas de 16 bytes."""


def timed(func, repeats: int = 1) -> tuple[object, float]:
    start = time.perf_counter()
    for _ in range(repeats):
        out = func()
    return out, (time.perf_counter() - start) / repeats


def sci(value: float, digits: int = 1) -> str:
    if value == 0 or not np.isfinite(value):
        return "0" if value == 0 else "--"
    exponent = int(np.floor(np.log10(abs(value))))
    return rf"${value / 10**exponent:.{digits}f}\times10^{{{exponent}}}$"


def seconds_text(value: float) -> str:
    if value < 1e-3:
        return rf"{value * 1e6:.1f}~$\mu$s"
    if value < 1:
        return f"{value * 1e3:.1f}~ms"
    if value < 3600:
        return f"{value:.1f}~s"
    return f"{value / 3600:.1f}~h"


def write_table(name: str, header: str, rows: list[str], columns: str) -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    (TABLES / f"{name}.tex").write_text(
        f"\\begin{{tabular}}{{{columns}}}\n\\toprule\n{header} \\\\\n\\midrule\n" + "\n".join(rows)
        + "\n\\bottomrule\n\\end{tabular}\n",
        encoding="utf-8",
    )


def power_fit(x, y) -> tuple[float, float]:
    """y ≈ c·x^α por mínimos cuadrados en log-log; devuelve (c, α)."""
    alpha, log_c = np.polyfit(np.log(x), np.log(y), 1)
    return float(np.exp(log_c)), float(alpha)


def pwe_seconds_per_frequency(spec, m, theta, nu, polarization, harmonics) -> float:
    omega = omega_from_nu_ghz(nu)
    solve_bloch(spec, m, omega[:1], theta, polarization, n_max=n_max_for(m, harmonics))
    _, seconds = timed(lambda: solve_bloch(spec, m, omega, theta, polarization, n_max=n_max_for(m, harmonics)))
    return seconds / nu.size


# ---------------------------------------------------------------- (a) escalamiento
def scaling_study(spec, cfg, polarization) -> dict:
    theta = cfg["theta_rad"]
    nu_vec = np.linspace(*cfg["window_ghz"], cfg["tmm_frequencies"])
    nu_scalar = np.linspace(*cfg["window_ghz"], cfg["tmm_scalar_frequencies"])
    omega_vec, omega_scalar = omega_from_nu_ghz(nu_vec), omega_from_nu_ghz(nu_scalar)
    tmm = []
    with np.errstate(all="ignore"):
        for m in cfg["tmm_orders"]:
            _, t_rec = timed(lambda: semitrace_by_recurrence(spec, m, omega_vec, theta, polarization), repeats=20)
            _, t_prod = timed(lambda: semitrace_by_product(spec, m, omega_vec, theta, polarization), repeats=3)
            _, t_scalar = timed(lambda: [semitrace_by_recurrence(spec, m, w, theta, polarization) for w in omega_scalar])
            tmm.append({"m": m, "layers": n_layers_total(m),
                        "recurrence_vectorized_s": t_rec / nu_vec.size,
                        "product_vectorized_s": t_prod / nu_vec.size,
                        "recurrence_scalar_s": t_scalar / nu_scalar.size})
            print(f"  TMM m={m:2d} F_m={n_layers_total(m):6d}: recurrencia {tmm[-1]['recurrence_vectorized_s'] * 1e6:.2f} µs, "
                  f"producto {tmm[-1]['product_vectorized_s'] * 1e6:.2f} µs, escalar {tmm[-1]['recurrence_scalar_s'] * 1e6:.1f} µs")

    h = cfg["pwe_harmonics_per_layer"]
    counts = cfg["pwe_frequencies"]
    pwe = []
    for m in cfg["pwe_orders"]:
        nu = np.linspace(*cfg["window_ghz"], counts.get(m, counts["default"]))
        seconds = pwe_seconds_per_frequency(spec, m, theta, nu, polarization, h)
        n_pw = 2 * n_max_for(m, h) + 1
        pwe.append({"m": m, "layers": n_layers_total(m), "harmonics_per_layer": h, "n_plane_waves": n_pw,
                    "seconds_per_frequency": seconds, "companion_mb": COMPANION_BYTES_PER_N2 * n_pw**2 / 1e6})
        print(f"  PWE m={m:2d} N={n_pw:5d}: {seconds:.4f} s/frecuencia")

    big = [r for r in pwe if r["n_plane_waves"] >= 200]
    c_pwe, alpha_pwe = power_fit([r["n_plane_waves"] for r in big], [r["seconds_per_frequency"] for r in big])
    big_tmm = [r for r in tmm if r["layers"] >= 50]
    c_prod, beta_prod = power_fit([r["layers"] for r in big_tmm], [r["product_vectorized_s"] for r in big_tmm])
    extrapolated = []
    for m in cfg["extrapolate_orders"]:
        n_pw = 2 * n_max_for(m, h) + 1
        extrapolated.append({"m": m, "layers": n_layers_total(m), "n_plane_waves": n_pw,
                             "seconds_per_frequency": c_pwe * n_pw**alpha_pwe,
                             "companion_mb": COMPANION_BYTES_PER_N2 * n_pw**2 / 1e6})
    return {"theta": theta, "tmm": tmm, "pwe": pwe, "extrapolated": extrapolated,
            "fit_pwe": {"c": c_pwe, "alpha": alpha_pwe, "variable": "N"},
            "fit_tmm_product": {"c": c_prod, "alpha": beta_prod, "variable": "F_m"}}


# ---------------------------------------------------------- (b) trabajo-precisión
def semitrace_extended(spec, m: int, nu_ghz, theta: float) -> np.ndarray:
    """R_m en TE con aritmética de precisión extendida (long double).

    R_m es par en cada Q_j, así que no importa la rama de la raíz cuadrada.
    """
    ld = np.longdouble
    omega = ld(2) * ld(np.pi) * np.asarray(nu_ghz, dtype=ld) * ld(1e9)
    kappa2 = (omega / ld(spec.speed_of_light)) ** 2
    drude = spec.medium_b
    eps_a, mu_a = ld(np.real(spec.medium_a.epsilon)), ld(np.real(spec.medium_a.mu))
    eps_b = ld(drude.epsilon_0) - ld(drude.omega_e) ** 2 / omega**2
    mu_b = ld(drude.mu_0) - ld(drude.omega_m) ** 2 / omega**2
    q2 = kappa2 * eps_a * mu_a * np.sin(ld(theta)) ** 2
    q_a = np.sqrt((kappa2 * eps_a * mu_a - q2).astype(np.clongdouble))
    q_b = np.sqrt((kappa2 * eps_b * mu_b - q2).astype(np.clongdouble))
    a, b = ld(spec.thickness_a), ld(spec.thickness_b)
    r0, r1 = np.cos(q_b * b), np.cos(q_a * a)
    ratio = (q_a * mu_b / (q_b * mu_a) + q_b * mu_a / (q_a * mu_b)) / 2
    r2 = np.cos(q_a * a) * np.cos(q_b * b) - ratio * np.sin(q_a * a) * np.sin(q_b * b)
    seq = [r0, r1, r2]
    for _ in range(3, m + 1):
        seq.append(2 * seq[-1] * seq[-2] - seq[-3])
    return seq[m].real


def work_precision_study(spec, cfg, polarization) -> dict:
    m, theta = cfg["fibonacci_order"], cfg["theta_rad"]
    nu = np.linspace(cfg["frequency_min_ghz"], cfg["frequency_max_ghz"], cfg["frequency_points"])
    omega = omega_from_nu_ghz(nu)
    reference = semitrace_extended(spec, m, nu, theta).astype(np.float64)
    band = np.abs(reference) <= 1

    def error(r):
        return float(np.max(np.abs(np.asarray(r).real - reference)[band]))

    rows = []
    r, t = timed(lambda: semitrace_by_recurrence(spec, m, omega, theta, polarization), repeats=200)
    rows.append({"method": "TMM, recurrencia (vectorizada)", "n_plane_waves": None,
                 "seconds_per_frequency": t / nu.size, "max_abs_error": error(r)})
    r, t = timed(lambda: semitrace_by_product(spec, m, omega, theta, polarization), repeats=100)
    rows.append({"method": "TMM, producto de capas (vectorizado)", "n_plane_waves": None,
                 "seconds_per_frequency": t / nu.size, "max_abs_error": error(r)})
    r, t = timed(lambda: np.array([semitrace_by_recurrence(spec, m, w, theta, polarization) for w in omega]), repeats=5)
    rows.append({"method": "TMM, recurrencia (una frecuencia por llamada)", "n_plane_waves": None,
                 "seconds_per_frequency": t / nu.size, "max_abs_error": error(r)})
    for h in cfg["harmonics_per_layer"]:
        n_max = n_max_for(m, h)
        solve_bloch(spec, m, omega[:1], theta, polarization, n_max=n_max)
        sol, t = timed(lambda: solve_bloch(spec, m, omega, theta, polarization, n_max=n_max))
        rows.append({"method": "PWE", "harmonics_per_layer": h, "n_plane_waves": 2 * n_max + 1,
                     "seconds_per_frequency": t / nu.size, "max_abs_error": error(sol.r)})
    for row in rows:
        print(f"  {row['method']:<46s} N={row['n_plane_waves']!s:>5s}: {row['seconds_per_frequency']:.2e} s/frec, "
              f"error {row['max_abs_error']:.1e}")
    return {"m": m, "theta": theta, "n_frequencies": int(nu.size), "n_in_bands": int(band.sum()), "rows": rows}


# ------------------------------------------------------------- (c) costo por figura
def figure_costs(cfg) -> dict:
    per_n: dict[tuple[int, int], float] = {}
    rows = []
    for figure_id in ("figure_01", "figure_02", "figure_03", "figure_04"):
        spec, physics, numerics, _ = load_pwe_config(figure_id)
        polarization = Polarization[physics["polarization"]]
        summary = load_json(COMPARISON / figure_id / "data" / "summary.json")
        nu = np.linspace(0.5, 5.0, cfg["sample_frequencies"])
        evaluations, tmm_s, wall_s, single_s = 0, 0.0, 0.0, 0.0
        for scan in summary["scans"]:
            value = numerics["harmonics_per_layer"]
            h = int(value[scan["m"]]) if isinstance(value, dict) else int(value)
            key = (scan["m"], h)
            if key not in per_n:
                per_n[key] = pwe_seconds_per_frequency(spec, scan["m"], scan["theta"], nu, polarization, h)
            evaluations += scan["n_points"]
            tmm_s += scan["tmm_seconds"]
            wall_s += scan["pwe_seconds"]
            single_s += scan["n_points"] * per_n[key]
        rows.append({"figure": figure_id, "pwe_evaluations": evaluations, "tmm_seconds": tmm_s,
                     "pwe_single_core_seconds": single_s, "pwe_single_core_estimated": True,
                     "pwe_wall_seconds": wall_s, "pwe_workers": cfg["parallel_workers_used"]})
    for figure_id in ("figure_05", "figure_06"):
        cases = load_json(COMPARISON / figure_id / "data" / "summary.json")["cases"]
        flat = [c for group in cases.values() for c in group]
        rows.append({"figure": figure_id, "pwe_evaluations": sum(c["pwe_evaluations"] for c in flat),
                     "tmm_seconds": sum(c["tmm_seconds"] for c in flat),
                     "pwe_single_core_seconds": sum(c["pwe_seconds"] for c in flat),
                     "pwe_single_core_estimated": False, "pwe_wall_seconds": None,
                     "pwe_workers": cfg["parallel_workers_used"]})
    for row in rows:
        row["ratio"] = row["pwe_single_core_seconds"] / row["tmm_seconds"]
        print(f"  {row['figure']}: {row['pwe_evaluations']} evaluaciones, TMM {row['tmm_seconds']:.3g} s, "
              f"PWE 1 núcleo {row['pwe_single_core_seconds']:.3g} s (×{row['ratio']:.0f})")
    return {"rows": rows, "per_frequency_single_core": [
        {"m": m, "harmonics_per_layer": h, "seconds": s} for (m, h), s in sorted(per_n.items())]}


# --------------------------------------------------- (d) caso favorable al PWE
def book_study(cfg) -> dict:
    book = yaml.safe_load((ROOT / "configs" / cfg["config"]).read_text(encoding="utf-8"))
    l1, l2 = book["thickness_1_um"], book["thickness_2_um"]
    eps = (book["epsilon_1"], book["epsilon_2"])
    n_bands = book["n_bands"]
    period = (l1 + l2) * 1e-6
    spec = spec_from_mapping({
        "layer_a_thickness_mm": l1 * 1e-3, "layer_b_thickness_mm": l2 * 1e-3, "epsilon_a": eps[0],
        "omega_e_over_2pi_ghz": 0.0, "omega_m_over_2pi_ghz": 0.0, "epsilon_0": eps[1],
        "source": "Sukhoivanov & Guryev, cap. 4 (1D)"})
    to_omega = spec.speed_of_light / period / 1e9

    def r_of_f(f):
        return semitrace_by_recurrence(spec, 2, omega_from_nu_ghz(np.asarray(f) * to_omega), 0.0,
                                       Polarization.TM).real

    ks = np.linspace(0.025, 0.975, cfg["k_points"]) * np.pi
    cell = bilayer_fourier(l1 / (l1 + l2), l2 / (l1 + l2), book["n_max"])
    nondispersive_bands(cell, eps, (1, 1), ks[:1], rule="inverse", n_bands=n_bands)
    pwe, t_pwe = timed(lambda: nondispersive_bands(cell, eps, (1, 1), ks, rule="inverse", n_bands=n_bands), repeats=3)
    f_max = float(pwe[:, -1].max()) * 1.05

    def tmm_roots():
        grid = np.linspace(1e-6, f_max, cfg["grid_points"])
        r_grid = r_of_f(grid)
        out = np.full((ks.size, n_bands), np.nan)
        for i, k in enumerate(ks):
            g = r_grid - np.cos(k)
            idx = np.nonzero(np.sign(g[:-1]) * np.sign(g[1:]) < 0)[0][:n_bands]
            out[i, : idx.size] = [brentq(lambda f: float(r_of_f(f)) - np.cos(k), grid[j], grid[j + 1],
                                         xtol=1e-14, rtol=1e-15) for j in idx]
        return out

    tmm, t_tmm = timed(tmm_roots, repeats=3)
    grid = np.linspace(1e-6, f_max, 60001)
    _, t_tmm_grid = timed(lambda: np.arccos(np.clip(r_of_f(grid), -1, 1)), repeats=20)
    rel = np.abs(pwe - tmm) / tmm
    out = {"n_plane_waves": 2 * book["n_max"] + 1, "k_points": int(ks.size), "n_bands": n_bands,
           "pwe_seconds": t_pwe, "tmm_root_seconds": t_tmm, "tmm_grid_seconds": t_tmm_grid,
           "tmm_grid_points": int(grid.size), "pwe_max_rel_error": float(np.nanmax(rel)),
           "tmm_missing_roots": int(np.isnan(tmm).sum())}
    print(f"  libro: PWE {t_pwe * 1e3:.1f} ms (error rel. {out['pwe_max_rel_error']:.1e}), "
          f"TMM raíces {t_tmm * 1e3:.1f} ms, TMM malla {t_tmm_grid * 1e3:.1f} ms")
    return out


# ------------------------------------------------------------------- figura y tablas
def plot(scaling, wp, dirs) -> dict:
    apply_prb_style()
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.4))

    ax = axes[0]
    m_tmm = [r["m"] for r in scaling["tmm"]]
    ax.semilogy(m_tmm, [r["recurrence_vectorized_s"] for r in scaling["tmm"]], "o-", color="#2c4d8c",
                ms=3, label="TMM, recurrencia")
    ax.semilogy(m_tmm, [r["recurrence_scalar_s"] for r in scaling["tmm"]], "^-", color="#6f8fc9",
                ms=3, label="TMM, recurrencia (escalar)")
    ax.semilogy(m_tmm, [r["product_vectorized_s"] for r in scaling["tmm"]], "s-", color="#2e8b4a",
                ms=3, label="TMM, producto de capas")
    m_pwe = [r["m"] for r in scaling["pwe"]]
    ax.semilogy(m_pwe, [r["seconds_per_frequency"] for r in scaling["pwe"]], "D-", color="#c23b22",
                ms=3, label=f"PWE, $h={scaling['pwe'][0]['harmonics_per_layer']}$")
    ext = scaling["extrapolated"]
    ax.semilogy([m_pwe[-1]] + [r["m"] for r in ext],
                [scaling["pwe"][-1]["seconds_per_frequency"]] + [r["seconds_per_frequency"] for r in ext],
                "D--", color="#c23b22", mfc="white", ms=3, label="PWE extrapolado")
    ax.set_xlabel("orden de Fibonacci $m$")
    ax.set_ylabel("tiempo por frecuencia (s)")
    ax.set_title("(a) escalamiento, un núcleo", fontsize=8)
    ax.legend(fontsize=5.8, loc="upper left")

    ax = axes[1]
    pwe_rows = [r for r in wp["rows"] if r["method"] == "PWE"]
    ax.loglog([r["seconds_per_frequency"] for r in pwe_rows], [r["max_abs_error"] for r in pwe_rows],
              "D-", color="#c23b22", ms=3, label="PWE ($h$ creciente)")
    for r in pwe_rows:
        ax.annotate(f"N={r['n_plane_waves']}", (r["seconds_per_frequency"], r["max_abs_error"]),
                    fontsize=5.5, xytext=(3, 2), textcoords="offset points")
    markers = {"TMM, recurrencia (vectorizada)": ("o", "#2c4d8c"),
               "TMM, producto de capas (vectorizado)": ("s", "#2e8b4a"),
               "TMM, recurrencia (una frecuencia por llamada)": ("^", "#6f8fc9")}
    for r in wp["rows"]:
        if r["method"] in markers:
            marker, color = markers[r["method"]]
            ax.loglog(r["seconds_per_frequency"], max(r["max_abs_error"], 1e-17), marker, color=color, ms=5,
                      label=r["method"])
    ax.set_xlabel("tiempo por frecuencia (s)")
    ax.set_ylabel(r"$\max|R-R_{\rm ref}|$ en bandas")
    ax.set_title(rf"(b) trabajo-precisión, $m={wp['m']}$, $\theta=\pi/6$", fontsize=8)
    ax.legend(fontsize=5.5, loc="center right")

    ax = axes[2]
    all_pwe = scaling["pwe"] + ext
    ax.semilogy([r["m"] for r in scaling["pwe"]], [r["companion_mb"] for r in scaling["pwe"]], "D-",
                color="#c23b22", ms=3, label="PWE, matriz compañera")
    ax.semilogy([scaling["pwe"][-1]["m"]] + [r["m"] for r in ext],
                [scaling["pwe"][-1]["companion_mb"]] + [r["companion_mb"] for r in ext], "D--",
                color="#c23b22", mfc="white", ms=3)
    tmm_mb = 4 * 16 / 1e6
    ax.semilogy([all_pwe[0]["m"], all_pwe[-1]["m"]], [tmm_mb, tmm_mb], "-", color="#2c4d8c",
                label="TMM, matriz $2\\times2$")
    ax.axhline(16e3, color="gray", linestyle=":", linewidth=0.8)
    ax.text(all_pwe[0]["m"], 16e3 / 4, "16 GB", fontsize=6, color="gray")
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlabel("orden de Fibonacci $m$")
    ax.set_ylabel("memoria por frecuencia (MB)")
    ax.set_title("(c) memoria", fontsize=8)
    ax.legend(fontsize=6, loc="center right")
    fig.tight_layout()
    return save_figure(fig, dirs["output"], "efficiency", formats=("pdf", "png"))


def tables(scaling, wp, figures, book) -> None:
    tmm_by_m = {r["m"]: r for r in scaling["tmm"]}
    rows = []
    for r in scaling["pwe"] + scaling["extrapolated"]:
        t = tmm_by_m[r["m"]]
        mark = "$^*$" if r in scaling["extrapolated"] else ""
        rows.append(f"{r['m']}{mark} & {r['layers']} & {r['n_plane_waves']} & "
                    f"{t['recurrence_vectorized_s'] * 1e6:.2f} & {t['recurrence_scalar_s'] * 1e6:.1f} & "
                    f"{seconds_text(r['seconds_per_frequency'])} & "
                    f"{sci(r['seconds_per_frequency'] / t['recurrence_scalar_s'], 0)} & "
                    f"{r['companion_mb']:.{2 if r['companion_mb'] < 10 else 0}f} \\\\")
    write_table("efficiency_scaling",
                r"$m$ & capas $F_m$ & $N$ & TMM vect. ($\mu$s) & TMM escalar ($\mu$s) & PWE & PWE/TMM escalar"
                r" & memoria PWE (MB)", rows, "rrrrrrrr")

    rows = []
    for r in wp["rows"]:
        label = r["method"] if r["method"] != "PWE" else rf"PWE, $h={r['harmonics_per_layer']}$, $N={r['n_plane_waves']}$"
        rows.append(f"{label} & {seconds_text(r['seconds_per_frequency'])} & {sci(max(r['max_abs_error'], 1e-17))} \\\\")
    write_table("efficiency_work_precision", r"Método & tiempo por frecuencia & $\max|R-R_{\rm ref}|$ en bandas",
                rows, "lrr")

    rows = []
    for r in figures["rows"]:
        single = seconds_text(r["pwe_single_core_seconds"]) + ("$^\\dagger$" if r["pwe_single_core_estimated"] else "")
        wall = seconds_text(r["pwe_wall_seconds"]) if r["pwe_wall_seconds"] is not None else "--"
        rows.append(f"{r['figure'][-1]} & {r['pwe_evaluations']} & {seconds_text(r['tmm_seconds'])} & {single} & "
                    f"{wall} & {sci(r['ratio'], 0)} \\\\")
    write_table("efficiency_figures",
                r"Fig. & evaluaciones PWE & TMM & PWE, 1 núcleo & PWE, pared (16 proc.) & PWE/TMM", rows, "rrrrrr")

    rows = [
        rf"PWE $\omega(k)$, regla inversa, $N={book['n_plane_waves']}$ & {seconds_text(book['pwe_seconds'])} & "
        rf"{sci(book['pwe_max_rel_error'])} \\",
        rf"TMM, raíces de $R(\omega)=\cos ka$ (Brent) & {seconds_text(book['tmm_root_seconds'])} & $<10^{{-13}}$ \\",
        rf"TMM, $k=\arccos R$ en {book['tmm_grid_points']} frecuencias & {seconds_text(book['tmm_grid_seconds'])}"
        r" & exacta en la malla \\",
    ]
    write_table("efficiency_book", rf"Método ({book['n_bands']} bandas, {book['k_points']} valores de $k$) & tiempo"
                r" & error relativo en $\omega$", rows, "lrr")


def sci_plain(value: float, digits: int = 0) -> str:
    exponent = int(np.floor(np.log10(abs(value))))
    mantissa = round(value / 10**exponent, digits)
    if mantissa >= 10:
        mantissa, exponent = mantissa / 10, exponent + 1
    return rf"{mantissa:.{digits}f}\times10^{{{exponent}}}"


def human_time(seconds: float) -> str:
    if seconds < 1e-3:
        return rf"{seconds * 1e6:.2g}~$\mu$s"
    if seconds < 1:
        return f"{seconds * 1e3:.{1 if seconds < 1e-2 else 0}f}~ms"
    if seconds < 60:
        return f"{seconds:.1f}~s"
    if seconds < 3600:
        return f"{seconds / 60:.0f}~min"
    return f"{seconds / 3600:.{0 if seconds >= 36000 else 1}f}~h"


def human_bytes(n_bytes: float) -> str:
    if n_bytes < 1e9:
        return f"{n_bytes / 1e6:.0f}~MB"
    return f"{n_bytes / 1e9:.{0 if n_bytes >= 1e10 else 1}f}~GB"


def write_macros(scaling, wp, figures, book) -> dict[str, str]:
    """Cifras del texto del anexo y del paper como macros LaTeX (``efficiency_macros.tex``)."""
    tmm, pwe, ext = scaling["tmm"], scaling["pwe"], scaling["extrapolated"]
    alpha = scaling["fit_pwe"]["alpha"]
    growth = [b["seconds_per_frequency"] / a["seconds_per_frequency"] for a, b in zip(pwe, pwe[1:]) if a["m"] >= 6]
    pwe_rows = [r for r in wp["rows"] if r["method"] == "PWE"]
    tmm_rows = [r for r in wp["rows"] if r["method"].startswith("TMM")]
    n_pw = np.array([r["n_plane_waves"] for r in pwe_rows], dtype=float)
    err_exp = -np.polyfit(np.log(n_pw), np.log([r["max_abs_error"] for r in pwe_rows]), 1)[0]
    time_exp = np.polyfit(np.log(n_pw), np.log([r["seconds_per_frequency"] for r in pwe_rows]), 1)[0]
    best = min(pwe_rows, key=lambda r: r["max_abs_error"])
    tmm_scalar = next(r for r in tmm_rows if "una frecuencia" in r["method"])
    target_n = best["n_plane_waves"] * (best["max_abs_error"] / 1e-12) ** (1 / err_exp)
    rows = figures["rows"]
    walls = [r["pwe_wall_seconds"] for r in rows if r["pwe_wall_seconds"] is not None]
    disp_tmm = [r["tmm_seconds"] for r in rows if r["pwe_wall_seconds"] is not None]
    ratios = [r["ratio"] for r in rows]
    by_m = {r["m"]: r for r in pwe + ext}
    macros = {
        "EffAlpha": f"{alpha:.2f}",
        "EffAlphaShort": f"{alpha:.1f}",
        "EffProdBeta": f"{scaling['fit_tmm_product']['alpha']:.2f}",
        "EffTmmRecMin": f"{min(r['recurrence_vectorized_s'] for r in tmm) * 1e6:.2f}",
        "EffTmmRecMax": f"{max(r['recurrence_vectorized_s'] for r in tmm) * 1e6:.2f}",
        "EffTmmScalar": f"{np.median([r['recurrence_scalar_s'] for r in tmm]) * 1e6:.0f}",
        "EffTmmMaxM": str(tmm[-1]["m"]),
        "EffTmmMaxLayers": f"{tmm[-1]['layers']:,}".replace(",", r"\,"),
        "EffProdFirst": f"{tmm[0]['product_vectorized_s'] * 1e6:.1f}",
        "EffProdLast": f"{tmm[-1]['product_vectorized_s'] * 1e6:.0f}",
        "EffPweFirstM": str(pwe[0]["m"]),
        "EffPweFirst": human_time(pwe[0]["seconds_per_frequency"]),
        "EffPweLastM": str(pwe[-1]["m"]),
        "EffPweLast": human_time(pwe[-1]["seconds_per_frequency"]),
        "EffPweGrowth": f"{np.median(growth):.1f}",
        "EffExtrapAM": str(ext[0]["m"]),
        "EffExtrapA": human_time(ext[0]["seconds_per_frequency"]),
        "EffExtrapBM": str(ext[-1]["m"]),
        "EffExtrapB": human_time(ext[-1]["seconds_per_frequency"]),
        "EffExtrapBMem": human_bytes(ext[-1]["companion_mb"] * 1e6),
        "EffMemEight": human_bytes(by_m[8]["companion_mb"] * 1e6) if 8 in by_m else "--",
        "EffWpM": str(wp["m"]),
        "EffWpTmmErr": sci_plain(max(r["max_abs_error"] for r in tmm_rows)),
        "EffWpTmmMin": f"{min(r['seconds_per_frequency'] for r in tmm_rows) * 1e6:.2g}",
        "EffWpTmmMax": f"{max(r['seconds_per_frequency'] for r in tmm_rows) * 1e6:.2g}",
        "EffWpErrExp": f"{err_exp:.1f}",
        "EffWpTimeExp": f"{time_exp:.1f}",
        "EffWpTradeExp": f"{err_exp / time_exp:.1f}",
        "EffWpBestErr": sci_plain(best["max_abs_error"]),
        "EffWpBestN": str(best["n_plane_waves"]),
        "EffWpBestTime": human_time(best["seconds_per_frequency"]),
        "EffWpSpeedRatio": sci_plain(best["seconds_per_frequency"] / tmm_scalar["seconds_per_frequency"]),
        "EffWpErrRatio": sci_plain(best["max_abs_error"] / tmm_scalar["max_abs_error"]),
        "EffWpTargetN": sci_plain(target_n),
        "EffWpTargetMem": human_bytes(COMPANION_BYTES_PER_N2 * target_n**2),
        "EffWpTargetTime": human_time(scaling["fit_pwe"]["c"] * target_n**alpha),
        "EffFigTmmTotal": f"{sum(r['tmm_seconds'] for r in rows):.2f}~s",
        "EffFigPweTotal": human_time(sum(r["pwe_single_core_seconds"] for r in rows)),
        "EffFigRatioMin": sci_plain(min(ratios)),
        "EffFigRatioMax": sci_plain(max(ratios)),
        "EffFigWallMin": f"{min(walls):.0f}",
        "EffFigWallMax": f"{max(walls):.0f}",
        "EffFigTmmMin": f"{min(disp_tmm) * 1e3:.1f}",
        "EffFigTmmMax": f"{max(disp_tmm) * 1e3:.1f}",
        "EffBookN": str(book["n_plane_waves"]),
        "EffBookPwe": human_time(book["pwe_seconds"]),
        "EffBookTmm": human_time(book["tmm_root_seconds"]),
        "EffBookGrid": human_time(book["tmm_grid_seconds"]),
        "EffBookGridPoints": f"{book['tmm_grid_points']:,}".replace(",", r"\,"),
        "EffBookErr": sci_plain(book["pwe_max_rel_error"]),
    }
    TABLES.mkdir(parents=True, exist_ok=True)
    (TABLES / "efficiency_macros.tex").write_text(
        "% Generado por scripts/comparison/efficiency.py: cifras de tiempo del anexo y del paper.\n"
        + "".join(f"\\providecommand{{\\{k}}}{{{v}}}\n" for k, v in macros.items()),
        encoding="utf-8",
    )
    return macros


def main() -> None:
    cfg = yaml.safe_load((ROOT / "configs" / "comparison" / "efficiency.yaml").read_text(encoding="utf-8"))
    spec, _, _ = load_figure_config(cfg["physics"])
    polarization = Polarization[cfg["polarization"]]
    dirs = result_dirs("comparison", "efficiency")
    if "--plot-only" in sys.argv:
        data = load_json(dirs["data"] / "summary.json")
        plot(data["scaling"], data["work_precision"], dirs)
        tables(data["scaling"], data["work_precision"], data["figures"], data["book"])
        write_macros(data["scaling"], data["work_precision"], data["figures"], data["book"])
        return
    print("(a) escalamiento con m")
    scaling = scaling_study(spec, cfg["scaling"], polarization)
    print("(b) trabajo-precisión")
    wp = work_precision_study(spec, cfg["work_precision"], polarization)
    print("(c) costo por figura")
    figures = figure_costs(cfg["figures"])
    print("(d) cristal del libro")
    book = book_study(cfg["book"])
    paths = plot(scaling, wp, dirs)
    tables(scaling, wp, figures, book)
    macros = write_macros(scaling, wp, figures, book)
    dump_json(dirs["data"] / "summary.json", {
        "config": "configs/comparison/efficiency.yaml", "threads": 1,
        "scaling": scaling, "work_precision": wp, "figures": figures, "book": book,
        "outputs": {k: relative(v) for k, v in paths.items()}, "macros": macros,
    })
    print("  ajuste PWE: t ∝ N^%.2f; TMM producto: t ∝ F_m^%.2f"
          % (scaling["fit_pwe"]["alpha"], scaling["fit_tmm_product"]["alpha"]))


if __name__ == "__main__":
    main()
