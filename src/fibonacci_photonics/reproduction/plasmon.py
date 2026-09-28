"""Figuras 5 y 6: subbandas plasmon-polaritón frente al orden m y al ángulo θ.

Cada método usa la estrategia que su costo permite (``analysis.plasmon_bands``):
la TMM, una malla uniforme densa con bordes refinados; el PWE, una malla
logarítmica en ν_m − ν con localización adaptativa y validación de cada banda
con 1.5× ondas planas. El PWE se limita a m ≤ 6 por su costo O(N³).
"""

from __future__ import annotations

import time
from collections.abc import Sequence
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from numpy.typing import NDArray

from fibonacci_photonics.analysis.band_edges import plasmon_grid
from fibonacci_photonics.analysis.bandwidth import bandwidth_versus_angle
from fibonacci_photonics.analysis.plasmon_bands import adaptive_bands, dense_grid_bands
from fibonacci_photonics.parallel import parallel_map
from fibonacci_photonics.physics.fibonacci import n_layers_b
from fibonacci_photonics.reproduction.common import FigureRun
from fibonacci_photonics.solvers.pwe import PWESolver
from fibonacci_photonics.solvers.scan import FrequencyInterval
from fibonacci_photonics.viz.style import BRANCH_COLORS, apply_prb_style, latex_theta

PLASMON_WINDOW_TOP_OFFSET_GHZ = 1e-4
"""La malla densa de la Fig. 5 (TMM) termina 0.1 MHz bajo ν_m."""
ANGLE_TICKS = ([0.0, np.pi / 12, np.pi / 6, np.pi / 3], [r"$0$", r"$\pi/12$", r"$\pi/6$", r"$\pi/3$"])


def band_records(intervals: Sequence[FrequencyInterval]) -> list[dict[str, float]]:
    return [{"min": i.nu_min_ghz, "max": i.nu_max_ghz} for i in intervals]


def pwe_plasmon_task(task: dict[str, Any]) -> dict[str, Any]:
    """Subbandas de un (m, θ) con el PWE: malla logarítmica, localización y validación."""
    solver: PWESolver = task["solver"]
    m, theta = task["m"], task["theta"]
    grid = plasmon_grid(task["nu_min"], solver.spec.nu_m_ghz(), task["grid_points"], task["min_offset"])
    start = time.perf_counter()
    result = adaptive_bands(solver, solver.refined(task["check_factor"]), m, theta, task["polarization"], grid)
    return {
        "m": m,
        "theta": theta,
        "harmonics_per_layer": solver.harmonics(m),
        "n_max": solver.n_max(m),
        "bands": band_records(result.bands),
        "rejected_bands": band_records(result.rejected),
        "evaluations": result.evaluations,
        "seconds": time.perf_counter() - start,
    }


def pwe_tasks(run: FigureRun, cases: list[tuple[int, float]], nu_min: float) -> list[dict[str, Any]]:
    numerics = run.numerics
    solver = run.solver(workers=1)
    return [
        {
            "solver": solver,
            "m": m,
            "theta": theta,
            "nu_min": nu_min,
            "grid_points": numerics.grid_points,
            "min_offset": numerics.nu_m_min_offset_ghz,
            "check_factor": numerics.check_factor,
            "polarization": run.polarization,
        }
        for m, theta in cases
    ]


def _label_panel(ax: Axes, panel: str, text: str, *, y: float) -> None:
    ax.text(0.04, y, rf"({panel})", transform=ax.transAxes, va="top")
    ax.text(0.96, y, text, transform=ax.transAxes, va="top", ha="right")


# ------------------------------------------------------------------ Figura 5
def _order_bands_tmm(run: FigureRun, label: str, theta: float) -> list[dict[str, Any]]:
    cfg = run.config
    window = cfg.frequency_windows_ghz[label]
    top = min(window[1], run.spec.nu_m_ghz() - PLASMON_WINDOW_TOP_OFFSET_GHZ)
    nu = np.linspace(window[0], top, cfg.frequency_points)
    solver = run.solver()
    rows = []
    for m in cfg.fibonacci_orders:
        # Sin fusionar: a θ = π/12 hay gaps físicos de apenas ~0.4 kHz (m = 7).
        intervals = dense_grid_bands(solver, m, theta, run.polarization, nu, window, min_points=3)
        rows.append({"m": m, "n": len(intervals), "expected": n_layers_b(m), "bands": band_records(intervals)})
    return rows


def order_bars(ax: Axes, rows: Sequence[dict[str, Any]], *, offset: float = 0.0, linewidth: float = 4.0) -> None:
    """Una barra vertical por subbanda en x = m (+ ``offset``), un color por subbanda."""
    for row in rows:
        for idx, band in enumerate(sorted(row["bands"], key=lambda b: b["min"])):
            ax.plot([row["m"] + offset] * 2, [band["min"], band["max"]],
                    color=BRANCH_COLORS[idx % len(BRANCH_COLORS)], linewidth=linewidth, solid_capstyle="butt")


def format_order_panel(
    ax: Axes, panel: str, label: str, window: Sequence[float], orders: Sequence[int]
) -> None:
    ax.set_xlim(1.5, 8.5)
    ax.set_ylim(window[0], window[1])
    ax.set_xticks(list(orders))
    ax.set_xlabel("Fibonacci order")
    ax.set_ylabel("bandwidth (GHz)")
    _label_panel(ax, panel, rf"$\theta = {latex_theta(label)}$", y=0.95)


def figure_05(run: FigureRun) -> dict[str, Any]:
    """Fig. 5: segmentos de frecuencia permitida frente al orden de Fibonacci."""
    cfg = run.config
    rows_by_label: dict[str, list[dict[str, Any]]] = {}
    if run.is_pwe:
        tasks = []
        for label, theta in cfg.angles():
            cases = [(m, theta) for m in run.numerics.fibonacci_orders]
            nu_min = cfg.frequency_windows_ghz[label][0]
            tasks += [task | {"label": label} for task in pwe_tasks(run, cases, nu_min)]
        results = parallel_map(pwe_plasmon_task, tasks, run.workers)
        for task, res in zip(tasks, results, strict=True):
            window = cfg.frequency_windows_ghz[task["label"]]
            m = res["m"]
            bands = [b for b in res["bands"] if b["max"] >= window[0]]
            rows_by_label.setdefault(task["label"], []).append(
                {"m": m, "n": len(bands), "expected": n_layers_b(m), "bands": bands,
                 **{k: res[k] for k in ("harmonics_per_layer", "n_max", "rejected_bands", "evaluations", "seconds")}}
            )
            print(f"  θ={task['label']} m={m}: {len(bands)} subbandas (esperadas {n_layers_b(m)}), "
                  f"{len(res['rejected_bands'])} espurias descartadas, "
                  f"{res['evaluations']} evaluaciones, {res['seconds']:.0f} s")
    else:
        for label, theta in cfg.angles():
            rows_by_label[label] = _order_bands_tmm(run, label, theta)

    apply_prb_style()
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 4.4))
    for ax, (label, _), panel in zip(axes, cfg.angles(), "ab", strict=False):
        order_bars(ax, rows_by_label[label])
        format_order_panel(ax, panel, label, cfg.frequency_windows_ghz[label], cfg.fibonacci_orders)
        if run.is_pwe:
            orders = run.numerics.fibonacci_orders
            ax.text(0.96, 0.05, f"PWE: m ≤ {max(orders)}", transform=ax.transAxes, ha="right", fontsize=7,
                    color="#555555")
    fig.tight_layout()
    description = "Estructura de bandas plasmon-polaritón frente al orden de Fibonacci."
    if run.is_pwe:
        description += ("\nBordes refinados con Brent sobre |R_PWE| − 1; cada banda se valida con 1.5× ondas planas."
                        "\nm = 7, 8 se omiten por costo O(N³).")
    return run.finish(fig, {"intervals": rows_by_label}, title="Figura 5", description=description)


# ------------------------------------------------------------------ Figura 6
def angle_grid(run: FigureRun) -> NDArray[np.float64]:
    """Malla de la TMM: logarítmica en ν_m − ν más uniforme fina, sin duplicados.

    Cerca de θ = 0 las subbandas quedan a ~1e-5 GHz de ν_m con anchos de
    ~5e-8 GHz, y a θ ~ π/12 hay gaps físicos de ~4e-7 GHz.
    """
    cfg = run.config
    nu_m = run.spec.nu_m_ghz()
    offset = cfg.nu_m_min_offset_ghz
    span = nu_m - cfg.plasmon_window_ghz[0]
    log_part = nu_m - np.logspace(np.log10(offset), np.log10(span), cfg.frequency_points_log)
    uniform_part = np.linspace(cfg.frequency_uniform_min_ghz, nu_m - offset, cfg.frequency_points_uniform)
    return np.unique(np.concatenate([log_part, uniform_part]))


def stack_by_order(
    rows: Sequence[dict[str, Any]], n_expected: int, nu_m_ghz: float
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """Bordes de cada subbanda, numeradas de menor a mayor frecuencia.

    ``rows`` son registros ``{"theta": θ, "bands": [{"min": ν, "max": ν}, ...]}``
    como los de ``summary.json``. Para θ > 0 hay exactamente F_{m−2} subbandas,
    así que la subbanda i es la misma rama física en todos los ángulos; un ángulo
    con otro número de intervalos se deja vacío. En θ = 0 todas colapsan a ν_m.
    """
    thetas = np.array([row["theta"] for row in rows], dtype=float)
    lo = np.full((n_expected, thetas.size), np.nan)
    hi = np.full((n_expected, thetas.size), np.nan)
    for j, row in enumerate(rows):
        if len(row["bands"]) != n_expected:
            continue
        for i, band in enumerate(sorted(row["bands"], key=lambda item: item["min"])):
            lo[i, j] = band["min"]
            hi[i, j] = band["max"]
    if thetas.size and thetas[0] == 0.0:
        lo[:, 0] = nu_m_ghz
        hi[:, 0] = nu_m_ghz
    return thetas, lo, hi


def fill_subbands(ax: Axes, rows: Sequence[dict[str, Any]], m: int, nu_m_ghz: float) -> None:
    """Regiones permitidas frente a θ, un color por subbanda (panel de la Fig. 6 de la TMM)."""
    theta, lo, hi = stack_by_order(rows, n_layers_b(m), nu_m_ghz)
    for mode_idx in range(lo.shape[0]):
        _fill_tracked(ax, theta, lo[mode_idx], hi[mode_idx], BRANCH_COLORS[mode_idx % len(BRANCH_COLORS)])


def _fill_tracked(ax: Axes, theta: NDArray[np.float64], lo: NDArray[np.float64], hi: NDArray[np.float64],
                  color: str) -> None:
    """Rellena por tramos contiguos; no une huecos con una diagonal."""
    finite = np.isfinite(lo) & np.isfinite(hi)
    if not np.any(finite):
        return
    padded = np.concatenate(([False], finite, [False]))
    starts = np.where(~padded[:-1] & padded[1:])[0]
    ends = np.where(padded[:-1] & ~padded[1:])[0]
    for s, e in zip(starts, ends, strict=True):
        ax.fill_between(theta[s:e], lo[s:e], hi[s:e], color=color, alpha=0.95, linewidth=0)


def _warn_incomplete(m: int, rows: list[tuple[float, int]]) -> None:
    expected = n_layers_b(m)
    incomplete = [round(float(np.degrees(theta)), 2) for theta, n in rows if theta > 0 and n != expected]
    if incomplete:
        print(f"Aviso m={m}: sin F_(m-2)={expected} subbandas en θ(°) = {incomplete}")


def format_angle_panel(ax: Axes, panel: str, m: int) -> None:
    ax.set_ylim(0.95, 1.002)
    ax.set_xticks(ANGLE_TICKS[0])
    ax.set_xticklabels(ANGLE_TICKS[1])
    _label_panel(ax, panel, rf"$m={m}$", y=0.92)


def figure_06(run: FigureRun) -> dict[str, Any]:
    """Fig. 6: regiones de frecuencia permitida frente a θ, una por subbanda."""
    cfg = run.config
    apply_prb_style()
    fig, axes = plt.subplots(3, 2, figsize=(7.0, 7.6), sharex=True, sharey=True)
    payload: dict[str, list[dict[str, Any]]] = {}
    if run.is_pwe:
        numerics = run.numerics
        thetas = np.linspace(numerics.theta_min_rad, numerics.theta_max_rad, numerics.theta_points)
        orders = list(numerics.fibonacci_orders)
        # Primero las tareas caras para repartir mejor la carga.
        cases = [(m, float(theta)) for m in sorted(orders, reverse=True) for theta in thetas]
        results = parallel_map(pwe_plasmon_task, pwe_tasks(run, cases, cfg.plasmon_window_ghz[0]), run.workers)
        width = 0.6 * (thetas[1] - thetas[0])
        for ax, m, panel in zip(axes.ravel(), orders, "abcdef", strict=False):
            rows = sorted((r for r in results if r["m"] == m), key=lambda r: r["theta"])
            payload[str(m)] = rows
            _warn_incomplete(m, [(r["theta"], len(r["bands"])) for r in rows])
            for row in rows:
                for idx, band in enumerate(sorted(row["bands"], key=lambda b: b["min"])):
                    ax.bar(row["theta"], band["max"] - band["min"], bottom=band["min"], width=width,
                           color=BRANCH_COLORS[idx % len(BRANCH_COLORS)], linewidth=0)
            ax.set_xlim(0.0, float(cfg.theta_max_rad) + width)
            format_angle_panel(ax, panel, m)
        for ax in axes.ravel()[len(orders):]:
            ax.axis("off")
            ax.text(0.5, 0.5, "m = 7: omitido en PWE\n(costo O(N³), N ∝ F$_m$)", transform=ax.transAxes,
                    ha="center", va="center", fontsize=8, color="#555555")
        description = ("Subbandas plasmon-polaritón frente al ángulo con ondas planas (m = 2..6, θ = 5°..60°).\n"
                       "Cada barra va de ν_min a ν_max de una subbanda, con bordes refinados por Brent.")
    else:
        thetas = np.linspace(cfg.theta_min_rad, cfg.theta_max_rad, cfg.theta_points)
        nu = angle_grid(run)
        nu_m = run.spec.nu_m_ghz()
        solver = run.solver()
        for ax, m, panel in zip(axes.ravel(), cfg.fibonacci_orders, "abcdef", strict=False):
            points = bandwidth_versus_angle(
                solver, m, thetas, nu, run.polarization, window=tuple(cfg.plasmon_window_ghz), min_points=2,
                merge_gap_ghz=None, refine=True,
            )
            _warn_incomplete(m, [(p.theta, p.n_modes) for p in points])
            payload[str(m)] = [
                {"m": m, "theta": p.theta, "n_modes": p.n_modes, "bands": band_records(p.intervals),
                 "bandwidths_ghz": list(p.bandwidths_ghz), "centers_ghz": list(p.centers_ghz)}
                for p in points
            ]
            fill_subbands(ax, payload[str(m)], m, nu_m)
            ax.set_xlim(float(thetas[0]), float(thetas[-1]))
            format_angle_panel(ax, panel, m)
        description = ("Bandas plasmon-polaritón vs ángulo: ν_min(θ)–ν_max(θ) relleno por subbanda.\n"
                       "A θ=0 las bandas colapsan hacia ν_m = 1 GHz; al aumentar θ se abren hacia abajo.")
    for ax in axes[:, 0]:
        ax.set_ylabel("frequency (GHz)")
    for ax in axes[-1, :]:
        ax.set_xlabel(r"$\theta$")
    fig.tight_layout()
    payload_key = {"note": "Eje vertical: frecuencia de las subbandas. El espesor de cada region es el bandwidth.",
                   "bands": payload}
    return run.finish(fig, payload_key, title="Figura 6", description=description)
