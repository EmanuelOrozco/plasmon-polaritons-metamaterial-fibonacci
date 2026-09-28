"""Figuras 1–4: relaciones de dispersión ν(k), con cualquiera de los dos métodos.

- Figs. 1 y 3: S3 y S4 en cuatro ángulos (malla uniforme de 0.15 a 5 GHz).
- Fig. 2: m = 3..6 cerca de ν_m, con el conteo de subbandas F_{m−2}.
- Fig. 4: zoom de las subbandas de m = 3 y 4 en ventanas por ángulo.

Cada barrido lleva además el cierre de las curvas en los bordes |R| = 1 y en
las puntas, refinados con el mismo método (``analysis.band_closure``).
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from numpy.typing import NDArray

from fibonacci_photonics.analysis.band_closure import BandEdges, close_band_edges
from fibonacci_photonics.analysis.plasmon_modes import detect_plasmon_modes, merge_touching_intervals
from fibonacci_photonics.physics.effective_medium import zero_average_index_frequency_ghz
from fibonacci_photonics.physics.electromagnetics import Polarization
from fibonacci_photonics.physics.fibonacci import n_layers_b
from fibonacci_photonics.physics.superlattice import SuperlatticeSpec
from fibonacci_photonics.reproduction.common import FigureRun
from fibonacci_photonics.solvers.base import DispersionSolver
from fibonacci_photonics.solvers.scan import DispersionScan, FrequencyInterval, scan_from_semitrace
from fibonacci_photonics.viz.dispersion import format_dispersion_axes, plot_dispersion_branches
from fibonacci_photonics.viz.style import M_COLORS, ORDER_STYLES, apply_prb_style, latex_theta

PANEL_LETTERS = "abcdef"
MODE_COUNT_MERGE_GAP_GHZ = {"figure_02": 3e-3, "figure_04": 2e-5}
"""Huecos de malla que se funden al contar subbandas en cada figura."""
MODE_COUNT_MIN_POINTS = 4


@dataclass(frozen=True)
class SolvedScan:
    """Barrido de un (m, θ) con sus tiempos y, en el PWE, los puntos de cierre."""

    scan: DispersionScan
    seconds: float
    edges: BandEdges | None = None
    edge_seconds: float = 0.0


def solve_scan(run: FigureRun, solver: DispersionSolver, m: int, theta: float, nu: NDArray[np.float64]) -> SolvedScan:
    """Barrido cronometrado y, si la figura lo pide, bordes |R| = 1 con el mismo solver."""
    start = time.perf_counter()
    scan = solver.scan(m, theta, nu, run.polarization)
    seconds = time.perf_counter() - start
    if not run.close_band_edges:
        return SolvedScan(scan, seconds)
    start = time.perf_counter()
    workers = run.workers if run.is_pwe else 1
    edges = close_band_edges(solver, scan, workers=workers, xtol_fraction=run.edge_xtol_fraction)
    return SolvedScan(scan, seconds, edges, time.perf_counter() - start)


def load_solved(path: Path, spec: SuperlatticeSpec, polarization: Polarization, method: str) -> SolvedScan:
    """Lee un barrido guardado por ``save_scan``, con sus puntos de cierre si los tiene."""
    with np.load(path) as data:
        scan = scan_from_semitrace(
            spec, int(data["m"]), float(data["theta"]), data["nu_ghz"],
            data["r_real"] + 1j * data["r_imag"], polarization, method,
        )
        if "edge_nu_ghz" not in data:
            return SolvedScan(scan, float(data["seconds"]))
        edges = BandEdges(
            nu_ghz=data["edge_nu_ghz"],
            r=data["edge_r"],
            k_lm_over_pi=data["edge_k_lm_over_pi"],
            kind=data["edge_kind"].astype(np.int8),
            evaluations=int(data["edge_evaluations"]),
        )
        return SolvedScan(scan, float(data["seconds"]), edges, float(data["edge_seconds"]))


def branch_style(figure_id: str, m: int) -> dict[str, Any]:
    """Estilo de las curvas de cada figura; lo comparten TMM, PWE y la superposición."""
    if figure_id in ("figure_01", "figure_03"):
        return dict(ORDER_STYLES[m])
    if figure_id == "figure_02":
        return {"color": "#2c4d8c", "linestyle": "-", "linewidth": 0.9}
    return {"color": M_COLORS.get(m, "black"), "linestyle": "-", "linewidth": 1.0}


def scan_arrays(solved: SolvedScan) -> dict[str, Any]:
    """Contenido del ``.npz`` de un barrido (mismo esquema en ambos métodos)."""
    scan = solved.scan
    arrays: dict[str, NDArray[Any]] = {
        "nu_ghz": scan.nu_ghz,
        "r_real": np.real(scan.r),
        "r_imag": np.imag(scan.r),
        "abs_r": scan.abs_r,
        "allowed": scan.allowed.astype(np.uint8),
        "k_lm_over_pi": scan.k_lm_over_pi,
        "m": np.array(scan.m),
        "theta": np.array(scan.theta),
        "seconds": np.array(solved.seconds),
    }
    if solved.edges is not None:
        edges = solved.edges
        arrays |= {
            "edge_nu_ghz": edges.nu_ghz,
            "edge_r": edges.r,
            "edge_k_lm_over_pi": edges.k_lm_over_pi,
            "edge_kind": edges.kind,
            "edge_evaluations": np.array(edges.evaluations),
            "edge_seconds": np.array(solved.edge_seconds),
        }
    return arrays


def scan_record(solved: SolvedScan) -> dict[str, Any]:
    """Resumen de un barrido para ``summary.json``."""
    scan = solved.scan
    record: dict[str, Any] = {
        "seconds": solved.seconds,
        "n_allowed": int(scan.allowed.sum()),
        "clipped_near_unit": scan.clipped_near_unit,
        "max_abs_imag_r": float(np.max(np.abs(np.imag(scan.r)))),
    }
    if solved.edges is not None:
        edges = solved.edges
        record |= {
            "n_edges": edges.n_edges,
            "n_extrema": edges.n_extrema,
            "n_hidden_gaps": edges.n_hidden_gaps,
            "edge_evaluations": edges.evaluations,
            "edge_seconds": solved.edge_seconds,
        }
    return record


def describe(solved: SolvedScan) -> str:
    text = f"{solved.seconds:.2f} s"
    if solved.edges is not None:
        e = solved.edges
        text += (f", {e.n_edges} bordes, {e.n_extrema} puntas, {e.n_hidden_gaps} gaps ocultos "
                 f"en {solved.edge_seconds:.1f} s")
    return text


def save_scan(run: FigureRun, stem: str, solved: SolvedScan) -> None:
    np.savez_compressed(run.dirs["data"] / f"{stem}.npz", **scan_arrays(solved))


def plot_solved(ax: Axes, solved: SolvedScan, **style: Any) -> None:
    edges = solved.edges
    if edges is None:
        plot_dispersion_branches(ax, solved.scan, **style)
    else:
        plot_dispersion_branches(
            ax, solved.scan, edge_nu_ghz=edges.nu_ghz, edge_k_lm_over_pi=edges.k_lm_over_pi, **style
        )


def frequency_points(run: FigureRun) -> int:
    """Puntos de malla: los de la física (TMM) o los de la numérica del PWE."""
    return int(run.numerics.frequency_points if run.is_pwe else run.config.frequency_points)


def interval_records(intervals: list[FrequencyInterval]) -> list[dict[str, float]]:
    return [{"nu_min": i.nu_min_ghz, "nu_max": i.nu_max_ghz, "width": i.bandwidth_ghz} for i in intervals]


def counted_modes(run: FigureRun, scan: DispersionScan, window: tuple[float, float]) -> list[FrequencyInterval]:
    report = detect_plasmon_modes(
        scan, nu_min_ghz=window[0], nu_max_ghz=window[1], min_points=MODE_COUNT_MIN_POINTS
    )
    return merge_touching_intervals(list(report.intervals), gap_ghz=MODE_COUNT_MERGE_GAP_GHZ[run.figure_id])


def dispersion_panels(run: FigureRun, *, title: str, description: str) -> dict[str, Any]:
    """Figs. 1 y 3: un panel por ángulo con S3 (discontinua) y S4 (continua)."""
    cfg = run.config
    nu = np.linspace(cfg.frequency_min_ghz, cfg.frequency_max_ghz, frequency_points(run))
    solver = run.solver()
    apply_prb_style()
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 6.2), sharex=True, sharey=True)
    runs = []
    for ax, (label, theta), panel in zip(axes.ravel(), cfg.angles(), PANEL_LETTERS, strict=False):
        for m in cfg.fibonacci_orders:
            solved = solve_scan(run, solver, m, theta, nu)
            save_scan(run, f"m{m}_theta_{label.replace('/', '_')}", solved)
            plot_solved(ax, solved, **branch_style(run.figure_id, m))
            runs.append({"m": m, "theta": theta, "theta_label": label, **scan_record(solved)})
            print(f"  m={m} θ={label}: {describe(solved)}")
        format_dispersion_axes(ax, nu_min=0.0, nu_max=5.0, panel=panel, theta_label=label)
    fig.tight_layout()
    payload = {
        "n0_gap_closed_ghz": {str(m): zero_average_index_frequency_ghz(run.spec, m) for m in cfg.fibonacci_orders},
        "runs": runs,
    }
    return run.finish(fig, payload, title=title, description=description)


def mode_count_figure(run: FigureRun) -> dict[str, Any]:
    """Fig. 2: m = 3..6 a θ = π/3 cerca de ν_m, con F_{m−2} subbandas bajo ν_m."""
    cfg = run.config
    nu = np.linspace(cfg.frequency_min_ghz, cfg.frequency_max_ghz, frequency_points(run))
    theta = float(cfg.thetas_rad[0])
    window = tuple(cfg.plasmon_window_ghz)
    solver = run.solver()
    apply_prb_style()
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 6.0), sharex=True, sharey=True)
    reports = []
    for ax, m, panel in zip(axes.ravel(), cfg.fibonacci_orders, PANEL_LETTERS, strict=False):
        solved = solve_scan(run, solver, m, theta, nu)
        save_scan(run, f"m{m}", solved)
        plot_solved(ax, solved, **branch_style(run.figure_id, m))
        ax.axhline(run.spec.nu_m_ghz(), color="#c23b22", linestyle="--", linewidth=0.8)
        format_dispersion_axes(ax, nu_min=2.0, nu_max=4.0, panel=panel, m_label=str(m))
        intervals = counted_modes(run, solved.scan, window)
        reports.append({"m": m, "expected": n_layers_b(m), "detected": len(intervals), **scan_record(solved),
                        "intervals": interval_records(intervals)})
        print(f"  m={m}: {len(intervals)} subbandas (esperadas {n_layers_b(m)}), {describe(solved)}")
    fig.tight_layout()
    return run.finish(
        fig,
        {"mode_counts": reports},
        title="Figura 2",
        description="Dispersión TE cerca de ν_m = 3 GHz. El número de subbandas debe ser F(m-2).",
    )


def zoom_figure(run: FigureRun) -> dict[str, Any]:
    """Fig. 4: subbandas plasmon-polaritón de m = 3 y 4 en ventanas por ángulo."""
    cfg = run.config
    solver = run.solver()
    panels = [(m, label, theta) for m in cfg.fibonacci_orders for label, theta in cfg.angles()]
    apply_prb_style()
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 6.0))
    summaries = []
    for ax, (m, label, theta), panel in zip(axes.ravel(), panels, PANEL_LETTERS, strict=False):
        window = cfg.frequency_windows_ghz[label]
        nu = np.linspace(window[0], window[1], frequency_points(run))
        solved = solve_scan(run, solver, m, theta, nu)
        save_scan(run, f"m{m}_{label.replace('/', '_')}", solved)
        plot_solved(ax, solved, **branch_style(run.figure_id, m))
        format_dispersion_axes(
            ax, nu_min=window[0], nu_max=window[1], panel=panel, theta_label=latex_theta(label), m_label=str(m)
        )
        intervals = counted_modes(run, solved.scan, window)
        summaries.append({"m": m, "theta": label, "detected": len(intervals), **scan_record(solved),
                          "intervals": interval_records(intervals)})
        print(f"  m={m} θ={label}: {len(intervals)} subbandas, {describe(solved)}")
    fig.tight_layout()
    return run.finish(
        fig,
        {"modes": summaries},
        title="Figura 4",
        description="Zoom de modos plasmon-polaritón. m=3: 1 subbanda; m=4: 2.",
    )


def figure_01(run: FigureRun) -> dict[str, Any]:
    return dispersion_panels(
        run,
        title="Figura 1",
        description=(
            "Relación de dispersión TE ν(k) para celdas S3 (discontinua) y S4 (continua).\n\n"
            "Parámetros del paper: a = b = 12 mm, εA = μA = 1, ωe/2π = ωm/2π = 3 GHz.\n"
            "Ángulos: 0, π/12, π/6, π/3."
        ),
    )


def figure_03(run: FigureRun) -> dict[str, Any]:
    return dispersion_panels(
        run, title="Figura 3", description="Como Fig. 1 con ν_m = 1 GHz dentro del gap ⟨n⟩=0."
    )
