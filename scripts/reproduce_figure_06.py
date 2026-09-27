#!/usr/bin/env python3
"""Figura 6: bandas plasmon-polaritón frente al ángulo (ν vs θ).

El original muestra regiones de frecuencia permitida que nacen cerca de ν_m = 1 GHz
en θ = 0 y se ensanchan hacia abajo al aumentar θ. El espesor de cada región es el
bandwidth; el eje vertical es frecuencia, no Δν.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from _common import dump_json, figure_dirs, load_figure_config, write_figure_readme, write_run_sidecar
from fibonacci_tmm.bandwidth import bandwidth_versus_angle
from fibonacci_tmm.electromagnetics import Polarization
from fibonacci_tmm.fibonacci import n_layers_b
from fibonacci_tmm.plotting import BRANCH_COLORS, apply_prb_style, save_figure


def _stack_by_order(points, n_expected: int, nu_m_ghz: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Bordes de cada subbanda, numeradas de menor a mayor frecuencia.

    Para θ > 0 hay exactamente F_{m-2} subbandas, así que la subbanda i (negro abajo,
    luego rojo, verde, …) es la misma rama física en todos los ángulos. Un ángulo con
    otro número de intervalos se deja vacío. En θ = 0 todas colapsan a ν_m.
    """
    thetas = np.array([p.theta for p in points], dtype=float)
    lo = np.full((n_expected, thetas.size), np.nan)
    hi = np.full((n_expected, thetas.size), np.nan)
    for j, point in enumerate(points):
        if point.n_modes != n_expected:
            continue
        pairs = sorted(
            (center - 0.5 * width, center + 0.5 * width)
            for center, width in zip(point.centers_ghz, point.bandwidths_ghz, strict=True)
        )
        for i, (a, b) in enumerate(pairs):
            lo[i, j] = a
            hi[i, j] = b
    if thetas.size and thetas[0] == 0.0:
        lo[:, 0] = nu_m_ghz
        hi[:, 0] = nu_m_ghz
    return thetas, lo, hi


def _frequency_grid(cfg: dict, nu_m_ghz: float) -> np.ndarray:
    window_min = cfg["plasmon_window_ghz"][0]
    offset = cfg["nu_m_min_offset_ghz"]
    log_part = nu_m_ghz - np.logspace(
        np.log10(offset), np.log10(nu_m_ghz - window_min), cfg["frequency_points_log"]
    )
    uniform_part = np.linspace(
        cfg["frequency_uniform_min_ghz"], nu_m_ghz - offset, cfg["frequency_points_uniform"]
    )
    return np.unique(np.concatenate([log_part, uniform_part]))


def _fill_tracked(ax, theta: np.ndarray, lo: np.ndarray, hi: np.ndarray, color: str) -> None:
    """Rellena por tramos contiguos; no une huecos con una diagonal."""
    finite = np.isfinite(lo) & np.isfinite(hi)
    if not np.any(finite):
        return
    padded = np.concatenate(([False], finite, [False]))
    starts = np.where(~padded[:-1] & padded[1:])[0]
    ends = np.where(padded[:-1] & ~padded[1:])[0]
    for s, e in zip(starts, ends, strict=True):
        ax.fill_between(theta[s:e], lo[s:e], hi[s:e], color=color, alpha=0.95, linewidth=0)


def main() -> None:
    spec, cfg, cfg_path = load_figure_config("figure_06.yaml")
    dirs = figure_dirs("figure_06")
    polarization = Polarization[cfg["polarization"]]
    orders = list(cfg["fibonacci_orders"])
    thetas = np.linspace(cfg["theta_min_rad"], cfg["theta_max_rad"], cfg["theta_points"])
    window = tuple(cfg["plasmon_window_ghz"])
    nu_m = spec.nu_m_ghz()
    nu = _frequency_grid(cfg, nu_m)

    apply_prb_style()
    fig, axes = plt.subplots(3, 2, figsize=(7.0, 7.6), sharex=True, sharey=True)
    payload = {}
    panel_ids = ["a", "b", "c", "d", "e", "f"]
    for ax, m, panel in zip(axes.ravel(), orders, panel_ids, strict=True):
        points = bandwidth_versus_angle(
            spec, m, thetas, nu, polarization, window=window, min_points=2, merge_gap_ghz=None
        )
        expected = n_layers_b(m)
        incomplete = [round(float(np.degrees(p.theta)), 2) for p in points if p.theta > 0 and p.n_modes != expected]
        if incomplete:
            print(f"Aviso m={m}: sin F_(m-2)={expected} subbandas en θ(°) = {incomplete}")
        th, lo, hi = _stack_by_order(points, expected, nu_m)
        payload[str(m)] = [
            {
                "theta": p.theta,
                "n_modes": p.n_modes,
                "bandwidths_ghz": list(p.bandwidths_ghz),
                "centers_ghz": list(p.centers_ghz),
            }
            for p in points
        ]
        n_tracks = lo.shape[0]
        for mode_idx in range(n_tracks):
            color = BRANCH_COLORS[mode_idx % len(BRANCH_COLORS)]
            _fill_tracked(ax, th, lo[mode_idx], hi[mode_idx], color)
        ax.set_xlim(float(thetas[0]), float(thetas[-1]))
        ax.set_ylim(0.95, 1.002)
        ax.set_xticks([0.0, np.pi / 12, np.pi / 6, np.pi / 3])
        ax.set_xticklabels([r"$0$", r"$\pi/12$", r"$\pi/6$", r"$\pi/3$"])
        ax.text(0.04, 0.92, rf"({panel})", transform=ax.transAxes, va="top")
        ax.text(0.96, 0.92, rf"$m={m}$", transform=ax.transAxes, va="top", ha="right")
    for ax in axes[:, 0]:
        ax.set_ylabel("frequency (GHz)")
    for ax in axes[-1, :]:
        ax.set_xlabel(r"$\theta$")
    fig.tight_layout()
    paths = save_figure(fig, dirs["output"], "figure_06")
    for dest in paths.values():
        (dirs["reproduced"] / dest.name).write_bytes(dest.read_bytes())
    extra = {
        "config": str(cfg_path),
        "note": "Eje vertical: frecuencia de las subbandas. El espesor de cada region es el bandwidth.",
        "bandwidths": payload,
        "outputs": {k: str(p) for k, p in paths.items()},
    }
    write_run_sidecar(dirs["data"] / "metadata.json", spec, extra)
    dump_json(dirs["data"] / "summary.json", extra)
    write_figure_readme(
        dirs["root"] / "README.md",
        "Figura 6",
        "Bandas plasmon-polaritón vs ángulo: ν_min(θ)–ν_max(θ) relleno por subbanda.\n"
        "A θ=0 las bandas colapsan hacia ν_m = 1 GHz; al aumentar θ se abren hacia abajo.\n\n"
        "Script: `python scripts/reproduce_figure_06.py`",
    )
    print("Figura 6 escrita")


if __name__ == "__main__":
    main()
