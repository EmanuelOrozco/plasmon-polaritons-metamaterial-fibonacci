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
from fibonacci_tmm.plotting import BRANCH_COLORS, apply_prb_style, save_figure


def _pairs_at_point(point) -> list[tuple[float, float]]:
    pairs = [
        (center - 0.5 * width, center + 0.5 * width)
        for center, width in zip(point.centers_ghz, point.bandwidths_ghz, strict=True)
    ]
    pairs.sort(key=lambda item: item[0])
    return pairs


def _track_edges(points, nu_m_ghz: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Sigue subbandas por solapamiento en ν (no por el centro).

    Así, al fragmentarse una banda, el trozo más ancho conserva el color y el
    nuevo aparece encima, como en el PRB: negro abajo, luego rojo, verde, …
    En θ = 0 las ramas que existen a θ > 0 colapsan a ν_m.
    """
    thetas = np.array([p.theta for p in points], dtype=float)
    all_pairs = [_pairs_at_point(p) for p in points]
    n_theta = thetas.size
    n_tracks = max((len(pairs) for pairs in all_pairs), default=0)
    if n_tracks == 0:
        return thetas, np.empty((0, n_theta)), np.empty((0, n_theta))

    lo = np.full((n_tracks, n_theta), np.nan)
    hi = np.full((n_tracks, n_theta), np.nan)
    last_lo = np.full(n_tracks, np.nan)
    last_hi = np.full(n_tracks, np.nan)

    start = next((j for j, pairs in enumerate(all_pairs) if pairs), 0)
    for i, (a, b) in enumerate(all_pairs[start]):
        lo[i, start] = a
        hi[i, start] = b
        last_lo[i] = a
        last_hi[i] = b

    for j in range(start + 1, n_theta):
        pairs = all_pairs[j]
        if not pairs:
            continue
        used_tracks: set[int] = set()
        used_new: set[int] = set()
        candidates: list[tuple[float, int, int]] = []
        for ni, (a, b) in enumerate(pairs):
            for ti in range(n_tracks):
                if not np.isfinite(last_lo[ti]):
                    continue
                overlap = min(b, last_hi[ti]) - max(a, last_lo[ti])
                if overlap > 0.0:
                    candidates.append((-overlap, ni, ti))
        candidates.sort()
        assignment: dict[int, int] = {}
        for _, ni, ti in candidates:
            if ni in used_new or ti in used_tracks:
                continue
            assignment[ni] = ti
            used_new.add(ni)
            used_tracks.add(ti)
        for ni, (a, b) in enumerate(pairs):
            if ni in assignment:
                ti = assignment[ni]
            else:
                empty = [
                    t
                    for t in range(n_tracks)
                    if t not in used_tracks and not np.isfinite(last_lo[t])
                ]
                if not empty:
                    empty = [t for t in range(n_tracks) if t not in used_tracks]
                if not empty:
                    continue
                ti = empty[0]
                used_tracks.add(ti)
            lo[ti, j] = a
            hi[ti, j] = b
        last_lo = lo[:, j].copy()
        last_hi = hi[:, j].copy()
        # Conserva el último estado de tracks no vistos en este θ (bandas muy finas).
        for ti in range(n_tracks):
            if not np.isfinite(last_lo[ti]) and np.isfinite(lo[ti, j - 1]):
                last_lo[ti] = lo[ti, j - 1]
                last_hi[ti] = hi[ti, j - 1]

    if n_theta and thetas[0] == 0.0:
        for ti in range(n_tracks):
            if np.any(np.isfinite(lo[ti])):
                lo[ti, 0] = nu_m_ghz
                hi[ti, 0] = nu_m_ghz

    return thetas, lo, hi


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
    nu = np.linspace(window[0], min(window[1], spec.nu_m_ghz() - 1e-4), cfg["frequency_points"])
    nu_m = spec.nu_m_ghz()

    apply_prb_style()
    fig, axes = plt.subplots(3, 2, figsize=(7.0, 7.6), sharex=True, sharey=True)
    payload = {}
    panel_ids = ["a", "b", "c", "d", "e", "f"]
    for ax, m, panel in zip(axes.ravel(), orders, panel_ids, strict=True):
        points = bandwidth_versus_angle(
            spec, m, thetas, nu, polarization, window=window, min_points=2, merge_gap_ghz=4e-5
        )
        th, lo, hi = _track_edges(points, nu_m)
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
