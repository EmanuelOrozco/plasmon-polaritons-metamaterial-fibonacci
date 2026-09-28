#!/usr/bin/env python3
"""Genera las imágenes conceptuales de la presentación en presentacion/imagenes/."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Arc, FancyArrowPatch, Rectangle

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from fibonacci_photonics.analysis.band_closure import close_band_edges  # noqa: E402
from fibonacci_photonics.config import load_figure, load_spec  # noqa: E402
from fibonacci_photonics.physics.effective_medium import average_epsilon_mu  # noqa: E402
from fibonacci_photonics.physics.fibonacci import n_layers_b, sequence  # noqa: E402
from fibonacci_photonics.physics.units import omega_from_nu_ghz  # noqa: E402
from fibonacci_photonics.solvers import TMMSolver  # noqa: E402
from fibonacci_photonics.viz.dispersion import plot_dispersion_branches  # noqa: E402

OUT = ROOT / "presentacion" / "imagenes"
TMM_RESULTS = ROOT / "results" / "tmm"

COLOR_A = "#cfe6fb"
COLOR_B = "#f4a261"
EDGE = "#333333"
BLUE = "#2c4a8c"
RED = "#c8412c"
GREEN = "#2e8b57"
PURPLE = "#7b4fa0"

SPEC_FIG1 = load_spec("figure_01.yaml")
SPEC_FIG3 = load_spec("figure_03.yaml")
TMM_FIG1 = TMMSolver(SPEC_FIG1)


def load_summary(name: str) -> dict:
    """``results/tmm/<name>/data/summary.json`` (ejecute antes ``make tmm``)."""
    return json.loads((TMM_RESULTS / name / "data" / "summary.json").read_text(encoding="utf-8"))


def style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.size": 12,
            "axes.labelsize": 12,
            "axes.titlesize": 13,
            "legend.fontsize": 10,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "savefig.dpi": 200,
            "savefig.bbox": "tight",
            "mathtext.fontset": "dejavusans",
        }
    )


def save(fig: plt.Figure, name: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{name}.png")
    fig.savefig(OUT / f"{name}.pdf")
    plt.close(fig)
    print("  ", name)


def draw_word(
    ax, word: str, x0: float, y0: float, w: float = 1.0, h: float = 1.0, labels: bool = True, fontsize: float = 11
) -> None:
    for i, layer in enumerate(word):
        color = COLOR_A if layer == "A" else COLOR_B
        ax.add_patch(Rectangle((x0 + i * w, y0), w, h, facecolor=color, edgecolor=EDGE, lw=1.0))
        if labels:
            ax.text(x0 + (i + 0.5) * w, y0 + h / 2, layer, ha="center", va="center", fontsize=fontsize, weight="bold")


def fig_superred() -> None:
    """Superred: celda S4 = ABAAB repetida, onda TE incidente a ángulo θ."""
    fig, ax = plt.subplots(figsize=(11, 4.2))
    word = sequence(4)
    x0 = 3.0
    h = 3.0
    for rep in range(2):
        draw_word(ax, word, x0 + rep * len(word), 0.0, h=h)
    ax.text(x0 - 0.5, h / 2, "···", ha="center", va="center", fontsize=20)
    ax.text(x0 + 2 * len(word) + 0.5, h / 2, "···", ha="center", va="center", fontsize=20)

    # Llave de la celda
    xa, xb = x0, x0 + len(word)
    ax.annotate("", xy=(xa, h + 0.35), xytext=(xb, h + 0.35), arrowprops=dict(arrowstyle="<->", lw=1.4))
    ax.text((xa + xb) / 2, h + 0.55, "celda $S_4 = ABAAB$\nperiodo $L_4 = 3a + 2b$", ha="center", va="bottom", fontsize=11.5)
    ax.annotate("", xy=(xb, h + 0.35), xytext=(xb + len(word), h + 0.35), arrowprops=dict(arrowstyle="<->", lw=1.0, color="gray"))
    ax.text(xb + len(word) / 2, h + 0.55, "se repite periódicamente", ha="center", va="bottom", fontsize=11, color="gray")

    # Onda incidente
    theta = np.deg2rad(35)
    start = np.array([0.0, 0.35])
    length = 2.6
    end = start + length * np.array([np.cos(theta), np.sin(theta)])
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=22, lw=2.2, color=RED))
    ax.text(end[0] - 0.25, end[1] + 0.2, r"$\mathbf{k}$", color=RED, fontsize=15)
    ax.plot([start[0], start[0] + 2.8], [start[1], start[1]], ls="--", color="gray", lw=1)
    ax.add_patch(Arc(start, 2.0, 2.0, theta1=0, theta2=np.rad2deg(theta), color=RED, lw=1.4))
    ax.text(start[0] + 1.1, start[1] + 0.3, r"$\theta$", color=RED, fontsize=15)
    ax.text(0.2, 2.75, "onda TE\nincidente", color=RED, fontsize=11, ha="left")

    # Ejes
    ox, oy = x0 + 3.0, -1.1
    ax.add_patch(FancyArrowPatch((ox, oy), (ox + 1.8, oy), arrowstyle="-|>", mutation_scale=15, lw=1.3))
    ax.text(ox + 1.95, oy, r"$z$ (dirección de crecimiento)", va="center", fontsize=11)
    ax.add_patch(FancyArrowPatch((ox, oy), (ox, oy + 0.9), arrowstyle="-|>", mutation_scale=15, lw=1.3))
    ax.text(ox - 0.4, oy + 0.7, r"$x$", fontsize=11)

    # Leyenda de materiales
    lx = x0 + 2 * len(word) + 1.4
    ax.add_patch(Rectangle((lx, 1.9), 0.7, 0.7, facecolor=COLOR_A, edgecolor=EDGE))
    ax.text(lx + 0.9, 2.25, "A: aire\n" + r"$\varepsilon_A=\mu_A=1$, $a=12$ mm", va="center", fontsize=10.5)
    ax.add_patch(Rectangle((lx, 0.6), 0.7, 0.7, facecolor=COLOR_B, edgecolor=EDGE))
    ax.text(lx + 0.9, 0.95, "B: metamaterial Drude\n" + r"$\varepsilon_B(\omega)$, $\mu_B(\omega)$, $b=12$ mm", va="center", fontsize=10.5)

    ax.set_xlim(-0.3, lx + 5.2)
    ax.set_ylim(-1.5, h + 1.8)
    ax.set_aspect("equal")
    ax.axis("off")
    save(fig, "superred_geometria")


def fig_fibonacci() -> None:
    """Construcción S_m = S_{m-1} S_{m-2}."""
    fig, ax = plt.subplots(figsize=(8.0, 6.0))
    x0, w, n_col = 1.1, 0.8, 1.1 + 13 * 0.8 + 0.4
    for m in range(7):
        y = -m * 1.25
        word = sequence(m)
        draw_word(ax, word, x0, y, w=w, h=0.95, fontsize=13)
        ax.text(0.0, y + 0.47, rf"$S_{m}$", va="center", fontsize=17, weight="bold")
        if m >= 2:
            xs = x0 + len(sequence(m - 1)) * w
            ax.plot([xs, xs], [y - 0.15, y + 1.1], color=RED, lw=3)
        ax.text(n_col, y + 0.47, rf"$N_B={n_layers_b(m)}$", va="center", fontsize=14)
    ax.text(x0, 1.3, "celda $S_m$ (la línea roja separa $S_{m-1}$ de $S_{m-2}$)", fontsize=13, color="#555555")
    ax.text(n_col, 1.3, "capas B", fontsize=13, color="#555555")
    ax.set_xlim(-0.1, n_col + 2.0)
    ax.set_ylim(-6 * 1.25 - 0.3, 1.9)
    ax.axis("off")
    save(fig, "fibonacci_construccion")


def fig_drude() -> None:
    """ε_B(ν) y μ_B(ν) para los dos juegos de parámetros del paper."""
    nu = np.linspace(0.25, 5.0, 2000)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    cases = [
        (axes[0], 3.0, 3.0, r"Figs. 1–2: $\nu_e=\nu_m=3$ GHz"),
        (axes[1], 3.0, 1.0, r"Figs. 3–6: $\nu_e=3$ GHz, $\nu_m=1$ GHz"),
    ]
    for ax, nu_e, nu_m, title in cases:
        eps = 1 - (nu_e / nu) ** 2
        mu = 1 - (nu_m / nu) ** 2
        ax.axvspan(nu[0], min(nu_e, nu_m), color=PURPLE, alpha=0.12)
        if nu_m < nu_e:
            ax.text(0.62, -3.05, "índice\nnegativo", ha="center", va="center", fontsize=9.5, color=PURPLE)
            ax.axvspan(nu_m, nu_e, color="gray", alpha=0.10)
            ax.text(2.35, -1.7, r"$\varepsilon<0,\ \mu>0$" + "\n(un solo signo\nnegativo)", ha="center", va="top", fontsize=10, color="#444444")
        else:
            ax.text(0.5 * (nu[0] + min(nu_e, nu_m)), 0.75, "índice\nnegativo\n" + r"($\varepsilon<0,\ \mu<0$)", ha="center", va="top", fontsize=10, color=PURPLE)
        ax.plot(nu, eps, color=BLUE, lw=2.4, label=r"$\varepsilon_B(\nu)=1-\nu_e^2/\nu^2$")
        ax.plot(nu, mu, color=RED, lw=2.0, ls="--", label=r"$\mu_B(\nu)=1-\nu_m^2/\nu^2$")
        ax.axhline(0, color="black", lw=0.8)
        ax.axvline(nu_m, color=RED, lw=1, ls=":")
        ax.axvline(nu_e, color=BLUE, lw=1, ls=":")
        ax.text(nu_m, -3.85, r" $\nu_m$", color=RED, fontsize=12)
        ax.text(nu_e, -3.45, r" $\nu_e$", color=BLUE, fontsize=12)
        ax.set_xlim(nu[0], nu[-1])
        ax.set_ylim(-4, 1.0)
        ax.set_xlabel("frecuencia $\\nu$ (GHz)")
        ax.set_title(title)
    axes[0].set_ylabel("valor de $\\varepsilon_B$, $\\mu_B$")
    axes[1].legend(loc="lower right", frameon=True)
    fig.tight_layout()
    save(fig, "drude_epsilon_mu")


def _wave_panel(ax, theta_deg: float, title: str) -> None:
    for xi in (0.0, 1.6, 3.2, 4.8):
        ax.add_patch(Rectangle((xi, -2.0), 0.8, 4.0, facecolor=COLOR_B, edgecolor=EDGE, alpha=0.55, lw=0.8))
    theta = np.deg2rad(theta_deg)
    o = np.array([-2.4, 0.0])
    kdir = np.array([np.cos(theta), np.sin(theta)])
    hdir = np.array([-np.sin(theta), np.cos(theta)])
    ax.add_patch(FancyArrowPatch(o, o + 2.2 * kdir, arrowstyle="-|>", mutation_scale=20, lw=2.4, color=RED))
    ax.text(*(o + 2.3 * kdir + np.array([0.05, 0.1])), r"$\mathbf{k}$", color=RED, fontsize=15)
    hvec = 1.5 * hdir
    ax.add_patch(FancyArrowPatch(o, o + hvec, arrowstyle="-|>", mutation_scale=18, lw=2.0, color=BLUE))
    ax.text(*(o + hvec + np.array([-0.45, 0.1])), r"$\mathbf{H}$", color=BLUE, fontsize=15)
    if theta_deg > 0:
        ax.add_patch(FancyArrowPatch(o, o + np.array([hvec[0], 0]), arrowstyle="-|>", mutation_scale=14, lw=1.6, color=GREEN))
        ax.text(o[0] + hvec[0] - 0.05, -0.45, r"$H_z$", color=GREEN, fontsize=14, ha="center")
        ax.plot([o[0] + hvec[0], o[0] + hvec[0]], [0, hvec[1]], ls=":", color=GREEN)
    ax.plot(*o, marker="o", ms=16, mfc="white", mec="black", mew=1.5)
    ax.plot(*o, marker=".", ms=8, color="black")
    ax.text(o[0] + 0.3, o[1] - 0.45, r"$\mathbf{E}\parallel y$" + "\n(sale del plano)", fontsize=10, ha="center", va="top")
    ax.add_patch(FancyArrowPatch((-4.1, -2.1), (-3.2, -2.1), arrowstyle="-|>", mutation_scale=12, lw=1.2))
    ax.text(-3.1, -2.1, "$z$", va="center")
    ax.add_patch(FancyArrowPatch((-4.1, -2.1), (-4.1, -1.2), arrowstyle="-|>", mutation_scale=12, lw=1.2))
    ax.text(-4.1, -1.05, "$x$", ha="center")
    ax.set_title(title, fontsize=13)
    ax.set_xlim(-4.2, 5.8)
    ax.set_ylim(-2.3, 2.4)
    ax.set_aspect("equal")
    ax.axis("off")


def fig_polarizacion() -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    _wave_panel(axes[0], 35, "Incidencia oblicua ($\\theta\\neq0$)\n$H_z\\neq0$ → excita el plasmón magnético")
    _wave_panel(axes[1], 0, "Incidencia normal ($\\theta=0$)\n$H_z=0$ → no hay acoplamiento")
    fig.tight_layout()
    save(fig, "polarizacion_TE")


def fig_osciladores() -> None:
    """Cada capa B es un oscilador: N_B capas → N_B modos (datos reales de la Fig. 2)."""
    counts = {entry["m"]: entry for entry in load_summary("figure_02")["mode_counts"]}
    fig, (ax_w, ax_l) = plt.subplots(1, 2, figsize=(11, 4.4), gridspec_kw={"width_ratios": [1.25, 1]})
    for row, m in enumerate((3, 4, 5, 6)):
        y = -row * 1.2
        word = sequence(m)
        draw_word(ax_w, word, 1.6, y, w=0.55, h=0.8, labels=True)
        ax_w.text(0.0, y + 0.4, rf"$S_{m}$", fontsize=14, va="center", weight="bold")
        ax_w.text(1.6 + 13 * 0.55 + 0.3, y + 0.4, rf"$N_B={n_layers_b(m)}$", fontsize=12, va="center", color="#b5541c")
    ax_w.set_xlim(-0.2, 1.6 + 13 * 0.55 + 2.2)
    ax_w.set_ylim(-3 * 1.2 - 0.3, 1.2)
    ax_w.set_title("capas B (naranja) = osciladores magnéticos", fontsize=12)
    ax_w.axis("off")

    for i, m in enumerate((3, 4, 5, 6)):
        for interval in counts[m]["intervals"]:
            ax_l.add_patch(
                Rectangle((i - 0.3, interval["nu_min"]), 0.6, interval["width"], facecolor=BLUE, edgecolor=BLUE, alpha=0.85)
            )
        n = counts[m]["detected"]
        ax_l.text(i, 2.93, f"{n} modo" + ("s" if n > 1 else ""), ha="center", fontsize=10.5)
    ax_l.axhline(3.0, color=RED, ls="--", lw=1.4)
    ax_l.text(-0.55, 3.02, r"$\nu_m = 3$ GHz", color=RED, va="bottom", fontsize=11)
    ax_l.set_xticks(range(4), [f"m={m}" for m in (3, 4, 5, 6)])
    ax_l.set_xlim(-0.6, 3.6)
    ax_l.set_ylim(2.0, 3.12)
    ax_l.set_ylabel("frecuencia (GHz)")
    ax_l.set_title(r"subbandas bajo $\nu_m$ ($\theta=\pi/3$, datos de la Fig. 2)", fontsize=12)
    fig.tight_layout()
    save(fig, "osciladores_modos")


def fig_semitraza() -> None:
    """R_m(ν) → bandas (|R|≤1) y gaps; relación cos(kL)=R."""
    nu = np.linspace(0.15, 5.0, 8000)
    scan = TMM_FIG1.scan(3, np.pi / 6, nu)
    r = np.real(scan.r)
    r_plot = np.where(np.abs(r) < 3.5, r, np.nan)
    fig, (ax_r, ax_k) = plt.subplots(1, 2, figsize=(11, 4.6), sharey=True, gridspec_kw={"width_ratios": [1, 1]})
    allowed = scan.allowed
    edges = np.flatnonzero(np.diff(allowed.astype(int)) != 0)
    starts = np.r_[0, edges + 1]
    ends = np.r_[edges, allowed.size - 1]
    for s, e in zip(starts, ends, strict=True):
        color = GREEN if allowed[s] else "gray"
        alpha = 0.13 if allowed[s] else 0.18
        for ax in (ax_r, ax_k):
            ax.axhspan(nu[s], nu[e], color=color, alpha=alpha, lw=0)
    ax_r.plot(r_plot, nu, color="black", lw=1.4)
    ax_r.axvline(1, color=RED, ls="--", lw=1.2)
    ax_r.axvline(-1, color=RED, ls="--", lw=1.2)
    ax_r.set_xlim(-3.5, 3.5)
    ax_r.set_xlabel(r"semitraza $R_m(\nu)=\frac{1}{2}\mathrm{Tr}\,T_m$")
    ax_r.set_ylabel("frecuencia (GHz)")
    ax_r.set_title(r"paso 1: calcular $R_m$ con matrices $2\times2$", fontsize=12)
    ax_r.text(0, 4.75, r"$|R_m|\leq 1$", ha="center", color=RED, fontsize=12, bbox=dict(fc="white", ec="none", alpha=0.8))

    plot_closed(ax_k, scan, color=BLUE, linestyle="-", linewidth=1.6)
    ax_k.set_xlim(-1, 1)
    ax_k.set_xlabel(r"$kL_m/\pi = \pm\arccos(R_m)/\pi$")
    ax_k.set_title(r"paso 2: $\cos(kL_m)=R_m$ → bandas $\nu(k)$", fontsize=12)
    ax_k.set_ylim(nu[0], nu[-1])
    ax_k.text(1.03, 0.62, "verde: banda\n(la onda se propaga)\n\ngris: gap\n(la onda se atenúa)", transform=ax_k.transAxes, fontsize=10.5, va="center")
    fig.suptitle(r"Ejemplo: $S_3$, $\nu_e=\nu_m=3$ GHz, $\theta=\pi/6$", fontsize=12, y=1.0)
    fig.tight_layout()
    save(fig, "semitraza_bandas")


def plot_closed(ax, scan, **style) -> None:
    edges = close_band_edges(TMM_FIG1, scan)
    plot_dispersion_branches(ax, scan, edge_nu_ghz=edges.nu_ghz, edge_k_lm_over_pi=edges.k_lm_over_pi, **style)


def fig_leer_dispersion() -> None:
    nu = np.linspace(0.15, 5.0, 10000)
    scan = TMM_FIG1.scan(4, np.pi / 6, nu)
    fig, ax = plt.subplots(figsize=(6.4, 5.6))
    plot_closed(ax, scan, color=BLUE, linestyle="-", linewidth=1.8)

    def point_on_band(nu_lo: float, nu_hi: float, k_target: float) -> tuple[float, float]:
        sel = scan.allowed & (nu > nu_lo) & (nu < nu_hi)
        idx = np.flatnonzero(sel)[np.argmin(np.abs(scan.k_lm_over_pi[sel] - k_target))]
        return float(scan.k_lm_over_pi[idx]), float(nu[idx])

    ax.set_xlim(-1, 1)
    ax.set_ylim(0, 5)
    ax.set_xlabel(r"número de onda de Bloch  $kL_m/\pi$")
    ax.set_ylabel("frecuencia (GHz)")
    gaps = []
    allowed = scan.allowed
    edges = np.flatnonzero(np.diff(allowed.astype(int)) != 0)
    starts = np.r_[0, edges + 1]
    ends = np.r_[edges, allowed.size - 1]
    for s, e in zip(starts, ends, strict=True):
        if not allowed[s] and nu[e] - nu[s] > 0.08:
            gaps.append((nu[s], nu[e]))
            ax.axhspan(nu[s], nu[e], color="gray", alpha=0.2, lw=0)
    top_gap = max(gaps, key=lambda g: g[0])
    arrow = dict(arrowstyle="->", lw=1.2)
    ax.annotate(
        "gap: no hay curva\n(frecuencia prohibida)",
        xy=(0.8, 0.5 * (top_gap[0] + top_gap[1])),
        xytext=(1.05, 4.55),
        fontsize=12.5,
        annotation_clip=False,
        arrowprops=arrow,
    )
    ax.annotate(
        "banda: cada punto es un\nmodo que se propaga",
        xy=point_on_band(3.3, 4.0, 0.55),
        xytext=(1.05, 3.55),
        fontsize=12.5,
        annotation_clip=False,
        arrowprops=arrow,
    )
    ax.annotate(
        "banda casi plana cerca de $\\nu_m$:\nmodo plasmon-polaritón\n(casi no se propaga)",
        xy=point_on_band(2.88, 2.97, 0.75),
        xytext=(1.05, 2.2),
        fontsize=12.5,
        annotation_clip=False,
        arrowprops=arrow,
    )
    ax.set_title(r"Ejemplo: $S_4$, $\nu_e=\nu_m=3$ GHz, $\theta=\pi/6$", fontsize=12)
    save(fig, "leer_dispersion")


def fig_gap_n0() -> None:
    nu = np.linspace(0.3, 3.0, 2000)
    omega = omega_from_nu_ghz(nu)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)

    ax = axes[0]
    for m, color in ((3, RED), (4, BLUE)):
        eps, mu = average_epsilon_mu(SPEC_FIG1, m, omega)
        ax.plot(nu, np.real(eps), color=color, lw=2.2, label=rf"$\langle\varepsilon\rangle_{m}=\langle\mu\rangle_{m}$")
        zero = 3.0 * np.sqrt(n_layers_b(m) / len(sequence(m)))
        ax.plot(zero, 0, "o", color=color, ms=7)
        ax.text(zero, 0.25 + 0.25 * (m - 3), f"{zero:.3f} GHz", color=color, ha="center", fontsize=10)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_title(r"Figs. 1–2 ($\nu_e=\nu_m$): ambos promedios se anulan juntos" + "\n" + r"→ a $\theta=0$ el gap $\langle n\rangle=0$ está cerrado", fontsize=11.5)
    ax.set_xlabel("frecuencia (GHz)")
    ax.set_ylabel("promedio ponderado por espesor")
    ax.legend(loc="lower right")

    ax = axes[1]
    eps, mu = average_epsilon_mu(SPEC_FIG3, 3, omega)
    z_mu = 1.0 * np.sqrt(1 / 3)
    z_eps = 3.0 * np.sqrt(1 / 3)
    ax.axvspan(z_mu, z_eps, color="gray", alpha=0.18)
    ax.plot(nu, np.real(eps), color=BLUE, lw=2.2, label=r"$\langle\varepsilon\rangle_3$")
    ax.plot(nu, np.real(mu), color=RED, lw=2.2, ls="--", label=r"$\langle\mu\rangle_3$")
    ax.axhline(0, color="black", lw=0.8)
    ax.axvline(1.0, color=PURPLE, lw=2)
    ax.text(1.02, -2.6, r"$\nu_m=1$ GHz" + "\ncae dentro\ndel gap", color=PURPLE, fontsize=10.5)
    ax.text(0.86, -0.45, r"gap $\langle n\rangle=0$", ha="center", fontsize=11, bbox=dict(fc="white", ec="none", alpha=0.9))
    ax.set_title(r"Figs. 3–6 ($\nu_m=1$ GHz): los ceros se separan" + "\n" + r"→ se abre un gap y el plasmón queda dentro", fontsize=11.5)
    ax.set_xlabel("frecuencia (GHz)")
    ax.legend(loc="lower right")
    for a in axes:
        a.set_ylim(-3.2, 1.0)
        a.set_xlim(nu[0], nu[-1])
    fig.tight_layout()
    save(fig, "gap_indice_nulo")


def fig_convergencia() -> None:
    path = TMM_RESULTS / "convergence" / "convergence_figure04.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    n = np.array([d["n_points"] for d in data])
    w = np.array([d["width_ghz"] for d in data]) * 1e3
    fig, ax = plt.subplots(figsize=(6.0, 4.6))
    ax.semilogx(n, w, "o-", color=BLUE, lw=2, ms=7)
    ax.axhline(w[-1], color="gray", ls="--", lw=1)
    ax.fill_between([n[0] * 0.8, n[-1] * 1.25], w[-1] * 0.998, w[-1] * 1.002, color=GREEN, alpha=0.18, label=r"$\pm 0.2\%$")
    ax.set_xlim(n[0] * 0.8, n[-1] * 1.25)
    ax.set_xlabel("puntos en el barrido de frecuencia")
    ax.set_ylabel(r"ancho de banda $\Delta\nu$ (MHz)")
    ax.set_title(r"Convergencia: $m=3$, $\theta=\pi/3$ (Fig. 4)", fontsize=12)
    ax.legend(loc="lower right")
    save(fig, "convergencia")


def fig_conteo_modos() -> None:
    """Subbandas detectadas en la Fig. 5 (TMM) frente a F_{m−2}."""
    orders = list(load_figure("figure_05").config.fibonacci_orders)
    expected = [n_layers_b(m) for m in orders]
    rows = load_summary("figure_05")["intervals"]
    detected = {label: {row["m"]: row["n"] for row in rows[label]} for label in rows}
    det_12 = [detected["pi/12"][m] for m in orders]
    det_3 = [detected["pi/3"][m] for m in orders]
    x = np.arange(len(orders))
    fig, ax = plt.subplots(figsize=(8.0, 4.2))
    ax.bar(x - 0.27, expected, 0.27, color="#bbbbbb", label=r"esperado $F_{m-2}$")
    ax.bar(x, det_12, 0.27, color=BLUE, label=r"detectado, $\theta=\pi/12$")
    ax.bar(x + 0.27, det_3, 0.27, color=RED, label=r"detectado, $\theta=\pi/3$")
    for xi, e in zip(x, expected, strict=True):
        ax.text(xi, e + 0.25, str(e), ha="center", fontsize=11, weight="bold")
    ax.set_xticks(x, [f"m={m}" for m in orders])
    ax.set_ylabel("número de subbandas")
    ax.set_title("Número de modos plasmon-polaritón vs orden de Fibonacci", fontsize=12)
    ax.legend(loc="upper left")
    save(fig, "conteo_modos")


def main() -> None:
    style()
    print("Generando imágenes en", OUT)
    fig_superred()
    fig_fibonacci()
    fig_drude()
    fig_polarizacion()
    fig_osciladores()
    fig_semitraza()
    fig_leer_dispersion()
    fig_gap_n0()
    fig_convergencia()
    fig_conteo_modos()


if __name__ == "__main__":
    main()
