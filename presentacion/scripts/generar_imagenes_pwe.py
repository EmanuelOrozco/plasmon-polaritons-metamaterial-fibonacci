#!/usr/bin/env python3
"""Genera las imágenes conceptuales de la presentación del PWE en presentacion/imagenes_pwe/."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from fibonacci_photonics.config import load_spec  # noqa: E402
from fibonacci_photonics.physics.electromagnetics import Polarization  # noqa: E402
from fibonacci_photonics.physics.units import omega_from_nu_ghz  # noqa: E402
from fibonacci_photonics.solvers import TMMSolver  # noqa: E402
from fibonacci_photonics.solvers.pwe.bloch import (  # noqa: E402
    BRILLOUIN_WINDOW,
    bloch_eigenvalues,
    n_max_for,
    select_bloch_wavevector,
)
from fibonacci_photonics.solvers.pwe.eigenfrequency import mode_profile, nondispersive_bands  # noqa: E402
from fibonacci_photonics.solvers.pwe.fourier import bilayer_fourier, cell_fourier, synthesize_profile  # noqa: E402

OUT = ROOT / "presentacion" / "imagenes_pwe"

COLOR_A = "#cfe6fb"
COLOR_B = "#f4a261"
BLUE = "#2c4a8c"
RED = "#c8412c"
GREEN = "#2e8b57"
GRAY = "#888888"

BOOK_FRACTIONS = (0.2, 0.8)
BOOK_EPS = (1.0, 9.0)


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


def shade_layers(ax, periods: int, fraction_1: float, color: str = COLOR_B) -> None:
    """Sombrea la capa 2 de la bicapa (z/a en [f1, 1) de cada celda)."""
    for j in range(-periods, periods):
        ax.axvspan(j + fraction_1, j + 1, color=color, alpha=0.35, lw=0)


def bilayer_semitrace(f: np.ndarray, eps: tuple[float, float] = BOOK_EPS) -> np.ndarray:
    """cos(ka) de la bicapa a incidencia normal; f = ωa/2πc."""
    q1 = 2 * np.pi * f * np.sqrt(eps[0])
    q2 = 2 * np.pi * f * np.sqrt(eps[1])
    d1, d2 = BOOK_FRACTIONS
    ratio = 0.5 * (q1 * eps[1] / (q2 * eps[0]) + q2 * eps[0] / (q1 * eps[1]))
    return np.cos(q1 * d1) * np.cos(q2 * d2) - ratio * np.sin(q1 * d1) * np.sin(q2 * d2)


def fig_ondas_planas() -> None:
    """Un perfil escalonado como suma de ondas planas."""
    cell = bilayer_fourier(*BOOK_FRACTIONS, 60)
    coeffs = cell.profile_coefficients(1 / BOOK_EPS[0], 1 / BOOK_EPS[1])
    z = np.linspace(-1.0, 1.0, 2001)
    exact = np.where(np.mod(z, 1.0) < BOOK_FRACTIONS[0], 1 / BOOK_EPS[0], 1 / BOOK_EPS[1])

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 3.8))
    ax = axes[0]
    for k, n in enumerate(range(0, 4)):
        mask = np.abs(cell.orders) == n
        term = synthesize_profile(coeffs[mask], cell.orders[mask], 1.0, z).real
        ax.plot(z, term - 0.55 * k, lw=1.4, color=[BLUE, GREEN, RED, "#7b4fa0"][k])
        ax.text(1.03, term.mean() - 0.55 * k, rf"$|n|={n}$", va="center", fontsize=10)
    ax.set_yticks([])
    ax.set_xlabel(r"$z/a$")
    ax.set_title(r"Ondas planas $f_n\,e^{iG_nz}+f_{-n}\,e^{-iG_nz}$", fontsize=12)
    ax.spines["left"].set_visible(False)

    ax = axes[1]
    shade_layers(ax, 1, BOOK_FRACTIONS[0])
    ax.plot(z, exact, color="black", lw=1.6, label="perfil exacto")
    for n_max, color in ((1, GREEN), (3, BLUE), (25, RED)):
        mask = np.abs(cell.orders) <= n_max
        ax.plot(
            z,
            synthesize_profile(coeffs[mask], cell.orders[mask], 1.0, z).real,
            color=color,
            lw=1.1,
            label=f"N = {2 * n_max + 1}",
        )
    ax.set_xlabel(r"$z/a$")
    ax.set_ylabel(r"$1/\varepsilon(z)$")
    ax.set_title("La suma reconstruye el perfil (Gibbs en los saltos)", fontsize=12)
    ax.legend(loc="center right", fontsize=9)
    fig.tight_layout()
    save(fig, "ondas_planas_suma")


def fig_bloch() -> None:
    """Modo de Bloch u(z) = h(z) e^{ikz} del cristal del libro."""
    cell = bilayer_fourier(*BOOK_FRACTIONS, 30)
    k = 0.4 * np.pi
    freqs, modes = nondispersive_bands(cell, BOOK_EPS, (1, 1), [k], rule="inverse", n_bands=2, return_modes=True)
    coeffs = modes[0][:, 1]
    z = np.linspace(0.0, 4.0, 4001)
    u = mode_profile(cell, coeffs, k, z)
    u = u / np.abs(u).max()
    h = u * np.exp(-1j * k * z)
    phase = np.exp(-1j * np.angle(h[0]))
    u, h = u * phase, h * phase

    fig, axes = plt.subplots(2, 1, figsize=(7.6, 5.4), sharex=True)
    for ax in axes:
        for j in range(4):
            ax.axvspan(j + BOOK_FRACTIONS[0], j + 1, color=COLOR_B, alpha=0.35, lw=0)
            ax.axvline(j, color=GRAY, lw=0.6, ls=":")
    axes[0].plot(z, h.real, color=GREEN, lw=1.5, label=r"$\mathrm{Re}\,h(z)$: periódica, se repite en cada celda")
    axes[0].legend(loc="upper center", fontsize=9, ncol=1, bbox_to_anchor=(0.5, 1.22), frameon=False)
    axes[0].set_yticks([])
    axes[1].plot(z, u.real, color=BLUE, lw=1.5, label=r"$\mathrm{Re}\,u(z)=\mathrm{Re}\,[h(z)\,e^{ikz}]$")
    axes[1].plot(z, np.cos(k * z) * np.abs(h).max(), color=RED, lw=1.0, ls="--", label=r"$\cos(kz)$")
    axes[1].legend(loc="upper center", fontsize=9, ncol=2, bbox_to_anchor=(0.5, 1.2), frameon=False)
    axes[1].set_yticks([])
    axes[1].set_xlabel(r"$z/a$  (sombreado: capa $\varepsilon_2=9$)")
    fig.suptitle(rf"Modo de Bloch, banda 2, $ka=0.4\pi$, $\omega a/2\pi c={freqs[0, 1]:.3f}$", fontsize=12)
    fig.tight_layout()
    save(fig, "bloch_onda")


def fig_plegado() -> None:
    """Zona extendida, plegado a la primera zona y apertura de gaps."""
    eps_avg = BOOK_FRACTIONS[0] * BOOK_EPS[0] + BOOK_FRACTIONS[1] * BOOK_EPS[1]
    n_eff = np.sqrt(eps_avg)
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), sharey=True)

    ax = axes[0]
    k = np.linspace(-3 * np.pi, 3 * np.pi, 600)
    ax.plot(k / np.pi, np.abs(k) / (2 * np.pi * n_eff), color=BLUE, lw=1.6)
    for j in (-3, -1, 1, 3):
        ax.axvline(j, color=GRAY, lw=0.7, ls=":")
    ax.axvspan(-1, 1, color=COLOR_A, alpha=0.6, lw=0)
    ax.set_xlabel(r"$ka/\pi$")
    ax.set_ylabel(r"$\omega a/2\pi c$")
    ax.set_title(r"Medio homogéneo: $\omega=ck/\bar n$", fontsize=12)
    ax.text(0, 0.9, "1.ª zona", ha="center", fontsize=10)

    ax = axes[1]
    kz = np.linspace(-np.pi, np.pi, 301)
    for g in range(-3, 4):
        ax.plot(kz / np.pi, np.abs(kz + 2 * np.pi * g) / (2 * np.pi * n_eff), color=GRAY, lw=0.9, ls="--")
    cell = bilayer_fourier(*BOOK_FRACTIONS, 40)
    bands = nondispersive_bands(cell, BOOK_EPS, (1, 1), kz, rule="inverse", n_bands=5)
    for j in range(bands.shape[1]):
        ax.plot(kz / np.pi, bands[:, j], color=RED, lw=1.6)
    for j in range(bands.shape[1] - 1):
        lo, hi = bands[:, j].max(), bands[:, j + 1].min()
        if hi > lo:
            ax.axhspan(lo, hi, color=COLOR_B, alpha=0.35, lw=0)
    ax.set_xlim(-1, 1)
    ax.set_ylim(0, 1.0)
    ax.set_xlabel(r"$ka/\pi$")
    ax.set_title("Plegado (gris) y cristal real (rojo): gaps", fontsize=12)
    fig.tight_layout()
    save(fig, "plegado_bandas")


def fig_toeplitz() -> None:
    """Estructura de la matriz de Toeplitz [[μ]] de la celda S4."""
    spec = load_spec("figure_03.yaml")
    n_max = n_max_for(4, 2)
    cell = cell_fourier(4, spec.thickness_a, spec.thickness_b, n_max)
    omega = float(omega_from_nu_ghz(0.9))
    mu_b = complex(spec.medium_b.mu(omega))
    t = cell.toeplitz(1.0, mu_b)
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.6), gridspec_kw={"width_ratios": [1.0, 1.2]})
    ax = axes[0]
    im = ax.imshow(np.log10(np.abs(t) + 1e-6), cmap="viridis", vmin=-3, vmax=0)
    ax.set_title(rf"$\log_{{10}}|[\![\mu]\!]_{{nn'}}|$, $S_4$, $\nu=0.9$ GHz, $N={cell.size}$", fontsize=12)
    ticks = [0, n_max, 2 * n_max]
    ax.set_xticks(ticks, [f"{-n_max}", "0", f"{n_max}"])
    ax.set_yticks(ticks, [f"{-n_max}", "0", f"{n_max}"])
    ax.set_xlabel(r"$n'$")
    ax.set_ylabel(r"$n$")
    fig.colorbar(im, ax=ax, fraction=0.046)
    ax = axes[1]
    orders = cell.orders
    coeffs = cell.profile_coefficients(1.0, mu_b)
    pos = (orders >= 0) & (np.abs(coeffs) > 1e-10)
    ax.semilogy(orders[pos], np.abs(coeffs[pos]), "o", color=BLUE, ms=4, label=r"$|\mu_n|$")
    nz = np.arange(1, orders.max() + 1)
    ax.semilogy(nz, 0.25 / nz, color=RED, lw=1.0, ls="--", label=r"$\propto 1/n$")
    ax.set_xlabel(r"orden $n$")
    ax.set_title(r"Cada diagonal repite un coeficiente $\mu_{n-n'}$", fontsize=12)
    ax.legend()
    fig.tight_layout()
    save(fig, "toeplitz")


def fig_omega_k_vs_k_omega() -> None:
    """Dos formas de cortar el diagrama de bandas: ω(k) y k(ω)."""
    f = np.linspace(1e-4, 0.75, 6000)
    r = bilayer_semitrace(f)
    allowed = np.abs(r) <= 1
    k = np.where(allowed, np.arccos(np.clip(r, -1, 1)) / np.pi, np.nan)
    cell = bilayer_fourier(*BOOK_FRACTIONS, 40)
    k_cut = 0.55
    omegas = nondispersive_bands(cell, BOOK_EPS, (1, 1), [k_cut * np.pi], rule="inverse", n_bands=4)[0]

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), sharey=True)
    for ax in axes:
        ax.plot(k, f, color=BLUE, lw=1.6)
        ax.plot(-k, f, color=BLUE, lw=1.6)
        ax.fill_betweenx(f, -1, 1, where=~allowed, color=COLOR_B, alpha=0.35, lw=0)
        ax.set_xlim(-1, 1)
        ax.set_ylim(0, 0.75)
        ax.set_xlabel(r"$ka/\pi$")
    ax = axes[0]
    ax.axvline(k_cut, color=RED, lw=1.2)
    ax.plot([k_cut] * len(omegas), omegas, "o", color=RED, ms=7)
    ax.set_ylabel(r"$\omega a/2\pi c$")
    ax.set_title(r"$\omega(k)$: fijo $k$, busco los $\omega$ (libro)", fontsize=12)
    ax = axes[1]
    for f_cut, label in ((0.27, "banda: $k$ real"), (0.37, "gap: $k$ complejo")):
        ax.axhline(f_cut, color=GREEN, lw=1.2)
        rr = bilayer_semitrace(np.array([f_cut]))[0]
        if abs(rr) <= 1:
            kk = np.arccos(rr) / np.pi
            ax.plot([kk, -kk], [f_cut, f_cut], "o", color=GREEN, ms=7)
        ax.text(-0.95, f_cut + 0.012, label, fontsize=10, color=GREEN)
    ax.set_title(r"$k(\omega)$: fijo $\omega$, busco $k$ (proyecto)", fontsize=12)
    fig.tight_layout()
    save(fig, "omega_k_vs_k_omega")


def fig_regla_li() -> None:
    """1/χ y u' saltan en las interfaces; su producto D = u'/χ es continuo."""
    eps = (1.0, 9.0)
    d = BOOK_FRACTIONS
    f = 0.3
    kappa = 2 * np.pi * f
    z_all, u_all, du_all, chi_all = [], [], [], []
    state = np.array([1.0, 0.0])
    z0 = 0.0
    for _ in range(2):
        for eps_j, d_j in zip(eps, d, strict=True):
            chi = eps_j
            qz = kappa * np.sqrt(eps_j)
            s = np.linspace(0.0, d_j, 300)
            u = state[0] * np.cos(qz * s) + state[1] * chi * np.sin(qz * s) / qz
            flux = -state[0] * qz * np.sin(qz * s) / chi + state[1] * np.cos(qz * s)
            z_all.append(z0 + s)
            u_all.append(u)
            du_all.append(chi * flux)
            chi_all.append(np.full_like(s, chi))
            state = np.array([u[-1], flux[-1]])
            z0 += d_j
    z = np.concatenate(z_all)
    du = np.concatenate(du_all)
    chi = np.concatenate(chi_all)
    flux = du / chi

    fig, axes = plt.subplots(3, 1, figsize=(7.6, 5.8), sharex=True)
    for ax in axes:
        for j in range(2):
            ax.axvspan(j + d[0], j + 1, color=COLOR_B, alpha=0.35, lw=0)
        ax.set_yticks([])
    axes[0].plot(z, 1 / chi, color="black", lw=1.5)
    axes[0].set_ylabel(r"$1/\chi$", rotation=0, labelpad=18)
    axes[1].plot(z, du, color=RED, lw=1.5)
    axes[1].set_ylabel(r"$u'$", rotation=0, labelpad=18)
    axes[2].plot(z, flux, color=GREEN, lw=1.8)
    axes[2].set_ylabel(r"$D=\frac{u'}{\chi}$", rotation=0, labelpad=24)
    axes[2].set_xlabel(r"$z/a$  (sombreado: capa $\chi=9$)")
    axes[0].set_title(r"Dos factores discontinuos cuyo producto es continuo", fontsize=12)
    fig.tight_layout()
    save(fig, "regla_li")


def fig_plano_k() -> None:
    """Autovalores k̃ de la forma compañera en una banda y en un gap."""
    spec = load_spec("figure_03.yaml")
    m, theta, h = 3, np.pi / 12, 8
    tmm = TMMSolver(spec)
    nu = np.linspace(1.2, 3.0, 4000)
    r = tmm.semitrace(m, theta, nu, Polarization.TE).real
    nu_band = float(nu[np.argmin(np.abs(r))])
    nu_gap = float(nu[np.argmax(np.abs(r) * (np.abs(r) < 3))])
    cell = cell_fourier(m, spec.thickness_a, spec.thickness_b, n_max_for(m, h))

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.4), sharey=True)
    for ax, nu0, title in ((axes[0], nu_band, "banda"), (axes[1], nu_gap, "gap")):
        values = bloch_eigenvalues(spec, cell, float(omega_from_nu_ghz(nu0)), theta, Polarization.TE)
        chosen = select_bloch_wavevector(values)
        r_tmm = tmm.semitrace(m, theta, [nu0], Polarization.TE).real[0]
        ax.axvspan(
            -BRILLOUIN_WINDOW,
            BRILLOUIN_WINDOW,
            color=COLOR_A,
            alpha=0.6,
            lw=0,
            label=r"ventana $|\mathrm{Re}\,\tilde k|\leq1.25\pi$",
        )
        ax.plot(values.real / np.pi, values.imag, ".", color=GRAY, ms=5, label=r"$2N$ autovalores $\tilde k$")
        ax.plot(
            chosen.real / np.pi,
            chosen.imag,
            "o",
            color=RED,
            ms=10,
            mfc="none",
            mew=2,
            label=r"elegido: mínimo $|\mathrm{Im}\,\tilde k|$",
        )
        ax.axhline(0, color="black", lw=0.6)
        ax.set_xlim(-4, 4)
        ax.set_ylim(-3, 3)
        ax.set_xlabel(r"$\mathrm{Re}\,\tilde k/\pi$")
        ax.set_title(
            rf"{title}: $\nu={nu0:.3f}$ GHz, $R_{{PWE}}={np.cos(chosen).real:.4f}$"
            rf" ($R_{{TMM}}={r_tmm:.4f}$)",
            fontsize=11,
        )
    axes[0].set_ylabel(r"$\mathrm{Im}\,\tilde k$")
    axes[0].legend(loc="lower left", fontsize=9)
    fig.suptitle(rf"$S_3$, TE, $\theta=\pi/12$, $\nu_m=1$ GHz, $N={cell.size}$", fontsize=12)
    fig.tight_layout()
    save(fig, "plano_k_complejo")


def fig_celda_s4() -> None:
    """Celda S4 con capas A/B y perfil μ(z) escalonado."""
    word = "ABAAB"
    fig, ax = plt.subplots(figsize=(8.0, 1.6))
    for i, layer in enumerate(word):
        ax.add_patch(Rectangle((i, 0), 1, 1, facecolor=COLOR_A if layer == "A" else COLOR_B, edgecolor="#333", lw=1.0))
        ax.text(i + 0.5, 0.5, layer, ha="center", va="center", fontsize=13, weight="bold")
    ax.annotate("", xy=(5, -0.25), xytext=(0, -0.25), arrowprops={"arrowstyle": "<->", "lw": 1.2})
    ax.text(2.5, -0.55, r"$L_4 = 3a + 2b$ (período)", ha="center", fontsize=11)
    ax.set_xlim(-0.2, 5.2)
    ax.set_ylim(-0.8, 1.1)
    ax.axis("off")
    save(fig, "celda_s4")


def main() -> None:
    style()
    print("Generando imágenes en", OUT)
    fig_ondas_planas()
    fig_bloch()
    fig_plegado()
    fig_toeplitz()
    fig_omega_k_vs_k_omega()
    fig_regla_li()
    fig_plano_k()
    fig_celda_s4()


if __name__ == "__main__":
    main()
