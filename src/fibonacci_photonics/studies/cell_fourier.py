"""Perfiles ε(z), μ(z) de la celda S4 y su serie de Fourier truncada.

La celda de Fibonacci es un perfil escalonado con F_m capas; sus coeficientes
decaen como 1/|n| y la síntesis truncada muestra el fenómeno de Gibbs en cada
interfaz.
"""

from __future__ import annotations

from typing import Any

import matplotlib.pyplot as plt
import numpy as np

from fibonacci_photonics.config import load_spec
from fibonacci_photonics.io.paths import relative, result_dirs
from fibonacci_photonics.io.results import dump_json
from fibonacci_photonics.physics.fibonacci import n_layers_total
from fibonacci_photonics.physics.units import omega_from_nu_ghz
from fibonacci_photonics.solvers.pwe.fourier import cell_fourier, layer_edges, synthesize_profile
from fibonacci_photonics.viz.style import PWE_FORMATS, apply_prb_style, save_figure

M = 4
NU_GHZ = 2.0


def main() -> dict[str, Any]:
    spec = load_spec("figure_03.yaml")
    dirs = result_dirs("pwe", "cell_fourier")
    omega = float(omega_from_nu_ghz(NU_GHZ))
    eps_b = float(spec.medium_b.epsilon(omega).real)
    mu_b = float(spec.medium_b.mu(omega).real)
    word, edges = layer_edges(M, 1.0, 1.0)
    bounds = list(zip(word, edges[:-1] / edges[-1], edges[1:] / edges[-1], strict=True))
    z = np.linspace(0.0, 1.0, 4001)
    exact_mu = np.ones_like(z)
    for symbol, lo, hi in bounds:
        if symbol == "B":
            exact_mu[(z >= lo) & (z < hi)] = mu_b

    apply_prb_style()
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.1), gridspec_kw={"width_ratios": [1.5, 1]})
    ax = axes[0]
    for symbol, lo, hi in bounds:
        ax.axvspan(lo, hi, color="#dde6f3" if symbol == "B" else "white", linewidth=0)
        ax.text(0.5 * (lo + hi), 1.12, symbol, ha="center", fontsize=8)
    ax.plot(z, exact_mu, color="black", linewidth=1.1, label=r"$\mu(z)$ exacto")
    for harmonics, color in ((2, "#c23b22"), (8, "#2c4d8c")):
        n_max = harmonics * n_layers_total(M)
        cell = cell_fourier(M, 1.0, 1.0, n_max)
        mask = np.abs(cell.orders) <= n_max
        values = synthesize_profile(cell.profile_coefficients(1.0, mu_b)[mask], cell.orders[mask], 1.0, z).real
        ax.plot(z, values, color=color, linewidth=0.7, label=f"N = {2 * n_max + 1}")
    ax.set_xlim(0, 1)
    ax.set_ylim(min(mu_b, 1) - 0.25, 1.22)
    ax.set_xlabel(r"$z/L_4$")
    ax.set_ylabel(r"$\mu(z)$")
    ax.set_title(rf"celda $S_4$ = {word}, $\nu$ = {NU_GHZ:g} GHz ($\mu_B$ = {mu_b:.2f}, $\varepsilon_B$ = {eps_b:.2f})",
                 fontsize=8)
    ax.legend(fontsize=6.5, loc="lower right")

    ax = axes[1]
    cell = cell_fourier(M, 1.0, 1.0, 200)
    coeffs = np.abs(cell.profile_coefficients(1.0, mu_b))
    positive = cell.orders > 0
    n = cell.orders[positive]
    ax.loglog(n, np.maximum(coeffs[positive], 1e-18), ".", markersize=2, color="#2c4d8c", label=r"$|\mu_n|$")
    ax.loglog(n, (1 - mu_b) / (np.pi * n) * 2, color="gray", linestyle=":", linewidth=0.8, label=r"$\propto 1/n$")
    ax.set_ylim(1e-5, 1)
    ax.set_xlabel(r"orden $n$")
    ax.set_ylabel(r"$|\mu_n|$")
    ax.legend(fontsize=7)
    fig.tight_layout()
    paths = save_figure(fig, dirs["output"], "cell_fourier", formats=PWE_FORMATS)
    summary = {"m": M, "nu_ghz": NU_GHZ, "epsilon_b": eps_b, "mu_b": mu_b, "word": word,
               "outputs": {k: relative(v) for k, v in paths.items()}}
    dump_json(dirs["data"] / "summary.json", summary)
    print("Celda S4 escrita en", relative(dirs["output"]))
    return summary
