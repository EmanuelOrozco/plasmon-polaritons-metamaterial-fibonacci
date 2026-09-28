"""Benchmark del capítulo 1D de Sukhoivanov y Guryev con PWE y TMM.

Reproduce, para la bicapa ε1 = 1 (0.2 µm) / ε2 = 9 (0.8 µm):

(a) la estructura de bandas ω(k) (fig. 4.5 del libro, notebook PWE-1d),
(b) la síntesis de Fourier del perfil 1/ε(x) (fig. 4.6),
(c) el perfil de campo de un modo (notebook) y
(d) el diagrama de bandas proyectado fuera del eje (fig. 4.9), TE y TM.

La TMM da la referencia exacta en todos los casos.
"""

from __future__ import annotations

from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import ArrayLike, NDArray

from fibonacci_photonics.config import load_model, spec_from_mapping
from fibonacci_photonics.config.schema import BookBenchmarkConfig
from fibonacci_photonics.io.paths import relative, result_dirs
from fibonacci_photonics.io.results import dump_json, save_npz
from fibonacci_photonics.physics.electromagnetics import Polarization
from fibonacci_photonics.physics.superlattice import SuperlatticeSpec
from fibonacci_photonics.physics.units import omega_from_nu_ghz
from fibonacci_photonics.solvers.pwe.eigenfrequency import mode_profile, nondispersive_bands
from fibonacci_photonics.solvers.pwe.fourier import bilayer_fourier, synthesize_profile
from fibonacci_photonics.solvers.tmm.transfer_matrix import semitrace_by_recurrence
from fibonacci_photonics.viz.style import PWE_FORMATS, PWE_STYLE, TMM_STYLE, apply_prb_style, save_figure

CONFIG = "pwe/book_benchmark.yaml"


def book_spec(cfg: BookBenchmarkConfig) -> SuperlatticeSpec:
    """La bicapa como superred m = 2 (S2 = AB) con B 'Drude' de ωe = ωm = 0 y ε0 = ε2."""
    return spec_from_mapping(
        {
            "layer_a_thickness_mm": cfg.thickness_1_um * 1e-3,
            "layer_b_thickness_mm": cfg.thickness_2_um * 1e-3,
            "epsilon_a": cfg.epsilon_1,
            "omega_e_over_2pi_ghz": 0.0,
            "omega_m_over_2pi_ghz": 0.0,
            "epsilon_0": cfg.epsilon_2,
            "source": "Sukhoivanov & Guryev, cap. 4 (1D)",
        }
    )


def bilayer_semitrace_fixed_q(
    f: ArrayLike, q: ArrayLike, fractions: tuple[float, float], eps: tuple[float, float], polarization: Polarization
) -> NDArray[np.float64]:
    """R = ½ Tr T de la bicapa con q fijo (ec. (9) del paper con Q_j = √(κ²ε_j − q²)).

    ``f`` = ωa/2πc y ``q`` = q a/2π; vale también fuera del cono de luz del aire.
    """
    kappa = 2 * np.pi * np.asarray(f, dtype=complex)
    qa = 2 * np.pi * np.asarray(q, dtype=complex)
    q1 = np.sqrt(kappa**2 * eps[0] - qa**2)
    q2 = np.sqrt(kappa**2 * eps[1] - qa**2)
    chi1, chi2 = eps if polarization == Polarization.TM else (1.0, 1.0)
    d1, d2 = fractions
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = 0.5 * (q1 * chi2 / (q2 * chi1) + q2 * chi1 / (q1 * chi2))
        r = np.cos(q1 * d1) * np.cos(q2 * d2) - ratio * np.sin(q1 * d1) * np.sin(q2 * d2)
    return np.real(r)


def main() -> dict[str, Any]:
    loaded = load_model(CONFIG, BookBenchmarkConfig)
    cfg = loaded.model
    dirs = result_dirs("pwe", "book_benchmark")
    spec = book_spec(cfg)
    l1, l2 = cfg.thickness_1_um, cfg.thickness_2_um
    eps = (cfg.epsilon_1, cfg.epsilon_2)
    period_m = (l1 + l2) * 1e-6
    n_bands = cfg.n_bands
    ks = np.linspace(-np.pi, np.pi, cfg.k_points)
    fractions = (l1 / (l1 + l2), l2 / (l1 + l2))
    cell = bilayer_fourier(*fractions, cfg.n_max)

    def tmm_r(f: ArrayLike) -> NDArray[np.float64]:
        omega = omega_from_nu_ghz(np.asarray(f) * spec.speed_of_light / period_m / 1e9)
        return np.real(semitrace_by_recurrence(spec, 2, omega, 0.0, Polarization.TM))

    bands = {rule: nondispersive_bands(cell, eps, (1, 1), ks, rule=rule, n_bands=n_bands)
             for rule in ("laurent", "inverse")}

    # Referencia TMM: ν(k) exacta desde cos(k a) = R(ν) en una malla fina.
    f_grid = np.linspace(1e-4, float(bands["inverse"][:, -1].max()) * 1.02, 60001)
    r_grid = tmm_r(f_grid)
    allowed = np.abs(r_grid) <= 1
    k_tmm = np.where(allowed, np.arccos(np.clip(r_grid, -1, 1)), np.nan)

    errors = {}
    for rule, data in bands.items():
        residual = [np.abs(tmm_r(row) - np.cos(k)) for k, row in zip(ks, data, strict=True)]
        errors[rule] = float(np.nanmax(np.array(residual)[:, 1:]))
    save_npz(dirs["data"] / "bands.npz", k=ks, laurent=bands["laurent"], inverse=bands["inverse"],
             f_tmm=f_grid, k_tmm=k_tmm)

    apply_prb_style()
    # (a) estructura de bandas
    fig, ax = plt.subplots(figsize=(4.2, 4.0))
    ax.plot(k_tmm / np.pi, f_grid, **TMM_STYLE, label="TMM (exacta)")
    ax.plot(-k_tmm / np.pi, f_grid, **TMM_STYLE)
    step = max(1, ks.size // 40)
    for j in range(n_bands):
        ax.plot(ks[::step] / np.pi, bands["laurent"][::step, j], **{**PWE_STYLE, "color": "#2e8b4a", "marker": "s"},
                label="PWE, Laurent (libro)" if j == 0 else None)
        ax.plot(ks[::step] / np.pi, bands["inverse"][::step, j], **PWE_STYLE,
                label="PWE, regla inversa" if j == 0 else None)
    ax.set_xlim(-1, 1)
    ax.set_ylim(0, float(bands["inverse"][:, -1].max()))
    ax.set_xlabel(r"$k a/\pi$")
    ax.set_ylabel(r"$\omega a / 2\pi c$")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=3, fontsize=6.5, frameon=False)
    ax.set_title(rf"$\varepsilon_1={eps[0]:g}$, $\varepsilon_2={eps[1]:g}$, $N={2 * cfg.n_max + 1}$")
    fig.tight_layout()
    out = {"band_structure": save_figure(fig, dirs["output"], "band_structure", formats=PWE_FORMATS)}

    # (b) síntesis de Fourier de 1/ε
    fig, ax = plt.subplots(figsize=(5.0, 3.2))
    x = np.linspace(-1.0, 1.0, 4001)
    exact = np.where(np.mod(x, 1.0) < fractions[0], 1 / eps[0], 1 / eps[1])
    ax.plot(x, exact, color="black", linewidth=1.2, label="perfil exacto")
    colors = ["#c23b22", "#2e8b4a", "#2c4d8c"]
    for n_max, color in zip(cfg.synthesis_n_max, colors, strict=False):
        c = bilayer_fourier(*fractions, n_max)
        coeffs = c.profile_coefficients(1 / eps[0], 1 / eps[1])
        mask = np.abs(c.orders) <= n_max
        values = synthesize_profile(coeffs[mask], c.orders[mask], 1.0, x).real
        ax.plot(x, values, color=color, linewidth=0.7, label=f"{2 * n_max + 1} ondas planas")
    ax.set_xlabel(r"$x/a$")
    ax.set_ylabel(r"$1/\varepsilon(x)$")
    ax.legend(fontsize=7, loc="center right")
    fig.tight_layout()
    out["fourier_synthesis"] = save_figure(fig, dirs["output"], "fourier_synthesis", formats=PWE_FORMATS)

    # (c) perfil de campo del modo (como el notebook del libro)
    k_field = cfg.field_k_over_pi * np.pi
    band_index = cfg.field_band
    freqs, modes = nondispersive_bands(cell, eps, (1, 1), [k_field], rule="inverse",
                                       n_bands=band_index + 1, return_modes=True)
    z = np.linspace(-1.5, 1.5, 3001)
    h = mode_profile(cell, modes[0][:, band_index], k_field, z)
    fig, ax = plt.subplots(figsize=(5.0, 3.0))
    ax.fill_between(z, 0, np.where(np.mod(z, 1.0) < fractions[0], 0, 1), color="#dddddd", step="mid",
                    label=r"capa $\varepsilon_2$")
    intensity = np.abs(h) ** 2
    ax.plot(z, intensity / intensity.max(), color="#2c4d8c", linewidth=0.9, label=r"$|H(x)|^2$")
    ax.set_xlabel(r"$x/a$")
    ax.set_ylabel("intensidad normalizada")
    ax.set_title(rf"banda {band_index + 1}, $k = \pi/(15a)$, $\omega a/2\pi c = {freqs[0, band_index]:.4f}$")
    ax.legend(fontsize=7, loc="upper right")
    fig.tight_layout()
    out["field_profile"] = save_figure(fig, dirs["output"], "field_profile", formats=PWE_FORMATS)

    # (d) diagrama proyectado fuera del eje: TE y TM
    q_values = np.linspace(0.0, cfg.off_axis_q_max, cfg.off_axis_q_points)
    cell_off = bilayer_fourier(*fractions, cfg.off_axis_n_max)
    f_map = np.linspace(1e-3, 1.0, 700)
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.6), sharey=True)
    projected = {}
    for ax, pol in zip(axes, (Polarization.TM, Polarization.TE), strict=True):
        def edge_band(k: float, pol: Polarization = pol) -> NDArray[np.float64]:
            return np.array([nondispersive_bands(cell_off, eps, (1, 1), [k], q_tilde=2 * np.pi * q,
                                                 polarization=pol, rule="inverse", n_bands=6)[0] for q in q_values])

        e0, epi = edge_band(0.0), edge_band(np.pi)
        projected[pol.value] = {"k0": e0.tolist(), "kpi": epi.tolist()}
        # Fondo TMM: regiones |R| ≤ 1 en el plano (q, ω).
        qq, ff = np.meshgrid(q_values, f_map)
        r_map = bilayer_semitrace_fixed_q(ff, qq, fractions, eps, pol)
        ax.contourf(qq, ff, (np.abs(r_map) <= 1).astype(float), levels=[0.5, 1.5], colors=["#c9d6ea"])
        for j in range(e0.shape[1]):
            ax.plot(q_values, np.minimum(e0[:, j], epi[:, j]), color="#c23b22", linewidth=0.6)
            ax.plot(q_values, np.maximum(e0[:, j], epi[:, j]), color="#c23b22", linewidth=0.6, linestyle="--")
        ax.plot(q_values, q_values / np.sqrt(eps[0]), color="gray", linewidth=0.7, linestyle=":",
                label=r"línea de luz, medio 1")
        if pol == Polarization.TM:
            # Ec. (4.45): q = (ω/c) n1 n2 / sqrt(n1² + n2²), sin reflexión de Fresnel TM.
            brewster = np.sqrt((eps[0] + eps[1]) / (eps[0] * eps[1]))
            ax.plot(q_values, brewster * q_values, color="black", linewidth=0.8, label="Brewster, ec. (4.45)")
        ax.legend(fontsize=6.5, loc="lower right")
        ax.set_xlim(0, cfg.off_axis_q_max)
        ax.set_ylim(0, 1.0)
        ax.set_xlabel(r"$q a/2\pi$")
        ax.set_title(f"{pol.value}: bandas (azul, TMM) y bordes PWE (rojo)")
    axes[0].set_ylabel(r"$\omega a/2\pi c$")
    fig.tight_layout()
    out["off_axis"] = save_figure(fig, dirs["output"], "off_axis_projected_bands", formats=PWE_FORMATS)

    summary = {
        "config": relative(loaded.path),
        "n_plane_waves": 2 * cfg.n_max + 1,
        "max_abs_tmm_residual_laurent": errors["laurent"],
        "max_abs_tmm_residual_inverse": errors["inverse"],
        "band_edges_k0_first4_inverse": bands["inverse"][ks.size // 2, :4].tolist(),
        "field_mode_frequency": float(freqs[0, band_index]),
        "projected_band_edges": projected,
        "outputs": {name: {k: relative(p) for k, p in paths.items()} for name, paths in out.items()},
    }
    dump_json(dirs["data"] / "summary.json", summary)
    print("Benchmark del libro: residuo TMM máx. Laurent =", errors["laurent"], " inversa =", errors["inverse"])
    return summary
