"""Forma ω(k) del libro extendida a Drude: contaminación espectral.

Se resuelve el QEP en λ = ω²/c² (regla de Laurent) para k fijo y distintos
truncamientos, y se comparan las frecuencias propias con las raíces exactas
de R_TMM(ν) = cos(k L_m). Los modos espurios cambian con N y se acumulan en
ν_m; por eso la superred se resuelve con la forma k(ω) y la regla inversa.
"""

from __future__ import annotations

from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray
from scipy.optimize import brentq

from fibonacci_photonics.config import load_model, load_spec
from fibonacci_photonics.config.schema import SpectralPollutionConfig
from fibonacci_photonics.io.paths import relative, result_dirs
from fibonacci_photonics.io.results import dump_json
from fibonacci_photonics.parallel import single_threaded_blas
from fibonacci_photonics.physics.electromagnetics import Polarization
from fibonacci_photonics.physics.superlattice import SuperlatticeSpec
from fibonacci_photonics.solvers.pwe.eigenfrequency import drude_qep_frequencies_ghz
from fibonacci_photonics.solvers.tmm import TMMSolver
from fibonacci_photonics.viz.style import PWE_FORMATS, apply_prb_style, save_figure

CONFIG = "pwe/spectral_pollution.yaml"
ROOT_RESIDUAL = 1e-6
"""Un cambio de signo a través de un polo converge con residuo enorme; se descarta."""


def tmm_roots(
    spec: SuperlatticeSpec, m: int, theta: float, k_tilde: float, window: tuple[float, float], nu_m: float
) -> NDArray[np.float64]:
    """Raíces exactas de R_TMM(ν) = cos k̃ en la ventana (GHz)."""
    evaluator = TMMSolver(spec).evaluator(m, theta, Polarization.TE)

    def f(nu: float) -> float:
        return evaluator(nu) - np.cos(k_tilde)

    # Malla más densa cerca de ν_m, donde R oscila sin cota.
    lower = np.linspace(window[0], nu_m - 1e-3, 20000)
    near = nu_m - np.logspace(-3, -9, 4000)
    upper = np.linspace(nu_m + 1e-6, window[1], 20000)
    roots = []
    for grid in (np.concatenate([lower, near]), upper):
        values = np.array([f(x) for x in grid])
        for i in np.flatnonzero(np.sign(values[1:]) != np.sign(values[:-1])):
            if not (np.isfinite(values[i]) and np.isfinite(values[i + 1])):
                continue
            root = brentq(f, grid[i], grid[i + 1], xtol=1e-15)
            if abs(f(root)) < ROOT_RESIDUAL:
                roots.append(root)
    return np.array(sorted(roots))


def main() -> dict[str, Any]:
    single_threaded_blas()
    loaded = load_model(CONFIG, SpectralPollutionConfig)
    cfg = loaded.model
    spec = load_spec(cfg.physics)
    dirs = result_dirs("pwe", "spectral_pollution")
    m, theta, window = cfg.fibonacci_order, cfg.theta_rad, cfg.window_ghz
    k_tilde = cfg.k_over_pi * np.pi
    nu_m = spec.nu_m_ghz()

    exact = tmm_roots(spec, m, theta, k_tilde, window, nu_m)
    spectra = {}
    for n_max in cfg.n_max_values:
        nu = drude_qep_frequencies_ghz(spec, m, k_tilde, theta, n_max)
        spectra[n_max] = nu[(nu >= window[0]) & (nu <= window[1])]
        print(f"  n_max={n_max}: {spectra[n_max].size} autovalores en la ventana")

    def classify(values: NDArray[np.float64]) -> NDArray[np.bool_]:
        if exact.size == 0:
            return np.zeros(values.size, dtype=bool)
        distance = np.min(np.abs(values[:, None] - exact[None, :]), axis=1)
        return distance < 1e-3 * np.maximum(values, 1e-9)

    apply_prb_style()
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.6), gridspec_kw={"width_ratios": [1.3, 1]})
    for i, (n_max, values) in enumerate(spectra.items()):
        good = classify(values)
        axes[0].plot(np.full(good.sum(), i), values[good], "o", color="#2c4d8c", markersize=3,
                     label="coincide con TMM" if i == 0 else None)
        axes[0].plot(np.full((~good).sum(), i), values[~good], "x", color="#c23b22", markersize=3,
                     label="espurio" if i == 0 else None)
        below = values[(values < nu_m) & (values > nu_m - 0.05)]
        axes[1].semilogy(np.full(below.size, 2 * n_max + 1), nu_m - below, "_", color="#c23b22", markersize=8)
    for root in exact:
        axes[0].axhline(root, color="#2c4d8c", linewidth=0.4, alpha=0.6)
    axes[0].axhline(nu_m, color="gray", linestyle="--", linewidth=0.7)
    axes[0].set_xticks(range(len(spectra)))
    axes[0].set_xticklabels([str(2 * n + 1) for n in spectra])
    axes[0].set_xlabel("número de ondas planas N")
    axes[0].set_ylabel(r"$\nu$ (GHz)")
    axes[0].set_ylim(*window)
    axes[0].legend(fontsize=7, loc="center right")
    axes[0].set_title(rf"$\omega(k)$ Drude, m={m}, $kL_m=\pi/2$, $\theta=\pi/12$", fontsize=8)
    axes[1].set_xlabel("número de ondas planas N")
    axes[1].set_ylabel(r"$\nu_m-\nu$ de los autovalores (GHz)")
    axes[1].set_title(r"acumulación en $\nu_m$", fontsize=8)
    fig.tight_layout()
    paths = save_figure(fig, dirs["output"], "spectral_pollution", formats=PWE_FORMATS)

    summary = {
        "config": relative(loaded.path),
        "tmm_roots_ghz": exact.tolist(),
        "spectra_ghz": {str(2 * n + 1): v.tolist() for n, v in spectra.items()},
        "spurious_counts": {str(2 * n + 1): int((~classify(v)).sum()) for n, v in spectra.items()},
        "outputs": {k: relative(p) for k, p in paths.items()},
    }
    dump_json(dirs["data"] / "summary.json", summary)
    print("TMM:", exact)
    return summary
