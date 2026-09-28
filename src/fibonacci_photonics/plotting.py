"""Estilo gráfico alineado con el Brief Report de Physical Review B."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from fibonacci_photonics.core.bands import DispersionScan

# Paleta del original: S3 discontinuo rojo, S4 continuo azul.
# Figs. 5–6, de menor a mayor ν: negro, rojo, verde, azul, magenta, cian.
S3_STYLE = {"color": "#c23b22", "linestyle": "--", "linewidth": 0.65}
S4_STYLE = {"color": "#2c4d8c", "linestyle": "-", "linewidth": 0.7}
BRANCH_COLORS = ["#1a1a1a", "#c23b22", "#2e8b4a", "#2c4d8c", "#c44ec0", "#17becf", "#d35400", "#7f7f7f"]

# Comparación de métodos: TMM en línea continua, PWE en marcadores huecos.
TMM_STYLE = {"color": "#2c4d8c", "linestyle": "-", "linewidth": 0.8}
PWE_STYLE = {"color": "#1a1a1a", "marker": "o", "markersize": 2.2, "markerfacecolor": "none",
             "markeredgewidth": 0.5, "linestyle": "none"}


def latex_theta(label: str) -> str:
    """Convierte 'pi/12' → '\\pi/12' si el caller aún no escapó π."""
    text = str(label)
    if r"\pi" in text:
        return text
    return text.replace("pi", r"\pi")


def apply_prb_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 9,
            "axes.labelsize": 9,
            "axes.titlesize": 9,
            "legend.fontsize": 8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "axes.linewidth": 0.7,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "pdf.fonttype": 42,
        }
    )


def _fill_isolated_nans(k: np.ndarray, max_dk: float) -> np.ndarray:
    """Cierra huecos de un solo punto (artefacto de malla en k = 0)."""
    k = np.asarray(k, dtype=float).copy()
    n = k.size
    for i in range(1, n - 1):
        if np.isfinite(k[i]) or not np.isfinite(k[i - 1]) or not np.isfinite(k[i + 1]):
            continue
        if abs(k[i - 1] - k[i + 1]) < max_dk:
            k[i] = 0.5 * (k[i - 1] + k[i + 1])
    return k


def _polyline_with_breaks(k: np.ndarray, nu: np.ndarray, max_dk: float) -> tuple[np.ndarray, np.ndarray]:
    """Inserta NaN cuando k salta, para no cruzar la zona de Brillouin con una recta."""
    k = np.asarray(k, dtype=float)
    nu = np.asarray(nu, dtype=float)
    if k.size < 2:
        return k, nu
    k_out = k.copy()
    finite = np.isfinite(k)
    jump = finite & np.concatenate(([False], np.abs(np.diff(k)) > max_dk))
    k_out[jump] = np.nan
    return k_out, nu


def plot_dispersion_branches(
    ax: Axes,
    scan: DispersionScan,
    *,
    color: str,
    linestyle: str,
    linewidth: float = 0.75,
    label: str | None = None,
    max_dk: float = 0.12,
) -> None:
    """Curvas ±k(ω) cortadas en gaps y en saltos de k (como en el PRB)."""
    nu = scan.nu_ghz
    k = scan.k_lm_over_pi
    allowed = scan.allowed
    if not np.any(allowed):
        return
    k_pos = np.where(allowed, k, np.nan)
    k_pos = _fill_isolated_nans(k_pos, max_dk)
    k_pos, nu_p = _polyline_with_breaks(k_pos, nu, max_dk)
    k_neg = -k_pos
    ax.plot(k_neg, nu_p, color=color, linestyle=linestyle, linewidth=linewidth, label=label)
    ax.plot(k_pos, nu_p, color=color, linestyle=linestyle, linewidth=linewidth)


def format_dispersion_axes(
    ax: Axes,
    *,
    nu_min: float,
    nu_max: float,
    panel: str | None = None,
    theta_label: str | None = None,
    m_label: str | None = None,
) -> None:
    ax.set_xlim(-1.0, 1.0)
    ax.set_ylim(nu_min, nu_max)
    ax.set_xlabel(r"$k L_m / \pi$")
    ax.set_ylabel("frequency (GHz)")
    ax.set_xticks([-1.0, -0.5, 0.0, 0.5, 1.0])
    if panel:
        ax.text(0.03, 0.95, rf"({panel})", transform=ax.transAxes, va="top", ha="left")
    if theta_label:
        ax.text(
            0.97,
            0.55,
            rf"$\theta = {latex_theta(theta_label)}$",
            transform=ax.transAxes,
            va="center",
            ha="right",
        )
    if m_label:
        ax.text(0.97, 0.92, rf"$m={m_label}$", transform=ax.transAxes, va="top", ha="right")


def save_figure(
    fig: Figure,
    output_dir: Path,
    stem: str,
    formats: tuple[str, ...] = ("pdf", "png", "svg"),
) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {fmt: output_dir / f"{stem}.{fmt}" for fmt in formats}
    for fmt, path in paths.items():
        fig.savefig(path, dpi=300 if fmt == "png" else None)
    plt.close(fig)
    return paths
