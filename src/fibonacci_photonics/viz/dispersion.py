"""Curvas de dispersión ±k(ν) como en el PRB, cortadas en gaps y saltos de k."""

from __future__ import annotations

from typing import Any

import numpy as np
from matplotlib.axes import Axes
from numpy.typing import ArrayLike, NDArray

from fibonacci_photonics.solvers.scan import DispersionScan
from fibonacci_photonics.viz.style import latex_theta

FloatArray = NDArray[np.float64]
BoolArray = NDArray[np.bool_]


def _fill_isolated_nans(k: ArrayLike, max_dk: float, protected: BoolArray | None = None) -> FloatArray:
    """Cierra huecos de un solo punto (artefacto de malla en k = 0).

    Las muestras ``protected`` no se rellenan: son gaps reales delimitados por bordes.
    """
    k = np.asarray(k, dtype=float).copy()
    n = k.size
    for i in range(1, n - 1):
        if protected is not None and protected[i]:
            continue
        if np.isfinite(k[i]) or not np.isfinite(k[i - 1]) or not np.isfinite(k[i + 1]):
            continue
        if abs(k[i - 1] - k[i + 1]) < max_dk:
            k[i] = 0.5 * (k[i - 1] + k[i + 1])
    return k


def _monotonic_run_to_anchor(k: FloatArray, anchors: BoolArray, left: int, right: int) -> bool:
    """¿El salto k[left] → k[right] sigue monótono, sin huecos, hasta un borde exacto?"""
    step = np.sign(k[right] - k[left])
    if step == 0.0:
        return True
    for start, direction in ((right, 1), (left, -1)):
        i = start
        while True:
            if anchors[i]:
                return True
            j = i + direction
            if j < 0 or j >= k.size or not np.isfinite(k[j]):
                break
            if np.sign(k[max(i, j)] - k[min(i, j)]) not in (step, 0.0):
                break
            i = j
    return False


def _polyline_with_breaks(
    k: ArrayLike, nu: ArrayLike, max_dk: float, anchors: ArrayLike | None = None
) -> tuple[FloatArray, FloatArray]:
    """Inserta NaN cuando k salta, para no cruzar la zona de Brillouin con una recta.

    Con ``anchors`` (bordes de banda exactos) no se corta un salto que forma parte
    de un tramo monótono que termina en un borde: cerca del borde k ∝ √|ν − ν_borde|
    varía rápido y, dentro de una banda, k es monótono hasta su borde.
    """
    k = np.asarray(k, dtype=float)
    nu = np.asarray(nu, dtype=float)
    if k.size < 2:
        return k, nu
    k_out = k.copy()
    finite = np.isfinite(k)
    jump = finite & np.concatenate(([False], np.abs(np.diff(k)) > max_dk))
    if anchors is not None:
        anchor_mask = np.asarray(anchors, dtype=bool)
        for i in np.flatnonzero(jump):
            if _monotonic_run_to_anchor(k, anchor_mask, int(i) - 1, int(i)):
                jump[i] = False
    k_out[jump] = np.nan
    return k_out, nu


def _merge_band_edges(
    nu: FloatArray,
    k: FloatArray,
    edge_nu: FloatArray,
    edge_k: FloatArray,
    max_dk: float,
) -> tuple[FloatArray, FloatArray, BoolArray]:
    """Intercala los puntos de cierre entre las muestras de la malla.

    Un punto con k finito (borde o punta) es un ancla de la curva; uno con k = NaN
    marca un gap oculto entre dos muestras y corta la curva.
    """
    protected = np.zeros(nu.size, dtype=bool)
    idx = np.searchsorted(nu, edge_nu)
    protected[np.clip(idx - 1, 0, nu.size - 1)] = True
    protected[np.clip(idx, 0, nu.size - 1)] = True
    protected &= ~np.isfinite(k)
    k = _fill_isolated_nans(k, max_dk, protected)
    nu_all = np.concatenate((nu, edge_nu))
    k_all = np.concatenate((k, edge_k))
    anchors = np.concatenate((np.zeros(nu.size, dtype=bool), np.isfinite(edge_k)))
    order = np.argsort(nu_all, kind="stable")
    return nu_all[order], k_all[order], anchors[order]


def plot_dispersion_branches(
    ax: Axes,
    scan: DispersionScan,
    *,
    color: str | tuple[float, float, float],
    linestyle: str,
    linewidth: float = 0.75,
    alpha: float = 1.0,
    label: str | None = None,
    max_dk: float = 0.12,
    edge_nu_ghz: ArrayLike | None = None,
    edge_k_lm_over_pi: ArrayLike | None = None,
) -> None:
    """Curvas ±k(ν) cortadas en gaps y en saltos de k (como en el PRB).

    Con ``edge_nu_ghz``/``edge_k_lm_over_pi`` (bordes |R| = 1 y puntas refinadas)
    cada banda se cierra en su borde exacto, k = 0 o k = ±1, en lugar de terminar
    en la última muestra de la malla.
    """
    nu = scan.nu_ghz
    k = scan.k_lm_over_pi
    allowed = scan.allowed
    if not np.any(allowed):
        return
    k_pos = np.where(allowed, k, np.nan)
    if edge_nu_ghz is not None and np.size(edge_nu_ghz):
        nu, k_pos, anchors = _merge_band_edges(
            nu,
            k_pos,
            np.asarray(edge_nu_ghz, dtype=float),
            np.asarray(edge_k_lm_over_pi, dtype=float),
            max_dk,
        )
        k_pos, nu_p = _polyline_with_breaks(k_pos, nu, max_dk, anchors)
    else:
        k_pos = _fill_isolated_nans(k_pos, max_dk)
        k_pos, nu_p = _polyline_with_breaks(k_pos, nu, max_dk)
    k_neg = -k_pos
    line: dict[str, Any] = {"color": color, "linestyle": linestyle, "linewidth": linewidth, "alpha": alpha}
    ax.plot(k_neg, nu_p, label=label, **line)
    ax.plot(k_pos, nu_p, **line)


def format_dispersion_axes(
    ax: Axes,
    *,
    nu_min: float,
    nu_max: float,
    panel: str | None = None,
    theta_label: str | None = None,
    m_label: str | None = None,
) -> None:
    """Ejes k L_m/π ∈ [−1, 1] y ν ∈ [nu_min, nu_max] GHz, con letra de panel, θ y m."""
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
