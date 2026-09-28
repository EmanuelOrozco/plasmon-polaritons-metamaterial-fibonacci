"""Raíces acotadas con el método de Brent."""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
from scipy.optimize import brentq

DEFAULT_XTOL_GHZ = 1.0e-13
"""Tolerancia absoluta de Brent sobre ν, en GHz."""


def refine_crossing(
    func: Callable[[float], float],
    lo: float,
    hi: float,
    *,
    xtol: float = DEFAULT_XTOL_GHZ,
) -> float | None:
    """Raíz de ``func`` en [lo, hi] por Brent; ``None`` si no hay cambio de signo.

    Brent converge de forma garantizada (superlineal cerca de la raíz) siempre
    que f(lo)·f(hi) ≤ 0; valores no finitos en los extremos anulan el intervalo.
    """
    f_lo, f_hi = func(lo), func(hi)
    if not (np.isfinite(f_lo) and np.isfinite(f_hi)) or f_lo * f_hi > 0.0:
        return None
    if f_lo == 0.0:
        return lo
    if f_hi == 0.0:
        return hi
    return float(brentq(func, lo, hi, xtol=xtol, rtol=4.0 * np.finfo(float).eps))
