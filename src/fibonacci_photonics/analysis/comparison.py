"""Comparación cuantitativa entre dos barridos de dispersión (p. ej. TMM y PWE)."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import asdict, dataclass

import numpy as np
from scipy.optimize import brentq

from fibonacci_photonics.core.bands import DispersionScan, FrequencyInterval


@dataclass(frozen=True)
class ScanAgreement:
    """Métricas punto a punto sobre la misma malla de frecuencias."""

    n_points: int
    classification_agreement: float
    n_classification_mismatch: int
    max_abs_delta_r_in_bands: float
    rms_delta_r_in_bands: float
    max_abs_delta_k_over_pi: float
    max_abs_r_reference: float

    def as_dict(self) -> dict[str, float | int]:
        return asdict(self)


def compare_scans(reference: DispersionScan, candidate: DispersionScan) -> ScanAgreement:
    if reference.nu_ghz.shape != candidate.nu_ghz.shape or not np.allclose(
        reference.nu_ghz, candidate.nu_ghz, rtol=0.0, atol=1e-12
    ):
        raise ValueError("los barridos deben compartir la malla de frecuencias")
    both = reference.allowed & candidate.allowed
    delta_r = np.abs(reference.r.real - candidate.r.real)
    delta_k = np.abs(reference.k_lm_over_pi - candidate.k_lm_over_pi)
    in_bands = delta_r[both]
    return ScanAgreement(
        n_points=int(reference.nu_ghz.size),
        classification_agreement=float(np.mean(reference.allowed == candidate.allowed)),
        n_classification_mismatch=int(np.count_nonzero(reference.allowed != candidate.allowed)),
        max_abs_delta_r_in_bands=float(in_bands.max()) if in_bands.size else 0.0,
        rms_delta_r_in_bands=float(np.sqrt(np.mean(in_bands**2))) if in_bands.size else 0.0,
        max_abs_delta_k_over_pi=float(np.nanmax(delta_k[both])) if np.any(both) else 0.0,
        max_abs_r_reference=float(np.nanmax(reference.abs_r[np.isfinite(reference.abs_r)])),
    )


def match_intervals(
    reference: list[FrequencyInterval],
    candidate: list[FrequencyInterval],
) -> list[tuple[FrequencyInterval, FrequencyInterval | None]]:
    """Empareja cada banda de referencia con la banda candidata de mayor solape."""
    pairs: list[tuple[FrequencyInterval, FrequencyInterval | None]] = []
    for ref in reference:
        best, best_overlap = None, 0.0
        for cand in candidate:
            overlap = min(ref.nu_max_ghz, cand.nu_max_ghz) - max(ref.nu_min_ghz, cand.nu_min_ghz)
            if overlap >= best_overlap and overlap >= 0.0:
                best, best_overlap = cand, overlap
        pairs.append((ref, best))
    return pairs


def refine_crossing(
    func: Callable[[float], float],
    lo: float,
    hi: float,
    *,
    xtol: float = 1.0e-13,
) -> float | None:
    """Raíz de ``func`` en [lo, hi] por Brent; None si no hay cambio de signo."""
    f_lo, f_hi = func(lo), func(hi)
    if not (np.isfinite(f_lo) and np.isfinite(f_hi)) or f_lo * f_hi > 0.0:
        return None
    if f_lo == 0.0:
        return lo
    if f_hi == 0.0:
        return hi
    return float(brentq(func, lo, hi, xtol=xtol, rtol=4.0 * np.finfo(float).eps))


def refine_band_edges(
    r_of_nu: Callable[[float], float],
    interval: FrequencyInterval,
    step_ghz: float,
    *,
    xtol: float = 1.0e-13,
) -> tuple[float, float]:
    """Bordes exactos de una banda donde |R| cruza 1, a partir de la malla gruesa.

    ``interval`` son los puntos de malla extremos con |R| ≤ 1; los cruces están
    a menos de un paso ``step_ghz`` hacia fuera.
    """

    def excess(nu: float) -> float:
        return abs(r_of_nu(nu)) - 1.0

    lo = refine_crossing(excess, interval.nu_min_ghz - step_ghz, interval.nu_min_ghz, xtol=xtol)
    hi = refine_crossing(excess, interval.nu_max_ghz, interval.nu_max_ghz + step_ghz, xtol=xtol)
    return (
        interval.nu_min_ghz if lo is None else lo,
        interval.nu_max_ghz if hi is None else hi,
    )
