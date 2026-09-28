"""Comparación cuantitativa entre métodos sobre la misma malla o las mismas bandas."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

import numpy as np

from fibonacci_photonics.solvers.scan import DispersionScan, FrequencyInterval

GRID_ATOL_GHZ = 1e-12


@dataclass(frozen=True)
class ScanAgreement:
    """Métricas punto a punto entre un barrido de referencia y uno candidato.

    Attributes
    ----------
    n_points : int
        Frecuencias comparadas.
    classification_agreement : float
        Fracción de frecuencias con la misma clasificación banda/gap.
    n_classification_mismatch : int
        Frecuencias con clasificación distinta.
    max_abs_delta_r_in_bands, rms_delta_r_in_bands : float
        max y RMS de |ΔR| donde ambos métodos dan banda.
    max_abs_delta_k_over_pi : float
        max |Δ(k L_m/π)| donde ambos dan banda.
    max_abs_r_reference : float
        max |R| finito de la referencia (escala del error en los gaps).
    """

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
    """Acuerdo entre dos barridos de la misma malla (p. ej. TMM de referencia y PWE)."""
    if reference.nu_ghz.shape != candidate.nu_ghz.shape or not np.allclose(
        reference.nu_ghz, candidate.nu_ghz, rtol=0.0, atol=GRID_ATOL_GHZ
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


def edge_errors(reference: list[dict[str, float]], candidate: list[dict[str, float]]) -> list[dict[str, Any]]:
    """Errores de borde y de ancho entre bandas ``{"min", "max"}`` emparejadas por orden en ν."""
    rows = []
    pairs = zip(sorted(reference, key=lambda b: b["min"]), sorted(candidate, key=lambda b: b["min"]), strict=False)
    for t, p in pairs:
        width = t["max"] - t["min"]
        rows.append({
            "tmm": t, "pwe": p,
            "delta_min_ghz": p["min"] - t["min"],
            "delta_max_ghz": p["max"] - t["max"],
            "relative_width_error": ((p["max"] - p["min"]) - width) / width if width > 0 else float("nan"),
        })
    return rows
