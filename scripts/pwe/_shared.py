"""Piezas comunes de los scripts PWE (barridos cronometrados y subbandas refinadas)."""

from __future__ import annotations

import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
for path in (ROOT / "src", ROOT / "scripts"):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from _common import single_threaded_blas, worker_count  # noqa: E402
from fibonacci_photonics.analysis.band_edges import locate_bands, plasmon_grid, validate_bands  # noqa: E402
from fibonacci_photonics.core.bands import DispersionScan  # noqa: E402
from fibonacci_photonics.core.electromagnetics import Polarization  # noqa: E402
from fibonacci_photonics.core.params import SuperlatticeSpec  # noqa: E402
from fibonacci_photonics.pwe.bloch import n_max_for  # noqa: E402
from fibonacci_photonics.pwe.dispersion import SemitraceEvaluator, scan_dispersion  # noqa: E402


CHECK_FACTOR = 1.5
"""Las bandas candidatas se validan con 1.5 veces más ondas planas."""


def harmonics_for(numerics: dict[str, Any], m: int) -> int:
    value = numerics["harmonics_per_layer"]
    if isinstance(value, dict):
        return int(value[m])
    return int(value)


def timed_pwe_scan(
    spec: SuperlatticeSpec,
    m: int,
    theta: float,
    nu: np.ndarray,
    polarization: Polarization,
    numerics: dict[str, Any],
) -> tuple[DispersionScan, float]:
    single_threaded_blas()
    start = time.perf_counter()
    scan = scan_dispersion(
        spec,
        m,
        theta,
        nu,
        polarization,
        harmonics_per_layer=harmonics_for(numerics, m),
        rule=numerics.get("rule", "inverse"),
        workers=worker_count(),
    )
    return scan, time.perf_counter() - start


def scan_arrays(scan: DispersionScan, seconds: float) -> dict[str, np.ndarray]:
    return {
        "nu_ghz": scan.nu_ghz,
        "r_real": np.real(scan.r),
        "abs_r": scan.abs_r,
        "allowed": scan.allowed.astype(np.uint8),
        "k_lm_over_pi": scan.k_lm_over_pi,
        "m": np.array(scan.m),
        "theta": np.array(scan.theta),
        "seconds": np.array(seconds),
    }


def plasmon_bands_task(task: dict[str, Any]) -> dict[str, Any]:
    """Subbandas de un (m, θ) con bordes refinados por Brent sobre R_PWE."""
    spec: SuperlatticeSpec = task["spec"]
    m, theta = task["m"], task["theta"]
    grid = plasmon_grid(task["nu_min"], spec.nu_m_ghz(), task["grid_points"], task["min_offset"])
    evaluator = SemitraceEvaluator(
        spec, m, theta, task["polarization"], n_max=n_max_for(m, task["harmonics"]), rule=task["rule"]
    )
    check = SemitraceEvaluator(
        spec, m, theta, task["polarization"],
        n_max=int(round(CHECK_FACTOR * n_max_for(m, task["harmonics"]))), rule=task["rule"],
    )
    start = time.perf_counter()
    r_grid = evaluator.many(grid)
    candidates = locate_bands(evaluator, grid, r_grid)
    bands, rejected = validate_bands(candidates, evaluator, check)
    return {
        "m": m,
        "theta": theta,
        "harmonics_per_layer": task["harmonics"],
        "n_max": n_max_for(m, task["harmonics"]),
        "bands": [{"min": b.nu_min_ghz, "max": b.nu_max_ghz} for b in bands],
        "rejected_bands": [{"min": b.nu_min_ghz, "max": b.nu_max_ghz} for b in rejected],
        "evaluations": evaluator.calls + check.calls,
        "seconds": time.perf_counter() - start,
    }
