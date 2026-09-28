"""Métodos numéricos que calculan la semitraza R_m(ν) = cos(k L_m).

Ambos implementan el protocolo ``DispersionSolver``:

- ``TMMSolver``: matriz de transferencia, exacta y O(m) por frecuencia.
- ``PWESolver``: ondas planas en la forma k(ω), con truncamiento N = 2 h F_m + 1.

``get_solver("tmm" | "pwe", spec, **opciones)`` los construye por nombre.
"""

from __future__ import annotations

from typing import Any

from fibonacci_photonics.errors import InvalidParameterError
from fibonacci_photonics.physics.superlattice import SuperlatticeSpec
from fibonacci_photonics.solvers.base import (
    DispersionSolver,
    SemitraceFunction,
    validate_angle,
    validate_frequencies,
)
from fibonacci_photonics.solvers.pwe.solver import PWESolver
from fibonacci_photonics.solvers.scan import (
    BAND_ABS_R_TOLERANCE,
    DispersionScan,
    FrequencyInterval,
    allowed_intervals,
    scan_from_semitrace,
)
from fibonacci_photonics.solvers.tmm.solver import TMMSolver

METHODS = ("tmm", "pwe")


def get_solver(method: str, spec: SuperlatticeSpec, **options: Any) -> TMMSolver | PWESolver:
    """Construye el método ``"tmm"`` o ``"pwe"`` para ``spec`` con sus opciones."""
    if method == "tmm":
        return TMMSolver(spec, **options)
    if method == "pwe":
        return PWESolver(spec, **options)
    raise InvalidParameterError(f"método desconocido: {method!r} (use 'tmm' o 'pwe')")


__all__ = [
    "BAND_ABS_R_TOLERANCE",
    "METHODS",
    "DispersionScan",
    "DispersionSolver",
    "FrequencyInterval",
    "PWESolver",
    "SemitraceFunction",
    "TMMSolver",
    "allowed_intervals",
    "get_solver",
    "scan_from_semitrace",
    "validate_angle",
    "validate_frequencies",
]
