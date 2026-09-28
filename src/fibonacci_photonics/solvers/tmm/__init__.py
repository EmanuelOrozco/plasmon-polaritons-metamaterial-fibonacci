"""Método de la matriz de transferencia (TMM).

- ``transfer_matrix``: matrices de capa, producto de celda y recurrencia de trazas.
- ``solver``: ``TMMSolver``, que implementa ``DispersionSolver`` en GHz.
"""

from fibonacci_photonics.solvers.tmm.solver import TMMSemitraceEvaluator, TMMSolver, scan_dispersion
from fibonacci_photonics.solvers.tmm.transfer_matrix import (
    analytic_r0_r1_r2,
    cell_transfer_matrix,
    semitrace_by_product,
    semitrace_by_recurrence,
)

__all__ = [
    "TMMSemitraceEvaluator",
    "TMMSolver",
    "analytic_r0_r1_r2",
    "cell_transfer_matrix",
    "scan_dispersion",
    "semitrace_by_product",
    "semitrace_by_recurrence",
]
