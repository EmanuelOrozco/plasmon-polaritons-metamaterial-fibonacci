"""Método de expansión en ondas planas (PWE).

- ``fourier``: coeficientes analíticos de perfiles escalonados en la celda S_m.
- ``bloch``: forma k(ω) con la regla inversa de factorización (método principal).
- ``semitrace``: R_PWE(ν) escalar para Brent.
- ``convergence``: verificación del truncamiento con más ondas planas.
- ``solver``: ``PWESolver``, que implementa ``DispersionSolver`` en GHz.
- ``eigenfrequency``: forma ω(k) del libro (benchmark no dispersivo y problema de Drude).
"""

from fibonacci_photonics.solvers.pwe.bloch import BlochSolution, n_max_for, select_bloch_wavevector, solve_bloch
from fibonacci_photonics.solvers.pwe.convergence import HarmonicConvergenceReport, check_harmonic_convergence
from fibonacci_photonics.solvers.pwe.semitrace import SemitraceEvaluator
from fibonacci_photonics.solvers.pwe.solver import (
    DEFAULT_HARMONICS_PER_LAYER,
    PWESolver,
    pwe_semitrace,
    scan_dispersion,
)

__all__ = [
    "DEFAULT_HARMONICS_PER_LAYER",
    "BlochSolution",
    "HarmonicConvergenceReport",
    "PWESolver",
    "SemitraceEvaluator",
    "check_harmonic_convergence",
    "n_max_for",
    "pwe_semitrace",
    "scan_dispersion",
    "select_bloch_wavevector",
    "solve_bloch",
]
