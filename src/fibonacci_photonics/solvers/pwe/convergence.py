"""Convergencia del PWE respecto al número de ondas planas.

El PWE trunca la serie de Fourier en N = 2 n_max + 1 términos; el error en R
decae algebraicamente con N (los coeficientes del perfil escalonado caen como
1/|n|). La prueba práctica es repetir el cálculo con más ondas planas: si R no
cambia dentro de la tolerancia en las bandas y la clasificación banda/gap es la
misma, el truncamiento es suficiente. Cerca de ν_m, donde |Im k L_m| es grande,
el PWE deja de converger y aquí se detecta.
"""

from __future__ import annotations

import warnings
from dataclasses import asdict, dataclass
from typing import TYPE_CHECKING, Any

import numpy as np
from numpy.typing import ArrayLike

from fibonacci_photonics.errors import ConvergenceWarning, InvalidParameterError
from fibonacci_photonics.physics.electromagnetics import Polarization

if TYPE_CHECKING:
    from fibonacci_photonics.solvers.pwe.solver import PWESolver

DEFAULT_REFINEMENT_FACTOR = 1.5
"""Las bandas se validan con 1.5 veces más ondas planas."""

DEFAULT_TOLERANCE = 1.0e-2
"""Diferencia máxima admitida en R entre los dos truncamientos."""


@dataclass(frozen=True)
class HarmonicConvergenceReport:
    """Resultado de comparar dos truncamientos N y N' ≈ factor · N.

    Attributes
    ----------
    m : int
        Orden de Fibonacci.
    n_max, n_max_check : int
        Truncamientos comparados.
    max_abs_delta_r_in_bands : float
        max |R_N − R_N'| donde |R_N'| ≤ 1.
    classification_mismatch : int
        Frecuencias donde los dos truncamientos no coinciden en banda/gap.
    tolerance : float
        Tolerancia usada.
    """

    m: int
    n_max: int
    n_max_check: int
    max_abs_delta_r_in_bands: float
    classification_mismatch: int
    tolerance: float

    @property
    def converged(self) -> bool:
        """R convergido en bandas y misma clasificación banda/gap."""
        return self.max_abs_delta_r_in_bands <= self.tolerance and self.classification_mismatch == 0

    def as_dict(self) -> dict[str, Any]:
        return {**asdict(self), "converged": self.converged}


def check_harmonic_convergence(
    solver: PWESolver,
    m: int,
    theta: float,
    nu_ghz: ArrayLike,
    polarization: Polarization | str = Polarization.TE,
    *,
    factor: float = DEFAULT_REFINEMENT_FACTOR,
    tolerance: float = DEFAULT_TOLERANCE,
    warn: bool = True,
) -> HarmonicConvergenceReport:
    """Compara R_PWE con ``solver`` y con ``factor`` veces más ondas planas.

    Parameters
    ----------
    solver : PWESolver
        Configuración a verificar.
    m, theta, nu_ghz, polarization
        Caso a evaluar; conviene una muestra pequeña de frecuencias.
    factor : float
        Razón N'/N (> 1).
    tolerance : float
        Diferencia máxima admitida en R dentro de las bandas.
    warn : bool
        Emite ``ConvergenceWarning`` si no hay convergencia.
    """
    if not factor > 1.0:
        raise InvalidParameterError(f"el factor de refinamiento debe ser > 1: {factor}")
    check = solver.refined(factor)
    r = solver.semitrace(m, theta, nu_ghz, polarization).real
    r_check = check.semitrace(m, theta, nu_ghz, polarization).real
    band = np.abs(r_check) <= 1.0
    finite = np.isfinite(r) & np.isfinite(r_check)
    delta = np.abs(r - r_check)[band & finite]
    report = HarmonicConvergenceReport(
        m=m,
        n_max=solver.n_max(m),
        n_max_check=check.n_max(m),
        max_abs_delta_r_in_bands=float(delta.max()) if delta.size else 0.0,
        classification_mismatch=int(np.count_nonzero(((np.abs(r) <= 1.0) != band) & finite)),
        tolerance=tolerance,
    )
    if warn and not report.converged:
        warnings.warn(
            f"PWE sin convergir en m={m}: max|ΔR| = {report.max_abs_delta_r_in_bands:.2e} "
            f"(tolerancia {tolerance:.1e}), {report.classification_mismatch} frecuencias cambian "
            f"de banda a gap entre n_max = {report.n_max} y {report.n_max_check}",
            ConvergenceWarning,
            stacklevel=2,
        )
    return report
