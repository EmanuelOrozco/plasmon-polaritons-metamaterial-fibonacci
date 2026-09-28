"""Semitraza del PWE, R_PWE(ν) = Re cos(k̃), con k̃ el vector de Bloch físico."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from fibonacci_photonics.physics.electromagnetics import Polarization
from fibonacci_photonics.physics.superlattice import SuperlatticeSpec
from fibonacci_photonics.physics.units import omega_from_nu_ghz
from fibonacci_photonics.solvers.pwe.bloch import (
    FactorizationRule,
    bloch_eigenvalues,
    select_bloch_wavevector,
    validate_rule,
)
from fibonacci_photonics.solvers.pwe.fourier import cell_fourier


class SemitraceEvaluator:
    """R_PWE(ν) escalar con la celda de Fourier precalculada (para Brent).

    Parameters
    ----------
    spec : SuperlatticeSpec
        Superred.
    m : int
        Orden de Fibonacci.
    theta : float
        Ángulo de incidencia en rad.
    polarization : Polarization
        TE o TM.
    n_max : int
        Truncamiento n = −n_max … n_max.
    rule : {"inverse", "laurent"}
        Regla de factorización de 1/χ.

    Attributes
    ----------
    calls : int
        Evaluaciones hechas (cada una resuelve un problema de 2N × 2N).
    """

    def __init__(
        self,
        spec: SuperlatticeSpec,
        m: int,
        theta: float,
        polarization: Polarization = Polarization.TE,
        *,
        n_max: int,
        rule: FactorizationRule = "inverse",
    ) -> None:
        self.spec = spec
        self.theta = theta
        self.polarization = polarization
        self.rule = validate_rule(rule)
        self.cell = cell_fourier(m, spec.thickness_a, spec.thickness_b, n_max)
        self.calls = 0

    def __call__(self, nu_ghz: float) -> float:
        """Re cos k̃ a la frecuencia ν (GHz); NaN si el problema es singular."""
        self.calls += 1
        omega = float(omega_from_nu_ghz(nu_ghz))
        try:
            eigenvalues = bloch_eigenvalues(self.spec, self.cell, omega, self.theta, self.polarization, self.rule)
        except np.linalg.LinAlgError:
            return float("nan")
        eigenvalues = eigenvalues[np.isfinite(eigenvalues)]
        if eigenvalues.size == 0:
            return float("nan")
        return float(np.cos(select_bloch_wavevector(eigenvalues)).real)

    def many(self, nu_ghz: ArrayLike) -> NDArray[np.float64]:
        """R_PWE en varias frecuencias, en serie."""
        return np.array([self(float(nu)) for nu in np.atleast_1d(nu_ghz)], dtype=np.float64)
