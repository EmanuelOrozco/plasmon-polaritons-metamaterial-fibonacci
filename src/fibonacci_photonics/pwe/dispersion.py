"""Relación de dispersión cos(k Lm) = Rm obtenida con el método de ondas planas."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike

from fibonacci_photonics.core.bands import DispersionScan, omega_from_nu_ghz, scan_from_semitrace
from fibonacci_photonics.core.constants import BAND_ABS_R_TOLERANCE
from fibonacci_photonics.core.electromagnetics import Polarization
from fibonacci_photonics.core.params import SuperlatticeSpec
from fibonacci_photonics.pwe.bloch import (
    FactorizationRule,
    bloch_eigenvalues,
    n_max_for,
    select_bloch_wavevector,
    solve_bloch,
)
from fibonacci_photonics.pwe.fourier import cell_fourier

DEFAULT_HARMONICS_PER_LAYER = 16


class SemitraceEvaluator:
    """R_PWE(ν) escalar con la celda de Fourier precalculada (para Brent)."""

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
        self.rule = rule
        self.cell = cell_fourier(m, spec.thickness_a, spec.thickness_b, n_max)
        self.calls = 0

    def __call__(self, nu_ghz: float) -> float:
        self.calls += 1
        omega = float(omega_from_nu_ghz(nu_ghz))
        try:
            eigenvalues = bloch_eigenvalues(
                self.spec, self.cell, omega, self.theta, self.polarization, self.rule
            )
        except np.linalg.LinAlgError:
            return float("nan")
        eigenvalues = eigenvalues[np.isfinite(eigenvalues)]
        if eigenvalues.size == 0:
            return float("nan")
        return float(np.cos(select_bloch_wavevector(eigenvalues)).real)

    def many(self, nu_ghz: ArrayLike) -> np.ndarray:
        return np.array([self(float(nu)) for nu in np.atleast_1d(nu_ghz)])


def pwe_semitrace(
    spec: SuperlatticeSpec,
    m: int,
    theta: float,
    nu_ghz: ArrayLike,
    polarization: Polarization = Polarization.TE,
    *,
    harmonics_per_layer: int = DEFAULT_HARMONICS_PER_LAYER,
    n_max: int | None = None,
    rule: FactorizationRule = "inverse",
    workers: int = 1,
) -> np.ndarray:
    """R_PWE(ν) = Re cos(k̃), con k̃ el vector de Bloch físico del QEP."""
    omega = omega_from_nu_ghz(np.asarray(nu_ghz, dtype=np.float64))
    n_max = n_max_for(m, harmonics_per_layer) if n_max is None else int(n_max)
    return solve_bloch(
        spec, m, omega, theta, polarization, n_max=n_max, rule=rule, workers=workers
    ).r


def scan_dispersion(
    spec: SuperlatticeSpec,
    m: int,
    theta: float,
    nu_ghz: ArrayLike,
    polarization: Polarization = Polarization.TE,
    *,
    harmonics_per_layer: int = DEFAULT_HARMONICS_PER_LAYER,
    n_max: int | None = None,
    rule: FactorizationRule = "inverse",
    workers: int = 1,
    atol: float = BAND_ABS_R_TOLERANCE,
) -> DispersionScan:
    """Mismo contrato que ``fibonacci_photonics.tmm.dispersion.scan_dispersion``."""
    nu = np.asarray(nu_ghz, dtype=np.float64)
    r = pwe_semitrace(
        spec,
        m,
        theta,
        nu,
        polarization,
        harmonics_per_layer=harmonics_per_layer,
        n_max=n_max,
        rule=rule,
        workers=workers,
    )
    return scan_from_semitrace(spec, m, theta, nu, r, polarization, f"pwe-{rule}", atol=atol)
