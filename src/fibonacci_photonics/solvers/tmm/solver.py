"""Relación de dispersión cos(k L_m) = R_m con la matriz de transferencia."""

from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar, Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray

from fibonacci_photonics.errors import InvalidParameterError
from fibonacci_photonics.physics.electromagnetics import Polarization, as_polarization
from fibonacci_photonics.physics.fibonacci import validate_order
from fibonacci_photonics.physics.superlattice import SuperlatticeSpec
from fibonacci_photonics.physics.units import omega_from_nu_ghz
from fibonacci_photonics.solvers.base import validate_angle, validate_frequencies
from fibonacci_photonics.solvers.scan import BAND_ABS_R_TOLERANCE, DispersionScan, scan_from_semitrace
from fibonacci_photonics.solvers.tmm.transfer_matrix import semitrace_by_product, semitrace_by_recurrence

TMMAlgorithm = Literal["recurrence", "product"]


@dataclass
class TMMSemitraceEvaluator:
    """Re R_m(ν) escalar por la recurrencia, con contador de evaluaciones."""

    spec: SuperlatticeSpec
    m: int
    theta: float
    polarization: Polarization
    calls: int = 0

    def __call__(self, nu_ghz: float) -> float:
        self.calls += 1
        omega = omega_from_nu_ghz(np.array([nu_ghz]))
        return float(semitrace_by_recurrence(self.spec, self.m, omega, self.theta, self.polarization).real[0])


@dataclass(frozen=True)
class TMMSolver:
    """Método de la matriz de transferencia.

    Attributes
    ----------
    spec : SuperlatticeSpec
        Superred.
    algorithm : {"recurrence", "product"}
        ``"recurrence"`` usa R_m = 2 R_{m−1} R_{m−2} − R_{m−3}, O(m) por
        frecuencia; ``"product"`` multiplica las F_m matrices de capa y sirve
        de referencia independiente.
    atol : float
        Tolerancia de |R| ≤ 1 + atol para clasificar bandas.

    Notes
    -----
    Costo por frecuencia: O(m) con la recurrencia y O(F_m) con el producto;
    memoria O(1) por frecuencia. Es exacta salvo redondeo: no hay truncamiento.
    """

    spec: SuperlatticeSpec
    algorithm: TMMAlgorithm = "recurrence"
    atol: float = BAND_ABS_R_TOLERANCE
    method: ClassVar[str] = "tmm"

    def __post_init__(self) -> None:
        if self.algorithm not in ("recurrence", "product"):
            raise InvalidParameterError(f"algoritmo TMM desconocido: {self.algorithm!r}")

    def semitrace(
        self, m: int, theta: float, nu_ghz: ArrayLike, polarization: Polarization | str = Polarization.TE
    ) -> NDArray[np.complex128]:
        """R_m(ν) para frecuencias ν en GHz de cualquier forma."""
        m = validate_order(m)
        theta = validate_angle(theta)
        omega = omega_from_nu_ghz(validate_frequencies(nu_ghz))
        pol = as_polarization(polarization)
        if self.algorithm == "recurrence":
            return semitrace_by_recurrence(self.spec, m, omega, theta, pol)
        return semitrace_by_product(self.spec, m, omega, theta, pol)

    def scan(
        self, m: int, theta: float, nu_ghz: ArrayLike, polarization: Polarization | str = Polarization.TE
    ) -> DispersionScan:
        """R_m en una malla creciente ``nu_ghz`` (GHz), clasificada en bandas y gaps."""
        nu = validate_frequencies(nu_ghz, grid=True)
        pol = as_polarization(polarization)
        r = self.semitrace(m, theta, nu, pol)
        return scan_from_semitrace(self.spec, m, theta, nu, r, pol, self.algorithm, atol=self.atol)

    def evaluator(
        self, m: int, theta: float, polarization: Polarization | str = Polarization.TE
    ) -> TMMSemitraceEvaluator:
        """Re R_m(ν) escalar (recurrencia) para refinar bordes."""
        return TMMSemitraceEvaluator(self.spec, validate_order(m), validate_angle(theta), as_polarization(polarization))


def scan_dispersion(
    spec: SuperlatticeSpec,
    m: int,
    theta: float,
    nu_ghz: ArrayLike,
    polarization: Polarization | str = Polarization.TE,
    *,
    method: TMMAlgorithm = "recurrence",
    atol: float = BAND_ABS_R_TOLERANCE,
) -> DispersionScan:
    """Atajo funcional de ``TMMSolver(spec, method, atol).scan(m, theta, nu_ghz, polarization)``."""
    return TMMSolver(spec, algorithm=method, atol=atol).scan(m, theta, nu_ghz, polarization)
