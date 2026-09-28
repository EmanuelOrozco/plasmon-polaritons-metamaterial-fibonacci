"""Relación de dispersión cos(k L_m) = R_m con el método de ondas planas."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field, replace
from typing import ClassVar

import numpy as np
from numpy.typing import ArrayLike, NDArray

from fibonacci_photonics.errors import InvalidParameterError
from fibonacci_photonics.physics.electromagnetics import Polarization, as_polarization
from fibonacci_photonics.physics.fibonacci import validate_order
from fibonacci_photonics.physics.superlattice import SuperlatticeSpec
from fibonacci_photonics.physics.units import omega_from_nu_ghz
from fibonacci_photonics.solvers.base import validate_angle, validate_frequencies, validate_positive_int
from fibonacci_photonics.solvers.pwe.bloch import FactorizationRule, n_max_for, solve_bloch, validate_rule
from fibonacci_photonics.solvers.pwe.convergence import (
    DEFAULT_REFINEMENT_FACTOR,
    DEFAULT_TOLERANCE,
    HarmonicConvergenceReport,
    check_harmonic_convergence,
)
from fibonacci_photonics.solvers.pwe.semitrace import SemitraceEvaluator
from fibonacci_photonics.solvers.scan import BAND_ABS_R_TOLERANCE, DispersionScan, scan_from_semitrace

DEFAULT_HARMONICS_PER_LAYER = 16

Harmonics = int | Mapping[int, int]
"""Armónicos por capa h: uno para todos los órdenes o uno por orden m."""


@dataclass(frozen=True)
class PWESolver:
    """Método de ondas planas en la forma k(ω).

    Attributes
    ----------
    spec : SuperlatticeSpec
        Superred.
    harmonics_per_layer : int or Mapping[int, int]
        h ≥ 1; el truncamiento es n_max = h · F_m, es decir N = 2 h F_m + 1
        ondas planas. Con un diccionario, h depende de m.
    rule : {"inverse", "laurent"}
        Regla de factorización de 1/χ; la inversa (Li) converge mucho más rápido.
    workers : int
        Procesos para los barridos.
    n_max_scale : float
        Factor sobre n_max (1.5 en el truncamiento de control).
    atol : float
        Tolerancia de |R| ≤ 1 + atol.

    Notes
    -----
    Costo por frecuencia O((2N)³) y memoria O((2N)²): con N ∝ F_m el costo
    crece como τ^{3m}. El error de truncamiento se estima con
    ``convergence_report``.
    """

    spec: SuperlatticeSpec
    harmonics_per_layer: Harmonics = DEFAULT_HARMONICS_PER_LAYER
    rule: FactorizationRule = "inverse"
    workers: int = 1
    n_max_scale: float = 1.0
    atol: float = BAND_ABS_R_TOLERANCE
    method: ClassVar[str] = "pwe"
    _harmonics: dict[int, int] | int = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        validate_rule(self.rule)
        validate_positive_int(self.workers, "workers")
        if not (np.isfinite(self.n_max_scale) and self.n_max_scale > 0.0):
            raise InvalidParameterError(f"n_max_scale debe ser > 0: {self.n_max_scale}")
        if isinstance(self.harmonics_per_layer, Mapping):
            value: dict[int, int] | int = {
                validate_order(int(m)): validate_positive_int(h, f"harmonics_per_layer[{m}]")
                for m, h in self.harmonics_per_layer.items()
            }
        else:
            value = validate_positive_int(self.harmonics_per_layer, "harmonics_per_layer")
        object.__setattr__(self, "_harmonics", value)

    @property
    def label(self) -> str:
        """Etiqueta del método en los barridos, p. ej. ``"pwe-inverse"``."""
        return f"pwe-{self.rule}"

    def harmonics(self, m: int) -> int:
        """Armónicos por capa h para el orden m."""
        if isinstance(self._harmonics, int):
            return self._harmonics
        try:
            return self._harmonics[validate_order(m)]
        except KeyError as exc:
            raise InvalidParameterError(f"harmonics_per_layer no define el orden m = {m}") from exc

    def n_max(self, m: int) -> int:
        """Truncamiento n_max = round(n_max_scale · h · F_m)."""
        base = n_max_for(m, self.harmonics(m))
        return base if self.n_max_scale == 1.0 else round(self.n_max_scale * base)

    def n_plane_waves(self, m: int) -> int:
        """N = 2 n_max + 1."""
        return 2 * self.n_max(m) + 1

    def refined(self, factor: float = DEFAULT_REFINEMENT_FACTOR) -> PWESolver:
        """El mismo método con ``factor`` veces más ondas planas (truncamiento de control)."""
        return replace(self, n_max_scale=self.n_max_scale * factor)

    def with_workers(self, workers: int) -> PWESolver:
        return replace(self, workers=workers)

    def semitrace(
        self, m: int, theta: float, nu_ghz: ArrayLike, polarization: Polarization | str = Polarization.TE
    ) -> NDArray[np.complex128]:
        """R_PWE(ν) = Re cos k̃ para frecuencias en GHz (NaN donde el problema es singular)."""
        m = validate_order(m)
        theta = validate_angle(theta)
        nu = validate_frequencies(nu_ghz)
        solution = solve_bloch(
            self.spec,
            m,
            omega_from_nu_ghz(nu.ravel()),
            theta,
            as_polarization(polarization),
            n_max=self.n_max(m),
            rule=self.rule,
            workers=self.workers,
        )
        return solution.r.reshape(nu.shape).astype(np.complex128)

    def scan(
        self, m: int, theta: float, nu_ghz: ArrayLike, polarization: Polarization | str = Polarization.TE
    ) -> DispersionScan:
        """R_PWE en una malla creciente (GHz), con el mismo contrato que la TMM."""
        nu = validate_frequencies(nu_ghz, grid=True)
        pol = as_polarization(polarization)
        r = self.semitrace(m, theta, nu, pol)
        return scan_from_semitrace(self.spec, m, theta, nu, r, pol, self.label, atol=self.atol)

    def evaluator(
        self, m: int, theta: float, polarization: Polarization | str = Polarization.TE
    ) -> SemitraceEvaluator:
        """R_PWE(ν) escalar con el mismo truncamiento y regla (para Brent)."""
        return SemitraceEvaluator(
            self.spec,
            validate_order(m),
            validate_angle(theta),
            as_polarization(polarization),
            n_max=self.n_max(m),
            rule=self.rule,
        )

    def convergence_report(
        self,
        m: int,
        theta: float,
        nu_ghz: ArrayLike,
        polarization: Polarization | str = Polarization.TE,
        *,
        factor: float = DEFAULT_REFINEMENT_FACTOR,
        tolerance: float = DEFAULT_TOLERANCE,
        warn: bool = True,
    ) -> HarmonicConvergenceReport:
        """Compara con ``factor`` veces más ondas planas; ver ``check_harmonic_convergence``."""
        return check_harmonic_convergence(
            self, m, theta, nu_ghz, polarization, factor=factor, tolerance=tolerance, warn=warn
        )


def pwe_semitrace(
    spec: SuperlatticeSpec,
    m: int,
    theta: float,
    nu_ghz: ArrayLike,
    polarization: Polarization | str = Polarization.TE,
    *,
    harmonics_per_layer: Harmonics = DEFAULT_HARMONICS_PER_LAYER,
    rule: FactorizationRule = "inverse",
    workers: int = 1,
) -> NDArray[np.float64]:
    """Atajo funcional: R_PWE(ν) real."""
    solver = PWESolver(spec, harmonics_per_layer, rule, workers)
    return solver.semitrace(m, theta, nu_ghz, polarization).real


def scan_dispersion(
    spec: SuperlatticeSpec,
    m: int,
    theta: float,
    nu_ghz: ArrayLike,
    polarization: Polarization | str = Polarization.TE,
    *,
    harmonics_per_layer: Harmonics = DEFAULT_HARMONICS_PER_LAYER,
    rule: FactorizationRule = "inverse",
    workers: int = 1,
    atol: float = BAND_ABS_R_TOLERANCE,
) -> DispersionScan:
    """Atajo funcional de ``PWESolver(...).scan(m, theta, nu_ghz, polarization)``."""
    solver = PWESolver(spec, harmonics_per_layer, rule, workers, atol=atol)
    return solver.scan(m, theta, nu_ghz, polarization)
