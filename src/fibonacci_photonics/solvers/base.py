"""Contrato común de los métodos numéricos y validación de sus entradas.

Un ``DispersionSolver`` calcula la semitraza R_m(ν) = cos(k L_m) de la celda
S_m para una superred fija. Todo el post-proceso (clasificación banda/gap,
bordes, subbandas, comparación) se escribe contra este contrato, sin saber si
R viene de la matriz de transferencia o de las ondas planas.
"""

from __future__ import annotations

import math
from typing import Protocol, runtime_checkable

import numpy as np
from numpy.typing import ArrayLike, NDArray

from fibonacci_photonics.errors import InvalidParameterError
from fibonacci_photonics.physics.electromagnetics import Polarization
from fibonacci_photonics.physics.superlattice import SuperlatticeSpec
from fibonacci_photonics.solvers.scan import DispersionScan

MAX_ANGLE = 0.5 * math.pi
"""Incidencia rasante θ = π/2 es el límite físico del ángulo en el medio A."""


class SemitraceFunction(Protocol):
    """R(ν) escalar, real, con contador de evaluaciones (para Brent)."""

    calls: int

    def __call__(self, nu_ghz: float) -> float: ...


@runtime_checkable
class DispersionSolver(Protocol):
    """Método que calcula R_m(ν) para una superred ``spec``.

    Attributes
    ----------
    method : str
        ``"tmm"`` o ``"pwe"``.
    spec : SuperlatticeSpec
        Superred que se resuelve.
    """

    @property
    def method(self) -> str: ...

    @property
    def spec(self) -> SuperlatticeSpec: ...

    def semitrace(
        self, m: int, theta: float, nu_ghz: ArrayLike, polarization: Polarization | str = ...
    ) -> NDArray[np.complex128]:
        """R_m(ν) vectorizada sobre ``nu_ghz`` (GHz)."""
        ...

    def scan(
        self, m: int, theta: float, nu_ghz: ArrayLike, polarization: Polarization | str = ...
    ) -> DispersionScan:
        """R_m en una malla creciente de frecuencias, clasificada en bandas y gaps."""
        ...

    def evaluator(self, m: int, theta: float, polarization: Polarization | str = ...) -> SemitraceFunction:
        """Re R_m(ν) escalar para refinar bordes con Brent."""
        ...


def validate_angle(theta: float) -> float:
    """Ángulo de incidencia θ ∈ [0, π/2] en rad.

    Raises
    ------
    InvalidParameterError
        Si θ no es finito o cae fuera de [0, π/2].
    """
    value = float(theta)
    if not math.isfinite(value) or value < 0.0 or value > MAX_ANGLE + 1.0e-12:
        raise InvalidParameterError(f"θ debe estar en [0, π/2] rad: {theta!r}")
    return value


def validate_frequencies(nu_ghz: ArrayLike, *, grid: bool = False) -> NDArray[np.float64]:
    """Frecuencias en GHz finitas y positivas.

    Parameters
    ----------
    nu_ghz : array_like
        Frecuencias ν en GHz; con ``grid=False`` se acepta cualquier forma.
    grid : bool
        Exige además una malla 1-D, no vacía y estrictamente creciente.

    Raises
    ------
    InvalidParameterError
        Si hay valores no finitos o ≤ 0 (el modelo de Drude diverge en ν = 0),
        o si ``grid`` y la malla no es 1-D creciente.
    """
    nu = np.asarray(nu_ghz, dtype=np.float64)
    if nu.size == 0:
        raise InvalidParameterError("la malla de frecuencias está vacía")
    if not np.all(np.isfinite(nu)):
        raise InvalidParameterError("las frecuencias deben ser finitas")
    if np.any(nu <= 0.0):
        raise InvalidParameterError(f"las frecuencias deben ser > 0 GHz (mínimo {float(nu.min())!r})")
    if grid:
        if nu.ndim != 1:
            raise InvalidParameterError(f"la malla de frecuencias debe ser 1-D, no de forma {nu.shape}")
        if nu.size > 1 and np.any(np.diff(nu) <= 0.0):
            raise InvalidParameterError("la malla de frecuencias debe ser estrictamente creciente")
    return nu


def validate_positive_int(value: int, name: str) -> int:
    """Entero ≥ 1 (armónicos por capa, procesos, puntos de malla)."""
    if isinstance(value, bool) or int(value) != value or int(value) < 1:
        raise InvalidParameterError(f"{name} debe ser un entero ≥ 1: {value!r}")
    return int(value)
