"""Esquemas Pydantic v2 de las configuraciones YAML.

Todos los modelos son inmutables y rechazan claves desconocidas
(``extra="forbid"``): una errata en un YAML falla al cargarlo, no en mitad de
una corrida de horas. Unidades: espesores en mm, frecuencias en GHz
(ν = ω/2π), ángulos en rad, espesores del cristal del libro en µm.
"""

from __future__ import annotations

import math
from typing import Annotated, Literal, Union

from pydantic import (
    AfterValidator,
    BaseModel,
    ConfigDict,
    Field,
    NonNegativeFloat,
    NonNegativeInt,
    PositiveFloat,
    PositiveInt,
    model_validator,
)

from fibonacci_photonics.physics.constants import SPEED_OF_LIGHT, SPEED_OF_LIGHT_ROUNDED
from fibonacci_photonics.physics.electromagnetics import Polarization
from fibonacci_photonics.physics.materials import DrudeMetamaterial, HomogeneousMedium
from fibonacci_photonics.physics.superlattice import SuperlatticeSpec
from fibonacci_photonics.physics.units import mm_to_m, omega_from_nu_ghz

ANGLE_SLACK = 1.0e-12


def _check_window(window: tuple[float, float]) -> tuple[float, float]:
    if not window[0] < window[1]:
        raise ValueError(f"la ventana debe cumplir min < max: {window}")
    return window


def _check_nonzero(value: float) -> float:
    if value == 0.0:
        raise ValueError("debe ser distinto de cero")
    return value


Angle = Annotated[float, Field(ge=0.0, le=math.pi / 2 + ANGLE_SLACK)]
"""Ángulo de incidencia θ ∈ [0, π/2] rad."""

Window = Annotated[tuple[PositiveFloat, PositiveFloat], AfterValidator(_check_window)]
"""Ventana (ν_min, ν_max) en GHz, con 0 < ν_min < ν_max."""

GridPoints = Annotated[int, Field(ge=2)]
FibonacciOrders = Annotated[list[NonNegativeInt], Field(min_length=1)]
NonZeroFloat = Annotated[float, AfterValidator(_check_nonzero)]
Rule = Literal["inverse", "laurent"]
Harmonics = Union[PositiveInt, dict[NonNegativeInt, PositiveInt]]
"""Armónicos por capa h (n_max = h F_m): un entero o un valor por orden m."""


class StrictModel(BaseModel):
    """Base inmutable que rechaza claves desconocidas."""

    model_config = ConfigDict(extra="forbid", frozen=True)


class SuperlatticeConfig(StrictModel):
    """Parámetros físicos de la superred (bloque común de ``configs/physics/figure_0N.yaml``)."""

    layer_a_thickness_mm: PositiveFloat
    layer_b_thickness_mm: PositiveFloat
    epsilon_a: NonZeroFloat = 1.0
    mu_a: NonZeroFloat = 1.0
    epsilon_0: PositiveFloat = 1.0
    mu_0: PositiveFloat = 1.0
    omega_e_over_2pi_ghz: NonNegativeFloat
    omega_m_over_2pi_ghz: NonNegativeFloat
    speed_of_light: Union[Literal["si", "SI", "rounded", "3e8"], PositiveFloat] = "si"
    source: str = "paper"
    notes: tuple[str, ...] = ()

    def speed_of_light_si(self) -> float:
        """c en m/s según la opción del YAML."""
        if self.speed_of_light in ("si", "SI"):
            return SPEED_OF_LIGHT
        if self.speed_of_light in ("rounded", "3e8"):
            return SPEED_OF_LIGHT_ROUNDED
        return float(self.speed_of_light)

    def to_spec(self) -> SuperlatticeSpec:
        """Especificación en SI (espesores en m, ω en rad/s)."""
        return SuperlatticeSpec(
            thickness_a=mm_to_m(self.layer_a_thickness_mm),
            thickness_b=mm_to_m(self.layer_b_thickness_mm),
            medium_a=HomogeneousMedium(epsilon=complex(self.epsilon_a), mu=complex(self.mu_a), label="A"),
            medium_b=DrudeMetamaterial(
                omega_e=float(omega_from_nu_ghz(self.omega_e_over_2pi_ghz)),
                omega_m=float(omega_from_nu_ghz(self.omega_m_over_2pi_ghz)),
                epsilon_0=self.epsilon_0,
                mu_0=self.mu_0,
                label="B",
            ),
            speed_of_light=self.speed_of_light_si(),
            source=self.source,
            notes=self.notes,
        )


class FigureConfig(SuperlatticeConfig):
    """Física más la polarización y los órdenes de Fibonacci de una figura."""

    polarization: Polarization = Polarization.TE
    fibonacci_orders: FibonacciOrders


class _LabelledAngles(FigureConfig):
    thetas_rad: Annotated[list[Angle], Field(min_length=1)]
    theta_labels: Annotated[list[str], Field(min_length=1)]

    @model_validator(mode="after")
    def _labels_match_angles(self) -> _LabelledAngles:
        if len(self.theta_labels) != len(self.thetas_rad):
            raise ValueError("theta_labels y thetas_rad deben tener la misma longitud")
        if len(set(self.theta_labels)) != len(self.theta_labels):
            raise ValueError("theta_labels no puede repetir etiquetas")
        return self

    def angles(self) -> list[tuple[str, float]]:
        """Pares (etiqueta, θ en rad) en el orden del YAML."""
        return list(zip(self.theta_labels, self.thetas_rad, strict=True))


class DispersionConfig(_LabelledAngles):
    """Figs. 1 y 3: paneles ν(k) para cada θ en una malla uniforme."""

    frequency_min_ghz: PositiveFloat
    frequency_max_ghz: PositiveFloat
    frequency_points: GridPoints

    @model_validator(mode="after")
    def _ordered_band(self) -> DispersionConfig:
        if not self.frequency_min_ghz < self.frequency_max_ghz:
            raise ValueError("frequency_min_ghz debe ser menor que frequency_max_ghz")
        return self


class ModeCountConfig(DispersionConfig):
    """Fig. 2: dispersión y conteo de subbandas bajo ν_m."""

    plasmon_window_ghz: Window


class WindowedConfig(_LabelledAngles):
    """Figs. 4 y 5: una ventana de frecuencia por ángulo."""

    frequency_windows_ghz: dict[str, Window]
    frequency_points: GridPoints

    @model_validator(mode="after")
    def _window_per_angle(self) -> WindowedConfig:
        missing = set(self.theta_labels) - set(self.frequency_windows_ghz)
        if missing:
            raise ValueError(f"falta frequency_windows_ghz para {sorted(missing)}")
        return self


class AngleSweepConfig(FigureConfig):
    """Fig. 6: anchos de banda frente a θ con malla logarítmica + uniforme."""

    theta_min_rad: Angle
    theta_max_rad: Angle
    theta_points: GridPoints
    plasmon_window_ghz: Window
    frequency_points_log: PositiveInt
    frequency_points_uniform: PositiveInt
    frequency_uniform_min_ghz: PositiveFloat
    nu_m_min_offset_ghz: PositiveFloat

    @model_validator(mode="after")
    def _ordered_angles(self) -> AngleSweepConfig:
        if not self.theta_min_rad < self.theta_max_rad:
            raise ValueError("theta_min_rad debe ser menor que theta_max_rad")
        return self


class _PWEBase(StrictModel):
    rule: Rule = "inverse"
    harmonics_per_layer: Harmonics

    def check_orders(self, orders: list[int]) -> None:
        """Error si ``harmonics_per_layer`` es un mapeo que no cubre ``orders``."""
        if isinstance(self.harmonics_per_layer, dict):
            missing = sorted(set(orders) - set(self.harmonics_per_layer))
            if missing:
                raise ValueError(f"harmonics_per_layer no define los órdenes {missing}")


class PWEDispersionNumerics(_PWEBase):
    """Numérica del PWE para las Figs. 1-4 (``configs/pwe/figure_0[1-4].yaml``)."""

    frequency_points: GridPoints
    close_band_edges: bool = True
    edge_xtol_fraction: PositiveFloat = 1.0e-6


class PWEPlasmonNumerics(_PWEBase):
    """Numérica del PWE para la Fig. 5: malla logarítmica en ν_m − ν y validación."""

    fibonacci_orders: FibonacciOrders
    grid_points: GridPoints
    nu_m_min_offset_ghz: PositiveFloat
    check_factor: Annotated[float, Field(gt=1.0)] = 1.5

    @model_validator(mode="after")
    def _orders_covered(self) -> PWEPlasmonNumerics:
        self.check_orders(self.fibonacci_orders)
        return self


class PWEAngleNumerics(PWEPlasmonNumerics):
    """Numérica del PWE para la Fig. 6: además, el barrido angular propio."""

    theta_min_rad: Angle
    theta_max_rad: Angle
    theta_points: GridPoints

    @model_validator(mode="after")
    def _ordered_angles(self) -> PWEAngleNumerics:
        if not self.theta_min_rad < self.theta_max_rad:
            raise ValueError("theta_min_rad debe ser menor que theta_max_rad")
        return self


class BookBenchmarkConfig(StrictModel):
    """Cristal no dispersivo del libro (Sukhoivanov y Guryev), espesores en µm."""

    thickness_1_um: PositiveFloat
    thickness_2_um: PositiveFloat
    epsilon_1: PositiveFloat
    epsilon_2: PositiveFloat
    n_max: PositiveInt
    k_points: GridPoints
    n_bands: PositiveInt
    field_k_over_pi: float
    field_band: NonNegativeInt
    synthesis_n_max: Annotated[list[PositiveInt], Field(min_length=1)]
    off_axis_q_max: PositiveFloat
    off_axis_q_points: GridPoints
    off_axis_n_max: PositiveInt

    @model_validator(mode="after")
    def _bands_fit(self) -> BookBenchmarkConfig:
        size = 2 * self.n_max + 1
        if self.n_bands > size or self.field_band >= size:
            raise ValueError(f"n_bands y field_band deben caber en 2 n_max + 1 = {size} ondas planas")
        return self


class BookConvergence(StrictModel):
    thickness_1_um: PositiveFloat
    thickness_2_um: PositiveFloat
    epsilon_1: PositiveFloat
    epsilon_2: PositiveFloat
    n_max_values: Annotated[list[PositiveInt], Field(min_length=1)]
    k_over_pi: Annotated[list[float], Field(min_length=1)]
    n_bands: PositiveInt


class DrudeConvergence(StrictModel):
    physics: str
    polarizations: Annotated[list[Polarization], Field(min_length=1)]
    theta_rad: Angle
    fibonacci_orders: FibonacciOrders
    harmonics_per_layer: Annotated[list[PositiveInt], Field(min_length=2)]
    frequency_min_ghz: PositiveFloat
    frequency_max_ghz: PositiveFloat
    frequency_points: GridPoints


class CostStudy(StrictModel):
    physics: str
    theta_rad: Angle
    window_ghz: Window
    frequency_points: GridPoints
    cases: Annotated[list[tuple[NonNegativeInt, PositiveInt]], Field(min_length=1)]


class NearPoleStudy(StrictModel):
    physics: str
    theta_rad: Angle
    fibonacci_order: NonNegativeInt
    harmonics_per_layer: PositiveInt
    offsets_ghz: Window
    points: GridPoints


class ConvergenceConfig(StrictModel):
    """Estudio de convergencia y costo del PWE (``configs/pwe/convergence.yaml``)."""

    book: BookConvergence
    drude: DrudeConvergence
    cost: CostStudy
    near_pole: NearPoleStudy


class SpectralPollutionConfig(StrictModel):
    """Forma ω(k) con regla de Laurent (``configs/pwe/spectral_pollution.yaml``)."""

    physics: str
    fibonacci_order: NonNegativeInt
    theta_rad: Angle
    k_over_pi: float
    n_max_values: Annotated[list[PositiveInt], Field(min_length=1)]
    window_ghz: Window


class ScalingStudy(StrictModel):
    theta_rad: Angle
    window_ghz: Window
    tmm_orders: FibonacciOrders
    tmm_frequencies: PositiveInt
    tmm_scalar_frequencies: PositiveInt
    pwe_orders: FibonacciOrders
    pwe_harmonics_per_layer: PositiveInt
    pwe_frequencies: dict[Union[int, str], PositiveInt]
    extrapolate_orders: list[NonNegativeInt]

    @model_validator(mode="after")
    def _default_count(self) -> ScalingStudy:
        if "default" not in self.pwe_frequencies:
            raise ValueError("pwe_frequencies necesita la clave 'default'")
        return self

    def pwe_frequency_count(self, m: int) -> int:
        return self.pwe_frequencies.get(m, self.pwe_frequencies["default"])


class WorkPrecisionStudy(StrictModel):
    fibonacci_order: NonNegativeInt
    theta_rad: Angle
    frequency_min_ghz: PositiveFloat
    frequency_max_ghz: PositiveFloat
    frequency_points: GridPoints
    harmonics_per_layer: Annotated[list[PositiveInt], Field(min_length=1)]


class FigureCostStudy(StrictModel):
    sample_frequencies: PositiveInt
    parallel_workers_used: PositiveInt


class BookCostStudy(StrictModel):
    config: str
    k_points: GridPoints
    grid_points: GridPoints


class EfficiencyConfig(StrictModel):
    """Estudio de eficiencia TMM frente a PWE (``configs/comparison/efficiency.yaml``)."""

    physics: str
    polarization: Polarization = Polarization.TE
    scaling: ScalingStudy
    work_precision: WorkPrecisionStudy
    figures: FigureCostStudy
    book: BookCostStudy
