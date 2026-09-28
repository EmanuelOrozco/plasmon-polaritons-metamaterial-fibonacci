"""Especificación física de la superred de Fibonacci, en SI."""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass, field
from typing import Any

from fibonacci_photonics.errors import InvalidParameterError
from fibonacci_photonics.physics.constants import GHZ, SPEED_OF_LIGHT
from fibonacci_photonics.physics.materials import DrudeMetamaterial, HomogeneousMedium


@dataclass(frozen=True)
class SuperlatticeSpec:
    """Parámetros físicos de la superred A/B, validados al construirse.

    Attributes
    ----------
    thickness_a, thickness_b : float
        Espesores a y b de las capas A y B, en m (> 0).
    medium_a : HomogeneousMedium
        Medio A (aire en el paper).
    medium_b : DrudeMetamaterial
        Metamaterial B de Drude.
    speed_of_light : float
        c en m/s (> 0).
    source : str
        Procedencia de los parámetros (``"paper"``, figura inferida, ...).
    notes : tuple of str
        Notas de implementación que acompañan a los resultados.
    """

    thickness_a: float
    thickness_b: float
    medium_a: HomogeneousMedium
    medium_b: DrudeMetamaterial
    speed_of_light: float = SPEED_OF_LIGHT
    source: str = "paper"
    notes: tuple[str, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        for name in ("thickness_a", "thickness_b", "speed_of_light"):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0:
                raise InvalidParameterError(f"{name} debe ser finito y > 0: {value}")

    @property
    def omega_e(self) -> float:
        """Frecuencia de plasma eléctrica ω_e del medio B, en rad/s."""
        return self.medium_b.omega_e

    @property
    def omega_m(self) -> float:
        """Frecuencia de plasma magnética ω_m del medio B, en rad/s."""
        return self.medium_b.omega_m

    def nu_e_ghz(self) -> float:
        """ν_e en GHz, donde ε_B = 0."""
        return self.medium_b.electric_plasmon_frequency_hz() / GHZ

    def nu_m_ghz(self) -> float:
        """ν_m en GHz, donde μ_B = 0 (plasmón magnético)."""
        return self.medium_b.magnetic_plasmon_frequency_hz() / GHZ


def spec_as_dict(spec: SuperlatticeSpec) -> dict[str, Any]:
    """Diccionario serializable de la especificación (para metadatos)."""
    payload = asdict(spec)
    payload["medium_a"] = {
        "epsilon": complex(spec.medium_a.epsilon),
        "mu": complex(spec.medium_a.mu),
    }
    payload["medium_b"] = {
        "omega_e": spec.medium_b.omega_e,
        "omega_m": spec.medium_b.omega_m,
        "epsilon_0": spec.medium_b.epsilon_0,
        "mu_0": spec.medium_b.mu_0,
    }
    return payload
