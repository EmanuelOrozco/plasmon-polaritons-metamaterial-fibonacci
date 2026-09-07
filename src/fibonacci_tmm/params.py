"""Especificación física de la superred y carga de configuraciones YAML."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import yaml

from fibonacci_tmm.constants import GHZ, MM, SPEED_OF_LIGHT
from fibonacci_tmm.electromagnetics import Polarization
from fibonacci_tmm.materials import DrudeMetamaterial, HomogeneousMedium


def ghz_to_omega(nu_ghz: float) -> float:
    return 2.0 * 3.141592653589793 * nu_ghz * GHZ


def omega_to_ghz(omega: float) -> float:
    return omega / (2.0 * 3.141592653589793 * GHZ)


def mm_to_m(length_mm: float) -> float:
    return length_mm * MM


@dataclass(frozen=True)
class SuperlatticeSpec:
    """Todos los parámetros físicos en SI, más metadatos de procedencia."""

    thickness_a: float
    thickness_b: float
    medium_a: HomogeneousMedium
    medium_b: DrudeMetamaterial
    speed_of_light: float = SPEED_OF_LIGHT
    source: str = "paper"
    notes: tuple[str, ...] = field(default_factory=tuple)

    @property
    def omega_e(self) -> float:
        return self.medium_b.omega_e

    @property
    def omega_m(self) -> float:
        return self.medium_b.omega_m

    def nu_e_ghz(self) -> float:
        return self.medium_b.electric_plasmon_frequency_hz() / GHZ

    def nu_m_ghz(self) -> float:
        return self.medium_b.magnetic_plasmon_frequency_hz() / GHZ


def spec_from_mapping(data: dict[str, Any]) -> SuperlatticeSpec:
    """Interpreta un YAML científico. Frecuencias de plasma en GHz, espesores en mm."""
    notes = tuple(data.get("notes", ()))
    source = str(data.get("source", "paper"))
    c_choice = str(data.get("speed_of_light", "si"))
    if c_choice in {"si", "SI"}:
        c = SPEED_OF_LIGHT
    elif c_choice in {"rounded", "3e8"}:
        from fibonacci_tmm.constants import SPEED_OF_LIGHT_ROUNDED

        c = SPEED_OF_LIGHT_ROUNDED
    else:
        c = float(c_choice)

    medium_a = HomogeneousMedium(
        epsilon=complex(data.get("epsilon_a", 1.0)),
        mu=complex(data.get("mu_a", 1.0)),
        label="A",
    )
    medium_b = DrudeMetamaterial(
        omega_e=ghz_to_omega(float(data["omega_e_over_2pi_ghz"])),
        omega_m=ghz_to_omega(float(data["omega_m_over_2pi_ghz"])),
        epsilon_0=float(data.get("epsilon_0", 1.0)),
        mu_0=float(data.get("mu_0", 1.0)),
        label="B",
    )
    return SuperlatticeSpec(
        thickness_a=mm_to_m(float(data["layer_a_thickness_mm"])),
        thickness_b=mm_to_m(float(data["layer_b_thickness_mm"])),
        medium_a=medium_a,
        medium_b=medium_b,
        speed_of_light=c,
        source=source,
        notes=notes,
    )


def load_spec(path: str | Path) -> SuperlatticeSpec:
    with Path(path).open(encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    return spec_from_mapping(data)


def polarization_from_name(name: str) -> Polarization:
    return Polarization[name.strip().upper()]


def spec_as_dict(spec: SuperlatticeSpec) -> dict[str, Any]:
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
