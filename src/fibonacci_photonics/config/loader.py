"""Carga y validación de los YAML de ``configs/``.

- ``configs/figure_0N.yaml``: física de la superred más la definición de la figura.
- ``configs/pwe/figure_0N.yaml``: parámetros numéricos del PWE para esa figura.
- ``configs/pwe/<estudio>.yaml`` y ``configs/comparison/efficiency.yaml``: estudios.

Un YAML inválido lanza ``ConfigError`` con la ruta y los campos que fallan.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Generic, TypeVar

import yaml
from pydantic import BaseModel, ValidationError

from fibonacci_photonics.config.schema import (
    AngleSweepConfig,
    DispersionConfig,
    FigureConfig,
    ModeCountConfig,
    PWEAngleNumerics,
    PWEDispersionNumerics,
    PWEPlasmonNumerics,
    SuperlatticeConfig,
    WindowedConfig,
)
from fibonacci_photonics.errors import InvalidParameterError
from fibonacci_photonics.io.paths import configs_dir
from fibonacci_photonics.physics.superlattice import SuperlatticeSpec

Model = TypeVar("Model", bound=BaseModel)
FigureModel = TypeVar("FigureModel", bound=FigureConfig)

FIGURE_SCHEMAS: dict[str, type[FigureConfig]] = {
    "figure_01": DispersionConfig,
    "figure_02": ModeCountConfig,
    "figure_03": DispersionConfig,
    "figure_04": WindowedConfig,
    "figure_05": WindowedConfig,
    "figure_06": AngleSweepConfig,
}
"""Esquema de la definición física de cada figura."""

PWE_SCHEMAS: dict[str, type[BaseModel]] = {
    "figure_01": PWEDispersionNumerics,
    "figure_02": PWEDispersionNumerics,
    "figure_03": PWEDispersionNumerics,
    "figure_04": PWEDispersionNumerics,
    "figure_05": PWEPlasmonNumerics,
    "figure_06": PWEAngleNumerics,
}
"""Esquema de la numérica del PWE de cada figura."""

FIGURE_IDS = tuple(FIGURE_SCHEMAS)


class ConfigError(InvalidParameterError):
    """Archivo de configuración ausente, ilegible o con valores inválidos."""


@dataclass(frozen=True)
class LoadedConfig(Generic[Model]):
    """Modelo validado junto con el archivo del que salió."""

    model: Model
    path: Path


@dataclass(frozen=True)
class FigureSetup(Generic[FigureModel]):
    """Definición validada de una figura y su especificación física en SI."""

    figure_id: str
    config: FigureModel
    spec: SuperlatticeSpec
    path: Path


def resolve(path: str | Path) -> Path:
    """Ruta absoluta; las relativas se toman respecto de ``configs/``."""
    path = Path(path)
    return path if path.is_absolute() else configs_dir() / path


def load_yaml(path: str | Path) -> dict[str, Any]:
    """Contenido de un YAML como diccionario."""
    path = resolve(path)
    try:
        data = yaml.safe_load(path.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise ConfigError(f"no existe el archivo de configuración {path}") from exc
    except yaml.YAMLError as exc:
        raise ConfigError(f"YAML inválido en {path}: {exc}") from exc
    if not isinstance(data, dict):
        raise ConfigError(f"{path} debe contener un mapeo clave: valor")
    return data


def validate(model: type[Model], data: dict[str, Any], source: str | Path = "<dict>") -> Model:
    """``model.model_validate(data)`` con un mensaje que indica el origen."""
    try:
        return model.model_validate(data)
    except ValidationError as exc:
        raise ConfigError(f"configuración inválida en {source}:\n{exc}") from exc


def load_model(path: str | Path, model: type[Model]) -> LoadedConfig[Model]:
    """YAML validado con ``model``."""
    resolved = resolve(path)
    return LoadedConfig(validate(model, load_yaml(resolved), resolved), resolved)


def spec_from_mapping(data: dict[str, Any]) -> SuperlatticeSpec:
    """Especificación física a partir de un mapeo; ignora las claves que no son física.

    Espesores en mm y frecuencias de plasma como ω/2π en GHz; ver ``SuperlatticeConfig``.
    """
    physics = {key: value for key, value in data.items() if key in SuperlatticeConfig.model_fields}
    return validate(SuperlatticeConfig, physics).to_spec()


def load_spec(path: str | Path = "figure_01.yaml") -> SuperlatticeSpec:
    """Especificación física de un YAML de figura (las claves de figura se ignoran).

    Por defecto, la física del artículo (``configs/figure_01.yaml``: a = b = 12 mm,
    ω_e/2π = ω_m/2π = 3 GHz), común a todas las figuras.
    """
    return spec_from_mapping(load_yaml(path))


def load_figure(figure_id: str) -> FigureSetup[Any]:
    """Definición validada de ``configs/<figure_id>.yaml``."""
    if figure_id not in FIGURE_SCHEMAS:
        raise ConfigError(f"figura desconocida {figure_id!r}; opciones: {', '.join(FIGURE_IDS)}")
    loaded = load_model(f"{figure_id}.yaml", FIGURE_SCHEMAS[figure_id])
    config = loaded.model
    return FigureSetup(figure_id=figure_id, config=config, spec=config.to_spec(), path=loaded.path)


def load_pwe_numerics(figure_id: str) -> LoadedConfig[Any]:
    """Numérica validada de ``configs/pwe/<figure_id>.yaml``."""
    if figure_id not in PWE_SCHEMAS:
        raise ConfigError(f"figura desconocida {figure_id!r}; opciones: {', '.join(FIGURE_IDS)}")
    return load_model(Path("pwe") / f"{figure_id}.yaml", PWE_SCHEMAS[figure_id])
