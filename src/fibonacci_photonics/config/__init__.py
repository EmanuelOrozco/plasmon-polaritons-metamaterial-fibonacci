"""Configuración tipada (Pydantic v2) de figuras, numérica del PWE y estudios."""

from fibonacci_photonics.config.loader import (
    FIGURE_IDS,
    ConfigError,
    FigureSetup,
    LoadedConfig,
    load_figure,
    load_model,
    load_pwe_numerics,
    load_spec,
    load_yaml,
    spec_from_mapping,
)
from fibonacci_photonics.config.schema import SuperlatticeConfig

__all__ = [
    "FIGURE_IDS",
    "ConfigError",
    "FigureSetup",
    "LoadedConfig",
    "SuperlatticeConfig",
    "load_figure",
    "load_model",
    "load_pwe_numerics",
    "load_spec",
    "load_yaml",
    "spec_from_mapping",
]
