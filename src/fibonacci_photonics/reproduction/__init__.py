"""Reproducción de las seis figuras de PRB 81, 153101 con la TMM o el PWE.

Uso::

    from fibonacci_photonics.reproduction import run_figure
    run_figure("figure_02", "pwe")

o desde la terminal: ``python scripts/run_figure.py figure_02 --method pwe``.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

from fibonacci_photonics.config import FIGURE_IDS, ConfigError
from fibonacci_photonics.reproduction import dispersion, plasmon
from fibonacci_photonics.reproduction.common import METHODS, FigureRun, Method

FIGURES: dict[str, Callable[[FigureRun], dict[str, Any]]] = {
    "figure_01": dispersion.figure_01,
    "figure_02": dispersion.mode_count_figure,
    "figure_03": dispersion.figure_03,
    "figure_04": dispersion.zoom_figure,
    "figure_05": plasmon.figure_05,
    "figure_06": plasmon.figure_06,
}


def run_figure(figure_id: str, method: str, *, workers: int | None = None) -> dict[str, Any]:
    """Calcula, dibuja y guarda una figura; devuelve su ``summary.json``.

    Parameters
    ----------
    figure_id : str
        ``"figure_01"`` … ``"figure_06"`` (se aceptan ``"1"`` o ``"01"``).
    method : {"tmm", "pwe"}
        Método que calcula R(ν).
    workers : int, optional
        Procesos del PWE; por defecto ``PWE_WORKERS`` o todos los núcleos.
    """
    figure_id = normalize_figure_id(figure_id)
    return FIGURES[figure_id](FigureRun.open(figure_id, method, workers=workers))


def normalize_figure_id(value: str) -> str:
    """``"2"``, ``"02"`` o ``"figure_02"`` → ``"figure_02"``."""
    text = str(value).removeprefix("figure_").removeprefix("fig")
    figure_id = f"figure_{int(text):02d}" if text.isdigit() else str(value)
    if figure_id not in FIGURES:
        raise ConfigError(f"figura desconocida {value!r}; opciones: {', '.join(FIGURE_IDS)}")
    return figure_id


__all__ = ["FIGURES", "METHODS", "FigureRun", "Method", "normalize_figure_id", "run_figure"]
