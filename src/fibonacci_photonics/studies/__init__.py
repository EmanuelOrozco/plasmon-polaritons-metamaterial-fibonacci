"""Estudios numéricos complementarios a las figuras del paper, separados por método.

- ``tmm``: convergencia de la TMM respecto a la malla.
- ``pwe``: benchmark del libro, celda de Fourier, contaminación espectral y
  convergencia, costo y validez del PWE.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

from fibonacci_photonics.studies.pwe import book_benchmark, cell_fourier, spectral_pollution
from fibonacci_photonics.studies.pwe import convergence as pwe_convergence
from fibonacci_photonics.studies.tmm import convergence as tmm_convergence

STUDIES: dict[str, Callable[..., Any]] = {
    "tmm/convergence": tmm_convergence.main,
    "pwe/book_benchmark": book_benchmark.main,
    "pwe/cell_fourier": cell_fourier.main,
    "pwe/spectral_pollution": spectral_pollution.main,
    "pwe/convergence": pwe_convergence.main,
}
"""Estudio por nombre ``<método>/<estudio>`` (carpeta de ``results/`` y de ``studies/``)."""

__all__ = ["STUDIES"]
