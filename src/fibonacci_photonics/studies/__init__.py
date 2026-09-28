"""Estudios numéricos complementarios a las figuras del paper.

- ``tmm_convergence``: ancho de banda de la Fig. 4 frente a la malla (TMM).
- ``book_benchmark``: bicapa no dispersiva del libro de Sukhoivanov y Guryev.
- ``cell_fourier``: serie de Fourier del perfil escalonado de la celda S4.
- ``spectral_pollution``: autovalores espurios de la forma ω(k) con Drude.
- ``pwe_convergence``: convergencia, costo y validez del PWE frente a la TMM.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

from fibonacci_photonics.studies import (
    book_benchmark,
    cell_fourier,
    pwe_convergence,
    spectral_pollution,
    tmm_convergence,
)

STUDIES: dict[str, Callable[..., Any]] = {
    "tmm/convergence": tmm_convergence.main,
    "pwe/book_benchmark": book_benchmark.main,
    "pwe/cell_fourier": cell_fourier.main,
    "pwe/spectral_pollution": spectral_pollution.main,
    "pwe/convergence": pwe_convergence.main,
}
"""Estudio por nombre ``<método>/<estudio>`` (carpeta de ``results/``)."""

__all__ = ["STUDIES"]
