"""Contexto compartido por los pipelines de las figuras: configuración, solver y salidas."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

from matplotlib.figure import Figure

from fibonacci_photonics.analysis.band_closure import DEFAULT_EDGE_XTOL_FRACTION
from fibonacci_photonics.config import FigureSetup, load_figure, load_pwe_numerics
from fibonacci_photonics.errors import InvalidParameterError
from fibonacci_photonics.io.paths import relative, result_dirs
from fibonacci_photonics.io.results import dump_json, write_readme, write_run_metadata
from fibonacci_photonics.parallel import single_threaded_blas, worker_count
from fibonacci_photonics.physics.electromagnetics import Polarization
from fibonacci_photonics.physics.superlattice import SuperlatticeSpec
from fibonacci_photonics.solvers.base import DispersionSolver
from fibonacci_photonics.solvers.pwe import PWESolver
from fibonacci_photonics.solvers.tmm import TMMSolver
from fibonacci_photonics.viz.style import PWE_FORMATS, TMM_FORMATS, save_figure

Method = Literal["tmm", "pwe"]
METHODS: tuple[Method, ...] = ("tmm", "pwe")
METHOD_NAMES = {"tmm": "TMM", "pwe": "PWE"}


@dataclass(frozen=True)
class FigureRun:
    """Una figura del paper reproducida con un método.

    Attributes
    ----------
    figure_id : str
        ``"figure_01"`` … ``"figure_06"``.
    method : {"tmm", "pwe"}
        Método que calcula R(ν).
    setup : FigureSetup
        Física (``configs/physics/<figura>.yaml``) y mallas de la TMM (``configs/tmm/<figura>.yaml``).
    numerics : pydantic model or None
        Numérica del PWE (``configs/pwe/<figura>.yaml``); ``None`` en la TMM.
    numerics_path : Path or None
        Archivo de ``numerics``.
    dirs : dict
        ``root``, ``data`` y ``output`` de ``results/<método>/<figura>``.
    workers : int
        Procesos para los barridos del PWE.
    """

    figure_id: str
    method: Method
    setup: FigureSetup[Any]
    numerics: Any
    numerics_path: Path | None
    dirs: dict[str, Path]
    workers: int

    @classmethod
    def open(cls, figure_id: str, method: str, *, workers: int | None = None) -> FigureRun:
        """Carga y valida las configuraciones y crea las carpetas de salida."""
        if method not in METHODS:
            raise InvalidParameterError(f"método desconocido {method!r}; opciones: {', '.join(METHODS)}")
        setup = load_figure(figure_id)
        numerics, numerics_path = None, None
        if method == "pwe":
            single_threaded_blas()
            loaded = load_pwe_numerics(figure_id)
            numerics, numerics_path = loaded.model, loaded.path
        return cls(
            figure_id=figure_id,
            method=method,
            setup=setup,
            numerics=numerics,
            numerics_path=numerics_path,
            dirs=result_dirs(method, figure_id),
            workers=worker_count(workers),
        )

    @property
    def spec(self) -> SuperlatticeSpec:
        return self.setup.spec

    @property
    def config(self) -> Any:
        return self.setup.config

    @property
    def polarization(self) -> Polarization:
        return Polarization(self.config.polarization)

    @property
    def is_pwe(self) -> bool:
        return self.method == "pwe"

    @property
    def name(self) -> str:
        return METHOD_NAMES[self.method]

    @property
    def formats(self) -> tuple[str, ...]:
        return PWE_FORMATS if self.is_pwe else TMM_FORMATS

    @property
    def close_band_edges(self) -> bool:
        """Cierre de las curvas en los bordes exactos: siempre en la TMM, configurable en el PWE.

        Aun con la malla densa de la TMM, las puntas de las bandas planas junto a ν_m
        caen entre muestras; refinarlas cuesta milisegundos.
        """
        if not self.is_pwe:
            return True
        return bool(getattr(self.numerics, "close_band_edges", False))

    @property
    def edge_xtol_fraction(self) -> float:
        if self.is_pwe:
            return float(self.numerics.edge_xtol_fraction)
        return DEFAULT_EDGE_XTOL_FRACTION

    def solver(self, *, workers: int | None = None) -> DispersionSolver:
        """Solver configurado para esta figura; ``workers`` solo afecta al PWE."""
        if not self.is_pwe:
            return TMMSolver(self.spec)
        return PWESolver(
            self.spec,
            harmonics_per_layer=self.numerics.harmonics_per_layer,
            rule=self.numerics.rule,
            workers=self.workers if workers is None else workers,
        )

    def finish(self, fig: Figure, payload: dict[str, Any], *, title: str, description: str) -> dict[str, Any]:
        """Guarda la figura, ``summary.json``, ``metadata.json`` y el README de la carpeta."""
        paths = save_figure(fig, self.dirs["output"], self.figure_id, formats=self.formats)
        config_path = self.numerics_path if self.is_pwe and self.numerics_path else self.setup.tmm_path
        extra: dict[str, Any] = {
            "config": relative(config_path),
            "physics": relative(self.setup.path),
            "method": self.method,
        }
        if self.is_pwe:
            extra["numerics"] = self.numerics.model_dump(mode="json")
        extra |= payload
        extra["outputs"] = {key: relative(path) for key, path in paths.items()}
        write_run_metadata(self.dirs["data"] / "metadata.json", self.spec, extra)
        dump_json(self.dirs["data"] / "summary.json", extra)
        suffix = " (PWE)" if self.is_pwe else ""
        write_readme(
            self.dirs["root"] / "README.md",
            f"{title}{suffix}",
            f"{description}\n\nScript: `python scripts/run_figure.py {self.figure_id} --method {self.method}`",
        )
        print(f"{title} ({self.name}) escrita en {relative(self.dirs['output'])}")
        return extra
