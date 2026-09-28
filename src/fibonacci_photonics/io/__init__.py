"""Entrada/salida: rutas del proyecto, resultados, metadatos y tablas LaTeX."""

from fibonacci_photonics.io.paths import configs_dir, project_root, relative, result_dirs, results_dir
from fibonacci_photonics.io.results import dump_json, load_json, save_npz, write_readme, write_run_metadata

__all__ = [
    "configs_dir",
    "dump_json",
    "load_json",
    "project_root",
    "relative",
    "result_dirs",
    "results_dir",
    "save_npz",
    "write_readme",
    "write_run_metadata",
]
