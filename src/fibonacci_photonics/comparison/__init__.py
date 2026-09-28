"""Comparación de los métodos TMM y PWE (salidas en ``results/comparison``).

- ``compare``: acuerdo punto a punto entre barridos y errores de borde entre bandas.
- ``timing``: tiempos, ajustes de potencias y ``benchmark_methods`` sobre cualquier
  conjunto de ``DispersionSolver``.
- ``dispersion_report``, ``plasmon_report``, ``efficiency`` y ``summary``: los
  informes de ``results/comparison`` (figuras, JSON y tablas LaTeX).
"""

from fibonacci_photonics.comparison.compare import ScanAgreement, compare_scans, edge_errors, match_intervals
from fibonacci_photonics.comparison.timing import MethodRun, benchmark_methods, power_fit, stopwatch, timed

__all__ = [
    "MethodRun",
    "ScanAgreement",
    "benchmark_methods",
    "compare_scans",
    "edge_errors",
    "match_intervals",
    "power_fit",
    "stopwatch",
    "timed",
]
