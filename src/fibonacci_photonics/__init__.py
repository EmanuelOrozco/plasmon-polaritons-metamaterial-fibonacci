"""Plasmon-polaritones en superredes de Fibonacci (Reyes-Gómez et al., PRB 81, 153101, 2010).

Dos métodos independientes calculan la misma semitraza R_m(ν) = cos(k L_m):

- ``solvers.tmm``: matriz de transferencia (``TMMSolver``), exacta por capas.
- ``solvers.pwe``: ondas planas en la forma k(ω) con la regla inversa (``PWESolver``).

Capas del paquete:

- ``physics``: materiales, polarización, palabras de Fibonacci y ``SuperlatticeSpec``.
- ``solvers``: protocolo ``DispersionSolver``, barridos y ambos métodos.
- ``analysis``: bordes de banda, subbandas plasmónicas y anchos de banda.
- ``config``: esquemas Pydantic de los YAML de ``configs/`` (``physics/``, ``tmm/``, ``pwe/``).
- ``comparison``: comparación y tiempos TMM vs PWE.
- ``reproduction``: figuras del artículo con cualquiera de los dos métodos.
- ``studies.tmm`` y ``studies.pwe``: estudios numéricos de cada método.
- ``io`` y ``viz``: rutas, resultados, tablas LaTeX y gráficas.

Ejemplo::

    import numpy as np
    from fibonacci_photonics import TMMSolver, load_spec

    solver = TMMSolver(load_spec())
    scan = solver.scan(m=4, theta=0.0, nu_ghz=np.linspace(0.5, 3.0, 2001), polarization="TE")
"""

from fibonacci_photonics.analysis.bandwidth import BandwidthPoint, bandwidth_versus_angle
from fibonacci_photonics.analysis.plasmon_modes import PlasmonModeReport, detect_plasmon_modes
from fibonacci_photonics.config import load_figure, load_spec
from fibonacci_photonics.errors import ConvergenceWarning, InvalidParameterError
from fibonacci_photonics.physics.constants import GOLDEN_RATIO, SPEED_OF_LIGHT
from fibonacci_photonics.physics.effective_medium import (
    average_epsilon_mu,
    zero_average_index_frequency_ghz,
)
from fibonacci_photonics.physics.electromagnetics import Polarization
from fibonacci_photonics.physics.fibonacci import (
    FibonacciCell,
    build_cell,
    cell_length,
    n_layers_a,
    n_layers_b,
    paper_fibonacci,
    sequence,
)
from fibonacci_photonics.physics.superlattice import SuperlatticeSpec
from fibonacci_photonics.solvers import (
    BAND_ABS_R_TOLERANCE,
    DispersionScan,
    DispersionSolver,
    FrequencyInterval,
    PWESolver,
    TMMSolver,
    allowed_intervals,
    get_solver,
)
from fibonacci_photonics.solvers.tmm.transfer_matrix import (
    analytic_r0_r1_r2,
    cell_transfer_matrix,
    semitrace_by_product,
    semitrace_by_recurrence,
)

__version__ = "2.0.0"

__all__ = [
    "BAND_ABS_R_TOLERANCE",
    "GOLDEN_RATIO",
    "SPEED_OF_LIGHT",
    "BandwidthPoint",
    "ConvergenceWarning",
    "DispersionScan",
    "DispersionSolver",
    "FibonacciCell",
    "FrequencyInterval",
    "InvalidParameterError",
    "PWESolver",
    "PlasmonModeReport",
    "Polarization",
    "SuperlatticeSpec",
    "TMMSolver",
    "__version__",
    "allowed_intervals",
    "analytic_r0_r1_r2",
    "average_epsilon_mu",
    "bandwidth_versus_angle",
    "build_cell",
    "cell_length",
    "cell_transfer_matrix",
    "detect_plasmon_modes",
    "get_solver",
    "load_figure",
    "load_spec",
    "n_layers_a",
    "n_layers_b",
    "paper_fibonacci",
    "semitrace_by_product",
    "semitrace_by_recurrence",
    "sequence",
    "zero_average_index_frequency_ghz",
]
