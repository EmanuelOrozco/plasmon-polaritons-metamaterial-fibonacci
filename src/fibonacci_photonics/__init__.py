"""Plasmon-polaritones en superredes de Fibonacci (Reyes-Gómez et al., PRB 81, 153101, 2010).

Dos métodos numéricos independientes calculan la misma semitraza R(ν) = cos(k Lm):

- ``fibonacci_photonics.tmm``: matriz de transferencia (analítica por capas).
- ``fibonacci_photonics.pwe``: expansión en ondas planas (k(ω), regla inversa).

``core`` contiene la física compartida y ``analysis`` el post-proceso, que no
depende del método que produjo el barrido.
"""

from fibonacci_photonics.analysis.bandwidth import BandwidthPoint, bandwidth_versus_angle
from fibonacci_photonics.analysis.plasmon_modes import PlasmonModeReport, detect_plasmon_modes
from fibonacci_photonics.core.bands import DispersionScan, FrequencyInterval, allowed_intervals
from fibonacci_photonics.core.constants import BAND_ABS_R_TOLERANCE, GOLDEN_RATIO, SPEED_OF_LIGHT
from fibonacci_photonics.core.effective_medium import (
    average_epsilon_mu,
    zero_average_index_frequency_ghz,
)
from fibonacci_photonics.core.electromagnetics import Polarization
from fibonacci_photonics.core.fibonacci import (
    FibonacciCell,
    build_cell,
    cell_length,
    n_layers_a,
    n_layers_b,
    paper_fibonacci,
    sequence,
)
from fibonacci_photonics.core.params import SuperlatticeSpec, load_spec
from fibonacci_photonics.pwe.dispersion import scan_dispersion as scan_dispersion_pwe
from fibonacci_photonics.tmm.dispersion import scan_dispersion as scan_dispersion_tmm
from fibonacci_photonics.tmm.transfer_matrix import (
    analytic_r0_r1_r2,
    cell_transfer_matrix,
    semitrace_by_product,
    semitrace_by_recurrence,
)

__all__ = [
    "BAND_ABS_R_TOLERANCE",
    "BandwidthPoint",
    "DispersionScan",
    "FibonacciCell",
    "FrequencyInterval",
    "GOLDEN_RATIO",
    "PlasmonModeReport",
    "Polarization",
    "SPEED_OF_LIGHT",
    "SuperlatticeSpec",
    "allowed_intervals",
    "analytic_r0_r1_r2",
    "average_epsilon_mu",
    "bandwidth_versus_angle",
    "build_cell",
    "cell_length",
    "cell_transfer_matrix",
    "detect_plasmon_modes",
    "load_spec",
    "n_layers_a",
    "n_layers_b",
    "paper_fibonacci",
    "scan_dispersion_pwe",
    "scan_dispersion_tmm",
    "semitrace_by_product",
    "semitrace_by_recurrence",
    "sequence",
    "zero_average_index_frequency_ghz",
]
