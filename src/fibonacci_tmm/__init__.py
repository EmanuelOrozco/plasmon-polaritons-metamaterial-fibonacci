"""Reimplementación TMM de Reyes-Gómez et al., Phys. Rev. B 81, 153101 (2010)."""

from fibonacci_tmm.bandwidth import BandwidthPoint, bandwidth_versus_angle
from fibonacci_tmm.constants import BAND_ABS_R_TOLERANCE, GOLDEN_RATIO, SPEED_OF_LIGHT
from fibonacci_tmm.dispersion import (
    DispersionScan,
    FrequencyInterval,
    allowed_intervals,
    scan_dispersion,
    zero_average_index_frequency_ghz,
)
from fibonacci_tmm.electromagnetics import Polarization
from fibonacci_tmm.fibonacci import (
    FibonacciCell,
    build_cell,
    cell_length,
    n_layers_a,
    n_layers_b,
    paper_fibonacci,
    sequence,
)
from fibonacci_tmm.params import SuperlatticeSpec, load_spec
from fibonacci_tmm.plasmon_modes import PlasmonModeReport, detect_plasmon_modes
from fibonacci_tmm.transfer_matrix import (
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
    "bandwidth_versus_angle",
    "build_cell",
    "cell_length",
    "cell_transfer_matrix",
    "detect_plasmon_modes",
    "load_spec",
    "n_layers_a",
    "n_layers_b",
    "paper_fibonacci",
    "scan_dispersion",
    "semitrace_by_product",
    "semitrace_by_recurrence",
    "sequence",
    "zero_average_index_frequency_ghz",
]
