"""Alias de compatibilidad: API del antiguo ``fibonacci_tmm.dispersion``.

Reparte entre ``core.bands``, ``core.effective_medium`` y ``tmm.dispersion``.
"""

from fibonacci_photonics.core.bands import (  # noqa: F401
    DispersionScan,
    FrequencyInterval,
    allowed_intervals,
    allowed_mask,
    near_drude_pole,
    nu_ghz_from_omega,
    omega_from_nu_ghz,
    real_semitrace,
    reduced_wavevector,
)
from fibonacci_photonics.core.effective_medium import (  # noqa: F401
    average_epsilon_mu,
    zero_average_index_frequency_ghz,
)
from fibonacci_photonics.tmm.dispersion import scan_dispersion  # noqa: F401
