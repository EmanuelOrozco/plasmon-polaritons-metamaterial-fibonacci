"""Relación de dispersión cos(k Lm) = Rm con la semitraza de la matriz de transferencia."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike

from fibonacci_photonics.core.bands import DispersionScan, omega_from_nu_ghz, scan_from_semitrace
from fibonacci_photonics.core.constants import BAND_ABS_R_TOLERANCE
from fibonacci_photonics.core.electromagnetics import Polarization
from fibonacci_photonics.core.params import SuperlatticeSpec
from fibonacci_photonics.tmm.transfer_matrix import semitrace_by_product, semitrace_by_recurrence


def scan_dispersion(
    spec: SuperlatticeSpec,
    m: int,
    theta: float,
    nu_ghz: ArrayLike,
    polarization: Polarization = Polarization.TE,
    *,
    method: str = "recurrence",
    atol: float = BAND_ABS_R_TOLERANCE,
) -> DispersionScan:
    nu = np.asarray(nu_ghz, dtype=np.float64)
    omega = omega_from_nu_ghz(nu)
    if method == "recurrence":
        r = semitrace_by_recurrence(spec, m, omega, theta, polarization)
    elif method == "product":
        r = semitrace_by_product(spec, m, omega, theta, polarization)
    else:
        raise ValueError("method debe ser 'recurrence' o 'product'")
    return scan_from_semitrace(spec, m, theta, nu, r, polarization, method, atol=atol)
