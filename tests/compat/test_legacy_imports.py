"""El alias ``fibonacci_tmm`` expone la API antigua (presentación y notebooks)."""

import importlib

import numpy as np
import pytest

LEGACY = {
    "constants": "core.constants",
    "electromagnetics": "core.electromagnetics",
    "materials": "core.materials",
    "fibonacci": "core.fibonacci",
    "params": "core.params",
    "provenance": "core.provenance",
    "transfer_matrix": "tmm.transfer_matrix",
    "plasmon_modes": "analysis.plasmon_modes",
    "bandwidth": "analysis.bandwidth",
    "plotting": "plotting",
}


@pytest.mark.parametrize("old,new", sorted(LEGACY.items()))
def test_legacy_module_is_the_new_module(old, new):
    legacy = importlib.import_module(f"fibonacci_tmm.{old}")
    assert legacy is importlib.import_module(f"fibonacci_photonics.{new}")


def test_legacy_dispersion_api(paper_equal_plasma):
    from fibonacci_tmm.dispersion import (
        allowed_intervals,
        average_epsilon_mu,
        omega_from_nu_ghz,
        scan_dispersion,
        zero_average_index_frequency_ghz,
    )

    nu = np.linspace(0.5, 5.0, 200)
    scan = scan_dispersion(paper_equal_plasma, 3, np.pi / 6, nu)
    assert scan.method == "recurrence"
    assert allowed_intervals(scan)
    eps, mu = average_epsilon_mu(paper_equal_plasma, 3, omega_from_nu_ghz(2.0))
    assert np.isfinite(eps) and np.isfinite(mu)
    assert zero_average_index_frequency_ghz(paper_equal_plasma, 3) is not None


def test_legacy_package_exports_tmm_scan():
    import fibonacci_tmm
    from fibonacci_photonics.tmm.dispersion import scan_dispersion

    assert fibonacci_tmm.scan_dispersion is scan_dispersion
