"""Los métodos reproducen los valores guardados en ``reference_values.json``.

El JSON se generó con la versión anterior a la refactorización; cualquier cambio
numérico (no solo de estructura) en la TMM, el PWE, el cierre de bandas o la
búsqueda de subbandas rompe estas pruebas.
"""

import json
from pathlib import Path

import numpy as np
import pytest

from fibonacci_photonics.analysis.band_closure import close_band_edges
from fibonacci_photonics.analysis.band_edges import locate_bands, plasmon_grid
from fibonacci_photonics.config import spec_from_mapping
from fibonacci_photonics.solvers import PWESolver, TMMSolver

REFERENCE = json.loads((Path(__file__).with_name("reference_values.json")).read_text(encoding="utf-8"))
SPECS = {name: spec_from_mapping(physics) for name, physics in REFERENCE["physics"].items()}
NU = np.array(REFERENCE["nu_ghz"])


def _case_id(case):
    return f"{case['spec']}-m{case['m']}-{case['polarization']}-th{case['theta']:.3f}"


@pytest.mark.parametrize("case", REFERENCE["tmm"], ids=_case_id)
def test_tmm_semitrace(case):
    r = TMMSolver(SPECS[case["spec"]]).semitrace(case["m"], case["theta"], NU, case["polarization"])
    reference = np.array(case["r_real"]) + 1j * np.array(case["r_imag"])
    np.testing.assert_allclose(r, reference, rtol=1e-12, atol=1e-12)


@pytest.mark.parametrize("case", REFERENCE["pwe"], ids=_case_id)
def test_pwe_semitrace(case):
    solver = PWESolver(SPECS[case["spec"]], harmonics_per_layer=case["harmonics_per_layer"])
    r = solver.semitrace(case["m"], case["theta"], NU, case["polarization"]).real
    reference = np.array(case["r_real"], dtype=float)
    np.testing.assert_array_equal(np.isnan(r), np.isnan(reference))
    np.testing.assert_allclose(r, reference, rtol=1e-9, atol=1e-9, equal_nan=True)


def test_pwe_band_edges():
    case = REFERENCE["pwe_edges"]
    solver = PWESolver(SPECS[case["spec"]], harmonics_per_layer=case["harmonics_per_layer"])
    scan = solver.scan(case["m"], case["theta"], np.linspace(*case["grid"]), "TE")
    edges = close_band_edges(solver, scan)
    np.testing.assert_allclose(edges.nu_ghz, case["nu_ghz"], rtol=0, atol=1e-10)
    np.testing.assert_array_equal(edges.kind, case["kind"])
    np.testing.assert_allclose(edges.k_lm_over_pi, case["k_lm_over_pi"], rtol=0, atol=1e-8)


def test_tmm_plasmon_bands():
    case = REFERENCE["tmm_plasmon_bands"]
    grid = plasmon_grid(*case["grid"])
    evaluator = TMMSolver(SPECS[case["spec"]]).evaluator(case["m"], case["theta"], "TE")
    bands = locate_bands(evaluator, grid, np.array([evaluator(float(nu)) for nu in grid]))
    found = [[band.nu_min_ghz, band.nu_max_ghz] for band in bands]
    np.testing.assert_allclose(found, case["bands"], rtol=0, atol=1e-12)
