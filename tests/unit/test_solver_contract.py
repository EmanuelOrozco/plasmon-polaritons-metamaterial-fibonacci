"""Contrato ``DispersionSolver`` y validación de entradas de ambos métodos."""

import warnings

import numpy as np
import pytest

from fibonacci_photonics.errors import ConvergenceWarning, InvalidParameterError
from fibonacci_photonics.solvers import (
    DispersionSolver,
    PWESolver,
    TMMSolver,
    get_solver,
    validate_frequencies,
)

SOLVERS = {
    "tmm": lambda spec: TMMSolver(spec),
    "pwe": lambda spec: PWESolver(spec, harmonics_per_layer=4),
}


@pytest.fixture(params=sorted(SOLVERS))
def solver(request, paper_split_plasma):
    return SOLVERS[request.param](paper_split_plasma)


def test_both_methods_satisfy_the_protocol(solver):
    assert isinstance(solver, DispersionSolver)
    assert solver.method in ("tmm", "pwe")


def test_get_solver_by_name(paper_split_plasma):
    assert isinstance(get_solver("tmm", paper_split_plasma), TMMSolver)
    pwe = get_solver("pwe", paper_split_plasma, harmonics_per_layer=3)
    assert isinstance(pwe, PWESolver) and pwe.n_max(3) == 9
    with pytest.raises(InvalidParameterError):
        get_solver("fdtd", paper_split_plasma)


def test_semitrace_keeps_the_input_shape(solver):
    nu = np.linspace(0.5, 4.5, 6).reshape(2, 3)
    r = solver.semitrace(3, np.pi / 6, nu, "TE")
    assert r.shape == (2, 3) and r.dtype == np.complex128


def test_scan_contract(solver):
    nu = np.linspace(0.5, 4.5, 40)
    scan = solver.scan(3, np.pi / 3, nu, "TM")
    assert scan.nu_ghz.shape == scan.r.shape == scan.allowed.shape == nu.shape
    assert scan.polarization == "TM" and scan.m == 3
    k = scan.k_lm_over_pi
    assert np.all((k[scan.allowed] >= 0.0) & (k[scan.allowed] <= 1.0))
    assert np.all(np.isnan(k[~scan.allowed]))


def test_evaluator_matches_semitrace_and_counts_calls(solver):
    evaluator = solver.evaluator(3, np.pi / 6, "TE")
    values = [evaluator(nu) for nu in (0.7, 2.2)]
    np.testing.assert_allclose(values, solver.semitrace(3, np.pi / 6, [0.7, 2.2], "TE").real, atol=1e-12)
    assert evaluator.calls == 2


@pytest.mark.parametrize(
    "nu",
    [
        pytest.param([0.0, 1.0], id="cero"),
        pytest.param([-1.0, 1.0], id="negativa"),
        pytest.param([1.0, np.nan], id="nan"),
        pytest.param([1.0, np.inf], id="inf"),
        pytest.param([], id="vacia"),
    ],
)
def test_non_positive_or_non_finite_frequencies_are_rejected(solver, nu):
    with pytest.raises(InvalidParameterError):
        solver.semitrace(3, 0.1, nu)


@pytest.mark.parametrize("nu", [[1.0, 1.0, 2.0], [2.0, 1.0], [[1.0, 2.0]]], ids=["repetida", "decreciente", "2d"])
def test_scan_requires_increasing_1d_grid(solver, nu):
    with pytest.raises(InvalidParameterError):
        solver.scan(3, 0.1, nu)


def test_validate_frequencies_accepts_any_shape_without_grid():
    assert validate_frequencies([[2.0, 1.0]]).shape == (1, 2)


@pytest.mark.parametrize("theta", [-0.1, 1.6, np.nan])
def test_invalid_angles_are_rejected(solver, theta):
    with pytest.raises(InvalidParameterError):
        solver.semitrace(3, theta, [1.0])


@pytest.mark.parametrize("m", [-1, 1.5])
def test_invalid_orders_are_rejected(solver, m):
    with pytest.raises(InvalidParameterError):
        solver.semitrace(m, 0.1, [1.0])


def test_invalid_polarization_is_rejected(solver):
    with pytest.raises((InvalidParameterError, ValueError)):
        solver.semitrace(3, 0.1, [1.0], "XY")


@pytest.mark.parametrize(
    "kwargs",
    [
        {"harmonics_per_layer": 0},
        {"harmonics_per_layer": 2.5},
        {"harmonics_per_layer": {3: 0}},
        {"rule": "fft"},
        {"workers": 0},
        {"n_max_scale": 0.0},
    ],
)
def test_pwe_rejects_invalid_numerics(paper_split_plasma, kwargs):
    with pytest.raises(InvalidParameterError):
        PWESolver(paper_split_plasma, **kwargs)


def test_tmm_rejects_unknown_algorithm(paper_split_plasma):
    with pytest.raises(InvalidParameterError):
        TMMSolver(paper_split_plasma, algorithm="fast")


def test_pwe_harmonics_per_order(paper_split_plasma):
    solver = PWESolver(paper_split_plasma, harmonics_per_layer={3: 4, 5: 2})
    assert (solver.n_max(3), solver.n_max(5)) == (12, 16)
    assert solver.n_plane_waves(3) == 25
    assert solver.refined(1.5).n_max(3) == 18
    with pytest.raises(InvalidParameterError):
        solver.harmonics(4)


def test_convergence_report_passes_with_enough_harmonics(paper_split_plasma):
    solver = PWESolver(paper_split_plasma, harmonics_per_layer=16)
    with warnings.catch_warnings():
        warnings.simplefilter("error", ConvergenceWarning)
        report = solver.convergence_report(3, np.pi / 6, np.linspace(1.5, 4.5, 8))
    assert report.converged
    assert report.n_max_check == 72
    assert report.as_dict()["converged"] is True


def test_convergence_warning_near_the_magnetic_plasmon(paper_split_plasma):
    """Pegado a ν_m, |Im k L| es enorme y un truncamiento corto no converge."""
    solver = PWESolver(paper_split_plasma, harmonics_per_layer=1)
    nu = 1.0 - np.logspace(-6, -2, 12)[::-1]
    with pytest.warns(ConvergenceWarning):
        report = solver.convergence_report(6, np.pi / 3, nu, factor=3.0, tolerance=1e-6)
    assert not report.converged


def test_convergence_factor_must_exceed_one(paper_split_plasma):
    with pytest.raises(InvalidParameterError):
        PWESolver(paper_split_plasma).convergence_report(3, 0.1, [1.5], factor=1.0)
