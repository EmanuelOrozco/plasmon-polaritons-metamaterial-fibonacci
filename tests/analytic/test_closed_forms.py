"""Casos con solución cerrada que ambos métodos deben reproducir.

- Límite homogéneo: si B = A, R_m = cos(Q_A L_m) con Q_A = (ω/c) cos θ.
- Bicapa de cuarto de onda (n_A a = n_B b, φ = π/2): R = −(η + 1/η)/2 = −5/3 para η = 3.
- Unimodularidad det T = 1, invariante de Kohmoto, TE = TM en θ = 0 y dualidad ε ↔ μ.
- Parseval para los indicadores de la celda y razón áurea de las palabras de Fibonacci.
"""

from itertools import pairwise

import numpy as np
import pytest

from fibonacci_photonics.config import spec_from_mapping
from fibonacci_photonics.physics.constants import GOLDEN_RATIO, SPEED_OF_LIGHT
from fibonacci_photonics.physics.effective_medium import zero_average_index_frequency_ghz
from fibonacci_photonics.physics.electromagnetics import Polarization
from fibonacci_photonics.physics.fibonacci import cell_length, paper_fibonacci, ratio_tau
from fibonacci_photonics.physics.units import omega_from_nu_ghz
from fibonacci_photonics.solvers import PWESolver, TMMSolver
from fibonacci_photonics.solvers.pwe.fourier import cell_fourier
from fibonacci_photonics.solvers.tmm.transfer_matrix import (
    cell_transfer_matrix,
    determinant,
    semitrace_by_product,
)

POLARIZATIONS = [Polarization.TE, Polarization.TM]


def _dielectric_superlattice(a_mm: float, b_mm: float, epsilon_b: float):
    """A = aire y B dieléctrico no dispersivo (Drude con ω_e = ω_m = 0)."""
    return spec_from_mapping(
        {
            "layer_a_thickness_mm": a_mm,
            "layer_b_thickness_mm": b_mm,
            "omega_e_over_2pi_ghz": 0.0,
            "omega_m_over_2pi_ghz": 0.0,
            "epsilon_0": epsilon_b,
        }
    )


@pytest.fixture
def homogeneous():
    """B idéntico al aire: la superred es vacío con período L_m."""
    return _dielectric_superlattice(12.0, 7.0, 1.0)


@pytest.mark.parametrize("polarization", POLARIZATIONS)
@pytest.mark.parametrize("m", [2, 3, 5, 7])
@pytest.mark.parametrize("theta", [0.0, np.pi / 5])
def test_homogeneous_limit_tmm(homogeneous, polarization, m, theta):
    nu = np.linspace(0.3, 20.0, 57)
    length = cell_length(m, homogeneous.thickness_a, homogeneous.thickness_b)
    exact = np.cos(omega_from_nu_ghz(nu) / SPEED_OF_LIGHT * np.cos(theta) * length)
    for algorithm in ("recurrence", "product"):
        r = TMMSolver(homogeneous, algorithm=algorithm).semitrace(m, theta, nu, polarization)
        np.testing.assert_allclose(r.real, exact, rtol=0, atol=1e-10)
        np.testing.assert_allclose(r.imag, 0.0, atol=1e-10)


@pytest.mark.parametrize("polarization", POLARIZATIONS)
def test_homogeneous_limit_pwe(homogeneous, polarization):
    """Sin contraste el operador de ondas planas es diagonal: exacto con cualquier truncamiento."""
    m, theta = 4, np.pi / 7
    nu = np.linspace(0.3, 10.0, 17)
    length = cell_length(m, homogeneous.thickness_a, homogeneous.thickness_b)
    exact = np.cos(omega_from_nu_ghz(nu) / SPEED_OF_LIGHT * np.cos(theta) * length)
    r = PWESolver(homogeneous, harmonics_per_layer=2).semitrace(m, theta, nu, polarization)
    np.testing.assert_allclose(r.real, exact, rtol=0, atol=1e-9)


def test_homogeneous_limit_has_no_gaps(homogeneous):
    scan = TMMSolver(homogeneous).scan(5, np.pi / 4, np.linspace(0.1, 30.0, 3000))
    assert scan.allowed.all()


@pytest.fixture
def quarter_wave():
    """Bicapa AB con n_A a = n_B b: a = 3 mm de aire y b = 1 mm con ε = 9 (η = 3)."""
    return _dielectric_superlattice(3.0, 1.0, 9.0)


def _quarter_wave_nu_ghz(spec) -> float:
    """ν con φ_A = ω a / c = π/2."""
    return spec.speed_of_light / (4.0 * spec.thickness_a) / 1e9


@pytest.mark.parametrize("polarization", POLARIZATIONS)
def test_quarter_wave_bilayer_midgap(quarter_wave, polarization):
    nu = np.array([_quarter_wave_nu_ghz(quarter_wave)])
    r = TMMSolver(quarter_wave).semitrace(2, 0.0, nu, polarization)
    assert r.real[0] == pytest.approx(-5.0 / 3.0, abs=1e-12)
    scan = TMMSolver(quarter_wave).scan(2, 0.0, nu, polarization)
    assert not scan.allowed[0]


def test_quarter_wave_bilayer_half_wave_is_band_edge(quarter_wave):
    """A doble frecuencia φ_A = φ_B = π y R = 1: borde de banda en k = 0."""
    nu = np.array([2.0 * _quarter_wave_nu_ghz(quarter_wave)])
    assert TMMSolver(quarter_wave).semitrace(2, 0.0, nu).real[0] == pytest.approx(1.0, abs=1e-12)


def test_quarter_wave_bilayer_pwe(quarter_wave):
    """El PWE resuelve k(ω) complejo en el gap y recupera R = −cosh(Im k L) = −5/3."""
    nu = np.array([_quarter_wave_nu_ghz(quarter_wave)])
    r = PWESolver(quarter_wave, harmonics_per_layer=16).semitrace(2, 0.0, nu)
    assert r.real[0] == pytest.approx(-5.0 / 3.0, abs=1e-4)


@pytest.mark.parametrize("polarization", POLARIZATIONS)
@pytest.mark.parametrize("m", [0, 1, 2, 3, 4, 5, 6])
def test_transfer_matrix_is_unimodular(paper_equal_plasma, polarization, m):
    omega = omega_from_nu_ghz(np.array([0.5, 1.2, 2.1, 2.7, 3.4, 4.5]))
    t = cell_transfer_matrix(paper_equal_plasma, m, omega, np.pi / 6, polarization)
    np.testing.assert_allclose(determinant(t), 1.0, atol=1e-8)


def _kohmoto(r_prev, r_mid, r_next):
    return r_next**2 + r_mid**2 + r_prev**2 - 2.0 * r_next * r_mid * r_prev - 1.0


@pytest.mark.parametrize("polarization", POLARIZATIONS)
def test_kohmoto_invariant_is_constant(paper_split_plasma, polarization):
    """I = R²_{m+1} + R²_m + R²_{m−1} − 2 R_{m+1} R_m R_{m−1} − 1 no depende de m.

    Las semitrazas se calculan con el producto de matrices, no con la recurrencia,
    para que la prueba sea independiente de la identidad que se verifica.
    """
    omega = omega_from_nu_ghz(np.array([0.45, 1.3, 2.2, 3.7, 4.6]))
    theta = np.pi / 5
    r = [semitrace_by_product(paper_split_plasma, m, omega, theta, polarization) for m in range(9)]
    invariants = np.array([_kohmoto(r[m - 1], r[m], r[m + 1]) for m in range(1, 8)])
    scale = np.maximum(1.0, np.max(np.abs(r), axis=0) ** 2)
    np.testing.assert_allclose(invariants - invariants[0], 0.0, atol=1e-9 * scale.max())


def test_kohmoto_invariant_of_dielectric_bilayer(quarter_wave):
    """A θ = 0: I = ¼ (η − 1/η)² sin² φ_A sin² φ_B."""
    nu = np.linspace(1.0, 60.0, 23)
    omega = omega_from_nu_ghz(nu)
    r = [semitrace_by_product(quarter_wave, m, omega, 0.0, Polarization.TE) for m in range(3)]
    phi_a = omega / quarter_wave.speed_of_light * quarter_wave.thickness_a
    phi_b = 3.0 * omega / quarter_wave.speed_of_light * quarter_wave.thickness_b
    expected = 0.25 * (3.0 - 1.0 / 3.0) ** 2 * np.sin(phi_a) ** 2 * np.sin(phi_b) ** 2
    np.testing.assert_allclose(_kohmoto(r[0], r[1], r[2]).real, expected, atol=1e-12)


@pytest.mark.parametrize("m", [3, 4, 6])
def test_normal_incidence_te_equals_tm(paper_split_plasma, m):
    nu = np.linspace(0.3, 5.0, 31)
    solver = TMMSolver(paper_split_plasma)
    np.testing.assert_allclose(solver.semitrace(m, 0.0, nu, "TE"), solver.semitrace(m, 0.0, nu, "TM"), atol=1e-10)


def test_duality_swaps_te_and_tm():
    """Intercambiar ε ↔ μ en ambos medios convierte TE en TM."""
    base = {"layer_a_thickness_mm": 12.0, "layer_b_thickness_mm": 12.0}
    spec = spec_from_mapping({**base, "epsilon_a": 2.0, "mu_a": 1.0, "epsilon_0": 1.0, "mu_0": 1.5,
                              "omega_e_over_2pi_ghz": 3.0, "omega_m_over_2pi_ghz": 1.0})
    dual = spec_from_mapping({**base, "epsilon_a": 1.0, "mu_a": 2.0, "epsilon_0": 1.5, "mu_0": 1.0,
                              "omega_e_over_2pi_ghz": 1.0, "omega_m_over_2pi_ghz": 3.0})
    nu = np.linspace(0.3, 5.0, 37)
    r_te = TMMSolver(spec).semitrace(5, np.pi / 5, nu, "TE")
    r_tm = TMMSolver(dual).semitrace(5, np.pi / 5, nu, "TM")
    np.testing.assert_allclose(r_te, r_tm, rtol=1e-10, atol=1e-10)


@pytest.mark.parametrize("m", [3, 5])
def test_parseval_for_cell_indicators(m):
    """Σ|c_n(1_A)|² = f_A y Σ c_n(1_A) c̄_n(1_B) = 0; la cola decae como 1/N."""
    cell = cell_fourier(m, 1.0, 1.3, 3000)
    fill_a = paper_fibonacci(m - 1) * 1.0 / cell_length(m, 1.0, 1.3)
    assert np.sum(np.abs(cell.indicator_a) ** 2) == pytest.approx(fill_a, abs=2e-4)
    assert abs(np.vdot(cell.indicator_b, cell.indicator_a)) < 2e-4


def test_fibonacci_ratio_tends_to_golden_ratio():
    errors = [abs(ratio_tau(m) - GOLDEN_RATIO) for m in range(3, 16)]
    assert errors[-1] < 1e-6
    assert all(later < earlier for earlier, later in pairwise(errors))
    np.testing.assert_allclose(errors[-1] / errors[-2], GOLDEN_RATIO**-2, rtol=1e-3)


def test_zero_average_index_frequency_formula(paper_equal_plasma):
    """⟨ε⟩ = ⟨μ⟩ = 0 en ν = ν_p √(N_B b / L_m)."""
    assert zero_average_index_frequency_ghz(paper_equal_plasma, 3) == pytest.approx(3.0 * np.sqrt(1 / 3), rel=1e-12)
    assert zero_average_index_frequency_ghz(paper_equal_plasma, 4) == pytest.approx(3.0 * np.sqrt(2 / 5), rel=1e-12)


def test_n0_gap_is_closed_at_normal_incidence(paper_equal_plasma):
    """Con ν_e = ν_m el gap ⟨n⟩ = 0 se cierra en θ = 0 y se abre al inclinar."""
    nu_gap = zero_average_index_frequency_ghz(paper_equal_plasma, 3)
    nu = np.linspace(nu_gap - 0.05, nu_gap + 0.05, 400)
    solver = TMMSolver(paper_equal_plasma)
    normal, oblique = solver.scan(3, 0.0, nu), solver.scan(3, np.pi / 3, nu)
    assert normal.allowed.all()
    assert oblique.allowed.mean() < normal.allowed.mean()
