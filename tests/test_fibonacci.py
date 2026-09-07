"""Test 1–3: secuencias, número de capas y longitud de celda."""

import pytest

from fibonacci_tmm.constants import GOLDEN_RATIO
from fibonacci_tmm.fibonacci import (
    build_cell,
    cell_length,
    n_layers_a,
    n_layers_b,
    n_layers_total,
    paper_fibonacci,
    ratio_tau,
    sequence,
    sequence_by_inflation,
    substitution_identity,
)

EXPECTED_SEQUENCES = {
    0: "B",
    1: "A",
    2: "AB",
    3: "ABA",
    4: "ABAAB",
    5: "ABAABABA",
    6: "ABAABABAABAAB",
}

EXPECTED_F = {0: 1, 1: 1, 2: 2, 3: 3, 4: 5, 5: 8, 6: 13, 7: 21, 8: 34}


def test_paper_fibonacci_initial_conditions():
    assert paper_fibonacci(0) == 1
    assert paper_fibonacci(1) == 1
    assert paper_fibonacci(-1) == 0
    assert paper_fibonacci(-2) == 1


@pytest.mark.parametrize("m, value", EXPECTED_F.items())
def test_paper_fibonacci_values(m, value):
    assert paper_fibonacci(m) == value


@pytest.mark.parametrize("m, seq", EXPECTED_SEQUENCES.items())
def test_concatenation_sequences(m, seq):
    assert sequence(m) == seq


@pytest.mark.parametrize("m", range(1, 9))
def test_inflation_matches_concatenation(m):
    assert sequence_by_inflation(m) == sequence(m)


@pytest.mark.parametrize("m", range(2, 9))
def test_layer_counts(m):
    seq = sequence(m)
    assert seq.count("A") == n_layers_a(m) == paper_fibonacci(m - 1)
    assert seq.count("B") == n_layers_b(m) == paper_fibonacci(m - 2)
    assert len(seq) == n_layers_total(m) == paper_fibonacci(m)


@pytest.mark.parametrize("m, expected", [(3, 1), (4, 2), (5, 3), (6, 5), (7, 8)])
def test_number_of_b_layers_is_f_m_minus_2(m, expected):
    assert n_layers_b(m) == expected


def test_cell_length_formula():
    a, b = 0.012, 0.012
    for m in range(2, 9):
        expected = paper_fibonacci(m - 1) * a + paper_fibonacci(m - 2) * b
        assert cell_length(m, a, b) == pytest.approx(expected)
        cell = build_cell(m, a, b)
        assert cell.length == pytest.approx(expected)
        assert cell.n_a == paper_fibonacci(m - 1)
        assert cell.n_b == paper_fibonacci(m - 2)


def test_golden_ratio_convergence():
    assert ratio_tau(12) == pytest.approx(GOLDEN_RATIO, rel=1e-4)


@pytest.mark.parametrize("m", range(2, 8))
@pytest.mark.parametrize("k", range(0, 5))
def test_substitution_identity(m, k):
    if k > m:
        pytest.skip("k > m")
    assert substitution_identity(m, k) == sequence(m)
