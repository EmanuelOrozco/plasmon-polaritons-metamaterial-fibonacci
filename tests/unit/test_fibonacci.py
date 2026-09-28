"""Palabras de Fibonacci, número de capas y longitud de la celda."""

import pytest

from fibonacci_photonics.errors import InvalidParameterError
from fibonacci_photonics.physics.fibonacci import (
    build_cell,
    cell_length,
    n_layers_a,
    n_layers_b,
    n_layers_total,
    paper_fibonacci,
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


def test_cell_length_formula():
    a, b = 0.012, 0.007
    for m in range(2, 9):
        expected = paper_fibonacci(m - 1) * a + paper_fibonacci(m - 2) * b
        assert cell_length(m, a, b) == pytest.approx(expected)
        cell = build_cell(m, a, b)
        assert cell.length == pytest.approx(expected)
        assert (cell.n_a, cell.n_b) == (paper_fibonacci(m - 1), paper_fibonacci(m - 2))


@pytest.mark.parametrize("m", range(2, 8))
@pytest.mark.parametrize("k", range(0, 5))
def test_substitution_identity(m, k):
    if k > m:
        pytest.skip("k > m")
    assert substitution_identity(m, k) == sequence(m)


@pytest.mark.parametrize("m", [-1, 2.5, True])
def test_invalid_orders_are_rejected(m):
    with pytest.raises(InvalidParameterError):
        sequence(m)
