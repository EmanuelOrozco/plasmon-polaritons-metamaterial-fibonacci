"""Secuencias de Fibonacci según Reyes-Gómez et al., PRB 81, 153101 (2010).

Indexación del paper, no la clásica F_0 = 0:

    F_0 = F_1 = 1,   F_m = F_{m−1} + F_{m−2}

    S_0 = B,   S_1 = A,   S_m = S_{m−1} S_{m−2}

La celda S_m tiene N_A = F_{m−1} capas A, N_B = F_{m−2} capas B y longitud

    L_m = F_{m−1} a + F_{m−2} b.
"""

from __future__ import annotations

import operator
from dataclasses import dataclass
from functools import cache

from fibonacci_photonics.errors import InvalidParameterError
from fibonacci_photonics.physics.constants import GOLDEN_RATIO

LAYER_A = "A"
LAYER_B = "B"


def validate_order(m: int, *, minimum: int = 0) -> int:
    """Orden de Fibonacci entero ≥ ``minimum``.

    Raises
    ------
    InvalidParameterError
        Si ``m`` no es entero o es menor que ``minimum``.
    """
    try:
        if isinstance(m, bool):
            raise TypeError
        value = operator.index(m)
    except TypeError as exc:
        raise InvalidParameterError(f"el orden de Fibonacci debe ser entero: {m!r}") from exc
    if value < minimum:
        raise InvalidParameterError(f"el orden de Fibonacci debe ser ≥ {minimum}: {value}")
    return value


@dataclass(frozen=True)
class FibonacciCell:
    """Celda elemental S_m.

    Attributes
    ----------
    order : int
        Orden m.
    sequence : str
        Palabra de capas, p. ej. ``"ABAAB"`` para m = 4.
    n_a, n_b, n_total : int
        Número de capas A, B y total (F_{m−1}, F_{m−2}, F_m).
    length : float
        L_m en las unidades de los espesores dados.
    """

    order: int
    sequence: str
    n_a: int
    n_b: int
    n_total: int
    length: float


def paper_fibonacci(m: int) -> int:
    """Número de Fibonacci del paper, F_0 = F_1 = 1.

    Para m = −1 y m = −2 se usa la extensión F_{−1} = 0, F_{−2} = 1, de modo que
    N_A = F_{m−1} y N_B = F_{m−2} valgan también en m = 0, 1.
    """
    m = validate_order(m, minimum=-2)
    if m == -2:
        return 1
    if m == -1:
        return 0
    if m in (0, 1):
        return 1
    prev_prev, prev = 1, 1
    for _ in range(2, m + 1):
        prev_prev, prev = prev, prev + prev_prev
    return prev


@cache
def sequence(m: int) -> str:
    """Palabra S_m por concatenación S_m = S_{m−1} S_{m−2}."""
    m = validate_order(m)
    if m == 0:
        return LAYER_B
    if m == 1:
        return LAYER_A
    return sequence(m - 1) + sequence(m - 2)


@cache
def sequence_by_inflation(m: int) -> str:
    """Palabra S_m por la regla de inflación B → A, A → AB."""
    m = validate_order(m)
    if m == 0:
        return LAYER_B
    if m == 1:
        return LAYER_A
    prev = sequence_by_inflation(m - 1)
    return prev.replace(LAYER_B, "X").replace(LAYER_A, "AB").replace("X", LAYER_A)


def n_layers_a(m: int) -> int:
    """Número de capas A en S_m: F_{m−1}."""
    return paper_fibonacci(validate_order(m) - 1)


def n_layers_b(m: int) -> int:
    """Número de capas B en S_m: F_{m−2}. Es también el número de subbandas plasmon-polaritón."""
    return paper_fibonacci(validate_order(m) - 2)


def n_layers_total(m: int) -> int:
    """Número total de capas de S_m: F_m."""
    return paper_fibonacci(validate_order(m))


def cell_length(m: int, thickness_a: float, thickness_b: float) -> float:
    """L_m = F_{m−1} a + F_{m−2} b, en las unidades de ``thickness_a`` y ``thickness_b``."""
    return n_layers_a(m) * thickness_a + n_layers_b(m) * thickness_b


def build_cell(m: int, thickness_a: float, thickness_b: float) -> FibonacciCell:
    """Celda S_m con sus conteos de capas y su longitud."""
    seq = sequence(m)
    return FibonacciCell(
        order=m,
        sequence=seq,
        n_a=seq.count(LAYER_A),
        n_b=seq.count(LAYER_B),
        n_total=len(seq),
        length=cell_length(m, thickness_a, thickness_b),
    )


def golden_ratio() -> float:
    """τ = (1 + √5)/2."""
    return GOLDEN_RATIO


def ratio_tau(m: int) -> float:
    """τ_m = F_{m+1}/F_m, que converge a τ."""
    return paper_fibonacci(m + 1) / paper_fibonacci(m)


def substitution_identity(m: int, k: int) -> str:
    """S_m(A, B) = S_{m−k}[S_{k+1}(A, B), S_k(A, B)], con 0 ≤ k ≤ m."""
    if k < 0 or k > m:
        raise InvalidParameterError("se requiere 0 <= k <= m")
    if m < 2:
        return sequence(m)
    block_plus = sequence(k + 1)
    block = sequence(k)
    pattern = sequence(m - k)
    return "".join(block_plus if symbol == LAYER_A else block for symbol in pattern)
