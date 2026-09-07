"""Secuencias de Fibonacci según Reyes-Gómez et al., PRB 81, 153101 (2010).

Indexación del paper, no la clásica F0 = 0:

    F0 = F1 = 1,  Fm = F(m-1) + F(m-2)

    S0 = B,  S1 = A,  Sm = S(m-1) S(m-2)

    Lm = F(m-1) a + F(m-2) b
    NB = F(m-2)
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from math import sqrt

from fibonacci_tmm.constants import GOLDEN_RATIO

LAYER_A = "A"
LAYER_B = "B"


@dataclass(frozen=True)
class FibonacciCell:
    """Celda elemental Sm: secuencia, conteos y longitud."""

    order: int
    sequence: str
    n_a: int
    n_b: int
    n_total: int
    length: float


def paper_fibonacci(m: int) -> int:
    """Número de Fibonacci del paper: F0 = F1 = 1.

    Para m = -1 y m = -2 se usa la extensión F(-1) = 0, F(-2) = 1,
    de modo que NA = F(m-1) y NB = F(m-2) valgan también en m = 0, 1.
    """
    if m < -2:
        raise ValueError(f"orden de Fibonacci no soportado: {m}")
    if m == -2:
        return 1
    if m == -1:
        return 0
    if m == 0 or m == 1:
        return 1
    prev_prev, prev = 1, 1
    for _ in range(2, m + 1):
        prev_prev, prev = prev, prev + prev_prev
    return prev


@lru_cache(maxsize=None)
def sequence(m: int) -> str:
    """Sm por concatenación Sm = S(m-1) S(m-2)."""
    if m < 0:
        raise ValueError("el orden m de la secuencia debe ser >= 0")
    if m == 0:
        return LAYER_B
    if m == 1:
        return LAYER_A
    return sequence(m - 1) + sequence(m - 2)


@lru_cache(maxsize=None)
def sequence_by_inflation(m: int) -> str:
    """Sm por la ley de inflación B → A, A → AB."""
    if m < 0:
        raise ValueError("el orden m de la secuencia debe ser >= 0")
    if m == 0:
        return LAYER_B
    if m == 1:
        return LAYER_A
    prev = sequence_by_inflation(m - 1)
    return prev.replace(LAYER_B, "X").replace(LAYER_A, "AB").replace("X", LAYER_A)


def n_layers_a(m: int) -> int:
    return paper_fibonacci(m - 1)


def n_layers_b(m: int) -> int:
    """Número de capas metamaterial B en Sm. Igual a F(m-2) y al número de modos."""
    return paper_fibonacci(m - 2)


def n_layers_total(m: int) -> int:
    return paper_fibonacci(m)


def cell_length(m: int, thickness_a: float, thickness_b: float) -> float:
    """Lm = F(m-1) a + F(m-2) b."""
    return n_layers_a(m) * thickness_a + n_layers_b(m) * thickness_b


def build_cell(m: int, thickness_a: float, thickness_b: float) -> FibonacciCell:
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
    return GOLDEN_RATIO


def ratio_tau(m: int) -> float:
    """τ_m = F(m+1) / F(m). Converge a (1+√5)/2."""
    return paper_fibonacci(m + 1) / paper_fibonacci(m)


def substitution_identity(m: int, k: int) -> str:
    """Sm(A,B) = S_{m-k}[S_{k+1}(A,B), S_k(A,B)]."""
    if k < 0 or k > m:
        raise ValueError("se requiere 0 <= k <= m")
    if m < 2:
        return sequence(m)
    block_plus = sequence(k + 1)
    block = sequence(k)
    pattern = sequence(m - k)
    return "".join(block_plus if symbol == LAYER_A else block for symbol in pattern)
