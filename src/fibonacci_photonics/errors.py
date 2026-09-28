"""Excepciones y advertencias propias del paquete."""

from __future__ import annotations


class InvalidParameterError(ValueError):
    """Parámetro físico o numérico fuera de su dominio (ν ≤ 0, θ ∉ [0, π/2], m < 0, ...)."""


class ConvergenceWarning(UserWarning):
    """El resultado no cambió menos que la tolerancia al refinar la discretización."""
