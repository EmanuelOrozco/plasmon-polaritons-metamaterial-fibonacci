"""Formato de cifras y tablas LaTeX generadas a partir de los datos."""

from __future__ import annotations

from pathlib import Path

import numpy as np

THETA_TEX = {0.0: "0", round(np.pi / 12, 6): r"\pi/12", round(np.pi / 6, 6): r"\pi/6", round(np.pi / 3, 6): r"\pi/3"}


def theta_tex(theta: float, *, degrees_fallback: bool = False) -> str:
    """θ en rad como ``\\pi/12``, ``\\pi/6``, ...; si no es uno de ellos, en grados o decimal."""
    name = THETA_TEX.get(round(theta, 6))
    if name is not None:
        return name
    return f"{np.degrees(theta):.0f}^\\circ" if degrees_fallback else f"{theta:.3f}"


def sci(value: float, digits: int = 1) -> str:
    """``$a.b\\times10^{n}$``; ``0`` y ``--`` para cero y no finitos."""
    if value == 0 or not np.isfinite(value):
        return "0" if value == 0 else "--"
    exponent = int(np.floor(np.log10(abs(value))))
    return rf"${value / 10**exponent:.{digits}f}\times10^{{{exponent}}}$"


def sci_plain(value: float, digits: int = 0) -> str:
    """``a\\times10^{n}`` sin delimitadores de modo matemático, con mantisa redondeada."""
    exponent = int(np.floor(np.log10(abs(value))))
    mantissa = round(value / 10**exponent, digits)
    if mantissa >= 10:
        mantissa, exponent = mantissa / 10, exponent + 1
    return rf"{mantissa:.{digits}f}\times10^{{{exponent}}}"


def seconds_text(value: float) -> str:
    """Duración con la unidad más legible (µs, ms, s, h)."""
    if value < 1e-3:
        return rf"{value * 1e6:.1f}~$\mu$s"
    if value < 1:
        return f"{value * 1e3:.1f}~ms"
    if value < 3600:
        return f"{value:.1f}~s"
    return f"{value / 3600:.1f}~h"


def human_time(seconds: float) -> str:
    """Duración con cifras redondeadas para el texto (µs, ms, s, min, h)."""
    if seconds < 1e-3:
        return rf"{seconds * 1e6:.2g}~$\mu$s"
    if seconds < 1:
        return f"{seconds * 1e3:.{1 if seconds < 1e-2 else 0}f}~ms"
    if seconds < 60:
        return f"{seconds:.1f}~s"
    if seconds < 3600:
        return f"{seconds / 60:.0f}~min"
    return f"{seconds / 3600:.{0 if seconds >= 36000 else 1}f}~h"


def human_bytes(n_bytes: float) -> str:
    """Tamaño en MB o GB."""
    if n_bytes < 1e9:
        return f"{n_bytes / 1e6:.0f}~MB"
    return f"{n_bytes / 1e9:.{0 if n_bytes >= 1e10 else 1}f}~GB"


def write_table(directory: Path, name: str, header: str, rows: list[str], columns: str) -> Path:
    """Tabla ``booktabs`` en ``directory/name.tex``; se incluye tal cual en los documentos."""
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"{name}.tex"
    body = "\n".join(rows)
    path.write_text(
        f"\\begin{{tabular}}{{{columns}}}\n\\toprule\n{header} \\\\\n\\midrule\n"
        f"{body}\n\\bottomrule\n\\end{{tabular}}\n",
        encoding="utf-8",
    )
    return path
