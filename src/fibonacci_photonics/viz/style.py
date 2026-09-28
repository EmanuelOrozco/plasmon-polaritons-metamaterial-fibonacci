"""Estilo gráfico del Brief Report de Physical Review B y guardado de figuras."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
from matplotlib import patheffects
from matplotlib.colors import to_rgb
from matplotlib.figure import Figure

# Paleta del original: S3 discontinuo rojo, S4 continuo azul.
# Figs. 5–6, de menor a mayor ν: negro, rojo, verde, azul, magenta, cian.
S3_STYLE: dict[str, Any] = {"color": "#c23b22", "linestyle": "--", "linewidth": 0.65}
S4_STYLE: dict[str, Any] = {"color": "#2c4d8c", "linestyle": "-", "linewidth": 0.7}
ORDER_STYLES: dict[int, dict[str, Any]] = {3: S3_STYLE, 4: S4_STYLE}
BRANCH_COLORS = ["#1a1a1a", "#c23b22", "#2e8b4a", "#2c4d8c", "#c44ec0", "#17becf", "#d35400", "#7f7f7f"]
M_COLORS = {3: "#c23b22", 4: "#2c4d8c", 5: "#2e8b4a", 6: "#c44ec0"}

# Comparación de métodos: TMM en línea continua, PWE en marcadores huecos.
TMM_STYLE: dict[str, Any] = {"color": "#2c4d8c", "linestyle": "-", "linewidth": 0.8}
PWE_STYLE: dict[str, Any] = {
    "color": "#1a1a1a",
    "marker": "o",
    "markersize": 2.2,
    "markerfacecolor": "none",
    "markeredgewidth": 0.5,
    "linestyle": "none",
}

# Superposición de las figuras: la curva TMM de la figura como trazo ancho y
# claro (opaco: las mitades ±k se solapan en k = 0), y encima la curva PWE con
# el estilo de la figura. Las subbandas PWE de las Figs. 5–6 son trazos blancos
# con borde negro, visibles sobre cualquier relleno.
TMM_UNDERLAY: dict[str, Any] = {"linestyle": "-", "linewidth": 3.0}
TMM_UNDERLAY_LIGHTNESS = 0.7
PWE_BAR_STYLE: dict[str, Any] = {
    "color": "white",
    "linewidth": 1.0,
    "path_effects": [patheffects.withStroke(linewidth=2.4, foreground="black")],
}


def lighten(color: str, amount: float = TMM_UNDERLAY_LIGHTNESS) -> tuple[float, float, float]:
    """Mezcla ``color`` con blanco: 0 lo deja igual, 1 da blanco."""
    r, g, b = to_rgb(color)
    return (r + (1 - r) * amount, g + (1 - g) * amount, b + (1 - b) * amount)


TMM_FORMATS = ("pdf", "png", "svg")
PWE_FORMATS = ("pdf", "png")


def latex_theta(label: str) -> str:
    """Convierte ``'pi/12'`` en ``'\\pi/12'`` si π aún no está escapado."""
    text = str(label)
    if r"\pi" in text:
        return text
    return text.replace("pi", r"\pi")


def apply_prb_style() -> None:
    """Tipografía serif de 9 pt, marcas hacia dentro y 300 dpi."""
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 9,
            "axes.labelsize": 9,
            "axes.titlesize": 9,
            "legend.fontsize": 8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "axes.linewidth": 0.7,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "pdf.fonttype": 42,
        }
    )


def save_figure(
    fig: Figure,
    output_dir: Path,
    stem: str,
    formats: tuple[str, ...] = TMM_FORMATS,
) -> dict[str, Path]:
    """Guarda ``fig`` como ``output_dir/stem.<fmt>`` en cada formato y la cierra."""
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {fmt: output_dir / f"{stem}.{fmt}" for fmt in formats}
    for fmt, path in paths.items():
        fig.savefig(path, dpi=300 if fmt == "png" else None)
    plt.close(fig)
    return paths
