"""Visualización: estilo PRB, curvas de dispersión y guardado de figuras."""

from fibonacci_photonics.viz.dispersion import format_dispersion_axes, plot_dispersion_branches
from fibonacci_photonics.viz.style import (
    BRANCH_COLORS,
    M_COLORS,
    ORDER_STYLES,
    PWE_FORMATS,
    PWE_STYLE,
    S3_STYLE,
    S4_STYLE,
    TMM_FORMATS,
    TMM_STYLE,
    apply_prb_style,
    latex_theta,
    save_figure,
)

__all__ = [
    "BRANCH_COLORS",
    "M_COLORS",
    "ORDER_STYLES",
    "PWE_FORMATS",
    "PWE_STYLE",
    "S3_STYLE",
    "S4_STYLE",
    "TMM_FORMATS",
    "TMM_STYLE",
    "apply_prb_style",
    "format_dispersion_axes",
    "latex_theta",
    "plot_dispersion_branches",
    "save_figure",
]
