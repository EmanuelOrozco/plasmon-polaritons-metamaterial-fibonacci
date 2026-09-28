"""Constantes físicas y factores de unidad.

Las figuras del paper usan ν en GHz y espesores en mm; el núcleo trabaja en SI.
El valor de c no aparece en el artículo y se documenta como decisión de
implementación.
"""

from __future__ import annotations

import math

SPEED_OF_LIGHT = 299_792_458.0
"""c en m/s (SI exacto). [IMPL] El paper no especifica c."""

SPEED_OF_LIGHT_ROUNDED = 3.0e8
"""c redondeada, en m/s; alternativa habitual en microondas (estudios de sensibilidad)."""

GHZ = 1.0e9
"""1 GHz en Hz."""

MM = 1.0e-3
"""1 mm en m."""

GOLDEN_RATIO = (1.0 + math.sqrt(5.0)) / 2.0
"""τ = (1 + √5)/2, límite de F_{m+1}/F_m. [PAPER]"""
