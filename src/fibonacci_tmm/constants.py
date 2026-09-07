"""Constantes físicas y conversiones de unidades.

Las figuras del paper usan ν en GHz y espesores en mm. El núcleo trabaja en SI.
El valor de c no aparece en el artículo; se documenta como decisión de implementación.
"""

from __future__ import annotations

import math

# [IMPL] SI exacto. El paper no especifica c.
SPEED_OF_LIGHT = 299_792_458.0  # m/s

# [IMPL] Alternativa usada a veces en microondas; para análisis de convergencia.
SPEED_OF_LIGHT_ROUNDED = 3.0e8  # m/s

GHZ = 1.0e9  # Hz
MM = 1.0e-3  # m

# [PAPER] τ = (1 + √5) / 2
GOLDEN_RATIO = (1.0 + math.sqrt(5.0)) / 2.0

# [IMPL] Tolerancia para |R| ≤ 1 (redondeo, no clipping silencioso).
BAND_ABS_R_TOLERANCE = 1.0e-12

# [IMPL] Umbral para la forma límite de sin(Qd)/Q cuando Q → 0.
SMALL_Q_ABS = 1.0e-14
