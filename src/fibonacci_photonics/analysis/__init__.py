"""Post-proceso de barridos R(ν), independiente del método que los produjo.

- ``roots``: Brent acotado.
- ``band_edges``: bandas y gaps más estrechos que la malla, con bordes refinados.
- ``band_closure``: bordes |R| = 1 y puntas para cerrar las curvas ν(k).
- ``plasmon_modes``: conteo de subbandas F_{m−2} en una ventana.
- ``plasmon_bands``: estrategias densa (TMM) y adaptativa validada (PWE).
- ``bandwidth``: anchos de banda frente al ángulo.
"""
