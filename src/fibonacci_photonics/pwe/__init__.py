"""Método de expansión en ondas planas (PWE) para la superred de Fibonacci.

- ``fourier``: coeficientes analíticos de perfiles escalonados en la celda Sm.
- ``bloch``: forma k(ω) con regla inversa de factorización (método principal).
- ``dispersion``: R_PWE(ν) = Re cos(k Lm) empaquetado como ``DispersionScan``.
- ``eigenfrequency``: forma ω(k) del libro (benchmark no dispersivo y QEP de Drude).
"""
