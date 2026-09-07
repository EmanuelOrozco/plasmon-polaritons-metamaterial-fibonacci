# Informe final de validación

**Autor:** Emanuel Orozco Gallego
**Fecha:** 31 de agosto de 2026
**Artículo:** Reyes-Gómez et al., Phys. Rev. B 81, 153101 (2010)

## Tabla de figuras

| Figura | Reproducida | Concordancia | Diferencia principal | Causa |
| ------ | ----------- | ------------ | -------------------- | ----- |
| Fig. 1 | Sí | Buena concordancia | Sin gaps a θ=0 | Impedancia adaptada si ωe=ωm; las curvas ya no se rellenan a baja ν (cortes en saltos de k) |
| Fig. 2 | Sí | Excelente concordancia | — | N_modos = F_{m-2} = 1,2,3,5 exacto |
| Fig. 3 | Sí | Buena concordancia | Comparación visual | ν_m dentro del gap ⟨n⟩=0; el original no lleva línea de referencia en 1 GHz |
| Fig. 4 | Sí | Excelente concordancia | Δν ~20 % vs lectura visual del impreso | Lectura a ojo, no dato digital |
| Fig. 5 | Sí | Excelente concordancia | El PRB rota el eje como “bandwidth (GHz)” | Los valores son frecuencias de las subbandas (0.994–1.000 y 0.95–1.00 GHz) |
| Fig. 6 | Sí (parámetros inferidos) | Buena concordancia | Caption incompleto | Misma gráfica que el original: ν(θ) relleno, no Δν en líneas; m=2…7 inferido |

## Correcciones visuales respecto a una primera versión

- **Fig. 1:** `plot()` unía puntos consecutivos en ω y cruzaba la zona de Brillouin con rectas, lo que manchaba el panel a baja frecuencia. Ahora se inserta NaN si \|Δk\| es grande.
- **Fig. 6:** la primera versión dibujaba Δν(θ) como líneas desde 0. El original son seis paneles de regiones rellenas de frecuencia vs θ, que nacen en ν_m=1 GHz a θ=0 y se abren hacia abajo. Colores de menor a mayor ν: negro, rojo, verde, azul, magenta, cian.

## Entorno

- Python 3.14.6
- numpy 2.5.2, scipy 1.18.1, matplotlib 3.11.1, pyyaml 6.0.3, pytest 9.1.1
- Tolerancia de banda: |R| ≤ 1 + 10^{-12}
- c = 299792458 m/s (IMPL; el paper no lo da)
- q = (ω/c) n_A sin θ (ω/c restaurado)
- Resolución Fig. 1: 12000 puntos en 0.15–5 GHz
- Resolución Fig. 4: 16000 puntos en cada ventana
- Convergencia Fig. 4, m=3, θ=π/3: Δν = 0.03310 GHz (16000 pts); estable desde 2000 pts (<0.2 %)
- Tests: todos los ejecutados pasan (3 skipped por k>m)

## Parámetros del paper usados

- a = b = 12 mm, εA = μA = 1, ε0 = μ0 = 1
- Figs. 1–2: ωe/2π = ωm/2π = 3 GHz
- Figs. 3–5: ωe/2π = 3 GHz, ωm/2π = 1 GHz
- Fig. 6: mismos que 3–5, m=2…7 (inferido)
