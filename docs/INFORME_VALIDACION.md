# Informe final de validación

**Autor:** Emanuel Orozco Gallego  
**Fecha:** 31 de agosto de 2026  
**Artículo:** Reyes-Gómez et al., Phys. Rev. B 81, 153101 (2010)

## Tabla de figuras

| Figura | Reproducida | Concordancia | Diferencia principal | Causa |
|--------|-------------|--------------|----------------------|-------|
| Fig. 1 | Sí | Buena | Sin gaps a $\theta=0$ | Impedancia adaptada si $\omega_e=\omega_m$; las curvas ya no se rellenan a baja $\nu$ (cortes en saltos de $k$) |
| Fig. 2 | Sí | Excelente | — | $N_{\mathrm{modos}} = F_{m-2} = 1,2,3,5$ exacto |
| Fig. 3 | Sí | Buena | Comparación visual | $\nu_m$ dentro del gap $\langle n\rangle=0$; el original no lleva línea de referencia en 1 GHz |
| Fig. 4 | Sí | Excelente | $\Delta\nu$ ~20 % vs lectura visual del impreso | Lectura a ojo, no dato digital |
| Fig. 5 | Sí | Excelente | El PRB rotula el eje como “bandwidth (GHz)” | Los valores son frecuencias de las subbandas (0.994–1.000 y 0.95–1.00 GHz) |
| Fig. 6 | Sí (parámetros inferidos) | Buena | Caption incompleto | Misma gráfica que el original: $\nu(\theta)$ relleno, no $\Delta\nu$ en líneas; $m=2\ldots 7$ inferido |

## Correcciones visuales respecto a una primera versión

- **Fig. 1:** `plot()` unía puntos consecutivos en $\omega$ y cruzaba la zona de Brillouin con rectas, lo que manchaba el panel a baja frecuencia. Ahora se inserta NaN si $|\Delta k|$ es grande.
- **Fig. 6:** la primera versión dibujaba $\Delta\nu(\theta)$ como líneas desde 0. El original son seis paneles de regiones rellenas de frecuencia vs $\theta$, que nacen en $\nu_m=1$ GHz a $\theta=0$ y se abren hacia abajo. Colores de menor a mayor $\nu$: negro, rojo, verde, azul, magenta, cian.

## Entorno

- Python 3.14.6
- numpy 2.5.2, scipy 1.18.1, matplotlib 3.11.1, pyyaml 6.0.3, pytest 9.1.1
- Tolerancia de banda: $|R| \le 1 + 10^{-12}$
- $c = 299792458$ m/s (`IMPL`; el paper no lo da)
- $q = (\omega/c)\, n_A \sin\theta$ ($\omega/c$ restaurado)
- Resolución Fig. 1: 12000 puntos en 0.15–5 GHz
- Resolución Fig. 4: 16000 puntos en cada ventana
- Convergencia Fig. 4, $m=3$, $\theta=\pi/3$: $\Delta\nu = 0.03310$ GHz (16000 pts); estable desde 2000 pts ($<0.2$ %)
- Tests: todos los ejecutados pasan (3 skipped por $k>m$)

## Parámetros del paper usados

- $a = b = 12$ mm, $\varepsilon_A = \mu_A = 1$, $\varepsilon_0 = \mu_0 = 1$
- Figs. 1–2: $\omega_e/2\pi = \omega_m/2\pi = 3$ GHz
- Figs. 3–5: $\omega_e/2\pi = 3$ GHz, $\omega_m/2\pi = 1$ GHz
- Fig. 6: mismos que 3–5, $m=2\ldots 7$ (inferido)
