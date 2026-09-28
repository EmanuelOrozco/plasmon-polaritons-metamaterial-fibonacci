# Presentaciones

Dos presentaciones Beamer (16:9, tema `metropolis`):

| Presentación | Fuente | PDF | Enfoque |
|---|---|---|---|
| Reimplementación del paper | [`presentacion.tex`](presentacion.tex) | [`presentacion.pdf`](presentacion.pdf) | Física, matriz de transferencia (TMM) y las seis figuras de Reyes-Gómez *et al.*, Phys. Rev. B **81**, 153101 (2010) |
| Método de ondas planas | [`presentacion_pwe.tex`](presentacion_pwe.tex) | [`presentacion_pwe.pdf`](presentacion_pwe.pdf) | Conceptos del PWE, por qué falla la receta del libro con Drude, la forma $k(\omega)$ con la regla inversa de Li y su implementación en el paquete |

## Contenido: reimplementación del paper

| Sección | Temas |
|---|---|
| El problema | Artículo original, qué significa reimplementar, geometría de la superred |
| Conceptos físicos | $\varepsilon$ y $\mu$, metamateriales, modelo de Drude, polarización TE, plasmon-polaritones, bandas y gaps, gap $\langle n\rangle = 0$, secuencia de Fibonacci, espectro tipo Cantor |
| Cómo se reimplementó | Metodología, auditoría del paper, ecuación de onda, matriz de transferencia, semitraza y bandas, validación cruzada, arquitectura, parámetros, tests y convergencia |
| Resultados | Figuras 1 a 6 explicadas, conteo de modos $N = F_{m-2}$, comparación con el original, dos regímenes físicos |
| Cierre | Limitaciones, conclusiones, cómo reproducir, referencias |

## Contenido: método de ondas planas

| Sección | Temas |
|---|---|
| Motivación | Por qué un segundo método, hoja de ruta |
| Conceptos necesarios | Ecuación TE/TM unificada, teorema de Bloch, red recíproca y zona de Brillouin, series de Fourier, coeficientes de la celda de Fibonacci, matrices de Toeplitz, ecuación maestra, algoritmo $\omega(k)$ del libro |
| Validación con el libro | Bicapa $\varepsilon_1=1$ / $\varepsilon_2=9$, convergencia Laurent vs inversa, perfil de campo y bandas fuera del eje |
| El reto: Drude | Problema cuadrático en $\omega^2$ y contaminación espectral junto a $\nu_m$ |
| La solución | Forma $k(\omega)$, reglas de factorización de Li, QEP en $k$ y matriz compañera, selección del modo, convergencia y límite de validez |
| Implementación | Arquitectura de `solvers/pwe/`, contrato común con la TMM, núcleo del solver, localización de bordes, parámetros y pruebas |
| Resultados y cierre | Figuras 1 a 6 con el PWE, superposición con la TMM, conteo de subbandas, costo, cuándo usar cada método |

## Estructura

```
presentacion/
├── presentacion.tex             # diapositivas del paper (TMM)
├── presentacion.pdf
├── presentacion_pwe.tex         # diapositivas del método de ondas planas
├── presentacion_pwe.pdf
├── imagenes/                    # imágenes conceptuales de presentacion.tex (PNG y PDF)
├── imagenes_pwe/                # imágenes conceptuales de presentacion_pwe.tex (PNG y PDF)
└── scripts/
    ├── generar_imagenes.py      # genera imagenes/
    └── generar_imagenes_pwe.py  # genera imagenes_pwe/
```

Las figuras de resultados se toman directamente de `../results/`: `presentacion.tex` usa `results/tmm/figure_0N/output/`, y `presentacion_pwe.tex` usa `results/pwe/` y `results/comparison/` (figuras, tablas y macros de eficiencia). Conviene regenerarlas (`make tmm`, `make pwe`, `make comparison`) antes de compilar si se cambió algún YAML de `configs/`.

## Regenerar

Desde la raíz del repositorio, con el entorno virtual activo y una distribución LaTeX con Beamer y el tema `metropolis`:

```bash
make presentation       # presentacion.pdf
make presentation-pwe   # presentacion_pwe.pdf
```

O paso a paso:

```bash
MPLCONFIGDIR=/tmp/mpl .venv/bin/python presentacion/scripts/generar_imagenes.py
MPLCONFIGDIR=/tmp/mpl .venv/bin/python presentacion/scripts/generar_imagenes_pwe.py
cd presentacion
pdflatex presentacion.tex && pdflatex presentacion.tex
pdflatex presentacion_pwe.tex && pdflatex presentacion_pwe.tex
```

La segunda pasada de `pdflatex` es necesaria para el índice y la numeración de diapositivas. `make clean` borra los archivos auxiliares.

`imagenes/metamaterial_ilustracion.png` es una ilustración esquemática fija; el script no la regenera.
