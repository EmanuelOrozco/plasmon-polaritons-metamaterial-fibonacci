# Presentación

Diapositivas (Beamer, 16:9) que explican la reimplementación de Reyes-Gómez *et al.*, Phys. Rev. B **81**, 153101 (2010): los conceptos físicos, el método numérico y los resultados de las seis figuras.

- PDF: [`presentacion.pdf`](presentacion.pdf)
- Fuente: [`presentacion.tex`](presentacion.tex)

## Contenido

| Sección | Temas |
|---|---|
| El problema | Artículo original, qué significa reimplementar, geometría de la superred |
| Conceptos físicos | $\varepsilon$ y $\mu$, metamateriales, modelo de Drude, polarización TE, plasmon-polaritones, bandas y gaps, gap $\langle n\rangle = 0$, secuencia de Fibonacci, espectro tipo Cantor |
| Cómo se reimplementó | Metodología, auditoría del paper, ecuación de onda, matriz de transferencia, semitraza y bandas, validación cruzada, arquitectura, parámetros, tests y convergencia |
| Resultados | Figuras 1 a 6 explicadas, conteo de modos $N = F_{m-2}$, comparación con el original, dos regímenes físicos |
| Cierre | Limitaciones, conclusiones, cómo reproducir, referencias |

## Estructura

```
presentacion/
├── presentacion.tex          # diapositivas
├── presentacion.pdf          # resultado compilado
├── imagenes/                 # imágenes conceptuales (PNG y PDF)
└── scripts/
    └── generar_imagenes.py   # genera todo lo que hay en imagenes/
```

Las figuras de resultados 1 a 6 se toman directamente de `../results/tmm/figure_0N/output/`, así que conviene regenerarlas (`make tmm`) antes de compilar si se cambió algún YAML de `configs/`.

## Regenerar

Desde la raíz del repositorio, con el entorno virtual activo y una distribución LaTeX con Beamer y el tema `metropolis`:

```bash
make presentation
```

O paso a paso:

```bash
MPLCONFIGDIR=/tmp/mpl .venv/bin/python presentacion/scripts/generar_imagenes.py
cd presentacion
pdflatex presentacion.tex
pdflatex presentacion.tex
```

La segunda pasada de `pdflatex` es necesaria para el índice y la numeración de diapositivas. `make clean` borra los archivos auxiliares.

`imagenes/metamaterial_ilustracion.png` es una ilustración esquemática fija; el script no la regenera.
