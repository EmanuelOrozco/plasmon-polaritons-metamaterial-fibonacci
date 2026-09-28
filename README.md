# Plasmon-polaritones en superredes Fibonacci (TMM y PWE)

Reimplementación computacional en Python del artículo de Reyes-Gómez *et al.*,
[Phys. Rev. B **81**, 153101 (2010)](https://doi.org/10.1103/PhysRevB.81.153101),
a cargo de **Emanuel Orozco Gallego**.

## Cómo funciona

El sistema es una superred unidimensional que alterna aire (A) y un
metamaterial de Drude (B), con $\varepsilon_B$ y $\mu_B$ negativos por debajo de sus
frecuencias de plasma, apilados según la secuencia de Fibonacci
$S_m=S_{m-1}S_{m-2}$. La luz incide con ángulo $\theta$ y polarización TE. Una
frecuencia $\nu$ pertenece a una banda permitida cuando la semitraza de la celda
cumple $\lvert R_m(\nu)\rvert\le 1$, y entonces $\cos(kL_m)=R_m$ da la relación de
dispersión $\nu(k)$. Las seis figuras del paper salen de este criterio: curvas
$\nu(k)$, subbandas plasmon-polaritón cerca de $\nu_m$ (hay $F_{m-2}$) y su
evolución con el orden $m$ y con $\theta$.

$R_m$ se calcula con dos métodos independientes que implementan el mismo contrato
(`DispersionSolver`), de modo que el post-proceso y las figuras no dependen del
método:

- **Matriz de transferencia (TMM)**: producto exacto de matrices de capa, con la
  recurrencia de trazas de Fibonacci. Cuesta microsegundos por frecuencia y es el
  método de producción.
- **Expansión en ondas planas (PWE)**: sigue el capítulo 1D de Sukhoivanov y
  Guryev, reformulado como problema cuadrático en $k(\omega)$ con la regla inversa
  de Li para tratar el Drude. Es la verificación independiente.

Flujo de una corrida:

```text
configs/physics/figure_0N.yaml  +  configs/tmm/ o configs/pwe/figure_0N.yaml
        │   (validación con Pydantic)
        ▼
SuperlatticeSpec ──► TMMSolver | PWESolver ──► barrido R_m(ν) ──► bordes exactos |R|=1
        ▼
results/<tmm|pwe>/figure_0N/{output,data}  ──►  results/comparison/  ──►  PDF en docs/
```

Los PDF no recalculan física: incluyen las figuras y tablas ya generadas en `results/`.

## Documentación en PDF

| Documento | Contenido |
|-----------|-----------|
| [Paper de reimplementación](docs/scientific_paper/paper_reimplementation.pdf) | Artículo de la reproducción: modelo, TMM, las seis figuras, concordancia con el PRB, convergencia y limitaciones |
| [Guía conceptual](docs/guia_conceptual/guia_conceptual.pdf) | Física desde cero: $\varepsilon$ y $\mu$, metamateriales, plasmon-polaritones, gap $\langle n\rangle=0$, Fibonacci, qué muestra cada figura |
| [Solución numérica](docs/solucion_numerica/solucion_numerica.pdf) | De las ecuaciones al código de la TMM: unidades, matriz de capa, semitraza, barridos, detección de modos |
| [Método de ondas planas](docs/metodo_ondas_planas/metodo_ondas_planas.pdf) | Formulación del PWE, validación con el cristal del libro, contaminación espectral con Drude, forma $k(\omega)$ y resultados |
| [Comparación TMM vs PWE](docs/comparacion_tmm_pwe/comparacion_tmm_pwe.pdf) | Superposiciones, error de semitraza, bordes, conteos de subbandas, tiempos y eficiencia |
| [Arquitectura y ejecución](docs/arquitectura/arquitectura.pdf) | Árbol del proyecto, capas, configuración, flujo, pruebas, dependencias y problemas frecuentes |
| [Documentación del software](docs/software_documentation/software_documentation.pdf) | Manual corto de la API de `fibonacci_photonics` |

**Presentaciones** (Beamer 16:9, tema `metropolis`):

| Presentación | Enfoque |
|--------------|---------|
| [Reimplementación del paper](presentacion/presentacion.pdf) | Conceptos físicos, TMM, las seis figuras y el conteo $N=F_{m-2}$ |
| [Método de ondas planas](presentacion/presentacion_pwe.pdf) | Bloch, Fourier y Toeplitz, por qué falla la receta del libro con Drude, forma $k(\omega)$, implementación y comparación con la TMM |

**En Markdown:** [resultados y validación](docs/resultados.md) (cifras de cada
figura, comparación entre métodos, convergencia, decisiones propias y entorno) y
[auditoría matemática](docs/auditoria_matematica.md) del artículo, previa al código.

Orden de lectura sugerido: guía conceptual → auditoría → solución numérica →
arquitectura → paper de reimplementación → método de ondas planas → comparación.

## Resultados

| Figura | Contenido | TMM | PWE | Comparación |
|--------|-----------|-----|-----|-------------|
| 1 | Dispersión TE $\nu(k)$ para $S_3$ y $S_4$, $\nu_e=\nu_m=3$ GHz, varios $\theta$ | [PDF](results/tmm/figure_01/output/figure_01.pdf) | [PDF](results/pwe/figure_01/output/figure_01.pdf) | [superposición](results/comparison/figure_01/output/figure_01_overlay.pdf) · [error](results/comparison/figure_01/output/figure_01_error.pdf) |
| 2 | Zoom cerca de $\nu_m=3$ GHz; subbandas $=F_{m-2}$ | [PDF](results/tmm/figure_02/output/figure_02.pdf) | [PDF](results/pwe/figure_02/output/figure_02.pdf) | [superposición](results/comparison/figure_02/output/figure_02_overlay.pdf) · [error](results/comparison/figure_02/output/figure_02_error.pdf) |
| 3 | Como Fig. 1 con $\nu_m=1$ GHz dentro del gap $\langle n\rangle=0$ | [PDF](results/tmm/figure_03/output/figure_03.pdf) | [PDF](results/pwe/figure_03/output/figure_03.pdf) | [superposición](results/comparison/figure_03/output/figure_03_overlay.pdf) · [error](results/comparison/figure_03/output/figure_03_error.pdf) |
| 4 | Zoom de los plasmon-polaritones ($m=3$: 1 subbanda; $m=4$: 2) | [PDF](results/tmm/figure_04/output/figure_04.pdf) | [PDF](results/pwe/figure_04/output/figure_04.pdf) | [superposición](results/comparison/figure_04/output/figure_04_overlay.pdf) · [error](results/comparison/figure_04/output/figure_04_error.pdf) |
| 5 | Fragmentación de bandas frente al orden de Fibonacci | [PDF](results/tmm/figure_05/output/figure_05.pdf) | [PDF](results/pwe/figure_05/output/figure_05.pdf) | [barras](results/comparison/figure_05/output/figure_05_comparison.pdf) |
| 6 | Bandas frente al ángulo $\theta$ (parámetros inferidos) | [PDF](results/tmm/figure_06/output/figure_06.pdf) | [PDF](results/pwe/figure_06/output/figure_06.pdf) | [superposición](results/comparison/figure_06/output/figure_06_overlay.pdf) · [bordes](results/comparison/figure_06/output/figure_06_edge_error.pdf) |

Las dos implementaciones coinciden: la semitraza difiere menos de $2\times10^{-2}$
en el peor punto (extremo de 0.15 GHz de la Fig. 1) y entre $10^{-7}$ y $10^{-5}$ en
casi todo el rango, los bordes de banda difieren menos de $10^{-4}$ GHz y los conteos
$F_{m-2}$ son iguales. La TMM hace las seis figuras en 0.31 s y el PWE en ≈5.2 h
([eficiencia](results/comparison/efficiency/output/efficiency.pdf)), así que la TMM
es el método óptimo en 1D. Detalle en [docs/resultados.md](docs/resultados.md).

## Requisitos

- Python **3.10** o superior, con `pip` y `venv`.
- Opcional: `make` y LaTeX (`pdflatex`, con Beamer y `metropolis` para las
  presentaciones) para regenerar los PDF.

## Instalación

```bash
git clone https://github.com/EmanuelOrozco/plasmon-polaritons-metamaterial-fibonacci.git
cd plasmon-polaritons-metamaterial-fibonacci

python3 -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install --upgrade pip
pip install -r requirements.txt    # = pip install -e ".[dev]"
```

Dependencias: `numpy`, `scipy`, `matplotlib`, `pyyaml`, `pydantic`; en desarrollo,
`pytest`, `ruff` y `mypy`.

## Cómo ejecutar

### Todo de una vez

```bash
make all            # pruebas + resultados (TMM, PWE, comparación) + PDF + verificación
make reproduce      # lo mismo sin las pruebas
```

### Por partes

```bash
make tests          # pytest (unit/, analytic/, regression/)
make lint           # ruff + mypy
make tmm            # results/tmm/        Figs. 1–6 y convergencia (≈1 min)
make pwe            # results/pwe/        Figs. 1–6 y estudios del PWE (horas; PWE_WORKERS=n limita procesos)
make comparison     # results/comparison/ superposiciones, errores, tablas y eficiencia (requiere tmm y pwe)
make documents      # los siete PDF de docs/
make presentation presentation-pwe   # las dos presentaciones
make verify         # comprueba salidas, acuerdo TMM/PWE y que los PDF estén al día
make clean          # borra auxiliares de LaTeX
```

`make verify` (`scripts/verify_results.py`) exige los conteos $F_{m-2}$ y cotas de
error de semitraza, bordes y anchos; con `--reference DIR` compara además los datos
con una corrida anterior. Cada PDF se recompila solo si cambió su `.tex` o alguna
figura o tabla que incluye.

### Una figura

```bash
.venv/bin/python scripts/run_figure.py 3 --method tmm     # TMM
.venv/bin/python scripts/run_figure.py 3 --method pwe     # PWE
.venv/bin/python scripts/run_figure.py all --method both  # las seis, con los dos métodos
```

Cada figura escribe `results/<método>/figure_0N/output/` (PDF y PNG; la TMM
también SVG) y `results/<método>/figure_0N/data/` (`.npz`, `summary.json`,
`metadata.json`).

### Cambiar parámetros

La física de cada figura está en `configs/physics/figure_0N.yaml` y la leen los dos
métodos; la resolución de cada método está aparte. Por ejemplo, para la Fig. 3:

```yaml
# configs/physics/figure_03.yaml
omega_m_over_2pi_ghz: 1.5          # mueve el plasmón magnético
layer_a_thickness_mm: 10.0         # espesor del bloque A
fibonacci_orders: [3, 4, 5]        # más órdenes de Fibonacci

# configs/tmm/figure_03.yaml
frequency_points: 20000            # más resolución en frecuencia
```

Después se vuelve a ejecutar la figura. Un valor inválido o una clave desconocida
detiene la carga con un `ConfigError` que indica el archivo y el campo.

## Estructura del repositorio

| | TMM | PWE | Común / comparación |
|---|---|---|---|
| Solver | `src/.../solvers/tmm/` | `src/.../solvers/pwe/` | `solvers/base.py`, `solvers/scan.py` (contrato y barrido) |
| Estudios | `src/.../studies/tmm/` | `src/.../studies/pwe/` | `src/.../comparison/` |
| Configuración | `configs/tmm/` | `configs/pwe/` | `configs/physics/`, `configs/comparison/` |
| Scripts | `scripts/tmm/` | `scripts/pwe/` | `scripts/comparison/`, `scripts/run_figure.py` |
| Resultados | `results/tmm/` | `results/pwe/` | `results/comparison/` |
| Pruebas | `tests/unit/test_tmm_*.py` | `tests/unit/test_pwe_*.py` | resto de `tests/` |
| Presentación | `presentacion/presentacion.tex` | `presentacion/presentacion_pwe.tex` | |

(`src/...` es `src/fibonacci_photonics`.)

```text
configs/
  physics/                 física de cada figura (común a TMM y PWE)
  tmm/                     mallas de frecuencia y ángulo de la TMM
  pwe/                     armónicos y mallas del PWE; estudios del PWE
  comparison/              estudio de eficiencia
src/fibonacci_photonics/
  physics/                 constantes, unidades, Drude, TE/TM, Fibonacci, SuperlatticeSpec
  solvers/                 contrato DispersionSolver y DispersionScan
    tmm/                   matriz de transferencia y TMMSolver
    pwe/                   Fourier de la celda, QEP de Bloch, ω(k), convergencia, PWESolver
  analysis/                raíces, bordes y cierre de bandas, subbandas, anchos vs θ
  config/                  esquemas Pydantic v2 y carga de YAML
  reproduction/            Figs. 1–6 con cualquiera de los dos métodos y CLI
  studies/tmm/             convergencia de la TMM
  studies/pwe/             benchmark del libro, celda de Fourier, contaminación espectral, convergencia
  comparison/              TMM vs PWE: errores, bordes, tiempos, eficiencia, tablas
  io/, viz/                rutas, resultados, procedencia, LaTeX; estilo de figuras
scripts/                   envoltorios delgados de la biblioteca y verify_results.py
tests/                     unit/, analytic/ (casos cerrados), regression/ (valores de referencia)
docs/                      PDF y fuentes LaTeX, resultados.md, auditoria_matematica.md
presentacion/              presentaciones; scripts/ genera imagenes/ e imagenes_pwe/
notebooks/                 exploración (opcional)
referencias/               solo local: paper original y material del PWE del libro
```

### Qué se sube al repositorio

- Se sube: código (`src/`, `scripts/`, `tests/`), configuración (`configs/`),
  figuras y resúmenes (`results/**/output/`, `summary.json`, `results/comparison/tables/`),
  los barridos del PWE (`results/pwe/**/data/*.npz`, horas de cómputo que usa la
  comparación) y las fuentes y PDF de `docs/` y `presentacion/`.
- No se sube: entorno y cachés, auxiliares de LaTeX, barridos de la TMM
  (`results/tmm/**/data/*.npz`, segundos de cómputo), metadatos de corrida y el
  material externo de `referencias/` (paper de APS, capítulo, cuaderno y `PWE_2D.m`
  del libro, con copyright de sus autores).

## Convención física importante

El paper escribe $q = n_A\sin\theta$, dimensionalmente inconsistente con
$Q=\sqrt{(\omega/c)^2 n^2-q^2}$. Aquí se usa

$$
q=\frac{\omega}{c}\,n_A\sin\theta
$$

como en el artículo hermano del mismo grupo (EPL **88**, 24002, 2009).
Velocidad de la luz: valor SI $c=299\,792\,458\,\mathrm{m/s}$ (el PRB 2010 no la especifica).

## Licencia y artículo original

Código y documentación propia: **MIT** (ver [`LICENSE`](LICENSE)).

El artículo de *Physical Review B* es de la American Physical Society.
**No** se incluye el PDF ni capturas del original en este repositorio;
consúltalo vía [DOI](https://doi.org/10.1103/PhysRevB.81.153101) o acceso institucional.
