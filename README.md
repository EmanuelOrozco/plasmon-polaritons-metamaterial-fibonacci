# Plasmon-polaritones en superredes Fibonacci (TMM y PWE)

Reimplementación computacional en Python del artículo de Reyes-Gómez *et al.*,
[Phys. Rev. B **81**, 153101 (2010)](https://doi.org/10.1103/PhysRevB.81.153101),
a cargo de **Emanuel Orozco Gallego**.

Las figuras del paper se calculan con dos métodos independientes:

- **Matriz de transferencia (TMM)** a incidencia oblicua: exacta para medios
  estratificados; es el método de producción.
- **Expansión en ondas planas (PWE)**, siguiendo el capítulo 1D de Sukhoivanov
  y Guryev (copia local en `referencias/`), en la forma $k(\omega)$ con la regla inversa de
  Li para tratar el metamaterial de Drude; sirve de verificación independiente.

No se usa FDTD ni elementos finitos.

## Documentación

| Documento | Descripción |
|-----------|-------------|
| [docs/ANEXO_DOCUMENTACION.md](docs/ANEXO_DOCUMENTACION.md) | Resumen de la guía conceptual, solución numérica, arquitectura y paper |
| [docs/RESULTADOS.md](docs/RESULTADOS.md) | Resultados frente al paper original (figuras 1–6) |
| [docs/INFORME_VALIDACION.md](docs/INFORME_VALIDACION.md) | Informe breve de validación |
| [docs/fase_A_auditoria_matematica.md](docs/fase_A_auditoria_matematica.md) | Auditoría matemática previa a la implementación |

### PDF incluidos en el repositorio

**Resultados de la simulación (TMM):**

- [figure_01.pdf](results/tmm/figure_01/output/figure_01.pdf)
- [figure_02.pdf](results/tmm/figure_02/output/figure_02.pdf)
- [figure_03.pdf](results/tmm/figure_03/output/figure_03.pdf)
- [figure_04.pdf](results/tmm/figure_04/output/figure_04.pdf)
- [figure_05.pdf](results/tmm/figure_05/output/figure_05.pdf)
- [figure_06.pdf](results/tmm/figure_06/output/figure_06.pdf)

Las mismas figuras con ondas planas están en `results/pwe/figure_0N/output/` y
las superposiciones y errores TMM vs PWE en `results/comparison/`.

**Documentos del proyecto:**

- [paper_reimplementation.pdf](docs/scientific_paper/paper_reimplementation.pdf)
- [software_documentation.pdf](docs/software_documentation/software_documentation.pdf)
- [guia_conceptual.pdf](docs/guia_conceptual/guia_conceptual.pdf)
- [solucion_numerica.pdf](docs/solucion_numerica/solucion_numerica.pdf)
- [arquitectura.pdf](docs/arquitectura/arquitectura.pdf)
- [metodo_ondas_planas.pdf](docs/metodo_ondas_planas/metodo_ondas_planas.pdf): paper del método de ondas planas
- [comparacion_tmm_pwe.pdf](docs/comparacion_tmm_pwe/comparacion_tmm_pwe.pdf): anexo de comparación TMM vs PWE

---

## Requisitos

- Python **3.10** o superior
- `pip` y `venv`
- Opcional: LaTeX (`pdflatex`) solo si quieres regenerar los PDF desde las fuentes `.tex`

## Instalación

```bash
git clone https://github.com/EmanuelOrozco/plasmon-polaritons-metamaterial-fibonacci.git
cd plasmon-polaritons-metamaterial-fibonacci

python3 -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install --upgrade pip
pip install -r requirements.txt
pip install -e .                   # opcional: instala el paquete fibonacci_photonics
```

Dependencias: `numpy`, `scipy`, `matplotlib`, `pyyaml`, `pytest`.

## Cómo ejecutar

### Todo de una vez

```bash
make all            # tests + resultados (TMM, PWE, comparación) + PDF + verificación
make reproduce      # lo mismo sin los tests
make verify         # solo la verificación de resultados y documentos
```

`make verify` (`scripts/verify_results.py`) comprueba que existan todas las
salidas, que TMM y PWE cumplan los criterios de acuerdo (conteos $F_{m-2}$,
error de semitraza, error de borde y de ancho) y que cada PDF esté al día con
sus figuras y tablas y compile sin errores. Con `--reference DIR` compara además
los datos con una corrida anterior.

### Tests

```bash
make tests          # o: .venv/bin/python -m pytest
```

### Todos los resultados

```bash
make results        # TMM + PWE + comparación
```

Por método:

```bash
.venv/bin/python scripts/tmm/run_all.py          # results/tmm/        (segundos)
.venv/bin/python scripts/pwe/run_all.py          # results/pwe/        (horas; PWE_WORKERS=n limita procesos)
.venv/bin/python scripts/comparison/run_all.py   # results/comparison/ (requiere los dos anteriores)
```

`make figures` sigue reproduciendo las figuras TMM (incluye la convergencia).

### Una figura por aparte

```bash
.venv/bin/python scripts/run_figure.py 3 --method tmm     # TMM
.venv/bin/python scripts/run_figure.py 3 --method pwe     # PWE
.venv/bin/python scripts/run_figure.py all --method both  # las seis, ambos métodos
```

`scripts/tmm/figure_0N.py` y `scripts/pwe/figure_0N.py` siguen disponibles como envoltorios
de una línea sobre el mismo runner.

Salidas:

- Imágenes: `results/<método>/figure_0N/output/` (PDF y PNG; la TMM también SVG)
- Datos y resumen: `results/<método>/figure_0N/data/`
- Parámetros físicos: `configs/figure_0N.yaml` (comunes a los dos métodos)
- Parámetros numéricos del PWE: `configs/pwe/`
- Estudio de tiempos TMM vs PWE (un núcleo): `configs/comparison/efficiency.yaml`,
  `scripts/comparison/efficiency.py` → `results/comparison/efficiency/`. Conclusión: la TMM es
  el método óptimo en 1D (exacta; 0.27 s frente a ≈5 h del PWE para las seis figuras).

### Regenerar documentación PDF (opcional)

```bash
make paper      # artículo de la reimplementación
make docs       # documentación del software
make guides     # guía conceptual, solución numérica, arquitectura
make pwe-docs   # paper del PWE y anexo de comparación
make documents  # todos los anteriores
```

Cada PDF depende de su `.tex` y de las figuras y tablas que incluye: `make` solo
recompila lo que cambió y se detiene en el primer error de LaTeX.

## Qué hace cada figura

| Figura | Contenido |
|--------|-----------|
| 1 | Dispersión TE $\nu(k)$ para $S_3$ y $S_4$, $\nu_e=\nu_m=3$ GHz, varios $\theta$ |
| 2 | Zoom cerca de $\nu_m=3$ GHz; número de subbandas $= F_{m-2}$ |
| 3 | Como Fig. 1 con $\nu_m=1$ GHz dentro del gap $\langle n\rangle=0$ |
| 4 | Zoom de modos plasmon-polaritón ($m=3$: 1 subbanda; $m=4$: 2) |
| 5 | Fragmentación de bandas vs orden de Fibonacci |
| 6 | Bandas vs ángulo $\theta$ (parámetros inferidos del resto del paper) |

Detalle y métricas: [docs/RESULTADOS.md](docs/RESULTADOS.md).

## Cambiar parámetros

Edita el YAML de la figura y vuelve a ejecutar su script. Ejemplo en `configs/figure_03.yaml`:

```yaml
omega_m_over_2pi_ghz: 1.5          # mueve el plasmón magnético
layer_a_thickness_mm: 10.0         # espesor del bloque A
fibonacci_orders: [3, 4, 5]        # más órdenes Fibonacci
frequency_points: 20000            # más resolución en frecuencia
```

## Estructura del repositorio

```text
configs/                     YAML por figura (física) + configs/pwe/ (numérica del PWE)
                             + configs/comparison/ (estudio de eficiencia)
src/fibonacci_photonics/     biblioteca (v2.0)
  physics/                   constantes, unidades, materiales de Drude, Fibonacci, SuperlatticeSpec
  solvers/                   contrato DispersionSolver, DispersionScan, intervalos permitidos
    tmm/                     matriz de transferencia y TMMSolver
    pwe/                     ondas planas: Fourier de la celda, solver k(ω), ω(k), convergencia, PWESolver
  config/                    esquema Pydantic v2 y carga de YAML
  analysis/                  modos de plasmón, anchos de banda, bordes, cierre de bandas, raíces
  benchmark/                 comparación TMM vs PWE, tiempos y eficiencia
  reproduction/              pipeline de figuras y CLI de scripts/run_figure.py
  studies/                   convergencia, benchmark del libro, contaminación espectral
  io/, viz/                  rutas, resultados, procedencia, LaTeX; estilo de las figuras
scripts/run_figure.py        runner único: figuras 1–6 con --method tmm|pwe|both
scripts/tmm/                 figuras 1–6 y convergencia con la TMM
scripts/pwe/                 benchmark del libro, estudios numéricos y figuras 1–6 con el PWE
scripts/comparison/          superposiciones, errores, tablas LaTeX y eficiencia (tiempo, memoria)
scripts/verify_results.py    verificación automática (make verify)
results/tmm/                 salidas TMM
results/pwe/                 salidas PWE
results/comparison/          comparación TMM vs PWE y tablas
tests/                       pytest: unit/, analytic/ (casos cerrados), regression/ (valores de referencia)
docs/                        LaTeX, PDF y anexos Markdown
presentacion/                presentación (independiente; no se regenera con make results)
notebooks/                   exploración / validación (opcional)
referencias/                 solo local (no se sube): paper original y material del PWE del libro
```

### Qué se sube al repositorio

El `.gitignore` deja fuera todo lo que se regenera o es local:

- Se sube: código (`src/`, `scripts/`, `tests/`), configuración (`configs/`),
  figuras y resúmenes (`results/**/output/`, `summary.json`, `results/comparison/tables/`),
  los barridos PWE (`results/pwe/**/data/*.npz`, horas de cómputo que usa la comparación),
  las fuentes y los PDF de `docs/` y `presentacion/`.
- No se sube: entorno y cachés, auxiliares de LaTeX,
  barridos TMM (`results/tmm/**/data/*.npz`, segundos de cómputo), metadatos de
  corrida y el material de partida externo (`referencias/`: paper de APS, capítulo,
  cuaderno y `PWE_2D.m` del libro, con copyright de sus autores).

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
