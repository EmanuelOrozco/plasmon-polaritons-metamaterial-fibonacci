# Plasmon-polaritones en superredes Fibonacci (TMM)

Reimplementación computacional en Python del artículo de Reyes-Gómez *et al.*,
[Phys. Rev. B **81**, 153101 (2010)](https://doi.org/10.1103/PhysRevB.81.153101),
a cargo de **Emanuel Orozco Gallego**.

El método es **matriz de transferencia (TMM)** a incidencia oblicua.
No se usa FDTD ni elementos finitos.

## Documentación

| Documento | Descripción |
|-----------|-------------|
| [docs/ANEXO_DOCUMENTACION.md](docs/ANEXO_DOCUMENTACION.md) | Resumen de la guía conceptual, solución numérica, arquitectura y paper |
| [docs/RESULTADOS.md](docs/RESULTADOS.md) | Resultados frente al paper original (figuras 1–6) |
| [docs/INFORME_VALIDACION.md](docs/INFORME_VALIDACION.md) | Informe breve de validación |
| [docs/fase_A_auditoria_matematica.md](docs/fase_A_auditoria_matematica.md) | Auditoría matemática previa a la implementación |

### PDF incluidos en el repositorio

**Resultados de la simulación (figuras):**

- [figure_01.pdf](figures/figure_01/output/figure_01.pdf)
- [figure_02.pdf](figures/figure_02/output/figure_02.pdf)
- [figure_03.pdf](figures/figure_03/output/figure_03.pdf)
- [figure_04.pdf](figures/figure_04/output/figure_04.pdf)
- [figure_05.pdf](figures/figure_05/output/figure_05.pdf)
- [figure_06.pdf](figures/figure_06/output/figure_06.pdf)

**Documentos del proyecto:**

- [paper_reimplementation.pdf](docs/scientific_paper/paper_reimplementation.pdf)
- [software_documentation.pdf](docs/software_documentation/software_documentation.pdf)
- [guia_conceptual.pdf](docs/guia_conceptual/guia_conceptual.pdf)
- [solucion_numerica.pdf](docs/solucion_numerica/solucion_numerica.pdf)
- [arquitectura.pdf](docs/arquitectura/arquitectura.pdf)

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
pip install -e .                   # opcional: instala el paquete fibonacci_tmm
```

Dependencias: `numpy`, `scipy`, `matplotlib`, `pyyaml`, `pytest`.

## Cómo ejecutar

### Tests

```bash
.venv/bin/pytest
```

### Todas las figuras

```bash
.venv/bin/python scripts/reproduce_all.py
```

También: `make figures` (incluye el análisis de convergencia).

### Una figura por aparte

```bash
.venv/bin/python scripts/reproduce_figure_01.py
.venv/bin/python scripts/reproduce_figure_02.py
.venv/bin/python scripts/reproduce_figure_03.py
.venv/bin/python scripts/reproduce_figure_04.py
.venv/bin/python scripts/reproduce_figure_05.py
.venv/bin/python scripts/reproduce_figure_06.py
```

Salidas:

- Imágenes: `figures/figure_0N/output/` (PNG, SVG, PDF)
- Datos y resumen: `figures/figure_0N/data/`
- Parámetros: `configs/figure_0N.yaml`

### Regenerar documentación PDF (opcional)

```bash
make paper    # artículo de la reimplementación
make docs     # documentación del software
make guides   # guía conceptual, solución numérica, arquitectura
make all      # tests + figuras + todos los PDF
```

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
configs/             YAML por figura (+ base)
src/fibonacci_tmm/   biblioteca TMM
scripts/             reproducción de figuras
figures/figure_0N/   output/ (PNG, SVG, PDF), data/, README
tests/               pytest
docs/                LaTeX, PDF y anexos Markdown
notebooks/           exploración / validación (opcional)
```

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
