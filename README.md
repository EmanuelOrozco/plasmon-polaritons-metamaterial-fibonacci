# Plasmon-polaritones en superredes Fibonacci (TMM)

Reimplementación computacional en Python del artículo de Reyes-Gómez *et al.*,
[Phys. Rev. B **81**, 153101 (2010)](https://doi.org/10.1103/PhysRevB.81.153101),
a cargo de **Emanuel Orozco Gallego**.

El método es **matriz de transferencia (TMM)** a incidencia oblicua.
No se usa FDTD ni elementos finitos.

| Documento | Descripción |
|-----------|-------------|
| [docs/ANEXO_DOCUMENTACION.md](docs/ANEXO_DOCUMENTACION.md) | Resumen detallado de la guía conceptual, solución numérica, arquitectura y paper de reimplementación |
| [docs/RESULTADOS.md](docs/RESULTADOS.md) | Resultados obtenidos frente al paper original (figuras 1–6) |
| [docs/INFORME_VALIDACION.md](docs/INFORME_VALIDACION.md) | Informe breve de validación |
| [docs/fase_A_auditoria_matematica.md](docs/fase_A_auditoria_matematica.md) | Auditoría matemática previa a la implementación |

---

## Requisitos

- Python **≥ 3.10**
- `pip` y `venv`
- Opcional: LaTeX (`pdflatex`) solo si quieres regenerar los PDF en `docs/`

## Instalación (otro computador)

```bash
git clone <URL-DE-ESTE-REPO>.git
cd analisis_numerico

python3 -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install --upgrade pip
pip install -r requirements.txt
pip install -e .                   # opcional: instala el paquete `fibonacci_tmm`
```

Dependencias principales: `numpy`, `scipy`, `matplotlib`, `pyyaml`, `pytest`.

## Cómo ejecutar

### Tests

```bash
.venv/bin/pytest
# o, con el venv activo:
pytest
```

### Todas las figuras

```bash
.venv/bin/python scripts/reproduce_all.py
```

Equivale a `make figures` (también corre el análisis de convergencia).

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

### Documentación PDF (opcional)

```bash
make paper    # artículo de la reimplementación
make docs     # documentación del software
make guides   # guía conceptual, solución numérica, arquitectura
make all      # tests + figuras + todos los PDF
```

Los PDF se generan desde las fuentes `.tex` en `docs/`. El contenido narrativo también está resumido en Markdown (ver tabla al inicio).

## Qué hace cada figura

| Figura | Contenido |
|--------|-----------|
| 1 | Dispersión TE \(ν(k)\) para \(S_3\) y \(S_4\), \(\nu_e=\nu_m=3\) GHz, varios \(\theta\) |
| 2 | Zoom cerca de \(\nu_m=3\) GHz; número de subbandas \(= F_{m-2}\) |
| 3 | Como Fig. 1 con \(\nu_m=1\) GHz dentro del gap \(\langle n\rangle=0\) |
| 4 | Zoom de modos plasmon-polaritón (\(m=3\): 1 subbanda; \(m=4\): 2) |
| 5 | Fragmentación de bandas vs orden de Fibonacci |
| 6 | Bandas vs ángulo \(\theta\) (parámetros inferidos del resto del paper) |

Detalle y métricas: [docs/RESULTADOS.md](docs/RESULTADOS.md).

## Cambiar parámetros

Edita el YAML de la figura y vuelve a ejecutar su script. Ejemplos en `configs/figure_03.yaml`:

```yaml
omega_m_over_2pi_ghz: 1.5          # mueve el plasmón magnético
layer_a_thickness_mm: 10.0         # espesor del bloque A
fibonacci_orders: [3, 4, 5]        # más órdenes Fibonacci
frequency_points: 20000            # más resolución en frecuencia
```

## Estructura del repositorio

```
configs/             YAML por figura (+ base)
src/fibonacci_tmm/   biblioteca TMM
scripts/             reproducción de figuras
figures/figure_0N/   output/, data/, README
tests/               pytest
docs/                LaTeX + anexos Markdown
notebooks/           exploración / validación (opcional)
```

## Convención física importante

El paper escribe \(q = n_A\sin\theta\), dimensionalmente inconsistente con
\(Q=\sqrt{(\omega/c)^2 n^2-q^2}\). Aquí se usa

\[
q=\frac{\omega}{c}\,n_A\sin\theta
\]

como en el artículo hermano del mismo grupo (EPL **88**, 24002, 2009).
Velocidad de la luz: valor SI \(c=299\,792\,458\,\mathrm{m/s}\) (el PRB 2010 no la especifica).

## Licencia y artículo original

Código y documentación propia: **MIT** (ver `LICENSE`).

El artículo de *Physical Review B* es de la American Physical Society.
**No** se incluye el PDF ni capturas del original en este repositorio;
consúltalo vía DOI o acceso institucional.
