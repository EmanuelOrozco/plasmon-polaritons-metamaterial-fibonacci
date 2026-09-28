# Anexo: documentación del proyecto

Este anexo resume el contenido de los documentos PDF/LaTeX del repositorio.
Las fuentes canónicas están en `docs/*/*.tex`. Los PDF ya están versionados;
también se pueden regenerar con `make paper`, `make docs`, `make guides` y `make pwe-docs`.

| Documento | Fuente | PDF |
|-----------|--------|-----|
| Guía conceptual | `docs/guia_conceptual/guia_conceptual.tex` | [guia_conceptual.pdf](guia_conceptual/guia_conceptual.pdf) |
| Solución numérica | `docs/solucion_numerica/solucion_numerica.tex` | [solucion_numerica.pdf](solucion_numerica/solucion_numerica.pdf) |
| Arquitectura y ejecución | `docs/arquitectura/arquitectura.tex` | [arquitectura.pdf](arquitectura/arquitectura.pdf) |
| Documentación del software | `docs/software_documentation/software_documentation.tex` | [software_documentation.pdf](software_documentation/software_documentation.pdf) |
| Paper de reimplementación | `docs/scientific_paper/paper_reimplementation.tex` | [paper_reimplementation.pdf](scientific_paper/paper_reimplementation.pdf) |
| Método de ondas planas (PWE) | `docs/metodo_ondas_planas/metodo_ondas_planas.tex` | [metodo_ondas_planas.pdf](metodo_ondas_planas/metodo_ondas_planas.pdf) |
| Comparación TMM vs PWE | `docs/comparacion_tmm_pwe/comparacion_tmm_pwe.tex` | [comparacion_tmm_pwe.pdf](comparacion_tmm_pwe/comparacion_tmm_pwe.pdf) |
| Auditoría matemática | [fase_A_auditoria_matematica.md](fase_A_auditoria_matematica.md) | solo Markdown |
| Informe de validación | [INFORME_VALIDACION.md](INFORME_VALIDACION.md) | solo Markdown |
| Resultados (figs. 1–6) | [RESULTADOS.md](RESULTADOS.md) | PDF por figura en `results/tmm/*/output/` y `results/pwe/*/output/` |

---

## 1. Guía conceptual (`guia_conceptual`)

Texto de física desde cero para entender qué se está replicando, sin asumir
experiencia previa en cristales fotónicos.

**Temas principales**

1. **Qué replica el proyecto** — Brief Report PRB 2010; polarización TE; TMM.
2. **Ondas, $\varepsilon$ y $\mu$** — respuesta electromagnética de los medios.
3. **Polarizaciones TE y TM** — en el paper las figuras son TE; TM se menciona por dualidad.
4. **Metamateriales e índice negativo** — régimen $\varepsilon_B<0$, $\mu_B<0$.
5. **Plasmones y plasmon-polaritones** — acoplamiento luz–plasmón magnético a $\theta\neq 0$.
6. **Cristales fotónicos y dispersión** — bandas permitidas vía $\lvert R_m\rvert\le 1$.
7. **Gap de índice promedio nulo $\langle n\rangle=0$** — gap no Bragg, casi invariante de escala.
8. **Secuencias de Fibonacci** — $S_m=S_{m-1}S_{m-2}$, $S_0=B$, $S_1=A$, $F_0=F_1=1$.
9. **Geometría numérica** — $a=b=12\,\mathrm{mm}$; A = aire; B = Drude.
10. **$q$ y $Q$** — componente transversal y vector de onda en cada capa.
11. **Helmholtz y vector continuo $\Psi$** — continuidad de campo e impedancia.
12. **Matriz de transferencia** — producto de matrices de capa → semitraza $R_m=\cos(kL_m)$.
13. **Qué muestra cada figura (1–6)** — interpretación física panel a panel.
14. **Dos regímenes** — $\nu_m$ en banda permitida (Figs. 1–2) vs en el gap (Figs. 3–6).
15. **Lo que el paper no especifica** — $c$, forma explícita de la matriz de capa, caption de Fig. 6.
16. **Glosario** — términos usados en el resto de la documentación.

**Para quién:** lectores que necesitan el marco físico antes de tocar el código.

---

## 2. Solución numérica e implementación (`solucion_numerica`)

Puente entre las ecuaciones del paper y el código TMM en `src/fibonacci_photonics/`.

**Temas principales**

1. **Flujo global** — YAML → `SuperlatticeSpec` → barrido en $\nu$ → detección de modos → plot.
2. **Arquitectura del código** — módulos: `physics/`, `solvers/tmm/`, `solvers/scan.py`, `analysis/`, `viz/`.
3. **Unidades** — espesores en mm → SI; frecuencias de plasma como $\omega/2\pi$ en GHz.
4. **YAML → especificación** — `config.load_figure` / `config.load_spec` (Pydantic v2).
5. **Fibonacci en código** — recurrencia de secuencias y longitud de celda $L_m$.
6. **Materiales** — A homogéneo; B Drude sin pérdidas $\varepsilon=1-\omega_e^2/\omega^2$, $\mu=1-\omega_m^2/\omega^2$.
7. **Cinemática** — $n$, $q=(\omega/c)n_A\sin\theta$, $Q=\sqrt{(\omega/c)^2 n^2-q^2}$, $\chi$.
8. **Matriz de una capa** — deducida y contrastada con las eqs. (7)–(9) del paper.
9. **Semitraza** — producto directo vs recurrencia de Fibonacci (equivalentes).
10. **Barrido de dispersión** — grilla en $\nu$; criterio de banda $\lvert R\rvert\le 1+\varepsilon_{\mathrm{tol}}$.
11. **Detección de modos** — intervalos conexos de frecuencia permitida cerca del plasmón.
12. **Dibujo de $k(\nu)$** — cortes con NaN en saltos de $k$ para evitar “manchas”.
13. **Pipeline por figura** — qué ventanas, órdenes $m$ y ángulos usa cada script.
14. **Regularizaciones** — polos Drude, clipping cerca de $\lvert R\rvert=1$, tolerancias.
15. **Tests** — qué valida cada archivo en `tests/`.
16. **Cómo reproducir** — comandos de scripts y Makefile.
17. **Decisiones `[IMPL]`** — $c$ SI, restauración de $\omega/c$ en $q$, $m=2..7$ en Fig. 6.

**Para quién:** quien vaya a modificar algoritmos, resoluciones o criterios numéricos.

---

## 3. Arquitectura, flujo y ejecución (`arquitectura`)

Manual operativo del repositorio: árbol de directorios, dependencias y problemas frecuentes.

**Temas principales**

1. **Qué es el repositorio** — reproducción científica, no un demonio ni una GUI.
2. **Árbol del proyecto** — `src/`, `scripts/{tmm,pwe,comparison}/`, `configs/`, `results/{tmm,pwe,comparison}/`, `tests/`, `docs/`.
3. **Módulos y responsabilidades** — mapa de `fibonacci_photonics/{physics,solvers,config,analysis,benchmark,reproduction,studies,io,viz}`.
4. **Objetos que cruzan fronteras** — `SuperlatticeSpec`, resultados de scan, metadatos.
5. **Capa de configuración** — un YAML por figura validado con Pydantic (`extra="forbid"`).
6. **Capa de scripts** — `scripts/tmm/figure_0N.py`, `scripts/pwe/figure_0N.py`, `scripts/comparison/` y el runner único `scripts/run_figure.py`.
7. **Flujo de una corrida** — de config a PNG/SVG/PDF + `summary.json`.
8. **Tests** — `pytest` en `unit/`, `analytic/` y `regression/`; algunos casos skipped por $k>m$.
9. **Dependencias** — runtime Python; LaTeX solo para PDF; notebooks opcionales.
10. **Cómo ejecutar** — venv, tests, figuras, docs, `make all`.
11. **Cambiar un parámetro** — editar YAML y re-ejecutar el script de esa figura.
12. **Problemas frecuentes** — path del venv, resolución insuficiente, polo Drude.
13. **Dónde está cada documento** — índice de PDFs y Markdown.

**Para quién:** quien clona el repo en un computador nuevo y necesita ponerlo en marcha.

---

## 4. Documentación del software (`software_documentation`)

Referencia más corta orientada a la API del paquete `fibonacci_photonics` (v2.0).

**Temas principales**

- Arquitectura en capas (física → TMM → barridos → I/O).
- Instalación (`pip install -r requirements.txt` / editable).
- Uso mínimo desde Python (cargar spec, escanear dispersión).
- API: materiales, Fibonacci, transferencia, dispersión, modos, anchos de banda.
- Modelo físico resumido y notas de TMM.
- Tests y reproducción de figuras.
- Solución de problemas y posibles extensiones (pérdidas, TM, perfiles de campo).

**Para quién:** desarrolladores que integren o extiendan la biblioteca.

---

## 5. Paper de reimplementación (`paper_reimplementation`)

Artículo propio que documenta la reproducción científica: modelo, método, figuras
generadas y comparación con el PRB 2010.

**Estructura**

1. Introducción y objetivos de reproducibilidad.
2. Fundamentos: gap $\langle n\rangle=0$, metamateriales, polaritones, Fibonacci.
3. Modelo: geometría, $S_m$, Drude.
4. Formulación electromagnética y TMM ($\cos(kL_m)=R_m$).
5. Implementación computacional.
6. **Reproducción de resultados** (Figs. 1–6) — ver también [RESULTADOS.md](RESULTADOS.md).
7. Tabla de concordancia con el original.
8. Análisis físico (acoplamiento fuerte vs casi-plasmónico).
9. Convergencia numérica.
10. Limitaciones y discrepancias documentadas.
11. Conclusiones y apéndices (unidades, cómo reproducir).

**Para quién:** evaluación académica o lectura tipo paper de lo hecho.

---

## 5b. Método de ondas planas (`metodo_ondas_planas`)

Paper que desarrolla el PWE para la superred siguiendo el capítulo 1D de
Sukhoivanov y Guryev (`referencias/metodo_ondas_planas/PWE-Method-1D.pdf`, cuaderno
`PWE-1d.ipynb`).

1. Planteamiento, celda unidad $S_m$, red recíproca y zona de Brillouin.
2. Coeficientes de Fourier analíticos (generalización de la ec. 4.38) y síntesis.
3. Ecuación maestra, hermiticidad, campo E vs H y propagación fuera del eje.
4. Validación con el cristal del libro (bandas, campo, Brewster).
5. Por qué falla la forma $\omega(k)$ con Drude (contaminación espectral).
6. Forma $k(\omega)$: problema cuadrático en $k$ con la regla inversa de Li.
7. Figs. 1–6 con el PWE, convergencia, costo y límite de validez cerca de $\nu_m$.

## 5c. Anexo de comparación (`comparacion_tmm_pwe`)

Superposiciones TMM/PWE, error punto a punto de $R_m$, conteos de subbandas,
errores de bordes y anchos (Figs. 5–6), tiempos y tablas generadas desde
`results/comparison/`.

---

## 6. Auditoría matemática (`fase_A_auditoria_matematica.md`)

Documento **previo** a escribir código. Etiqueta cada afirmación:

| Etiqueta | Significado |
|----------|-------------|
| `[PAPER]` | Dicho explícitamente en el PRB 2010 |
| `[DEDUC]` | Consecuencia unívoca de las ecuaciones |
| `[REF16]` | Tomado del EPL 2009 (mismo grupo) |
| `[AMB]` | Ambiguo o incompleto en el paper |
| `[IMPL]` | Decisión de implementación propia |

Incluye geometría, Fibonacci, materiales, Helmholtz, TMM, $q$, parámetros de cada
figura y lista de ambigüedades (entre ellas la omisión de $\omega/c$ en $q$).

---

## 7. Informe de validación (`INFORME_VALIDACION.md`)

Tabla condensada figura por figura: reproducida / concordancia / diferencia / causa.
Lista entorno (versiones de librerías), tolerancia de banda, resolución usada y
parámetros del paper.

---

## Relación recomendada de lectura

1. Este anexo (mapa) → **guía conceptual** (física).
2. **Auditoría matemática** (qué es del paper vs inferido).
3. **Solución numérica** + código en `src/`.
4. **Arquitectura** (cómo correr).
5. **RESULTADOS** / paper de reimplementación (qué se obtuvo).
