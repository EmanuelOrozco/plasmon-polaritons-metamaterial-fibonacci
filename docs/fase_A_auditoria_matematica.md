# Fase A — Auditoría matemática del artículo

**Artículo:** E. Reyes-Gómez, N. Raigoza, S. B. Cavalcanti, C. A. A. de Carvalho y L. E. Oliveira,
*Plasmon polaritons in photonic metamaterial Fibonacci superlattices*,
Phys. Rev. B **81**, 153101 (2010). DOI: 10.1103/PhysRevB.81.153101

**Auditor:** Emanuel Orozco Gallego
**Fecha:** 31 de agosto de 2026
**Fuente primaria:** `Plasmon_polaritons_in_photonic_metamaterial_Fibonacci_superlattices.pdf` (4 páginas, Brief Report)

**Estado:** auditoría previa a cualquier implementación. No se ha escrito código científico.

Convención de etiquetas usada en todo el documento:

| Etiqueta | Significado |
| --- | --- |
| `[PAPER]` | Afirmación explícita en el artículo de 2010 |
| `[DEDUC]` | Consecuencia matemática unívoca de las ecuaciones del paper |
| `[REF16]` | Información tomada del artículo hermano EPL **88**, 24002 (2009), citado como Ref. 16 |
| `[AMB]` | Ambiguo, incompleto o no especificado en el paper |
| `[IMPL]` | Decisión de implementación propuesta, **no** es un parámetro original |

---

## 1. Identificación del artículo

- Tipo: Brief Report de 4 páginas.
- Polarización principal de las figuras: TE `[PAPER]`.
- TM: se afirma que las ecuaciones (6)–(9) valen sustituyendo \(\mu_A,\mu_B\) por \(\varepsilon_A,\varepsilon_B\) en (9); no hay figuras TM `[PAPER]`.
- Método: formalismo de matriz de transferencia a incidencia oblicua `[PAPER]`.
- PACS: 42.70.Qs, 41.20.Jb, 42.70.Gi, 78.20.Bh.

---

## 2. Geometría y secuencia de Fibonacci

### 2.1 Bloques constitutivos `[PAPER]`

El sistema es periódico a lo largo de \(z\). La celda elemental es la generación \(m\)-ésima de Fibonacci, construida con dos bloques:

| Bloque | Espesor | Permitividad | Permeabilidad | Papel |
| --- | --- | --- | --- | --- |
| A | \(a\) | \(\varepsilon_A\) | \(\mu_A\) | índice positivo (aire) |
| B | \(b\) | \(\varepsilon_B(\omega)\) | \(\mu_B(\omega)\) | metamaterial (Drude) |

### 2.2 Recurrencia, inflación y relación general `[PAPER]`

\[
S_m(A,B)=S_{m-1}(A,B)\,S_{m-2}(A,B),\qquad m\ge 2,
\]

con

\[
S_0(A,B)=B,\qquad S_1(A,B)=A.
\]

Ley de inflación: \(B\to A\), \(A\to AB\), o equivalentemente \(S_m(A,B)=S_{m-1}(AB,A)\).

Relación general `[PAPER]`:

\[
S_m(A,B)=S_{m-k}\bigl[S_{k+1}(A,B),\,S_k(A,B)\bigr],
\qquad m\ge 2,\quad 0\le k\le m.
\]

### 2.3 Números de Fibonacci usados en el paper `[PAPER]`

\[
F_m=F_{m-1}+F_{m-2},\qquad F_0=F_1=1.
\]

Esta indexación **no** es la de \(F_0=0,\,F_1=1\). Equivale a la sucesión clásica desplazada: \(F_m^{\text{(paper)}}=F_{m+1}^{\text{(clásica)}}\).

Razón áurea `[PAPER]`:

\[
\tau_m=\frac{F_{m+1}}{F_m}\xrightarrow{m\to\infty}\tau=\frac{1+\sqrt{5}}{2}.
\]

(El texto extraído por OCR llegó a leer \(\sqrt{2}\); la página impresa dice \(\sqrt{5}\).)

Longitud de la celda `[PAPER]`:

\[
L_m=F_{m-1}\,a+F_{m-2}\,b.
\]

Número de capas B (metamaterial) `[PAPER]`:

\[
N_B(S_m)=F_{m-2}.
\]

Número total de capas `[DEDUC]`: \(N(S_m)=F_m\).
Número de capas A `[DEDUC]`: \(N_A(S_m)=F_{m-1}\).

Extensión hacia índices negativos `[DEDUC]`, para que las fórmulas cubran \(m=0,1\):

\[
F_{-1}=0,\qquad F_{-2}=1.
\]

### 2.4 Secuencias explícitas `[DEDUC]`

| \(m\) | \(S_m\) | \(N_A=F_{m-1}\) | \(N_B=F_{m-2}\) | \(N=F_m\) | \(L_m\) si \(a=b=d\) |
| --- | --- | --- | --- | --- | --- |
| 0 | B | 0 | 1 | 1 | \(d\) |
| 1 | A | 1 | 0 | 1 | \(d\) |
| 2 | AB | 1 | 1 | 2 | \(2d\) |
| 3 | ABA | 2 | 1 | 3 | \(3d\) |
| 4 | ABAAB | 3 | 2 | 5 | \(5d\) |
| 5 | ABAABABA | 5 | 3 | 8 | \(8d\) |
| 6 | ABAABABAABAAB | 8 | 5 | 13 | \(13d\) |
| 7 | \(S_6S_5\) | 13 | 8 | 21 | \(21d\) |
| 8 | \(S_7S_6\) | 21 | 13 | 34 | \(34d\) |

Con \(d=12\,\mathrm{mm}\): \(L_3=36\,\mathrm{mm}\), \(L_4=60\,\mathrm{mm}\), \(L_5=96\,\mathrm{mm}\), \(L_6=156\,\mathrm{mm}\), \(L_7=252\,\mathrm{mm}\), \(L_8=408\,\mathrm{mm}\).

Número de modos plasmon-polaritón `[PAPER]`: \(N_{\text{modos}}=F_{m-2}\).

| \(m\) | \(F_{m-2}\) | Esperado |
| --- | --- | --- |
| 3 | 1 | 1 |
| 4 | 2 | 2 |
| 5 | 3 | 3 |
| 6 | 5 | 5 |
| 7 | 8 | 8 |
| 8 | 13 | 13 |

Verificación de inflación `[DEDUC]`:

- \(S_1=A\)
- \(S_2=AB\)
- \(S_3=ABA\)
- \(S_4=ABAAB\)
- \(S_5=ABAABABA\)

Coincide con la concatenación \(S_m=S_{m-1}S_{m-2}\).

---

## 3. Modelo de materiales

### 3.1 Medio A `[PAPER]`

\[
\varepsilon_A=1,\qquad \mu_A=1
\]

(aire). Espesor en todas las figuras: \(a=12\,\mathrm{mm}\).

### 3.2 Medio B — Drude sin pérdidas `[PAPER]`

Ecuación (10):

\[
\varepsilon_B(\omega)=\varepsilon_0-\frac{\omega_e^2}{\omega^2},\qquad
\mu_B(\omega)=\mu_0-\frac{\omega_m^2}{\omega^2}.
\]

Valores numéricos usados `[PAPER]`: \(\varepsilon_0=1\), \(\mu_0=1\).

No hay término de amortiguamiento \(\gamma\). El modelo es **sin pérdidas** `[PAPER]` (y `[REF16]` lo confirma para el mismo grupo).

### 3.3 Frecuencias de plasmón `[PAPER]`

\[
\nu_e=\frac{\omega_e}{2\pi\sqrt{\varepsilon_0}},\qquad
\nu_m=\frac{\omega_m}{2\pi\sqrt{\mu_0}}.
\]

Con \(\varepsilon_0=\mu_0=1\):

\[
\nu_e=\frac{\omega_e}{2\pi},\qquad \nu_m=\frac{\omega_m}{2\pi}.
\]

Ceros de las funciones de respuesta `[DEDUC]`:

- \(\varepsilon_B(\omega)=0\) \(\Leftrightarrow\) \(\omega=\omega_e/\sqrt{\varepsilon_0}\) \(\Leftrightarrow\) \(\nu=\nu_e\)
- \(\mu_B(\omega)=0\) \(\Leftrightarrow\) \(\omega=\omega_m/\sqrt{\mu_0}\) \(\Leftrightarrow\) \(\nu=\nu_m\)

Índice negativo (LHM) cuando \(\varepsilon_B<0\) y \(\mu_B<0\), es decir `[DEDUC]`:

\[
\nu<\min(\nu_e,\nu_m)\qquad\text{(con }\varepsilon_0=\mu_0=1\text{)}.
\]

### 3.4 Unidades de frecuencia — no mezclar \(f\) y \(\omega\) `[PAPER]`

El paper es explícito:

- \(\omega\): frecuencia angular.
- \(\nu=\omega/(2\pi)\): frecuencia cíclica.
- Ejes verticales de las figuras: **frequency (GHz)** = \(\nu\), no \(\omega\).
- Los parámetros de plasma se dan como \(\omega_e/2\pi\) y \(\omega_m/2\pi\) en GHz.

Por tanto `[DEDUC]`:

\[
\omega=2\pi\nu,\qquad
\frac{\omega_e}{2\pi}=3\,\mathrm{GHz}\ \Rightarrow\ \omega_e=2\pi\times 3\times 10^9\,\mathrm{rad/s}.
\]

**No** debe usarse \(\omega_e=3\,\mathrm{GHz}\) sin el \(2\pi\).

### 3.5 Índice de refracción `[PAPER]`

\[
n(z)=\sqrt{\mu(z)}\,\sqrt{\varepsilon(z)}.
\]

Dos raíces, no \(\sqrt{\mu\varepsilon}\) como un solo radical. Esto fija el signo en medios LHM: si ambas raíces principales de números negativos se toman como \(+i\sqrt{|\cdot|}\), entonces \(\sqrt{\mu}\sqrt{\varepsilon}=-\sqrt{|\mu\varepsilon|}\) cuando \(\varepsilon<0\) y \(\mu<0\) `[DEDUC]`.

Para el vector de onda longitudinal solo entra \(n^2=\mu\varepsilon\), de modo que el signo de \(n\) no cambia la relación de dispersión \(\cos(kL_m)=R_m\) `[DEDUC]`.

### 3.6 Promedios y gap \(\langle n\rangle_m=0\) `[PAPER]`

\[
\langle\varepsilon\rangle_m=\frac{F_{m-1}\varepsilon_A a+F_{m-2}\varepsilon_B b}{L_m},
\qquad
\langle\mu\rangle_m=\frac{F_{m-1}\mu_A a+F_{m-2}\mu_B b}{L_m}.
\]

En la Fig. 1 (\(\omega_e=\omega_m\)) el paper afirma que a \(\theta=0\) los gaps no Bragg \(\langle n\rangle_m=0\) son nulos porque \(\langle\varepsilon\rangle_m=\langle\mu\rangle_m=0\). Eso ocurre **a una frecuencia concreta**, no idénticamente `[DEDUC]`. Con \(a=b\) y \(\varepsilon_A=\mu_A=\varepsilon_0=\mu_0=1\), \(\varepsilon_B=\mu_B=1-\omega_p^2/\omega^2\):

\[
\langle\varepsilon\rangle_m=0 \quad\Leftrightarrow\quad
\omega=\omega_p\sqrt{\frac{F_{m-2}}{F_m}}.
\]

| \(m\) | \(F_{m-2}/F_m\) | \(\nu_{\langle n\rangle=0}\) si \(\nu_p=3\,\mathrm{GHz}\) |
| --- | --- | --- |
| 3 | \(1/3\) | \(1.732\,\mathrm{GHz}\) |
| 4 | \(2/5\) | \(1.897\,\mathrm{GHz}\) |
| 5 | \(3/8\) | \(1.837\,\mathrm{GHz}\) |
| 6 | \(5/13\) | \(1.861\,\mathrm{GHz}\) |

En la Fig. 3, \(\nu_m=1\,\mathrm{GHz}\) se coloca **dentro** del gap \(\langle n\rangle_m=0\) a incidencia normal `[PAPER]`.

---

## 4. Electromagnetismo

### 4.1 Polarización TE `[PAPER]` + campo `[REF16]`

El paper 2010 escribe la ecuación (1) para la amplitud \(E(z)\) de modos TE (campo eléctrico paralelo a las interfaces). La dependencia transversal y temporal no se escribe. El artículo hermano `[REF16]`, mismo grupo y misma ecuación de Helmholtz, usa

\[
\mathbf{E}(\mathbf{r},t)=E(z)\,e^{i(qx-\omega t)}\,\mathbf{e}_y.
\]

Convención temporal: \(e^{-i\omega t}\) `[REF16]`.

Ecuación (1) `[PAPER]`:

\[
\frac{d}{dz}\left[\frac{1}{\mu(z)}\frac{dE}{dz}\right]
=-\varepsilon(z)\left[\frac{\omega^2}{c^2}-\frac{q^2}{n^2(z)}\right]E(z).
\]

Equivalente a `[DEDUC]`:

\[
\frac{d}{dz}\left[\frac{1}{\mu}\frac{dE}{dz}\right]
+\left(\varepsilon\frac{\omega^2}{c^2}-\frac{q^2}{\mu}\right)E=0.
\]

### 4.2 Polarización TM `[PAPER]` + `[REF16]`

El paper 2010 no escribe la ecuación para \(H\). Indica sustituir \(\mu\to\varepsilon\) en (9). `[REF16]` da

\[
\frac{d}{dz}\left[\frac{1}{\varepsilon(z)}\frac{dH}{dz}\right]
=-\mu(z)\left[\frac{\omega^2}{c^2}-\frac{q^2}{n^2(z)}\right]H(z).
\]

TE \(\leftrightarrow\) plasmón **magnético** (\(\nu_m\)); TM \(\leftrightarrow\) plasmón **eléctrico** (\(\nu_e\)) `[PAPER]`.

### 4.3 Vector de estado y continuidad `[PAPER]`

\[
\Psi(z)=\begin{pmatrix} E(z) \\ \dfrac{1}{\mu(z)}\dfrac{dE}{dz} \end{pmatrix}
\qquad\text{es continuo en las interfaces.}
\]

Para TM, el análogo es \(\bigl(H,\;\varepsilon^{-1}dH/dz\bigr)\) `[DEDUC]` / `[REF16]`.

### 4.4 Componente transversal \(q\) — ambigüedad resuelta

El paper escribe `[PAPER]`:

\[
q=n_A\sin\theta_A=n_B\sin\theta_B,\qquad \theta_A=\theta.
\]

Esto es **dimensionalmente inconsistente** con

\[
Q_{A,B}=\sqrt{\frac{\omega^2}{c^2}n_{A,B}^2-q^2},
\]

porque \(Q\) tiene unidades de inverso de longitud y \(n\sin\theta\) es adimensional `[AMB]`.

Interpretación dimensionally consistente `[DEDUC]`, confirmada por el mismo grupo `[REF16]`:

\[
q=\frac{\omega}{c}\,n_A\sin\theta=\frac{\omega}{c}\,n_B\sin\theta_B.
\]

La notación del paper 2010 omite el factor \(\omega/c\) (ley de Snell escrita para índices, no para el vector de onda). **Se implementará la forma con \(\omega/c\)**. Usar \(q=n\sin\theta\) sin \(\omega/c\) haría \(Q\) complejo de forma espuria y destruiría las escalas en GHz/mm.

Luego `[DEDUC]`:

\[
Q_{A,B}=\frac{\omega}{c}\sqrt{n_{A,B}^2-n_A^2\sin^2\theta}.
\]

Cuando el radicando es negativo, \(Q\) es imaginario (ondas evanescentes / hiperbólicas). `[REF16]` escribe explícitamente las formas \(\cos/\cosh\) según el signo de \(n_B^2-n_A^2\sin^2\theta\). Una implementación unificada con \(Q\) complejo (continuación analítica de senos y cosenos) reproduce ambas.

### 4.5 Velocidad de la luz `[AMB]`

El paper no da el valor numérico de \(c\). Con \(a=12\,\mathrm{mm}\) y \(\nu\) en GHz,

\[
\frac{\omega}{c}=\frac{2\pi\nu}{c}.
\]

Propuesta `[IMPL]`: \(c=299\,792\,458\,\mathrm{m/s}\). Alternativa común: \(c=3\times 10^8\,\mathrm{m/s}\) (diferencia relativa \(\approx 0.07\%\), inferior a la resolución visual de las figuras).

---

## 5. Método de matriz de transferencia

### 5.1 Definiciones `[PAPER]`

Operador de traslación sobre la celda \(S_m\):

\[
\hat T_{L_m}\Psi(0)=\Psi(L_m)=e^{ik L_m}\Psi(0).
\]

Matriz de transferencia \(2\times 2\):

\[
\Psi(L_m)=T_m\Psi(0).
\]

Autovalores: \(e^{\pm ik L_m}\) (matriz unimodular). Relación de dispersión, ecuación (6):

\[
\cos(k L_m)=R_m,\qquad R_m=\frac12\mathrm{Tr}(T_m).
\]

Bandas permitidas: \(|R_m|\le 1\). Gaps: \(|R_m|>1\) `[DEDUC]`.

Eje horizontal de las figuras: \(k L_m/\pi\in[-1,1]\) (primera zona de Brillouin de la superred de periodo \(L_m\)) `[PAPER]`.

Como \(\cos\) es par, \(\nu(k)=\nu(-k)\) `[DEDUC]`.

### 5.2 Orden de multiplicación `[PAPER]`

Concatenación: \(S_m=S_{m-1}S_{m-2}\) (primero \(S_{m-1}\), luego \(S_{m-2}\), al crecer \(z\)).

Matrices: \(T_m=T_{m-2}T_{m-1}\).

Esto es coherente con \(\Psi(L)=M_{\text{última}}\cdots M_{\text{primera}}\Psi(0)\) `[DEDUC]`. Ejemplo: \(S_2=AB\) \(\Rightarrow\) \(T_2=T_0 T_1=M_B M_A\).

La semitraza no depende del orden porque \(\mathrm{Tr}(AB)=\mathrm{Tr}(BA)\) `[DEDUC]`.

### 5.3 Recurrencia de la semitraza `[PAPER]`

\[
R_m=2 R_{m-1} R_{m-2}-R_{m-3},\qquad m>2.
\]

Basta conocer \(R_0,R_1,R_2\). Es el mapa de Kohmoto (Ref. 17).

### 5.4 Semitrazas iniciales `[PAPER]`

\[
R_0=\cos(Q_B b),
\]

\[
R_1=\cos(Q_A a),
\]

\[
R_2=\cos(Q_A a)\cos(Q_B b)
-\frac12\left[\frac{Q_A\mu_B}{Q_B\mu_A}+\frac{Q_B\mu_A}{Q_A\mu_B}\right]
\sin(Q_A a)\sin(Q_B b).
\]

\[
Q_{A,B}=\sqrt{\frac{\omega^2}{c^2}n_{A,B}^2-q^2},\qquad
n_{A,B}=\sqrt{\mu_{A,B}}\,\sqrt{\varepsilon_{A,B}}.
\]

TM: reemplazar \(\mu_A,\mu_B\) por \(\varepsilon_A,\varepsilon_B\) **solo en (9)** `[PAPER]`. \(R_0\) y \(R_1\) no cambian de forma.

### 5.5 Matriz de una capa — no está escrita; se reconstruye de forma unívoca

El paper **no** da \(M_j\) explícita `[AMB]`. Se deduce de la continuidad de \(\Psi\) en un medio homogéneo `[DEDUC]`.

En una capa homogénea de espesor \(d\), permeabilidad \(\mu\) y \(Q=\sqrt{(\omega/c)^2 n^2-q^2}\),

\[
E(z)=A\cos(Qz)+B\sin(Qz)
\]

respecto a un origen local. Imponiendo \(\Psi(z+d)=M\Psi(z)\) se obtiene la matriz unimodular (\(\det M=1\)):

\[
M(Q,d,\mu)
=\begin{pmatrix}
\cos(Qd) & (\mu/Q)\sin(Qd)\\
-(Q/\mu)\sin(Qd) & \cos(Qd)
\end{pmatrix}.
\]

Comprobaciones contra el paper `[DEDUC]`:

- \(S_0=B\): \(\frac12\mathrm{Tr}(M_B)=\cos(Q_B b)=R_0\). OK.
- \(S_1=A\): \(\frac12\mathrm{Tr}(M_A)=\cos(Q_A a)=R_1\). OK.
- \(S_2=AB\), \(T_2=M_B M_A\): \(\frac12\mathrm{Tr}(T_2)\) reproduce exactamente (9). OK.

Para TM, \(\mu\to\varepsilon\) y \(E\to H\).

Signo de \(Q\) en LHM `[AMB]`: si \(n<0\), algunos autores toman \(Q<0\). \(R_m\) es **par** en \(Q_A\) y \(Q_B\), así que la dispersión \(\omega(k)\) no cambia `[DEDUC]`. El signo sí afectaría perfiles de campo, que el paper no publica.

Valores singulares `[AMB]` / `[IMPL]`:

- \(Q=0\) (umbral de corte): el límite \(M\to\begin{pmatrix}1& \mu d\\ 0&1\end{pmatrix}\).
- \(\omega\to 0\) o polos de Drude \(\omega\to 0^+\): \(\varepsilon_B,\mu_B\to-\infty\).
- \(\omega=\omega_e\) o \(\omega=\omega_m\): \(\varepsilon_B=0\) o \(\mu_B=0\); \(n=0\); hay que regularizar \(Q/\mu\) y \(\mu/Q\).

---

## 6. Condición de banda y \(k(\omega)\)

Dado un barrido en \(\nu\) (o \(\omega\)) a \(\theta\) fijo `[IMPL]`:

1. Evaluar \(\varepsilon_B(\omega),\mu_B(\omega)\).
2. Calcular \(q\), \(Q_A\), \(Q_B\).
3. Obtener \(R_m\) por producto de matrices **o** por recurrencia.
4. Si \(|R_m|\le 1\), \(k L_m=\pm\arccos(R_m)+2\pi n\). En la primera ZB se grafica \(k L_m/\pi=\pm\arccos(R_m)/\pi\in[-1,1]\).
5. Si \(|R_m|>1\), gap: no hay \(k\) real.

Clipping cerca de \(\pm 1\) `[AMB]`: el paper no lo menciona. `[IMPL]`: tolerancia del tipo \(|R|-1<10^{-12}\) se atribuirá a redondeo y se documentará; no se recortará de forma silenciosa.

---

## 7. Parámetros numéricos por figura

Parámetros globales en **todas** las figuras 1–6 `[PAPER]` salvo plasma:

| Magnitud | Valor | Fuente |
| --- | --- | --- |
| Polarización | TE | todas las figuras |
| \(\varepsilon_A\), \(\mu_A\) | 1, 1 | Fig. 1, texto |
| \(a\) | \(12\,\mathrm{mm}\) | todas |
| \(b\) | \(12\,\mathrm{mm}\) | todas |
| \(\varepsilon_0\), \(\mu_0\) | 1, 1 | texto, “sake of simplicity” |
| Eje \(x\) | \(k L_m/\pi\in[-1,1]\) | Figs. 1–4 |
| Frecuencia graficada | \(\nu=\omega/2\pi\) en GHz | todas |

### Figura 1 — dispersión TE, \(\nu_e=\nu_m=3\,\mathrm{GHz}\)

| Ítem | Valor |
| --- | --- |
| \(m\) | 3 (líneas discontinuas) y 4 (líneas continuas) |
| \(\omega_e/2\pi\) | \(3\,\mathrm{GHz}\) |
| \(\omega_m/2\pi\) | \(3\,\mathrm{GHz}\) |
| \(\theta\) | (a) \(0\); (b) \(\pi/12\); (c) \(\pi/6\); (d) \(\pi/3\) |
| Eje \(y\) | \(0\)–\(5\,\mathrm{GHz}\) |
| Observables | a \(\theta=0\) no hay polaritones longitudinales ni gap \(\langle n\rangle=0\) abierto; a \(\theta\neq 0\) aparece 1 modo (\(m=3\)) y 2 modos (\(m=4\)) |

### Figura 2 — vecindad del plasmón magnético, mismos parámetros que Fig. 1

| Ítem | Valor |
| --- | --- |
| \(m\) | (a) 3, (b) 4, (c) 5, (d) 6 |
| \(\theta\) | \(\pi/3\) |
| \(\omega_e/2\pi=\omega_m/2\pi\) | \(3\,\mathrm{GHz}\) |
| Eje \(y\) | \(2\)–\(4\,\mathrm{GHz}\) |
| Línea de referencia | \(\nu_m=3\,\mathrm{GHz}\) (discontinua) |
| Test cuantitativo | número de subbandas \(=F_{m-2}=1,2,3,5\) |

### Figura 3 — \(\nu_m\) dentro del gap \(\langle n\rangle=0\)

| Ítem | Valor |
| --- | --- |
| Como Fig. 1 en geometría y \(\theta\) | \(m=3\) (discontinuo) y \(m=4\) (continuo) |
| \(\omega_e/2\pi\) | \(3\,\mathrm{GHz}\) |
| \(\omega_m/2\pi\) | \(1\,\mathrm{GHz}\) |
| Eje \(y\) | \(0\)–\(5\,\mathrm{GHz}\) |
| Física | a \(\theta=0\), \(\nu_m\) cae en el gap no Bragg; a \(\theta\neq 0\) aparecen bandas casi planas de carácter plasmónico |

### Figura 4 — zoom de los polaritones de la Fig. 3

| Panel | \(m\) | \(\theta\) | Eje \(y\) (GHz) | N. de curvas |
| --- | --- | --- | --- | --- |
| (a) | 3 | \(\pi/12\) | \(0.994\)–\(1.002\) | 1 |
| (b) | 3 | \(\pi/3\) | \(0.94\)–\(1.02\) | 1 |
| (c) | 4 | \(\pi/12\) | \(0.994\)–\(1.002\) | 2 |
| (d) | 4 | \(\pi/3\) | \(0.94\)–\(1.02\) | 2 |

Parámetros: \(a=b=12\,\mathrm{mm}\), \(\omega_e/2\pi=3\,\mathrm{GHz}\), \(\omega_m/2\pi=1\,\mathrm{GHz}\).

Anchos de banda leídos visualmente (orden de magnitud, **no** datos digitales del paper):

- \(m=3\), \(\theta=\pi/12\): \(\Delta\nu\sim 5\times 10^{-3}\,\mathrm{GHz}\)
- \(m=3\), \(\theta=\pi/3\): \(\Delta\nu\sim 4\times 10^{-2}\,\mathrm{GHz}\)

El ancho **crece** con \(\theta\). Las bandas son cóncavas hacia arriba (mínimo en \(k=0\)).

### Figura 5 — estructura de bandas vs orden de Fibonacci

| Ítem | Valor |
| --- | --- |
| \(m\) | \(2,3,4,5,6,7,8\) (eje \(x\)) |
| \(\theta\) | (a) \(\pi/12\); (b) \(\pi/3\) |
| \(\omega_e/2\pi\) | \(3\,\mathrm{GHz}\) |
| \(\omega_m/2\pi\) | \(1\,\mathrm{GHz}\) |
| Representación | segmentos verticales en cada \(m\): intervalos de frecuencia permitidos (subbandas) |
| Eje \(y\) (a) | \(\approx 0.994\)–\(1.000\,\mathrm{GHz}\) |
| Eje \(y\) (b) | \(\approx 0.95\)–\(1.00\,\mathrm{GHz}\) |

El eje vertical es **frecuencia** de las subbandas (la longitud de cada segmento es el ancho de banda). El OCR del PDF a veces etiqueta este eje como “bandwidth”; el caption dice “band structures” y los valores numéricos están centrados en \(\nu_m=1\,\mathrm{GHz}\).

A \(m=2\) y \(m=3\) hay un solo segmento (coherente con \(F_{0}=1\) y \(F_{1}=1\)). Luego hay splitting tipo Cantor.

### Figura 6 — anchos de banda vs ángulo

Texto del cuerpo `[PAPER]`: *“Fig. 6 exhibits the various plasmon-polariton bandwidths of photonic superlattices as functions of the incidence angle, for different values of the order \(m\)”*.

Caption en APS (recorte público): *“FIG. 6. (Color online) TE plasmon-polariton bandwidth of photonic superlattices, as …”* (caption truncado en las fuentes disponibles).

| Ítem | Estado |
| --- | --- |
| Cantidad | anchos de banda plasmon-polaritón vs \(\theta\) |
| Polarización | TE |
| Parámetros de plasma y \(a,b\) | **no reiterados** en el párrafo de Fig. 6; el contexto es el de Figs. 3–5 `[AMB]` |
| Inferencia `[IMPL]` | \(a=b=12\,\mathrm{mm}\), \(\omega_e/2\pi=3\,\mathrm{GHz}\), \(\omega_m/2\pi=1\,\mathrm{GHz}\) |
| Valores de \(m\) | no listados explícitamente `[AMB]`; layout de 6 paneles en la columna derecha de la p. 4, compatible con \(m=2,\ldots,7\) |
| Cálculo del ancho | no descrito; debe obtenerse de los datos como \(\nu_{\max}-\nu_{\min}\) de cada subbanda con \(|R_m|\le 1\) `[IMPL]` |

---

## 8. Inventario de ecuaciones (paper 2010)

| Eq. | Contenido |
| --- | --- |
| (1) | Helmholtz TE para \(E(z)\) |
| (2) | Vector \(\Psi(z)\) |
| (3) | Bloch / traslación \(e^{ik L_m}\) |
| (4) | \(\Psi(L_m)=T_m\Psi(0)\) |
| (5) | Problema de autovalores de \(T_m\) |
| (6) | \(\cos(k L_m)=R_m\) |
| (7) | \(R_0=\cos(Q_B b)\) |
| (8) | \(R_1=\cos(Q_A a)\) |
| (9) | \(R_2\) (celda AB) |
| (10) | Drude \(\varepsilon_B\), \(\mu_B\) |
| — | \(T_m=T_{m-2}T_{m-1}\) |
| — | \(R_m=2R_{m-1}R_{m-2}-R_{m-3}\) |
| — | \(Q_{A,B}=\sqrt{(\omega/c)^2 n_{A,B}^2-q^2}\) |
| — | Snell (notación abreviada) |
| — | \(L_m=F_{m-1}a+F_{m-2}b\) |
| — | \(\langle\varepsilon\rangle_m\), \(\langle\mu\rangle_m\) |

No aparecen: matriz \(2\times 2\) de una capa, valor de \(c\), \(\gamma\) de Drude, receta numérica de barrido, TM explícita más allá del reemplazo \(\mu\to\varepsilon\).

---

## 9. Información que el paper no especifica

1. Factor \(\omega/c\) en \(q\) (notación de Snell incompleta).
2. Valor numérico de \(c\).
3. Matriz de transferencia de una sola capa.
4. Rama de las raíces cuadradas de \(\varepsilon\), \(\mu\), \(Q\) (signo en LHM; \(Q\) evanescente).
5. Convención \(e^{\pm i\omega t}\) (sí está en `[REF16]`).
6. Pérdidas / \(\gamma\).
7. Resolución en \(\nu\) y en \(k\); número de puntos.
8. Tratamiento de \(|R_m|\approx 1\) (redondeo).
9. Regularización cuando \(Q=0\), \(\mu=0\) o \(\varepsilon=0\).
10. Caption y lista explícita de \(m\) y \(\theta\) de la Fig. 6.
11. Parámetros de la Fig. 6 (se infieren).
12. Definición operativa de “bandwidth” (¿\(\nu_{\max}-\nu_{\min}\) por subbanda? ¿a \(k=0\)?).
13. Figuras TM (solo afirmación cualitativa).
14. Perfiles de campo de los modos.

Ninguno de estos huecos impide implementar el TMM para Figs. 1–5; la Fig. 6 requiere una inferencia documentada.

---

## 10. Decisiones de implementación propuestas (aún no ejecutadas)

Estas **no** son parámetros del paper. Se usarán en fases posteriores y quedarán etiquetadas como `[IMPL]`:

| Decisión | Valor propuesto | Motivo |
| --- | --- | --- |
| \(q\) | \((\omega/c)n_A\sin\theta\) | consistencia dimensional + `[REF16]` |
| \(c\) | \(299792458\,\mathrm{m/s}\) | SI; comparar con \(3\times 10^8\) en convergencia |
| Matriz de capa | \(M(Q,d,\mu)\) unimodular de la Sec. 5.5 | reproduce (7)–(9) |
| \(Q\) | `sqrt` complejo principal; trigonométricas complejas | cubre propagante y evanescente |
| Signo de \(Q\) | \(\mathrm{Re}\,Q\ge 0\) o \(\mathrm{Im}\,Q\ge 0\) | la dispersión es par en \(Q\) |
| Recurrencia vs producto | ambos; test de igualdad de \(R_m\) | validación |
| Banda | \(\lvert R_m\rvert\le 1+10^{-12}\) | redondeo documentado |
| Ancho de banda | \(\nu_{\max}-\nu_{\min}\) sobre \(k\) real por subbanda conexa | Fig. 6 |
| Fig. 6 params | iguales a Figs. 3–5; \(m=2\ldots 7\) | inferencia; se contrastará con la figura impresa |

---

## 11. Predicciones cuantitativas para validar el código (antes de mirar las figuras)

1. Secuencias y \(N_B=F_{m-2}\), \(L_m=F_{m-1}a+F_{m-2}b\).
2. \(\det T_m=1\) (unimodular, lossless).
3. \(R_m\) por producto = \(R_m\) por recurrencia.
4. \(R_2\) analítico (9) = \(\frac12\mathrm{Tr}(M_B M_A)\).
5. \(\theta=0\): no hay polaritones longitudinales (TE y TM coinciden).
6. Fig. 1, \(\theta=0\): gap \(\langle n\rangle_m=0\) cerrado en \(\nu=\nu_p\sqrt{F_{m-2}/F_m}\).
7. Fig. 2: \(1,2,3,5\) subbandas cerca de \(3\,\mathrm{GHz}\).
8. Fig. 3–4: \(1\) y \(2\) bandas casi planas cerca de \(1\,\mathrm{GHz}\); ancho mayor a \(\theta=\pi/3\) que a \(\pi/12\).
9. Fig. 5: splitting sucesivo; \(m=2\) y \(m=3\) degenerados a un solo intervalo (\(F_{0}=F_{1}=1\)).

---

## 12. Referencias usadas en esta auditoría

1. E. Reyes-Gómez et al., Phys. Rev. B **81**, 153101 (2010). **Fuente primaria.**
2. E. Reyes-Gómez et al., EPL **88**, 24002 (2009) / arXiv:0907.1644. **Ref. 16 del paper.** Usada solo para: (i) \(q=(\omega/c)n_A\sin\theta\); (ii) convención \(e^{i(qx-\omega t)}\); (iii) formas \(\cosh\) de la dispersión periódica AB, que es \(S_2\); (iv) ecuación TM. No se tomaron de allí \(\varepsilon_0=1.21\) ni otros parámetros numéricos distintos de los del PRB 2010.
3. M. Kohmoto, L. P. Kadanoff y C. Tang, Phys. Rev. Lett. **50**, 1870 (1983). **Ref. 17.** Mapa de la traza de Fibonacci.

No se han inventado otras referencias.

---

## 13. Conclusión de la Fase A

El modelo físico, las secuencias, el Drude, la relación \(\cos(kL_m)=R_m\), las semitrazas (7)–(9) y los parámetros de las **Figuras 1–5** están suficientemente especificados para una reimplementación TMM fiel.

Las únicas decisiones no unívocas que afectan a las figuras son:

- restaurar \(\omega/c\) en \(q\) (obligatorio; de lo contrario las escalas no tienen sentido);
- el valor de \(c\) (efecto subvisual);
- la receta de la Fig. 6 (caption incompleto; parámetros inferidos).

**La Fase A queda cerrada.** Siguiente paso autorizado por el plan: Fase B (ya cubierta en este documento como extracción estructurada) y, con el visto bueno de esta auditoría, Fase C–D de diseño e implementación.
