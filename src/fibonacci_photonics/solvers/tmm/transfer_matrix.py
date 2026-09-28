"""Matriz de transferencia de una celda de Fibonacci.

Matriz de una capa homogénea de espesor d (reconstruida de Ψ; reproduce las
ecs. 7–9 del paper):

.. math::

    M = \\begin{pmatrix} \\cos Qd & (\\chi/Q)\\sin Qd \\\\
                         -(Q/\\chi)\\sin Qd & \\cos Qd \\end{pmatrix},
    \\qquad \\det M = 1,

con χ = μ en TE y χ = ε en TM. Si S_m = X_1 … X_N de izquierda a derecha en z,
la celda es T_m = M_N ⋯ M_1, coherente con T_m = T_{m−2} T_{m−1}.

Como det T = 1, la semitraza R_m = ½ Tr T_m cumple la recurrencia de trazas
de Fibonacci

    R_m = 2 R_{m−1} R_{m−2} − R_{m−3},

que calcula R_m en O(m) operaciones por frecuencia, frente a O(F_m) del
producto de capas. Las entradas son frecuencias angulares ω en rad/s.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from fibonacci_photonics.physics.electromagnetics import (
    Polarization,
    characteristic_parameter,
    in_plane_wavevector,
    longitudinal_wavevector,
    refractive_index,
)
from fibonacci_photonics.physics.fibonacci import sequence, validate_order
from fibonacci_photonics.physics.superlattice import SuperlatticeSpec

SMALL_Q_ABS = 1.0e-14
"""Por debajo de |Q| (1/m) se usa sin(Qd)/Q ≈ d(1 − (Qd)²/6). [IMPL]"""

Stack = NDArray[np.complex128]
"""Pila de matrices 2×2 con forma (2, 2, *forma_de_ω)."""


def layer_matrix(
    Q: ArrayLike,
    thickness: float,
    chi: ArrayLike,
    *,
    small_q: float = SMALL_Q_ABS,
) -> Stack:
    """Matriz de una capa homogénea, apilada sobre las frecuencias.

    Parameters
    ----------
    Q : array_like
        Vector de onda longitudinal en la capa, en 1/m.
    thickness : float
        Espesor d en m.
    chi : array_like
        Parámetro característico χ (μ en TE, ε en TM).
    small_q : float
        Umbral de |Q| para la forma límite de sin(Qd)/Q.

    Returns
    -------
    ndarray
        Forma (2, 2, *broadcast(Q, chi).shape).
    """
    Q_arr = np.asarray(Q, dtype=np.complex128)
    chi_arr = np.asarray(chi, dtype=np.complex128)
    Q_arr, chi_arr = np.broadcast_arrays(Q_arr, chi_arr)
    delta = Q_arr * thickness
    cosine = np.cos(delta)
    sine = np.sin(delta)
    with np.errstate(divide="ignore", invalid="ignore"):
        sinc_q = np.where(np.abs(Q_arr) < small_q, thickness * (1.0 - delta**2 / 6.0), sine / Q_arr)
        m21 = -Q_arr * sine / chi_arr
    stacked = np.stack((np.stack((cosine, chi_arr * sinc_q)), np.stack((m21, cosine))), axis=0)
    return np.asarray(stacked, dtype=np.complex128)


def matmul_2x2(a: Stack, b: Stack) -> Stack:
    """Producto a·b de dos pilas de matrices 2×2, elemento a elemento en frecuencia."""
    out = np.empty(np.broadcast_shapes(a.shape, b.shape), dtype=np.complex128)
    out[0, 0] = a[0, 0] * b[0, 0] + a[0, 1] * b[1, 0]
    out[0, 1] = a[0, 0] * b[0, 1] + a[0, 1] * b[1, 1]
    out[1, 0] = a[1, 0] * b[0, 0] + a[1, 1] * b[1, 0]
    out[1, 1] = a[1, 0] * b[0, 1] + a[1, 1] * b[1, 1]
    return out


def identity_stack(shape_tail: tuple[int, ...]) -> Stack:
    """Identidad 2×2 apilada con forma (2, 2, *shape_tail)."""
    eye = np.zeros((2, 2, *shape_tail), dtype=np.complex128)
    eye[0, 0] = 1.0
    eye[1, 1] = 1.0
    return eye


def _layer_parameters(
    spec: SuperlatticeSpec,
    omega: ArrayLike,
    theta: float,
    polarization: Polarization,
) -> tuple[NDArray[np.complex128], NDArray[np.complex128], NDArray[np.complex128], NDArray[np.complex128]]:
    """(Q_A, Q_B, χ_A, χ_B) a las frecuencias ``omega`` (rad/s)."""
    omega_arr = np.asarray(omega, dtype=np.float64)
    eps_a, mu_a = spec.medium_a.epsilon_mu(omega_arr)
    eps_b, mu_b = spec.medium_b.epsilon_mu(omega_arr)
    n_a = refractive_index(eps_a, mu_a)
    q = in_plane_wavevector(omega_arr, theta, n_a, spec.speed_of_light)
    q_a = longitudinal_wavevector(omega_arr, eps_a, mu_a, q, spec.speed_of_light)
    q_b = longitudinal_wavevector(omega_arr, eps_b, mu_b, q, spec.speed_of_light)
    chi_a = characteristic_parameter(polarization, eps_a, mu_a)
    chi_b = characteristic_parameter(polarization, eps_b, mu_b)
    return q_a, q_b, chi_a, chi_b


def layer_matrices_ab(
    spec: SuperlatticeSpec,
    omega: ArrayLike,
    theta: float,
    polarization: Polarization,
) -> tuple[Stack, Stack]:
    """Matrices M_A y M_B de una capa de cada material."""
    q_a, q_b, chi_a, chi_b = _layer_parameters(spec, omega, theta, polarization)
    return layer_matrix(q_a, spec.thickness_a, chi_a), layer_matrix(q_b, spec.thickness_b, chi_b)


def cell_transfer_matrix(
    spec: SuperlatticeSpec,
    m: int,
    omega: ArrayLike,
    theta: float,
    polarization: Polarization,
) -> Stack:
    """T_m = M_N ⋯ M_1 por producto directo de las F_m matrices de capa (referencia)."""
    omega_arr = np.asarray(omega, dtype=np.float64)
    m_a, m_b = layer_matrices_ab(spec, omega_arr, theta, polarization)
    matrices = {"A": m_a, "B": m_b}
    t = identity_stack(omega_arr.shape)
    for symbol in sequence(m):
        t = matmul_2x2(matrices[symbol], t)
    return t


def cell_matrix_by_recurrence_blocks(
    spec: SuperlatticeSpec,
    m: int,
    omega: ArrayLike,
    theta: float,
    polarization: Polarization,
) -> Stack:
    """T_m = T_{m−2} T_{m−1} con T_0 = M_B y T_1 = M_A: O(m) productos 2×2."""
    m = validate_order(m)
    m_a, m_b = layer_matrices_ab(spec, omega, theta, polarization)
    if m == 0:
        return m_b
    t_prev2, t_prev1 = m_b, m_a
    for _ in range(2, m + 1):
        t_prev2, t_prev1 = t_prev1, matmul_2x2(t_prev2, t_prev1)
    return t_prev1


def semitrace(matrix: Stack) -> NDArray[np.complex128]:
    """R = ½ Tr T."""
    return 0.5 * (matrix[0, 0] + matrix[1, 1])


def determinant(matrix: Stack) -> NDArray[np.complex128]:
    """det T; vale 1 para toda matriz de transferencia sin pérdidas."""
    return matrix[0, 0] * matrix[1, 1] - matrix[0, 1] * matrix[1, 0]


def analytic_r0_r1_r2(
    spec: SuperlatticeSpec,
    omega: ArrayLike,
    theta: float,
    polarization: Polarization,
) -> tuple[NDArray[np.complex128], NDArray[np.complex128], NDArray[np.complex128]]:
    """Semitrazas iniciales, ecs. (7)–(9) del paper (en TM, χ = ε en la (9)).

    .. math::

        R_0 = \\cos Q_B b, \\quad R_1 = \\cos Q_A a, \\quad
        R_2 = \\cos Q_A a \\cos Q_B b - \\tfrac12\\Big(\\frac{Q_A\\chi_B}{Q_B\\chi_A}
              + \\frac{Q_B\\chi_A}{Q_A\\chi_B}\\Big)\\sin Q_A a \\sin Q_B b.
    """
    q_a, q_b, chi_a, chi_b = _layer_parameters(spec, omega, theta, polarization)
    a = spec.thickness_a
    b = spec.thickness_b
    r0 = np.cos(q_b * b)
    r1 = np.cos(q_a * a)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = 0.5 * (q_a * chi_b / (q_b * chi_a) + q_b * chi_a / (q_a * chi_b))
    r2 = np.cos(q_a * a) * np.cos(q_b * b) - ratio * np.sin(q_a * a) * np.sin(q_b * b)
    return r0, r1, r2


def semitrace_by_recurrence(
    spec: SuperlatticeSpec,
    m: int,
    omega: ArrayLike,
    theta: float,
    polarization: Polarization,
) -> NDArray[np.complex128]:
    """R_m a partir de R_0, R_1, R_2 y R_m = 2 R_{m−1} R_{m−2} − R_{m−3}.

    Cerca de ν_m, |R| crece sin cota (gap evanescente) y para m grande puede
    desbordar a ±inf; se deja propagar sin advertencias, porque esas
    frecuencias son gap de todos modos.
    """
    m = validate_order(m)
    r0, r1, r2 = analytic_r0_r1_r2(spec, omega, theta, polarization)
    if m == 0:
        return r0
    if m == 1:
        return r1
    prev3, prev2, prev1 = r0, r1, r2
    with np.errstate(over="ignore", invalid="ignore"):
        for _ in range(3, m + 1):
            prev3, prev2, prev1 = prev2, prev1, 2.0 * prev1 * prev2 - prev3
    return prev1


def semitrace_by_product(
    spec: SuperlatticeSpec,
    m: int,
    omega: ArrayLike,
    theta: float,
    polarization: Polarization,
) -> NDArray[np.complex128]:
    """R_m = ½ Tr(M_N ⋯ M_1): referencia independiente de la recurrencia, O(F_m)."""
    return semitrace(cell_transfer_matrix(spec, m, omega, theta, polarization))
