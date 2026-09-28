"""Método de matriz de transferencia para una celda de Fibonacci.

Matriz de una capa homogénea (reconstruida de Ψ; reproduce las eqs. 7–9):

    M = [[ cos(Qd),  (χ/Q) sin(Qd) ],
         [ -(Q/χ) sin(Qd), cos(Qd) ]]

χ = μ (TE) o χ = ε (TM). det M = 1.

Orden: si Sm = X1...XN de izquierda a derecha en z,
T = M_N ... M_1, coherente con Tm = T(m-2) T(m-1).
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from fibonacci_photonics.core.constants import SMALL_Q_ABS, SPEED_OF_LIGHT
from fibonacci_photonics.core.electromagnetics import (
    Polarization,
    characteristic_parameter,
    in_plane_wavevector,
    longitudinal_wavevector,
    refractive_index,
)
from fibonacci_photonics.core.fibonacci import LAYER_A, sequence
from fibonacci_photonics.core.params import SuperlatticeSpec


def _as_omega(omega: ArrayLike) -> NDArray[np.float64]:
    return np.asarray(omega, dtype=np.float64)


def layer_matrix(
    Q: ArrayLike,
    thickness: float,
    chi: ArrayLike,
    *,
    small_q: float = SMALL_Q_ABS,
) -> NDArray[np.complex128]:
    """Matriz 2×2 (o pila 2×2×Nω) de una capa homogénea."""
    Q_arr = np.asarray(Q, dtype=np.complex128)
    chi_arr = np.asarray(chi, dtype=np.complex128)
    Q_arr, chi_arr = np.broadcast_arrays(Q_arr, chi_arr)
    delta = Q_arr * thickness
    cosine = np.cos(delta)
    sine = np.sin(delta)
    # sin(Qd)/Q → d cuando Q → 0
    sinc_q = np.where(np.abs(Q_arr) < small_q, thickness * (1.0 - delta**2 / 6.0), sine / Q_arr)
    m11 = cosine
    m12 = chi_arr * sinc_q
    m21 = -Q_arr * sine / chi_arr
    m22 = cosine
    stacked = np.stack((np.stack((m11, m12)), np.stack((m21, m22))), axis=0)
    return np.asarray(stacked, dtype=np.complex128)


def _matmul(a: NDArray[np.complex128], b: NDArray[np.complex128]) -> NDArray[np.complex128]:
    return np.einsum("ij...,jk...->ik...", a, b)


def identity_stack(shape_tail: tuple[int, ...]) -> NDArray[np.complex128]:
    eye = np.zeros((2, 2, *shape_tail), dtype=np.complex128)
    eye[0, 0] = 1.0
    eye[1, 1] = 1.0
    return eye


def material_response(
    spec: SuperlatticeSpec,
    omega: ArrayLike,
) -> dict[str, tuple[NDArray[np.complex128], NDArray[np.complex128]]]:
    omega_arr = _as_omega(omega)
    eps_a, mu_a = spec.medium_a.epsilon_mu(omega_arr)
    eps_b, mu_b = spec.medium_b.epsilon_mu(omega_arr)
    return {"A": (eps_a, mu_a), "B": (eps_b, mu_b)}


def wavevectors(
    spec: SuperlatticeSpec,
    omega: ArrayLike,
    theta: float,
    polarization: Polarization,
) -> dict[str, NDArray[np.complex128]]:
    del polarization  # Q no depende de la polarización; χ sí.
    omega_arr = _as_omega(omega)
    responses = material_response(spec, omega_arr)
    eps_a, mu_a = responses["A"]
    n_a = refractive_index(eps_a, mu_a)
    q = in_plane_wavevector(omega_arr, theta, n_a, spec.speed_of_light)
    out: dict[str, NDArray[np.complex128]] = {"q": q}
    for label, (eps, mu) in responses.items():
        out[f"Q_{label}"] = longitudinal_wavevector(
            omega_arr, eps, mu, q, spec.speed_of_light
        )
        out[f"chi_{label}"] = None  # placeholder, filled by caller if needed
    return out


def layer_matrices_ab(
    spec: SuperlatticeSpec,
    omega: ArrayLike,
    theta: float,
    polarization: Polarization,
) -> tuple[NDArray[np.complex128], NDArray[np.complex128]]:
    omega_arr = _as_omega(omega)
    responses = material_response(spec, omega_arr)
    eps_a, mu_a = responses["A"]
    eps_b, mu_b = responses["B"]
    n_a = refractive_index(eps_a, mu_a)
    q = in_plane_wavevector(omega_arr, theta, n_a, spec.speed_of_light)
    q_a = longitudinal_wavevector(omega_arr, eps_a, mu_a, q, spec.speed_of_light)
    q_b = longitudinal_wavevector(omega_arr, eps_b, mu_b, q, spec.speed_of_light)
    chi_a = characteristic_parameter(polarization, eps_a, mu_a)
    chi_b = characteristic_parameter(polarization, eps_b, mu_b)
    m_a = layer_matrix(q_a, spec.thickness_a, chi_a)
    m_b = layer_matrix(q_b, spec.thickness_b, chi_b)
    return m_a, m_b


def cell_transfer_matrix(
    spec: SuperlatticeSpec,
    m: int,
    omega: ArrayLike,
    theta: float,
    polarization: Polarization,
) -> NDArray[np.complex128]:
    """Tm por producto directo de matrices de capa."""
    omega_arr = _as_omega(omega)
    m_a, m_b = layer_matrices_ab(spec, omega_arr, theta, polarization)
    matrices = {"A": m_a, "B": m_b}
    t = identity_stack(omega_arr.shape)
    for symbol in sequence(m):
        t = _matmul(matrices[symbol], t)
    return t


def semitrace(matrix: NDArray[np.complex128]) -> NDArray[np.complex128]:
    """Rm = (1/2) Tr(Tm)."""
    return 0.5 * (matrix[0, 0] + matrix[1, 1])


def determinant(matrix: NDArray[np.complex128]) -> NDArray[np.complex128]:
    return matrix[0, 0] * matrix[1, 1] - matrix[0, 1] * matrix[1, 0]


def analytic_r0_r1_r2(
    spec: SuperlatticeSpec,
    omega: ArrayLike,
    theta: float,
    polarization: Polarization,
) -> tuple[NDArray[np.complex128], NDArray[np.complex128], NDArray[np.complex128]]:
    """Ecuaciones (7)–(9) del paper. En TM, χ = ε en (9)."""
    omega_arr = _as_omega(omega)
    responses = material_response(spec, omega_arr)
    eps_a, mu_a = responses["A"]
    eps_b, mu_b = responses["B"]
    n_a = refractive_index(eps_a, mu_a)
    q = in_plane_wavevector(omega_arr, theta, n_a, spec.speed_of_light)
    q_a = longitudinal_wavevector(omega_arr, eps_a, mu_a, q, spec.speed_of_light)
    q_b = longitudinal_wavevector(omega_arr, eps_b, mu_b, q, spec.speed_of_light)
    chi_a = characteristic_parameter(polarization, eps_a, mu_a)
    chi_b = characteristic_parameter(polarization, eps_b, mu_b)
    a = spec.thickness_a
    b = spec.thickness_b
    r0 = np.cos(q_b * b)
    r1 = np.cos(q_a * a)
    # (1/2)[QA χB /(QB χA) + QB χA /(QA χB)]
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
    """Rm vía R0, R1, R2 y Rm = 2 R(m-1) R(m-2) − R(m-3)."""
    r0, r1, r2 = analytic_r0_r1_r2(spec, omega, theta, polarization)
    if m == 0:
        return r0
    if m == 1:
        return r1
    if m == 2:
        return r2
    prev3, prev2, prev1 = r0, r1, r2
    current = r2
    for _ in range(3, m + 1):
        current = 2.0 * prev1 * prev2 - prev3
        prev3, prev2, prev1 = prev2, prev1, current
    return current


def semitrace_by_product(
    spec: SuperlatticeSpec,
    m: int,
    omega: ArrayLike,
    theta: float,
    polarization: Polarization,
) -> NDArray[np.complex128]:
    return semitrace(cell_transfer_matrix(spec, m, omega, theta, polarization))


def cell_matrix_by_recurrence_blocks(
    spec: SuperlatticeSpec,
    m: int,
    omega: ArrayLike,
    theta: float,
    polarization: Polarization,
) -> NDArray[np.complex128]:
    """Tm = T(m-2) T(m-1), para comparar con el producto de capas."""
    omega_arr = _as_omega(omega)
    m_a, m_b = layer_matrices_ab(spec, omega_arr, theta, polarization)
    if m == 0:
        return m_b
    if m == 1:
        return m_a
    t_prev2 = m_b
    t_prev1 = m_a
    t_current = m_a
    for order in range(2, m + 1):
        t_current = _matmul(t_prev2, t_prev1)
        t_prev2, t_prev1 = t_prev1, t_current
        del order
    return t_current
