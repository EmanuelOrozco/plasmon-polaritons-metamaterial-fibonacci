"""Física compartida por todos los métodos numéricos.

- ``constants`` y ``units``: constantes SI y conversiones GHz ↔ rad/s, mm ↔ m.
- ``materials``: medio A homogéneo y metamaterial B de Drude sin pérdidas.
- ``electromagnetics``: polarización TE/TM y vectores de onda q, Q.
- ``fibonacci``: palabras S_m, números F_m y longitud de la celda.
- ``superlattice``: ``SuperlatticeSpec`` validada (todo en SI).
- ``effective_medium``: ⟨ε⟩_m, ⟨μ⟩_m y la frecuencia del gap ⟨n⟩ = 0.
"""

from fibonacci_photonics.physics.electromagnetics import Polarization, as_polarization
from fibonacci_photonics.physics.materials import DrudeMetamaterial, HomogeneousMedium
from fibonacci_photonics.physics.superlattice import SuperlatticeSpec, spec_as_dict
from fibonacci_photonics.physics.units import mm_to_m, nu_ghz_from_omega, omega_from_nu_ghz

__all__ = [
    "DrudeMetamaterial",
    "HomogeneousMedium",
    "Polarization",
    "SuperlatticeSpec",
    "as_polarization",
    "mm_to_m",
    "nu_ghz_from_omega",
    "omega_from_nu_ghz",
    "spec_as_dict",
]
