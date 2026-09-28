"""Alias de compatibilidad del paquete anterior (solo TMM).

El código nuevo debe importar ``fibonacci_photonics``. Este paquete se mantiene
para que la presentación y los notebooks existentes sigan funcionando sin cambios.
"""

from fibonacci_photonics import *  # noqa: F401,F403
from fibonacci_photonics import __all__ as _all
from fibonacci_photonics.tmm.dispersion import scan_dispersion  # noqa: F401

__all__ = [*_all, "scan_dispersion"]
