"""Alias de compatibilidad de ``fibonacci_photonics.tmm.transfer_matrix``."""

import importlib
import sys

sys.modules[__name__] = importlib.import_module("fibonacci_photonics.tmm.transfer_matrix")
