"""Alias de compatibilidad de ``fibonacci_photonics.analysis.plasmon_modes``."""

import importlib
import sys

sys.modules[__name__] = importlib.import_module("fibonacci_photonics.analysis.plasmon_modes")
