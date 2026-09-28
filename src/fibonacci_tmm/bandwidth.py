"""Alias de compatibilidad de ``fibonacci_photonics.analysis.bandwidth``."""

import importlib
import sys

sys.modules[__name__] = importlib.import_module("fibonacci_photonics.analysis.bandwidth")
