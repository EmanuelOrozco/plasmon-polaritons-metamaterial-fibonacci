"""Alias de compatibilidad de ``fibonacci_photonics.plotting``."""

import importlib
import sys

sys.modules[__name__] = importlib.import_module("fibonacci_photonics.plotting")
