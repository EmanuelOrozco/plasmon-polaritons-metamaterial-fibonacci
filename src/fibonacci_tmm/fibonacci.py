"""Alias de compatibilidad de ``fibonacci_photonics.core.fibonacci``."""

import importlib
import sys

sys.modules[__name__] = importlib.import_module("fibonacci_photonics.core.fibonacci")
