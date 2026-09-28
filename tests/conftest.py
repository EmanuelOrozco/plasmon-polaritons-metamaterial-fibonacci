import os

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_name, "1")
os.environ.setdefault("MPLBACKEND", "Agg")

from pathlib import Path  # noqa: E402

import pytest  # noqa: E402

from fibonacci_photonics.config import spec_from_mapping  # noqa: E402
from fibonacci_photonics.physics.superlattice import SuperlatticeSpec  # noqa: E402


def paper_spec(nu_e_ghz: float, nu_m_ghz: float) -> SuperlatticeSpec:
    """Superred del PRB: a = b = 12 mm, aire y Drude con ν_e, ν_m en GHz."""
    return spec_from_mapping(
        {
            "layer_a_thickness_mm": 12.0,
            "layer_b_thickness_mm": 12.0,
            "omega_e_over_2pi_ghz": nu_e_ghz,
            "omega_m_over_2pi_ghz": nu_m_ghz,
        }
    )


@pytest.fixture
def paper_equal_plasma() -> SuperlatticeSpec:
    """Figs. 1–2: ν_e = ν_m = 3 GHz."""
    return paper_spec(3.0, 3.0)


@pytest.fixture
def paper_split_plasma() -> SuperlatticeSpec:
    """Figs. 3–6: ν_e = 3 GHz, ν_m = 1 GHz."""
    return paper_spec(3.0, 1.0)


@pytest.fixture
def configs_dir() -> Path:
    return Path(__file__).resolve().parents[1] / "configs"
