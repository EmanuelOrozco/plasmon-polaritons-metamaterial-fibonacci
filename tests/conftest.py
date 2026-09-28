import os

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_name, "1")

from pathlib import Path  # noqa: E402

import pytest  # noqa: E402
import yaml

from fibonacci_photonics.core.params import SuperlatticeSpec, spec_from_mapping


@pytest.fixture
def paper_equal_plasma() -> SuperlatticeSpec:
    return spec_from_mapping(
        {
            "layer_a_thickness_mm": 12.0,
            "layer_b_thickness_mm": 12.0,
            "omega_e_over_2pi_ghz": 3.0,
            "omega_m_over_2pi_ghz": 3.0,
        }
    )


@pytest.fixture
def paper_split_plasma() -> SuperlatticeSpec:
    return spec_from_mapping(
        {
            "layer_a_thickness_mm": 12.0,
            "layer_b_thickness_mm": 12.0,
            "omega_e_over_2pi_ghz": 3.0,
            "omega_m_over_2pi_ghz": 1.0,
        }
    )


@pytest.fixture
def configs_dir() -> Path:
    return Path(__file__).resolve().parents[1] / "configs"
