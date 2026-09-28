"""Los YAML de las figuras cargan y coinciden con el paper."""

from fibonacci_photonics.core.params import load_spec, spec_from_mapping


def test_figure_yamls_load(configs_dir):
    for name in [f"figure_0{i}.yaml" for i in range(1, 7)]:
        spec = load_spec(configs_dir / name)
        assert spec.thickness_a == spec.thickness_b
        assert spec.thickness_a == 0.012


def test_figure_01_plasma(configs_dir):
    spec = load_spec(configs_dir / "figure_01.yaml")
    assert spec.nu_e_ghz() == spec.nu_m_ghz() == 3.0


def test_figure_03_plasma(configs_dir):
    spec = load_spec(configs_dir / "figure_03.yaml")
    assert spec.nu_e_ghz() == 3.0
    assert spec.nu_m_ghz() == 1.0
