"""Los YAML validan con sus esquemas Pydantic y los errores se detectan al cargar."""

import copy

import pytest

from fibonacci_photonics.config import (
    FIGURE_IDS,
    ConfigError,
    load_figure,
    load_model,
    load_pwe_numerics,
    load_spec,
    load_yaml,
    spec_from_mapping,
)
from fibonacci_photonics.config.schema import (
    BookBenchmarkConfig,
    ConvergenceConfig,
    EfficiencyConfig,
    SpectralPollutionConfig,
)
from fibonacci_photonics.physics.constants import SPEED_OF_LIGHT, SPEED_OF_LIGHT_ROUNDED


@pytest.mark.parametrize("figure_id", FIGURE_IDS)
def test_every_figure_and_its_pwe_numerics_validate(figure_id):
    setup = load_figure(figure_id)
    assert setup.spec.thickness_a == setup.spec.thickness_b == pytest.approx(0.012)
    assert (setup.path.parent.name, setup.tmm_path.parent.name) == ("physics", "tmm")
    assert setup.path.name == setup.tmm_path.name == f"{figure_id}.yaml"
    numerics = load_pwe_numerics(figure_id).model
    orders = getattr(numerics, "fibonacci_orders", None) or setup.config.fibonacci_orders
    assert set(orders) <= set(setup.config.fibonacci_orders)
    numerics.check_orders(list(orders))


@pytest.mark.parametrize(
    "path, model",
    [
        ("pwe/convergence.yaml", ConvergenceConfig),
        ("pwe/book_benchmark.yaml", BookBenchmarkConfig),
        ("pwe/spectral_pollution.yaml", SpectralPollutionConfig),
        ("comparison/efficiency.yaml", EfficiencyConfig),
    ],
)
def test_study_configs_validate(path, model):
    assert isinstance(load_model(path, model).model, model)


def test_plasma_frequencies_of_the_paper():
    assert load_spec("figure_01.yaml").nu_m_ghz() == pytest.approx(3.0)
    split = load_spec("figure_03.yaml")
    assert split.nu_e_ghz() == pytest.approx(3.0)
    assert split.nu_m_ghz() == pytest.approx(1.0)


def test_default_spec_is_the_paper_superlattice():
    assert load_spec() == load_spec("figure_01.yaml")


@pytest.mark.parametrize(
    "option, expected", [("si", SPEED_OF_LIGHT), ("rounded", SPEED_OF_LIGHT_ROUNDED), (2.5e8, 2.5e8)]
)
def test_speed_of_light_option(option, expected):
    spec = spec_from_mapping(
        {"layer_a_thickness_mm": 1, "layer_b_thickness_mm": 1, "omega_e_over_2pi_ghz": 1,
         "omega_m_over_2pi_ghz": 1, "speed_of_light": option}
    )
    assert spec.speed_of_light == expected


@pytest.mark.parametrize("figure_id", FIGURE_IDS)
def test_tmm_grids_are_not_physics(figure_id):
    """La física es común a los dos métodos; las mallas de la TMM viven aparte."""
    physics, grids = load_yaml(f"physics/{figure_id}.yaml"), load_yaml(f"tmm/{figure_id}.yaml")
    assert not set(physics) & set(grids)
    assert all("points" in key or key.startswith(("frequency_", "nu_m_")) for key in grids)


def test_tmm_grid_cannot_override_physics(tmp_path, monkeypatch):
    from fibonacci_photonics.config import loader

    for folder in ("physics", "tmm"):
        (tmp_path / folder).mkdir()
    (tmp_path / "physics" / "figure_01.yaml").write_text(
        loader.resolve("physics/figure_01.yaml").read_text(encoding="utf-8"), encoding="utf-8"
    )
    (tmp_path / "tmm" / "figure_01.yaml").write_text("frequency_points: 100\nfrequency_max_ghz: 9.0\n")
    monkeypatch.setattr(loader, "configs_dir", lambda: tmp_path)
    with pytest.raises(ConfigError, match="frequency_max_ghz"):
        load_figure("figure_01")


def _figure_05() -> dict:
    return load_yaml("physics/figure_05.yaml") | load_yaml("tmm/figure_05.yaml")


@pytest.mark.parametrize(
    "mutate",
    [
        pytest.param(lambda d: d.update(unknown_key=1), id="clave-desconocida"),
        pytest.param(lambda d: d.update(layer_a_thickness_mm=-12.0), id="espesor-negativo"),
        pytest.param(lambda d: d.update(omega_m_over_2pi_ghz=-1.0), id="plasma-negativo"),
        pytest.param(lambda d: d.update(epsilon_a=0.0), id="epsilon-nula"),
        pytest.param(lambda d: d.update(frequency_points=1), id="malla-de-un-punto"),
        pytest.param(lambda d: d.update(thetas_rad=[0.1, 2.0]), id="angulo-mayor-que-pi/2"),
        pytest.param(lambda d: d.update(theta_labels=["pi/12"]), id="etiquetas-sin-angulo"),
        pytest.param(lambda d: d["frequency_windows_ghz"].update({"pi/3": [1.0, 0.9]}), id="ventana-invertida"),
        pytest.param(lambda d: d["frequency_windows_ghz"].pop("pi/3"), id="ventana-faltante"),
        pytest.param(lambda d: d.update(polarization="XY"), id="polarizacion"),
        pytest.param(lambda d: d.update(fibonacci_orders=[]), id="sin-ordenes"),
    ],
)
def test_invalid_figure_config_is_rejected(mutate):
    from fibonacci_photonics.config.loader import FIGURE_SCHEMAS, validate

    data = copy.deepcopy(_figure_05())
    mutate(data)
    with pytest.raises(ConfigError):
        validate(FIGURE_SCHEMAS["figure_05"], data)


def test_unknown_figure_and_missing_file():
    with pytest.raises(ConfigError):
        load_figure("figure_99")
    with pytest.raises(ConfigError):
        load_yaml("no_existe.yaml")


def test_harmonics_mapping_must_cover_orders():
    numerics = load_pwe_numerics("figure_05").model
    with pytest.raises(ValueError):
        type(numerics).model_validate(
            {**numerics.model_dump(), "harmonics_per_layer": {2: 16}}
        ).check_orders([2, 3])
