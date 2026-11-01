import numpy as np
import pytest

from rlfea import config as C
from rlfea.units import UnitError, quantity, to_internal


def test_unit_conversions():
    assert to_internal(1, "m", "length") == 1000
    assert to_internal(1, "in", "length") == 25.4
    assert to_internal(5.22, "kPa", "stress") == pytest.approx(0.00522)
    assert to_internal(1, "lbf", "force") == pytest.approx(4.4482216)


def test_units_must_be_explicit():
    with pytest.raises(UnitError):
        quantity(75, "force")
    with pytest.raises(UnitError):
        quantity({"value": 1, "unit": "furlong"}, "length")


def test_rated_tension_requires_source_and_meaning(tmp_path):
    raw = C.load("configs/synthetic_human.yaml")
    raw["cable"]["rated_tension"] = {"value": 100, "unit": "N"}
    with pytest.raises(C.ConfigError):
        C.resolve(raw)
    raw["cable"]["rated_tension"] = {"value": 100, "unit": "N", "source": "x", "meaning": "dial torque"}
    with pytest.raises(C.ConfigError):
        C.resolve(raw)


def test_example_tension_is_flagged():
    cfg, warns = C.resolve(C.load("configs/synthetic_human.yaml"))
    assert any("EXAMPLE" in w for w in warns)
    assert any("NOT run" in w for w in warns)


def test_suction_requires_justification():
    raw = C.load("configs/synthetic_human.yaml")
    raw["load_cases"]["suction"] = {"enabled": True, "pressure": {"value": 10, "unit": "kPa"}}
    with pytest.raises(C.ConfigError):
        C.resolve(raw)


def test_applied_load_force_control():
    raw = C.load("configs/synthetic_human.yaml")
    raw["load_cases"]["pullout"] = {"enabled": False}
    raw["load_cases"]["applied_load"] = {"enabled": True, "axial_force": {"value": 1.2, "unit": "kN"},
                                         "moment_x": {"value": 61, "unit": "N*m"}}
    with pytest.raises(C.ConfigError):          # a source is required
        C.resolve(raw)
    raw["load_cases"]["applied_load"]["source"] = "test"
    cfg, _ = C.resolve(raw)
    al = cfg["load_cases"]["applied_load"]
    assert al["axial_force"] == pytest.approx(1200.0)
    assert al["moment_x"] == pytest.approx(61000.0)     # N*mm internally
    assert al["moment_y"] == 0.0
