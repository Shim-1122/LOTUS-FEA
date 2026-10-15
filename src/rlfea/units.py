"""Unit system.

All internal quantities use a consistent mm-N-s system:

    length  mm
    force   N
    stress  N/mm^2 = MPa
    energy  N*mm = mJ
    density tonne/mm^3 (not used: quasi-static analyses)

Every user-facing input in a config file carries an explicit unit string and
is converted here.  Unknown units raise instead of being guessed.
"""
from __future__ import annotations

LENGTH = {"m": 1000.0, "mm": 1.0, "cm": 10.0, "in": 25.4, "inch": 25.4, "um": 1e-3}
FORCE = {"N": 1.0, "kN": 1000.0, "lbf": 4.4482216152605}
STRESS = {"MPa": 1.0, "kPa": 1e-3, "Pa": 1e-6, "GPa": 1e3, "N/mm2": 1.0, "psi": 6.894757293168e-3}
ANGLE = {"deg": 1.0, "rad": 57.29577951308232}
MOMENT = {"N*mm": 1.0, "N*m": 1000.0}
STIFFNESS = {"N/mm": 1.0, "N/m": 1e-3, "kN/m": 1.0}
DIMENSIONLESS = {"-": 1.0, "1": 1.0, "": 1.0}

_TABLES = {"length": LENGTH, "force": FORCE, "stress": STRESS, "angle": ANGLE,
           "moment": MOMENT, "stiffness": STIFFNESS, "dimensionless": DIMENSIONLESS}

INTERNAL = {"length": "mm", "force": "N", "stress": "MPa", "angle": "deg",
            "moment": "N*mm", "stiffness": "N/mm", "dimensionless": "-"}


class UnitError(ValueError):
    pass


def to_internal(value: float, unit: str, kind: str) -> float:
    """Convert ``value`` given in ``unit`` to the internal unit of ``kind``."""
    table = _TABLES.get(kind)
    if table is None:
        raise UnitError(f"unknown quantity kind {kind!r}")
    if unit not in table:
        raise UnitError(f"unit {unit!r} is not a recognised {kind} unit; allowed: {sorted(table)}")
    return float(value) * table[unit]


def quantity(spec, kind: str) -> float:
    """Parse ``{value: x, unit: u}`` or ``[x, u]`` into internal units.

    Bare numbers are rejected for dimensional kinds so that units are always
    explicit in configuration files.
    """
    if isinstance(spec, dict):
        return to_internal(spec["value"], spec["unit"], kind)
    if isinstance(spec, (list, tuple)) and len(spec) == 2:
        return to_internal(spec[0], spec[1], kind)
    if isinstance(spec, (int, float)) and kind == "dimensionless":
        return float(spec)
    raise UnitError(f"quantity of kind {kind!r} must carry an explicit unit, got {spec!r}")
