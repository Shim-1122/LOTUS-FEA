import numpy as np
import pytest

from rlfea.geometry import repair as R
from rlfea.geometry.surface import TriSurface, audit_assembly, audit_body, winding_number
from rlfea.geometry.synthetic import synthetic_limb


def _box(s=10.0, c=(0, 0, 0)):
    import trimesh
    m = trimesh.creation.box(extents=(s, s, s))
    m.apply_translation(c)
    return TriSurface("box", np.asarray(m.vertices), np.asarray(m.faces))


def test_transform_logging_and_units():
    b = _box()
    M = np.eye(4); M[:3, 3] = [1, 2, 3]
    t = R.transform(b, M, "test shift")
    assert np.allclose(t.vertices - b.vertices, [1, 2, 3])
    assert t.history[-1]["op"] == "transform" and "before" in t.history[-1]
    assert b.history == []                      # source object untouched
    # mirror flips orientation -> faces reversed so volume stays positive
    Mm = np.diag([-1.0, 1, 1, 1])
    assert R.transform(b, Mm, "mirror").volume() == pytest.approx(b.volume())


def test_audit_detects_defects():
    b = _box()
    a = audit_body(b)
    assert a["watertight"] and a["status"] == "VALID_CLOSED_SURFACE"
    holed = TriSurface("holed", b.vertices, b.faces[:-1])
    assert not audit_body(holed)["watertight"]
    inv = TriSurface("inv", b.vertices, b.faces[:, ::-1])
    assert "inward oriented (negative volume)" in audit_body(inv)["defects"]
    empty = TriSurface("empty", np.zeros((0, 3)), np.zeros((0, 3), int))
    assert audit_body(empty)["status"] == "EMPTY"


def test_self_intersection_detected():
    a, b = _box(), _box(c=(3, 0, 0))
    merged = TriSurface("m", np.vstack([a.vertices, b.vertices]), np.vstack([a.faces, b.faces + len(a.vertices)]))
    assert audit_body(merged)["self_intersecting_pairs"] > 0


def test_assembly_containment_and_duplicates():
    outer, inner, outside = _box(20.0), _box(4.0), _box(4.0, c=(30, 0, 0))
    rep = audit_assembly(outer, [inner, outside])
    assert rep["bodies"]["box"]["fraction_vertices_inside_envelope"] in (0.0, 1.0)
    w = winding_number(outer, np.array([[0, 0, 0], [30, 0, 0]], float))
    assert w[0] == pytest.approx(1, abs=1e-9) and abs(w[1]) < 1e-9


def test_synthetic_limb_is_valid_and_nested():
    g = synthetic_limb()
    for s in g.values():
        a = audit_body(R.basic_cleanup(s))
        assert a["status"] == "VALID_CLOSED_SURFACE", (s.name, a["defects"])
    rep = audit_assembly(g["skin"], [g["tibia"], g["fibula"]])
    for k in ("tibia", "fibula"):
        assert rep["bodies"][k]["fraction_vertices_inside_envelope"] == 1.0
        assert rep["bodies"][k]["intersects_envelope_tri_pairs"] == 0
