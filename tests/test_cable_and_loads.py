import numpy as np
import pytest

from rlfea.loads import cable as CB


def _ring(n=6, r=60.0, z=100.0):
    pts = []
    for i in range(n):
        th = 2 * np.pi * i / n
        pts.append(CB.RoutePoint(f"p{i}", np.array([r * np.cos(th - 0.2), r * np.sin(th - 0.2), z]), f"p{i}a"))
        pts.append(CB.RoutePoint(f"p{i}", np.array([r * np.cos(th + 0.2), r * np.sin(th + 0.2), z]), f"p{i}b"))
    return pts + [CB.RoutePoint("p0", pts[0].xyz + [0, 0, 1e-3], "end")]


def test_internal_balance_and_no_multiplier():
    pts = _ring()
    coms = {f"p{i}": np.zeros(3) for i in range(6)}
    bl = CB.body_loads(pts, 75.0, coms)
    assert np.allclose(bl["sum_force_N"], 0, atol=1e-9)
    assert np.allclose(bl["sum_moment_about_origin_Nmm"], 0, atol=1e-7)
    # each span carries exactly the strand tension (no factor)
    assert all(abs(s["tension_N"] - 75.0) < 1e-12 for s in bl["spans"])
    # regular hexagonal lace: inward resultant on a non-dial petal = 2 T sin(half exterior angle)
    F = np.linalg.norm(bl["bodies"]["p2"]["F"])
    assert F > 0 and F < 2 * 75.0


def test_capstan_friction_decays_toward_middle():
    pts = _ring()
    T = CB.lace_tensions(pts, 75.0, mu_guide=0.2)
    assert T[0] == pytest.approx(75.0) and T[-1] == pytest.approx(75.0)
    assert T.min() < 75.0
    bl = CB.body_loads(pts, 75.0, {f"p{i}": np.zeros(3) for i in range(6)}, 0.2)
    assert np.allclose(bl["sum_force_N"], 0, atol=1e-9)


def test_petal_moment_balance_sign():
    class P:  # minimal petal: hinge at base on +x side, axis = tangent (+y)
        hinge_origin = np.array([60.0, 0, 0]); hinge_axis = np.array([0, 1.0, 0])
    load = {"point_forces": [{"xyz": [62.0, 0, 100.0], "force_N": [-10.0, 0, 0]}]}   # cable pulls petal top inward
    N = CB.petal_pad_force_moment_balance(P, load, np.array([60.0, 0, 50.0]), np.array([1.0, 0, 0]))
    assert N == pytest.approx(20.0)   # lever 100 vs 50 -> compressive pad force 20 N
