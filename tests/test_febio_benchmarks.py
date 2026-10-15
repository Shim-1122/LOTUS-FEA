"""Light FEBio benchmarks (run in CI when FEBio is available)."""
import numpy as np
import pytest

from conftest import requires_febio
from rlfea.verification import benchmarks as B
from rlfea.verification import reference as R

pytestmark = [requires_febio, pytest.mark.febio]
HGO = dict(c=0.2014, k1=24.53, k2=0.1327, kappa=0.1404, gamma=41.0, k=20.14)


@pytest.mark.parametrize("F", [np.diag([1.25, 1.25 ** -0.5, 1.25 ** -0.5]), np.array([[1, 0.3, 0], [0, 1, 0], [0, 0, 1.0]])])
def test_hgo_matches_independent_reference(tmp_path, F):
    e1, e2 = np.eye(3)[0], np.eye(3)[1]
    ref = lambda F_: R.cauchy_from_W(lambda C: R.W_hgo(C, *[HGO[k] for k in ("c", "k1", "k2", "kappa", "gamma", "k")], e1, e2), F_)
    r = B.homogeneous_case(F, {"febio_type": "Holzapfel-Gasser-Ogden", "params": HGO}, ref, str(tmp_path))
    assert r["status"] == "converged" and r["max_rel_error"] < 1e-6


def test_ogden_matches_reference(tmp_path):
    p = dict(c1=0.00522, m1=4.79, c2=0.00522, m2=-4.79, k=0.522)
    F = np.diag([0.8, 0.8 ** -0.5, 0.8 ** -0.5])
    r = B.homogeneous_case(F, {"febio_type": "Ogden", "params": p}, lambda F_: R.cauchy_ogden(F_, [p["c1"], p["c2"]], [p["m1"], p["m2"]], p["k"]), str(tmp_path))
    assert r["max_rel_error"] < 1e-6


def test_rigid_rotation_is_stress_free(tmp_path):
    th = np.deg2rad(60)
    Rz = np.array([[np.cos(th), -np.sin(th), 0], [np.sin(th), np.cos(th), 0], [0, 0, 1.0]])
    r = B.homogeneous_case(Rz, {"febio_type": "Mooney-Rivlin", "params": {"c1": 0.05, "c2": 0, "k": 10}}, lambda F_: np.zeros((3, 3)), str(tmp_path))
    assert np.abs(np.array(r["sigma_febio"])).max() < 1e-9


def test_coulomb_sliding(tmp_path):
    r = B.contact_slide_case(str(tmp_path), mu=0.3)
    assert r["status"] == "converged" and abs(r["ratio_end"] - 0.3) / 0.3 < 0.03


def test_contact_opening(tmp_path):
    r = B.contact_opening_case(str(tmp_path))
    assert r["status"] == "converged" and r["Fz_end_over_Fz_max"] < 1e-6


def test_hinge_spring(tmp_path):
    r = B.rigid_hinge_case(str(tmp_path))
    assert r["rel_error"] < 0.01


def test_energy_balance(tmp_path):
    r = B.energy_balance_case(str(tmp_path))
    assert r["rel_diff"] < 0.02
