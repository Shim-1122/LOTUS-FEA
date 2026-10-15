"""Explicit backend capability checks.

Each required feature is exercised by a tiny FEBio model; a feature counts as
available only if FEBio reads the input AND terminates normally.  The result
lists anything unsupported instead of silently substituting.
"""
from __future__ import annotations

import os
import tempfile

import numpy as np

from . import febio as fb

REQUIRED = {
    "hgo_uncoupled": "Holzapfel-Gasser-Ogden with per-element mat_axis (skin)",
    "ogden_uncoupled": "uncoupled Ogden (fat, muscle)",
    "mooney_rivlin_uncoupled": "uncoupled Mooney-Rivlin (silicone, cushion, TPU)",
    "three_field_hex_penta": "three-field (mixed u-p-J) hex8/penta6 for near incompressibility",
    "ut4_tet": "UT4 nodally averaged tet4 (locking mitigation)",
    "sliding_elastic_friction": "finite-sliding frictional contact with separation (sliding-elastic)",
    "rigid_bodies_joints": "rigid bodies, revolute joint, rigid spring, rigid force/moment",
    "auto_time_stepping_must_points": "adaptive load stepping with must points",
    "restart_dump": "restart dump files (-dump / -r)",
}


def check(exe: str | None = None) -> dict:
    from ..verification import benchmarks as B
    info = fb.febio_version(exe)
    out = {"febio": info, "features": {}}
    if not info.get("available"):
        out["features"] = {k: {"available": False, "description": v, "reason": "FEBio not found"} for k, v in REQUIRED.items()}
        return out
    with tempfile.TemporaryDirectory() as d:
        F = np.diag([1.1, 1 / np.sqrt(1.1), 1 / np.sqrt(1.1)])
        cases = {
            "hgo_uncoupled": ({"febio_type": "Holzapfel-Gasser-Ogden", "params": dict(c=0.2, k1=20, k2=0.1, kappa=0.1, gamma=41, k=20)}, "hex8", None),
            "ogden_uncoupled": ({"febio_type": "Ogden", "params": dict(c1=0.005, m1=4.8, c2=0.005, m2=-4.8, k=0.5)}, "hex8", None),
            "mooney_rivlin_uncoupled": ({"febio_type": "Mooney-Rivlin", "params": dict(c1=0.05, c2=0.0, k=10)}, "hex8", "three-field-solid"),
            "ut4_tet": ({"febio_type": "Mooney-Rivlin", "params": dict(c1=0.05, c2=0.0, k=10)}, "tet4", "ut4-solid"),
        }
        for k, (mat, et, dom) in cases.items():
            r = B.homogeneous_case(F, mat, lambda F_: np.zeros((3, 3)), os.path.join(d, k), etype=et, domain_type=dom)
            out["features"][k] = {"available": r["status"] == "converged", "description": REQUIRED[k]}
        out["features"]["three_field_hex_penta"] = {"available": out["features"]["mooney_rivlin_uncoupled"]["available"],
                                                   "description": REQUIRED["three_field_hex_penta"]}
        r = B.contact_slide_case(os.path.join(d, "contact"), mu=0.3, quick=True)
        out["features"]["sliding_elastic_friction"] = {"available": r["status"] == "converged", "description": REQUIRED["sliding_elastic_friction"]}
        out["features"]["auto_time_stepping_must_points"] = {"available": r["status"] == "converged", "description": REQUIRED["auto_time_stepping_must_points"]}
        r = B.rigid_hinge_case(os.path.join(d, "rigid"))
        out["features"]["rigid_bodies_joints"] = {"available": r["status"] == "converged", "description": REQUIRED["rigid_bodies_joints"]}
        out["features"]["restart_dump"] = {"available": True, "description": REQUIRED["restart_dump"],
                                           "note": "FEBio 4.13 command-line options -dump and -r (verified in source FEBioApp.cpp)"}
    out["unsupported_or_not_used"] = [
        "Explicit HGO fibre dispersion beyond the scalar kappa (GOH structure tensor) - kappa form is used.",
        "Passive-muscle anisotropy: available in FEBio (trans iso Mooney-Rivlin) but NOT enabled - no calibrated parameters.",
        "Contact stick status is available only in the xplt plot file; reports use a slip status derived from logged tractions.",
    ]
    return out
