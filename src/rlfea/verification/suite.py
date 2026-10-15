"""Verification suite: analytical / independent references for every solver
feature the limb model relies on.  Writes verification.json + a markdown
summary.  These are NUMERICAL verification results, not validation.
"""
from __future__ import annotations

import json
import os

import numpy as np

from ..loads import cable as CB
from . import benchmarks as B
from . import reference as R

HGO = dict(c=0.2014, k1=24.53, k2=0.1327, kappa=0.1404, gamma=41.0, k=20.14)
OGD = dict(c1=0.00522, m1=4.79, c2=0.00522, m2=-4.79, k=0.522)


def _deformations():
    lam = 1.25
    th = np.deg2rad(60)
    Rz = np.array([[np.cos(th), -np.sin(th), 0], [np.sin(th), np.cos(th), 0], [0, 0, 1.0]])
    U = np.array([[1.2, 0.05, 0.0], [0.05, 0.9, 0.02], [0.0, 0.02, 0.95]])
    return {"uniaxial_tension": np.diag([lam, lam ** -0.5, lam ** -0.5]),
            "uniaxial_compression": np.diag([0.8, 0.8 ** -0.5, 0.8 ** -0.5]),
            "simple_shear": np.array([[1, 0.3, 0], [0, 1, 0], [0, 0, 1.0]]),
            "volumetric": np.diag([1.02, 1.02, 1.02]),
            "rotated_stretch_RU": Rz @ U,
            "rigid_rotation_60deg": Rz}


def run_suite(out_dir: str, quick: bool = False) -> dict:
    os.makedirs(out_dir, exist_ok=True)
    checks = []

    def add(name, metric, value, tol, ok, detail=None):
        checks.append({"name": name, "metric": metric, "value": value, "tolerance": tol, "pass": bool(ok), "detail": detail})

    e1, e2 = np.array([1.0, 0, 0]), np.array([0, 1.0, 0])
    mats = {
        "HGO skin": ({"febio_type": "Holzapfel-Gasser-Ogden", "params": HGO},
                     lambda F: R.cauchy_from_W(lambda C: R.W_hgo(C, HGO["c"], HGO["k1"], HGO["k2"], HGO["kappa"], HGO["gamma"], HGO["k"], e1, e2), F)),
        "Ogden fat/muscle": ({"febio_type": "Ogden", "params": OGD},
                             lambda F: R.cauchy_ogden(F, [OGD["c1"], OGD["c2"]], [OGD["m1"], OGD["m2"]], OGD["k"])),
        "Mooney-Rivlin silicone": ({"febio_type": "Mooney-Rivlin", "params": {"c1": 0.05, "c2": 0.01, "k": 10.0}},
                                   lambda F: R.cauchy_from_W(lambda C: R.W_mooney_rivlin(C, 0.05, 0.01, 10.0), F)),
    }
    defs = _deformations()
    for mname, (mat, ref) in mats.items():
        for dname, F in defs.items():
            for et, dom in (("hex8", None), ("tet4", "ut4-solid")) if not quick else (("hex8", None),):
                r = B.homogeneous_case(F, mat, ref, os.path.join(out_dir, f"mat_{mname.split()[0]}_{dname}_{et}"), etype=et, domain_type=dom)
                if dname == "rigid_rotation_60deg":
                    s = np.abs(np.array(r.get("sigma_febio", [[np.inf]]))).max()
                    mu = mat["params"].get("c", mat["params"].get("c1", 1.0))
                    add(f"{mname}: rigid rotation 60 deg ({et}{'/' + dom if dom else ''})", "max |sigma| / mu", s / mu, 1e-8,
                        r["status"] == "converged" and s / mu < 1e-8)
                else:
                    add(f"{mname}: {dname} ({et}{'/' + dom if dom else ''})", "max rel. error vs independent W(C)",
                        r.get("max_rel_error", np.inf), 1e-6, r["status"] == "converged" and r.get("max_rel_error", 1) < 1e-6)
    # contact
    r = B.contact_slide_case(os.path.join(out_dir, "contact_slide"), mu=0.3)
    add("Frictional sliding (Coulomb), rigid plate on block", "|Fx|/|Fz| at full slip (target mu = 0.3)",
        r.get("ratio_end", np.nan), 0.03, r["status"] == "converged" and abs(r.get("ratio_end", 0) - 0.3) / 0.3 < 0.03, r)
    r = B.contact_opening_case(os.path.join(out_dir, "contact_opening"))
    add("Contact opening (compression only, separation)", "|Fz| end / |Fz| max after lifting past contact",
        r.get("Fz_end_over_Fz_max", np.nan), 1e-6, r["status"] == "converged" and r.get("Fz_end_over_Fz_max", 1) < 1e-6
        and r.get("max_pressure_end", 1) <= 1e-12, r)
    # rigid connectors
    r = B.rigid_hinge_case(os.path.join(out_dir, "rigid_hinge"))
    add("Rigid revolute joint + rigid spring + rigid moment", "relative error of hinge rotation",
        r.get("rel_error", np.nan), 0.01, r["status"] == "converged" and r.get("rel_error", 1) < 0.01, r)
    # cable routing force balance (no FEBio)
    rng = np.random.default_rng(1)
    pts = [CB.RoutePoint(f"b{k % 4}", rng.normal(size=3) * 50, f"p{k}") for k in range(12)]
    coms = {f"b{k}": rng.normal(size=3) for k in range(4)}
    for mu_g in (0.0, 0.2):
        bl = CB.body_loads(pts, 75.0, coms, mu_g)
        sF = np.linalg.norm(bl["sum_force_N"]); sM = np.linalg.norm(bl["sum_moment_about_origin_Nmm"])
        add(f"Cable routing force/moment balance (guide friction {mu_g})", "|sum F| N, |sum M| N mm", max(sF, sM), 1e-9, sF < 1e-9 and sM < 1e-7)
    # energy balance
    r = B.energy_balance_case(os.path.join(out_dir, "energy"))
    add("Energy balance, frictionless indentation", "|W_ext - U| / U (trapezoid over 10 must points)",
        r.get("rel_diff", np.nan), 0.02, r["status"] == "converged" and r.get("rel_diff", 1) < 0.02, r)
    # locking / near incompressibility: compression-dominated indentation at the
    # bulk/shear ratio used for the tissues (k = 100 mu); hex8 three-field mesh
    # refinement gives the reference.  Bending cantilever results are kept as
    # information (linear tets are poor in thin bending regardless of locking).
    lock = {}
    refs = {}
    for n in ((4, 8) if quick else (4, 8, 12)):
        refs[n] = B.indentation_locking_case(os.path.join(out_dir, f"ind_hex8_{n}"), "hex8", None, None, n, 100.0)
        lock[f"hex8 three-field n={n}"] = refs[n]
    nref = max(refs)
    F0 = refs[nref]["patch_force_N"]
    if len(refs) >= 3:
        ns = sorted(refs)
        f = [refs[k]["patch_force_N"] for k in ns]
        add("Indentation mesh refinement (hex8 three-field, k = 100 mu)", "relative change between the two finest meshes",
            abs(f[-1] - f[-2]) / abs(f[-1]), 0.05, abs(f[-1] - f[-2]) / abs(f[-1]) < 0.05, {"n": ns, "force": f})
    cases = [("tet4", None, None, 4, True), ("tet4", "ut4-solid", {"alpha": 0.05, "iso_stab": 0}, 8, False),
             ("tet10", None, None, 4, False)]
    for et, dom, par, n, expect_lock in cases:
        r = B.indentation_locking_case(os.path.join(out_dir, f"ind_{et}_{dom}_{n}"), et, dom, par, n, 100.0)
        key = f"{et}/{dom or 'standard'} n={n}"
        lock[key] = r
        ratio = r.get("patch_force_N", np.nan) / F0
        add(f"Near-incompressible indentation (k = 100 mu): {key}",
            f"patch force / hex8 n={nref} reference" + (" (locking EXPECTED for standard tet4)" if expect_lock else ""),
            ratio, None if expect_lock else 0.10, (ratio > 1.5) if expect_lock else abs(ratio - 1) < 0.10, r)
    r = B.indentation_locking_case(os.path.join(out_dir, "ind_ut4_k1000"), "tet4", "ut4-solid", {"alpha": 0.05, "iso_stab": 0}, 8, 1000.0)
    rh = B.indentation_locking_case(os.path.join(out_dir, "ind_hex8_k1000"), "hex8", None, None, 8, 1000.0)
    lock["INFO ut4 n=8 at k = 1000 mu (ratio to hex8 n=8)"] = r.get("patch_force_N", np.nan) / rh.get("patch_force_N", np.nan)
    rep = {"checks": checks, "locking_cases": lock, "scope": "numerical verification only"}
    with open(os.path.join(out_dir, "verification.json"), "w") as fh:
        json.dump(rep, fh, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    md = ["| check | metric | value | tolerance | verdict |", "|---|---|---|---|---|"]
    for c in checks:
        v = c["value"]
        md.append(f"| {c['name']} | {c['metric']} | {v:.3g} | {c['tolerance']} | {'PASS' if c['pass'] else 'FAIL'} |"
                  if isinstance(v, float) else f"| {c['name']} | {c['metric']} | {v} | {c['tolerance']} | {'PASS' if c['pass'] else 'FAIL'} |")
    open(os.path.join(out_dir, "verification.md"), "w").write("\n".join(md) + "\n")
    return rep
