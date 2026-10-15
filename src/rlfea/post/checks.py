"""Run-level numerical checks: equilibrium, inversion, distortion, energy,
contact penetration and stabilisation contributions.

These are NUMERICAL VERIFICATION checks only.  None of them is evidence of
clinical validity, grip success or tissue safety.
"""
from __future__ import annotations

import numpy as np

from ..mesh.model import FEModel


def rigid_table(spec):
    """FEBio 4.13 rigid_body_data rows are keyed by the rigid MATERIAL id
    (1-based position in the Material section), verified against the logged
    rigid-body positions."""
    return {i + 1: n for i, (n, m) in enumerate(spec.materials.items()) if m.get("rigid")}


def equilibrium(model: FEModel, spec, state, rigid_ids) -> dict:
    """Sum of support reactions (constrained nodes + fixed/prescribed rigid
    bodies) and of applied loads.  For an internally balanced lace the applied
    loads sum to zero, so the support reactions must also sum to zero."""
    X = model.nodes
    out = {}
    R = state.get("R")
    supported = set()
    for bc in spec.bcs:
        supported.update(model.nodesets[bc["node_set"]].tolist())
    sup = np.array(sorted(supported), int)
    # Sign conventions (verified with the indentation/slide benchmarks):
    #   FEBio node "Rx" for a prescribed dof = force the structure exerts ON the support (= -constraint force);
    #   FEBio rigid_body_data "Fx" for a constrained rigid body = constraint force acting ON the body.
    # Both are converted to "external force acting on the model" before summing.
    Fn = -R[sup].sum(0) if (R is not None and len(sup)) else np.zeros(3)
    Mn = np.cross(X[sup] + np.nan_to_num(state["U"][sup]), -R[sup]).sum(0) if (R is not None and len(sup)) else np.zeros(3)
    constrained_rb = {r["rb"] for r in spec.rigid_bcs} | {r["rb"] for st in spec.steps for r in st.rigid_bcs}
    Fr, Mr = np.zeros(3), np.zeros(3)
    rows = {}
    if "rigid" in state:
        for rid, d in zip(state["rigid"]["ids"], state["rigid"]["data"]):
            name = rigid_ids.get(int(rid), str(rid))
            pos, F, M = d[0:3], d[6:9], d[9:12]
            rows[name] = {"position": pos.tolist(), "F": F.tolist(), "M": M.tolist()}
            if name in constrained_rb:
                Fr += F
                Mr += M + np.cross(pos, F)
    out["node_support_force_N"] = Fn.tolist()
    out["rigid_support_force_N"] = Fr.tolist()
    total = Fn + Fr
    # normalise by the largest force magnitude acting on any rigid body or support group
    mags = [np.linalg.norm(Fn), np.linalg.norm(Fr)] + [np.linalg.norm(v["F"]) for v in rows.values()]
    scale = max(mags + [1e-12])
    out["sum_support_force_N"] = total.tolist()
    out["force_scale_N"] = float(scale)
    out["force_imbalance_relative"] = float(np.linalg.norm(total) / scale)
    out["sum_support_moment_Nmm"] = (Mn + Mr).tolist()
    out["rigid_bodies"] = rows
    return out


def deformation_checks(J, stretch_max, stretch_min) -> dict:
    ok = np.isfinite(J)
    return {"min_J": float(np.nanmin(J)), "max_J": float(np.nanmax(J)), "inverted_elements": int((J[ok] <= 0).sum()),
            "elements_J_below_0p5": int((J[ok] < 0.5).sum()), "max_principal_stretch": float(np.nanmax(stretch_max)),
            "min_principal_stretch": float(np.nanmin(stretch_min)),
            "excessive_distortion_flag": bool((J[ok] < 0.5).any() or np.nanmax(stretch_max) > 2.0 or np.nanmin(stretch_min) < 0.4)}


def strain_energy(sed, V0, block_slices) -> dict:
    out = {"total_Nmm": float(np.nansum(sed * V0))}
    for name, (a, b) in block_slices.items():
        out[name] = float(np.nansum(sed[a:b] * V0[a:b]))
    return out
