"""Standalone, shareable HTML viewer of a run's actual solver fields.

    rlfea export-web RUN_DIR --out viewer.html

Embeds the skin surface (every written output time, true-scale node
positions), bone and socket surfaces, nodal |u| and per-facet contact fields
(pressure, tangential traction, derived stick/slip).  Nothing is smoothed or
interpolated; facet values are drawn flat.  Never export patient-derived runs
to a public location.
"""
from __future__ import annotations

import json
import os
import pickle

import numpy as np


def _boundary_tris(conn_list):
    """Boundary triangles of a set of tet/wedge/hex blocks (faces used once)."""
    from collections import Counter
    faces = []
    for etype, conn in conn_list:
        if etype == "tet4":
            fl = [(0, 2, 1), (0, 1, 3), (1, 2, 3), (0, 3, 2)]
            for f in fl:
                faces.append(conn[:, f])
        elif etype == "hex8":
            q = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
            for a, b, c, d in q:
                faces.append(conn[:, [a, b, c]]); faces.append(conn[:, [a, c, d]])
        elif etype == "penta6":
            for f in [(0, 2, 1), (3, 4, 5)]:
                faces.append(conn[:, f])
            for a, b, c, d in [(0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)]:
                faces.append(conn[:, [a, b, c]]); faces.append(conn[:, [a, c, d]])
    F = np.vstack(faces)
    key = np.sort(F, axis=1)
    cnt = Counter(map(tuple, key))
    keep = np.array([cnt[tuple(k)] == 1 for k in key])
    return F[keep]


def export(run_path: str, out_html: str, max_states: int = 12) -> str:
    with open(os.path.join(run_path, "mesh", "model.pkl"), "rb") as fh:
        pk = pickle.load(fh)
    model, spec = pk["model"], pk["spec"]
    man = json.load(open(os.path.join(run_path, "manifest.json")))
    st = json.load(open(os.path.join(run_path, "status.json")))
    summ = json.load(open(os.path.join(run_path, "results", "summary.json")))
    load = json.load(open(os.path.join(run_path, "solver", "load_report.json")))
    cfg = json.load(open(os.path.join(run_path, "config.resolved.json")))
    import pyvista as pv
    X0 = model.nodes
    skin_tris = np.vstack([model.surfaces[s].conn for s in spec.log_surfaces] + ([model.surfaces["skin_free"].conn] if "skin_free" in model.surfaces else []))
    n_contact = sum(len(model.surfaces[s].conn) for s in spec.log_surfaces)
    bones = [(b.etype, b.conn) for b in model.blocks if b.region.startswith("bone")]
    bone_tris = _boundary_tris(bones) if bones else np.zeros((0, 3), int)
    sock = [(b.etype, b.conn) for b in model.blocks if b.region.startswith("socket") and b.formulation == "rigid"]
    sock_tris = _boundary_tris(sock) if sock else np.zeros((0, 3), int)
    used = np.unique(np.concatenate([skin_tris.ravel(), bone_tris.ravel(), sock_tris.ravel()]))
    remap = -np.ones(len(X0), int); remap[used] = np.arange(len(used))
    times = summ["times"]
    if len(times) > max_states:
        idx = np.unique(np.r_[np.linspace(0, len(times) - 1, max_states).round().astype(int)])
        times = [times[i] for i in idx]
    states = []
    T = load["tension"]["applied_N"] or 0
    steps = [s.name for s in spec.steps]
    for t in times:
        g = pv.read(os.path.join(run_path, "results", f"volume_t{t:.4f}.vtu"))
        U = np.asarray(g.point_data["displacement"])[used]
        rec = {"t": t, "label": (f"lace tension {T * t:.1f} N" if t <= 1 + 1e-9 else f"{steps[1] if len(steps) > 1 else ''} {100 * (t - 1):.0f}% of {cfg['load_cases'].get(steps[1], {}).get('displacement', '?')} mm, lace {T:.0f} N"),
               "u": np.round(U, 4).ravel().tolist()}
        vtp = os.path.join(run_path, "results", f"contact_t{t:.4f}.vtp")
        if os.path.exists(vtp):
            c = pv.read(vtp)
            rec["p"] = np.round(np.asarray(c.cell_data["contact_pressure"]) * 1000, 3).tolist()          # kPa
            rec["tt"] = np.round(np.asarray(c.cell_data["traction_tangential"]) * 1000, 3).tolist()       # kPa
            rec["ss"] = np.asarray(c.cell_data["stick_slip_status"]).astype(int).tolist()
        srec = next(s for s in summ["states"] if abs(s["time"] - t) < 1e-9)
        rec["umax"] = srec["max_tissue_displacement_mm"]
        rec["minJ"] = srec.get("deformation", {}).get("min_J")
        rec["N"] = srec.get("contact", {}).get("normal_force_N")
        rec["pmax"] = 1000 * srec.get("contact", {}).get("max_pressure_MPa", 0)
        rec["slip"] = srec.get("contact", {}).get("slip_fraction_of_contact_area")
        rec["imb"] = srec["equilibrium"]["force_imbalance_relative"]
        cup = srec["equilibrium"]["rigid_bodies"].get("rigid_cup", {}).get("F")
        rec["cupFz"] = cup[2] if cup else None
        states.append(rec)
    data = {"run_id": man["run_id"], "run_name": man.get("run_name"), "state": st["state"], "last_converged_time": st.get("last_converged_time"),
            "febio": man["solver"].get("version"), "synthetic": not man.get("patient_geometry"),
            "x0": np.round(X0[used], 3).ravel().tolist(), "skin": remap[skin_tris].ravel().tolist(), "n_contact": int(n_contact),
            "bone": remap[bone_tris].ravel().tolist(), "socket": remap[sock_tris].ravel().tolist(),
            "friction": cfg["contact"]["friction"], "warnings": man.get("warnings", []), "states": states}
    tpl = open(os.path.join(os.path.dirname(__file__), "web_template.html")).read()
    html = tpl.replace("/*__DATA__*/null", json.dumps(data, separators=(",", ":")))
    with open(out_html, "w") as fh:
        fh.write(html)
    return out_html
