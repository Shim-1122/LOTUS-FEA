"""Export FEA surface fields for the LOTUS Gait Tester (site/gait-tester.html).

Reads the result VTU/VTP files of the stance runs and writes, for a set of
output times, the outer skin surface, the four petal pads and shells, the cup
wall and cushion, with per-node displacement, contact pressure, tangential
traction and stick/slip status. The data is quantised and base64-encoded and
inlined into the page between the FIELDS_BEGIN / FIELDS_END markers, so the
page stays a single file. FEA figures are inlined the same way.

usage: python scripts/export_gait_fields.py [--html site/gait-tester.html]
"""
from __future__ import annotations

import argparse
import base64
import io
import json
import os
import re

import numpy as np
import pyvista as pv

RUNS = {
    # key: (run dir, lace tension at t = 1 in lbf, output times to export)
    "base": ("runs/lotus_v5_stance_25lbf", 25.0, [0.25, 0.5, 0.75, 1.0, 1.1, 1.2, 1.3]),
    "relief": ("runs/lotus_v5_stance_60lbf_relief", 60.0, [0.1, 0.2, 0.3, 0.4, 0.5]),
}
FIGS = [("sections", "Mid-limb sections at 25 lbf, end of tightening: displacement, max principal log strain, hydrostatic pressure, max principal stress"),
        ("contact", "Skin contact at 25 lbf: contact pressure, tangential traction, stick/slip, gap (two views)")]
BLOCKS = {"pad": [8, 10, 12, 14], "shell": [9, 11, 13, 15], "cup": [16], "cushion": [17]}


def surf_faces(vol, ids):
    """Triangulated boundary faces of the cells whose block_id is in ids (global node ids)."""
    sub = vol.extract_cells(np.isin(vol.cell_data["block_id"], ids))
    s = sub.extract_surface(pass_pointid=True).triangulate()
    gid = s.point_data["vtkOriginalPointIds"]
    gid = sub.point_data["vtkOriginalPointIds"][gid] if "vtkOriginalPointIds" in sub.point_data else gid
    return gid[s.faces.reshape(-1, 4)[:, 1:]]


def b64(a):
    return base64.b64encode(np.ascontiguousarray(a).tobytes()).decode()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--html", default="site/gait-tester.html")
    a = ap.parse_args()
    out = {"note": "LOTUS-FEA surface fields; displacement in units of 0.005 mm (int16), pressure and traction in 0.1 kPa (uint16), slip 0 open / 1 stick / 2 slip",
           "sets": {}}
    for key, (run, lace, times) in RUNS.items():
        v0 = pv.read(f"{run}/results/volume_t{times[0]:.4f}.vtu")
        v0.point_data["vtkOriginalPointIds"] = np.arange(v0.n_points)
        # outer skin = skin-block boundary faces that are not shared with the fat block
        skin = surf_faces(v0, [0]); fat = set(map(tuple, np.sort(surf_faces(v0, [1]), 1)))
        skin = skin[[tuple(f) not in fat for f in np.sort(skin, 1)]]
        parts = {"skin": skin}
        for name, ids in BLOCKS.items():
            parts[name] = surf_faces(v0, ids)
        used = np.unique(np.concatenate([f.ravel() for f in parts.values()]))
        remap = -np.ones(v0.n_points, int); remap[used] = np.arange(len(used))
        X = v0.points[used]
        sset = {"run": os.path.basename(run), "lace_lbf_at_t1": lace, "n": int(len(used)),
                "X": b64(X.astype(np.float32)),
                "parts": {k: b64(remap[f].astype(np.uint16)) for k, f in parts.items()}, "states": []}
        summ = {round(s["time"], 4): s for s in json.load(open(f"{run}/results/summary.json"))["states"]}
        skin_nodes = np.unique(skin)
        for t in times:
            vol = pv.read(f"{run}/results/volume_t{t:.4f}.vtu")
            con = pv.read(f"{run}/results/contact_t{t:.4f}.vtp")
            U = vol.point_data["displacement"][used]
            # contact fields: face -> node
            F = con.faces.reshape(-1, 4)[:, 1:] if con.faces.size else con.cells.reshape(-1, 4)[:, 1:]
            p = np.zeros(v0.n_points); tt = np.zeros(v0.n_points); sl = np.zeros(v0.n_points); c = np.zeros(v0.n_points)
            cp, ct, ss = (con.cell_data[k] for k in ("contact_pressure", "traction_tangential", "stick_slip_status"))
            # node value = maximum over the contact faces touching the node, so the peak matches FEBio's face peak
            for j in range(3):
                np.maximum.at(p, F[:, j], cp); np.maximum.at(tt, F[:, j], ct); np.add.at(sl, F[:, j], ss); np.add.at(c, F[:, j], 1)
            c[c == 0] = 1
            # stick/slip status: rounded mean of the adjacent faces (majority), so slip area is not overstated
            P, Tt, SL = p[used], tt[used], np.rint(sl / c)[used]
            st = summ[round(t, 4)]
            sset["states"].append({
                "t": t, "lace_lbf": lace * min(t, 1.0), "push_mm": max(0.0, t - 1.0) * 10.0,
                "axial_N": abs(st["equilibrium"]["rigid_bodies"].get("rigid_cup", {}).get("F", [0, 0, 0])[2]),
                "peak_kPa": 1000 * st["contact"]["max_pressure_MPa"], "min_J": st["deformation"]["min_J"],
                "umax_skin_mm": float(np.linalg.norm(vol.point_data["displacement"][skin_nodes], axis=1).max()),
                "U": b64(np.clip(np.rint(U / 0.005), -32767, 32767).astype(np.int16)),
                "P": b64(np.clip(np.rint(P * 1e4), 0, 65535).astype(np.uint16)),
                "Tt": b64(np.clip(np.rint(Tt * 1e4), 0, 65535).astype(np.uint16)),
                "S": b64(SL.astype(np.uint8)),
            })
        out["sets"][key] = sset
        print(key, sset["n"], "nodes", {k: len(f) for k, f in parts.items()}, "faces")
    figs = []
    from PIL import Image
    for name, cap in FIGS:
        im = Image.open(f"{RUNS['base'][0]}/figures/{name}.png").convert("RGB")
        im.thumbnail((1400, 1400)); buf = io.BytesIO(); im.save(buf, "WEBP", quality=72)
        figs.append({"name": name, "caption": cap, "src": "data:image/webp;base64," + base64.b64encode(buf.getvalue()).decode()})
    out["figures"] = figs
    payload = "const FIELDS = " + json.dumps(out, separators=(",", ":")) + ";"
    html = open(a.html).read()
    new = re.sub(r"/\*FIELDS_BEGIN\*/.*?/\*FIELDS_END\*/", lambda m: "/*FIELDS_BEGIN*/" + payload + "/*FIELDS_END*/", html, flags=re.S)
    assert new != html or payload in html, "markers not found"
    open(a.html, "w").write(new)
    print("inlined", round(len(payload) / 1024), "KB into", a.html)


if __name__ == "__main__":
    main()
