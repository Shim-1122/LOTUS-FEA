"""Socket built from a supplied CAD assembly (STEP), e.g. LOTUS V5.

The CAD bodies replace the parametric petal shapes, but the FE representation
is kept the same as the parametric socket, so the mechanics, contact and
cable model are unchanged:

* each petal = deformable silicone pad (hex8, three-field) between the CAD
  pad's inner and outer faces + a rigid shell layer behind it, sampled on a
  structured (angle x height) grid by rays from the socket axis;
* each petal hinges on the CAD axle (axis = principal direction of the axle
  body, origin = its centroid), as a revolute joint to the rigid cup;
* the rigid cup wall = the CAD cup + inner hinge shield cavity, sampled the
  same way (contact with the skin between the distal cap and the pads);
* the distal cushion is NOT in the CAD: when enabled it is kept as the
  parametric layer under the distal skin cap, backed by the rigid cup
  (standing in for the Anclote puck seat).

Frame: the model frame has the socket axis on x = y = 0 (+z proximal) and
+y anterior; the caller supplies the 4x4 transform from CAD to model.
"""
from __future__ import annotations

import hashlib
import os
import pickle

import numpy as np

from ..mesh.model import ElementBlock, FEModel, Surface
from .socket import Petal, _quad_area
from .surface import TriSurface

PETAL_PARTS = ("anatomical_petal", "silicone_pad", "metal_moving_leaf", "inner_steel_strap", "smooth_axle")
CUP_PARTS = ("continuous_socket_rim_and_cup", "flexible_inner_hinge_shield")


def tessellate_step(path: str, size_mm: float = 3.0, cache_dir: str | None = None,
                    keep=PETAL_PARTS + CUP_PARTS) -> dict:
    """Tessellate the named solids of a STEP assembly (gmsh/OpenCASCADE, mm).
    Returns {part_name: (V, F)}; cached by file hash + size."""
    h = hashlib.sha256(open(path, "rb").read()).hexdigest()
    cache = os.path.join(cache_dir or os.path.dirname(os.path.abspath(path)), f".tess_{h[:16]}_{size_mm:g}.pkl")
    if os.path.exists(cache):
        with open(cache, "rb") as fh:
            out = pickle.load(fh)
        out["_meta"]["cached"] = True
        return out
    import gmsh
    import trimesh
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.option.setNumber("Geometry.OCCImportLabels", 1)
        gmsh.option.setString("Geometry.OCCTargetUnit", "MM")
        gmsh.model.occ.importShapes(path)
        gmsh.model.occ.synchronize()
        names = {t: gmsh.model.getEntityName(3, t).split("/")[-1] for _, t in gmsh.model.getEntities(3)}
        drop = [(3, t) for t, n in names.items() if not any(k in n for k in keep)]
        gmsh.model.occ.remove(drop, recursive=True)
        gmsh.model.occ.synchronize()
        gmsh.option.setNumber("Mesh.MeshSizeMax", size_mm)
        gmsh.option.setNumber("Mesh.MeshSizeMin", 0.5)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 12)
        gmsh.model.mesh.generate(2)
        tags, xyz, _ = gmsh.model.mesh.getNodes()
        X = np.asarray(xyz).reshape(-1, 3)
        pos = np.full(int(tags.max()) + 1, -1)
        pos[tags.astype(int)] = np.arange(len(tags))
        out = {"_meta": {"source": path, "sha256": h, "size_mm": size_mm, "cached": False, "parts": {}}}
        for _, t in gmsh.model.getEntities(3):
            tris = []
            for (_, f) in gmsh.model.getBoundary([(3, t)], oriented=False):
                et, _, en = gmsh.model.mesh.getElements(2, abs(f))
                for e, nn in zip(et, en):
                    if e == 2:
                        tris.append(np.asarray(nn, int).reshape(-1, 3))
            T = pos[np.vstack(tris)]
            m = trimesh.Trimesh(X, T, process=True)
            trimesh.repair.fix_normals(m)
            n = names[t]
            out[n] = (np.asarray(m.vertices), np.asarray(m.faces))
            out["_meta"]["parts"][n] = {"faces": int(len(m.faces)), "watertight": bool(m.is_watertight),
                                        "volume_mm3": float(abs(m.volume))}
    finally:
        gmsh.finalize()
    with open(cache, "wb") as fh:
        pickle.dump(out, fh)
    return out


def transform_parts(parts: dict, M: np.ndarray) -> dict:
    M = np.asarray(M, float)
    out = {"_meta": parts["_meta"]}
    for k, v in parts.items():
        if k == "_meta":
            continue
        V, F = v
        out[k] = (V @ M[:3, :3].T + M[:3, 3], F)
    return out


def _tm(parts, names):
    import trimesh
    ms = [trimesh.Trimesh(*parts[n], process=False) for n in names if n in parts]
    return trimesh.util.concatenate(ms) if len(ms) > 1 else ms[0]


def _ray_hits(mesh, z, th):
    """All hit radii (sorted) of horizontal rays from the axis at height z."""
    o = np.c_[np.zeros(len(th)), np.zeros(len(th)), np.full(len(th), z)]
    d = np.c_[np.cos(th), np.sin(th), np.zeros(len(th))]
    loc, ri, _ = mesh.ray.intersects_location(o, d, multiple_hits=True)
    hits = [[] for _ in th]
    r = np.linalg.norm(loc[:, :2], axis=1)
    for k, rr in zip(ri, r):
        hits[k].append(rr)
    return [sorted(h) for h in hits]


def _principal_axis(V):
    c = V.mean(0)
    _, _, vt = np.linalg.svd(V - c, full_matrices=False)
    return c, vt[0] / np.linalg.norm(vt[0])


def _rotate_about(P, origin, axis, ang):
    """Rodrigues rotation of points P about the line (origin, axis)."""
    k = axis / np.linalg.norm(axis)
    v = P - origin
    c, s_ = np.cos(ang), np.sin(ang)
    return origin + v * c + np.cross(k, v) * s_ + np.outer(v @ k, k) * (1 - c)


def _tube_rows(petal_mesh, thc_deg, zmin, zmax, dz=0.25, min_frac=0.6):
    """Hollow lace tubes running across a petal: horizontal rays from the axis
    that cross the petal wall twice (>= 4 hits) over most of its width."""
    th = np.deg2rad(thc_deg + np.linspace(-70, 70, 281))
    rows, cur = [], []
    for z in np.arange(zmin, zmax, dz):
        H = _ray_hits(petal_mesh, z, th)
        n4 = np.array([len(h) >= 4 for h in H]); n2 = np.array([len(h) >= 2 for h in H])
        full = n2.sum() > 0 and n4.sum() > min_frac * n2.sum()
        if full:
            cur.append((z, H, n4, n2))
        elif cur:
            rows.append(cur); cur = []
    if cur:
        rows.append(cur)
    out = []
    for r in rows:
        z, H, n4, n2 = r[len(r) // 2]                   # tube centre row
        k4 = np.where(n4)[0]; k2 = np.where(n2)[0]
        def opening(k):
            h = H[k]
            return float(0.5 * (h[1] + h[2]))
        out.append({"z": float(z), "diameter_mm": float(r[-1][0] - r[0][0] + dz),
                    "open_cw": (float(th[k4[0]]), opening(k4[0])), "open_ccw": (float(th[k4[-1]]), opening(k4[-1])),
                    "outer_cw": (float(th[k2[0]]), float(H[k2[0]][-1])), "outer_ccw": (float(th[k2[-1]]), float(H[k2[-1]][-1]))})
    return out


def _cad_lace_route(parts, petals, names, closing, rep):
    """Lace through the hollow tubes in the petals (CAD).  One strand from the
    dial on the first petal runs through the lower tubes of every petal, then
    crosses into the upper tubes and returns to the dial.  A row missing on a
    petal is passed OUTSIDE that petal (frictionless contact with its outer
    surface at the row height)."""
    def P(th, r, z):
        return np.array([r * np.cos(th), r * np.sin(th), z])
    tubes = {}
    for p, pn in zip(petals, names):
        m = _tm(parts, [f"{pn}_anatomical_petal"])
        tubes[p.name] = _tube_rows(m, p.theta_center_deg, m.bounds[0][2] + 1, m.bounds[1][2] - 1)
    zs = sorted({round(t["z"]) for v in tubes.values() for t in v})
    rows = []
    for zr in zs:
        if not rows or abs(zr - rows[-1]) > 3:
            rows.append(zr)
    dial = petals[0]
    m0 = _tm(parts, [f"{names[0]}_anatomical_petal"])
    z_dial = float(np.mean(rows))
    h0 = _ray_hits(m0, z_dial, np.deg2rad([dial.theta_center_deg]))[0]
    dial_xyz = P(np.deg2rad(dial.theta_center_deg), h0[-1] if h0 else 70.0, z_dial)
    from ..loads.cable import RoutePoint
    route, info = [], {"rows_z_mm": rows, "petals": {}}

    def openings(p, ri, zr):
        t = [x for x in tubes[p.name] if abs(x["z"] - zr) <= 3]
        if t:
            t = t[0]
            info["petals"].setdefault(p.name, {})[f"row{ri + 1}"] = "tube"
            return P(t["open_cw"][0], t["open_cw"][1], t["z"]), P(t["open_ccw"][0], t["open_ccw"][1], t["z"])
        # no tube at this height: the strand passes over the petal's outer surface
        m = _tm(parts, [f"{names[petals.index(p)]}_anatomical_petal"])
        th = np.deg2rad(p.theta_center_deg + np.linspace(-70, 70, 281))
        H = _ray_hits(m, zr, th)
        k2 = [i for i, h in enumerate(H) if len(h) >= 2]
        info["petals"].setdefault(p.name, {})[f"row{ri + 1}"] = "outside (no tube)"
        return P(th[k2[0]], H[k2[0]][-1] + 1.0, zr), P(th[k2[-1]], H[k2[-1]][-1] + 1.0, zr)

    # RevoFit path (Shyam, 2026-10-26): one strand leaves the dial, runs through the LOWER
    # tubes of every petal first, then crosses up into the UPPER tubes and returns to the
    # dial.  The petal without an upper tube is passed on the outside (by design).
    route.append(RoutePoint(dial.name, dial_xyz, "dial:out"))
    for ri, zr in enumerate(rows):
        for k, p in enumerate(petals):
            cw, ccw = openings(p, ri, zr)
            if not (ri == 0 and k == 0):                    # the dial feeds the first lower tube directly
                route.append(RoutePoint(p.name, cw, f"{p.name}:row{ri + 1}_cw"))
            route.append(RoutePoint(p.name, ccw, f"{p.name}:row{ri + 1}_ccw"))
    route.append(RoutePoint(dial.name, dial_xyz + np.array([0, 0, 1e-3]), "dial:in"))
    # move the route points with each petal's closing rotation
    for q in route:
        o, a, ang = closing.get(q.body, (None, None, 0.0))
        if ang:
            q.xyz = _rotate_about(q.xyz[None], o, a, ang)[0]
    info["n_points"] = len(route)
    info["note"] = ("one strand from the dial on the " f"{names[0]} petal through the lower tubes of every petal, "
                    "then crossing up through the upper tubes and back to the dial; frictionless tubes; a petal "
                    "without a tube at a row is passed on its outer surface")
    rep["lace_route"] = info
    return route


def build_socket_from_cad(model: FEModel, skin: TriSurface, cfg: dict, parts: dict) -> tuple[FEModel, dict]:
    """parts: CAD tessellation already in the model frame.  cfg as for the
    parametric socket, plus cad.petals = [names], cad.cap_height (mm)."""
    cad = cfg["cad"]
    h = float(cfg.get("element_size", 8.0))
    ts = float(cfg.get("shell_thickness", 4.0))
    margin_deg = float(cad.get("edge_margin_deg", 1.0))
    rep = {"status": f"socket from CAD {os.path.basename(parts['_meta']['source'])} (tentative; see notes)",
           "cad_sha256": parts["_meta"]["sha256"], "tessellation_mm": parts["_meta"]["size_mm"], "petals": {}}
    petals: list[Petal] = []
    names = cad.get("petals", ["anterior", "lateral", "posterior", "medial"])
    # order petals by angle starting at the first name (dial petal first)
    for i, pn in enumerate(names):
        pad = _tm(parts, [f"{pn}_silicone_pad"])
        shellparts = [f"{pn}_{k}" for k in ("anatomical_petal", "metal_moving_leaf", "inner_steel_strap")]
        c = pad.triangles_center
        thc = float(np.rad2deg(np.arctan2(*np.average(c[:, [1, 0]], axis=0, weights=pad.area_faces))) % 360)
        # angular and vertical extent of the pad's inner face
        ang = (np.rad2deg(np.arctan2(c[:, 1], c[:, 0])) - thc + 180) % 360 - 180
        lo, hi = np.percentile(ang, 1) + margin_deg, np.percentile(ang, 99) - margin_deg
        zlo, zhi = c[:, 2].min() + 0.5, c[:, 2].max() - 0.5
        # shrink until every ray of the grid hits the pad twice (inner + outer face)
        for _ in range(80):
            rmid = np.median(np.linalg.norm(c[:, :2], axis=1))
            nt = max(3, int(np.ceil(np.deg2rad(hi - lo) * rmid / h)) + 1)
            nz = max(3, int(np.ceil((zhi - zlo) / h)) + 1)
            th = np.deg2rad(thc + np.linspace(lo, hi, nt))
            Z = np.linspace(zlo, zhi, nz)
            H = [_ray_hits(pad, z, th) for z in Z]
            ok = np.array([[len(H[j][k]) >= 2 for j in range(nz)] for k in range(nt)])
            if ok.all():
                break
            bad = ~ok                                   # (nt, nz)
            inner = bad[1:-1, 1:-1]
            if bad[0].any():
                lo += 1.0
            if bad[-1].any():
                hi -= 1.0
            if not (bad[0].any() or bad[-1].any()):
                if bad[:, 0].any():
                    zlo += 1.0
                if bad[:, -1].any():
                    zhi -= 1.0
                if inner.any() and not (bad[:, 0].any() or bad[:, -1].any()):
                    raise RuntimeError(f"{pn}: pad has holes inside the sampling grid; refine tessellation")
        else:
            raise RuntimeError(f"{pn}: could not find a fully covered pad grid")
        R_in = np.array([[H[j][k][0] for j in range(nz)] for k in range(nt)])
        R_mid = np.array([[H[j][k][1] for j in range(nz)] for k in range(nt)])
        Zg = np.broadcast_to(Z, (nt, nz)).copy()
        trim = next((t for t in cad.get("trims", []) if t.get("petal") == pn), None)
        if trim:
            # swept top-edge trim (relief): per column, the top row is lowered to z_cut(theta) and
            # the column is re-sampled, so the pad and petal follow the trimmed profile smoothly
            u = (np.rad2deg(th) - (thc + lo)) / (hi - lo)              # 0 .. 1 across the petal
            if trim.get("from", "low_theta") == "high_theta":
                u = 1.0 - u
            f, r = float(trim["flat_frac"]), float(trim["ramp_frac"])
            w = np.clip((u - f) / max(r, 1e-9), 0.0, 1.0)
            w = 0.5 - 0.5 * np.cos(np.pi * w)                          # smooth ramp
            ztop_k = zhi - float(trim["depth_mm"]) * (1.0 - w)
            for k in range(nt):
                Zg[k] = np.linspace(zlo, ztop_k[k], nz)
                hk = [_ray_hits(pad, z, th[k:k + 1])[0] for z in Zg[k]]
                if any(len(x) < 2 for x in hk):
                    raise RuntimeError(f"{pn}: trimmed column {k} missed the pad")
                R_in[k] = [x[0] for x in hk]
                R_mid[k] = [x[1] for x in hk]
            rep["petals"].setdefault(f"petal{i + 1}", {})["trim"] = {
                **trim, "z_top_min_mm": float(ztop_k.min()), "z_top_max_mm": float(ztop_k.max()),
                "note": "approximate swept top-edge relief applied in the model (not from CAD)"}
        er = np.stack([np.cos(th), np.sin(th), np.zeros(nt)], axis=1)[:, None, :]
        base = np.concatenate([np.zeros((nt, nz, 2)), Zg[..., None]], axis=2)
        ids_in = model.add_nodes((base + R_in[..., None] * er).reshape(-1, 3)).reshape(nt, nz)
        ids_mid = model.add_nodes((base + R_mid[..., None] * er).reshape(-1, 3)).reshape(nt, nz)
        ids_out = model.add_nodes((base + (R_mid + ts)[..., None] * er).reshape(-1, 3)).reshape(nt, nz)

        def hexes(bot, top):
            return np.array([[bot[a, b], bot[a + 1, b], bot[a + 1, b + 1], bot[a, b + 1],
                              top[a, b], top[a + 1, b], top[a + 1, b + 1], top[a, b + 1]]
                             for a in range(nt - 1) for b in range(nz - 1)])

        name = f"petal{i + 1}"
        padc, shell = hexes(ids_in, ids_mid), hexes(ids_mid, ids_out)
        model.blocks.append(ElementBlock(f"{name}_pad", "hex8", padc, "silicone_pad", f"socket:{pn} petal:silicone pad (CAD)", "three-field"))
        model.blocks.append(ElementBlock(f"{name}_shell", "hex8", shell, f"rigid_{name}", f"socket:{pn} petal:rigid petal (CAD)", "rigid"))
        quads = np.array([[ids_in[a, b], ids_in[a, b + 1], ids_in[a + 1, b + 1], ids_in[a + 1, b]]
                          for a in range(nt - 1) for b in range(nz - 1)])
        model.surfaces[f"{name}_pad_contact"] = Surface(f"{name}_pad_contact", "quad4", quads)
        ax_c, ax_d = _principal_axis(parts[f"{pn}_smooth_axle"][0])
        # orient the hinge axis so that +rotation swings the petal top inward (right-hand rule, like the parametric model)
        ref = np.array([-np.sin(np.deg2rad(thc)), np.cos(np.deg2rad(thc)), 0.0])
        if ax_d @ ref < 0:
            ax_d = -ax_d
        X = model.nodes
        cad_shell = _tm(parts, [n for n in shellparts if n in parts])
        com = np.asarray(cad_shell.center_mass) if cad_shell.is_watertight else X[np.unique(shell)].mean(0)
        p = Petal(i, name, thc, 0.5 * (hi - lo), float(zlo), float(zhi), ids_in, ids_mid, ids_out, ax_c, ax_d, com)
        petals.append(p)
        Xn = model.nodes
        thick = np.linalg.norm(Xn[ids_mid] - Xn[ids_in], axis=2)
        pad_area = float(pad.area_faces[np.einsum("ij,ij->i", pad.face_normals, -np.c_[c[:, :2], 0 * c[:, 2]] /
                                                  np.linalg.norm(c[:, :2], axis=1)[:, None]) > 0.6].sum())
        rep["petals"].setdefault(name, {}).update({"cad_name": pn, "theta_center_deg": thc, "angular_span_deg": [thc + lo, thc + hi],
                               "z_range_mm": [float(zlo), float(zhi)], "grid": [int(nt), int(nz)],
                               "pad_thickness_mm": {"min": float(thick.min()), "median": float(np.median(thick)), "max": float(thick.max())},
                               "meshed_pad_contact_area_mm2": float(_quad_area(Xn, quads).sum()),
                               "cad_pad_inner_area_mm2": pad_area,
                               "hinge_origin": ax_c.tolist(), "hinge_axis": ax_d.tolist(), "com_cad": com.tolist()})
    # ---- closing rotation (start in contact) --------------------------------
    import trimesh
    sk = trimesh.Trimesh(skin.vertices, skin.faces, process=False)
    pq = trimesh.proximity.ProximityQuery(sk)
    snug = cad.get("snug_fit", {"percentile": 20.0, "interference": 0.1})
    closing = {}
    if snug:
        X = model.nodes
        for p in petals:
            ids = p.pad_nodes.ravel()
            clear = -pq.signed_distance(X[ids])                    # + = gap
            # rotation needed at each node to close its gap + interference (small-angle lever about the hinge)
            sgn = 1.0
            test = _rotate_about(X[ids], p.hinge_origin, p.hinge_axis, 1e-3)
            if np.mean(np.linalg.norm(test[:, :2], axis=1)) > np.mean(np.linalg.norm(X[ids][:, :2], axis=1)):
                sgn = -1.0                                         # this sign swings the pad outward; flip
            lever = np.linalg.norm(np.cross(X[ids] - p.hinge_origin, p.hinge_axis), axis=1)
            need = (clear + float(snug.get("interference", 0.1))) / np.maximum(lever, 1e-6)
            ang = sgn * float(np.percentile(need, float(snug.get("percentile", 20.0))))
            allp = np.unique(np.concatenate([p.pad_nodes.ravel(), p.pad_outer_nodes.ravel(), p.shell_outer_nodes.ravel()]))
            model.nodes[allp] = _rotate_about(X[allp], p.hinge_origin, p.hinge_axis, ang)
            p.com = _rotate_about(p.com[None], p.hinge_origin, p.hinge_axis, ang)[0]
            closing[p.name] = (p.hinge_origin, p.hinge_axis, ang)
            after = -pq.signed_distance(model.nodes[ids])
            rep["petals"][p.name]["closing_rotation_deg"] = float(np.degrees(abs(ang)))
            rep["petals"][p.name]["clearance_before_mm"] = {"min": float(clear.min()), "median": float(np.median(clear)), "max": float(clear.max())}
            rep["petals"][p.name]["clearance_after_mm"] = {"min": float(after.min()), "median": float(np.median(after)), "max": float(after.max())}
        rep["snug_fit"] = ("Each CAD petal (closed state) is rotated about its own hinge axle until the "
                           f"{snug.get('percentile', 20.0):g}th percentile of its pad nodes is {snug.get('interference', 0.1):g} mm into the skin, "
                           "so the solve starts with the pads in contact (the lace closes the remaining design clearance).")
    # ---- hinge stops: each petal may rotate at most N degrees inward ----------
    # A rigid tab hangs below each hinge (part of the petal). When the petal's top swings inward,
    # the tab's lower end swings outward and meets a stiff bumper fixed to the cup. The gap is set
    # so contact starts exactly at the rotation limit (measured from the CAD closed position by
    # default, so the snug-fit closing rotation counts toward it). Ordinary contact, frictionless.
    stop_pairs = []
    lim = cad.get("rotation_limit")
    if lim:
        lim_deg = float(lim.get("max_inward_deg", 2.0))
        L, wdt, tt = float(lim.get("lever_mm", 30.0)), float(lim.get("width_mm", 10.0)), float(lim.get("tab_thickness_mm", 6.0))
        bt, bk, bh = float(lim.get("bumper_mm", 2.0)), 3.0, float(lim.get("bumper_height_mm", 8.0))
        rep["rotation_limit"] = {"max_inward_deg": lim_deg, "measured_from": lim.get("measured_from", "vertical"), "lever_mm": L, "petals": {},
                                 "note": "petal lean = slope of the pad inner surface (radius vs height); the stop engages when the petal leans lim_deg inward of vertical"}

        def hex8(g, i, j, k):  # g: node-id grid [radial][axial][z]; hex between index i..i+1 etc, positive volume
            h = [g[i, j, k], g[i + 1, j, k], g[i + 1, j + 1, k], g[i, j + 1, k],
                 g[i, j, k + 1], g[i + 1, j, k + 1], g[i + 1, j + 1, k + 1], g[i, j + 1, k + 1]]
            P = model.nodes[h]
            if np.dot(np.cross(P[1] - P[0], P[3] - P[0]), P[4] - P[0]) < 0:
                h = [h[0], h[3], h[2], h[1], h[4], h[7], h[6], h[5]]
            return h

        def quad(g, i, toward):  # face at radial index i, normal along `toward`
            q = [g[i, 0, 0], g[i, 1, 0], g[i, 1, 1], g[i, 0, 1]]
            P = model.nodes[q]
            if np.dot(np.cross(P[1] - P[0], P[3] - P[0]), toward) < 0:
                q = q[::-1]
            return np.array([q])

        for p in petals:
            o = np.asarray(p.hinge_origin, float); a = np.asarray(p.hinge_axis, float); a /= np.linalg.norm(a)
            er = np.array([o[0], o[1], 0.0]); er /= np.linalg.norm(er)
            ez = np.array([0.0, 0.0, 1.0]); ea = np.cross(ez, er)        # horizontal, along the hinge
            ang0 = closing.get(p.name, (None, None, 0.0))[2]
            probe = o + 50.0 * ez                                        # a point above the hinge
            sgn = 1.0 if np.linalg.norm(_rotate_about(probe[None], o, a, 1e-3)[0, :2]) < np.linalg.norm(probe[:2]) else -1.0
            # Lean of the petal from vertical in its current (snug) pose: slope of the pad inner surface,
            # radius vs height, per column across the petal width; + = top leans outward.
            Pin = model.nodes[p.pad_nodes]                                   # (nt, nz, 3)
            rr = np.linalg.norm(Pin[..., :2], axis=-1); zz = Pin[..., 2]
            slopes = [np.polyfit(zz[k], rr[k], 1)[0] for k in range(rr.shape[0])]
            lean = float(np.degrees(np.arctan(np.median(slopes))))
            if lim.get("measured_from", "vertical") == "vertical":
                # the petal may tilt inward until it is lim_deg past vertical
                rem = max(lean + lim_deg, 0.0)
                used = abs(np.degrees(ang0))
            elif lim.get("measured_from") == "initial":
                # lim_deg of inward rotation from the petal's initial (simulation start, snug-fit) position
                used = abs(np.degrees(ang0))
                rem = lim_deg
            else:                                                            # "cad": lim_deg of rotation from the CAD pose
                used = abs(np.degrees(ang0))
                rem = max(lim_deg - used, 0.0)
            tab = np.array([[[o + r * er + s_ * ea + z * ez for z in (-L, 0.0)] for s_ in (-wdt / 2, wdt / 2)] for r in (-tt / 2, tt / 2)])
            if ang0:
                tab = _rotate_about(tab.reshape(-1, 3), o, a, ang0).reshape(tab.shape)
            pt = tab[1, :, 0].mean(0)                                    # outer lower edge of the tab
            moved = _rotate_about(pt[None], o, a, sgn * np.radians(rem))[0]
            gap = float(np.dot(moved - pt, er))
            r0 = float(np.dot(pt - o, er)) + gap
            zb = float(pt[2])
            stop = np.array([[[o + r * er + s_ * ea + (z - o[2]) * ez for z in (zb - 1.0, zb + bh)]
                              for s_ in (-wdt / 2 - 2, wdt / 2 + 2)] for r in (r0, r0 + bt, r0 + bt + bk)])
            tg = model.add_nodes(tab.reshape(-1, 3)).reshape(2, 2, 2)
            sg = model.add_nodes(stop.reshape(-1, 3)).reshape(3, 2, 2)
            model.blocks.append(ElementBlock(f"{p.name}_tab", "hex8", np.array([hex8(tg, 0, 0, 0)]), f"rigid_{p.name}",
                                             f"socket:{p.name} hinge-stop tab (rigid, part of petal)", "rigid"))
            model.blocks.append(ElementBlock(f"stopbump{p.index + 1}", "hex8", np.array([hex8(sg, 0, 0, 0)]), "petal_stop",
                                             f"socket:hinge-stop bumper for {p.name}", "standard"))
            model.blocks.append(ElementBlock(f"cupstop{p.index + 1}", "hex8", np.array([hex8(sg, 1, 0, 0)]), "rigid_cup",
                                             f"socket:hinge-stop backing for {p.name} (rigid cup)", "rigid"))
            model.surfaces[f"stop{p.index + 1}_bumper"] = Surface(f"stop{p.index + 1}_bumper", "quad4", quad(sg, 0, -er))
            model.surfaces[f"stop{p.index + 1}_tab"] = Surface(f"stop{p.index + 1}_tab", "quad4", quad(tg, 1, er))
            stop_pairs.append((f"stop{p.index + 1}_bumper", f"stop{p.index + 1}_tab", f"stop_{p.name}"))
            rep["rotation_limit"]["petals"][p.name] = {"closing_rotation_deg": used, "lean_from_vertical_deg": lean,
                                                       "lean_at_stop_deg": lean - rem, "remaining_inward_deg": rem,
                                                       "initial_gap_mm": gap, "inward_sign": sgn}
    # ---- rigid cup wall from the CAD cavity ---------------------------------
    cup = _tm(parts, CUP_PARTS)
    z_tip = float(skin.vertices[:, 2].min())
    cap_h = float(cad.get("cap_height", 8.0))
    z0 = z_tip + cap_h
    z1 = min(p.z_base for p in petals) - 0.5
    nth = max(12, int(np.ceil(2 * np.pi * 45.0 / h)))
    thg = np.linspace(0, 2 * np.pi, nth, endpoint=False)
    nzc = max(2, int(np.ceil((z1 - z0) / h)) + 1)
    Zc = np.linspace(z0, z1, nzc)
    Rc = np.array([[ (hh[0] if hh else np.nan) for hh in _ray_hits(cup, z, thg)] for z in Zc]).T   # (nth, nzc)
    if np.isnan(Rc).any():
        raise RuntimeError("cup cavity rays missed the CAD cup/shield; check the CAD-to-model transform")
    erc = np.stack([np.cos(thg), np.sin(thg), np.zeros(nth)], axis=1)[:, None, :]
    basec = np.concatenate([np.zeros((nth, nzc, 2)), np.broadcast_to(Zc, (nth, nzc))[..., None]], axis=2)
    push = np.zeros_like(Rc)
    for _ in range(8):                                           # radial push is not along the normal: iterate
        P_c = (basec + (Rc + push)[..., None] * erc).reshape(-1, 3)
        clear_c = -pq.signed_distance(P_c).reshape(nth, nzc)
        if clear_c.min() >= 0.099:
            break
        push = push + np.maximum(0.0, 0.1 - clear_c) * 1.2      # keep >= 0.1 mm initial gap
    rep["cup_wall_adjustment"] = {"nodes_moved_outward": int((push > 0).sum()), "max_mm": float(push.max()),
                                  "note": "cavity nodes where the CAD cup/hinge shield sat inside the skin were moved radially "
                                          "outward to 0.1 mm clearance (the hinge shield is flexible in reality)"}
    Rc = Rc + push
    c_in = model.add_nodes((basec + Rc[..., None] * erc).reshape(-1, 3)).reshape(nth, nzc)
    c_out = model.add_nodes((basec + (Rc + 3.0)[..., None] * erc).reshape(-1, 3)).reshape(nth, nzc)
    wall, wq = [], []
    for a in range(nth):
        a2 = (a + 1) % nth
        for b in range(nzc - 1):
            wall.append([c_in[a, b], c_in[a2, b], c_in[a2, b + 1], c_in[a, b + 1],
                         c_out[a, b], c_out[a2, b], c_out[a2, b + 1], c_out[a, b + 1]])
            wq.append([c_in[a, b], c_in[a, b + 1], c_in[a2, b + 1], c_in[a2, b]])
    model.blocks.append(ElementBlock("cup_wall", "hex8", np.array(wall), "rigid_cup", "socket:rigid cup wall + hinge shield (CAD)", "rigid"))
    model.surfaces["cup_wall_contact"] = Surface("cup_wall_contact", "quad4", np.array(wq))
    rep["cup_wall"] = {"z_range_mm": [z0, z1], "grid": [nth, nzc], "radius_mm": {"min": float(Rc.min()), "max": float(Rc.max())},
                       "note": "flexible_inner_hinge_shield is treated as rigid together with the cup (tentative)"}
    # ---- distal cushion + puck support (not in the CAD) ---------------------
    tc, tcup, g0 = float(cfg.get("cushion_thickness", 6.0)), float(cfg.get("cup_thickness", 4.0)), float(cfg.get("initial_gap", 0.0))
    V, F = skin.vertices, skin.faces
    fz = V[F][:, :, 2].max(axis=1)
    capf = fz <= z0
    if cad.get("distal_cushion", True) and capf.any():
        sub = F[capf]
        vids = np.unique(sub.ravel())
        remap = -np.ones(len(V), int); remap[vids] = np.arange(len(vids))
        loc = remap[sub]
        nrm = skin.vertex_normals(smooth_iters=3)[vids]
        Vc = V[vids]
        q_in = model.add_nodes(Vc + g0 * nrm)
        q_mid = model.add_nodes(Vc + (g0 + tc) * nrm)
        q_out = model.add_nodes(Vc + (g0 + tc + tcup) * nrm)
        model.blocks.append(ElementBlock("distal_cushion", "penta6", np.c_[q_in[loc], q_mid[loc]], "distal_cushion",
                                         "socket:distal cushion (not in CAD; parametric)", "three-field"))
        model.blocks.append(ElementBlock("distal_cup", "penta6", np.c_[q_mid[loc], q_out[loc]], "rigid_cup",
                                         "socket:puck seat behind cushion (rigid cup)", "rigid"))
        model.surfaces["cushion_contact"] = Surface("cushion_contact", "tri3", q_in[loc][:, [0, 2, 1]])
        rep["distal_cushion"] = {"faces": int(capf.sum()), "thickness_mm": tc, "cap_height_mm": cap_h,
                                 "note": "cushion is not in the CAD; modelled as a layer under the distal skin cap"}
    # ---- skin contact zones -------------------------------------------------
    skin_faces = model.surfaces["skin_outer"].conn
    Xs = model.nodes
    cen = Xs[skin_faces].mean(axis=1)
    ang = np.rad2deg(np.arctan2(cen[:, 1], cen[:, 0])) % 360
    zone = np.full(len(skin_faces), -1)
    zone[cen[:, 2] <= z0] = 0
    zone[(cen[:, 2] > z0) & (cen[:, 2] <= z1 + 0.5 * h)] = 99
    ctr = np.array([p.theta_center_deg for p in petals])
    nearest = (np.abs(((ang[:, None] - ctr[None, :]) + 180) % 360 - 180)).argmin(axis=1)
    zm = float(cfg.get("contact_zone_margin", 10.0))
    ztop = np.array([p.z_top_center for p in petals])
    inz = (cen[:, 2] > z1 + 0.5 * h) & (cen[:, 2] <= ztop[nearest] + zm)
    zone[inz] = nearest[inz] + 1
    model.surfaces["skin_cushion_zone"] = Surface("skin_cushion_zone", "tri3", skin_faces[zone == 0])
    model.surfaces["skin_cup_zone"] = Surface("skin_cup_zone", "tri3", skin_faces[zone == 99])
    for p in petals:
        model.surfaces[f"skin_{p.name}_zone"] = Surface(f"skin_{p.name}_zone", "tri3", skin_faces[zone == p.index + 1])
    model.surfaces["skin_free"] = Surface("skin_free", "tri3", skin_faces[zone < 0])
    model.meta["petals"] = petals
    model.meta["socket_cfg"] = dict(cfg)
    model.meta["extra_contact_pairs"] = [("skin_cup_zone", "cup_wall_contact", "contact_cup_wall")] + stop_pairs
    if cad.get("lace", "tubes") == "tubes":
        model.meta["cad_route"] = _cad_lace_route(parts, petals, names, closing, rep)
    # initial clearance of each socket surface to the skin (+ = gap)
    for nm in [f"{p.name}_pad_contact" for p in petals] + ["cup_wall_contact"]:
        q = Xs[np.unique(model.surfaces[nm].conn)]
        sd = -pq.signed_distance(q)   # trimesh: + inside -> flip so + = gap
        rep.setdefault("initial_clearance_mm", {})[nm] = {"min": float(sd.min()), "median": float(np.median(sd)), "max": float(sd.max())}
    return model, rep
