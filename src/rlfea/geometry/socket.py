"""Parametric petal-socket geometry built fresh around a limb surface.

PROVISIONAL DESIGN MODEL.  No current socket CAD was supplied (the V6 design
is explicitly excluded), so this module generates a parametric stand-in for
the described concept:

* N organic load-distributing petals (rigid shells) with silicone interface
  pads (deformable hexahedra) on their inner faces,
* posterior petals trimmed lower for hamstring clearance,
* a revolute hinge at each petal base (ribbed hinge supports),
* a rigid distal cup with a deformable distal cushion.

The Anclote suction puck, flat plate and Bulldog adapter sit below the
cushion and are represented only as the rigid support of the cup (they are
not meshed).  Every dimension is a configuration value; none is a verified
component dimension unless the configuration says so.

Frame: limb axis = +z (proximal), origin at the distal skin apex, +y anterior.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..mesh.model import ElementBlock, FEModel, Surface
from .surface import TriSurface


@dataclass
class Petal:
    index: int
    name: str
    theta_center_deg: float
    half_width_deg: float
    z_base: float
    z_top_center: float
    pad_nodes: np.ndarray        # (nt, nz) global node ids, pad inner surface
    pad_outer_nodes: np.ndarray  # (nt, nz)
    shell_outer_nodes: np.ndarray
    hinge_origin: np.ndarray
    hinge_axis: np.ndarray
    com: np.ndarray
    guide_points: dict = field(default_factory=dict)


def _ray_radius(skin_tm, z, thetas):
    """Distance from the limb axis to the skin along radial rays at height z."""
    o = np.c_[np.zeros(len(thetas)), np.zeros(len(thetas)), np.full(len(thetas), z)]
    d = np.c_[np.cos(thetas), np.sin(thetas), np.zeros(len(thetas))]
    loc, iray, _ = skin_tm.ray.intersects_location(o, d, multiple_hits=True)
    r = np.full(len(thetas), np.nan)
    if len(loc):
        dist = np.linalg.norm(loc[:, :2], axis=1)
        for k in range(len(thetas)):
            m = iray == k
            if m.any():
                r[k] = dist[m].max()   # outermost crossing = skin
    return r


def build_socket(model: FEModel, skin: TriSurface, cfg: dict) -> tuple[FEModel, dict]:
    """Append socket components to ``model`` (in place) and return a report.

    cfg (lengths mm, angles deg):
      n_petals, petal_centers_deg (optional list), gap_deg, z_base, z_top,
      posterior_clearance_drop, posterior_center_deg (default 270),
      posterior_clearance_halfwidth_deg, tip_drop, pad_thickness,
      shell_thickness, initial_gap, element_size, cushion_thickness,
      cup_thickness, stabilizer_k (N/mm)
    """
    import trimesh
    rep = {"status": "PROVISIONAL parametric socket (not verified design geometry)"}
    tm = trimesh.Trimesh(skin.vertices, skin.faces, process=False)
    n_p = int(cfg.get("n_petals", 6))
    centers = cfg.get("petal_centers_deg") or list((90.0 + np.arange(n_p) * 360.0 / n_p) % 360.0)
    gap = float(cfg.get("gap_deg", 8.0))
    half = 180.0 / n_p - gap / 2.0
    zb, zt = float(cfg["z_base"]), float(cfg["z_top"])
    post_c = float(cfg.get("posterior_center_deg", 270.0))
    drop = float(cfg.get("posterior_clearance_drop", 25.0))
    post_hw = float(cfg.get("posterior_clearance_halfwidth_deg", 50.0))
    tip = float(cfg.get("tip_drop", 12.0))
    tp, ts, g0 = float(cfg.get("pad_thickness", 3.0)), float(cfg.get("shell_thickness", 4.0)), float(cfg.get("initial_gap", 0.0))
    h = float(cfg.get("element_size", 5.0))

    petals: list[Petal] = []
    pad_blocks, shell_blocks = [], []
    for i, thc in enumerate(centers):
        dpost = abs(((thc - post_c) + 180) % 360 - 180)
        z_top_c = zt - (drop if dpost <= post_hw else 0.0)
        width = np.deg2rad(2 * half)
        r_mid = np.nanmedian(_ray_radius(tm, 0.5 * (zb + z_top_c), np.array([np.deg2rad(thc)])))
        nt = max(3, int(np.ceil(width * r_mid / h)) + 1)
        nz = max(3, int(np.ceil((z_top_c - zb) / h)) + 1)
        th = np.deg2rad(thc) + np.linspace(-width / 2, width / 2, nt)
        u = np.linspace(-1, 1, nt)
        z_top_th = z_top_c - tip * (1 - np.cos(np.pi * u / 2) ** 2) if tip > 0 else np.full(nt, z_top_c)
        s = np.linspace(0, 1, nz)
        Z = zb + s[None, :] * (z_top_th[:, None] - zb)       # (nt, nz)
        R = np.empty_like(Z)
        for j in range(nz):
            # rays at varying z per column
            for k in range(nt):
                R[k, j] = _ray_radius(tm, Z[k, j], th[k:k + 1])[0]
        if np.isnan(R).any():
            raise RuntimeError(f"petal {i}: radial rays missed the skin; check limb axis alignment")
        er = np.stack([np.cos(th), np.sin(th), np.zeros(nt)], axis=1)[:, None, :]
        base = np.concatenate([np.zeros((nt, nz, 2)), Z[..., None]], axis=2)
        P_in = base + (R + g0)[..., None] * er
        P_mid = base + (R + g0 + tp)[..., None] * er
        P_out = base + (R + g0 + tp + ts)[..., None] * er
        ids_in = model.add_nodes(P_in.reshape(-1, 3)).reshape(nt, nz)
        ids_mid = model.add_nodes(P_mid.reshape(-1, 3)).reshape(nt, nz)
        ids_out = model.add_nodes(P_out.reshape(-1, 3)).reshape(nt, nz)

        def hexes(bot, top):
            H = []
            for a in range(nt - 1):
                for b in range(nz - 1):
                    q = [bot[a, b], bot[a + 1, b], bot[a + 1, b + 1], bot[a, b + 1]]
                    r = [top[a, b], top[a + 1, b], top[a + 1, b + 1], top[a, b + 1]]
                    H.append(q + r)
            return np.array(H)

        pad = hexes(ids_in, ids_mid)
        shell = hexes(ids_mid, ids_out)
        name = f"petal{i + 1}"
        model.blocks.append(ElementBlock(f"{name}_pad", "hex8", pad, "silicone_pad", f"socket:{name}:silicone pad", "three-field"))
        model.blocks.append(ElementBlock(f"{name}_shell", "hex8", shell, f"rigid_{name}", f"socket:{name}:rigid petal", "rigid"))
        # contact face of pad (normal pointing toward the skin)
        quads = []
        for a in range(nt - 1):
            for b in range(nz - 1):
                quads.append([ids_in[a, b], ids_in[a, b + 1], ids_in[a + 1, b + 1], ids_in[a + 1, b]])
        model.surfaces[f"{name}_pad_contact"] = Surface(f"{name}_pad_contact", "quad4", np.array(quads))
        hinge_r = float(np.mean(R[:, 0])) + g0 + tp + ts / 2
        thr = np.deg2rad(thc)
        hinge_o = np.array([hinge_r * np.cos(thr), hinge_r * np.sin(thr), zb])
        hinge_ax = np.array([-np.sin(thr), np.cos(thr), 0.0])
        X = model.nodes
        shell_nodes = np.unique(shell.ravel())
        com = X[shell_nodes].mean(axis=0)
        p = Petal(i, name, thc, half, zb, z_top_c, ids_in, ids_mid, ids_out, hinge_o, hinge_ax, com)
        petals.append(p)
        rep[name] = {"theta_center_deg": thc, "angular_width_deg": 2 * half, "z_base": zb, "z_top_center": z_top_c,
                     "posterior_clearance_applied": bool(dpost <= post_hw), "grid": [int(nt), int(nz)],
                     "pad_elements": int(len(pad)), "pad_contact_area_mm2": float(_quad_area(X, np.array(quads)).sum())}

    # ---- distal cup + cushion (offset of the distal skin region) --------
    zc = float(cfg.get("z_cup_top", zb))
    tc, tcup = float(cfg.get("cushion_thickness", 6.0)), float(cfg.get("cup_thickness", 4.0))
    V, F = skin.vertices, skin.faces
    fz = V[F][:, :, 2].max(axis=1)
    cupf = fz <= zc
    sub = F[cupf]
    vids = np.unique(sub.ravel())
    remap = -np.ones(len(V), int); remap[vids] = np.arange(len(vids))
    loc = remap[sub]
    nrm = skin.vertex_normals(smooth_iters=3)[vids]
    Vc = V[vids]
    c_in = model.add_nodes(Vc + g0 * nrm)
    c_mid = model.add_nodes(Vc + (g0 + tc) * nrm)
    c_out = model.add_nodes(Vc + (g0 + tc + tcup) * nrm)
    cushion = np.c_[c_in[loc], c_mid[loc]]
    cup = np.c_[c_mid[loc], c_out[loc]]
    model.blocks.append(ElementBlock("distal_cushion", "penta6", cushion, "distal_cushion", "socket:distal cushion", "three-field"))
    model.blocks.append(ElementBlock("distal_cup", "penta6", cup, "rigid_cup", "socket:rigid distal cup (+ puck/plate/adapter as support)", "rigid"))
    model.surfaces["cushion_contact"] = Surface("cushion_contact", "tri3", c_in[loc][:, [0, 2, 1]])
    rep["distal_cup"] = {"z_top": zc, "faces": int(cupf.sum()), "cushion_thickness": tc,
                         "note": "Anclote puck, plate and Bulldog adapter are not meshed; they are the rigid support of the cup."}

    # ---- skin contact zones (disjoint, primary surfaces) ----------------
    skin_faces = model.surfaces["skin_outer"].conn
    Xs = model.nodes
    cen = Xs[skin_faces].mean(axis=1)
    ang = np.rad2deg(np.arctan2(cen[:, 1], cen[:, 0])) % 360
    zone = np.full(len(skin_faces), -1)
    zone[cen[:, 2] <= zc + 0.5 * h] = 0
    margin = float(cfg.get("contact_zone_margin", 10.0))
    ctr = np.array(centers)
    dang = np.abs(((ang[:, None] - ctr[None, :]) + 180) % 360 - 180)
    nearest = dang.argmin(axis=1)
    ztop_each = np.array([p.z_top_center for p in petals])
    inz = (cen[:, 2] > zc + 0.5 * h) & (cen[:, 2] <= ztop_each[nearest] + margin)
    zone[inz] = nearest[inz] + 1
    model.surfaces["skin_cushion_zone"] = Surface("skin_cushion_zone", "tri3", skin_faces[zone == 0])
    for p in petals:
        model.surfaces[f"skin_{p.name}_zone"] = Surface(f"skin_{p.name}_zone", "tri3", skin_faces[zone == p.index + 1])
    model.surfaces["skin_free"] = Surface("skin_free", "tri3", skin_faces[zone < 0])
    model.meta["petals"] = petals
    model.meta["socket_cfg"] = dict(cfg)
    return model, rep


def _quad_area(X, Q):
    a = X[Q[:, 0]]; b = X[Q[:, 1]]; c = X[Q[:, 2]]; d = X[Q[:, 3]]
    return 0.5 * np.linalg.norm(np.cross(c - a, d - b), axis=1)
