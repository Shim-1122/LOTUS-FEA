"""Conforming layered full-limb volume mesh.

Construction (all layer interfaces share nodes, so the mesh is conforming):

1. The closed skin envelope S0 (outward oriented) is offset inward along
   smoothed vertex normals to S1 (skin/fat interface) and S2 (fat/muscle
   interface).  Skin and fat are linear wedges (penta6) extruded between
   S0-S1 and S1-S2.
2. The muscle volume between S2 and the bone surfaces is tetrahedralised by
   gmsh with the S2 and bone surface meshes kept fixed.
3. Bone interiors are tetrahedralised (used as rigid bodies or deformable
   bone, depending on the configuration).

Layer thickness is an ASSUMPTION unless segmented layers are supplied.
Where bone lies close to the skin the nominal thicknesses are reduced to a
fraction of the local skin-to-bone distance so that layers never overlap
bone; every clamped vertex is reported.  Wedges generated under the proximal
cut (an artificial truncation of the limb, not skin) are assigned to muscle.
"""
from __future__ import annotations

import numpy as np

from ..geometry.surface import TriSurface, closest_distance, self_intersections, surfaces_intersect
from . import quality as Q
from .model import ElementBlock, FEModel, Surface


def _laplacian_smooth_scalar(surf: TriSurface, f: np.ndarray, iters: int) -> np.ndarray:
    nbr = surf.vertex_neighbors()
    for _ in range(iters):
        acc = np.zeros_like(f)
        cnt = np.zeros(len(f))
        np.add.at(acc, nbr[:, 0], f[nbr[:, 1]])
        np.add.at(cnt, nbr[:, 0], 1)
        f = 0.5 * f + 0.5 * acc / np.maximum(cnt, 1)
    return f


def detect_cut_faces(surf: TriSurface, axis: np.ndarray, angle_deg: float = 20.0, plane_tol: float = 1.0):
    """Faces of a flat proximal truncation: normal within angle_deg of +axis
    and centroid within plane_tol (mm) of the most proximal extent."""
    n = surf.face_normals()
    c = surf.vertices[surf.faces].mean(axis=1)
    h = c @ axis
    hmax = (surf.vertices @ axis).max()
    return (n @ axis > np.cos(np.deg2rad(angle_deg))) & (h > hmax - plane_tol)


def tangent_frames(normals: np.ndarray, axis: np.ndarray):
    """Per-face (e_long, e_circ): projection of the limb axis onto the tangent
    plane and the circumferential direction n x e_long."""
    el = axis[None, :] - (normals @ axis)[:, None] * normals
    ln = np.linalg.norm(el, axis=1)
    bad = ln < 0.2  # near the distal apex / proximal cut: axis ~ normal
    alt = np.array([1.0, 0.0, 0.0]) if abs(axis[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    el2 = alt[None, :] - (normals @ alt)[:, None] * normals
    el[bad] = el2[bad]
    el /= np.linalg.norm(el, axis=1)[:, None]
    ec = np.cross(normals, el)
    ec /= np.linalg.norm(ec, axis=1)[:, None]
    return el, ec, bad


def mesh_volume_gmsh(outer: tuple[np.ndarray, np.ndarray, np.ndarray],
                     holes: list[tuple[np.ndarray, np.ndarray, np.ndarray]],
                     mesh_holes: bool, optimize: bool = True, algorithm3d: int = 1,
                     size_fn=None):
    """Tetrahedralise the region between a closed outer surface and closed
    inner surfaces whose triangulations are kept.

    Each surface is (node_tags (n,), coords (n,3), tri (m,3) as node tags).
    Returns (new_node_tags, new_coords, tets_region (k,4) tags, [tets_hole_i]).
    """
    import gmsh
    init_here = not gmsh.isInitialized()
    if init_here:
        gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.clear()
    gmsh.model.add("volume")
    stags = []
    for tags, X, tri in [outer] + holes:
        s = gmsh.model.addDiscreteEntity(2)
        gmsh.model.mesh.addNodes(2, s, tags, X.ravel())
        gmsh.model.mesh.addElementsByType(s, 2, [], tri.ravel())
        stags.append(s)
    # volumes bounded by the fixed discrete surface meshes (geo kernel loops)
    loops = [gmsh.model.geo.addSurfaceLoop([s]) for s in stags]
    vol_main = gmsh.model.geo.addVolume(loops)
    vol_holes = []
    if mesh_holes:
        for lp in loops[1:]:
            vol_holes.append(gmsh.model.geo.addVolume([lp]))
    gmsh.model.geo.synchronize()
    gmsh.option.setNumber("Mesh.Algorithm3D", algorithm3d)
    gmsh.option.setNumber("Mesh.Optimize", 1 if optimize else 0)
    gmsh.option.setNumber("Mesh.OptimizeNetgen", 1 if optimize else 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 1)
    gmsh.option.setNumber("Mesh.Renumber", 0)  # keep our boundary node tags
    if size_fn is not None:
        gmsh.model.mesh.setSizeCallback(lambda dim, t, x, y, z, lc: min(lc, float(size_fn(x, y, z))))
    gmsh.model.mesh.generate(3)
    nt, xyz, _ = gmsh.model.mesh.getNodes(-1, -1, includeBoundary=False, returnParametricCoord=False)
    xyz = xyz.reshape(-1, 3)

    def tets(v):
        et, _, conn = gmsh.model.mesh.getElements(3, v)
        out = [np.asarray(c, np.int64).reshape(-1, 4) for t, c in zip(et, conn) if t == 4]
        return np.vstack(out) if out else np.zeros((0, 4), np.int64)

    main = tets(vol_main)
    hol = [tets(v) for v in vol_holes]
    gmsh.clear()
    if init_here:
        gmsh.finalize()
    return np.asarray(nt, np.int64), xyz, main, hol


def build_layered_model(skin: TriSurface, bones: dict[str, TriSurface], cfg: dict) -> tuple[FEModel, dict]:
    """Create the conforming skin/fat/muscle/bone model.

    cfg keys (mm unless noted): skin_thickness, fat_thickness,
    max_fraction_skin, max_fraction_fat, min_layer, normal_smoothing,
    thickness_smoothing, axis (3,), cut_angle_deg, mesh_bones (bool),
    skin_fiber_reference ('circumferential'|'longitudinal'),
    muscle_fiber (bool), optimize (bool).
    """
    rep: dict = {"assumption": "Layer thicknesses are configurable geometric assumptions, "
                               "NOT patient-specific segmentation."}
    axis = np.asarray(cfg.get("axis", [0, 0, 1.0]), float)
    axis /= np.linalg.norm(axis)
    V0 = skin.vertices.copy()
    F = skin.faces.copy()
    nV = len(V0)
    n = skin.vertex_normals(smooth_iters=int(cfg.get("normal_smoothing", 3)))
    cut = detect_cut_faces(skin, axis, cfg.get("cut_angle_deg", 20.0)) if cfg.get("proximal_cut", True) else np.zeros(len(F), bool)
    rep["proximal_cut_faces"] = int(cut.sum())
    # at cut-rim vertices keep the offset inside the cut plane
    cut_v = np.zeros(nV, bool); cut_v[F[cut].ravel()] = True
    skin_v = np.zeros(nV, bool); skin_v[F[~cut].ravel()] = True
    rim = cut_v & skin_v
    if rim.any():
        nr = n[rim] - (n[rim] @ axis)[:, None] * axis
        n[rim] = nr / np.linalg.norm(nr, axis=1)[:, None]

    # available depth to bone
    if bones:
        d = np.min([closest_distance(b, V0) for b in bones.values() if not b.empty], axis=0)
    else:
        d = np.full(nV, np.inf)
    ts_nom, tf_nom = float(cfg["skin_thickness"]), float(cfg["fat_thickness"])
    fs, ff = float(cfg.get("max_fraction_skin", 0.2)), float(cfg.get("max_fraction_fat", 0.2))
    tmin = float(cfg.get("min_layer", 0.15))
    ts = np.minimum(ts_nom, fs * d)
    tf = np.minimum(tf_nom, ff * d)
    clamped = (ts < ts_nom - 1e-9) | (tf < tf_nom - 1e-9)
    it = int(cfg.get("thickness_smoothing", 5))
    ts = np.minimum(_laplacian_smooth_scalar(skin, ts, it), ts_nom)
    tf = np.minimum(_laplacian_smooth_scalar(skin, tf, it), tf_nom)
    ts = np.minimum(ts, fs * d); tf = np.minimum(tf, ff * d)
    too_thin = (ts < tmin) | (tf < tmin)
    ts = np.maximum(ts, tmin); tf = np.maximum(tf, tmin)
    rep["skin_to_bone_distance_mm"] = Q.summarize(d[np.isfinite(d)]) if np.isfinite(d).any() else None
    rep["thickness_skin_mm"] = Q.summarize(ts)
    rep["thickness_fat_mm"] = Q.summarize(tf)
    rep["vertices_thickness_clamped"] = int(clamped.sum())
    rep["vertices_below_min_layer"] = int(too_thin.sum())
    if too_thin.any():
        rep["warning_min_layer"] = (f"{int(too_thin.sum())} skin vertices are so close to bone that the layer "
                                    f"thickness was raised to min_layer={tmin} mm; bone clearance there is "
                                    f"min {float(d[too_thin].min()):.3f} mm")

    # iterate: offset, check inner-surface validity, locally thin where invalid
    for attempt in range(6):
        V1 = V0 - ts[:, None] * n
        V2 = V1 - tf[:, None] * n
        S2 = TriSurface("S2", V2, F)
        W_skin = np.c_[F + nV, F]          # bottom = S1, top = S0
        W_fat = np.c_[F + 2 * nV, F + nV]  # bottom = S2, top = S1
        Xall = np.vstack([V0, V1, V2])
        jac = np.minimum(Q.penta6_corner_jacobians(Xall, W_skin).min(axis=1),
                         Q.penta6_corner_jacobians(Xall, W_fat).min(axis=1))
        bad_faces = jac <= 0.02
        si = self_intersections(S2)
        bone_hits = sum(surfaces_intersect(S2, b) for b in bones.values() if not b.empty)
        rep.setdefault("offset_attempts", []).append({"attempt": attempt, "bad_wedges": int(bad_faces.sum()),
                                                      "S2_self_intersections": si, "S2_bone_intersections": bone_hits})
        if not bad_faces.any() and si == 0 and bone_hits == 0:
            break
        vbad = np.unique(F[bad_faces].ravel()) if bad_faces.any() else np.zeros(0, int)
        if si or bone_hits:
            # thin everything slightly where curvature is high (fallback, logged)
            vbad = np.arange(nV)
        ts[vbad] *= 0.7; tf[vbad] *= 0.7
        n = TriSurface("tmp", V0, F).vertex_normals(smooth_iters=int(cfg.get("normal_smoothing", 3)) + 2 * (attempt + 1))
    else:
        raise RuntimeError(f"layer offset could not be made valid: {rep['offset_attempts'][-1]}")
    rep["min_wedge_scaled_jacobian"] = float(jac.min())

    # ---- muscle + bone tetrahedra -------------------------------------
    Xs = np.vstack([V0, V1, V2])
    tag_S2 = np.arange(2 * nV, 3 * nV) + 1
    holes, hole_names, off = [], [], 3 * nV
    for name, b in bones.items():
        if b.empty:
            continue
        tags = np.arange(off, off + len(b.vertices)) + 1
        holes.append((tags, b.vertices, tags[b.faces]))
        hole_names.append(name)
        Xs = np.vstack([Xs, b.vertices])
        off += len(b.vertices)
    new_tags, new_xyz, tet_main, tet_holes = mesh_volume_gmsh(
        (tag_S2, V2, tag_S2[F]), holes, mesh_holes=bool(cfg.get("mesh_bones", True)),
        optimize=bool(cfg.get("optimize", True)), algorithm3d=int(cfg.get("algorithm3d", 1)))
    # map gmsh tags -> global indices
    maxtag = int(max(off, new_tags.max() if len(new_tags) else 0)) + 1
    tag2idx = np.full(maxtag + 1, -1, np.int64)
    tag2idx[1:off + 1] = np.arange(off)
    extra = new_tags[new_tags > off]
    ex_xyz = new_xyz[new_tags > off]
    tag2idx[extra] = off + np.arange(len(extra))
    X = np.vstack([Xs, ex_xyz])
    # sanity: gmsh must not have moved or duplicated boundary nodes
    kept = new_tags[new_tags <= off]
    if len(kept):
        drift = np.abs(new_xyz[new_tags <= off] - Xs[kept - 1]).max()
        rep["boundary_node_drift_mm"] = float(drift)
        if drift > 1e-6:
            raise RuntimeError(f"gmsh moved boundary nodes by {drift} mm")

    def fix_orient(T):
        T = tag2idx[T]
        if (T < 0).any():
            raise RuntimeError("unmapped node tag in tetrahedra")
        v = Q.tet_volume(X, T)
        T[v < 0] = T[v < 0][:, [0, 2, 1, 3]]
        return T

    tet_main = fix_orient(tet_main)
    tet_holes = [fix_orient(t) for t in tet_holes]

    # ---- element blocks -------------------------------------------------
    model = FEModel(nodes=X)
    fn = skin.face_normals()
    el, ec, apex = tangent_frames(fn, axis)
    ref = cfg.get("skin_fiber_reference", "circumferential")
    a = ec if ref == "circumferential" else el
    dvec = el if ref == "circumferential" else ec
    ax = np.stack([a, dvec], axis=1)
    sk = ~cut
    model.blocks.append(ElementBlock("skin", "penta6", W_skin[sk], "skin", "skin", "three-field", mat_axis=ax[sk],
                                     meta={"fiber_reference_axis": ref, "apex_fallback_faces": int(apex[sk].sum())}))
    model.blocks.append(ElementBlock("fat", "penta6", W_fat[sk], "fat", "fat", "three-field"))
    if cut.any():
        model.blocks.append(ElementBlock("proximal_cap", "penta6", np.vstack([W_skin[cut], W_fat[cut]]), "muscle",
                                         "muscle (layers under artificial proximal cut)", "three-field"))
    mform = cfg.get("muscle_formulation", "ut4")
    mfib = None
    if cfg.get("muscle_fiber", False):
        mfib = np.repeat(axis[None, :], len(tet_main), axis=0)
    model.blocks.append(ElementBlock("muscle", "tet4", tet_main, "muscle", "muscle", mform, fiber=mfib))
    for name, t in zip(hole_names, tet_holes):
        model.blocks.append(ElementBlock(f"bone_{name}", "tet4", t, f"bone_{name}", f"bone:{name}", "bone"))
    # surfaces / node sets
    model.surfaces["skin_outer"] = Surface("skin_outer", "tri3", F[~cut])
    if cut.any():
        model.surfaces["proximal_cut"] = Surface("proximal_cut", "tri3", F[cut])
        model.nodesets["proximal_cut"] = np.unique(F[cut].ravel())
    for name, (tags, Xb, tri) in zip(hole_names, holes):
        model.surfaces[f"bone_{name}_surface"] = Surface(f"bone_{name}_surface", "tri3", tag2idx[tri])
    model.meta.update({"axis": axis.tolist(), "n_skin_vertices": nV, "layer_node_offsets": [0, nV, 2 * nV],
                       "skin_thickness_per_vertex": ts, "fat_thickness_per_vertex": tf,
                       "skin_to_bone_distance": d, "bone_names": hole_names})
    # ---- quality ----------------------------------------------------------
    tq = Q.tet_quality(X, tet_main)
    rep["muscle_tets"] = int(len(tet_main))
    rep["muscle_tet_quality"] = {"radius_ratio": Q.summarize(tq["radius_ratio"]),
                                 "scaled_jacobian": Q.summarize(tq["scaled_jacobian"]),
                                 "inverted": int((tq["volume"] <= 0).sum())}
    for name, t in zip(hole_names, tet_holes):
        if len(t):
            q = Q.tet_quality(X, t)
            rep[f"bone_{name}_tets"] = int(len(t))
            rep[f"bone_{name}_inverted"] = int((q["volume"] <= 0).sum())
    vol_parts = {b.name: float(_block_volume(X, b)) for b in model.blocks}
    rep["block_volumes_mm3"] = vol_parts
    rep["model_summary"] = model.summary()
    return model, rep


def _block_volume(X, b: ElementBlock) -> float:
    if b.etype == "tet4":
        return Q.tet_volume(X, b.conn).sum()
    if b.etype == "penta6":
        # split wedge into 3 tets
        c = b.conn
        t = np.vstack([c[:, [0, 1, 2, 3]], c[:, [1, 2, 3, 4]], c[:, [2, 3, 4, 5]]])
        return np.abs(Q.tet_volume(X, t)).sum()
    if b.etype == "hex8":
        c = b.conn
        tets = [(0, 1, 3, 4), (1, 2, 3, 6), (1, 3, 4, 6), (3, 4, 6, 7), (1, 4, 5, 6)]
        return sum(np.abs(Q.tet_volume(X, c[:, list(t)])).sum() for t in tets)
    return 0.0
