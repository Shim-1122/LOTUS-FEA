"""Triangulated surface container, import and audit.

Source files are never modified.  Every body carries a ``history`` list that
records how it was obtained (reader, unit conversion, transforms, repairs)
together with before/after audit summaries.
"""
from __future__ import annotations

import copy
import hashlib
import os
from dataclasses import dataclass, field

import numpy as np

from . import step_tessellated as stt

SURFACE_EXT = {".stl", ".obj", ".ply", ".off"}
STEP_EXT = {".step", ".stp"}


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


@dataclass
class TriSurface:
    name: str
    vertices: np.ndarray           # (n,3) float, mm
    faces: np.ndarray              # (m,3) int
    source: str = ""
    source_sha256: str = ""
    source_unit: str = ""
    history: list = field(default_factory=list)

    @property
    def empty(self) -> bool:
        return len(self.faces) == 0

    def copy(self) -> "TriSurface":
        return copy.deepcopy(self)

    def log(self, op: str, **info):
        self.history.append({"op": op, **info})

    # ---- basic geometry -------------------------------------------------
    def face_normals(self, unit=True):
        v = self.vertices
        f = self.faces
        n = np.cross(v[f[:, 1]] - v[f[:, 0]], v[f[:, 2]] - v[f[:, 0]])
        if unit:
            ln = np.linalg.norm(n, axis=1)
            n = n / np.where(ln > 0, ln, 1.0)[:, None]
        return n

    def areas(self):
        return 0.5 * np.linalg.norm(self.face_normals(unit=False), axis=1)

    def volume(self) -> float:
        v = self.vertices
        f = self.faces
        return float(np.einsum("ij,ij->i", v[f[:, 0]], np.cross(v[f[:, 1]], v[f[:, 2]])).sum() / 6.0)

    def vertex_normals(self, smooth_iters: int = 0) -> np.ndarray:
        fn = self.face_normals(unit=False)  # area weighted
        vn = np.zeros_like(self.vertices)
        for k in range(3):
            np.add.at(vn, self.faces[:, k], fn)
        vn /= np.linalg.norm(vn, axis=1)[:, None] + 1e-30
        if smooth_iters:
            nbr = self.vertex_neighbors()
            for _ in range(smooth_iters):
                acc = np.zeros_like(vn)
                cnt = np.zeros(len(vn))
                np.add.at(acc, nbr[:, 0], vn[nbr[:, 1]])
                np.add.at(cnt, nbr[:, 0], 1)
                vn = 0.5 * vn + 0.5 * acc / cnt[:, None]
                vn /= np.linalg.norm(vn, axis=1)[:, None]
        return vn

    def edges(self):
        f = self.faces
        e = np.vstack([f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]])
        return e

    def vertex_neighbors(self):
        e = self.edges()
        e = np.vstack([e, e[:, ::-1]])
        return np.unique(e, axis=0)

    def to_trimesh(self):
        import trimesh
        return trimesh.Trimesh(self.vertices, self.faces, process=False)


# --------------------------------------------------------------------------
# Import
# --------------------------------------------------------------------------

def _gmsh_step_brep(path: str, mesh_size: float):
    """Import a B-rep STEP with OpenCASCADE (via gmsh) and tessellate it.

    gmsh is asked to convert to millimetres (Geometry.OCCTargetUnit = MM) so
    the unit conversion is performed by OCC from the STEP unit context.
    Returns a list of (name, V, F, info).
    """
    import gmsh
    init_here = not gmsh.isInitialized()
    if init_here:
        gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.clear()
    gmsh.option.setString("Geometry.OCCTargetUnit", "MM")
    gmsh.model.occ.importShapes(path)
    gmsh.model.occ.synchronize()
    vols = gmsh.model.getEntities(3)
    out = []
    info = {"occ_volumes": len(vols), "occ_surfaces": len(gmsh.model.getEntities(2)),
            "tessellation_mesh_size_mm": mesh_size}
    if vols or gmsh.model.getEntities(2):
        gmsh.option.setNumber("Mesh.MeshSizeMax", mesh_size)
        gmsh.option.setNumber("Mesh.MeshSizeMin", mesh_size * 0.25)
        gmsh.option.setNumber("Mesh.Algorithm", 6)
        gmsh.model.mesh.generate(2)
        nt, xyz, _ = gmsh.model.mesh.getNodes()
        xyz = xyz.reshape(-1, 3)
        idx = {t: i for i, t in enumerate(nt)}
        F = []
        for d, s in gmsh.model.getEntities(2):
            et, _, conn = gmsh.model.mesh.getElements(2, s)
            for t, c in zip(et, conn):
                if t == 2:
                    F.append(np.array([idx[x] for x in c]).reshape(-1, 3))
        F = np.vstack(F) if F else np.zeros((0, 3), int)
        info["brep_volume_mm3"] = float(sum(gmsh.model.occ.getMass(3, v) for _, v in vols)) if vols else None
        out.append((os.path.splitext(os.path.basename(path))[0], xyz, F, info))
    else:
        out.append((os.path.splitext(os.path.basename(path))[0], np.zeros((0, 3)), np.zeros((0, 3), int), info))
    gmsh.clear()
    if init_here:
        gmsh.finalize()
    return out


def load_body(path: str, unit: str | None = None, name: str | None = None,
              brep_mesh_size: float = 2.0) -> list[TriSurface]:
    """Load one geometry file into one or more TriSurface bodies (mm).

    * STEP AP242 tessellated shells: parsed directly; unit from STEP context.
    * STEP B-rep: OpenCASCADE via gmsh, converted to mm by OCC.
    * STL/OBJ/PLY/OFF: unit must be given explicitly (files carry no units).
    """
    from ..units import LENGTH
    ext = os.path.splitext(path)[1].lower()
    digest = sha256_file(path)
    bodies: list[TriSurface] = []
    if ext in STEP_EXT:
        if stt.is_tessellated_step(path):
            for b in stt.read_tessellated_step(path):
                s = TriSurface(name or b.name, b.vertices * b.scale_to_mm, b.faces.copy(), path, digest, b.unit_name)
                s.log("import", reader="rlfea.step_tessellated", step_unit=b.unit_name,
                      scale_to_mm=b.scale_to_mm, face_entities=b.n_face_entities)
                bodies.append(s)
        elif stt.has_brep_solid(path):
            for nm, V, F, info in _gmsh_step_brep(path, brep_mesh_size):
                s = TriSurface(name or nm, V, F, path, digest, "from STEP context (OCC->mm)")
                s.log("import", reader="gmsh/OpenCASCADE", **info)
                bodies.append(s)
        else:
            s = TriSurface(name or os.path.splitext(os.path.basename(path))[0], np.zeros((0, 3)),
                           np.zeros((0, 3), int), path, digest, "n/a")
            s.log("import", reader="none", note="STEP file contains no B-rep or tessellated geometry (EMPTY BODY)")
            bodies.append(s)
    elif ext in SURFACE_EXT:
        if unit is None:
            raise ValueError(f"{path}: surface-mesh files carry no units; specify 'unit' explicitly")
        import trimesh
        m = trimesh.load(path, force="mesh", process=False)
        scale = LENGTH[unit]
        s = TriSurface(name or os.path.splitext(os.path.basename(path))[0],
                       np.asarray(m.vertices, float) * scale, np.asarray(m.faces, int), path, digest, unit)
        s.log("import", reader="trimesh", declared_unit=unit, scale_to_mm=scale)
        bodies.append(s)
    else:
        raise ValueError(f"unsupported geometry format {ext}")
    return bodies


# --------------------------------------------------------------------------
# Geometric queries
# --------------------------------------------------------------------------

def winding_number(surf: TriSurface, pts: np.ndarray, chunk: int | None = None) -> np.ndarray:
    """Generalised winding number (solid-angle sum / 4pi), Van Oosterom &
    Strackee formula.  ~1 inside a closed, outward-oriented surface, ~0
    outside; robust to small defects."""
    V = surf.vertices
    F = surf.faces
    A = V[F[:, 0]]
    B = V[F[:, 1]]
    C = V[F[:, 2]]
    out = np.empty(len(pts))
    if chunk is None:  # bound temporary memory to ~150 MB
        chunk = max(1, int(1_500_000 // max(len(F), 1)))
    for i in range(0, len(pts), chunk):
        p = pts[i:i + chunk, None, :]
        a = A[None] - p
        b = B[None] - p
        c = C[None] - p
        la = np.linalg.norm(a, axis=2)
        lb = np.linalg.norm(b, axis=2)
        lc = np.linalg.norm(c, axis=2)
        num = np.einsum("ijk,ijk->ij", a, np.cross(b, c))
        den = (la * lb * lc + np.einsum("ijk,ijk->ij", a, b) * lc
               + np.einsum("ijk,ijk->ij", b, c) * la + np.einsum("ijk,ijk->ij", c, a) * lb)
        out[i:i + chunk] = 2.0 * np.arctan2(num, den).sum(axis=1) / (4 * np.pi)
    return out


def closest_distance(surf: TriSurface, pts: np.ndarray) -> np.ndarray:
    """Exact unsigned distance from points to a triangle surface."""
    import trimesh
    m = trimesh.Trimesh(surf.vertices, surf.faces, process=False)
    _, d, _ = trimesh.proximity.closest_point(m, pts)
    return np.asarray(d)


def _tri_tri_intersect(P, Q):
    """Vectorised Moller-style triangle/triangle intersection test.
    P, Q: (k,3,3).  Returns bool (k,).  Coplanar pairs are reported False
    (they are handled by the duplicate-face check)."""
    def plane(T):
        n = np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0])
        d = -np.einsum("ij,ij->i", n, T[:, 0])
        return n, d

    nq, dq = plane(Q)
    dp = np.einsum("ij,ikj->ik", nq, P) + dq[:, None]
    np_, dp_ = plane(P)
    dqq = np.einsum("ij,ikj->ik", np_, Q) + dp_[:, None]
    eps = 1e-12 * (np.abs(dp).max(axis=1) + 1)
    sep = (np.all(dp > eps[:, None], axis=1) | np.all(dp < -eps[:, None], axis=1) |
           np.all(dqq > eps[:, None], axis=1) | np.all(dqq < -eps[:, None], axis=1))
    res = np.zeros(len(P), bool)
    cand = ~sep
    if not cand.any():
        return res
    # interval overlap on intersection line
    D = np.cross(np_[cand], nq[cand])
    Pc, Qc = P[cand], Q[cand]
    dpc, dqc = dp[cand], dqq[cand]

    def interval(T, dist):
        proj = np.einsum("ij,ikj->ik", D, T)
        lo = np.full(len(T), np.inf)
        hi = np.full(len(T), -np.inf)
        for a, b in ((0, 1), (1, 2), (2, 0)):
            da, db = dist[:, a], dist[:, b]
            crosses = (da * db) <= 0
            denom = np.where(np.abs(da - db) > 1e-30, da - db, 1e-30)
            t = da / denom
            x = proj[:, a] + (proj[:, b] - proj[:, a]) * t
            x = np.where(crosses, x, np.nan)
            lo = np.fmin(lo, x)
            hi = np.fmax(hi, x)
        return lo, hi

    l1, h1 = interval(Pc, dpc)
    l2, h2 = interval(Qc, dqc)
    coplanar = np.linalg.norm(D, axis=1) < 1e-12 * (np.linalg.norm(np_[cand], axis=1) * np.linalg.norm(nq[cand], axis=1) + 1e-300)
    tol = 1e-9 * (np.abs(np.r_[l1, h1, l2, h2][np.isfinite(np.r_[l1, h1, l2, h2])]).max() + 1) if np.isfinite(np.r_[l1, h1]).any() else 0
    ov = (np.minimum(h1, h2) - np.maximum(l1, l2)) > tol
    res[np.where(cand)[0]] = ov & ~coplanar
    return res


def self_intersections(surf: TriSurface, max_faces: int = 400_000) -> int:
    """Count intersecting pairs of non-adjacent triangles (AABB prefilter)."""
    if surf.empty or len(surf.faces) > max_faces:
        return -1
    from rtree import index
    V, F = surf.vertices, surf.faces
    T = V[F]
    lo, hi = T.min(axis=1), T.max(axis=1)
    pad = 1e-9 * (np.ptp(V, axis=0).max() + 1)
    p = index.Property(); p.dimension = 3
    idx = index.Index(((i, (*lo[i] - pad, *hi[i] + pad), None) for i in range(len(F))), properties=p)
    pairs = []
    for i in range(len(F)):
        for j in idx.intersection((*lo[i] - pad, *hi[i] + pad)):
            if j > i:
                pairs.append((i, j))
    if not pairs:
        return 0
    pairs = np.array(pairs)
    shared = (F[pairs[:, 0]][:, :, None] == F[pairs[:, 1]][:, None, :]).any(axis=2).any(axis=1)
    pairs = pairs[~shared]
    if len(pairs) == 0:
        return 0
    n = 0
    for k in range(0, len(pairs), 200_000):
        pp = pairs[k:k + 200_000]
        n += int(_tri_tri_intersect(T[pp[:, 0]], T[pp[:, 1]]).sum())
    return n


def surfaces_intersect(a: TriSurface, b: TriSurface) -> int:
    """Number of intersecting triangle pairs between two surfaces."""
    from rtree import index
    Ta, Tb = a.vertices[a.faces], b.vertices[b.faces]
    lo, hi = Tb.min(axis=1), Tb.max(axis=1)
    p = index.Property(); p.dimension = 3
    idx = index.Index(((i, (*lo[i], *hi[i]), None) for i in range(len(Tb))), properties=p)
    la, ha = Ta.min(axis=1), Ta.max(axis=1)
    pairs = [(i, j) for i in range(len(Ta)) for j in idx.intersection((*la[i], *ha[i]))]
    if not pairs:
        return 0
    pairs = np.array(pairs)
    return int(_tri_tri_intersect(Ta[pairs[:, 0]], Tb[pairs[:, 1]]).sum())


# --------------------------------------------------------------------------
# Audit
# --------------------------------------------------------------------------

def audit_body(s: TriSurface, check_self_intersection: bool = True) -> dict:
    rep = {"name": s.name, "source": s.source, "source_sha256": s.source_sha256,
           "source_unit": s.source_unit, "internal_unit": "mm",
           "n_vertices": int(len(s.vertices)), "n_faces": int(len(s.faces)),
           "empty": bool(s.empty), "history_ops": [h["op"] for h in s.history]}
    if s.empty:
        rep["status"] = "EMPTY"
        return rep
    V, F = s.vertices, s.faces
    bb0, bb1 = V.min(0), V.max(0)
    rep["bbox_min_mm"] = bb0.round(4).tolist()
    rep["bbox_max_mm"] = bb1.round(4).tolist()
    rep["dimensions_mm"] = (bb1 - bb0).round(4).tolist()
    c = V.mean(0)
    _, sv, vt = np.linalg.svd(V - c, full_matrices=False)
    rep["vertex_centroid_mm"] = c.round(4).tolist()
    rep["principal_axes"] = vt.round(6).tolist()
    rep["principal_extent_rms_mm"] = (sv / np.sqrt(len(V))).round(4).tolist()
    ar = s.areas()
    rep["area_mm2"] = float(ar.sum())
    rep["signed_volume_mm3"] = s.volume()
    rep["degenerate_faces"] = int((ar < 1e-10 * max(ar.max(), 1e-30)).sum())
    fs = np.sort(F, axis=1)
    rep["duplicate_faces"] = int(len(fs) - len(np.unique(fs, axis=0)))
    e = np.sort(s.edges(), axis=1)
    ue, cnt = np.unique(e, axis=0, return_counts=True)
    rep["boundary_edges"] = int((cnt == 1).sum())
    rep["nonmanifold_edges"] = int((cnt > 2).sum())
    # orientation consistency: each interior directed edge should appear once each way
    de = s.edges()
    _, dcnt = np.unique(de, axis=0, return_counts=True)
    rep["inconsistently_oriented_edges"] = int((dcnt > 1).sum())
    rep["watertight"] = rep["boundary_edges"] == 0 and rep["nonmanifold_edges"] == 0
    rep["euler_characteristic"] = int(len(np.unique(F)) - len(ue) + len(F))
    # connected components (face adjacency via shared vertices)
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    nf = len(F)
    rows = np.repeat(np.arange(nf), 3)
    g = coo_matrix((np.ones(3 * nf), (rows, F.ravel())), shape=(nf, len(V))).tocsr()
    adj = g @ g.T
    ncomp, lab = connected_components(adj, directed=False)
    rep["connected_components"] = int(ncomp)
    rep["outward_oriented"] = bool(rep["signed_volume_mm3"] > 0) if rep["watertight"] else None
    if check_self_intersection:
        rep["self_intersecting_pairs"] = self_intersections(s)
    defects = []
    if rep["degenerate_faces"]:
        defects.append("degenerate faces")
    if rep["duplicate_faces"]:
        defects.append("duplicate faces")
    if not rep["watertight"]:
        defects.append("not watertight")
    if rep["inconsistently_oriented_edges"]:
        defects.append("inconsistent orientation")
    if rep.get("self_intersecting_pairs", 0) and rep["self_intersecting_pairs"] > 0:
        defects.append("self-intersections")
    if rep["outward_oriented"] is False:
        defects.append("inward oriented (negative volume)")
    rep["defects"] = defects
    rep["status"] = "VALID_CLOSED_SURFACE" if not defects and ncomp == 1 else ("VALID_MULTI_COMPONENT" if not defects else "DEFECTIVE")
    return rep


def audit_assembly(envelope: TriSurface | None, inner: list[TriSurface], tol_mm: float = 1e-3) -> dict:
    """Inter-body checks: duplicates, containment in the limb envelope,
    pairwise intersections and minimum clearance to the envelope."""
    out = {"pairs": [], "bodies": {}}
    bodies = ([envelope] if envelope is not None else []) + list(inner)
    # duplicates (geometric fingerprint)
    fp = {}
    for b in bodies:
        if b.empty:
            continue
        key = (round(abs(b.volume()), 1), tuple(np.round(b.vertices.min(0), 2)), tuple(np.round(b.vertices.max(0), 2)))
        fp.setdefault(key, []).append(b.name)
    out["duplicate_groups"] = [v for v in fp.values() if len(v) > 1]
    for b in inner:
        info = {}
        if b.empty:
            info["status"] = "EMPTY"
            out["bodies"][b.name] = info
            continue
        if envelope is not None and not envelope.empty:
            w = winding_number(envelope, b.vertices)
            info["fraction_vertices_inside_envelope"] = float((w > 0.5).mean())
            d = closest_distance(envelope, b.vertices)
            inside = w > 0.5
            info["min_clearance_to_envelope_mm"] = float(d[inside].min()) if inside.any() else None
            info["outside_envelope"] = bool((~inside).any())
            info["intersects_envelope_tri_pairs"] = surfaces_intersect(b, envelope)
        out["bodies"][b.name] = info
    for i in range(len(inner)):
        for j in range(i + 1, len(inner)):
            a, b = inner[i], inner[j]
            if a.empty or b.empty:
                continue
            w_ab = winding_number(b, a.vertices)
            w_ba = winding_number(a, b.vertices)
            d = closest_distance(b, a.vertices).min()
            out["pairs"].append({"a": a.name, "b": b.name,
                                 "a_vertices_inside_b": int((w_ab > 0.5).sum()),
                                 "b_vertices_inside_a": int((w_ba > 0.5).sum()),
                                 "intersecting_tri_pairs": surfaces_intersect(a, b),
                                 "min_distance_mm": float(d)})
    return out
