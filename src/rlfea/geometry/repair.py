"""Logged geometry repairs and transformations.

Every operation returns a *new* TriSurface; the input object is never changed,
and each call appends a history record with before/after audit summaries so a
reviewer can see exactly what was done to the source geometry.
"""
from __future__ import annotations

import numpy as np

from .surface import TriSurface, audit_body, closest_distance

_SUMMARY_KEYS = ("n_vertices", "n_faces", "signed_volume_mm3", "area_mm2", "watertight",
                 "connected_components", "degenerate_faces", "duplicate_faces",
                 "boundary_edges", "nonmanifold_edges", "inconsistently_oriented_edges", "status")


def _summary(s: TriSurface) -> dict:
    a = audit_body(s, check_self_intersection=False)
    return {k: a.get(k) for k in _SUMMARY_KEYS}


def transform(s: TriSurface, matrix: np.ndarray, reason: str) -> TriSurface:
    """Apply a 4x4 homogeneous rigid/scaling transform (logged)."""
    M = np.asarray(matrix, float)
    out = s.copy()
    before = _summary(s)
    out.vertices = (np.c_[s.vertices, np.ones(len(s.vertices))] @ M.T)[:, :3]
    if np.linalg.det(M[:3, :3]) < 0:
        out.faces = out.faces[:, ::-1].copy()
    out.log("transform", matrix=M.round(12).tolist(), reason=reason, before=before, after=_summary(out))
    return out


def basic_cleanup(s: TriSurface, merge_tol_mm: float = 1e-6) -> TriSurface:
    """Merge coincident vertices, drop degenerate/duplicate faces and unused
    vertices, and make face orientation consistent and outward."""
    import trimesh
    before = _summary(s)
    m = trimesh.Trimesh(s.vertices.copy(), s.faces.copy(), process=False)
    nv0, nf0 = len(m.vertices), len(m.faces)
    m.merge_vertices(digits_vertex=int(max(0, -np.log10(merge_tol_mm))))
    m.update_faces(m.nondegenerate_faces())
    m.update_faces(m.unique_faces())
    m.remove_unreferenced_vertices()
    trimesh.repair.fix_winding(m)
    flipped = False
    if m.is_watertight and m.volume < 0:
        m.invert(); flipped = True
    out = s.copy()
    out.vertices = np.asarray(m.vertices, float)
    out.faces = np.asarray(m.faces, int)
    out.log("basic_cleanup", merge_tol_mm=merge_tol_mm, removed_vertices=nv0 - len(m.vertices),
            removed_faces=nf0 - len(m.faces), inverted_to_outward=flipped,
            before=before, after=_summary(out))
    return out


def remesh(s: TriSurface, size_mm: float, size_fn=None, angle_deg: float = 40.0,
           reason: str = "uniform FE surface resolution") -> TriSurface:
    """Isotropic remeshing of a closed triangulated surface with gmsh.

    The surface is classified into patches (dihedral threshold ``angle_deg``),
    reparametrised and remeshed at ``size_mm`` (or ``size_fn(x,y,z)``).  The
    deviation from the input is measured both ways and logged.
    """
    import gmsh
    before = _summary(s)
    init_here = not gmsh.isInitialized()
    if init_here:
        gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.clear()
    gmsh.model.add("remesh")
    tag = gmsh.model.addDiscreteEntity(2)
    nodes = np.arange(1, len(s.vertices) + 1)
    gmsh.model.mesh.addNodes(2, tag, nodes, s.vertices.ravel())
    gmsh.model.mesh.addElementsByType(tag, 2, [], (s.faces + 1).ravel())
    gmsh.model.mesh.classifySurfaces(np.deg2rad(angle_deg), True, True, np.pi)
    gmsh.model.mesh.createGeometry()
    gmsh.option.setNumber("Mesh.MeshSizeMax", size_mm)
    gmsh.option.setNumber("Mesh.MeshSizeMin", size_mm * (0.2 if size_fn is not None else 0.9))
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    if size_fn is not None:
        gmsh.model.mesh.setSizeCallback(lambda dim, t, x, y, z, lc: float(size_fn(x, y, z)))
    gmsh.model.mesh.generate(2)
    nt, xyz, _ = gmsh.model.mesh.getNodes()
    xyz = xyz.reshape(-1, 3)
    idx = np.zeros(int(nt.max()) + 1, int)
    idx[nt.astype(int)] = np.arange(len(nt))
    F = []
    for d, t in gmsh.model.getEntities(2):
        et, _, conn = gmsh.model.mesh.getElements(2, t)
        for typ, c in zip(et, conn):
            if typ == 2:
                F.append(idx[np.asarray(c, int)].reshape(-1, 3))
    gmsh.clear()
    if init_here:
        gmsh.finalize()
    out = s.copy()
    out.vertices = xyz
    out.faces = np.vstack(F)
    out = basic_cleanup(out)
    out.history.pop()  # fold the cleanup into the remesh record
    dev_new_to_old = closest_distance(s, out.vertices)
    rng = np.random.default_rng(0)
    sample = s.vertices[rng.choice(len(s.vertices), min(len(s.vertices), 20000), replace=False)]
    dev_old_to_new = closest_distance(out, sample)
    out.log("remesh", tool="gmsh classifySurfaces/createGeometry + 2D mesh", size_mm=size_mm,
            local_size_function=size_fn is not None, angle_deg=angle_deg, reason=reason,
            max_deviation_new_to_source_mm=float(dev_new_to_old.max()),
            max_deviation_source_to_new_mm=float(dev_old_to_new.max()),
            volume_change_percent=float(100 * (out.volume() - s.volume()) / s.volume()),
            before=before, after=_summary(out))
    return out
