"""Read FEBio raw output and derive solver-consistent fields.

Only quantities computed by FEBio (nodal displacement, element deformation
gradient, Cauchy stress, strain-energy density, contact gap/pressure/traction,
rigid-body and connector data, nodal reactions) are read.  Derived fields are
pointwise algebra on those values (principal values, invariants, pressure).
No smoothing or interpolation across regions is applied: element fields are
exported as cell data on the element where FEBio evaluated them.
"""
from __future__ import annotations

import os

import numpy as np

from ..mesh import quality as Q
from ..mesh.model import FEModel
from ..solvers.febio import read_logdata

VTK_ORDER = {"tet4": None, "hex8": None, "penta6": [0, 2, 1, 3, 5, 4]}  # FEBio -> VTK wedge orientation


def _by_time(recs):
    return {round(r["time"], 9): r for r in recs}


def load_raw(solver_dir: str, model: FEModel, surfaces: list[str]) -> dict:
    nd = _by_time(read_logdata(os.path.join(solver_dir, "out_node_disp.txt")))
    nr = _by_time(read_logdata(os.path.join(solver_dir, "out_node_react.txt")))
    eF = _by_time(read_logdata(os.path.join(solver_dir, "out_elem_F.txt")))
    eS = _by_time(read_logdata(os.path.join(solver_dir, "out_elem_stress.txt")))
    eW = _by_time(read_logdata(os.path.join(solver_dir, "out_elem_sed.txt")))
    rg = _by_time(read_logdata(os.path.join(solver_dir, "out_rigid.txt")))
    rc = _by_time(read_logdata(os.path.join(solver_dir, "out_connectors.txt")))
    fc = {s: _by_time(read_logdata(os.path.join(solver_dir, f"out_face_{s}.txt"))) for s in surfaces}
    times = sorted(nd)
    N, E = len(model.nodes), model.n_elements()
    states = []
    for t in times:
        st = {"time": t}
        U = np.full((N, 3), np.nan)
        U[nd[t]["ids"] - 1] = nd[t]["data"]
        st["U"] = U
        if t in nr:
            Rf = np.zeros((N, 3)); Rf[nr[t]["ids"] - 1] = nr[t]["data"]; st["R"] = Rf
        if t in eF:
            F = np.full((E, 9), np.nan); F[eF[t]["ids"] - 1] = eF[t]["data"]; st["F"] = F.reshape(E, 3, 3)
        if t in eS:
            S = np.full((E, 6), np.nan); S[eS[t]["ids"] - 1] = eS[t]["data"]; st["S"] = S
        if t in eW:
            W = np.full(E, np.nan); W[eW[t]["ids"] - 1] = eW[t]["data"][:, 0]; st["sed"] = W
        if t in rg:
            st["rigid"] = {"ids": rg[t]["ids"], "data": rg[t]["data"]}
        if t in rc:
            st["connectors"] = {"ids": rc[t]["ids"], "data": rc[t]["data"]}
        st["faces"] = {}
        for s, recs in fc.items():
            if t in recs:
                nf = len(model.surfaces[s].conn)
                D = np.full((nf, 5), np.nan); D[recs[t]["ids"] - 1] = recs[t]["data"]
                st["faces"][s] = D
        states.append(st)
    return {"times": times, "states": states}


def element_fields(F: np.ndarray, S6: np.ndarray) -> dict:
    """Per-element derived fields from FEBio's F and Cauchy stress."""
    J = np.linalg.det(F)
    C = np.einsum("eki,ekj->eij", F, F)
    lam2 = np.linalg.eigvalsh(C)                     # ascending
    lam2 = np.clip(lam2, 1e-300, None)
    GL = 0.5 * (lam2 - 1.0)                          # principal Green-Lagrange
    H = 0.5 * np.log(lam2)                           # principal Hencky (log) strain
    sig = np.zeros((len(S6), 3, 3))
    sig[:, 0, 0], sig[:, 1, 1], sig[:, 2, 2] = S6[:, 0], S6[:, 1], S6[:, 2]
    sig[:, 0, 1] = sig[:, 1, 0] = S6[:, 3]
    sig[:, 1, 2] = sig[:, 2, 1] = S6[:, 4]
    sig[:, 0, 2] = sig[:, 2, 0] = S6[:, 5]
    sp = np.linalg.eigvalsh(sig)
    p = -np.trace(sig, axis1=1, axis2=2) / 3.0
    vm = np.sqrt(0.5 * ((sp[:, 0] - sp[:, 1]) ** 2 + (sp[:, 1] - sp[:, 2]) ** 2 + (sp[:, 2] - sp[:, 0]) ** 2))
    return {"J": J, "GL_max": GL[:, 2], "GL_min": GL[:, 0], "logstrain_max": H[:, 2], "logstrain_min": H[:, 0],
            "stretch_max": np.sqrt(lam2[:, 2]), "stretch_min": np.sqrt(lam2[:, 0]),
            "sigma_1": sp[:, 2], "sigma_3": sp[:, 0], "pressure": p, "von_mises": vm,
            "shear_strain_max": 0.5 * (H[:, 2] - H[:, 0])}


def contact_fields(X0, U, tri, D, mu):
    """Normal/tangential decomposition on the DEFORMED primary facets.
    Status: 0 open (p<=0), 1 stick (|t_t| < 0.98 mu p), 2 slip (|t_t| >= 0.98 mu p).
    The slip classification is DERIVED from FEBio's traction output."""
    x = X0 + np.nan_to_num(U)
    n = np.cross(x[tri[:, 1]] - x[tri[:, 0]], x[tri[:, 2]] - x[tri[:, 0]])
    area = 0.5 * np.linalg.norm(n, axis=1)
    n = n / (np.linalg.norm(n, axis=1)[:, None] + 1e-300)
    gap, p, t = D[:, 0], D[:, 1], D[:, 2:5]
    tn = np.einsum("ij,ij->i", t, n)
    tt_vec = t - tn[:, None] * n
    tt = np.linalg.norm(tt_vec, axis=1)
    status = np.zeros(len(p), int)
    act = p > 1e-9
    status[act] = 1
    if mu > 0:
        status[act & (tt >= 0.98 * mu * p)] = 2
    else:
        status[act] = 2          # frictionless: every active point is free to slide
    return {"contact_pressure": p, "contact_gap": gap, "traction_normal": tn, "traction_tangential": tt,
            "friction_utilisation": np.where(act, tt / (mu * p + 1e-300), 0.0) if mu > 0 else np.zeros(len(p)),
            "stick_slip_status": status, "area_deformed": area, "traction_vector": t, "normal_deformed": n,
            "tangential_vector": tt_vec}


def to_vtk(model: FEModel, U=None, cell=None, point=None):
    """pyvista UnstructuredGrid with correct VTK node ordering."""
    import pyvista as pv
    import vtk
    vtk_type = {"tet4": vtk.VTK_TETRA, "penta6": vtk.VTK_WEDGE, "hex8": vtk.VTK_HEXAHEDRON}
    cells, types = [], []
    for b in model.blocks:
        c = b.conn if VTK_ORDER[b.etype] is None else b.conn[:, VTK_ORDER[b.etype]]
        cells.append(np.c_[np.full(len(c), c.shape[1]), c].ravel())
        types.append(np.full(len(c), vtk_type[b.etype], np.uint8))
    g = pv.UnstructuredGrid(np.concatenate(cells), np.concatenate(types), model.nodes.copy())
    bid = np.concatenate([np.full(len(b.conn), i) for i, b in enumerate(model.blocks)])
    g.cell_data["block_id"] = bid
    if U is not None:
        g.point_data["displacement"] = U
        g.point_data["displacement_magnitude"] = np.linalg.norm(U, axis=1)
    for k, v in (cell or {}).items():
        g.cell_data[k] = v
    for k, v in (point or {}).items():
        g.point_data[k] = v
    return g


def initial_quality(model: FEModel) -> dict:
    E = model.n_elements()
    rr = np.full(E, np.nan)
    sj = np.full(E, np.nan)
    off = 0
    for b in model.blocks:
        n = len(b.conn)
        if b.etype == "tet4":
            q = Q.tet_quality(model.nodes, b.conn)
            rr[off:off + n] = q["radius_ratio"]; sj[off:off + n] = q["scaled_jacobian"]
        elif b.etype == "penta6":
            sj[off:off + n] = Q.penta6_corner_jacobians(model.nodes, b.conn).min(axis=1)
        elif b.etype == "hex8":
            sj[off:off + n] = Q.hex8_corner_jacobians(model.nodes, b.conn).min(axis=1)
        off += n
    return {"radius_ratio": rr, "min_scaled_jacobian": sj}


def element_volumes(model: FEModel) -> np.ndarray:
    from ..mesh.layers import _block_volume
    from ..mesh.model import ElementBlock
    out = []
    for b in model.blocks:
        if b.etype == "tet4":
            out.append(np.abs(Q.tet_volume(model.nodes, b.conn)))
        elif b.etype == "penta6":
            c = b.conn
            v = sum(np.abs(Q.tet_volume(model.nodes, c[:, t])) for t in ([0, 1, 2, 3], [1, 2, 3, 4], [2, 3, 4, 5]))
            out.append(v)
        else:
            c = b.conn
            v = sum(np.abs(Q.tet_volume(model.nodes, c[:, t])) for t in
                    ([0, 1, 3, 4], [1, 2, 3, 6], [1, 3, 4, 6], [3, 4, 6, 7], [1, 4, 5, 6]))
            out.append(v)
    return np.concatenate(out)
