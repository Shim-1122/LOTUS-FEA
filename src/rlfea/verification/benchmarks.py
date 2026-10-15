"""Verification benchmarks run through the same FEBio adapter as the limb.

Each benchmark returns a dict with the computed quantity, the reference
(analytical or independently computed) value, the relative error and a
pass/fail verdict against a stated tolerance.
"""
from __future__ import annotations

import json
import os

import numpy as np

from ..mesh.model import ElementBlock, FEModel, Surface
from ..solvers import febio as fb
from . import reference as ref

HEX_CORNERS = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], float)
HEX_TO_TET = np.array([[0, 1, 3, 4], [1, 2, 3, 6], [1, 3, 4, 6], [3, 4, 6, 7], [1, 4, 5, 6]])  # 5-tet split (alternating)
HEX_TO_TET6 = np.array([[0, 1, 2, 6], [0, 2, 3, 6], [0, 3, 7, 6], [0, 7, 4, 6], [0, 4, 5, 6], [0, 5, 1, 6]])


def block_mesh(nx, ny, nz, L=(1.0, 1.0, 1.0), etype="hex8") -> FEModel:
    xs = np.linspace(0, L[0], nx + 1); ys = np.linspace(0, L[1], ny + 1); zs = np.linspace(0, L[2], nz + 1)
    X = np.array([[x, y, z] for z in zs for y in ys for x in xs])
    nid = lambda i, j, k: i + (nx + 1) * (j + (ny + 1) * k)
    H = []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                H.append([nid(i, j, k), nid(i + 1, j, k), nid(i + 1, j + 1, k), nid(i, j + 1, k),
                          nid(i, j, k + 1), nid(i + 1, j, k + 1), nid(i + 1, j + 1, k + 1), nid(i, j + 1, k + 1)])
    H = np.array(H)
    m = FEModel(nodes=X)
    if etype == "hex8":
        m.blocks.append(ElementBlock("block", "hex8", H, "mat", "block", "three-field"))
    elif etype in ("tet4", "tet10"):
        T = H[:, HEX_TO_TET6].reshape(-1, 4)
        if etype == "tet10":
            T, m.nodes = _tet10(T, X)
        m.blocks.append(ElementBlock("block", etype, T, "mat", "block"))
    tol = 1e-9
    Xn = m.nodes
    for ax, nm in enumerate("xyz"):
        m.nodesets[f"{nm}0"] = np.where(np.abs(Xn[:, ax]) < tol)[0]
        m.nodesets[f"{nm}1"] = np.where(np.abs(Xn[:, ax] - L[ax]) < tol)[0]
    # top surface (z = L) quads / tris with outward normal +z
    top = []
    for j in range(ny):
        for i in range(nx):
            top.append([nid(i, j, nz), nid(i + 1, j, nz), nid(i + 1, j + 1, nz), nid(i, j + 1, nz)])
    m.surfaces["top"] = Surface("top", "quad4", np.array(top))
    bot = []
    for j in range(ny):
        for i in range(nx):
            bot.append([nid(i, j, 0), nid(i, j + 1, 0), nid(i + 1, j + 1, 0), nid(i + 1, j, 0)])
    m.surfaces["bottom"] = Surface("bottom", "quad4", np.array(bot))
    return m


def _tet10(T, X):
    """Add mid-edge nodes (FEBio tet10 ordering: 4 corners, then edges
    01,12,20,03,13,23)."""
    edges = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)]
    key = {}
    newX = list(X)
    T10 = np.zeros((len(T), 10), int)
    T10[:, :4] = T
    for e, (a, b) in enumerate(edges):
        for t in range(len(T)):
            k = tuple(sorted((T[t, a], T[t, b])))
            if k not in key:
                key[k] = len(newX)
                newX.append(0.5 * (X[k[0]] + X[k[1]]))
            T10[t, 4 + e] = key[k]
    return T10, np.array(newX)


def _spec_for(model, mat, domain_type=None, domain_params=None, solver=None):
    spec = fb.AnalysisSpec(materials={"mat": mat}, block_material={"block": "mat"},
                           block_domain={"block": {"type": domain_type, "params": domain_params or {}} if domain_type else {}})
    spec.solver = solver or {"max_refs": 25, "dtol": 1e-6, "etol": 1e-8, "rtol": 1e-8, "lstol": 0.9,
                             "qn_method": {"type": "BFGS", "max_ups": 0}, "linear_solver": "pardiso"}
    spec.time_stepper = {"dtmin": 1e-4, "dtmax": 0.25, "max_retries": 10, "opt_iter": 8}
    spec.load_curves = {1: {"points": [(0, 0), (1, 1)]}, 2: {"points": [(0, 0), (1, 1)], "interpolate": "STEP"}}
    spec.plot_vars = ["displacement", "stress"]
    return spec


def run_case(model, spec, workdir, name="case"):
    os.makedirs(workdir, exist_ok=True)
    feb = os.path.join(workdir, f"{name}.feb")
    fb.write_feb(feb, model, spec)
    res = fb.run_febio(feb, workdir, dump=False)
    res["elem_stress"] = fb.read_logdata(os.path.join(workdir, "out_elem_stress.txt"))
    res["node_disp"] = fb.read_logdata(os.path.join(workdir, "out_node_disp.txt"))
    res["node_react"] = fb.read_logdata(os.path.join(workdir, "out_node_react.txt"))
    res["elem_F"] = fb.read_logdata(os.path.join(workdir, "out_elem_F.txt"))
    res["elem_sed"] = fb.read_logdata(os.path.join(workdir, "out_elem_sed.txt"))
    return res


def homogeneous_case(Fend, mat, wref, workdir, etype="hex8", domain_type=None, name="homog", rot=None):
    """Single cube; every node is prescribed u = (F - I) X (ramped).  The
    FEBio Cauchy stress is compared with the independent reference."""
    m = block_mesh(1, 1, 1, etype=etype)
    X = m.nodes
    U = X @ (np.asarray(Fend) - np.eye(3)).T
    spec = _spec_for(m, mat, domain_type)
    step = fb.StepSpec("homogeneous", 1.0, 4, [0.25, 0.5, 0.75, 1.0])
    # nodal displacement varies linearly in time (exact for each node's straight path only if F is
    # reached linearly; intermediate states are not compared)
    for i in range(len(X)):
        m.nodesets[f"n{i}"] = np.array([i])
        for d, dn in enumerate("xyz"):
            step.bcs.append({"type": "prescribed displacement", "node_set": f"n{i}", "dof": dn,
                             "value": float(U[i, d]), "lc": 1, "name": f"n{i}{dn}"})
    spec.steps = [step]
    step.mp_lc = 2
    spec.load_curves[2] = {"points": [(0.0, 0.25), (0.25, 0.25), (0.5, 0.25), (0.75, 0.25), (1.0, 0.25)], "interpolate": "STEP"}
    res = run_case(m, spec, workdir, name)
    out = {"status": res["status"], "workdir": workdir}
    if res["status"] != "converged" or not res["elem_stress"]:
        out["pass"] = False
        return out
    s = res["elem_stress"][-1]["data"]           # (ne, 6): sx sy sz sxy syz sxz
    sig = np.mean(s, axis=0)
    S_fe = np.array([[sig[0], sig[3], sig[5]], [sig[3], sig[1], sig[4]], [sig[5], sig[4], sig[2]]])
    S_ref = wref(np.asarray(Fend))
    scale = max(np.abs(S_ref).max(), 1e-12)
    err = np.abs(S_fe - S_ref).max() / scale
    out.update({"sigma_febio": S_fe.tolist(), "sigma_reference": S_ref.tolist(), "max_rel_error": float(err)})
    return out


# ---------------------------------------------------------------------------
# contact, rigid connectors, locking, energy
# ---------------------------------------------------------------------------
def _two_block(n=4, L=10.0, H=2.0, gap=0.0):
    """Deformable block (n x n x 2 hex) with a rigid plate (n x n x 1 hex) on top."""
    from ..mesh.model import ElementBlock, FEModel, Surface
    a = block_mesh(n, n, 2, L=(L, L, H))
    b = block_mesh(n, n, 1, L=(L, L, 1.0))
    Xb = b.nodes + np.array([0, 0, H + gap])
    m = FEModel(nodes=np.vstack([a.nodes, Xb]))
    off = len(a.nodes)
    m.blocks.append(ElementBlock("soft", "hex8", a.blocks[0].conn, "soft", "soft"))
    m.blocks.append(ElementBlock("plate", "hex8", b.blocks[0].conn + off, "plate", "plate", "rigid"))
    m.surfaces["soft_top"] = Surface("soft_top", "quad4", a.surfaces["top"].conn)
    m.surfaces["plate_bottom"] = Surface("plate_bottom", "quad4", b.surfaces["bottom"].conn + off)
    m.nodesets["base"] = a.nodesets["z0"]
    return m


def _rigid_rows(workdir):
    recs = fb.read_logdata(os.path.join(workdir, "out_rigid.txt"))
    return [(r["time"], r["data"]) for r in recs]


def contact_slide_case(workdir, mu=0.3, quick=False, slide=2.0, press=0.1, n=4, n_steps=10, frictionless_energy=False):
    """Rigid plate pressed into a compliant block, then slid tangentially.
    Coulomb: once gross sliding occurs |Fx| = mu |Fz|."""
    m = _two_block(n=n)
    mat = {"febio_type": "neo-Hookean", "params": {"E": 1.0, "v": 0.3}}
    spec = _spec_for(m, mat)
    spec.materials = {"soft": mat, "plate": {"rigid": True, "com": m.nodes[m.blocks[1].conn].reshape(-1, 3).mean(0)}}
    spec.block_material = {"soft": "soft", "plate": "plate"}
    spec.block_domain = {}
    spec.bcs = [{"type": "zero displacement", "node_set": "base", "dofs": "xyz"}]
    spec.rigid_bcs = [{"type": "rigid_fixed", "rb": "plate", "dofs": ["Ry", "Ru", "Rv", "Rw"]}]
    spec.contacts = [{"type": "sliding-elastic", "name": "c", "primary": "soft_top", "secondary": "plate_bottom",
                      "params": {"laugon": "PENALTY", "penalty": 1.0, "auto_penalty": 1, "two_pass": 0, "fric_coeff": mu,
                                 "symmetric_stiffness": 0, "tension": 0, "search_tol": 0.01}}]
    spec.log_surfaces = ["soft_top"]; spec.log_rigid = True
    spec.plot_vars = spec.plot_vars + ["contact traction", "contact stick", "contact pressure", "contact gap"]
    spec.solver["symmetric_stiffness"] = 0
    spec.load_curves = {1: {"points": [(0, 0), (1, 1)]}, 2: {"points": [(1, 0), (2, 1)]}}
    s1 = fb.StepSpec("press", 1.0, n_steps, [0.5, 1.0], mp_lc=3)
    s1.rigid_bcs = [{"type": "rigid_displacement", "rb": "plate", "dof": "z", "value": -press, "lc": 1, "relative": 0},
                    {"type": "rigid_fixed", "rb": "plate", "dofs": ["Rx"]}]
    s2 = fb.StepSpec("slide", 1.0, n_steps, [1.0 + k / 10 for k in range(1, 11)], mp_lc=4)
    s2.rigid_bcs = [{"type": "rigid_fixed", "rb": "plate", "dofs": ["Rz"]},
                    {"type": "rigid_displacement", "rb": "plate", "dof": "x", "value": slide, "lc": 2, "relative": 1}]
    spec.steps = [s1, s2]
    spec.load_curves[3] = {"points": [(0, 0.1), (0.5, 0.1), (1.0, 0.1)], "interpolate": "STEP"}
    spec.load_curves[4] = {"points": [(1.0, 0.1)] + [(t, 0.1) for t in s2.must_points], "interpolate": "STEP"}
    res = run_case(m, spec, workdir, "slide")
    out = {"status": res["status"]}
    rows = _rigid_rows(workdir)
    if rows:
        F = np.array([r[1][0][6:9] for r in rows]); t = np.array([r[0] for r in rows])
        out["time"] = t.tolist(); out["Fx"] = F[:, 0].tolist(); out["Fz"] = F[:, 2].tolist()
        out["ratio_end"] = float(abs(F[-1, 0]) / max(abs(F[-1, 2]), 1e-300))
        out["ratio_first_slide_output"] = float(abs(F[t > 1][0, 0]) / max(abs(F[t > 1][0, 2]), 1e-300)) if (t > 1).any() else None
    return out


def contact_opening_case(workdir):
    """Press then lift beyond the initial position: no tensile contact
    force may remain (compression-only, separation allowed)."""
    m = _two_block(n=4)
    mat = {"febio_type": "neo-Hookean", "params": {"E": 1.0, "v": 0.3}}
    spec = _spec_for(m, mat)
    spec.materials = {"soft": mat, "plate": {"rigid": True, "com": m.nodes[m.blocks[1].conn].reshape(-1, 3).mean(0)}}
    spec.block_material = {"soft": "soft", "plate": "plate"}
    spec.block_domain = {}
    spec.bcs = [{"type": "zero displacement", "node_set": "base", "dofs": "xyz"}]
    spec.rigid_bcs = [{"type": "rigid_fixed", "rb": "plate", "dofs": ["Rx", "Ry", "Ru", "Rv", "Rw"]}]
    spec.contacts = [{"type": "sliding-elastic", "name": "c", "primary": "soft_top", "secondary": "plate_bottom",
                      "params": {"laugon": "PENALTY", "penalty": 1.0, "auto_penalty": 1, "fric_coeff": 0.0, "tension": 0}}]
    spec.log_surfaces = ["soft_top"]; spec.log_rigid = True
    spec.load_curves = {1: {"points": [(0, 0), (1, 1), (2, -1.0)]}, 3: {"points": [(0, 0.1), (1.0, 0.1), (1.5, 0.1), (2.0, 0.1)], "interpolate": "STEP"}}
    s1 = fb.StepSpec("press_lift", 2.0, 20, [1.0, 1.5, 2.0], mp_lc=3)
    s1.rigid_bcs = [{"type": "rigid_displacement", "rb": "plate", "dof": "z", "value": -0.1, "lc": 1, "relative": 0}]
    spec.steps = [s1]
    res = run_case(m, spec, workdir, "opening")
    out = {"status": res["status"]}
    rows = _rigid_rows(workdir)
    faces = fb.read_logdata(os.path.join(workdir, "out_face_soft_top.txt"))
    if rows:
        Fz = np.array([r[1][0][8] for r in rows])
        out["Fz"] = Fz.tolist(); out["times"] = [r[0] for r in rows]
        out["Fz_end_over_Fz_max"] = float(abs(Fz[-1]) / max(np.abs(Fz).max(), 1e-300))
    if faces:
        out["max_pressure_end"] = float(faces[-1]["data"][:, 1].max())
        out["gap_end_min_max"] = [float(faces[-1]["data"][:, 0].min()), float(faces[-1]["data"][:, 0].max())]
    return out


def rigid_hinge_case(workdir, M=100.0, k=2.0, Lr=20.0):
    """Rigid body B on a revolute joint (axis z through origin) to fixed body
    A, loaded by a couple M about z and restrained by a rigid spring of
    stiffness k attached at radius Lr.  Small rotation: theta = M/(k Lr^2)."""
    from ..mesh.model import ElementBlock, FEModel
    X = np.vstack([HEX_CORNERS * [4, 4, 4] + [-2, -2, -10], HEX_CORNERS * [4, 4, 4] + [Lr - 2, -2, -2]])
    m = FEModel(nodes=X)
    m.blocks.append(ElementBlock("A", "hex8", np.array([[0, 1, 2, 3, 4, 5, 6, 7]]), "A", "A", "rigid"))
    m.blocks.append(ElementBlock("B", "hex8", np.array([[8, 9, 10, 11, 12, 13, 14, 15]]), "B", "B", "rigid"))
    spec = _spec_for(m, {})
    spec.materials = {"A": {"rigid": True, "com": [0, 0, -8.0]}, "B": {"rigid": True, "com": [Lr, 0, 0.0]}}
    spec.block_material = {"A": "A", "B": "B"}; spec.block_domain = {}
    spec.rigid_bcs = [{"type": "rigid_fixed", "rb": "A", "dofs": ["Rx", "Ry", "Rz", "Ru", "Rv", "Rw"]}]
    spec.connectors = [{"type": "rigid revolute joint", "name": "hinge", "params": {
        "body_a": "A", "body_b": "B", "laugon": "PENALTY", "tolerance": 0, "force_penalty": 1e5, "moment_penalty": 1e8,
        "auto_penalty": 0, "joint_origin": [0, 0, 0.0], "rotation_axis": [0, 0, 1.0], "transverse_axis": [1, 0, 0.0]}},
        {"type": "rigid spring", "name": "spring", "params": {"body_a": "B", "body_b": "A", "k": k,
                                                              "insertion_a": [Lr, 0, 0.0], "insertion_b": [Lr, -1000.0, 0.0], "free_length": 1000.0}}]
    spec.rigid_loads = [{"type": "rigid_moment", "rb": "B", "dof": "Rw", "value": M, "lc": 1}]
    spec.log_rigid = True; spec.log_connectors = True
    spec.load_curves = {1: {"points": [(0, 0), (1, 1)]}, 2: {"points": [(0, 0.25), (1, 0.25)], "interpolate": "STEP"}}
    s = fb.StepSpec("load", 1.0, 4, [1.0], mp_lc=2)
    spec.steps = [s]
    res = run_case(m, spec, workdir, "hinge")
    out = {"status": res["status"]}
    rows = _rigid_rows(workdir)
    if rows:
        d = rows[-1][1]
        b = d[1]                                   # body B (second rigid material)
        theta = float(b[5])                        # thz
        out["theta_rad"] = theta
        # exact for the finite rotation with a long (nearly parallel) spring: M = k Lr^2 sin(theta) cos(theta)
        out["theta_small_angle_ref"] = M / (k * Lr ** 2)
        from scipy.optimize import brentq
        out["theta_ref"] = brentq(lambda th: k * Lr ** 2 * np.sin(th) * np.cos(th) * (1 - 0) - M, 0, 1.2)
        out["rel_error"] = abs(theta - out["theta_ref"]) / out["theta_ref"]
    return out


def cantilever_locking_case(workdir, etype="hex8", domain=None, n=(10, 2, 2), delta=1.0, k_over_c=1e4):
    """Near-incompressible cantilever (40 x 4 x 4 mm, Mooney-Rivlin c1, k = k_over_c*c1);
    tip face prescribed vertical displacement; returns total reaction (stiffness).
    Locking shows up as an over-stiff reaction compared with the converged
    mixed-formulation reference."""
    m = block_mesh(n[0], n[1], n[2], L=(40.0, 4.0, 4.0), etype=etype)
    c1 = 0.1
    mat = {"febio_type": "Mooney-Rivlin", "params": {"c1": c1, "c2": 0.0, "k": k_over_c * c1}}
    spec = _spec_for(m, mat, domain)
    spec.bcs = [{"type": "zero displacement", "node_set": "x0", "dofs": "xyz"}]
    st = fb.StepSpec("bend", 1.0, 5, [1.0], mp_lc=2)
    st.bcs = [{"type": "prescribed displacement", "node_set": "x1", "dof": "z", "value": delta, "lc": 1}]
    spec.steps = [st]
    spec.load_curves[2] = {"points": [(0, 0.2), (1.0, 0.2)], "interpolate": "STEP"}
    res = run_case(m, spec, workdir, f"cant_{etype}_{domain}")
    out = {"status": res["status"], "etype": etype, "domain": domain, "n": list(n), "n_nodes": len(m.nodes)}
    if res["node_react"]:
        R = res["node_react"][-1]
        idx = np.isin(R["ids"] - 1, m.nodesets["x1"])
        out["tip_force_N"] = float(R["data"][idx, 2].sum())
        E = 6 * c1
        I = 4 * 4 ** 3 / 12
        out["euler_bernoulli_force_N"] = 3 * E * I * delta / 40 ** 3
    if res["elem_F"]:
        F = res["elem_F"][-1]["data"].reshape(-1, 3, 3)
        J = np.linalg.det(F)
        out["J_range"] = [float(J.min()), float(J.max())]
    return out


def energy_balance_case(workdir, n=4, n_mp=10):
    """Frictionless indentation: external work (trapezoid of reaction x
    displacement over must points) vs. stored strain energy (sum sed*V0)."""
    m = _two_block(n=n)
    mat = {"febio_type": "neo-Hookean", "params": {"E": 1.0, "v": 0.3}}
    spec = _spec_for(m, mat)
    spec.materials = {"soft": mat, "plate": {"rigid": True, "com": m.nodes[m.blocks[1].conn].reshape(-1, 3).mean(0)}}
    spec.block_material = {"soft": "soft", "plate": "plate"}; spec.block_domain = {}
    spec.bcs = [{"type": "zero displacement", "node_set": "base", "dofs": "xyz"}]
    spec.rigid_bcs = [{"type": "rigid_fixed", "rb": "plate", "dofs": ["Rx", "Ry", "Ru", "Rv", "Rw"]}]
    spec.contacts = [{"type": "sliding-elastic", "name": "c", "primary": "soft_top", "secondary": "plate_bottom",
                      "params": {"laugon": "AUGLAG", "tolerance": 0.01, "penalty": 1.0, "auto_penalty": 1, "fric_coeff": 0.0}}]
    spec.log_surfaces = ["soft_top"]; spec.log_rigid = True
    mps = [k / n_mp for k in range(1, n_mp + 1)]
    spec.load_curves = {1: {"points": [(0, 0), (1, 1)]}, 2: {"points": [(0, 0.1)] + [(t, 0.1) for t in mps], "interpolate": "STEP"}}
    s = fb.StepSpec("press", 1.0, n_mp, mps, mp_lc=2)
    s.rigid_bcs = [{"type": "rigid_displacement", "rb": "plate", "dof": "z", "value": -0.2, "lc": 1}]
    spec.steps = [s]
    res = run_case(m, spec, workdir, "energy")
    out = {"status": res["status"]}
    rows = _rigid_rows(workdir)
    if rows and res["elem_sed"]:
        t = np.r_[0, [r[0] for r in rows if r[0] > 0]]
        Fz = np.r_[0, [r[1][0][8] for r in rows if r[0] > 0]]
        z = -0.2 * t
        W = float(np.trapezoid(-Fz, -z)) if hasattr(np, "trapezoid") else float(np.trapz(-Fz, -z))
        from ..post.results import element_volumes
        V0 = element_volumes(m)
        sed = res["elem_sed"][-1]["data"][:, 0]
        U = float((sed * V0[: len(sed)]).sum())
        out.update({"external_work_Nmm": abs(W), "strain_energy_Nmm": U, "rel_diff": abs(abs(W) - U) / max(U, 1e-300)})
    return out


def indentation_locking_case(workdir, etype="hex8", domain=None, dom_params=None, n=8, k_over_mu=100.0, delta=1.0):
    """Compression-dominated locking check representative of soft tissue
    under a pad: 40 x 40 x 20 mm block, base fixed, central 10 x 10 mm top
    patch pushed down by delta.  Returns the patch reaction (stiffness)."""
    m = block_mesh(2 * n, 2 * n, n, L=(40.0, 40.0, 20.0), etype=etype)
    mu = 0.005
    c1 = mu / 2
    mat = {"febio_type": "Mooney-Rivlin", "params": {"c1": c1, "c2": 0.0, "k": k_over_mu * mu}}
    spec = _spec_for(m, mat, domain)
    if domain:
        spec.block_domain = {"block": {"type": domain, "params": dom_params or {}}}
    X = m.nodes
    patch = np.where((np.abs(X[:, 2] - 20) < 1e-9) & (np.abs(X[:, 0] - 20) <= 5 + 1e-9) & (np.abs(X[:, 1] - 20) <= 5 + 1e-9))[0]
    m.nodesets["patch"] = patch
    spec.bcs = [{"type": "zero displacement", "node_set": "z0", "dofs": "xyz"}]
    st = fb.StepSpec("indent", 1.0, 5, [1.0], mp_lc=2)
    st.bcs = [{"type": "prescribed displacement", "node_set": "patch", "dof": "z", "value": -delta, "lc": 1}]
    spec.steps = [st]
    spec.load_curves[2] = {"points": [(0, 0.2), (1.0, 0.2)], "interpolate": "STEP"}
    res = run_case(m, spec, workdir, f"indent_{etype}")
    out = {"status": res["status"], "etype": etype, "domain": domain, "params": dom_params, "n": n, "n_nodes": len(X)}
    if res["node_react"]:
        R = res["node_react"][-1]
        idx = np.isin(R["ids"] - 1, patch)
        out["patch_force_N"] = float(-R["data"][idx, 2].sum())
    return out


def restart_case(workdir, mu=0.3, kill_after_dumps=5):
    """Restart verification: the Coulomb slide model (rigid plate, frictional
    contact, two steps) is solved once uninterrupted, and once killed with
    SIGKILL after `kill_after_dumps` restart points and resumed with
    `febio4 -r`.  The resumed run must terminate normally and reproduce the
    uninterrupted rigid-body reaction history at every common output time."""
    import re
    import shutil
    import signal
    import subprocess
    import time
    ref_dir, rst_dir = os.path.join(workdir, "uninterrupted"), os.path.join(workdir, "killed_resumed")
    ref_out = contact_slide_case(ref_dir, mu=mu)
    ref_rows = dict((round(t, 6), np.array(d[0][6:9])) for t, d in _rigid_rows(ref_dir))

    orig = fb.run_febio
    state = {}

    def run_and_kill(feb, wd, **kw):
        exe = fb.find_febio()
        log = os.path.join(wd, "slide.log")
        p = subprocess.Popen([exe, "-nosplash", "-silent", "-noconfig", "-i", os.path.basename(feb), "-dump", "slide.dmp"],
                             cwd=wd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL)
        while p.poll() is None:
            if os.path.exists(log) and len(re.findall("Restart point created", open(log, errors="replace").read())) >= kill_after_dumps:
                p.send_signal(signal.SIGKILL)
                break
            time.sleep(0.002)
        p.wait()
        state["rc"] = p.returncode
        conv = fb.parse_log(log)
        state["killed_at"] = conv.get("last_converged_time")
        return {"status": fb.classify(p.returncode, conv), "returncode": p.returncode}

    shutil.rmtree(rst_dir, ignore_errors=True)
    fb.run_febio = run_and_kill
    try:
        first = contact_slide_case(rst_dir, mu=mu)
    finally:
        fb.run_febio = orig
    rst = fb.write_restart(os.path.join(rst_dir, "restart.feb"), os.path.join(rst_dir, "slide.dmp"))
    res = fb.run_febio(os.path.join(rst_dir, "slide.feb"), rst_dir, restart=rst)
    rows = dict((round(t, 6), np.array(d[0][6:9])) for t, d in _rigid_rows(rst_dir))
    common = sorted(set(rows) & set(ref_rows))
    scale = max(np.abs(np.array(list(ref_rows.values()))).max(), 1e-300)
    diff = max((float(np.abs(rows[t] - ref_rows[t]).max()) / scale for t in common), default=float("inf"))
    out = {"first_status": first["status"], "killed_returncode": state.get("rc"), "killed_at_time": state.get("killed_at"),
           "restart_status": res["status"], "restart_returncode": res["returncode"],
           "final_time": res["convergence"].get("last_converged_time"), "reference_final_time": ref_out.get("time", [None])[-1],
           "common_output_times": len(common), "max_rel_reaction_difference": diff}
    out["pass"] = bool(first["status"] == "interrupted" and res["status"] == "converged"
                       and abs((out["final_time"] or 0) - 2.0) < 1e-9 and len(common) == len(ref_rows) and diff < 1e-6)
    return out
