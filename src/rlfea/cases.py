"""Assemble a complete analysis from a resolved configuration.

geometry -> audit/repair log -> layered conforming mesh -> (socket) ->
materials -> supports -> contact -> cable loads -> load steps -> AnalysisSpec
"""
from __future__ import annotations

import copy

import numpy as np

from .geometry import repair as R
from .geometry.socket import build_socket
from .geometry.surface import TriSurface, audit_assembly, audit_body, load_body
from .geometry.synthetic import synthetic_limb
from .loads import cable as CB
from .materials.library import load_material_set
from .mesh import quality as Q
from .mesh.layers import build_layered_model, detect_cut_faces
from .solvers.febio import AnalysisSpec, StepSpec


# ---------------------------------------------------------------------------
# geometry
# ---------------------------------------------------------------------------
def prepare_geometry(cfg: dict) -> dict:
    g = cfg["geometry"]
    out = {"audit_source": {}, "audit_prepared": {}, "missing_anatomy": [], "roles": {}, "notes": []}
    if g["source"] == "synthetic":
        bodies = synthetic_limb(length=g["synthetic"]["length"])
        skin = R.basic_cleanup(bodies.pop("skin"))
        bones = {k: R.basic_cleanup(v) for k, v in bodies.items()}
        out["roles"] = {"tibia": "tibia (synthetic)", "fibula": "fibula (synthetic)"}
        out["notes"].append("SYNTHETIC geometry: not patient data, not anatomically measured.")
        for k, s in [("skin", skin)] + list(bones.items()):
            out["audit_source"][k] = audit_body(s)
    else:
        skin_spec = g["skin"]
        sk = load_body(skin_spec["path"], skin_spec.get("unit"), name="skin_envelope",
                       brep_mesh_size=g.get("brep_tessellation_size", 2.0))
        if len(sk) != 1:
            raise ValueError("skin file must contain exactly one body")
        skin = sk[0]
        out["audit_source"]["skin"] = audit_body(skin)
        bones = {}
        for b in g.get("bones", []):
            bl = load_body(b["path"], b.get("unit"), name=b["name"])
            for body in bl:
                a = audit_body(body)
                out["audit_source"][b["name"]] = a
                if body.empty:
                    out["notes"].append(f"{b['name']}: EMPTY body in {b['path']} - excluded (no geometry).")
                    continue
                if not b.get("include", True):
                    out["notes"].append(f"{b['name']}: excluded by configuration ({b.get('reason', 'no reason given')}).")
                    continue
                bones[b["name"]] = R.basic_cleanup(body)
                out["roles"][b["name"]] = b.get("role", "unidentified")
        for e in g.get("ignored_files", []):
            ign = load_body(e["path"], e.get("unit"), name=e["name"])
            for body in ign:
                out["audit_source"][e["name"]] = audit_body(body)
                out["notes"].append(f"{e['name']}: listed as ignored ({e.get('reason', '')}).")
        skin = R.basic_cleanup(skin)
    # expected anatomy check
    for exp in g.get("expected_bones", []):
        if not any(exp.lower() in str(r).lower() for r in out["roles"].values()):
            out["missing_anatomy"].append(exp)
    # alignment into the limb frame (+z proximal, origin at distal apex)
    al = g.get("alignment", {"method": "none"})
    if al.get("method") in ("pca", "axes"):
        M = _pca_frame(skin, np.asarray(al.get("proximal_hint", [0, 0, 1.0]), float),
                       np.asarray(al.get("anterior_hint", [0, 1.0, 0]), float), use_pca=al["method"] == "pca")
        skin = R.transform(skin, M, f"limb frame ({al['method']}): long axis -> +z, anterior -> +y, distal apex -> z=0")
        bones = {k: R.transform(v, M, "same limb-frame transform as skin") for k, v in bones.items()}
        out["frame_transform"] = M.tolist()
    elif al.get("method") == "translate_apex":
        apex_z = skin.vertices[:, 2].min()
        M = np.eye(4); M[2, 3] = -apex_z
        skin = R.transform(skin, M, "translate distal apex to z = 0")
        bones = {k: R.transform(v, M, "same translation") for k, v in bones.items()}
        out["frame_transform"] = M.tolist()
    # FE remeshing (logged)
    rm = g.get("remesh", {})
    if rm.get("skin_size"):
        skin = R.remesh(skin, rm["skin_size"])
    if rm.get("bone_size"):
        bones = {k: R.remesh(v, rm["bone_size"]) for k, v in bones.items()}
    out["audit_prepared"]["skin"] = audit_body(skin)
    for k, v in bones.items():
        out["audit_prepared"][k] = audit_body(v)
    out["assembly"] = audit_assembly(skin, list(bones.values()))
    out["skin"] = skin
    out["bones"] = bones
    out["history"] = {"skin": skin.history, **{k: v.history for k, v in bones.items()}}
    return out


def _pca_frame(skin: TriSurface, proximal_hint, anterior_hint, use_pca: bool = True) -> np.ndarray:
    """Rigid transform into the limb frame.  use_pca=False keeps the source
    long axis given by proximal_hint (preferred when the source model has a
    flat proximal cut perpendicular to a known axis)."""
    V = skin.vertices
    c = V.mean(0)
    if use_pca:
        _, _, vt = np.linalg.svd(V - c, full_matrices=False)
        ez = vt[0] if vt[0] @ proximal_hint > 0 else -vt[0]
    else:
        ez = proximal_hint / np.linalg.norm(proximal_hint)
    ey = anterior_hint - (anterior_hint @ ez) * ez
    ey /= np.linalg.norm(ey)
    ex = np.cross(ey, ez)
    Rm = np.vstack([ex, ey, ez])           # rows: new axes in old coordinates
    Vn = (V - c) @ Rm.T
    apex = Vn[np.argmin(Vn[:, 2])]
    # put the limb axis through the centroid of the distal half (x,y), apex at z=0
    lower = Vn[Vn[:, 2] < np.percentile(Vn[:, 2], 50)]
    shift = np.array([-lower[:, 0].mean(), -lower[:, 1].mean(), -apex[2]])
    M = np.eye(4)
    M[:3, :3] = Rm
    M[:3, 3] = -Rm @ c + shift
    return M


# ---------------------------------------------------------------------------
# model + analysis
# ---------------------------------------------------------------------------
def build_case(cfg: dict, geo: dict) -> dict:
    mset = load_material_set(cfg["materials"]["set"])
    L = cfg["layers"]
    lcfg = dict(skin_thickness=L["skin_thickness"], fat_thickness=L["fat_thickness"],
                max_fraction_skin=L["max_fraction_skin"], max_fraction_fat=L["max_fraction_fat"],
                min_layer=L["min_layer"], normal_smoothing=L["normal_smoothing"],
                thickness_smoothing=L["thickness_smoothing"], axis=[0, 0, 1.0],
                skin_fiber_reference=L["skin_fiber_reference"], muscle_formulation=L["muscle_formulation"],
                muscle_fiber=L["muscle_fiber"], optimize=L["optimize"], mesh_bones=True,
                proximal_cut=cfg["geometry"].get("proximal_cut", True))
    model, layer_rep = build_layered_model(geo["skin"], geo["bones"], lcfg)
    sock_rep = None
    S = cfg.get("socket", {"enabled": False})
    mode = cfg.get("load_model", "resolved_contact")
    traction_model = None
    if S.get("enabled"):
        if mode == "resolved_contact":
            model, sock_rep = build_socket(model, geo["skin"], S)
        elif mode == "prescribed_traction":
            probe = copy.deepcopy(model)
            probe, sock_rep = build_socket(probe, geo["skin"], S)
            traction_model = probe
            for k, v in probe.surfaces.items():
                if k.startswith("skin_"):
                    model.surfaces[k] = v
        else:
            raise ValueError(f"unknown load_model {mode}")
    spec, load_rep = _assemble(cfg, model, mset, traction_model)
    return {"model": model, "spec": spec, "material_set": mset, "layer_report": layer_rep,
            "socket_report": sock_rep, "load_report": load_rep, "load_model": mode}


def _rigid_com(X, conn, etype):
    """Volume-weighted centroid of a rigid block (uniform density)."""
    from .mesh.model import ElementBlock, FEModel
    from .post.results import element_volumes
    v = element_volumes(FEModel(nodes=X, blocks=[ElementBlock("t", etype, conn, "", "")]))
    c = X[conn].mean(axis=1)
    if not np.isfinite(v).all() or v.sum() <= 0:
        raise ValueError("rigid block with zero/invalid volume")
    return (c * v[:, None]).sum(0) / v.sum()


def _assemble(cfg, model, mset, traction_model=None):
    mats = mset["materials"]
    reg = cfg["materials"]["regions"]
    spec = AnalysisSpec(materials={}, block_material={}, block_domain={})
    rep: dict = {"assumptions": []}
    X = model.nodes
    bone_mode = cfg["materials"].get("bones", "rigid")
    supports = cfg["supports"]
    for b in model.blocks:
        if b.formulation == "rigid" or (b.region.startswith("bone:") and bone_mode == "rigid"):
            mname = f"rigid_{b.name}" if not b.name.startswith("petal") else f"rigid_{b.name.split('_')[0]}"
            if b.name == "distal_cup":
                mname = "rigid_cup"
            com = _rigid_com(X, b.conn, b.etype)
            if mname in spec.materials:
                # several blocks may form one rigid body (not used currently)
                pass
            spec.materials[mname] = {"rigid": True, "com": com}
            spec.block_material[b.name] = mname
            continue
        if b.region.startswith("bone:"):
            key = cfg["materials"].get("bone_material", "bone_deformable")
        elif b.name.endswith("_pad"):
            key = reg.get("silicone_pad", "silicone_pad")
        else:
            key = reg[b.name] if b.name in reg else reg.get(b.material, b.material)
        if key not in mats:
            raise KeyError(f"block {b.name}: material {key!r} not in set {mset['set_name']}")
        m = mats[key]
        spec.materials.setdefault(key, {"febio_type": m["febio_type"], "params": dict(m["params"])})
        spec.block_material[b.name] = key
        uncoupled = m["febio_type"] in ("Holzapfel-Gasser-Ogden", "Ogden", "Mooney-Rivlin", "trans iso Mooney-Rivlin")
        if b.etype == "tet4" and b.formulation == "ut4":
            spec.block_domain[b.name] = {"type": "ut4-solid", "params": {"alpha": cfg["layers"].get("ut4_alpha", 0.05), "iso_stab": 0}}
        elif b.etype in ("hex8", "penta6") and uncoupled:
            spec.block_domain[b.name] = {"type": "three-field-solid", "params": {}}
        else:
            spec.block_domain[b.name] = {}
    # --- tissue supports -------------------------------------------------
    pc = supports.get("proximal_cut", "axial_roller")
    if "proximal_cut" in model.nodesets and pc != "free":
        dofs = {"fixed": "xyz", "axial_roller": "z"}[pc]
        _support(spec, "proximal_cut", dofs, f"proximal_cut_{pc}")
    rep["supports"] = {"proximal_cut": pc}
    for bname, how in supports.get("bones", {}).items():
        mname = f"rigid_bone_{bname}"
        if bone_mode == "rigid":
            if mname not in spec.materials:
                raise KeyError(f"support for unknown bone {bname}")
            if how == "fixed":
                spec.rigid_bcs.append({"type": "rigid_fixed", "rb": mname, "dofs": ["Rx", "Ry", "Rz", "Ru", "Rv", "Rw"],
                                       "name": f"fix_{bname}"})
            elif how != "free":
                raise ValueError(f"bone support {how!r} not supported (fixed|free)")
        else:
            if how == "fixed":
                blk = model.block(f"bone_{bname}")
                nodes = np.unique(blk.conn.ravel())
                ztop = X[nodes, 2].max()
                sel = nodes[X[nodes, 2] > ztop - cfg["supports"].get("deformable_bone_fixed_length_mm", 10.0)]
                model.nodesets[f"bone_{bname}_proximal"] = sel
                _support(spec, f"bone_{bname}_proximal", "xyz", f"fix_{bname}")
        rep["supports"][bname] = how
    # --- socket: rigid cup, hinges, stabilizers, contact ---------------------
    S = cfg.get("socket", {})
    C = cfg["contact"]
    petals = model.meta.get("petals", []) if traction_model is None else traction_model.meta.get("petals", [])
    lc = {"cable": 1, "case2": 2}
    spec.load_curves[1] = {"name": "cable_tightening", "points": [(0.0, 0.0), (1.0, 1.0)]}
    nonsym = False
    if S.get("enabled") and traction_model is None:
        spec.rigid_bcs.append({"type": "rigid_fixed", "rb": "rigid_cup", "dofs": ["Rx", "Ry", "Ru", "Rv", "Rw"], "name": "cup_lateral_rot"})
        k_stab = float(S.get("stabilizer_k", 0.05))
        for p in petals:
            rb = f"rigid_{p.name}"
            spec.materials[rb]["com"] = _rigid_com(X, model.block(f"{p.name}_shell").conn, "hex8")
            radial = np.cross(p.hinge_axis, [0, 0, 1.0])
            spec.connectors.append({"type": "rigid revolute joint", "name": f"hinge_{p.name}", "params": {
                "body_a": "rigid_cup", "body_b": rb, "laugon": "PENALTY", "tolerance": 0.0,
                "force_penalty": 1e4, "moment_penalty": 1e7, "auto_penalty": 0,
                "joint_origin": p.hinge_origin, "rotation_axis": p.hinge_axis, "transverse_axis": radial}})
            if k_stab > 0:
                top = X[p.shell_outer_nodes[:, -1]].mean(0) * 0.5 + X[p.pad_outer_nodes[:, -1]].mean(0) * 0.5
                anchor = np.array([0.0, 0.0, top[2]])
                spec.connectors.append({"type": "rigid spring", "name": f"stabilizer_{p.name}", "params": {
                    "body_a": rb, "body_b": "rigid_cup", "k": k_stab, "insertion_a": top, "insertion_b": anchor,
                    "free_length": float(np.linalg.norm(top - anchor))}})
        rep["assumptions"].append(f"Each petal hinge is a frictionless revolute joint (penalty 1e4 N/mm, 1e7 N mm/rad). "
                                  f"A radial stabiliser spring k = {k_stab} N/mm per petal is NUMERICAL STABILISATION; "
                                  "its force and energy are reported.")
        ctype = "tied-elastic" if C.get("bonded") else C["type"]
        pairs = [(f"skin_{p.name}_zone", f"{p.name}_pad_contact", f"contact_{p.name}") for p in petals]
        pairs.append(("skin_cushion_zone", "cushion_contact", "contact_cushion"))
        for prim, sec, name in pairs:
            if len(model.surfaces[prim].conn) == 0:
                continue
            if ctype == "sliding-elastic":
                params = {"laugon": "AUGLAG" if C["augmented_lagrangian"] else "PENALTY", "tolerance": 0.2, "gaptol": 0,
                          "penalty": C["penalty_scale"], "auto_penalty": int(C["auto_penalty"]), "two_pass": int(C["two_pass"]),
                          "search_tol": C["search_tol"], "search_radius": C["search_radius"], "symmetric_stiffness": 0,
                          "tension": 0, "fric_coeff": C["friction"], "minaug": 0, "maxaug": 10}
                nonsym = nonsym or C["friction"] > 0
            else:
                params = {"laugon": "PENALTY", "tolerance": 0.2, "penalty": C["penalty_scale"], "auto_penalty": 1}
            spec.contacts.append({"type": ctype, "name": name, "primary": prim, "secondary": sec, "params": params})
        spec.log_surfaces = [p for p, _, _ in pairs if len(model.surfaces[p].conn)]
        spec.log_rigid = True
        spec.log_connectors = True
        rep["contact"] = {"type": ctype, "friction": C["friction"], "enforcement": "augmented Lagrangian" if C["augmented_lagrangian"] else "penalty",
                          "penalty_scale": C["penalty_scale"], "auto_penalty": bool(C["auto_penalty"]),
                          "separation_allowed": ctype == "sliding-elastic", "compression_only": ctype == "sliding-elastic"}
    elif any(spec.materials[m].get("rigid") for m in spec.materials):
        spec.log_rigid = True
    # --- cable loads ---------------------------------------------------------
    cb = cfg["cable"]
    T_example = cb.get("strand_tension")
    T_rated = cb.get("rated_tension_N")
    T_apply = T_rated if (cfg.get("tension_basis") == "rated" and T_rated) else T_example
    rep["tension"] = {"applied_N": T_apply, "basis": "rated" if T_apply == T_rated and T_rated else "example",
                      "example_N": T_example, "rated_N": T_rated,
                      "levels_N": [f * T_apply for f in cb["rating_fractions"]] if T_apply else []}
    if S.get("enabled") and T_apply:
        tm = traction_model if traction_model is not None else model
        Xs = tm.nodes

        def shell_r(p, th, z):
            nodes = p.shell_outer_nodes
            P = Xs[nodes.ravel()]
            thp = np.arctan2(P[:, 1], P[:, 0])
            d = np.abs(((thp - th) + np.pi) % (2 * np.pi) - np.pi) * 60 + np.abs(P[:, 2] - z)
            k = np.argsort(d)[:4]
            w = 1.0 / (d[k] + 1e-9)
            return float((np.linalg.norm(P[k, :2], axis=1) * w).sum() / w.sum())

        if cb.get("routing", "default") == "default":
            route = CB.default_route(petals, cb["route_height_fraction"], cb["guide_inset"], cb["dial_petal"], shell_r)
            rep["routing_status"] = "PROVISIONAL default single-row lace generated from the parametric petals"
        else:
            route = CB.route_from_config(cb["routing"])
            rep["routing_status"] = "explicit routing from configuration"
        coms = {f"{p.name}": (tm and _rigid_com(Xs, tm.block(f"{p.name}_shell").conn, "hex8")) for p in petals}
        bl = CB.body_loads(route, T_apply, coms, cb["guide_friction"])
        rep["cable"] = {"route": [{"body": r.body, "xyz": r.xyz.tolist(), "label": r.label} for r in route],
                        "body_loads": {k: {"F_N": v["F"].tolist(), "M_com_Nmm": v["M_com"].tolist(),
                                           "point_forces": v["point_forces"]} for k, v in bl["bodies"].items()},
                        "spans": bl["spans"], "sum_force_N": bl["sum_force_N"],
                        "sum_moment_about_origin_Nmm": bl["sum_moment_about_origin_Nmm"],
                        "segment_tensions_N": bl["segment_tensions_N"], "guide_friction": cb["guide_friction"],
                        "coms": {k: v.tolist() for k, v in coms.items()}}
        if traction_model is None:
            for pname, v in bl["bodies"].items():
                for i, dof in enumerate(("Rx", "Ry", "Rz")):
                    if abs(v["F"][i]) > 0:
                        spec.rigid_loads.append({"type": "rigid_force", "rb": f"rigid_{pname}", "dof": dof,
                                                 "value": float(v["F"][i]), "lc": 1, "name": f"cable_{pname}_{dof}"})
                for i, dof in enumerate(("Ru", "Rv", "Rw")):
                    if abs(v["M_com"][i]) > 0:
                        spec.rigid_loads.append({"type": "rigid_moment", "rb": f"rigid_{pname}", "dof": dof,
                                                 "value": float(v["M_com"][i]), "lc": 1, "name": f"cable_{pname}_{dof}"})
            rep["assumptions"].append("Cable span forces are applied to the rigid petals as fixed-direction resultant force "
                                      "and moment about each petal's centre of mass (directions do not follow petal rotation).")
        else:
            # prescribed traction experiment: pad resultant from petal moment balance, uniform pressure on skin zone
            loads = []
            rep["prescribed_traction"] = {}
            for p in petals:
                surf = tm.surfaces[f"{p.name}_pad_contact"].conn
                Pq = Xs[surf]
                cen = Pq.mean(axis=(0, 1))
                th = np.deg2rad(p.theta_center_deg)
                n_out = np.array([np.cos(th), np.sin(th), 0.0])
                N = CB.petal_pad_force_moment_balance(p, bl["bodies"][p.name], cen, n_out)
                zs = model.surfaces[f"skin_{p.name}_zone"]
                A = _tri_area(model.nodes, zs.conn).sum()
                pr = max(N, 0.0) / A
                loads.append({"surface": f"skin_{p.name}_zone", "value": pr, "lc": 1, "name": f"traction_{p.name}"})
                rep["prescribed_traction"][p.name] = {"pad_force_N": N, "zone_area_mm2": A, "pressure_MPa": pr}
            rep["assumptions"].append("PRESCRIBED TRACTION EXPERIMENT: uniform pressure on each petal's skin zone equal to the "
                                      "pad force from a rigid-petal / frictionless-hinge moment balance. No socket contact is "
                                      "solved; this cannot establish grip.")
            rep["_traction_loads"] = loads
    # --- steps -------------------------------------------------------------------
    sol = cfg["solver"]
    spec.solver = {"max_refs": sol["max_refs"], "dtol": sol["dtol"], "etol": sol["etol"], "rtol": sol["rtol"],
                   "lstol": sol["lstol"],
                   "qn_method": {"type": (("Broyden" if nonsym else "BFGS") if sol["qn_method"] == "auto" else sol["qn_method"]),
                                 "max_ups": sol["max_ups"]},
                   "symmetric_stiffness": 0 if nonsym else 1, "linear_solver": "pardiso", "diverge_reform": 1,
                   "reform_each_time_step": 1}
    spec.time_stepper = {"dtmin": sol["dtmin"], "dtmax": 0.25, "max_retries": sol["max_retries"], "opt_iter": sol["opt_iter"]}
    fr = cb["rating_fractions"]
    # must points: requested load fractions plus every 5% so that a failed run
    # always has a written state within 5% of its last accepted increment
    mp1 = sorted(set([round(f, 6) for f in fr] + [round(0.05 * k, 6) for k in range(1, 21)]))
    st1 = StepSpec("tightening", 1.0, sol["n_steps"], mp1)
    if S.get("enabled") and traction_model is None:
        st1.rigid_bcs.append({"type": "rigid_fixed", "rb": "rigid_cup", "dofs": ["Rz"], "name": "cup_axial_hold"})
    if traction_model is not None:
        st1.loads = rep.pop("_traction_loads")
    spec.steps.append(st1)
    lcs = cfg["load_cases"]
    second = [k for k in ("pullout", "body_weight", "bending", "torsion") if lcs[k].get("enabled")]
    if len(second) > 1:
        raise ValueError("enable at most one of pullout/body_weight/bending/torsion per run")
    if second:
        k = second[0]
        spec.load_curves[2] = {"name": k, "points": [(1.0, 0.0), (2.0, 1.0)]}
        st2 = StepSpec(k, 1.0, sol.get("n_steps_case2", sol["n_steps"]), [round(1.0 + 0.1 * i, 6) for i in range(1, 11)])
        if k in ("pullout", "body_weight"):
            d = lcs[k]["displacement"] * (-1 if k == "pullout" else 1)
            st2.rigid_bcs.append({"type": "rigid_displacement", "rb": "rigid_cup", "dof": "z", "value": d, "lc": 2,
                                  "relative": 1, "name": f"cup_{k}"})
        else:
            # bending: rotation about x through the cup reference point; torsion: about the limb axis
            dof = "Ru" if k == "bending" else "Rw"
            st2.rigid_bcs.append({"type": "rigid_fixed", "rb": "rigid_cup", "dofs": ["Rz"], "name": "cup_axial_hold2"})
            st2.rigid_bcs.append({"type": "rigid_rotation", "rb": "rigid_cup", "dof": dof,
                                  "value": float(np.deg2rad(lcs[k]["rotation"])), "lc": 2, "relative": 1, "name": f"cup_{k}"})
            rep["assumptions"].append(f"{k}: prescribed cup rotation about the cup centre of mass (not a gait load).")
        spec.steps.append(st2)
        # cup Rx,Ry,Ru,Rv,Rw remain fixed globally; for bending/torsion the rotated dof is released below
        if k in ("bending", "torsion"):
            dof = "Ru" if k == "bending" else "Rw"
            for r in spec.rigid_bcs:
                if r.get("name") == "cup_lateral_rot":
                    r["dofs"] = [x for x in r["dofs"] if x != dof]
                    st1.rigid_bcs.append({"type": "rigid_fixed", "rb": "rigid_cup", "dofs": [dof], "name": f"cup_{dof}_hold"})
    if lcs["suction"].get("enabled"):
        sc = lcs["suction"]
        rep["assumptions"].append(f"Suction: {sc['pressure']} MPa sub-atmospheric pressure on the distal skin zone; "
                                  f"justification: {sc['justification']}; sealing: {sc['sealing_assumption']}")
        for st in spec.steps:
            st.loads.append({"surface": "skin_cushion_zone", "value": -sc["pressure"], "lc": 1, "name": "suction"})
    # must-point curves (step-function dtmax)
    nid = max(spec.load_curves) + 1          # FEBio requires consecutive load-controller ids
    for st in spec.steps:
        t0 = st.must_points[0] - (st.must_points[0] - (0 if st is spec.steps[0] else 1.0))
        pts = [(t0, 0.25 * st.duration)] + [(t, 0.25 * st.duration) for t in st.must_points]
        spec.load_curves[nid] = {"name": f"must_points_{st.name}", "points": pts, "interpolate": "STEP"}
        st.mp_lc = nid
        nid += 1
    rep["steps"] = [{"name": s.name, "must_points": s.must_points} for s in spec.steps]
    return spec, rep


def _support(spec, node_set, dofs, name):
    """Displacement supports are written as PRESCRIBED zero displacements:
    FEBio 4.13 logs nodal reaction forces only for prescribed (not fixed)
    degrees of freedom, and the reactions are needed for equilibrium checks."""
    for d in dofs:
        spec.bcs.append({"type": "prescribed displacement", "node_set": node_set, "dof": d, "value": 0.0, "lc": 1,
                         "name": f"{name}_{d}", "support": True})


def _tri_area(X, T):
    return 0.5 * np.linalg.norm(np.cross(X[T[:, 1]] - X[T[:, 0]], X[T[:, 2]] - X[T[:, 0]]), axis=1)
