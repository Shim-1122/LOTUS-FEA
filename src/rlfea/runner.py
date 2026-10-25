"""End-to-end run orchestration with provenance.

prepare -> mesh -> (socket) -> write FEBio input -> run -> parse -> derive ->
export VTK -> checks -> figures -> report.  Partial outputs of failed or
interrupted runs are processed and labelled as such; nothing cached from a
different run is ever shown as this run's solution.
"""
from __future__ import annotations

import json
import os
import pickle
import shutil
import time

import numpy as np

from . import config as C
from . import provenance as P
from .cases import build_case, prepare_geometry
from .post import checks as CK
from .post import results as RS
from .solvers import febio as FB


def prepare(config_path: str, runs_root: str = "runs", run_id: str | None = None) -> P.RunDir:
    raw = C.load(config_path)
    cfg, warns = C.resolve(raw)
    rd = P.RunDir(runs_root, run_id or f"{raw.get('run_name', 'run')}_{time.strftime('%Y%m%dT%H%M%S')}")
    t0 = time.time()
    inputs = {"config": {"path": raw["_config_path"], "sha256": raw["_config_sha256"]},
              "material_set": {"path": cfg["materials"]["set"], "sha256": P.sha256(cfg["materials"]["set"])}}
    g = cfg["geometry"]
    if g["source"] == "files":
        for b in [g["skin"]] + g.get("bones", []) + g.get("ignored_files", []):
            inputs[f"geometry:{b.get('name', 'skin')}"] = {"path": b["path"], "sha256": P.sha256(b["path"]),
                                                          "declared_unit": b.get("unit")}
    shutil.copy(raw["_config_path"], rd.p("config.yaml"))
    P.dump(rd.p("config.resolved.json"), cfg)
    geo = prepare_geometry(cfg)
    P.dump(rd.p("geometry", "audit_source.json"), geo["audit_source"])
    P.dump(rd.p("geometry", "audit_prepared.json"), geo["audit_prepared"])
    P.dump(rd.p("geometry", "assembly_audit.json"), geo["assembly"])
    P.dump(rd.p("geometry", "history.json"), geo["history"])
    P.dump(rd.p("geometry", "anatomy.json"), {"roles": geo["roles"], "missing_anatomy": geo["missing_anatomy"],
                                              "notes": geo["notes"], "frame_transform": geo.get("frame_transform")})
    case = build_case(cfg, geo)
    model, spec = case["model"], case["spec"]
    used = set(spec.block_material.values())
    for k, m in case["material_set"]["materials"].items():
        if k not in used:
            continue
        if m.get("requires_user_calibration"):
            warns.append(f"Material '{k}' REQUIRES USER CALIBRATION (see material set); results are provisional.")
        if m.get("assumed_parameters"):
            warns.append(f"Material '{k}' uses assumed (unsourced) parameters: {', '.join(m['assumed_parameters'])}.")
        if m.get("approximated_parameters"):
            warns.append(f"Material '{k}' uses literature-based APPROXIMATIONS (not measured on this material): "
                         f"{', '.join(m['approximated_parameters'])}.")
    P.dump(rd.p("mesh", "layer_report.json"), case["layer_report"])
    P.dump(rd.p("mesh", "socket_report.json"), case["socket_report"])
    P.dump(rd.p("solver", "load_report.json"), case["load_report"])
    P.dump(rd.p("mesh", "model_summary.json"), model.summary())
    P.dump(rd.p("materials.json"), case["material_set"])
    # mesh export with region ids and initial quality
    q = RS.initial_quality(model)
    grid = RS.to_vtk(model, cell={"radius_ratio": q["radius_ratio"], "min_scaled_jacobian": q["min_scaled_jacobian"]})
    grid.save(rd.p("mesh", "model.vtu"))
    P.dump(rd.p("mesh", "blocks.json"), [{"id": i, "name": b.name, "etype": b.etype, "material": spec.block_material[b.name],
                                          "region": b.region, "formulation": spec.block_domain.get(b.name, {}).get("type") or b.formulation}
                                         for i, b in enumerate(model.blocks)])
    with open(rd.p("mesh", "model.pkl"), "wb") as fh:
        pickle.dump({"model": model, "spec": spec, "case_meta": {k: v for k, v in case.items() if k not in ("model", "spec")}}, fh)
    feb = rd.p("solver", "model.feb")
    info = FB.write_feb(feb, model, spec)
    ver = FB.febio_version()
    manifest = {"run_id": rd.run_id, "run_name": raw.get("run_name"), "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                "inputs": inputs, "repository": P.repo_commit(), "environment": P.environment(),
                "solver": {"backend": "FEBio", **{k: v for k, v in ver.items() if k != "raw"}},
                "solver_input": {"path": "solver/model.feb", "sha256": info["sha256"]},
                "warnings": warns, "load_model": case["load_model"], "prepare_time_s": time.time() - t0,
                "geometry_notes": geo["notes"], "missing_anatomy": geo["missing_anatomy"],
                "patient_geometry": cfg["geometry"]["source"] == "files"}
    P.dump(rd.p("manifest.json"), manifest)
    rd.set_status("prepared")
    return rd


def execute(rd: P.RunDir, threads: int | None = None, timeout: float | None = None, restart: str | None = None,
            overrides: dict | None = None) -> dict:
    """Run FEBio.  With `restart` (a .dmp written at the last converged step)
    the solve resumes from that state; `overrides` (solver tolerances, time
    stepper) are written into the restart file and recorded in the manifest."""
    if restart:
        overrides = overrides or {}
        solver = {k: v for k, v in overrides.items() if k in FB.RESTART_SOLVER_KEYS}
        stepper = {k: v for k, v in overrides.items() if k in FB.RESTART_STEPPER_KEYS}
        control = {k: v for k, v in overrides.items() if k in FB.RESTART_CONTROL_KEYS}
        unknown = set(overrides) - set(solver) - set(stepper) - set(control)
        if unknown:
            raise ValueError(f"unsupported restart override(s): {sorted(unknown)}")
        if restart.endswith(".dmp"):
            dump = os.path.abspath(restart)
            if os.path.dirname(dump) != os.path.abspath(rd.p("solver")):
                raise ValueError("restart dump must be inside the run's solver directory")
            restart = FB.write_restart(rd.p("solver", "restart.feb"), dump, solver, stepper, control)
        man = json.load(open(rd.p("manifest.json")))
        man.setdefault("restarts", []).append({
            "restart_file": os.path.basename(restart), "overrides": overrides,
            "previous_status": P.read_status(rd.path), "time": time.strftime("%Y-%m-%dT%H:%M:%S%z")})
        P.dump(rd.p("manifest.json"), man)
    rd.set_status("running")
    try:
        res = FB.run_febio(rd.p("solver", "model.feb"), rd.p("solver"), threads=threads, timeout=timeout, restart=restart)
    except KeyboardInterrupt:
        rd.set_status("interrupted", note="keyboard interrupt")
        raise
    man = json.load(open(rd.p("manifest.json")))
    man["execution"] = {k: v for k, v in res.items() if k != "convergence"}
    P.dump(rd.p("manifest.json"), man)
    P.dump(rd.p("solver", "convergence.json"), res["convergence"])
    rd.set_status(res["status"], last_converged_time=res["convergence"].get("last_converged_time", 0.0),
                  returncode=res["returncode"], wall_time_s=res["wall_time_s"])
    return res


def _convergence(run_path: str) -> dict:
    """convergence.json, or (if the wrapper itself was killed before writing
    it) a fresh parse of the FEBio log, which is then saved."""
    p = os.path.join(run_path, "solver", "convergence.json")
    if not os.path.exists(p):
        log = os.path.join(run_path, "solver", "model.log")
        conv = FB.parse_log(log) if os.path.exists(log) else {"last_converged_time": 0.0, "steps": []}
        P.dump(p, conv)
    return json.load(open(p))


def postprocess(run_path: str) -> dict:
    """Derive fields, export VTK and compute checks for every output time
    that FEBio actually wrote (converged must points)."""
    with open(os.path.join(run_path, "mesh", "model.pkl"), "rb") as fh:
        pk = pickle.load(fh)
    model, spec = pk["model"], pk["spec"]
    status = P.read_status(run_path)
    cfg = json.load(open(os.path.join(run_path, "config.resolved.json")))
    raw = RS.load_raw(os.path.join(run_path, "solver"), model, spec.log_surfaces)
    from .solvers.febio import spec_surfaces
    raw["_all_surfaces"] = list(spec_surfaces(model, spec).keys())
    contact_source = RS.attach_true_contact(raw, os.path.join(run_path, "solver"), model, spec.log_surfaces)
    V0 = RS.element_volumes(model)
    offs, slices = 0, {}
    for b in model.blocks:
        slices[b.name] = (offs, offs + len(b.conn)); offs += len(b.conn)
    rigid_ids = CK.rigid_table(spec)
    mu = cfg["contact"]["friction"]
    q = RS.initial_quality(model)
    region_id = np.concatenate([np.full(len(b.conn), i) for i, b in enumerate(model.blocks)])
    summary = {"status": status, "times": raw["times"], "states": [], "contact_source": contact_source}
    pvd = []
    skin_tris = np.vstack([model.surfaces[s].conn for s in spec.log_surfaces]) if spec.log_surfaces else np.zeros((0, 3), int)
    for k, st in enumerate(raw["states"]):
        t = st["time"]
        rec = {"time": t}
        U = st["U"]
        cell = {"radius_ratio": q["radius_ratio"], "min_scaled_jacobian_initial": q["min_scaled_jacobian"], "region_id": region_id}
        if "F" in st and "S" in st:
            ef = RS.element_fields(st["F"], st["S"])
            cell.update(ef)
            cell["sed"] = st.get("sed")
            rec["deformation"] = CK.deformation_checks(ef["J"], ef["stretch_max"], ef["stretch_min"])
            rec["strain_energy_Nmm"] = CK.strain_energy(st["sed"], V0, slices)
            # extremes per region (tissue only)
            rec["region_extremes"] = {}
            for b in model.blocks:
                a, e = slices[b.name]
                if b.formulation == "rigid" or spec.materials[spec.block_material[b.name]].get("rigid"):
                    continue
                rec["region_extremes"][b.name] = {
                    "max_disp_mm": float(np.nanmax(np.linalg.norm(U[np.unique(b.conn)], axis=1))),
                    "max_logstrain": float(np.nanmax(ef["logstrain_max"][a:e])), "min_logstrain": float(np.nanmin(ef["logstrain_min"][a:e])),
                    "max_pressure_MPa": float(np.nanmax(ef["pressure"][a:e])), "min_pressure_MPa": float(np.nanmin(ef["pressure"][a:e])),
                    "max_von_mises_MPa": float(np.nanmax(ef["von_mises"][a:e])), "min_J": float(np.nanmin(ef["J"][a:e]))}
        tissue_nodes = np.unique(np.concatenate([b.conn.ravel() for b in model.blocks if b.region in ("skin", "fat", "muscle") or b.name == "proximal_cap"]))
        rec["max_tissue_displacement_mm"] = float(np.nanmax(np.linalg.norm(U[tissue_nodes], axis=1)))
        rec["equilibrium"] = CK.equilibrium(model, spec, st, rigid_ids)
        grid = RS.to_vtk(model, U=np.nan_to_num(U), cell={k: v for k, v in cell.items() if v is not None})
        fn = f"volume_t{t:.4f}.vtu"
        grid.save(os.path.join(run_path, "results", fn))
        pvd.append((t, fn))
        # contact surface
        if st["faces"]:
            parts, zone = [], []
            for zi, s in enumerate(spec.log_surfaces):
                if s in st["faces"]:
                    tc = st.get("true_contact", {}).get(s)
                    cf = RS.contact_fields(model.nodes, U, model.surfaces[s].conn, st["faces"][s], mu,
                                           tc["traction"] if tc else None, tc["stick"] if tc else None)
                    cf["zone"] = np.full(len(cf["contact_pressure"]), zi)
                    parts.append(cf)
            cat = {kk: np.concatenate([p_[kk] for p_ in parts]) for kk in parts[0]}
            import pyvista as pv
            faces = np.c_[np.full(len(skin_tris), 3), skin_tris].ravel()
            surf = pv.PolyData(model.nodes.copy(), faces)
            surf.point_data["displacement"] = np.nan_to_num(U)
            for kk, v in cat.items():
                surf.cell_data[kk] = v
            surf.save(os.path.join(run_path, "results", f"contact_t{t:.4f}.vtp"))
            A = cat["area_deformed"]
            act = cat["contact_pressure"] > 1e-9
            rec["contact"] = {"contact_area_mm2": float(A[act].sum()), "max_pressure_MPa": float(cat["contact_pressure"].max()),
                              "mean_pressure_MPa": float((cat["contact_pressure"] * A).sum() / max(A[act].sum(), 1e-12)),
                              "max_tangential_MPa": float(cat["traction_tangential"].max()),
                              "slip_fraction_of_contact_area": float(A[cat["stick_slip_status"] == 2].sum() / max(A[act].sum(), 1e-12)),
                              "max_penetration_mm": float(np.nanmax(cat["contact_gap"])),
                              "normal_force_N": float((cat["contact_pressure"] * A).sum()),
                              "tangential_force_vector_N": (cat["tangential_vector"] * A[:, None]).sum(0).tolist(),
                              "traction_resultant_N": (cat["traction_vector"] * A[:, None]).sum(0).tolist(),
                              "per_zone": {}}
            for zi, s in enumerate(spec.log_surfaces):
                m_ = cat["zone"] == zi
                rec["contact"]["per_zone"][s] = {"normal_force_N": float((cat["contact_pressure"][m_] * A[m_]).sum()),
                                                 "traction_resultant_N": (cat["traction_vector"][m_] * A[m_, None]).sum(0).tolist(),
                                                 "max_pressure_MPa": float(cat["contact_pressure"][m_].max()) if m_.any() else 0.0,
                                                 "slip_area_fraction": float(A[m_ & (cat["stick_slip_status"] == 2)].sum() / max(A[m_ & act].sum(), 1e-12))}
        if "connectors" in st:
            names = [c["name"] for c in spec.connectors]
            rec["connectors"] = {names[int(i) - 1] if int(i) - 1 < len(names) else str(i): d.tolist()
                                 for i, d in zip(st["connectors"]["ids"], st["connectors"]["data"])}
        # per-petal free body: cable (applied) + hinge + stabiliser + pad contact (= -skin traction resultant)
        try:
            lr = json.load(open(os.path.join(run_path, "solver", "load_report.json")))
            if lr.get("cable") and "connectors" in rec and "contact" in rec:
                fbd = {}
                lam = min(t, 1.0)
                for pname, bl in lr["cable"]["body_loads"].items():
                    Fc = lam * np.array(bl["F_N"])
                    H = -np.array(rec["connectors"].get(f"hinge_{pname}", [0, 0, 0])[:3])          # connector force reported on body_a (cup)
                    Sk = np.array(rec["connectors"].get(f"stabilizer_{pname}", [0, 0, 0])[:3])
                    Pc = -np.array(rec["contact"]["per_zone"].get(f"skin_{pname}_zone", {}).get("traction_resultant_N", [0, 0, 0]))
                    r = Fc + H + Sk + Pc
                    fbd[pname] = {"cable_N": Fc.tolist(), "hinge_on_petal_N": H.tolist(), "stabiliser_N": Sk.tolist(),
                                  "pad_contact_N": Pc.tolist(), "residual_N": r.tolist(),
                                  "residual_over_cable": float(np.linalg.norm(r) / max(np.linalg.norm(Fc), 1e-12)),
                                  "stabiliser_share_of_cable": float(np.linalg.norm(Sk) / max(np.linalg.norm(Fc), 1e-12))}
                rec["petal_free_body"] = fbd
        except Exception as e:  # noqa: BLE001
            rec["petal_free_body_error"] = str(e)
        summary["states"].append(rec)
    with open(os.path.join(run_path, "results", "collection.pvd"), "w") as fh:
        fh.write('<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1">\n<Collection>\n')
        for t, f in pvd:
            fh.write(f'<DataSet timestep="{t}" part="0" file="{f}"/>\n')
        fh.write("</Collection>\n</VTKFile>\n")
    P.dump(os.path.join(run_path, "results", "summary.json"), summary)
    return summary


def make_figures(run_path: str, exaggeration: float | None = 5.0) -> dict:
    """Render report figures from the solved fields of THIS run."""
    from .post import figures as FG
    with open(os.path.join(run_path, "mesh", "model.pkl"), "rb") as fh:
        pk = pickle.load(fh)
    model, spec = pk["model"], pk["spec"]
    man = json.load(open(os.path.join(run_path, "manifest.json")))
    st = P.read_status(run_path)
    cfg = json.load(open(os.path.join(run_path, "config.resolved.json")))
    summ = json.load(open(os.path.join(run_path, "results", "summary.json")))
    load = json.load(open(os.path.join(run_path, "solver", "load_report.json")))
    figdir = os.path.join(run_path, "figures")
    idx = {}
    T = (load.get("tension") or {}).get("applied_N") or 0.0
    steps = [s.name for s in spec.steps]

    def load_text(t):
        if t <= 1.0 + 1e-9:
            return f"strand tension {T * t:.3g} N ({load['tension']['basis']})"
        return f"{steps[1] if len(steps) > 1 else ''} {100 * (t - 1):.0f}% + tension {T:.3g} N"

    meta = {"run_name": man.get("run_name"), "case": " + ".join(steps), "solver_version": man["solver"].get("version"),
            "state": st.get("state"), "load_text": load_text, "fiber_ref": cfg["layers"]["skin_fiber_reference"],
            "skin_gamma_deg": spec.materials.get("skin", {}).get("params", {}).get("gamma", 41.0)}
    res = os.path.join(run_path, "results")
    times = summ["times"]
    if not times:
        return idx
    vtu0 = os.path.join(res, f"volume_t{times[0]:.4f}.vtu")
    t_end1 = max([t for t in times if t <= 1.0 + 1e-9] or [times[-1]])
    vtu1 = os.path.join(res, f"volume_t{t_end1:.4f}.vtu")
    meta["time"] = t_end1
    idx["full_limb"] = FG.full_limb_views(run_path, vtu1, meta, os.path.join(figdir, "full_limb_1x.png"))
    if exaggeration:
        idx["full_limb_exaggerated"] = FG.full_limb_views(run_path, vtu1, meta, os.path.join(figdir, f"full_limb_x{exaggeration:g}.png"), scale=exaggeration)
    idx["sections"] = FG.section_fields(run_path, vtu1, meta, os.path.join(figdir, "sections.png"),
                                        [("displacement_magnitude", "|u|", "mm", "viridis"),
                                         ("logstrain_max", "max principal log strain", "-", "plasma"),
                                         ("pressure", "hydrostatic pressure p = -tr(sigma)/3", "MPa", "coolwarm"),
                                         ("sigma_1", "max principal Cauchy stress", "MPa", "coolwarm")])
    vtp = os.path.join(res, f"contact_t{t_end1:.4f}.vtp")
    if os.path.exists(vtp):
        idx["contact"] = FG.contact_maps(vtp, meta, os.path.join(figdir, "contact.png"), cfg["contact"]["friction"])
    idx["regions"] = FG.regions_and_fibers(run_path, vtu0, meta, os.path.join(figdir, "regions_fibers.png"), model)
    idx["quality"] = FG.mesh_quality(run_path, os.path.join(run_path, "mesh", "model.vtu"), meta, os.path.join(figdir, "mesh_quality.png"))
    conv = _convergence(run_path)
    idx["convergence"] = FG.convergence_plot(conv, os.path.join(figdir, "convergence.png"), title=f"{man.get('run_name')} | state {st.get('state')}")
    if load.get("cable"):
        petals = [{"name": p.name, "theta_center_deg": p.theta_center_deg, "half_width_deg": p.half_width_deg,
                   "r_ref": float(np.linalg.norm(model.nodes[p.shell_outer_nodes[:, 0]][:, :2], axis=1).mean())}
                  for p in model.meta.get("petals", [])] if model.meta.get("petals") else []
        idx["free_body"] = FG.free_body_diagram(load, petals, os.path.join(figdir, "free_body.png"), title=man.get("run_name", ""))
    # force-displacement / load-response
    S = summ["states"]
    t1 = [s for s in S if s["time"] <= 1.0 + 1e-9]
    ser = [{"x": [T * s["time"] for s in t1], "y": [s["max_tissue_displacement_mm"] for s in t1], "label": "max tissue |u| [mm]"}]
    if t1 and "contact" in t1[-1]:
        ser.append({"x": [T * s["time"] for s in t1], "y": [s.get("contact", {}).get("normal_force_N", 0) / 10 for s in t1],
                    "label": "total contact normal force / 10 [N]", "style": "s--"})
    idx["force_displacement"] = FG.xy_plot(ser, "strand tension [N]", "response", os.path.join(figdir, "force_displacement.png"),
                                           title=f"{man.get('run_name')} tightening | state {st.get('state')} | true values")
    t2 = [s for s in S if s["time"] > 1.0 + 1e-9]
    if t2 and len(steps) > 1 and steps[1] in ("pullout", "body_weight"):
        d = cfg["load_cases"][steps[1]]["displacement"]
        cup = []
        for s in [x for x in S if x["time"] >= 1.0 - 1e-9]:
            rb = s["equilibrium"]["rigid_bodies"].get("rigid_cup", {})
            cup.append((d * (s["time"] - 1.0), rb.get("F", [0, 0, 0])[2]))
        idx["pullout"] = FG.xy_plot([{"x": [c[0] for c in cup], "y": [c[1] for c in cup], "label": "cup axial reaction Fz [N]"}],
                                    f"prescribed socket displacement ({'distal' if steps[1] == 'pullout' else 'proximal'}) [mm]", "force [N]",
                                    os.path.join(figdir, "pullout.png"), title=f"{steps[1]} | state {st.get('state')}")
    rel = {k: os.path.relpath(v, os.path.join(run_path, "report")) for k, v in idx.items()}
    P.dump(os.path.join(figdir, "index.json"), rel)
    return rel
