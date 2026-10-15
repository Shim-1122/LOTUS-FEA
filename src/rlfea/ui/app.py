"""Interactive results/job interface (trame + PyVista, runs in a browser).

    rlfea gui --runs runs --port 8080

Features: run browser with job state; launch a configuration as a background
job; live solver log tail and convergence table; full-limb view of the ACTUAL
solver fields for any written output time; true-scale (1x) deformation by
default with an explicit, labelled exaggeration control; tissue opacity, bone
and socket visibility, mesh edges, clip plane; field selection (displacement
magnitude/components, strains, principal stresses, pressure, J, quality,
regions, contact pressure/traction/stick-slip); node probing with physical
coordinates and units; load-step animation by slider.
"""
from __future__ import annotations

import glob
import json
import os
import subprocess
import sys

import numpy as np

POINT_FIELDS = {"|u| (mm)": ("displacement_magnitude", None), "u_x (mm)": ("displacement", 0),
                "u_y (mm)": ("displacement", 1), "u_z (mm)": ("displacement", 2)}
CELL_FIELDS = {"J (-)": "J", "max principal log strain (-)": "logstrain_max", "min principal log strain (-)": "logstrain_min",
               "max principal Green-Lagrange (-)": "GL_max", "max principal Cauchy stress (MPa)": "sigma_1",
               "min principal Cauchy stress (MPa)": "sigma_3", "hydrostatic pressure (MPa)": "pressure",
               "von Mises (MPa)": "von_mises", "strain energy density (MPa)": "sed",
               "tet radius ratio (-)": "radius_ratio", "min scaled Jacobian, initial (-)": "min_scaled_jacobian_initial",
               "material region (block id)": "block_id"}
CONTACT_FIELDS = {"contact pressure (MPa)": "contact_pressure", "tangential traction (MPa)": "traction_tangential",
                  "stick/slip 0 open 1 stick 2 slip (FEBio flag)": "stick_slip_status", "contact gap (mm)": "contact_gap",
                  "friction utilisation (-)": "friction_utilisation"}


def list_runs(root):
    from ..provenance import read_status
    out = []
    for d in sorted(glob.glob(os.path.join(root, "*", "manifest.json"))):
        rp = os.path.dirname(d)
        st = read_status(rp)
        out.append({"path": rp, "name": os.path.basename(rp), "state": st.get("state")})
    return out


def main(runs_root="runs", port=8080):
    import pyvista as pv
    from pyvista.trame.ui import plotter_ui
    from trame.app import get_server
    from trame.ui.vuetify3 import SinglePageWithDrawerLayout
    from trame.widgets import html
    from trame.widgets import vuetify3 as v3

    pv.OFF_SCREEN = True
    server = get_server(client_type="vue3")
    state, ctrl = server.state, server.controller
    pl = pv.Plotter()
    pl.set_background("white")
    cache = {}

    runs = list_runs(runs_root)
    state.runs = [f"{r['name']}  [{r['state']}]" for r in runs]
    with_results = [i for i, r in enumerate(runs) if os.path.exists(os.path.join(r["path"], "results", "summary.json"))]
    state.run_sel = state.runs[with_results[-1]] if with_results else (state.runs[-1] if runs else None)
    state.times = []
    state.t_index = 0
    state.field = "|u| (mm)"
    state.fields = list(POINT_FIELDS) + list(CELL_FIELDS) + list(CONTACT_FIELDS)
    state.scale = 1.0
    state.opacity = 0.35
    state.show_bones = True
    state.show_socket = False
    state.edges = False
    state.clip = "none"
    state.clip_pos = 0.0
    state.banner = ""
    state.probe_node = ""
    state.probe_out = ""
    state.log_tail = ""
    state.config_path = "configs/synthetic_human.yaml"
    state.job_msg = ""

    def run_path():
        if not state.run_sel:
            return None
        name = state.run_sel.split("  [")[0]
        return os.path.join(runs_root, name)

    def load_run():
        rp = run_path()
        if rp is None:
            return
        summ_p = os.path.join(rp, "results", "summary.json")
        state.times = json.load(open(summ_p))["times"] if os.path.exists(summ_p) else []
        state.t_index = len(state.times) - 1 if state.times else 0
        cache.clear()
        cache["_camera_set"] = False
        cache["blocks"] = json.load(open(os.path.join(rp, "mesh", "blocks.json")))
        update()

    def grids(rp, t):
        key = (rp, t)
        if key not in cache:
            vtu = os.path.join(rp, "results", f"volume_t{t:.4f}.vtu")
            vtp = os.path.join(rp, "results", f"contact_t{t:.4f}.vtp")
            cache[key] = (pv.read(vtu), pv.read(vtp) if os.path.exists(vtp) else None)
        return cache[key]

    def update(**_):
        pl.clear()
        rp = run_path()
        if rp is None or not state.times:
            pl.add_text("No processed results for this run (prepared/running/failed before first output?)", font_size=10)
            ctrl.view_update()
            return
        from ..provenance import read_status
        st = read_status(rp)
        man = json.load(open(os.path.join(rp, "manifest.json")))
        t = state.times[int(state.t_index)]
        g, s = grids(rp, t)
        sc = float(state.scale)
        gw = g.warp_by_vector("displacement", factor=sc)
        blocks = cache["blocks"]
        tis = [b["id"] for b in blocks if b["name"] in ("skin", "fat", "muscle", "proximal_cap")]
        bone = [b["id"] for b in blocks if b["region"].startswith("bone")]
        sock = [b["id"] for b in blocks if b["region"].startswith("socket")]

        def sub(G, ids):
            return G.extract_cells(np.where(np.isin(G.cell_data["block_id"], ids))[0])

        T = sub(gw, tis)
        if state.clip != "none":
            T = T.clip(normal=state.clip, origin=[state.clip_pos if a == state.clip else 0 for a in "xyz"], invert=True, crinkle=True)
        f = state.field
        if f in POINT_FIELDS:
            name, comp = POINT_FIELDS[f]
            arr = T.point_data[name] if comp is None else T.point_data[name][:, comp]
            T.point_data["_show"] = arr
            pl.add_mesh(T, scalars="_show", opacity=1.0 if state.clip != "none" else float(state.opacity), show_edges=state.edges,
                        cmap="viridis", scalar_bar_args={"title": f})
        elif f in CELL_FIELDS:
            pl.add_mesh(T, scalars=CELL_FIELDS[f], preference="cell", opacity=1.0 if state.clip != "none" else float(state.opacity),
                        show_edges=state.edges, cmap="coolwarm" if "stress" in f or "pressure" in f else "plasma", scalar_bar_args={"title": f})
        else:
            if s is not None:
                sw = s.warp_by_vector("displacement", factor=sc)
                pl.add_mesh(sw, scalars=CONTACT_FIELDS[f], preference="cell", cmap="inferno", scalar_bar_args={"title": f})
            pl.add_mesh(T.extract_surface(), color="#e8c4a8", opacity=0.15)
        if state.show_bones and bone:
            pl.add_mesh(sub(gw, bone).extract_surface(), color="#f2efe6")
        if state.show_socket and sock:
            pl.add_mesh(sub(gw, sock).extract_surface(), color="#5a6b7c", opacity=0.3)
        lbl = "deformation scale 1x (TRUE SCALE)" if sc == 1 else f"EXAGGERATED deformation x{sc:g}"
        state.banner = f"{man.get('run_name')} | state: {st.get('state')} | FEBio {man['solver'].get('version')} | t = {t:g} | {lbl}"
        pl.add_text(state.banner, font_size=9, position="upper_left", color="red" if sc != 1 else "black")
        if not cache.get("_camera_set"):
            pl.view_vector((1.0, -1.0, 0.35))
            pl.reset_camera()
            cache["_camera_set"] = True
        ctrl.view_update()

    def probe(**_):
        rp = run_path()
        if rp is None or not state.times:
            return
        g, _ = grids(rp, state.times[int(state.t_index)])
        txt = str(state.probe_node).strip()
        try:
            if "," in txt:
                xyz = np.array([float(v) for v in txt.split(",")])
                nid = int(np.argmin(np.linalg.norm(g.points - xyz, axis=1)))
            else:
                nid = int(txt) - 1
            X = g.points[nid]; U = g.point_data["displacement"][nid]
            cells = g.extract_points([nid], adjacent_cells=True)
            vals = {k: float(np.nanmean(cells.cell_data[k])) for k in ("J", "logstrain_max", "pressure", "sigma_1") if k in cells.cell_data}
            state.probe_out = (f"node {nid + 1}: X = ({X[0]:.3f}, {X[1]:.3f}, {X[2]:.3f}) mm | u = ({U[0]:.4f}, {U[1]:.4f}, {U[2]:.4f}) mm, "
                               f"|u| = {np.linalg.norm(U):.4f} mm | adjacent-element means: " + ", ".join(f"{k}={v:.4g}" for k, v in vals.items())
                               + " (units: J,-; log strain,-; MPa)")
        except Exception as e:  # noqa: BLE001
            state.probe_out = f"probe error: {e}"

    def refresh_log(**_):
        rp = run_path()
        if rp is None:
            return
        p = os.path.join(rp, "solver", "model.log")
        if os.path.exists(p):
            lines = open(p, errors="replace").read().splitlines()
            keep = [l for l in lines if any(k in l for k in ("beginning time step", "converged at", "residual", "Retrying", "ERROR", "T E R M"))]
            state.log_tail = "\n".join(keep[-40:])
        conv = os.path.join(rp, "solver", "convergence.json")
        if os.path.exists(conv):
            c = json.load(open(conv))
            state.log_tail += f"\n\nparsed: {len(c.get('time_steps', []))} attempts, last converged t = {c.get('last_converged_time')}"

    def launch(**_):
        cfgp = state.config_path
        if not os.path.exists(cfgp):
            state.job_msg = f"config not found: {cfgp}"
            return
        subprocess.Popen([sys.executable, "-m", "rlfea.cli", "run", cfgp, "--runs", runs_root], stdout=subprocess.DEVNULL,
                         stderr=subprocess.DEVNULL, start_new_session=True)
        state.job_msg = f"launched {cfgp}; use Refresh runs to follow the job state"

    def refresh_runs(**_):
        r = list_runs(runs_root)
        state.runs = [f"{x['name']}  [{x['state']}]" for x in r]

    state.change("run_sel")(lambda **_: load_run())
    for k in ("t_index", "field", "scale", "opacity", "show_bones", "show_socket", "edges", "clip", "clip_pos"):
        state.change(k)(update)

    with SinglePageWithDrawerLayout(server) as layout:
        layout.title.set_text("Residual-limb FEA viewer (numerical results; not clinical validation)")
        with layout.drawer:
            v3.VSelect(label="Run", v_model=("run_sel",), items=("runs",), density="compact")
            v3.VBtn("Refresh runs", click=refresh_runs, size="small")
            v3.VTextField(label="Config to run", v_model=("config_path",), density="compact")
            v3.VBtn("Prepare + solve (background job)", click=launch, size="small")
            html.Div("{{ job_msg }}", style="font-size:11px")
            v3.VSlider(label="output time", v_model=("t_index",), min=0, max=("times.length - 1",), step=1, density="compact")
            html.Div("t = {{ times[t_index] }}", style="font-size:11px")
            v3.VSelect(label="Field", v_model=("field",), items=("fields",), density="compact")
            v3.VTextField(label="Deformation scale (1 = true)", v_model=("scale",), type="number", density="compact")
            html.Div("{{ scale == 1 ? 'TRUE SCALE (1x)' : 'EXAGGERATED x' + scale }}", style="font-weight:bold; color:{{ scale == 1 ? 'black' : 'red' }}")
            v3.VSlider(label="tissue opacity", v_model=("opacity",), min=0.05, max=1.0, step=0.05, density="compact")
            v3.VCheckbox(label="bones", v_model=("show_bones",), density="compact")
            v3.VCheckbox(label="socket", v_model=("show_socket",), density="compact")
            v3.VCheckbox(label="mesh edges", v_model=("edges",), density="compact")
            v3.VSelect(label="clip plane normal", v_model=("clip",), items=(["none", "x", "y", "z"],), density="compact")
            v3.VSlider(label="clip position (mm)", v_model=("clip_pos",), min=-120, max=250, step=1, density="compact")
            v3.VTextField(label="Probe: node id or x,y,z (mm)", v_model=("probe_node",), density="compact")
            v3.VBtn("Probe", click=probe, size="small")
            html.Div("{{ probe_out }}", style="font-size:11px")
            v3.VBtn("Solver log / convergence", click=refresh_log, size="small")
            html.Pre("{{ log_tail }}", style="font-size:10px; max-height:260px; overflow:auto")
        with layout.content:
            with v3.VContainer(fluid=True, classes="pa-0 fill-height"):
                view = plotter_ui(pl)
                ctrl.view_update = view.update
    load_run()
    server.start(port=port, open_browser=False)


if __name__ == "__main__":
    main(*(sys.argv[1:2] or ["runs"]))
