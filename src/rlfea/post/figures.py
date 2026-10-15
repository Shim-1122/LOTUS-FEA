"""Publication/report figures rendered from the actual solver fields.

Rules enforced here:
* default deformation scale is 1x and the scale is printed on every image;
  exaggerated views are generated only on request and say so in red;
* element fields are drawn as cell data (flat per element, no smoothing or
  interpolation across material regions);
* every image carries load case, time/load level, solver and run state;
* comparisons use a shared, explicit colour range.
"""
from __future__ import annotations

import json
import os

import warnings

import numpy as np

warnings.filterwarnings("ignore", message=".*extract_surface.*")
TISSUE = ("skin", "fat", "muscle", "proximal_cap")


def _pv():
    import pyvista as pv
    pv.OFF_SCREEN = True
    return pv


def _blocks(run_path):
    return json.load(open(os.path.join(run_path, "mesh", "blocks.json")))


def _ids(blocks, pred):
    return [b["id"] for b in blocks if pred(b)]


def _subset(grid, ids):
    m = np.isin(grid.cell_data["block_id"], ids)
    return grid.extract_cells(np.where(m)[0])


def label(meta: dict, t: float, scale: float) -> str:
    sc = "deformation scale 1x (true scale)" if scale == 1 else f"EXAGGERATED deformation x{scale:g}"
    return (f"{meta.get('run_name', '')} | {meta.get('case', '')} | t = {t:.3g} "
            f"({meta.get('load_text', lambda t: '')(t)}) | FEBio {meta.get('solver_version', '?')} | "
            f"state: {meta.get('state', '?')} | {sc}")


def full_limb_views(run_path, vtu_t, meta, out_png, field="displacement_magnitude", clim=None, scale=1.0,
                    show_socket=True, units="mm"):
    """Undeformed (left) and deformed (right) whole limb; skin translucent,
    bones opaque, plus a sagittal section of the deformed limb (right)."""
    pv = _pv()
    blocks = _blocks(run_path)
    g = pv.read(vtu_t)
    g0 = g.copy()
    gw = g.warp_by_vector("displacement", factor=scale)
    tis = _ids(blocks, lambda b: b["name"] in TISSUE)
    bone = _ids(blocks, lambda b: b["region"].startswith("bone"))
    sock = _ids(blocks, lambda b: b["region"].startswith("socket"))
    t = meta["time"]
    if clim is None:
        clim = [0.0, float(np.nanmax(g.point_data["displacement_magnitude"]))]
    p = pv.Plotter(off_screen=True, shape=(1, 3), window_size=(1800, 820), border=False)
    for col, (G, ttl) in enumerate([(g0, "Undeformed"), (gw, f"Deformed ({'1x' if scale == 1 else f'x{scale:g} EXAGGERATED'})")]):
        p.subplot(0, col)
        p.add_text(ttl, font_size=11, color="red" if (col == 1 and scale != 1) else "black")
        T = _subset(G, tis)
        if col == 0:
            p.add_mesh(T.extract_surface(), color="#e8c4a8", opacity=0.25)
        else:
            p.add_mesh(T.extract_surface(), scalars=field, clim=clim, cmap="viridis", opacity=0.35,
                       scalar_bar_args={"title": f"|u| [{units}]", "vertical": True, "position_x": 0.85})
        p.add_mesh(_subset(G, bone).extract_surface(), color="#f2efe6", opacity=1.0)
        if show_socket and sock:
            p.add_mesh(_subset(G, sock).extract_surface(), color="#5a6b7c", opacity=0.18)
        p.view_vector((1.0, -1.0, 0.35)); p.camera.zoom(1.25)
    p.subplot(0, 2)
    p.add_text("Sagittal section (deformed, cell data)", font_size=11)
    T = _subset(gw, tis + bone)
    cut = T.slice(normal="x", origin=(0, 0, 0))
    p.add_mesh(cut, scalars=field, clim=clim, cmap="viridis", show_scalar_bar=False, lighting=False)
    p.view_vector((1.0, 0.0, 0.0)); p.camera.zoom(1.2)
    p.add_text(label(meta, t, scale), position="lower_left", font_size=11, viewport=False)
    p.screenshot(out_png)
    p.close()
    return out_png


def section_fields(run_path, vtu_t, meta, out_png, fields, normal="x", scale=1.0, clims=None):
    """Grid of sections through the deformed full limb, one per field.
    fields: list of (array_name, title, unit, cmap)."""
    pv = _pv()
    blocks = _blocks(run_path)
    g = pv.read(vtu_t).warp_by_vector("displacement", factor=scale)
    tis = _ids(blocks, lambda b: b["name"] in TISSUE)
    bone = _ids(blocks, lambda b: b["region"].startswith("bone"))
    T = _subset(g, tis)
    B = _subset(g, bone)
    n = len(fields)
    p = pv.Plotter(off_screen=True, shape=(1, n), window_size=(560 * n, 820), border=False)
    for i, (name, title, unit, cmap) in enumerate(fields):
        p.subplot(0, i)
        cut = T.slice(normal=normal, origin=(0, 0, 0))
        arr = cut.cell_data[name] if name in cut.cell_data else cut.point_data[name]
        cl = (clims or {}).get(name) or [float(np.nanpercentile(arr, 0.5)), float(np.nanpercentile(arr, 99.5))]
        if cl[0] < 0 < cl[1] and cmap == "coolwarm":
            a = max(abs(cl[0]), abs(cl[1])); cl = [-a, a]
        p.add_mesh(cut, scalars=name, clim=cl, cmap=cmap, preference="cell", lighting=False,
                   scalar_bar_args={"title": f"{title} [{unit}]", "vertical": False, "position_y": 0.02, "width": 0.8, "position_x": 0.1})
        p.add_mesh(B.slice(normal=normal, origin=(0, 0, 0)), color="#d9d4c7", lighting=False)
        p.add_text(title, font_size=10)
        p.view_vector((1.0, 0, 0) if normal == "x" else (0, -1.0, 0)); p.camera.zoom(1.25)
    p.add_text(label(meta, meta["time"], scale) + " | percentile colour range shown on bar; peaks may be mesh-dependent",
               position="lower_left", font_size=11, viewport=False)
    p.screenshot(out_png)
    p.close()
    return out_png


def contact_maps(vtp_t, meta, out_png, mu, scale=1.0):
    pv = _pv()
    s = pv.read(vtp_t)
    s = s.warp_by_vector("displacement", factor=scale)
    s.cell_data["contact_pressure_kPa"] = s.cell_data["contact_pressure"] * 1000
    s.cell_data["traction_tangential_kPa"] = s.cell_data["traction_tangential"] * 1000
    panels = [("contact_pressure_kPa", "Contact pressure", "kPa", "inferno"),
              ("traction_tangential_kPa", "Tangential traction", "kPa", "magma"),
              ("stick_slip_status", "0 open / 1 stick / 2 slip (derived)", "-", "Set1"),
              ("contact_gap", "Gap (FEBio sign: >0 penetration)", "mm", "coolwarm")]
    p = pv.Plotter(off_screen=True, shape=(2, 4), window_size=(2000, 1250), border=False)
    for r, view in enumerate([(1.0, -1.0, 0.2), (-1.0, 1.0, 0.2)]):
        for c, (name, title, unit, cmap) in enumerate(panels):
            p.subplot(r, c)
            arr = s.cell_data[name]
            if name == "stick_slip_status":
                cl = [0, 2]
            elif name == "contact_gap":
                a = float(np.nanmax(np.abs(arr))) or 1e-6; cl = [-a, a]
            else:
                cl = [0.0, float(np.nanpercentile(arr, 99.5)) or 1e-9]
                if r == 0:
                    p.add_text(f"max {np.nanmax(arr):.3g} {unit} (bar clipped at 99.5th pct)", position="upper_right", font_size=8)
            p.add_mesh(s, scalars=name, clim=cl, cmap=cmap, preference="cell", lighting=False, n_colors=3 if name == "stick_slip_status" else 256,
                       scalar_bar_args={"title": f"{title} [{unit}]", "vertical": False, "position_y": 0.02, "width": 0.8, "position_x": 0.1})
            if r == 0:
                p.add_text(title, font_size=9)
            p.view_vector(view); p.camera.zoom(1.3)
    p.add_text(label(meta, meta["time"], scale) + f" | friction mu = {mu} | slip status derived from FEBio tractions (|t_t| >= 0.98 mu p)",
               position="lower_left", font_size=11, viewport=False)
    p.screenshot(out_png)
    p.close()
    return out_png


def regions_and_fibers(run_path, vtu0, meta, out_png, model):
    """Material regions on a section and skin fibre families on the surface."""
    pv = _pv()
    blocks = _blocks(run_path)
    g = pv.read(vtu0)
    tis = _ids(blocks, lambda b: b["name"] in TISSUE or b["region"].startswith("bone"))
    T = _subset(g, tis)
    cut = T.clip(normal="x", origin=(0, 0, 0), invert=True, crinkle=True)
    p = pv.Plotter(off_screen=True, shape=(1, 3), window_size=(1800, 820), border=False)
    p.subplot(0, 0)
    p.add_mesh(cut, scalars="block_id", cmap="tab10", preference="cell", show_edges=False,
               scalar_bar_args={"title": "region (block id; see table)", "vertical": False})
    p.add_text("Material regions (sagittal section)", font_size=10)
    p.view_vector((1, 0, 0)); p.camera.zoom(1.25)
    p.subplot(0, 1)
    zoom = cut.clip_box(bounds=_zoom_box(cut), invert=False)
    p.add_mesh(zoom, scalars="block_id", cmap="tab10", preference="cell", show_edges=True, show_scalar_bar=False)
    p.add_text("Zoom: skin / fat wedges and muscle tets", font_size=10)
    p.view_vector((1, 0, 0)); p.camera.zoom(1.1)
    p.subplot(0, 2)
    skin = model.block("skin")
    gam = meta.get("skin_gamma_deg", 41.0)
    X = model.nodes
    c = X[skin.conn].mean(axis=1)
    a, d = skin.mat_axis[:, 0], skin.mat_axis[:, 1]
    sel = np.arange(0, len(c), max(1, len(c) // 900))
    gr = np.deg2rad(gam)
    f1 = np.cos(gr) * a + np.sin(gr) * d
    f2 = np.cos(gr) * a - np.sin(gr) * d
    surf = _subset(g, _ids(blocks, lambda b: b["name"] == "skin")).extract_surface()
    p.add_mesh(surf, color="#e8c4a8", opacity=0.6)
    for f, col in ((f1, "#1f4e9c"), (f2, "#b2182b")):
        pts = pv.PolyData(c[sel]); pts["v"] = f[sel]
        p.add_mesh(pts.glyph(orient="v", scale=False, factor=5.0, geom=pv.Line(pointa=(-0.5, 0, 0), pointb=(0.5, 0, 0))), color=col, line_width=2)
    p.add_text(f"Skin HGO fibre families (+/-{gam:g} deg about e1 = {meta.get('fiber_ref', '')})", font_size=10)
    p.view_vector((1.0, -1.0, 0.3)); p.camera.zoom(1.3)
    p.add_text(meta.get("run_name", "") + " | reference configuration", position="lower_left", font_size=11, viewport=False)
    p.screenshot(out_png)
    p.close()
    return out_png


def _zoom_box(cut):
    b = cut.bounds
    zc = 0.5 * (b[4] + b[5])
    return (b[0], b[1], b[2], b[2] + 0.35 * (b[3] - b[2]), zc - 25, zc + 25)


def mesh_quality(run_path, vtu0, meta, out_png):
    pv = _pv()
    blocks = _blocks(run_path)
    g = pv.read(vtu0)
    tis = _ids(blocks, lambda b: b["name"] in TISSUE)
    T = _subset(g, tis)
    cut = T.clip(normal="x", origin=(0, 0, 0), invert=True, crinkle=True)
    p = pv.Plotter(off_screen=True, shape=(1, 2), window_size=(1400, 820), border=False)
    for i, (name, title) in enumerate([("radius_ratio", "Tet radius ratio (1 = regular)"),
                                       ("min_scaled_jacobian", "Min scaled Jacobian (initial)")]):
        p.subplot(0, i)
        p.add_mesh(cut, scalars=name, cmap="RdYlGn", clim=[0, 1], preference="cell", show_edges=True, edge_color="#333333",
                   line_width=0.3, scalar_bar_args={"title": title, "vertical": False})
        p.add_text(title, font_size=10)
        p.view_vector((1, 0, 0)); p.camera.zoom(1.25)
    p.add_text(meta.get("run_name", "") + " | mesh inspection (reference configuration)", position="lower_left", font_size=11, viewport=False)
    p.screenshot(out_png)
    p.close()
    return out_png


def convergence_plot(conv: dict, out_png, title=""):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
    steps = conv.get("time_steps", [])
    for s in steps:
        r = [it["residual"][1] / max(it["residual"][0], 1e-300) for it in s["norms"] if "residual" in it and it["residual"][0] > 0]
        e = [it["energy"][1] / max(it["energy"][0], 1e-300) for it in s["norms"] if "energy" in it and it["energy"][0] > 0]
        col = "tab:blue" if s["converged"] else "tab:red"
        if r:
            ax[0].semilogy(range(1, len(r) + 1), r, color=col, alpha=0.5)
        if e:
            ax[1].semilogy(range(1, len(e) + 1), e, color=col, alpha=0.5)
    ax[0].set_title("Residual norm / initial (blue = converged step, red = failed attempt)")
    ax[1].set_title("Energy norm / initial")
    for a in ax:
        a.set_xlabel("Newton iteration"); a.grid(alpha=0.3)
    fig.suptitle(title, fontsize=9)
    fig.tight_layout(); fig.savefig(out_png, dpi=130); plt.close(fig)
    return out_png


def xy_plot(series, xlabel, ylabel, out_png, title=""):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    for s in series:
        ax.plot(s["x"], s["y"], s.get("style", "o-"), label=s["label"])
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel); ax.grid(alpha=0.3); ax.legend(fontsize=8)
    ax.set_title(title, fontsize=9)
    fig.tight_layout(); fig.savefig(out_png, dpi=130); plt.close(fig)
    return out_png


def free_body_diagram(load_report: dict, petals_meta: list, out_png, title=""):
    """Top and oblique views of the lace route with per-guide forces."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cab = load_report["cable"]
    pts = np.array([r["xyz"] for r in cab["route"]])
    fig = plt.figure(figsize=(13, 5.6))
    ax = fig.add_subplot(1, 2, 1)
    ax.plot(pts[:, 0], pts[:, 1], "-", color="0.5", lw=1, label="lace (free spans + channels)")
    for p in petals_meta:
        th = np.deg2rad(np.linspace(p["theta_center_deg"] - p["half_width_deg"], p["theta_center_deg"] + p["half_width_deg"], 30))
        r = p.get("r_ref", 60)
        ax.plot(r * np.cos(th), r * np.sin(th), lw=6, alpha=0.35)
        ax.text(1.15 * r * np.cos(np.deg2rad(p["theta_center_deg"])), 1.15 * r * np.sin(np.deg2rad(p["theta_center_deg"])), p["name"], ha="center", fontsize=8)
    smax = max(np.linalg.norm(pf["force_N"]) for b in cab["body_loads"].values() for pf in b["point_forces"])
    for b, v in cab["body_loads"].items():
        for pf in v["point_forces"]:
            x, f = np.array(pf["xyz"]), np.array(pf["force_N"])
            ax.arrow(x[0], x[1], 25 * f[0] / smax, 25 * f[1] / smax, width=0.6, color="tab:red", length_includes_head=True)
        com = np.array(load_report["cable"]["coms"][b]); F = np.array(v["F_N"])
        ax.arrow(com[0], com[1], 25 * F[0] / smax, 25 * F[1] / smax, width=1.2, color="tab:blue", length_includes_head=True)
        ax.text(com[0], com[1], f"{np.linalg.norm(F):.1f} N", fontsize=7, color="tab:blue")
    ax.set_aspect("equal"); ax.set_title("Top view: guide forces (red) and petal resultants at COM (blue)", fontsize=9)
    ax.set_xlabel("x [mm]"); ax.set_ylabel("y (anterior) [mm]")
    ax2 = fig.add_subplot(1, 2, 2)
    rows = [(b, *np.round(v["F_N"], 2), *np.round(v["M_com_Nmm"], 1)) for b, v in sorted(cab["body_loads"].items())]
    ax2.axis("off")
    tbl = ax2.table(cellText=[[r[0]] + [f"{x:g}" for x in r[1:]] for r in rows] +
                    [["SUM F", *[f"{x:.2e}" for x in cab["sum_force_N"]], "", "", ""],
                     ["SUM M@O", "", "", "", *[f"{x:.2e}" for x in cab["sum_moment_about_origin_Nmm"]]]],
                    colLabels=["body", "Fx N", "Fy N", "Fz N", "Mx Nmm", "My Nmm", "Mz Nmm"], loc="center")
    tbl.auto_set_font_size(False); tbl.set_fontsize(8); tbl.scale(1, 1.3)
    ax2.set_title(f"Lace tension {load_report['tension']['applied_N']} N ({load_report['tension']['basis']}); "
                  f"guide friction mu_g = {cab['guide_friction']}", fontsize=9)
    fig.suptitle(title, fontsize=9)
    fig.tight_layout(); fig.savefig(out_png, dpi=130); plt.close(fig)
    return out_png
