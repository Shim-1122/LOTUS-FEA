"""FEBio 4 adapter: input generation, capability checks, execution, parsing.

Target: FEBio 4.13 (febio_spec 4.0), built from the official source tag v4.13
with Intel MKL PARDISO.  The adapter only emits constructs whose parameter
names were verified against that version's source and user manual (see
docs/solver_selection.md).  Nothing here implements a constitutive model or
a nonlinear solver: FEBio does all mechanics.
"""
from __future__ import annotations

import hashlib
import os
import re
import shutil
import subprocess
import time
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field

import numpy as np

from ..mesh.model import FEModel

FEBIO_EXE_ENV = "RLFEA_FEBIO"


def find_febio() -> str | None:
    exe = os.environ.get(FEBIO_EXE_ENV) or shutil.which("febio4")
    return exe


_TINY = """<?xml version="1.0" encoding="ISO-8859-1"?>
<febio_spec version="4.0"><Module type="solid"/>
<Material><material id="1" name="m" type="neo-Hookean"><E>1</E><v>0.3</v></material></Material>
<Mesh><Nodes name="n"><node id="1">0,0,0</node><node id="2">1,0,0</node><node id="3">0,1,0</node><node id="4">0,0,1</node></Nodes>
<Elements type="tet4" name="p"><elem id="1">1,2,3,4</elem></Elements><NodeSet name="s">1,2,3,4</NodeSet></Mesh>
<MeshDomains><SolidDomain name="p" mat="m"/></MeshDomains>
<Boundary><bc type="zero displacement" node_set="s"><x_dof>1</x_dof><y_dof>1</y_dof><z_dof>1</z_dof></bc></Boundary>
<Step><step id="1"><Control><time_steps>1</time_steps><step_size>1</step_size><solver type="solid"/></Control></step></Step>
</febio_spec>
"""


def febio_version(exe: str | None = None) -> dict:
    """Run a one-element model and read the version banner from its log."""
    exe = exe or find_febio()
    if not exe:
        return {"available": False}
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        open(os.path.join(d, "v.feb"), "w").write(_TINY)
        p = subprocess.run([exe, "-nosplash", "-silent", "-noconfig", "-i", "v.feb"], capture_output=True, text=True,
                           cwd=d, timeout=120, stdin=subprocess.DEVNULL)
        log = open(os.path.join(d, "v.log"), errors="replace").read() if os.path.exists(os.path.join(d, "v.log")) else ""
    m = re.search(r"[Vv]ersion\s+([0-9]+\.[0-9]+\.[0-9]+(?:\.[0-9a-f]+)?)", log + p.stdout)
    return {"available": True, "executable": os.path.realpath(exe), "version": m.group(1) if m else None,
            "smoke_test_normal_termination": "N O R M A L   T E R M I N A T I O N" in log}


# ---------------------------------------------------------------------------
# Analysis specification
# ---------------------------------------------------------------------------
@dataclass
class StepSpec:
    name: str
    duration: float
    n_steps: int
    must_points: list[float]               # absolute times
    rigid_bcs: list[dict] = field(default_factory=list)
    rigid_loads: list[dict] = field(default_factory=list)
    loads: list[dict] = field(default_factory=list)
    bcs: list[dict] = field(default_factory=list)
    mp_lc: int | None = None               # must-point load curve id (set by the case builder)


@dataclass
class AnalysisSpec:
    materials: dict                    # name -> {febio_type|rigid, params, com}
    block_material: dict               # block name -> material name
    block_domain: dict                 # block name -> {type, params}
    bcs: list = field(default_factory=list)
    rigid_bcs: list = field(default_factory=list)
    rigid_loads: list = field(default_factory=list)
    connectors: list = field(default_factory=list)
    contacts: list = field(default_factory=list)
    load_curves: dict = field(default_factory=dict)
    steps: list = field(default_factory=list)
    solver: dict = field(default_factory=dict)
    time_stepper: dict = field(default_factory=dict)
    log_surfaces: list = field(default_factory=list)
    log_rigid: list = field(default_factory=list)
    log_connectors: list = field(default_factory=list)
    plot_vars: list = field(default_factory=lambda: ["displacement", "stress", "Lagrange strain", "relative volume",
                                                      "contact pressure", "contact gap", "contact traction",
                                                      "contact stick", "reaction forces"])


def _fmt(v):
    if isinstance(v, (list, tuple, np.ndarray)):
        return ",".join(_fmt(x) for x in v)
    if isinstance(v, bool):
        return "1" if v else "0"
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    return repr(float(v))


def _sub(parent, tag, text=None, **attrib):
    e = ET.SubElement(parent, tag, {k: str(v) for k, v in attrib.items()})
    if text is not None:
        e.text = _fmt(text) if not isinstance(text, str) else text
    return e


def _ids(arr):
    return ",".join(str(int(x) + 1) for x in arr)


def _renumber_load_curves(spec: AnalysisSpec) -> AnalysisSpec:
    """FEBio requires load-controller ids 1..n; remap all references."""
    import copy
    sp = copy.deepcopy(spec)
    old = sorted(sp.load_curves)
    mp = {o: i + 1 for i, o in enumerate(old)}
    sp.load_curves = {mp[o]: sp.load_curves[o] for o in old}

    def fix(lst):
        for d in lst:
            if "lc" in d and d["lc"] is not None:
                d["lc"] = mp[d["lc"]]
    fix(sp.bcs); fix(sp.rigid_bcs); fix(sp.rigid_loads)
    for st in sp.steps:
        fix(st.bcs); fix(st.rigid_bcs); fix(st.rigid_loads); fix(st.loads)
        if st.mp_lc is not None:
            st.mp_lc = mp[st.mp_lc]
    return sp


def write_feb(path: str, model: FEModel, spec: AnalysisSpec) -> dict:
    """Write a febio_spec 4.0 input file.  Returns id maps for parsing."""
    spec = _renumber_load_curves(spec)
    root = ET.Element("febio_spec", version="4.0")
    _sub(root, "Module", type="solid")
    # NOTE: an empty <Globals><Constants/></Globals> block makes FEBio 4.13
    # discard the Material section (verified), so Globals is omitted.
    # ---- materials -----------------------------------------------------
    mats = _sub(root, "Material")
    mat_id = {}
    for i, (name, m) in enumerate(spec.materials.items(), start=1):
        mat_id[name] = i
        if m.get("rigid"):
            me = _sub(mats, "material", id=i, name=name, type="rigid body")
            _sub(me, "density", 1.0)
            if m.get("com") is not None:
                _sub(me, "center_of_mass", list(m["com"]))
            continue
        me = _sub(mats, "material", id=i, name=name, type=m["febio_type"])
        for p, v in m["params"].items():
            _sub(me, p, v)
        if m["febio_type"] == "trans iso Mooney-Rivlin":
            _sub(me, "fiber", "1,0,0", type="vector")
    # ---- mesh -----------------------------------------------------------
    mesh = _sub(root, "Mesh")
    nodes = _sub(mesh, "Nodes", name="all_nodes")
    X = model.nodes
    lines = [f'<node id="{i + 1}">{x[0]:.10g},{x[1]:.10g},{x[2]:.10g}</node>' for i, x in enumerate(X)]
    nodes.text = "__NODES__"
    elem_offset = {}
    off = 0
    for b in model.blocks:
        el = _sub(mesh, "Elements", type=b.etype, name=b.name)
        elem_offset[b.name] = off
        el.text = f"__ELEMS_{b.name}__"
        off += len(b.conn)
    for name, s in spec_surfaces(model, spec).items():
        se = _sub(mesh, "Surface", name=name)
        se.text = f"__SURF_{name}__"
    used_ns = {bc.get("node_set") for bc in spec.bcs} | {bc.get("node_set") for st in spec.steps for bc in st.bcs}
    for name, ns in model.nodesets.items():
        if name in used_ns:
            _sub(mesh, "NodeSet", _ids(ns), name=name)
    for c in spec.contacts:
        sp = _sub(mesh, "SurfacePair", name=c["name"])
        _sub(sp, "primary", c["primary"])
        _sub(sp, "secondary", c["secondary"])
    # ---- domains --------------------------------------------------------
    dom = _sub(root, "MeshDomains")
    for b in model.blocks:
        d = spec.block_domain.get(b.name, {})
        attrs = {"name": b.name, "mat": spec.block_material[b.name]}
        if d.get("type"):
            attrs["type"] = d["type"]
        de = _sub(dom, "SolidDomain", **attrs)
        for k, v in d.get("params", {}).items():
            _sub(de, k, v)
    # ---- mesh data: per-element material axes ---------------------------
    md = None
    for b in model.blocks:
        mname = spec.block_material[b.name]
        m = spec.materials[mname]
        if b.mat_axis is not None and not m.get("rigid") and m.get("febio_type") == "Holzapfel-Gasser-Ogden":
            md = md if md is not None else _sub(root, "MeshData")
            ed = _sub(md, "ElementData", type="mat_axis", elem_set=b.name)
            ed.text = f"__MATAXIS_{b.name}__"
    # ---- boundary -------------------------------------------------------
    if spec.bcs:
        bnd = _sub(root, "Boundary")
        for bc in spec.bcs:
            _bc(bnd, bc)
    # ---- rigid (global) ---------------------------------------------------
    if spec.rigid_bcs or spec.connectors or spec.rigid_loads:
        rg = _sub(root, "Rigid")
        for r in spec.rigid_bcs:
            _rigid_bc(rg, r)
        for r in spec.rigid_loads:
            _rigid_load(rg, r)
        for c in spec.connectors:
            e = _sub(rg, "rigid_connector", type=c["type"], name=c["name"])
            for k, v in c["params"].items():
                _sub(e, k, v)
    # ---- contact -------------------------------------------------------------
    if spec.contacts:
        ct = _sub(root, "Contact")
        for c in spec.contacts:
            e = _sub(ct, "contact", type=c["type"], name=c["name"], surface_pair=c["name"])
            for k, v in c["params"].items():
                _sub(e, k, v)
    # ---- load data -------------------------------------------------------
    ld = _sub(root, "LoadData")
    for lc_id, lc in spec.load_curves.items():
        e = _sub(ld, "load_controller", id=lc_id, name=lc.get("name", f"lc{lc_id}"), type="loadcurve")
        _sub(e, "interpolate", lc.get("interpolate", "LINEAR"))
        _sub(e, "extend", lc.get("extend", "CONSTANT"))
        pts = _sub(e, "points")
        for t, v in lc["points"]:
            _sub(pts, "pt", f"{t:.10g},{v:.10g}")
    # ---- output -----------------------------------------------------------
    out = _sub(root, "Output")
    pf = _sub(out, "plotfile", type="febio")
    for v in spec.plot_vars:
        _sub(pf, "var", type=v)
    lf = _sub(out, "logfile")
    _sub(lf, "node_data", data="ux;uy;uz", file="out_node_disp.txt", delim=",")
    _sub(lf, "node_data", data="Rx;Ry;Rz", file="out_node_react.txt", delim=",")
    _sub(lf, "element_data", data="Fxx;Fxy;Fxz;Fyx;Fyy;Fyz;Fzx;Fzy;Fzz", file="out_elem_F.txt", delim=",")
    _sub(lf, "element_data", data="sx;sy;sz;sxy;syz;sxz", file="out_elem_stress.txt", delim=",")
    _sub(lf, "element_data", data="sed", file="out_elem_sed.txt", delim=",")
    for s in spec.log_surfaces:
        _sub(lf, "face_data", data="contact gap;contact pressure;contact_traction.x;contact_traction.y;contact_traction.z",
             surface=s, file=f"out_face_{s}.txt", delim=",")
    if spec.log_rigid:
        _sub(lf, "rigid_body_data", data="x;y;z;thx;thy;thz;Fx;Fy;Fz;Mx;My;Mz", file="out_rigid.txt", delim=",")
    if spec.log_connectors:
        _sub(lf, "rigid_connector_data", data="RCFx;RCFy;RCFz;RCMx;RCMy;RCMz", file="out_connectors.txt", delim=",")
    # ---- steps -------------------------------------------------------------
    stp = _sub(root, "Step")
    for i, s in enumerate(spec.steps, start=1):
        se = _sub(stp, "step", id=i, name=s.name)
        ctl = _sub(se, "Control")
        _sub(ctl, "analysis", "STATIC")
        _sub(ctl, "time_steps", s.n_steps)
        _sub(ctl, "step_size", s.duration / s.n_steps)
        _sub(ctl, "plot_level", "PLOT_MUST_POINTS")
        _sub(ctl, "output_level", "OUTPUT_MUST_POINTS")
        _sub(ctl, "plot_zero_state", 1)
        sol = _sub(ctl, "solver", type="solid")
        for k, v in spec.solver.items():
            if k == "qn_method":
                q = _sub(sol, "qn_method", type=v["type"])
                for kk, vv in v.items():
                    if kk != "type":
                        _sub(q, kk, vv)
            elif k == "linear_solver":
                _sub(sol, "linear_solver", type=v)
            else:
                _sub(sol, k, v)
        ts = _sub(ctl, "time_stepper", type="default")
        for k, v in spec.time_stepper.items():
            if k == "dtmax":
                if s.mp_lc is not None:
                    _sub(ts, "dtmax", v, lc=s.mp_lc)
                else:
                    _sub(ts, "dtmax", v)
            else:
                _sub(ts, k, v)
        if s.bcs:
            bnd = _sub(se, "Boundary")
            for bc in s.bcs:
                _bc(bnd, bc)
        if s.rigid_bcs or s.rigid_loads:
            rg = _sub(se, "Rigid")
            for r in s.rigid_bcs:
                _rigid_bc(rg, r)
            for r in s.rigid_loads:
                _rigid_load(rg, r)
        if s.loads:
            lo = _sub(se, "Loads")
            for L in s.loads:
                e = _sub(lo, "surface_load", type="pressure", surface=L["surface"], name=L.get("name", ""))
                _sub(e, "pressure", L["value"], lc=L["lc"])
                _sub(e, "linear", 0)
                _sub(e, "symmetric_stiffness", 1)

    ET.indent(root, space=" ")
    txt = '<?xml version="1.0" encoding="ISO-8859-1"?>\n' + ET.tostring(root, encoding="unicode")
    # bulk sections are inserted as text for speed
    txt = txt.replace("__NODES__", "\n" + "\n".join(lines) + "\n")
    eid = 0
    for b in model.blocks:
        ls = []
        for c in b.conn:
            eid += 1
            ls.append(f'<elem id="{eid}">{",".join(str(int(n) + 1) for n in c)}</elem>')
        txt = txt.replace(f"__ELEMS_{b.name}__", "\n" + "\n".join(ls) + "\n")
    for name, s in spec_surfaces(model, spec).items():
        ls = [f'<{s.etype} id="{i + 1}">{",".join(str(int(n) + 1) for n in f)}</{s.etype}>' for i, f in enumerate(s.conn)]
        txt = txt.replace(f"__SURF_{name}__", "\n" + "\n".join(ls) + "\n")
    for b in model.blocks:
        tag = f"__MATAXIS_{b.name}__"
        if tag in txt:
            ls = [f'<elem lid="{i + 1}"><a>{a[0]:.8g},{a[1]:.8g},{a[2]:.8g}</a><d>{d[0]:.8g},{d[1]:.8g},{d[2]:.8g}</d></elem>'
                  for i, (a, d) in enumerate(b.mat_axis)]
            txt = txt.replace(tag, "\n" + "\n".join(ls) + "\n")
    with open(path, "w") as fh:
        fh.write(txt)
    return {"elem_offset": elem_offset, "mat_id": mat_id, "sha256": hashlib.sha256(txt.encode()).hexdigest()}


def spec_surfaces(model: FEModel, spec: AnalysisSpec):
    names = set(spec.log_surfaces)
    for c in spec.contacts:
        names.add(c["primary"]); names.add(c["secondary"])
    for s in spec.steps:
        for L in s.loads:
            names.add(L["surface"])
    return {n: model.surfaces[n] for n in sorted(names)}


def _rigid_load(parent, r):
    e = _sub(parent, "rigid_load", type=r["type"], name=r.get("name", ""))
    _sub(e, "rb", r["rb"])
    _sub(e, "dof", r["dof"])
    _sub(e, "value", r["value"], lc=r["lc"])
    if r["type"] == "rigid_force":
        _sub(e, "load_type", r.get("load_type", 0))
    _sub(e, "relative", 0)


def _bc(parent, bc):
    t = bc.get("type", "zero displacement")
    e = _sub(parent, "bc", type=t, node_set=bc["node_set"], name=bc.get("name", bc["node_set"] + "_" + t.split()[0]))
    if t == "zero displacement":
        for dof in ("x", "y", "z"):
            _sub(e, f"{dof}_dof", 1 if dof in bc["dofs"] else 0)
    elif t == "prescribed displacement":
        _sub(e, "dof", bc["dof"])
        _sub(e, "value", bc["value"], lc=bc["lc"])
        _sub(e, "relative", int(bc.get("relative", 0)))
    else:
        raise ValueError(t)


def _rigid_bc(parent, r):
    if r["type"] == "rigid_fixed":
        e = _sub(parent, "rigid_bc", type="rigid_fixed", name=r.get("name", ""))
        _sub(e, "rb", r["rb"])
        for d in ("Rx", "Ry", "Rz", "Ru", "Rv", "Rw"):
            if d in r["dofs"]:
                _sub(e, f"{d}_dof", 1)
    elif r["type"] in ("rigid_displacement", "rigid_rotation"):
        e = _sub(parent, "rigid_bc", type=r["type"], name=r.get("name", ""))
        _sub(e, "rb", r["rb"])
        _sub(e, "dof", r["dof"])
        _sub(e, "value", r["value"], lc=r["lc"])
        _sub(e, "relative", int(r.get("relative", 0)))
    else:
        raise ValueError(r["type"])


# ---------------------------------------------------------------------------
# Execution
# ---------------------------------------------------------------------------
def run_febio(feb_path: str, workdir: str, exe: str | None = None, dump: bool = True,
              restart: str | None = None, timeout: float | None = None, threads: int | None = None) -> dict:
    exe = exe or find_febio()
    if not exe:
        raise RuntimeError("FEBio executable not found (set RLFEA_FEBIO or put febio4 on PATH)")
    base = os.path.splitext(os.path.basename(feb_path))[0]
    cmd = [exe, "-nosplash", "-silent", "-noconfig"]
    if restart:
        cmd = [exe, "-nosplash", "-silent", "-noconfig", "-r", restart]
    else:
        cmd += ["-i", os.path.basename(feb_path)]
        if dump:
            cmd += ["-dump", f"{base}.dmp"]
    env = dict(os.environ)
    if threads:
        env["OMP_NUM_THREADS"] = str(threads)
        env["MKL_NUM_THREADS"] = str(threads)
    t0 = time.time()
    with open(os.path.join(workdir, "febio_stdout.txt"), "w") as so:
        try:
            p = subprocess.run(cmd, cwd=workdir, stdout=so, stderr=subprocess.STDOUT, timeout=timeout, env=env,
                               stdin=subprocess.DEVNULL)
            rc = p.returncode
            interrupted = False
        except subprocess.TimeoutExpired:
            rc, interrupted = None, True
    wall = time.time() - t0
    log = os.path.join(workdir, f"{base}.log")
    conv = parse_log(log) if os.path.exists(log) else {"terminated": "no log", "last_converged_time": 0.0,
                                                         "input_error": open(os.path.join(workdir, "febio_stdout.txt")).read()[-2000:]}
    status = classify(rc, conv, interrupted)
    return {"command": cmd, "returncode": rc, "wall_time_s": wall, "status": status, "convergence": conv,
            "log": log}


def classify(rc, conv, interrupted=False) -> str:
    if interrupted:
        return "interrupted"
    if conv.get("normal_termination") and rc == 0:
        return "converged"
    return "failed"


_TS_RE = re.compile(r"^\s*([0-9]+)\s+([0-9]+)\s+([0-9]+)\s+([0-9]+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)")


def parse_log(path: str) -> dict:
    """Extract per-time-step convergence information from an FEBio log."""
    txt = open(path, errors="replace").read()
    out = {"normal_termination": "N O R M A L   T E R M I N A T I O N" in txt,
           "error_termination": "E R R O R   T E R M I N A T I O N" in txt}
    m = re.search(r"[Vv]ersion\s+([0-9.]+[0-9a-f.]*)", txt)
    out["febio_version"] = m.group(1) if m else None
    steps = []
    # iteration blocks: "convergence norms :" followed by residual / energy / displacement lines
    for blk in re.finditer(r"===== beginning time step\s+(\d+)\s*:\s*([0-9.eE+-]+)\s*=====(.*?)(?=(?:===== beginning time step)|\Z)", txt, re.S):
        nstep, tval, body = int(blk.group(1)), float(blk.group(2)), blk.group(3)
        iters = []
        for it in re.finditer(r"convergence norms\s*:\s*INITIAL\s+CURRENT\s+REQUIRED(.*?)(?=\n\s*\n)", body, re.S):
            norms = {}
            for line in it.group(1).strip().splitlines():
                mm = re.match(r"\s*(residual|energy|displacement)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)", line)
                if mm:
                    norms[mm.group(1)] = [float(mm.group(2)), float(mm.group(3)), float(mm.group(4))]
            if norms:
                iters.append(norms)
        conv = "------- converged at time" in body
        fail = re.search(r"(failed to converge|negative jacobian|Negative jacobian|ERROR)", body)
        aug = len(re.findall(r"augmentation", body, re.I))
        steps.append({"step": nstep, "time": tval, "converged": conv, "iterations": len(iters), "norms": iters,
                      "failure": fail.group(1) if (fail and not conv) else None, "augmentations": aug})
    out["time_steps"] = steps
    acc = [s for s in steps if s["converged"]]
    out["last_converged_time"] = acc[-1]["time"] if acc else 0.0
    m2 = re.findall(r"converged at time\s*:\s*([0-9.eE+-]+)", txt)
    if m2:
        out["last_converged_time"] = float(m2[-1])
    out["n_retries"] = len(re.findall(r"Retrying time step", txt))
    out["negative_jacobian_events"] = len(re.findall(r"[Nn]egative jacobian", txt))
    m3 = re.search(r"Total elapsed time\s*\.+\s*:\s*([0-9:]+)", txt)
    out["elapsed"] = m3.group(1) if m3 else None
    return out


def read_logdata(path: str) -> list[dict]:
    """Parse a logfile data file written with file="..." (one record per
    output time).  Returns [{'step', 'time', 'ids', 'data'}]."""
    recs = []
    if not os.path.exists(path):
        return recs
    cur = None
    with open(path) as fh:
        for line in fh:
            if line.startswith("*Step"):
                if cur:
                    recs.append(cur)
                cur = {"step": int(line.split("=")[1]), "rows": []}
            elif line.startswith("*Time"):
                cur["time"] = float(line.split("=")[1])
            elif line.startswith("*"):
                continue
            elif line.strip() and cur is not None:
                cur["rows"].append(line.strip())
    if cur:
        recs.append(cur)
    out = []
    for r in recs:
        arr = np.array([[float(x) for x in row.split(",")] for row in r["rows"]]) if r["rows"] else np.zeros((0, 1))
        out.append({"step": r["step"], "time": r.get("time"), "ids": arr[:, 0].astype(int) if len(arr) else np.zeros(0, int),
                    "data": arr[:, 1:] if len(arr) else arr})
    return out
