"""Run directories, manifests and job status.

A run directory contains everything needed to audit or reproduce the run:

    manifest.json        inputs + hashes, repository commit, environment, solver
    status.json          prepared | running | converged | failed | interrupted
    config.resolved.json resolved configuration (internal units)
    geometry/            audits, transform/repair history
    mesh/                model mesh (VTU) with region ids and quality
    solver/              FEBio input (.feb), log, stdout, raw logfile data, xplt
    results/             VTU/VTP per output time, PVD collection, JSON checks
    figures/             PNG figures (labelled with case/solver/state/scale)
    report/              LaTeX source (+ PDF when compiled)
"""
from __future__ import annotations

import datetime as dt
import hashlib
import json
import os
import platform
import subprocess
import sys

STATES = ("prepared", "running", "converged", "failed", "interrupted")


def _json_default(o):
    import numpy as np
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if hasattr(o, "__dict__"):
        return {k: v for k, v in o.__dict__.items() if not k.startswith("_")}
    return str(o)


def dump(path, obj):
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=1, default=_json_default)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def repo_commit() -> dict:
    root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    try:
        c = subprocess.run(["git", "-C", root, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
        dirty = bool(subprocess.run(["git", "-C", root, "status", "--porcelain"], capture_output=True, text=True).stdout.strip())
        return {"commit": c or None, "dirty_worktree": dirty}
    except Exception:  # noqa: BLE001
        return {"commit": None, "dirty_worktree": None}


def environment() -> dict:
    pk = {}
    for name in ("numpy", "scipy", "gmsh", "meshio", "trimesh", "pyvista", "vtk", "matplotlib", "yaml"):
        try:
            mod = __import__(name)
            pk[name] = getattr(mod, "__version__", getattr(mod, "GMSH_API_VERSION", "?"))
        except Exception:  # noqa: BLE001
            pk[name] = None
    return {"python": sys.version.split()[0], "platform": platform.platform(), "packages": pk}


class RunDir:
    def __init__(self, root: str, run_id: str | None = None):
        stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        self.run_id = run_id or stamp
        self.path = os.path.join(root, self.run_id)
        for sub in ("geometry", "mesh", "solver", "results", "figures", "report"):
            os.makedirs(os.path.join(self.path, sub), exist_ok=True)

    def p(self, *parts):
        return os.path.join(self.path, *parts)

    def set_status(self, state: str, **info):
        assert state in STATES, state
        rec = {"state": state, "updated_utc": dt.datetime.now(dt.timezone.utc).isoformat(), **info}
        if state == "running":
            rec["pid"] = os.getpid()
        dump(self.p("status.json"), rec)


def read_status(run_path: str) -> dict:
    p = os.path.join(run_path, "status.json")
    if not os.path.exists(p):
        return {"state": "unknown"}
    st = json.load(open(p))
    if st.get("state") == "running" and st.get("pid"):
        try:
            os.kill(int(st["pid"]), 0)
        except OSError:
            st["state"] = "interrupted"
            st["note"] = "process no longer alive; partial outputs only"
    return st
