"""Aggregate sensitivity runs (same load case, one parameter changed each)."""
from __future__ import annotations

import json
import os

import numpy as np

METRICS = [("max_tissue_displacement_mm", "max tissue displacement [mm]"),
           ("contact.max_pressure_MPa", "max contact pressure [MPa]"),
           ("contact.normal_force_N", "total contact normal force [N]"),
           ("contact.slip_fraction_of_contact_area", "slip area fraction [-]"),
           ("contact.max_penetration_mm", "max penetration [mm]"),
           ("deformation.min_J", "min J [-]"),
           ("region_extremes.muscle.max_logstrain", "muscle max log strain [-]"),
           ("region_extremes.skin.max_logstrain", "skin max log strain [-]")]


def _get(d, path):
    for k in path.split("."):
        if d is None or k not in d:
            return None
        d = d[k]
    return d


def collect(runs: dict, t: float = 1.0) -> dict:
    """runs: {label: run_path}.  Uses the state at time t (or the last
    converged state if t was not reached, flagged)."""
    out = {}
    for label, rp in runs.items():
        p = os.path.join(rp, "results", "summary.json")
        st = json.load(open(os.path.join(rp, "status.json"))) if os.path.exists(os.path.join(rp, "status.json")) else {}
        if not os.path.exists(p):
            out[label] = {"state": st.get("state", "missing"), "reached_t": None}
            continue
        s = json.load(open(p))
        cand = [x for x in s["states"] if abs(x["time"] - t) < 1e-6]
        rec = cand[0] if cand else (s["states"][-1] if s["states"] else {})
        row = {"state": st.get("state"), "reached_t": rec.get("time")}
        for key, _ in METRICS:
            row[key] = _get(rec, key)
        man = json.load(open(os.path.join(rp, "manifest.json")))
        row["n_elements"] = json.load(open(os.path.join(rp, "mesh", "model_summary.json")))["n_elements"]
        row["wall_time_s"] = man.get("execution", {}).get("wall_time_s")
        out[label] = row
    return out


def table_md(rows: dict, baseline: str, t: float | None = None) -> str:
    """Relative changes are only given for runs that reached time t."""
    hdr = "| case | state | t reached | elements | " + " | ".join(m[1] for m in METRICS) + " |"
    lines = [hdr, "|" + "---|" * (4 + len(METRICS))]
    b = rows.get(baseline, {})
    for k, r in rows.items():
        cells = []
        for key, _ in METRICS:
            v = r.get(key)
            if v is None:
                cells.append("-")
                continue
            bv = b.get(key)
            reached = t is None or (r.get("reached_t") is not None and abs(r["reached_t"] - t) < 1e-6)
            rel = f" ({100 * (v - bv) / abs(bv):+.0f}%)" if (bv not in (None, 0) and k != baseline and reached) else ""
            cells.append(f"{v:.4g}{rel}")
        flag = "" if (t is None or (r.get("reached_t") is not None and abs(r["reached_t"] - t) < 1e-6)) else " (target t NOT reached)"
        lines.append(f"| {k} | {r.get('state')}{flag} | {r.get('reached_t')} | {r.get('n_elements', '-')} | " + " | ".join(cells) + " |")
    return "\n".join(lines)


def plot(rows: dict, baseline: str, out_png: str, title="", t: float = 1.0):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    keys = [m for m in METRICS[:4]]
    labels = [k for k in rows if k != baseline and rows[k].get("reached_t") is not None and abs(rows[k]["reached_t"] - t) < 1e-6]
    fig, ax = plt.subplots(1, len(keys), figsize=(4.2 * len(keys), 4.2))
    for a, (key, name) in zip(ax, keys):
        bv = rows[baseline].get(key)
        vals = [100 * (rows[l][key] - bv) / abs(bv) if (rows[l].get(key) is not None and bv) else np.nan for l in labels]
        a.barh(labels, vals, color=["tab:red" if v > 0 else "tab:blue" for v in vals])
        a.axvline(0, color="k", lw=0.8)
        a.set_title(name, fontsize=9); a.set_xlabel("% change vs baseline")
        a.tick_params(labelsize=8)
    fig.suptitle(title + f" | baseline: {baseline}", fontsize=9)
    fig.tight_layout(); fig.savefig(out_png, dpi=130); plt.close(fig)
    return out_png
