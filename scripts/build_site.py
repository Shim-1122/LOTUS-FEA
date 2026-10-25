"""Build the static website in site/ from one run directory.

usage: python scripts/build_site.py --run runs/knee_limb_segmented_T75 \
           --verification docs/verification/verification.json \
           --repo-url https://github.com/Shim-1122/knee-fea [--out site]

The page states the run's real solver state; a failed or interrupted step is
shown as such, and nothing is filled in from other runs.
"""
from __future__ import annotations

import argparse
import datetime as dt
import html
import json
import os
import shutil
import string
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))


def _state_at(states, t):
    c = [s for s in states if abs(s["time"] - t) < 1e-6]
    return c[0] if c else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--verification", required=True)
    ap.add_argument("--repo-url", required=True)
    ap.add_argument("--out", default=os.path.join(HERE, "..", "site"))
    a = ap.parse_args()
    run, out = a.run, a.out
    os.makedirs(os.path.join(out, "figures"), exist_ok=True)

    man = json.load(open(os.path.join(run, "manifest.json")))
    st = json.load(open(os.path.join(run, "status.json")))
    summ = json.load(open(os.path.join(run, "results", "summary.json")))
    lay = json.load(open(os.path.join(run, "mesh", "layer_report.json")))
    load = json.load(open(os.path.join(run, "solver", "load_report.json")))
    ver = json.load(open(a.verification))
    T = load["tension"]["applied_N"]
    states = summ["states"]
    s1 = _state_at(states, 1.0)
    if s1 is None:
        raise SystemExit("run never reached the end of tightening (t = 1); the site needs that state")

    def stat(v, k, n=""):
        return f'<div class="stat"><span class="v">{v}</span><span class="k">{html.escape(k)}</span><span class="n">{html.escape(n)}</span></div>'

    c1 = s1["contact"]
    stats = "".join([
        stat(f"{s1['max_tissue_displacement_mm']:.1f} mm", "Max tissue displacement", "anywhere in skin, fat or muscle"),
        stat(f"{1000 * c1['max_pressure_MPa']:.0f} kPa", "Peak contact pressure", "mesh-dependent, underestimated"),
        stat(f"{c1['normal_force_N']:.0f} N", "Total contact force", "skin on pads and cushion"),
        stat(f"{s1['deformation']['min_J']:.2f}", "Lowest volume ratio J", "1 = unchanged"),
        stat(f"{c1.get('max_penetration_mm', 0):.2f} mm", "Max penetration", "penalty contact regularisation"),
    ])
    rows = []
    for s in states:
        t = s["time"]
        if t > 1 + 1e-9 or abs(round(t * 4) - t * 4) < 1e-6:   # quarter steps of tightening + all pull-out outputs
            load_txt = f"lace {T * t:.1f} N" if t <= 1 + 1e-9 else f"lace {T:.0f} N + pull-out {2.0 * (t - 1):.2f} mm"
            c = s.get("contact", {})
            rows.append(f"<tr><td class='num'>t = {t:.2f}</td><td>{load_txt}</td>"
                        f"<td class='num'>{s['max_tissue_displacement_mm']:.2f}</td>"
                        f"<td class='num'>{1000 * c.get('max_pressure_MPa', 0):.0f}</td>"
                        f"<td class='num'>{c.get('normal_force_N', 0):.0f}</td>"
                        f"<td class='num'>{s['deformation']['min_J']:.3f}</td>"
                        f"<td class='num'>{s['equilibrium']['force_imbalance_relative']:.1e}</td></tr>")
    last = st.get("last_converged_time", 0) or 0
    state = st["state"]
    if state == "converged":
        label, cls = "solver converged (tightening + pull-out)", "ok"
        pull = "The pull-out step converged. Even so, a converged pull-out is not proof that the socket grips a real limb."
    else:
        label = f"tightening converged; run {state} at t = {last:.3f}" if last >= 1 - 1e-9 else f"run {state} at t = {last:.3f}"
        cls = "bad"
        pull = (f"The pull-out step did not complete: the run is {state} at t = {last:.3f} "
                f"(t = 1 is the end of tightening). {st.get('note', '')}").strip()
    fb = [v["residual_over_cable"] for v in (s1.get("petal_free_body") or {}).values()]
    vals = dict(
        febio_version=html.escape(str(man["solver"].get("version"))),
        state_label=html.escape(label), state_class=cls,
        repo_url=html.escape(a.repo_url.rstrip("/")),
        stats=stats, rows="".join(rows),
        verification_pass=sum(c["pass"] for c in ver["checks"]), verification_total=len(ver["checks"]),
        dev_skin=f"{lay.get('offset_vs_segmentation_mm', {}).get('skin_inner', {}).get('median', float('nan')):.2f}",
        dev_fat=f"{lay.get('offset_vs_segmentation_mm', {}).get('fat_inner', {}).get('median', float('nan')):.2f}",
        fb_resid=f"{100 * max(fb):.2f}" if fb else "n/a",
        pullout_note=html.escape(pull),
        run_id=html.escape(man["run_id"]),
        commit=html.escape(str(man.get("repository", {}).get("commit", "?"))[:10]),
        generated=dt.date.today().isoformat())
    tpl = string.Template(open(os.path.join(HERE, "site_template.html")).read())
    open(os.path.join(out, "index.html"), "w").write(tpl.substitute(vals))

    for f in ("full_limb_1x.png", "contact.png", "sections.png", "force_displacement.png", "free_body.png"):
        src = os.path.join(run, "figures", f)
        if os.path.exists(src):
            shutil.copy(src, os.path.join(out, "figures", f))
    shutil.copy(os.path.join(run, "report", "report.pdf"), os.path.join(out, "report.pdf"))
    from rlfea.ui.web_export import export
    v = os.path.join(out, "viewer.html")
    export(run, v)
    body = open(v).read()
    if not body.lstrip().lower().startswith("<!doctype"):
        open(v, "w").write("<!doctype html>\n" + body.replace("<title>", '<meta name="viewport" content="width=device-width, initial-scale=1">\n<title>', 1))
    print("site written to", os.path.abspath(out))


if __name__ == "__main__":
    main()
