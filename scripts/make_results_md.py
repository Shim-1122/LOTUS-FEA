"""Write docs/RESULTS.md from run directories (synthetic runs only; supplied
patient-derived runs are summarised privately, never in the repository).

usage: python scripts/make_results_md.py --verification DIR --main RUN --sens label=RUN ... --out docs/RESULTS.md
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from rlfea.post import sensitivity as SE  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--verification", required=True)
    ap.add_argument("--main", required=True)
    ap.add_argument("--sens", nargs="*", default=[])
    ap.add_argument("--baseline", default=None)
    ap.add_argument("--out", default="docs/RESULTS.md")
    ap.add_argument("--figdir", default="docs/figures")
    a = ap.parse_args()
    os.makedirs(a.figdir, exist_ok=True)
    L = ["# Results status (generated)", "",
         "All numbers below are NUMERICAL results for the SYNTHETIC limb and a PROVISIONAL socket, with the 75 N EXAMPLE strand tension.",
         "They are not clinical predictions, and they say nothing about grip or tissue safety. Supplied (patient-derived) runs stay private.", ""]
    v = json.load(open(os.path.join(a.verification, "verification.json")))
    npass = sum(c["pass"] for c in v["checks"])
    L += [f"## Verification ({npass}/{len(v['checks'])} checks passed)", "", open(os.path.join(a.verification, "verification.md")).read(), ""]
    man = json.load(open(os.path.join(a.main, "manifest.json")))
    st = json.load(open(os.path.join(a.main, "status.json")))
    s = json.load(open(os.path.join(a.main, "results", "summary.json")))
    L += [f"## Synthetic full-limb run `{man['run_id']}`", "",
          f"* state: **{st['state']}**, last converged time {st.get('last_converged_time')} (t = 1: tightening complete; t = 2: pull-out complete)",
          f"* solver: FEBio {man['solver'].get('version')}; wall time {man.get('execution', {}).get('wall_time_s', 0):.0f} s",
          f"* mesh: {json.load(open(os.path.join(a.main, 'mesh', 'model_summary.json')))['n_elements']} elements", "",
          "| t | load | max tissue \\|u\\| mm | min J | contact normal force N | max contact pressure kPa | slip area fraction | max penetration mm | cup axial constraint force N | global force imbalance |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    T = json.load(open(os.path.join(a.main, "solver", "load_report.json")))["tension"]["applied_N"]
    for r in s["states"]:
        c = r.get("contact", {})
        cup = r["equilibrium"]["rigid_bodies"].get("rigid_cup", {}).get("F", [0, 0, 0])
        load = f"{T * r['time']:.3g} N lace" if r["time"] <= 1 + 1e-9 else f"75 N + pull-out {2.0 * (r['time'] - 1):.2f} mm"
        L.append(f"| {r['time']:.3g} | {load} | {r['max_tissue_displacement_mm']:.3f} | {r['deformation']['min_J']:.3f} | "
                 f"{c.get('normal_force_N', 0):.1f} | {1000 * c.get('max_pressure_MPa', 0):.1f} | {c.get('slip_fraction_of_contact_area', 0):.3f} | "
                 f"{c.get('max_penetration_mm', 0):.3f} | {cup[2]:.2f} | {r['equilibrium']['force_imbalance_relative']:.1e} |")
    L.append("")
    for f in ("full_limb_1x.png", "sections.png", "contact.png", "pullout.png", "force_displacement.png", "free_body.png"):
        src = os.path.join(a.main, "figures", f)
        if os.path.exists(src):
            shutil.copy(src, os.path.join(a.figdir, f))
            L.append(f"![{f}](figures/{f})")
    L.append("")
    if a.sens:
        runs = dict(x.split("=", 1) for x in a.sens)
        base = a.baseline or list(runs)[0]
        rows = SE.collect(runs, 1.0)
        L += ["## Sensitivity (coarse synthetic mesh, end of tightening, t = 1)", "", SE.table_md(rows, base), ""]
        SE.plot(rows, base, os.path.join(a.figdir, "sensitivity.png"), "coarse synthetic model, 75 N example")
        L.append("![sensitivity](figures/sensitivity.png)")
        json.dump(rows, open(os.path.join(os.path.dirname(a.out), "sensitivity.json"), "w"), indent=1)
    open(a.out, "w").write("\n".join(L) + "\n")
    print("wrote", a.out)


if __name__ == "__main__":
    main()
