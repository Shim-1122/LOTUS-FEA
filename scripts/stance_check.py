"""Stance-grip verdict for a run with the body_weight load case.

Reads results/summary.json of a run whose step 2 pushes the socket
proximally (cup rigid displacement +z) and reports, per output time:
  * total axial load the socket transmits (cup axial constraint force),
  * the share carried by the distal end (cushion zone contact, axial part),
  * the share carried by the petals (petal-zone contact, axial part),
  * slip area fraction and peak contact pressure,
then interpolates the state at the target load (default 1.3 x 923 N) and
checks the design requirements:
  * the target load is reached (the grip holds it),
  * distal share < 15 %,
  * no element inversion (min J > 0) and the run converged that far.

usage: python scripts/stance_check.py RUN_DIR [--target 1846] [--distal-max 0.15] [--out docs/STANCE_CHECK.md]
"""
from __future__ import annotations

import argparse
import json
import os

import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--target", type=float, default=1200.0, help="axial load to hold, N (1.3 x 923 N stance peak)")
    ap.add_argument("--distal-max", type=float, default=0.15)
    ap.add_argument("--push-mm", type=float, default=None, help="step-2 cup displacement (default: from config)")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    s = json.load(open(os.path.join(a.run, "results", "summary.json")))
    st = json.load(open(os.path.join(a.run, "status.json")))
    cfg = json.load(open(os.path.join(a.run, "config.resolved.json")))
    T = json.load(open(os.path.join(a.run, "solver", "load_report.json")))["tension"]["applied_N"]
    push = a.push_mm or cfg["load_cases"]["body_weight"]["displacement"]
    rows = []
    for r in s["states"]:
        t = r["time"]
        z = r["contact"].get("per_zone", {})
        cup = r["equilibrium"]["rigid_bodies"].get("rigid_cup", {}).get("F", [0, 0, 0])[2]
        dist = z.get("skin_cushion_zone", {}).get("traction_resultant_N", [0, 0, 0])[2]
        pet = sum(v.get("traction_resultant_N", [0, 0, 0])[2] for k, v in z.items() if "petal" in k)
        wall = z.get("skin_cup_zone", {}).get("traction_resultant_N", [0, 0, 0])[2]
        rows.append({"t": t, "push_mm": max(0.0, (t - 1.0)) * push, "lace_N": T * min(t, 1.0),
                     "axial_N": abs(cup), "distal_N": abs(dist), "petals_N": abs(pet), "wall_N": abs(wall),
                     "slip": r["contact"].get("slip_fraction_of_contact_area", np.nan),
                     "pmax_kPa": 1000 * r["contact"].get("max_pressure_MPa", np.nan),
                     "min_J": r["deformation"]["min_J"], "imbalance": r["equilibrium"]["force_imbalance_relative"]})
    L = [f"# Stance-grip check: run `{os.path.basename(os.path.abspath(a.run))}`", "",
         f"Run state: **{st['state']}**, last converged time {st.get('last_converged_time')} "
         f"(t = 1 is the end of lace tightening, t = 2 is the full {push:g} mm push).",
         f"Lace tension: {T:.0f} N ({T / 4.4482:.0f} lbf). Target axial load: {a.target:.0f} N; distal share limit: {100 * a.distal_max:.0f}%.", "",
         "| t | push (mm) | axial load (N) | petals (N) | distal (N) | distal share | cup wall (N) | slip area | peak pressure (kPa) | min J |",
         "|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        if r["t"] >= 1 - 1e-9 or abs(r["t"] * 4 - round(r["t"] * 4)) < 1e-6:
            share = r["distal_N"] / r["axial_N"] if r["axial_N"] > 1 else np.nan
            L.append(f"| {r['t']:.2f} | {r['push_mm']:.1f} | {r['axial_N']:.0f} | {r['petals_N']:.0f} | {r['distal_N']:.0f} | "
                     f"{100 * share:.0f}% | {r['wall_N']:.0f} | {r['slip']:.2f} | {r['pmax_kPa']:.0f} | {r['min_J']:.2f} |")
    L.append("")
    push_rows = [r for r in rows if r["t"] > 1 + 1e-9]
    verdict = []
    if not push_rows:
        verdict.append("**No push state was reached**, so the grip was not tested.")
    else:
        F = np.array([r["axial_N"] for r in push_rows])
        reached = F.max() >= a.target
        if reached:
            k = int(np.argmax(F >= a.target))
            if k == 0:
                rA = rB = push_rows[0]; w = 0.0
            else:
                rA, rB = push_rows[k - 1], push_rows[k]
                w = (a.target - rA["axial_N"]) / max(rB["axial_N"] - rA["axial_N"], 1e-9)
            lerp = lambda key: rA[key] + w * (rB[key] - rA[key])
            share = lerp("distal_N") / a.target
            verdict.append(f"* **Target reached:** the socket carries {a.target:.0f} N after a {lerp('push_mm'):.1f} mm push.")
            verdict.append(f"* **Distal share at target:** {100 * share:.0f}% ({'PASS' if share < a.distal_max else 'FAIL'} "
                           f"against < {100 * a.distal_max:.0f}%). Petals: {lerp('petals_N'):.0f} N.")
            verdict.append(f"* **At target:** peak pressure ~{lerp('pmax_kPa'):.0f} kPa, slip area {lerp('slip'):.2f}, min J {lerp('min_J'):.2f}.")
        else:
            verdict.append(f"* **Target NOT reached:** the largest axial load is {F.max():.0f} N "
                           f"({100 * F.max() / a.target:.0f}% of the target) at a {push_rows[int(F.argmax())]['push_mm']:.1f} mm push.")
            if st["state"] != "converged":
                verdict.append("* The run stopped before the full push, so this is a lower bound on the grip, not its capacity.")
    L += ["## Verdict", ""] + verdict + ["",
          "Quasi-static axial push only: no gait bending moment, torsion or dynamics. "
          "Grip here means the model's frictional and geometric resistance to the socket sliding "
          "proximally, with penalty contact (which underestimates peak pressure)."]
    txt = "\n".join(L) + "\n"
    if a.out:
        open(a.out, "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
