"""Muscle incompressibility check: volume lost by the muscle material under confined indentation.

A 20 mm cube of the muscle material (uncoupled Ogden, Sengeh 2016 shear parameters) on UT4
tetrahedra (the muscle element used in the limb) has its sides held in their normal direction
and its base fixed axially. A 10 x 10 mm indenter in the middle of the top face displaces it 5 mm
(25% strain). The material cannot spread sideways, so it can keep its volume only by flowing up
around the indenter. This is the situation of muscle squeezed between a pad and bone.

For each bulk/shear ratio k/c it reports the total volume change (from the deformed node
positions). The material is incompressible in the sense required by the team when the volume
change is well under 1%.

usage: python scripts/muscle_incompressibility_check.py [--ratios 100 1000 10000] [--febio febio4]
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import tempfile

import numpy as np

C1, M1 = 0.0095, 5.0           # muscle: Ranger 2022 shear modulus, Mukhina 2020 exponent (MPa, -)
L, N, INDENT = 20.0, 8, 5.0    # cube edge (mm), divisions, indentation (mm)


def mesh():
    g = np.linspace(0, L, N + 1)
    X = np.array([[x, y, z] for z in g for y in g for x in g])
    idx = lambda i, j, k: k * (N + 1) ** 2 + j * (N + 1) + i
    tets = []
    for k in range(N):
        for j in range(N):
            for i in range(N):
                v = [idx(i, j, k), idx(i + 1, j, k), idx(i + 1, j + 1, k), idx(i, j + 1, k),
                     idx(i, j, k + 1), idx(i + 1, j, k + 1), idx(i + 1, j + 1, k + 1), idx(i, j + 1, k + 1)]
                for t in [(0, 1, 3, 4), (1, 2, 3, 6), (1, 4, 5, 6), (3, 4, 6, 7), (1, 3, 4, 6)]:
                    tets.append([v[a] for a in t])
    return X, np.array(tets)


def feb(X, T, k, out):
    sets = {"fx": [q for q, p in enumerate(X) if p[0] in (0, L)],
            "fy": [q for q, p in enumerate(X) if p[1] in (0, L)],
            "fz": [q for q, p in enumerate(X) if p[2] == 0],
            "ind": [q for q, p in enumerate(X) if p[2] == L and 5 <= p[0] <= 15 and 5 <= p[1] <= 15]}
    ns = "".join(f'<NodeSet name="{n}">' + ",".join(str(i + 1) for i in ids) + "</NodeSet>" for n, ids in sets.items())
    nodes = "".join(f'<node id="{i + 1}">{p[0]},{p[1]},{p[2]}</node>' for i, p in enumerate(X))
    els = "".join(f'<elem id="{i + 1}">{",".join(str(a + 1) for a in t)}</elem>' for i, t in enumerate(T))
    return f'''<?xml version="1.0"?><febio_spec version="4.0"><Module type="solid"/>
<Control><analysis>STATIC</analysis><time_steps>10</time_steps><step_size>0.1</step_size>
<time_stepper type="default"><dtmin>0.0001</dtmin><dtmax>0.1</dtmax><max_retries>10</max_retries><opt_iter>10</opt_iter></time_stepper>
<solver type="solid"><max_refs>25</max_refs><qn_method type="BFGS"><max_ups>0</max_ups></qn_method></solver></Control>
<Material><material id="1" name="muscle" type="Ogden"><c1>{C1}</c1><m1>{M1}</m1><c2>{C1}</c2><m2>{-M1}</m2><k>{k}</k></material></Material>
<Mesh><Nodes name="all">{nodes}</Nodes><Elements type="tet4" name="b">{els}</Elements>{ns}</Mesh>
<MeshDomains><SolidDomain name="b" mat="muscle" type="ut4-solid"><alpha>0.05</alpha><iso_stab>0</iso_stab></SolidDomain></MeshDomains>
<Boundary><bc type="zero displacement" node_set="fx"><x_dof>1</x_dof></bc><bc type="zero displacement" node_set="fy"><y_dof>1</y_dof></bc>
<bc type="zero displacement" node_set="fz"><z_dof>1</z_dof></bc>
<bc type="prescribed displacement" node_set="ind"><dof>z</dof><value lc="1">{-INDENT}</value></bc></Boundary>
<LoadData><load_controller id="1" type="loadcurve"><points><pt>0,0</pt><pt>1,1</pt></points></load_controller></LoadData>
<Output><logfile><node_data data="x;y;z" file="{out}"/></logfile></Output></febio_spec>'''


def volume(x, T):
    a, b, c, d = (x[T[:, i]] for i in range(4))
    return np.einsum("ij,ij->i", np.cross(b - a, c - a), d - a).sum() / 6


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ratios", type=float, nargs="+", default=[100, 1000, 10000])
    ap.add_argument("--febio", default=os.environ.get("RLFEA_FEBIO", shutil.which("febio4") or "febio4"))
    a = ap.parse_args()
    X, T = mesh()
    V0 = volume(X, T)
    print("| k / c | bulk modulus k (MPa) | Poisson ratio | reached | volume change |")
    print("|---|---|---|---|---|")
    for r in a.ratios:
        with tempfile.TemporaryDirectory() as d:
            open(os.path.join(d, "m.feb"), "w").write(feb(X, T, r * C1, "x.txt"))
            subprocess.run([a.febio, "-i", "m.feb", "-silent"], cwd=d, capture_output=True)
            ok = "N O R M A L" in open(os.path.join(d, "m.log")).read()
            blocks = open(os.path.join(d, "x.txt")).read().split("*Step")[1:]
            rows = [l.split() for l in blocks[-1].splitlines() if l and l[0].isdigit()]
            x = np.array([[float(v) for v in row[1:]] for row in rows])
            nu = (3 * r - 2) / (6 * r + 2)
            print(f"| {r:g} | {r * C1:.3g} | {nu:.4f} | {'full indent' if ok else 'stopped early'} | "
                  f"{100 * (volume(x, T) / V0 - 1):+.2f}% |")


if __name__ == "__main__":
    main()
