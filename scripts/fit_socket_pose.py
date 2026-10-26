"""Fit the pose of a socket CAD to the limb skin (rigid: x, y, z shift and a
rotation about the socket axis) and print the 4x4 transforms used in
configs/lotus_v5_socket.yaml.

The objective is the mean squared distance from points sampled on the inner
(skin-facing) faces of the silicone pads to the outer skin surface, i.e. the
"closed" CAD state should sit uniformly close to the skin.  The fitted
clearance distribution is printed so the result can be judged.

usage: python scripts/fit_socket_pose.py SOCKET.step data/limb/skin.stl [--guess DX DY]
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))


def Rz(a):
    c, s = np.cos(a), np.sin(a)
    M = np.eye(4)
    M[:2, :2] = [[c, -s], [s, c]]
    return M


def main():
    import trimesh
    from scipy.optimize import minimize
    from rlfea.geometry.cad_socket import tessellate_step

    ap = argparse.ArgumentParser()
    ap.add_argument("step")
    ap.add_argument("skin_stl")
    ap.add_argument("--guess", nargs=2, type=float, default=None, help="initial limb-centre offset (mm) in socket x, y")
    ap.add_argument("--samples", type=int, default=400)
    a = ap.parse_args()
    parts = tessellate_step(a.step, 3.0)
    skin = trimesh.load(a.skin_stl, force="mesh")
    skin.merge_vertices()
    skin = sorted(skin.split(only_watertight=False), key=lambda c: -abs(c.volume))[0]   # outer skin shell
    pq = trimesh.proximity.ProximityQuery(skin)
    rng = np.random.default_rng(0)
    P = []
    for k, (V, F) in parts.items():
        if k == "_meta" or not k.endswith("silicone_pad"):
            continue
        m = trimesh.Trimesh(V, F, process=False)
        c, n = m.triangles_center, m.face_normals
        rh = c[:, :2] / np.linalg.norm(c[:, :2], axis=1)[:, None]
        inner = (n[:, 0] * rh[:, 0] + n[:, 1] * rh[:, 1]) < -0.6
        w = m.area_faces[inner]
        P.append(c[inner][rng.choice(inner.sum(), a.samples, replace=False, p=w / w.sum())])
    P = np.vstack(P)
    sk_c = skin.vertices.mean(0)
    guess = np.array(a.guess) if a.guess else -sk_c[:2]

    def sd(x):
        dx, dy, dz, rz = x
        c, s = np.cos(rz), np.sin(rz)
        Q = P.copy()
        Q[:, 0], Q[:, 1] = c * P[:, 0] - s * P[:, 1], s * P[:, 0] + c * P[:, 1]
        return pq.signed_distance(Q - np.array([dx, dy, dz]))   # + = pad point inside the skin

    r = minimize(lambda x: np.mean(sd(x) ** 2), np.r_[-guess, 0.0, 0.0], method="Nelder-Mead",
                 options={"xatol": 0.05, "fatol": 1e-3, "maxiter": 600})
    dx, dy, dz, phi = r.x
    d = -sd(r.x)
    print(f"pose: limb = Rz({np.degrees(phi):.3f} deg) * socket - ({dx:.3f}, {dy:.3f}, {dz:.3f}) mm")
    print("pad-to-skin clearance (mm, + = gap): p5 %.2f  median %.2f  p95 %.2f  min %.2f  max %.2f"
          % (*np.percentile(d, [5, 50, 95]), d.min(), d.max()))
    T = np.eye(4); T[:3, 3] = [dx, dy, dz]
    src_to_sock = Rz(-phi) @ T
    V = np.asarray(skin.vertices)
    z0 = ((src_to_sock[:3, :3] @ V.T).T + src_to_sock[:3, 3])[:, 2].min()
    Z = np.eye(4); Z[2, 3] = -z0
    sock_to_model = Z @ Rz(np.pi)
    np.set_printoptions(precision=6, suppress=True)
    print("geometry.alignment.source_to_model =", (sock_to_model @ src_to_sock).round(6).tolist())
    print("socket.cad.cad_to_model            =", sock_to_model.round(6).tolist())


if __name__ == "__main__":
    main()
