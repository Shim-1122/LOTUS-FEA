"""Element quality and validity metrics (tet4, penta6, hex8)."""
from __future__ import annotations

import numpy as np


def tet_volume(X, T):
    a, b, c, d = (X[T[:, i]] for i in range(4))
    return np.einsum("ij,ij->i", b - a, np.cross(c - a, d - a)) / 6.0


def tet_quality(X, T) -> dict:
    """Signed volume, radius ratio (3 r_in / r_circ, 1 = regular) and the
    minimum scaled Jacobian (1 = regular tet)."""
    P = X[T]
    vol = tet_volume(X, T)
    # face areas
    faces = [(1, 2, 3), (0, 2, 3), (0, 1, 3), (0, 1, 2)]
    A = sum(0.5 * np.linalg.norm(np.cross(P[:, j] - P[:, i], P[:, k] - P[:, i]), axis=1) for i, j, k in faces)
    r_in = 3 * np.abs(vol) / A
    # circumradius
    a = P[:, 1] - P[:, 0]; b = P[:, 2] - P[:, 0]; c = P[:, 3] - P[:, 0]
    num = (np.einsum("ij,ij->i", a, a)[:, None] * np.cross(b, c) + np.einsum("ij,ij->i", b, b)[:, None] * np.cross(c, a)
           + np.einsum("ij,ij->i", c, c)[:, None] * np.cross(a, b))
    r_circ = np.linalg.norm(num, axis=1) / (12 * np.abs(vol) + 1e-300)
    rr = 3 * r_in / (r_circ + 1e-300)
    # scaled jacobian at each corner
    sj = np.full(len(T), np.inf)
    for i in range(4):
        o = P[:, i]
        others = [P[:, j] - o for j in range(4) if j != i]
        if i % 2 == 1:
            others = [others[0], others[2], others[1]]
        det = np.einsum("ij,ij->i", others[0], np.cross(others[1], others[2]))
        ln = np.prod([np.linalg.norm(e, axis=1) for e in others], axis=0)
        sj = np.minimum(sj, det / (ln + 1e-300))
    sj *= np.sqrt(2.0)  # normalise so that a regular tet scores 1
    return {"volume": vol, "radius_ratio": rr, "scaled_jacobian": sj}


def penta6_corner_jacobians(X, W):
    """Jacobian determinants at the 6 corners of linear wedges with nodes
    (1,2,3) bottom and (4,5,6) top (FEBio/Abaqus ordering).  Positive for a
    valid element.  Returned scaled by edge lengths (dimensionless)."""
    P = X[W]
    out = []
    for (i, j, k) in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
        for layer in (0, 3):
            o = P[:, i + layer]
            e1 = P[:, j + layer] - o
            e2 = P[:, k + layer] - o
            e3 = P[:, i + 3] - P[:, i]          # through-thickness edge
            det = np.einsum("ij,ij->i", np.cross(e1, e2), e3)
            ln = np.linalg.norm(e1, axis=1) * np.linalg.norm(e2, axis=1) * np.linalg.norm(e3, axis=1)
            out.append(det / (ln + 1e-300))
    return np.array(out).T  # (n, 6)


def hex8_corner_jacobians(X, H):
    P = X[H]
    # corner -> (neighbours along r, s, t) for standard hex ordering
    nb = {0: (1, 3, 4), 1: (2, 0, 5), 2: (3, 1, 6), 3: (0, 2, 7),
          4: (7, 5, 0), 5: (4, 6, 1), 6: (5, 7, 2), 7: (6, 4, 3)}
    out = []
    for c, (i, j, k) in nb.items():
        e1 = P[:, i] - P[:, c]; e2 = P[:, j] - P[:, c]; e3 = P[:, k] - P[:, c]
        det = np.einsum("ij,ij->i", np.cross(e1, e2), e3)
        if c >= 4:
            det = -det
        ln = np.linalg.norm(e1, axis=1) * np.linalg.norm(e2, axis=1) * np.linalg.norm(e3, axis=1)
        out.append(det / (ln + 1e-300))
    return np.array(out).T


def summarize(values: np.ndarray) -> dict:
    v = np.asarray(values)
    if v.size == 0:
        return {}
    return {"min": float(v.min()), "p01": float(np.percentile(v, 1)), "median": float(np.median(v)),
            "max": float(v.max())}
