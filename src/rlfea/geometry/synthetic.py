"""Synthetic transtibial residual-limb geometry (no patient data).

Frame: +z proximal (limb long axis), +y anterior, +x lateral; origin at the
distal skin apex.  All dimensions are illustrative round numbers chosen to
resemble a mid-length transtibial residuum; they are NOT anatomical
measurements and must never be presented as patient-specific.

Bodies: skin envelope (closed, flat proximal cut), tibia (with anterior crest
placed close to the skin to create a thin-tissue region) and fibula
(posterolateral, shorter).  The bones end a few millimetres below the
proximal cut so that all bodies are closed and nested.
"""
from __future__ import annotations

import numpy as np

from .surface import TriSurface


def _tube(zs, ring_fn, n_theta, cap_bottom=True, cap_top=True):
    """Closed tube from rings: ring_fn(z, theta) -> (x, y). Bottom/top capped
    with a single pole (bottom) or flat fan (top)."""
    th = np.linspace(0, 2 * np.pi, n_theta, endpoint=False)
    rings = []
    for z in zs:
        x, y = ring_fn(z, th)
        rings.append(np.c_[x, y, np.full(n_theta, z)])
    V = np.vstack(rings)
    F = []
    nz = len(zs)
    for i in range(nz - 1):
        a = i * n_theta
        b = (i + 1) * n_theta
        for j in range(n_theta):
            j1 = (j + 1) % n_theta
            F.append([a + j, a + j1, b + j1])
            F.append([a + j, b + j1, b + j])
    if cap_bottom:
        pole = len(V)
        x, y = ring_fn(zs[0], th)
        V = np.vstack([V, [np.mean(x), np.mean(y), zs[0]]])
        for j in range(n_theta):
            F.append([pole, (j + 1) % n_theta, j])
    if cap_top:
        top = len(V)
        last = (nz - 1) * n_theta
        V = np.vstack([V, [V[last:last + n_theta, 0].mean(), V[last:last + n_theta, 1].mean(), zs[-1]]])
        for j in range(n_theta):
            F.append([top, last + j, last + (j + 1) % n_theta])
    return V, np.array(F)


def _dome(z, z0, r_scale, n=2.0):
    """Smooth rounding factor approaching 0 at z=z0 (apex) over length r_scale."""
    t = np.clip((z - z0) / r_scale, 0, 1)
    return (1 - (1 - t) ** n) ** (1 / n)


def synthetic_limb(length=200.0, distal_round=45.0, n_theta=72, dz=2.5):
    """Return dict of TriSurface: skin, tibia, fibula."""
    # ---- skin envelope -------------------------------------------------
    def skin_ring(z, th):
        s = z / length
        a = 46 + 10 * s - 6 * s ** 2          # mediolateral semi-axis (mm)
        b = 44 + 12 * s - 5 * s ** 2          # anteroposterior semi-axis (mm)
        f = _dome(z, 0.0, distal_round)
        p = 2.3                               # superellipse exponent
        c, sn = np.cos(th), np.sin(th)
        r = (np.abs(c / a) ** p + np.abs(sn / b) ** p) ** (-1 / p)
        # slightly flattened anterior tibial face, fuller posterior calf
        r = r * (1 + 0.06 * np.clip(-sn, 0, 1) - 0.03 * np.clip(sn, 0, 1) ** 2)
        return f * r * c, f * r * sn - 4 * s * np.clip(-sn, 0, 1) * 0

    zs = np.r_[np.linspace(0.0, distal_round, int(distal_round / 1.0) + 1)[1:],
               np.arange(distal_round + dz, length + 1e-9, dz)]
    zs = np.unique(np.r_[0.4, zs])
    V, F = _tube(zs, skin_ring, n_theta, cap_bottom=True, cap_top=True)
    skin = TriSurface("skin_envelope", V, F, source="rlfea.geometry.synthetic", source_unit="mm")
    skin.log("generate", generator="synthetic_limb", length_mm=length, note="SYNTHETIC, not patient data")

    # ---- tibia ------------------------------------------------------------
    z_t0, z_t1 = 62.0, length - 4.0          # residual tibia from distal cut to just below proximal cut

    def tib_ring(z, th):
        s = (z - z_t0) / (z_t1 - z_t0)
        flare = np.clip((s - 0.72) / 0.28, 0, 1) ** 2   # proximal metaphyseal flare
        a = 11 + 22 * flare
        b = 10 + 13 * flare
        f = _dome(z, z_t0, 9.0)
        c, sn = np.cos(th), np.sin(th)
        # triangular-ish section with anterior crest (+y)
        r = 1.0 / np.sqrt((c / a) ** 2 + (sn / b) ** 2)
        r = r * (1 + 0.28 * np.clip(sn, 0, 1) ** 6)
        y_off = 14 - 6 * flare                               # anterior placement: thin tissue over crest
        return f * r * c - 4.0, f * r * sn + y_off

    zt = np.unique(np.r_[np.linspace(z_t0 + 0.3, z_t0 + 9, 10), np.arange(z_t0 + 10, z_t1 + 1e-9, 3.0), z_t1])
    Vt, Ft = _tube(zt, tib_ring, 48)
    tibia = TriSurface("tibia", Vt, Ft, source="rlfea.geometry.synthetic", source_unit="mm")
    tibia.log("generate", generator="synthetic_limb", note="SYNTHETIC tibia")

    # ---- fibula -----------------------------------------------------------
    z_f0, z_f1 = 80.0, length - 22.0

    def fib_ring(z, th):
        s = (z - z_f0) / (z_f1 - z_f0)
        a = 6.0 + 3.0 * np.clip((s - 0.85) / 0.15, 0, 1)
        f = _dome(z, z_f0, 5.0)
        return f * a * np.cos(th) + 27.0, f * a * 0.9 * np.sin(th) - 14.0

    zf = np.unique(np.r_[np.linspace(z_f0 + 0.2, z_f0 + 5, 6), np.arange(z_f0 + 6, z_f1 + 1e-9, 3.0), z_f1])
    Vf, Ff = _tube(zf, fib_ring, 24)
    fibula = TriSurface("fibula", Vf, Ff, source="rlfea.geometry.synthetic", source_unit="mm")
    fibula.log("generate", generator="synthetic_limb", note="SYNTHETIC fibula")
    return {"skin": skin, "tibia": tibia, "fibula": fibula}


def write_stl(surf: TriSurface, path: str):
    import trimesh
    trimesh.Trimesh(surf.vertices, surf.faces, process=False).export(path)
