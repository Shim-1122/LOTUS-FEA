"""Cable (lace) routing and guide-force derivation.

The lace is an ordered list of points, each attached to a body (a petal or
the cup).  Consecutive points on the SAME body are joined by a segment that
runs inside that body (a channel); its tension is an internal force of that
body.  Consecutive points on DIFFERENT bodies form a free span; the span
pulls each end point toward the other with the span tension.

Free-body consequences (verified in tests):
* sum of all span forces on all bodies = 0 (the lace is internal to the
  socket assembly) and sum of moments = 0 (each span's pair is collinear);
* the force and moment on body B are sum over its span ends of T*u and
  (r - c_B) x T*u.

Strand tension is applied at both lace ends by the reel (dial).  Optional
guide friction uses the capstan relation T_out = T_in exp(-mu * wrap) with the
wrap angle at each direction change; tension then decays from each end toward
the middle of the lace (T = max of the two end-driven values), which is the
quasi-static tightening state.  mu_guide = 0 is a frictionless assumption.

The strand tension is NOT multiplied by any factor and NOT applied
independently to each pad: pads receive load only through resolved contact
with the rigid petals (or, in the explicitly labelled prescribed-traction
experiment, through the petal moment balance).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class RoutePoint:
    body: str
    xyz: np.ndarray
    label: str = ""


def _unit(v):
    n = np.linalg.norm(v)
    if n == 0:
        raise ValueError("zero-length lace span (coincident guide points)")
    return v / n


def lace_tensions(points: list[RoutePoint], T0: float, mu_guide: float = 0.0) -> np.ndarray:
    """Tension in each segment k (between point k and k+1)."""
    n = len(points)
    if n < 2:
        raise ValueError("route needs at least two points")
    dirs = [_unit(points[k + 1].xyz - points[k].xyz) for k in range(n - 1)]
    wrap = np.zeros(n)
    for k in range(1, n - 1):
        c = np.clip(np.dot(dirs[k - 1], dirs[k]), -1, 1)
        wrap[k] = np.arccos(c)
    # forward from start, backward from end
    cum_f = np.concatenate([[0.0], np.cumsum(wrap[1:n - 1])])            # wrap before segment k
    cum_b = np.concatenate([np.cumsum(wrap[1:n - 1][::-1])[::-1], [0.0]])  # wrap after segment k
    Tf = T0 * np.exp(-mu_guide * cum_f)
    Tb = T0 * np.exp(-mu_guide * cum_b)
    return np.maximum(Tf, Tb)


def body_loads(points: list[RoutePoint], T0: float, com: dict[str, np.ndarray], mu_guide: float = 0.0) -> dict:
    """Return per-body point forces, resultant force and moment about the
    body's reference point (centre of mass used by the solver)."""
    T = lace_tensions(points, T0, mu_guide)
    out: dict = {b: {"point_forces": [], "F": np.zeros(3), "M_com": np.zeros(3)} for b in {p.body for p in points}}
    spans = []
    for k in range(len(points) - 1):
        a, b = points[k], points[k + 1]
        if a.body == b.body:
            continue
        u = _unit(b.xyz - a.xyz)
        fa = T[k] * u          # on body a at point a, toward b
        fb = -T[k] * u         # on body b at point b, toward a
        spans.append({"from": a.label, "to": b.label, "tension_N": float(T[k]),
                      "length_mm": float(np.linalg.norm(b.xyz - a.xyz))})
        for p, f in ((a, fa), (b, fb)):
            rec = out[p.body]
            rec["point_forces"].append({"label": p.label, "xyz": p.xyz.tolist(), "force_N": f.tolist()})
            rec["F"] = rec["F"] + f
            rec["M_com"] = rec["M_com"] + np.cross(p.xyz - com[p.body], f)
    total_F = sum((v["F"] for v in out.values()), np.zeros(3))
    total_M0 = np.zeros(3)
    for b, v in out.items():
        for pf in v["point_forces"]:
            total_M0 += np.cross(np.array(pf["xyz"]), np.array(pf["force_N"]))
    return {"bodies": out, "spans": spans, "segment_tensions_N": T.tolist(),
            "sum_force_N": total_F.tolist(), "sum_moment_about_origin_Nmm": total_M0.tolist(),
            "T0_N": T0, "mu_guide": mu_guide}


def default_route(petals, z_frac: float = 0.75, inset_deg: float = 4.0, dial_petal: int = 0,
                  shell_radius_fn=None) -> list[RoutePoint]:
    """Single-row circumferential lace crossing every inter-petal gap.

    dial (center of petal `dial_petal`) -> its right guide -> next petal left
    guide -> right guide -> ... -> dial petal left guide -> dial.
    PROVISIONAL: guide positions follow the parametric petals, not a verified
    RevoFit routing.
    """
    n = len(petals)
    order = [(dial_petal + k) % n for k in range(n)]

    # one common lace height (fraction of the LOWEST petal height) so the
    # default row is horizontal; posterior petals trimmed for hamstring
    # clearance would otherwise tilt every span.
    z_row = petals[0].z_base + z_frac * (min(q.z_top_center for q in petals) - petals[0].z_base)

    def pt(p, dth):
        th = np.deg2rad(p.theta_center_deg + dth)
        z = z_row
        r = shell_radius_fn(p, th, z)
        return np.array([r * np.cos(th), r * np.sin(th), z])

    dp = petals[order[0]]
    pts = [RoutePoint(dp.name, pt(dp, 0.0), f"{dp.name}:dial_out")]
    pts.append(RoutePoint(dp.name, pt(dp, dp.half_width_deg - inset_deg), f"{dp.name}:guide_ccw"))
    for idx in order[1:]:
        p = petals[idx]
        pts.append(RoutePoint(p.name, pt(p, -(p.half_width_deg - inset_deg)), f"{p.name}:guide_cw"))
        pts.append(RoutePoint(p.name, pt(p, p.half_width_deg - inset_deg), f"{p.name}:guide_ccw"))
    pts.append(RoutePoint(dp.name, pt(dp, -(dp.half_width_deg - inset_deg)), f"{dp.name}:guide_cw"))
    pts.append(RoutePoint(dp.name, pt(dp, 0.0) + np.array([0, 0, 1e-3]), f"{dp.name}:dial_in"))
    return pts


def route_from_config(spec: list[dict]) -> list[RoutePoint]:
    """Explicit routing: [{body, xyz_mm:[x,y,z], label}] in the model frame."""
    return [RoutePoint(s["body"], np.asarray(s["xyz_mm"], float), s.get("label", f"p{i}")) for i, s in enumerate(spec)]


def petal_pad_force_moment_balance(petal, body_load: dict, pad_centroid: np.ndarray, pad_normal: np.ndarray) -> float:
    """Simplified prescribed-traction model: rigid petal, frictionless pin
    hinge, pad resultant acting along -pad_normal at the pad centroid.
    Returns the pad normal force magnitude from moment balance about the
    hinge axis (N).  Used ONLY for the labelled prescribed-traction
    experiment; it cannot establish grip.

    pad_normal is the OUTWARD (radial) unit normal of the pad; the skin
    pushes the petal along +pad_normal with magnitude N (N > 0 = compression).
    Moment balance about the hinge axis a:  M_cable + N ((c - o) x n) . a = 0.
    """
    a = petal.hinge_axis / np.linalg.norm(petal.hinge_axis)
    M_cable = 0.0
    for pf in body_load["point_forces"]:
        M_cable += np.dot(np.cross(np.array(pf["xyz"]) - petal.hinge_origin, np.array(pf["force_N"])), a)
    lever = np.dot(np.cross(pad_centroid - petal.hinge_origin, pad_normal), a)
    if abs(lever) < 1e-9:
        raise ValueError("pad resultant passes through the hinge axis")
    return -M_cable / lever
