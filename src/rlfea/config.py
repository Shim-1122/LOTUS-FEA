"""Run configuration: loading, unit conversion and validation.

Configuration files are YAML.  Every dimensional input is written as
``{value: x, unit: u}``; bare numbers are only accepted for dimensionless
quantities.  ``resolve()`` returns a plain dict in internal units (mm, N, MPa,
deg) and a list of validation messages.  Errors stop the run; warnings are
recorded in the run manifest and in the report.
"""
from __future__ import annotations

import copy
import hashlib
import json
import os

import yaml

from .units import quantity

DEFAULTS = {
    "layers": {"skin_thickness": {"value": 1.5, "unit": "mm"}, "fat_thickness": {"value": 3.0, "unit": "mm"},
               "max_fraction_skin": 0.2, "max_fraction_fat": 0.2, "min_layer": {"value": 0.15, "unit": "mm"},
               "normal_smoothing": 3, "thickness_smoothing": 5, "skin_fiber_reference": "circumferential",
               "muscle_formulation": "ut4", "ut4_alpha": 0.05, "muscle_fiber": False, "optimize": True},
    "supports": {"proximal_cut": "axial_roller"},
    "contact": {"type": "sliding-elastic", "friction": 0.5, "augmented_lagrangian": False, "penalty_scale": 1.0,
                "auto_penalty": True, "two_pass": False, "search_tol": 0.01, "search_radius": {"value": 5.0, "unit": "mm"},
                "bonded": False},
    "solver": {"backend": "febio", "dtol": 0.01, "etol": 0.001, "rtol": 0.001, "lstol": 0.9, "max_refs": 50,
               "qn_method": "auto", "max_ups": 10, "n_steps": 20, "dtmin": 1e-5, "max_retries": 10, "opt_iter": 10,
               "threads": None, "timeout_s": None},
    "cable": {"guide_friction": 0.0, "route_height_fraction": 0.75, "guide_inset": {"value": 4.0, "unit": "deg"},
              "dial_petal": 0, "tension_is_example": True, "rated_tension": None,
              "rating_fractions": [0.25, 0.5, 0.75, 1.0]},
    "load_cases": {"tightening": {"enabled": True}, "pullout": {"enabled": False},
                   "body_weight": {"enabled": False}, "bending": {"enabled": False},
                   "torsion": {"enabled": False}, "suction": {"enabled": False}},
}


class ConfigError(ValueError):
    pass


def _merge(base, over):
    out = copy.deepcopy(base)
    for k, v in (over or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


def load(path: str) -> dict:
    with open(path) as fh:
        raw = yaml.safe_load(fh)
    raw["_config_path"] = os.path.abspath(path)
    raw["_config_sha256"] = hashlib.sha256(open(path, "rb").read()).hexdigest()
    return raw


def _q(d, key, kind, required=True, default=None):
    if key not in d or d[key] is None:
        if required and default is None:
            raise ConfigError(f"missing required quantity {key!r}")
        return default
    return quantity(d[key], kind)


def resolve(raw: dict) -> tuple[dict, list[str]]:
    """Merge defaults, convert units and validate.  Returns (cfg, warnings)."""
    warn: list[str] = []
    cfg = _merge(DEFAULTS, {k: v for k, v in raw.items()})
    base = os.path.dirname(raw.get("_config_path", "."))
    L = cfg["layers"]
    for k in ("skin_thickness", "fat_thickness", "min_layer"):
        L[k] = quantity(L[k], "length")
    if L["skin_fiber_reference"] not in ("circumferential", "longitudinal"):
        raise ConfigError("layers.skin_fiber_reference must be circumferential|longitudinal")
    if L["muscle_formulation"] not in ("ut4", "standard"):
        raise ConfigError("layers.muscle_formulation must be ut4|standard")
    g = cfg["geometry"]
    if g["source"] not in ("synthetic", "files"):
        raise ConfigError("geometry.source must be synthetic|files")
    for k in ("skin_size", "bone_size"):
        if k in g.get("remesh", {}):
            g["remesh"][k] = quantity(g["remesh"][k], "length")
    if g["source"] == "synthetic":
        g.setdefault("synthetic", {})
        g["synthetic"]["length"] = quantity(g["synthetic"].get("length", {"value": 200, "unit": "mm"}), "length")
    else:
        seg = g.get("segmented_layers") or {}
        for k in seg:
            if k not in ("skin_inner", "fat_inner", "note"):
                raise ConfigError(f"geometry.segmented_layers.{k}: allowed keys are skin_inner, fat_inner, note")
        if seg and not ("skin_inner" in seg and "fat_inner" in seg):
            raise ConfigError("geometry.segmented_layers needs both skin_inner and fat_inner")
        if seg:
            warn.append("Layer thicknesses are measured from a supplied segmentation whose provenance and accuracy "
                        "are not verified by this software.")
        for b in [g["skin"]] + g.get("bones", []) + g.get("ignored_files", []) + \
                [v for k, v in seg.items() if k != "note"]:
            b["path"] = b["path"] if os.path.isabs(b["path"]) else os.path.normpath(os.path.join(base, b["path"]))
            if not os.path.exists(b["path"]):
                raise ConfigError(f"geometry file not found: {b['path']}")
        if g.get("patient_data", True):
            warn.append("Supplied geometry is treated as private patient-derived data; it is never written into the repository.")
    m = cfg["materials"]
    m["set"] = m["set"] if os.path.isabs(m["set"]) else os.path.normpath(os.path.join(base, m["set"]))
    if m.get("bones", "rigid") not in ("rigid", "deformable"):
        raise ConfigError("materials.bones must be rigid|deformable")
    C = cfg["contact"]
    if not 0 <= float(C["friction"]) <= 1.5:
        raise ConfigError("contact.friction out of range")
    C["search_radius"] = quantity(C["search_radius"], "length")
    if C.get("bonded"):
        warn.append("contact.bonded = true: tissue-socket interfaces are TIED by explicit user choice (no sliding/separation).")
    S = cfg.get("socket", {"enabled": False})
    if S.get("enabled"):
        for k in ("z_base", "z_top", "posterior_clearance_drop", "tip_drop", "pad_thickness", "shell_thickness",
                  "initial_gap", "element_size", "cushion_thickness", "cup_thickness", "contact_zone_margin"):
            if k in S:
                S[k] = quantity(S[k], "length")
        for k in ("gap_deg", "posterior_center_deg", "posterior_clearance_halfwidth_deg"):
            if k in S and isinstance(S[k], dict):
                S[k] = quantity(S[k], "angle")
        S["stabilizer_k"] = quantity(S.get("stabilizer_k", {"value": 0.05, "unit": "N/mm"}), "stiffness")
        warn.append("Socket geometry is a PROVISIONAL parametric model generated from configuration values; "
                    "no verified current socket CAD was supplied (V6 excluded).")
    cb = cfg["cable"]
    if "strand_tension" in cb:
        cb["strand_tension"] = quantity(cb["strand_tension"], "force")
    cb["guide_inset"] = quantity(cb["guide_inset"], "angle")
    rt = cb.get("rated_tension")
    if rt:
        for req in ("value", "unit", "source", "meaning"):
            if req not in rt or rt[req] in (None, ""):
                raise ConfigError(f"cable.rated_tension requires '{req}' (verified manufacturer rating with source and meaning)")
        if rt["meaning"] not in ("strand_tension",):
            raise ConfigError("cable.rated_tension.meaning must be 'strand_tension' (the quantity the model applies); "
                              "other ratings (e.g. dial torque, lace break strength) must be converted with a documented method first")
        cb["rated_tension_N"] = quantity({"value": rt["value"], "unit": rt["unit"]}, "force")
    else:
        warn.append("No verified RevoFit rated strand tension supplied: rated-tension cases (25/50/75/100%) are NOT run.")
    if cb.get("tension_is_example", True) and "strand_tension" in cb:
        warn.append(f"Strand tension {cb['strand_tension']:.3g} N is an EXAMPLE value, not a verified rating.")
    lc = cfg["load_cases"]
    for k in ("pullout", "body_weight"):
        if lc[k].get("enabled"):
            lc[k]["displacement"] = quantity(lc[k]["displacement"], "length")
    for k in ("bending", "torsion"):
        if lc[k].get("enabled"):
            lc[k]["rotation"] = quantity(lc[k]["rotation"], "angle")
    if lc["suction"].get("enabled"):
        sc = lc["suction"]
        for req in ("pressure", "justification", "sealing_assumption"):
            if not sc.get(req):
                raise ConfigError(f"suction enabled but '{req}' is missing; suction is only modelled with a justified "
                                  "pressure boundary condition and an explicit sealing assumption")
        sc["pressure"] = quantity(sc["pressure"], "stress")
    if any(lc[k].get("enabled") for k in ("pullout", "body_weight", "bending", "torsion")) and not S.get("enabled"):
        raise ConfigError("pullout/body-weight/bending/torsion require the socket model (they act through the socket)")
    return cfg, warn


def config_digest(cfg: dict) -> str:
    return hashlib.sha256(json.dumps(cfg, sort_keys=True, default=str).encode()).hexdigest()
