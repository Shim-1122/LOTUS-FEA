"""Material parameter sets with units and provenance.

A material set is a YAML file (see configs/materials).  Parameters are
converted to internal units (MPa, deg, dimensionless) and every value keeps
its source key.  Scenario A (human tissue) and Scenario B (physical
prototype) are different files and are never mixed: a run uses exactly one
set, and the region->material mapping must reference keys of that set.

Supported FEBio constitutive models and the exact strain-energy functions
implemented by FEBio 4.13 (verified against FEBioMech source):

* "Holzapfel-Gasser-Ogden" (uncoupled)
      W = c/2 (I1~ - 3) + k1/(2 k2) sum_{a=1,2} [exp(k2 <E_a>^2) - 1] + U(J)
      E_a = kappa (I1~ - 3) + (1 - 3 kappa)(I4a~ - 1)  (fibres active if E_a >= 0)
      U(J) = k/2 [ (J^2 - 1)/2 - ln J ]            (pressure_model 3, set by HGO)
      fibre families at +/-gamma from local e1 in the e1-e2 plane.
* "Ogden" (uncoupled)
      W = sum_i c_i/m_i^2 (l1~^m_i + l2~^m_i + l3~^m_i - 3) + k/2 (ln J)^2
* "Mooney-Rivlin" (uncoupled)
      W = c1 (I1~ - 3) + c2 (I2~ - 3) + k/2 (ln J)^2
* "neo-Hookean" (compressible, E, v)
      W = mu/2 (I1 - 3) - mu ln J + lambda/2 (ln J)^2
* "trans iso Mooney-Rivlin" (uncoupled; passive fibre, optional)
* "rigid body"
"""
from __future__ import annotations

import copy

import yaml

from ..units import STRESS, quantity

ALLOWED = {
    "Holzapfel-Gasser-Ogden": {"c": "stress", "k1": "stress", "k2": "dimensionless", "kappa": "dimensionless",
                               "gamma": "angle", "k": "stress"},
    "Ogden": {"c1": "stress", "m1": "dimensionless", "c2": "stress", "m2": "dimensionless", "c3": "stress",
              "m3": "dimensionless", "k": "stress"},
    "Mooney-Rivlin": {"c1": "stress", "c2": "stress", "k": "stress"},
    "neo-Hookean": {"E": "stress", "v": "dimensionless"},
    "trans iso Mooney-Rivlin": {"c1": "stress", "c2": "stress", "c3": "stress", "c4": "dimensionless",
                                "c5": "stress", "lam_max": "dimensionless", "k": "stress"},
}


class MaterialError(ValueError):
    pass


def _conv(p: dict, kind: str) -> float:
    unit = p.get("unit", "-")
    if kind == "dimensionless":
        if unit not in ("-", "1", ""):
            raise MaterialError(f"dimensionless parameter given unit {unit!r}")
        return float(p["value"])
    return quantity({"value": p["value"], "unit": unit}, kind)


def load_material_set(path: str) -> dict:
    with open(path) as fh:
        raw = yaml.safe_load(fh)
    out = {"set_name": raw["set_name"], "scenario": raw["scenario"], "references": raw.get("references", {}),
           "source_file": path, "materials": {}}
    for key, m in raw["materials"].items():
        ftype = m["febio_type"]
        if ftype not in ALLOWED:
            raise MaterialError(f"{key}: unsupported FEBio material {ftype!r}")
        params, prov, assumed, approx = {}, {}, [], []
        for pname, p in m["parameters"].items():
            if pname not in ALLOWED[ftype]:
                raise MaterialError(f"{key}: parameter {pname!r} not valid for {ftype}")
            params[pname] = _conv(p, ALLOWED[ftype][pname])
            src = p.get("source")
            if not src:
                raise MaterialError(f"{key}.{pname}: missing source")
            if src != "assumption" and src not in out["references"]:
                raise MaterialError(f"{key}.{pname}: unknown source {src!r}")
            prov[pname] = {"value_input": p["value"], "unit_input": p.get("unit", "-"), "source": src,
                           "note": p.get("note", "")}
            if src in ("assumption", "PriorPackageAssumption"):
                assumed.append(pname)
            elif (out["references"].get(src) or {}).get("kind") == "approximation":
                approx.append(pname)
        pa = m.get("passive_anisotropy")
        if pa and pa.get("enabled"):
            raise MaterialError(f"{key}: passive anisotropy enabled but no calibrated parameters with provenance "
                                "are available; refusing to substitute a different model")
        unc = {}
        for pname, rng in (m.get("uncertainty") or {}).items():
            unit = m["parameters"][pname].get("unit", "-")
            kind = ALLOWED[ftype][pname]
            unc[pname] = [_conv({"value": v, "unit": unit}, kind) for v in rng]
        out["materials"][key] = {"febio_type": ftype, "params": params, "provenance": prov,
                                 "assumed_parameters": assumed, "approximated_parameters": approx,
                                 "calibration_range": m.get("calibration_range", ""),
                                 "uncertainty": unc, "region": m.get("region", key),
                                 "requires_user_calibration": bool(m.get("requires_user_calibration", False)),
                                 "model": m.get("model", "")}
    return out


def scaled(mset: dict, key: str, factor: float, params=None) -> dict:
    """Copy of a material set with stiffness parameters of one material
    scaled (used for sensitivity studies; recorded in the run manifest)."""
    new = copy.deepcopy(mset)
    m = new["materials"][key]
    stiff = params or [p for p, k in ALLOWED[m["febio_type"]].items() if k == "stress" and p in m["params"]]
    for p in stiff:
        m["params"][p] *= factor
    m.setdefault("sensitivity", []).append({"scaled": stiff, "factor": factor})
    return new


def small_strain_moduli(mat: dict) -> dict:
    """Initial shear and bulk moduli implied by a parameter set (MPa)."""
    t, p = mat["febio_type"], mat["params"]
    if t == "Holzapfel-Gasser-Ogden":
        return {"mu0_matrix": p["c"], "K": p["k"]}
    if t == "Ogden":
        mu = 0.5 * sum(p.get(f"c{i}", 0.0) for i in (1, 2, 3))  # FEBio convention c_i/m_i^2 -> mu0 = sum(c_i)/2
        return {"mu0": mu, "K": p["k"]}
    if t == "Mooney-Rivlin":
        return {"mu0": 2 * (p["c1"] + p["c2"]), "K": p["k"]}
    if t == "neo-Hookean":
        E, v = p["E"], p["v"]
        return {"mu0": E / (2 * (1 + v)), "K": E / (3 * (1 - 2 * v))}
    return {}


_ = STRESS  # unit table re-exported for documentation tools
