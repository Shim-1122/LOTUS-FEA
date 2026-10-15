"""Reader for STEP AP242 tessellated geometry (TESSELLATED_SHELL /
TRIANGULATED_FACE / COMPLEX_TRIANGULATED_FACE).

OpenCASCADE builds bundled with gmsh do not reliably import AP242
tessellated bodies (the supplied bone files crash the gmsh OCC importer), so
these entities are parsed directly.  Only the subset of ISO 10303-42/-242
needed for triangle meshes is implemented:

* COORDINATES_LIST('', npts, ((x,y,z),...))
* TRIANGULATED_FACE(name, coords, pnmax, normals, geometric_link, pnindex, triangles)
* COMPLEX_TRIANGULATED_FACE(name, coords, pnmax, normals, geometric_link,
  pnindex, triangle_strips, triangle_fans)
* the length unit from the GLOBAL_UNIT_ASSIGNED_CONTEXT (SI_UNIT with optional
  prefix, or CONVERSION_BASED_UNIT for inches/feet).

Anything else in the file is ignored.  The reader never modifies the source
file; it returns vertices in the file's own length unit plus the scale to mm.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

import numpy as np

_ENTITY_RE = re.compile(r"#(\d+)\s*=\s*(.*?);\s*(?=#\d+\s*=|ENDSEC)", re.S)

_SI_PREFIX = {
    None: 1.0, "MILLI": 1e-3, "CENTI": 1e-2, "DECI": 1e-1, "KILO": 1e3,
    "MICRO": 1e-6,
}


@dataclass
class TessellatedBody:
    name: str
    vertices: np.ndarray            # (n,3) in file units
    faces: np.ndarray               # (m,3) int
    unit_name: str
    scale_to_mm: float
    n_face_entities: int
    warnings: list = field(default_factory=list)


def _split_top(s: str) -> list[str]:
    """Split a STEP argument list at top-level commas."""
    out, depth, cur, in_str = [], 0, [], False
    for ch in s:
        if ch == "'" :
            in_str = not in_str
        if not in_str:
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
            elif ch == "," and depth == 0:
                out.append("".join(cur).strip()); cur = []
                continue
        cur.append(ch)
    if cur:
        out.append("".join(cur).strip())
    return out


def _args(body: str) -> list[str]:
    i = body.index("(")
    j = body.rindex(")")
    return _split_top(body[i + 1:j])


def _nums(s: str) -> np.ndarray:
    return np.array(re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", s), dtype=float)


def _nested_int_lists(s: str) -> list[list[int]]:
    """Parse '((1,2,3),(4,5,6,7))' -> [[1,2,3],[4,5,6,7]]."""
    s = s.strip()
    if s in ("$", "()", ""):
        return []
    inner = s[1:-1]
    return [[int(v) for v in re.findall(r"\d+", grp)] for grp in re.findall(r"\(([^()]*)\)", inner)]


def _length_unit(entities: dict[int, str]) -> tuple[str, float]:
    for body in entities.values():
        if "LENGTH_UNIT" in body and "SI_UNIT" in body:
            m = re.search(r"SI_UNIT\s*\(\s*(\.\w+\.|\$)\s*,\s*\.(\w+)\.\s*\)", body)
            if m:
                pre = m.group(1).strip(".") if m.group(1) != "$" else None
                base = m.group(2)
                if base != "METRE":
                    raise ValueError(f"unsupported SI length base {base}")
                f = _SI_PREFIX.get(pre)
                if f is None:
                    raise ValueError(f"unsupported SI prefix {pre}")
                return (f"{pre or ''}METRE", f * 1000.0)
        if "LENGTH_UNIT" in body and "CONVERSION_BASED_UNIT" in body:
            m = re.search(r"CONVERSION_BASED_UNIT\s*\(\s*'([^']*)'", body)
            name = (m.group(1) if m else "").upper()
            if name in ("INCH", "IN"):
                return ("INCH", 25.4)
            if name in ("FOOT", "FT"):
                return ("FOOT", 304.8)
            raise ValueError(f"unsupported conversion-based length unit {name!r}")
    raise ValueError("no length unit found in STEP header/context")


def read_tessellated_step(path: str) -> list[TessellatedBody]:
    text = open(path, "r", errors="replace").read()
    data = text[text.index("DATA;") + 5:]
    entities = {int(m.group(1)): " ".join(m.group(2).split()) for m in _ENTITY_RE.finditer(data + "ENDSEC")}
    unit_name, scale = _length_unit(entities)

    coords: dict[int, np.ndarray] = {}
    for eid, body in entities.items():
        if body.startswith("COORDINATES_LIST"):
            a = _args(body)
            pts = _nums(a[2]).reshape(-1, 3)
            if int(a[1]) != len(pts):
                raise ValueError(f"#{eid}: COORDINATES_LIST count mismatch {a[1]} vs {len(pts)}")
            coords[eid] = pts

    shells: dict[str, list[int]] = {}
    for eid, body in entities.items():
        if body.startswith("TESSELLATED_SHELL"):
            a = _args(body)
            name = a[0].strip("'")
            shells[name or f"shell_{eid}"] = [int(x) for x in re.findall(r"#(\d+)", a[1])]

    def face_tris(eid: int, warn: list) -> tuple[int, np.ndarray]:
        body = entities[eid]
        a = _args(body)
        cref = int(a[1].lstrip("#"))
        pnindex = _nested_int_lists("(" + a[5] + ")") if a[5] not in ("$", "()") else []
        pn = np.array(pnindex[0], dtype=int) if pnindex else None
        tris: list[list[int]] = []
        if body.startswith("TRIANGULATED_FACE"):
            tris = _nested_int_lists(a[6])
        elif body.startswith("COMPLEX_TRIANGULATED_FACE"):
            for strip in _nested_int_lists(a[6]):
                for k in range(len(strip) - 2):
                    t = [strip[k], strip[k + 1], strip[k + 2]] if k % 2 == 0 else [strip[k + 1], strip[k], strip[k + 2]]
                    if len(set(t)) == 3:
                        tris.append(t)
            for fan in _nested_int_lists(a[7]):
                for k in range(1, len(fan) - 1):
                    t = [fan[0], fan[k], fan[k + 1]]
                    if len(set(t)) == 3:
                        tris.append(t)
        else:
            raise ValueError(f"#{eid}: unsupported face entity")
        t = np.array(tris, dtype=int) - 1  # STEP indices are 1-based
        if pn is not None:
            t = pn[t] - 1
        return cref, t

    bodies = []
    for name, face_ids in shells.items():
        warn: list = []
        verts, faces, off = [], [], 0
        cache: dict[int, int] = {}
        for fid in face_ids:
            cref, t = face_tris(fid, warn)
            if cref not in cache:
                cache[cref] = off
                verts.append(coords[cref]); off += len(coords[cref])
            faces.append(t + cache[cref])
        V = np.vstack(verts) if verts else np.zeros((0, 3))
        F = np.vstack(faces) if faces else np.zeros((0, 3), int)
        bodies.append(TessellatedBody(name, V, F, unit_name, scale, len(face_ids), warn))
    return bodies


def is_tessellated_step(path: str) -> bool:
    with open(path, "r", errors="replace") as fh:
        txt = fh.read()
    return "TESSELLATED_SHELL" in txt or "TRIANGULATED_FACE" in txt


def has_brep_solid(path: str) -> bool:
    with open(path, "r", errors="replace") as fh:
        txt = fh.read()
    return any(k in txt for k in ("MANIFOLD_SOLID_BREP", "ADVANCED_FACE", "BREP_WITH_VOIDS"))
