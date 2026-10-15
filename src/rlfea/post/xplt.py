"""Minimal reader for FEBio .xplt plot files (format 0x0035, FEBio 4.x).

Only what the contact post-processing needs is decoded: the surface
dictionary, surface names and per-state FACE_DATA arrays (FMT_ITEM).  FEBio
writes the true contact traction (Lagrange multipliers, `m_tr`) and the stick
flag ONLY to the plot file; the logfile "contact_traction" is a stress-based
surface traction of the underlying solid element (verified in
FEContactSurface::GetSurfaceTraction), so it must not be used as contact
traction.

Layout: uint32 magic 0x00464542, then chunks [uint32 id][uint32 size][payload].
"""
from __future__ import annotations

import struct
import zlib

import numpy as np

MAGIC = 0x00464542
ROOT, HEADER, HDR_COMPRESSION = 0x01000000, 0x01010000, 0x01010004
DICTIONARY, DIC_ITEM, DIC_ITEM_TYPE, DIC_ITEM_FMT, DIC_ITEM_NAME = 0x01020000, 0x01020001, 0x01020002, 0x01020003, 0x01020004
DIC_GLOBAL, DIC_NODAL, DIC_DOMAIN, DIC_SURFACE, DIC_EDGE = 0x01021000, 0x01023000, 0x01024000, 0x01025000, 0x01026000
MESH, SURFACE_SECTION, SURFACE, SURFACE_HDR, SURFACE_ID, SURFACE_FACES, SURFACE_NAME = (
    0x01040000, 0x01043000, 0x01043100, 0x01043101, 0x01043102, 0x01043103, 0x01043104)
STATE, STATE_HEADER, STATE_HDR_TIME = 0x02000000, 0x02010000, 0x02010002
STATE_DATA, STATE_VARIABLE, STATE_VAR_ID, STATE_VAR_DATA, FACE_DATA = 0x02020000, 0x02020001, 0x02020002, 0x02020003, 0x02020500
VAR_SIZE = {0: 1, 1: 3, 2: 6, 3: 1, 4: 9, 5: 3, 6: 9}   # FLOAT, VEC3F, MAT3FS, MAT3FD, TENS4FS, MAT3F, ARRAY(approx)


def _chunks(buf, start, end):
    p = start
    while p + 8 <= end:
        cid, size = struct.unpack_from("<II", buf, p)
        yield cid, p + 8, p + 8 + size
        p += 8 + size


def _str(buf, a, b):
    raw = buf[a:b]
    if len(raw) >= 4:
        n = struct.unpack_from("<i", raw, 0)[0]
        if 0 <= n <= len(raw) - 4:
            return raw[4:4 + n].split(b"\0")[0].decode(errors="replace")
    return raw.split(b"\0")[0].decode(errors="replace")


def read_face_data(path: str, wanted=("contact traction", "contact stick", "contact pressure", "contact gap")) -> dict:
    """Returns {'surfaces': [names], 'times': [...], 'data': {var: [ {surface_name: array} per state ]}}."""
    raw = open(path, "rb").read()
    if struct.unpack_from("<I", raw, 0)[0] != MAGIC:
        raise ValueError("not an FEBio xplt file")
    buf = raw
    out = {"surfaces": [], "times": [], "data": {}}
    surf_dict = []
    compressed = 0
    p = 4
    # root chunk
    for cid, a, b in _chunks(buf, p, len(buf)):
        if cid == ROOT:
            for c2, a2, b2 in _chunks(buf, a, b):
                if c2 == HEADER:
                    for c3, a3, b3 in _chunks(buf, a2, b2):
                        if c3 == HDR_COMPRESSION:
                            compressed = struct.unpack_from("<I", buf, a3)[0]
                elif c2 == DICTIONARY:
                    for c3, a3, b3 in _chunks(buf, a2, b2):
                        if c3 == DIC_SURFACE:
                            for c4, a4, b4 in _chunks(buf, a3, b3):
                                if c4 == DIC_ITEM:
                                    item = {}
                                    for c5, a5, b5 in _chunks(buf, a4, b4):
                                        if c5 == DIC_ITEM_TYPE:
                                            item["type"] = struct.unpack_from("<I", buf, a5)[0]
                                        elif c5 == DIC_ITEM_FMT:
                                            item["fmt"] = struct.unpack_from("<I", buf, a5)[0]
                                        elif c5 == DIC_ITEM_NAME:
                                            item["name"] = _str(buf, a5, b5)
                                    surf_dict.append(item)
                elif c2 == MESH:
                    for c3, a3, b3 in _chunks(buf, a2, b2):
                        if c3 == SURFACE_SECTION:
                            for c4, a4, b4 in _chunks(buf, a3, b3):
                                if c4 == SURFACE:
                                    for c5, a5, b5 in _chunks(buf, a4, b4):
                                        if c5 == SURFACE_HDR:
                                            info = {}
                                            for c6, a6, b6 in _chunks(buf, a5, b5):
                                                if c6 == SURFACE_ID:
                                                    info["id"] = struct.unpack_from("<i", buf, a6)[0]
                                                elif c6 == SURFACE_FACES:
                                                    info["faces"] = struct.unpack_from("<i", buf, a6)[0]
                                                elif c6 == SURFACE_NAME:
                                                    info["name"] = _str(buf, a6, b6)
                                            out["surfaces"].append(info)
            break
    if compressed:
        raise NotImplementedError("compressed xplt not supported (FEBio built without zlib writes uncompressed files)")
    names = [s.get("name", f"surface{i + 1}") for i, s in enumerate(out["surfaces"])]
    for w in wanted:
        out["data"][w] = []
    for cid, a, b in _chunks(buf, 4, len(buf)):
        if cid != STATE:
            continue
        t = None
        rec = {w: {} for w in wanted}
        for c2, a2, b2 in _chunks(buf, a, b):
            if c2 == STATE_HEADER:
                for c3, a3, b3 in _chunks(buf, a2, b2):
                    if c3 == STATE_HDR_TIME:
                        t = struct.unpack_from("<f", buf, a3)[0]
            elif c2 == STATE_DATA:
                for c3, a3, b3 in _chunks(buf, a2, b2):
                    if c3 != FACE_DATA:
                        continue
                    for c4, a4, b4 in _chunks(buf, a3, b3):
                        if c4 != STATE_VARIABLE:
                            continue
                        vid, payload = None, None
                        for c5, a5, b5 in _chunks(buf, a4, b4):
                            if c5 == STATE_VAR_ID:
                                vid = struct.unpack_from("<I", buf, a5)[0]
                            elif c5 == STATE_VAR_DATA:
                                payload = (a5, b5)
                        if vid is None or payload is None or vid - 1 >= len(surf_dict):
                            continue
                        item = surf_dict[vid - 1]
                        if item.get("name") not in wanted:
                            continue
                        nc = VAR_SIZE.get(item.get("type", 0), 1)
                        for c6, a6, b6 in _chunks(buf, *payload):
                            arr = np.frombuffer(buf[a6:b6], dtype="<f4").astype(float)
                            region = c6 - 1
                            nm = names[region] if region < len(names) else f"surface{c6}"
                            rec[item["name"]][nm] = arr.reshape(-1, nc) if nc > 1 else arr
        out["times"].append(t)
        for w in wanted:
            out["data"][w].append(rec[w])
    out["dictionary"] = surf_dict
    return out
