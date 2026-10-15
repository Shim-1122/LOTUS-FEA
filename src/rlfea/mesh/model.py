"""Solver-neutral finite element model container."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

NODES_PER = {"tet4": 4, "penta6": 6, "hex8": 8, "tri3": 3, "quad4": 4}


@dataclass
class ElementBlock:
    name: str
    etype: str                      # tet4 | penta6 | hex8
    conn: np.ndarray                # (n, k) 0-based global node indices
    material: str                   # key into the material set
    region: str                     # anatomical / component label
    formulation: str = "default"    # ut4 | three-field | standard | rigid
    mat_axis: np.ndarray | None = None  # (n, 2, 3): rows a, d (FEBio vector convention)
    fiber: np.ndarray | None = None     # (n, 3) optional single fibre direction
    meta: dict = field(default_factory=dict)


@dataclass
class Surface:
    name: str
    etype: str                      # tri3 | quad4
    conn: np.ndarray                # 0-based, outward oriented (w.r.t. owning solid)
    meta: dict = field(default_factory=dict)


@dataclass
class FEModel:
    nodes: np.ndarray
    blocks: list[ElementBlock] = field(default_factory=list)
    surfaces: dict[str, Surface] = field(default_factory=dict)
    nodesets: dict[str, np.ndarray] = field(default_factory=dict)
    meta: dict = field(default_factory=dict)

    def block(self, name: str) -> ElementBlock:
        for b in self.blocks:
            if b.name == name:
                return b
        raise KeyError(name)

    def n_elements(self) -> int:
        return int(sum(len(b.conn) for b in self.blocks))

    def add_nodes(self, X: np.ndarray) -> np.ndarray:
        start = len(self.nodes)
        self.nodes = np.vstack([self.nodes, X])
        return np.arange(start, start + len(X))

    def element_offsets(self) -> dict[str, int]:
        off, out = 0, {}
        for b in self.blocks:
            out[b.name] = off
            off += len(b.conn)
        return out

    def to_meshio(self, point_data=None, cell_data=None):
        """meshio Mesh with one cell block per element block (VTK types)."""
        import meshio
        vtk = {"tet4": "tetra", "penta6": "wedge", "hex8": "hexahedron"}
        cells = [(vtk[b.etype], b.conn) for b in self.blocks]
        cd = {"block_id": [np.full(len(b.conn), i, np.int32) for i, b in enumerate(self.blocks)]}
        if cell_data:
            cd.update(cell_data)
        return meshio.Mesh(self.nodes, cells, point_data=point_data or {}, cell_data=cd)

    def summary(self) -> dict:
        return {"n_nodes": int(len(self.nodes)), "n_elements": self.n_elements(),
                "blocks": [{"name": b.name, "etype": b.etype, "n": int(len(b.conn)), "material": b.material,
                            "region": b.region, "formulation": b.formulation} for b in self.blocks],
                "surfaces": {k: {"etype": s.etype, "n": int(len(s.conn))} for k, s in self.surfaces.items()},
                "nodesets": {k: int(len(v)) for k, v in self.nodesets.items()}}
