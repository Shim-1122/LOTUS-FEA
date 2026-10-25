"""Measured (segmented) layer thickness: inner surfaces at known offsets
must be recovered by the layered mesher."""
import numpy as np

from rlfea.geometry import repair as R
from rlfea.geometry.surface import TriSurface
from rlfea.geometry.synthetic import synthetic_limb
from rlfea.mesh.layers import _chunked_distance, build_layered_model


def _offset(s: TriSurface, d: float, name: str) -> TriSurface:
    n = s.vertex_normals(smooth_iters=3)
    return TriSurface(name, s.vertices - d * n, s.faces.copy())


def test_measured_layers_recover_known_thickness():
    b = synthetic_limb(length=200.0)
    skin = R.remesh(R.basic_cleanup(b.pop("skin")), 10.0)
    bones = {k: R.remesh(R.basic_cleanup(v), 8.0) for k, v in b.items()}
    skin_inner = _offset(skin, 2.0, "skin_inner")
    fat_inner = _offset(skin, 5.0, "fat_inner")
    cfg = dict(skin_thickness=1.5, fat_thickness=3.0, min_layer=0.15, normal_smoothing=3,
               thickness_smoothing=5, axis=[0, 0, 1.0], skin_fiber_reference="circumferential",
               muscle_formulation="ut4", muscle_fiber=False, optimize=True, mesh_bones=True,
               proximal_cut=True, max_fraction_total=0.9,
               measured_layers={"skin_inner": skin_inner, "fat_inner": fat_inner})
    model, rep = build_layered_model(skin, bones, cfg)
    assert abs(rep["measured_skin_thickness_mm"]["median"] - 2.0) < 0.2
    assert abs(rep["measured_fat_thickness_mm"]["median"] - 3.0) < 0.3
    assert rep["offset_vs_segmentation_mm"]["skin_inner"]["median"] < 0.2
    assert rep["offset_vs_segmentation_mm"]["fat_inner"]["median"] < 0.3
    assert "MEASURED" in rep["assumption"]
    assert rep["min_wedge_scaled_jacobian"] > 0


def test_chunked_distance_matches_plane():
    V = np.array([[0, 0, 0], [10, 0, 0], [0, 10, 0], [10, 10, 0.0]])
    F = np.array([[0, 1, 2], [1, 3, 2]])
    s = TriSurface("plane", V, F)
    pts = np.array([[2, 2, 3.0], [5, 5, -1.5]])
    assert np.allclose(_chunked_distance(s, pts, chunk=1), [3.0, 1.5])
