"""Helpers of the CAD socket builder."""
import numpy as np
import trimesh

from rlfea.geometry.cad_socket import _principal_axis, _ray_hits, _rotate_about


def test_rotation_about_hinge_axis():
    P = np.array([[0.0, 10.0, 5.0]])
    Q = _rotate_about(P, np.array([0.0, 10.0, 0.0]), np.array([1.0, 0.0, 0.0]), np.pi / 2)
    assert np.allclose(Q, [[0.0, 5.0, 0.0]], atol=1e-12)          # 5 mm lever swung 90 deg about x


def test_rays_from_axis_hit_inner_and_outer_face_of_a_tube_segment():
    tube = trimesh.creation.annulus(r_min=40.0, r_max=43.0, height=50.0, sections=128)
    hits = _ray_hits(tube, 0.0, np.deg2rad([0.0, 90.0, 200.0]))
    for h in hits:
        assert len(h) == 2 and abs(h[0] - 40.0) < 0.1 and abs(h[1] - 43.0) < 0.1


def test_principal_axis_of_an_axle():
    ax = trimesh.creation.cylinder(radius=3.0, height=40.0, sections=32)
    ax.apply_transform(trimesh.transformations.rotation_matrix(np.pi / 2, [0, 1, 0]))
    c, d = _principal_axis(np.asarray(ax.vertices))
    assert np.allclose(np.abs(d), [1, 0, 0], atol=1e-6) and np.allclose(c, 0, atol=1e-9)
