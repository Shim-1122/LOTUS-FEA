# Segmented limb geometry

These are surface meshes (binary STL, **millimetres**) of one below-knee
residual limb with the knee, exported from a CAD model named
"Limb w- Bone". The team built the model from approximations of an MRI
scan, so it is MRI-informed CAD rather than a direct segmentation of image
voxels. The owner chose to publish it with this repository.
`SHA256SUMS` pins every file, and the run manifest records the hashes it
used.

| File | Shells (largest first) | Size (mm) | Used as |
|---|---|---|---|
| `skin.stl` | 2: outer skin (3.15e6 mm³ enclosed), inner skin (2.70e6 mm³) | 143.7 × 152.1 × 289.1 | outer = model envelope; inner = skin thickness |
| `fat.stl` | 1: the region between inner skin and outer muscle (3.83e5 mm³) | 136.5 × 144.5 × 274.6 | reference only (it is fully defined by the other two surfaces) |
| `muscle.stl` | 3: outer muscle surface (2.00e6 mm³), cavities for part_1 and part_3 | 129.7 × 137.3 × 260.9 | outer = fat thickness |
| `part_1.stl` | 1 (2.75e5 mm³) | 88.1 × 62.5 × 150.0 | bone, **candidate tibia** |
| `part_2.stl` | 1 (3.00e5 mm³) | 95.2 × 74.1 × 114.8 | bone, **candidate distal femur** |
| `part_3.stl` | 1 (2.55e4 mm³) | 47.7 × 23.7 × 47.5 | bone, **candidate patella** |
| `part_4.stl` | 1 (1.99e4 mm³) | 28.3 × 26.2 × 100.4 | bone, **candidate fibula** |

Measured from these files:

* Skin thickness ranges from 1.2 to 9.9 mm (median 3.8 mm).
* Fat thickness ranges from 0 to 11.7 mm (median 3.2 mm). Where it is 0, the
  inner skin touches the muscle; the model raises it to 0.15 mm there.
* Every surface is closed (watertight).
* Two cavities in `muscle.stl` coincide exactly with part_1 and part_3.
  part_2 and part_4 are not cut out of the muscle; the mesher carves them.
* part_2 extends up to 6.8 mm outside the muscle surface, near the proximal
  cut, and comes within 1.5 mm of the skin. The model thins the skin and fat
  layers there and reports each change.

**The bone labels are heuristic.** They come from size, shape and position
(for example, part_4 is slender and lies posterolateral to part_1), and
nobody has confirmed them. Because the layers and bones were approximated from MRI, their accuracy
against the scan has not been quantified, so measured does not mean
validated.
