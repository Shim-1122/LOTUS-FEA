# Knee model results (`configs/knee_limb_segmented.yaml`)

Run `knee_limb_segmented_T75` used FEBio 4.13 with the restart patch. Its
state is **converged**: tightening to the 75 N example lace tension
(t = 1), then a 2 mm pull-out (t = 2), took 1,382 s. These numbers come
from a model and have not been clinically validated. Peak contact pressure
depends on the mesh and is underestimated by penalty contact.

## Mesh

The mesh has 58,600 elements in the tissue, bones and socket:

* skin and fat as wedges with **measured** thickness;
* UT4 muscle tetrahedra;
* four rigid bones (candidate tibia, femur, patella and fibula);
* six rigid petals with silicone pads, a distal cup and a cushion.

Measured skin thickness ranges from 1.2 to 9.9 mm (median 3.8 mm). Measured
fat thickness ranges from 0 to 11.7 mm (median 3.2 mm). The median deviation
between the meshed layers and the segmented surfaces is 0.01 mm for inner
skin and 0.08 mm for inner fat.

## Results by load step

| t | load | max tissue displacement (mm) | peak contact pressure (kPa) | contact force (N) | slip area fraction | max penetration (mm) | min J | cup axial force (N) | force imbalance |
|---|---|---|---|---|---|---|---|---|---|
| 0.25 | lace 18.8 N | 2.06 | 22 | 165 | 0.211 | 0.35 | 0.884 | 23.8 | 1.9e-6 |
| 0.50 | lace 37.5 N | 4.14 | 36 | 306 | 0.115 | 0.66 | 0.794 | 43.3 | 1.4e-4 |
| 0.75 | lace 56.2 N | 6.08 | 59 | 448 | 0.082 | 0.88 | 0.733 | 63.1 | 2.9e-9 |
| 1.00 | lace 75 N | 7.85 | 103 | 591 | 0.072 | 1.58 | 0.683 | 83.3 | 1.2e-6 |
| 1.50 | + pull-out 1.0 mm | 8.14 | 106 | 566 | 0.077 | 1.61 | 0.647 | 55.3 | 2.2e-5 |
| 2.00 | + pull-out 2.0 mm | 8.42 | 113 | 536 | 0.095 | 1.47 | 0.576 | 26.7 | 2.6e-4 |

Every petal's free body (lace, hinge, stabiliser and pad contact) balances
to within 1% of the lace force. Over the 2 mm pull-out, the cup's axial
constraint force fell from 83 N to 27 N while only 9.5% of the contact area
slipped. That is consistent with the petals holding the limb by friction in
this model, but it is not proof of grip.

## Compared with the earlier STEP model

The earlier model had no tibia, a 1.5 mm skin and a 3 mm fat layer. At
75 N, the segmented model gives 7.9 mm max displacement (vs 9.9 mm), 103 kPa
peak pressure (vs 120 kPa) and 591 N contact force (vs 623 N). Pull-out now
converges.
