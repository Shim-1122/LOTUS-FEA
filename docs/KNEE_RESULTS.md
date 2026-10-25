# Knee model results (`configs/knee_limb_segmented.yaml`)

Run `knee_limb_4petal_T75` used FEBio 4.13 with the restart patch. It
**converged** through tightening to the 75 N example lace tension (t = 1)
and a 2 mm pull-out (t = 2), in 1,265 s. These numbers come from a model
and have not been clinically validated. Peak contact pressure depends on the
mesh and is underestimated by penalty contact.

## Model

The mesh has 58,612 elements:

* **Geometry:** CAD of the limb built from approximations of an MRI scan.
* **Skin and fat:** wedges with **measured** thickness. Skin is 1.2 to 9.9 mm
  thick (median 3.8 mm) and fat is 0 to 11.7 mm (median 3.2 mm). The median
  deviation between the meshed layers and the segmentation is 0.01 mm for
  inner skin and 0.08 mm for inner fat.
* **Muscle:** UT4 tetrahedra.
* **Bones:** four rigid bones (candidate tibia, femur, patella and fibula).
  The femur, tibia and fibula are fixed; the patella moves with the tissue.
* **Socket:** **four** rigid petals (anterior, posterior, medial and
  lateral) with silicone pads, a rigid distal cup and a TPU cushion.

| Material | Model | Values | Basis |
|---|---|---|---|
| Skin | HGO | c = 201 kPa, k1 = 24.5 MPa, k2 = 0.133, kappa = 0.140, gamma = 41° | Ni Annaidh 2012 (human skin); bulk modulus = 100 × shear (approximation) |
| Fat, muscle | 2-term Ogden | c = 5.2 kPa, m = 4.8 | Sengeh 2016 (in-vivo transtibial residuum) |
| Silicone pad | neo-Hookean | shear modulus 23 kPa | median of five silicone-gel liners, Sanders 2004 (range 19 to 53 kPa) |
| Distal cushion | compressible neo-Hookean | E = 71 kPa, Poisson ratio 0.3 | 85A TPU, 2 walls + 8% gyroid: Gibson-Ashby E·(0.08)² with E = 11.1 MPa (NinjaFlex 85A, Reppel & Weinberg 2018); range 36 to 252 kPa |

## Results by load step

| t | load | max tissue displacement (mm) | peak contact pressure (kPa) | contact force (N) | slip area fraction | max penetration (mm) | min J | cup axial force (N) | force imbalance |
|---|---|---|---|---|---|---|---|---|---|
| 0.25 | lace 18.8 N | 1.84 | 24 | 157 | 0.185 | 0.26 | 0.905 | 23.6 | 7.8e-6 |
| 0.50 | lace 37.5 N | 3.52 | 37 | 291 | 0.099 | 0.46 | 0.820 | 42.9 | 1.6e-10 |
| 0.75 | lace 56.2 N | 5.08 | 49 | 427 | 0.075 | 0.71 | 0.737 | 62.1 | 1.9e-6 |
| 1.00 | lace 75 N | 6.49 | 57 | 560 | 0.057 | 0.96 | 0.663 | 81.2 | 7.0e-7 |
| 1.50 | + pull-out 1.0 mm | 5.99 | 58 | 539 | 0.064 | 0.98 | 0.648 | 52.2 | 1.6e-6 |
| 2.00 | + pull-out 2.0 mm | 5.46 | 55 | 511 | 0.081 | 0.98 | 0.637 | 23.5 | 1.9e-9 |

Each petal's free body (lace, hinge, stabiliser and pad contact) balances to
within 0.9% of the lace force. Over the 2 mm pull-out, the cup's axial
constraint force fell from 81 N to 24 N while only 8% of the contact area
slipped. That is consistent with the petals holding the limb by friction in
this model, but it is not proof of grip.

At 75 N, 9 muscle elements (0.02%) have a volume ratio J below 0.8 (lowest
0.66), in one small patch 133 to 154 mm above the distal end. Muscle is
nearly incompressible, so this is flagged as a possible mesh or
local-loading artefact.

## Compared with earlier runs (all at 75 N)

| Run | Petals | Pad shear modulus | Cushion | Max displacement (mm) | Peak pressure (kPa) | Contact force (N) | Pull-out |
|---|---|---|---|---|---|---|---|
| knee_limb_4petal_T75 (this run) | 4 | 23 kPa (gel liner) | 85A TPU gyroid, E = 71 kPa | 6.49 | 57 | 560 | converged |
| knee_limb_segmented_T75 | 6 | 100 kPa (assumed) | assumed gel, shear 20 kPa | 7.85 | 103 | 591 | converged |
| STEP model, no tibia, assumed layers | 6 | 100 kPa (assumed) | assumed gel, shear 20 kPa | 9.89 | 120 | 623 | did not converge |

The run-to-run differences mix several changes at once (petal count, pad
and cushion stiffness, anatomy), so none of them is attributed to a single
cause. The softer gel pads are the most likely reason peak pressure
roughly halved: they spread the load over a larger area.
