# Results

> **Superseded tension.** The runs recorded in this file used a 75 N line tension. That value was a placeholder based on incomplete data, and the team has retired it. Current runs use the 60 lbf (267 N) design tension from [GRIP_ESTIMATE.md](GRIP_ESTIMATE.md). These results are kept as a record.


## LOTUS V5 CAD socket (`configs/lotus_v5_socket.yaml`, run `lotus_v5_tubes_T75`)

The socket comes from the LOTUS V5 (closed) CAD, used tentatively: petals,
silicone pads (about 2.7 mm thick), hinge axles and the cup/hinge-shield
cavity are all taken from the CAD. The CAD file itself is not in this
repository.

* **Fit to the limb:** fitted automatically (`scripts/fit_socket_pose.py`).
  The closed pads sit 0.7 to 2.5 mm off the skin (median 1.0 mm).
* **Starting state:** each petal is turned about 0.9 to 1.0° about its own
  hinge so the pads start just touching the skin. Seven cup-cavity points
  that sat inside the skin were moved out to 0.1 mm clearance.
* **Lace:** RevoFit routing through the hollow tubes in the petals. One
  strand leaves the dial on the anterior petal and runs through the lower
  tubes of all four petals. It then crosses up through the upper tubes
  (anterior, medial, posterior) and passes over the lateral petal, which
  has no upper tube by design, back to the dial. The tubes are modelled as
  frictionless.
* **Unchanged from the earlier knee runs:** the limb, materials, loads and
  solver settings. The distal cushion is still a parametric layer, because
  it is not in the CAD.

**State: tightening to 75 N converged. Pull-out FAILED at t = 1.2 (0.4 mm
of the 2 mm).** Two muscle elements inverted under the top edge of the
lateral pad. Their volume ratio had already fallen to 0.41 at 75 N.

| t | load | max tissue displacement (mm) | peak contact pressure (kPa) | contact force (N) | slip area fraction | max penetration (mm) | min J | force imbalance |
|---|---|---|---|---|---|---|---|---|
| 0.25 | lace 18.8 N | 4.18 | 56 | 141 | 0.409 | 0.47 | 0.706 | 1.7e-7 |
| 0.50 | lace 37.5 N | 7.04 | 93 | 300 | 0.231 | 1.07 | 0.579 | 4.0e-7 |
| 0.75 | lace 56.2 N | 9.74 | 176 | 463 | 0.176 | 3.11 | 0.475 | 7.8e-9 |
| 1.00 | lace 75 N | 12.3 | 196 | 627 | 0.128 | 2.38 | 0.411 | 1.9e-6 |

Contact at 75 N, by petal:

| Petal | Normal force (N) | Peak pressure (kPa) |
|---|---|---|
| anterior | 157 | 101 |
| medial | 176 | 196 |
| posterior | 105 | 58 |
| lateral | 171 | 185 |
| cushion | 18 | 17 |
| cup wall | 0 | 0 |

Each petal's free body balances to within 0.06% of the lace force.

**Findings to act on:**
* **Muscle crushed under the lateral pad's top edge.** About 30 to 60 muscle
  elements under that edge (θ ≈ 0°, 120 to 136 mm above the distal end)
  compress to J = 0.41 to 0.7, and this happens with both lace routings.
  That is a sharp edge loading nearly incompressible tissue. It is worth
  checking whether the real pad has a flared or rounded top lip, and the
  mesh there should be refined.
* **Penetration too large.** Penalty penetration reaches 2 to 3 mm at the
  pad edges, about the same as the pad thickness. Peak pressures are
  therefore underestimated, and contact at the edges is the weakest part of
  this model.
* **Cup wall unused.** The cup wall never touches the skin up to 75 N.

An earlier run (`lotus_v5_T75`) used the same CAD with a guessed single-row
lace. At 75 N it gave 8.7 mm max displacement, 141 kPa peak pressure and
485 N contact force. Its pull-out reached t = 1.987 before hitting the run
time limit.

## Earlier parametric four-petal socket (`configs/knee_limb_segmented.yaml`)

Run `knee_limb_4petal_T75` used FEBio 4.13 with the restart patch. It
**converged** through tightening to the 75 N example lace tension (t = 1)
and a 2 mm pull-out (t = 2), in 1,265 s. These numbers come from a model
and have not been clinically validated. Peak contact pressure depends on the
mesh and is underestimated by penalty contact.

### Model

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

### Results by load step

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

### Compared with earlier runs (all at 75 N)

| Run | Petals | Pad shear modulus | Cushion | Max displacement (mm) | Peak pressure (kPa) | Contact force (N) | Pull-out |
|---|---|---|---|---|---|---|---|
| knee_limb_4petal_T75 (this run) | 4 | 23 kPa (gel liner) | 85A TPU gyroid, E = 71 kPa | 6.49 | 57 | 560 | converged |
| knee_limb_segmented_T75 | 6 | 100 kPa (assumed) | assumed gel, shear 20 kPa | 7.85 | 103 | 591 | converged |
| STEP model, no tibia, assumed layers | 6 | 100 kPa (assumed) | assumed gel, shear 20 kPa | 9.89 | 120 | 623 | did not converge |

The run-to-run differences mix several changes at once (petal count, pad
and cushion stiffness, anatomy), so none of them is attributed to a single
cause. The softer gel pads are the most likely reason peak pressure
roughly halved: they spread the load over a larger area.
