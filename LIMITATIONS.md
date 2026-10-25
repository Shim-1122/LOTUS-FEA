# Limitations

This software produces **numerical** finite element results. Verification
(the equations are solved correctly) is separate from validation (the model
represents the real limb). No validation has been performed.

## Geometry and anatomy
* The segmented limb in `data/limb/` has four bones: a candidate tibia,
  distal femur, patella and fibula. The labels are heuristic (size, shape,
  position) and unconfirmed. Nobody has verified where the segmentation came
  from or how accurate it is.
* Skin and fat thickness are measured at every skin vertex as the
  closest-point distance to the segmented inner surfaces. Meshed layers
  follow those surfaces closely (median deviation 0.01 mm for inner skin and
  0.08 mm for inner fat). Where offsetting would self-intersect or cut into
  bone, the layers are thinned locally, and each change is logged with its
  deviation; the worst local deviation is reported in the run's
  `mesh/layer_report.json`. Fat that is 0 mm thick in the segmentation is
  raised to 0.15 mm.
* The synthetic limb (`configs/synthetic_human.yaml`) still uses assumed
  layer thicknesses.
* Bones are rigid by default, with explicit supports. In the knee model the
  femur, tibia and fibula are fixed (extended, locked knee) and the patella
  moves with the tissue. The deformable-bone option uses homogeneous
  cortical properties for whole bones.
* Tissue-bone interfaces are bonded (shared nodes). No ligaments, tendons or
  joint capsule are modelled.

## Socket and loads
* No current socket CAD was supplied, and V6 is excluded. The petal socket
  is a **provisional parametric model**: six rigid petals, silicone pads,
  hinges, a distal cup and a cushion. The Anclote suction puck, the plate and
  the Bulldog adapter are treated only as the rigid support of the cup.
* Lace routing is a provisional single-row circumferential lace unless an
  explicit routing is configured. The 75 N strand tension is an **example**.
  No RevoFit rating has been verified, so no rated-percentage cases are run.
* Cable span forces act on rigid petals as fixed-direction resultants. They
  do not follow petal rotation.
* Each petal has a small radial stabiliser spring (default 0.05 N/mm). It is
  numerical stabilisation. Its force and share of the cable moment are
  reported.
* Friction is Coulomb with a screening coefficient (0.5). No sweat, liner or
  adhesion models are included.
* Suction is not modelled unless a justified pressure and a sealing
  assumption are supplied.

## Materials
* Skin HGO parameters come from excised back skin of elderly donors (Ni
  Annaidh 2012), tested in uniaxial tension. Fat and muscle parameters come
  from inverse FE of in-vivo indentation of a single transtibial residuum
  (Sengeh 2016). The fat layer uses that study's combined skin-adipose
  parameters. Silicone, cushion and TPU-infill properties are uncalibrated
  assumptions.
* The behaviour is quasi-static hyperelastic. It does not include
  viscoelasticity, perfusion, active muscle, damage or residual stress.

## Numerical
* Peak stresses and strains near bony prominences and contact edges are
  mesh-dependent. They are not injury predictions.
* A converged solution does not establish grip, comfort or safety. A
  prescribed-traction experiment cannot establish grip.
* Contact traction and stick/slip are FEBio's own contact quantities, read
  from the `.xplt` plot file. If that file is missing, the code falls back to
  the logfile surface traction, which is stress-based rather than contact
  traction, and labels it as such.
* The FEBio build carries a small source patch (`scripts/patches/`) that fixes
  two restart serialization defects in FEBio 4.13. Without it, every `-r`
  restart of a model with rigid bodies or logfile face data crashed. The patch
  does not change constitutive, element, contact or solver code. The restart
  check in `rlfea verify` compares a killed-and-resumed run against an
  uninterrupted one. Dump files from an unpatched build cannot be resumed by
  the patched build.
* Energy convergence tolerance: the synthetic runs use `etol = 1e-3`. The
  knee and earlier supplied-limb runs use FEBio's default `etol = 0.01`, because with `1e-3`
  full Newton entered a friction stick/slip two-cycle at 48 N. That cycle
  kept the energy norm just above the 1e-3 tolerance, although the residual
  and displacement criteria were met. The setting is recorded in every
  run's resolved configuration.
* Penalty contact (scale 0.1 on FEBio's auto-penalty) leaves up to about
  1 mm of penetration at 75 N on the synthetic model, which is on the order
  of the 3 mm pad thickness. In the augmented-Lagrangian sensitivity run at
  26 N, penetration fell from 0.42 to 0.08 mm and peak contact pressure rose
  by 60%. Penalty peak pressures are therefore underestimates. The
  augmented-Lagrangian run did not reach 75 N (it stalled at the minimum
  step), so no penalty-free 75 N result exists.
* Pull-out on the earlier STEP geometry (no tibia, assumed layers) did not converge; on the segmented knee model it converged: Newton iterations
  entered a frictional stick/slip limit cycle at the first increment. The
  synthetic pull-out converged. Frictionless contact did not converge past
  about 49 N on the synthetic model.
* Knee model at 75 N: 15 of 58,600 elements (0.03%) have a volume ratio J
  below 0.8 (lowest 0.68, falling to 0.58 by the end of pull-out). They sit
  in one small patch of muscle, 126 to 154 mm above the distal end. Muscle
  is nearly incompressible, so this is flagged as a possible mesh or
  local-loading artefact rather than reported as a result.
