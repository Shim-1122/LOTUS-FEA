# Limitations

This software produces **numerical** finite element results. Verification
(the equations are solved correctly) is separate from validation (the model
represents the real limb). No validation has been performed.

## Geometry and anatomy
* The segmented limb in `data/limb/` is CAD built from approximations of an
  MRI scan. It has four bones: a candidate tibia, distal femur, patella and
  fibula. The labels are heuristic (size, shape, position) and unconfirmed,
  and nobody has measured how closely the layers and bones match the scan.
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
* With `configs/lotus_v5_socket.yaml`, the petals, pads, hinge axles and
  cup cavity come from the LOTUS V5 CAD (used tentatively). The petals are
  rotated about 1° about their hinges so the pads start in contact, cup
  points inside the skin are moved out, and the hinge shield is treated as
  rigid. The cushion is not in the CAD and stays parametric.
* Otherwise (`configs/knee_limb_segmented.yaml`), the petal socket
  is a **provisional parametric model**: four rigid petals (anterior, posterior, medial, lateral) with silicone pads,
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
* Skin uses an HGO fit to excised back skin from elderly donors (Ni Annaidh
  2012, uniaxial tension). Its bulk modulus is approximated as 100 x shear,
  the same ratio used for fat and muscle.
* Fat and muscle use Ogden fits from inverse FE of in-vivo indentation of a
  single transtibial residuum (Sengeh 2016). The fat layer uses that study's
  combined skin-adipose parameters.
* The silicone pads stand in for medical-grade gel-liner silicone. They use
  a neo-Hookean model with shear modulus 23 kPa, the median of five
  silicone-gel liner products in Sanders 2004 (range 19 to 53 kPa). This is
  a literature approximation, not a test of the team's silicone. It also
  misses the strong stiffening gels show in compression, so pad compression
  under high pressure is overestimated.
* The distal cushion is 85A TPU printed with 2 walls and 8% gyroid infill.
  It is modelled as one homogenised compressible solid: E = 71 kPa,
  Poisson ratio 0.3. That comes from Gibson-Ashby scaling, E* = Es (0.08)^2,
  with Es = 11.1 MPa from the NinjaFlex 85A fit; the plausible range is 36
  to 252 kPa. The model does not capture the lattice's plateau and
  densification at large compression, and the walls are not resolved.
* All materials are quasi-static hyperelastic: no viscoelasticity,
  perfusion, active muscle, damage or residual stress.

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
* Knee model at 75 N: 9 of 58,612 elements (0.02%) have a volume ratio J
  below 0.8 (lowest 0.66, then 0.64 at the end of pull-out). They sit in one
  small patch of muscle, 133 to 154 mm above the distal end. Muscle is nearly
  incompressible, so this is flagged as a possible mesh or local-loading
  artefact rather than reported as a result.
