# Limitations

This software produces **numerical** finite element results. Verification
(the equations are solved correctly) is separate from validation (the model
represents the real limb). No validation has been performed.

## Geometry and anatomy
* Supplied residual-limb model (private, not in this repository): the bone
  set is **incomplete**. It contains a distal femur, a patella and one slender
  unidentified bone (about 100 mm x 28 mm). **No tibia is present**, and
  `Limb w- Bone - Part 1 (1).step` is an empty body. Missing bones are never
  created or substituted. Results for that geometry cannot represent a real
  transtibial limb's load path, because the tibia normally carries most of
  the socket load.
* Skin, fat and muscle are **assumed layers** of configurable thickness offset
  from the skin envelope. They are not segmentation. Layer thickness is
  reduced automatically where bone approaches the skin, and every reduction
  is reported.
* Bones are rigid by default, with explicit supports. The deformable-bone
  option uses homogeneous cortical properties for whole bones.
* Tissue-bone interfaces are bonded (shared nodes).

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
* Stick or slip is derived from logged tractions. FEBio's own stick flag is
  available only in the `.xplt` file.
