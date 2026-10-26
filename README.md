# LOTUS FEA

LOTUS FEA (the finite element model for the LOTUS socket: Lighted Overnight Transtibial Universal Socket) runs full residual-limb finite element analysis of a transtibial
(below-knee) limb inside an adjustable petal socket. It was built for a
senior design project on a nighttime prosthesis. It covers:

- geometry audit;
- conforming skin, fat, muscle and bone meshing;
- nonlinear frictional-contact solves in **FEBio 4.13**;
- solver-field visualisation: static figures, VTK/ParaView and an
  interactive browser viewer;
- reproducible LaTeX reports.

> **Scope.** The results are *numerical*. Verification (the equations are
> solved correctly) is kept separate from validation (the model represents a
> real limb), and **no clinical validation has been done**. A converged run
> does not establish grip, comfort or tissue safety. See
> [LIMITATIONS.md](LIMITATIONS.md).

## What's in the repository

* `data/limb/` holds the segmented limb geometry as STL files in mm: skin
  (outer and inner shells), fat, muscle and four bones. The owner chose to
  publish this geometry. See [data/limb/README.md](data/limb/README.md) for
  the inventory, hashes and bone labels.
* `configs/knee_limb_segmented.yaml` is the main model. Skin and fat
  thickness are **measured** at every skin vertex from the segmented inner
  surfaces, and the four bones include a candidate tibia and fibula.
* `configs/lotus_v5_socket.yaml` puts the same limb inside the LOTUS V5
  socket CAD, with the lace through the hollow petal tubes. The STEP file is
  not included; place it at `data/socket/LOTUS_V5_Closed.step`.
* `configs/synthetic_human.yaml` is a fully synthetic limb used for
  verification and for the sensitivity study.
* `site/gait-tester.html` (LOTUS Gait Tester) is a single-file screening page:
  grip vs gait load at a chosen lace tension, built on the FEA results in its
  `RESULTS` object. It runs no FEA solve and is not clinically validated.
* `site/` is a static website (results page, 3D viewer and report) that is
  ready for Vercel (`vercel.json`).

## Pipeline

```
geometry (STEP B-rep | STEP AP242 tessellated | STL/OBJ/PLY)  -- immutable source, SHA-256 recorded
  -> audit (units, bbox, closure, components, orientation, self-intersection, containment, empty bodies)
  -> logged repairs/transforms (before/after audit per operation)
  -> layered conforming mesh: skin + fat wedges (penta6, three-field) offset from the skin
     (thickness measured per vertex from segmented surfaces, or assumed when none are given),
     muscle tetrahedra (UT4) to the bone surfaces, bone tetrahedra (rigid or deformable)
  -> provisional parametric petal socket: rigid petals + silicone pads (hex8), revolute hinges,
     rigid distal cup + cushion; lace routing -> guide forces -> petal force/moment
  -> FEBio input (febio_spec 4.0) -> FEBio 4.13 (PARDISO, adaptive steps, must points)
  -> raw logfile fields -> derived fields (J, principal log/Green strains, principal stress, pressure)
  -> VTU/VTP/PVD + equilibrium/inversion/energy/contact checks + figures + LaTeX/PDF report
```

## Installation

```bash
# 1) FEBio 4.13 from the official source tag (pinned commit) with MKL PARDISO
sudo PREFIX=/opt/febio-4.13 bash scripts/build_febio.sh
export RLFEA_FEBIO=/opt/febio-4.13/bin/febio4            # or put febio4 on PATH

# 2) Python 3.11+ environment (exact versions in requirements.lock.txt)
python3 -m venv .venv && . .venv/bin/activate
pip install -r requirements.lock.txt && pip install -e . --no-deps
# headless rendering needs OSMesa/GL libraries, e.g. on Ubuntu:
# sudo apt-get install libglu1-mesa libgl1 libosmesa6 libxrender1 libxcursor1 libxft2 libxinerama1
```

Alternatively, run `docker build -t lotus-fea .` for a complete pinned environment
(FEBio, Python and TeX).

## End-to-end commands

```bash
rlfea capabilities                                  # FEBio found? every required feature exercised
rlfea verify --out verification_runs                # analytical / independent-reference benchmarks
rlfea run configs/lotus_v5_stance.yaml --runs runs      # V5 socket: lace to 60 lbf, then stance push (grip check)
rlfea run configs/knee_limb_segmented.yaml --runs runs  # parametric socket: lace to 60 lbf + pull-out
rlfea run configs/synthetic_human.yaml --runs runs  # synthetic full limb: same load case
rlfea status runs/<run_id>                          # prepared | running | converged | failed | interrupted
rlfea post runs/<run_id>                            # re-derive fields, VTK, figures, report from raw output
rlfea solve runs/<run_id> --restart runs/<run_id>/solver/model.dmp [--set etol=0.01]
                                                    # resume from the last converged step (patched FEBio, see docs)
rlfea gui --runs runs --port 8080                   # interactive viewer at http://localhost:8080
rlfea audit path/to/limb.step path/to/bone.step     # geometry audit only (no modification)
pytest -q                                           # unit tests (+ FEBio benchmarks when febio4 is available)
```

### Other geometry

To run different geometry, keep private scans out of public repositories.

1. Copy `configs/supplied_limb_template.yaml` to a private folder and fill in
   the file paths, units (required for STL/OBJ/PLY), the bone roles you have
   confirmed, the expected anatomy and the supports.
2. Fill in `docs/PROVENANCE_TEMPLATE.md` next to the data.
3. Run `rlfea run /private/configs/my_limb.yaml --runs /private/runs`.

## Run directory (what is preserved)

`manifest.json` records the input hashes, repository commit, environment,
FEBio version, command and warnings. `status.json` holds the job state.
`config.yaml` and `config.resolved.json` hold the configuration as written and
in internal units. `geometry/` holds the source and prepared audits, the
transform/repair history and the anatomy roles, including any missing
anatomy. `mesh/` holds `model.vtu` with block/region ids and quality, plus
the layer, socket and blocks reports. `solver/` holds `model.feb`, the FEBio
log, stdout, raw logfile data (`out_*.txt`), `model.xplt`, the restart dump,
`convergence.json` and `load_report.json` (the free-body data). `results/`
holds `volume_t*.vtu`, `contact_t*.vtp`, `collection.pvd` and `summary.json`
(checks per output time). `figures/` and `report/` hold the figures and
`report.tex`/`report.pdf`. Partial outputs of failed or interrupted runs are
processed and labelled with the run state. Results are never shared between
runs.

Open `results/collection.pvd` in ParaView for independent inspection
(element fields are cell data).

## Units and frame

All internal values are in mm, N and MPa (energy in N mm). Every configuration
quantity carries an explicit unit. The limb frame puts the long axis along
+z (proximal), +y anterior, and the distal skin apex at z = 0. The transform
from the source frame is logged.

## LOTUS V5 model: assumptions and test parameters

This section collects every assumption the current LOTUS V5 runs rest on,
and the parameters of each test. Values marked *estimate* or
*approximation* should be replaced by measurements when they exist.

### Patient and loads

| Item | Value | Basis |
|---|---|---|
| Patient | 50 years, 175 cm (5'9"), 81.6 kg (180 lb) before amputation | team |
| Body mass after amputation | 78.1 kg (3.5 kg lost with the limb) | Winter segment fractions (estimate) |
| Prosthesis mass | 1.5 kg | **assumed**; weigh the device |
| Walking stance peak (governing load) | 923 N (207 lbf) into the socket | 1.2 × body weight − prosthesis weight (estimate) |
| Sit-to-stand peak | about 470 N | about 0.6 × body weight (estimate) |
| Swing pull-out | about 34 N | 1.5 kg at 6.5 rad/s, 0.30 m from the knee (estimate) |
| Design target | **1.2 kN** axial (1.3 × the stance peak, "a bit of headroom") | team requirement |
| Distal end share | carries load, but **less than 15%** of it | team requirement |
| Petal share | therefore at least 85% (at least 1.02 kN) | follows from the above |

The derivation is in [docs/GRIP_ESTIMATE.md](docs/GRIP_ESTIMATE.md).

### Lace (RevoFit)

* **Design tension: 60 lbf (267 N)** per strand. A friction-only estimate
  gives 35 to 57 lbf for the 1.2 kN target (skin-silicone friction 0.8 to
  0.5). The earlier 75 N value was a placeholder from incomplete data and
  is retired.
* **Rating:** 200 lbf line tension, as stated by the team; about 300 lbf
  reliably in external (third-party) testing, as reported by the team. 60 lbf is 30% of the rating.
* **Routing** (from the team's sketch): the strand leaves the dial on the
  anterior petal, runs through the **lower** tubes of all four petals, then
  crosses into the **upper** tubes (anterior, medial, posterior). The
  lateral petal has no upper tube by design, so the strand passes over its
  outer surface back to the dial. Tubes are modelled as frictionless guide
  points, and span forces act as fixed-direction resultants on the rigid
  petals.

### Socket geometry (LOTUS V5 CAD, used tentatively)

* Source: `LOTUS_V5_Closed (4).step`, all parts in mm. The STEP file is not
  published (`data/socket/` is git-ignored).
* **Four petals:** anterior (carries the dial), medial, posterior, lateral.
  Petals and hinge axles are rigid; each petal rotates about its own hinge
  axis, taken from the axle's principal axis.
* **Pose:** fitted automatically (`scripts/fit_socket_pose.py`). Each petal
  is then closed about its hinge until the 20th percentile of its pad nodes
  is 0.1 mm into the skin, so contact starts snug.
* **Silicone pads** (about 2.7 to 3 mm thick) are backed by **rigid nylon**:
  each pad's back face shares nodes with its rigid petal shell, so the pad
  can only compress through its thickness (and bulge at its free edges); it
  cannot stretch or slide relative to the petal.
* **Distal cup** wall from the CAD cavity (rigid); any cup point that started
  inside the skin is moved out to 0.1 mm clearance.
* **Distal cushion** (not in the CAD): a parametric 6 mm layer under the
  distal 8 mm of skin.
* **Fibular relief (approximate, under test):** a transverse scallop off the
  lateral petal's top edge, traced from the team's swept-profile sketch. It
  lowers the top edge by 16 mm over the posterior 45% of the petal's width,
  then ramps back to full height over the next 25% (cosine ramp). The pad
  follows the new edge. Set in `socket.cad.trims`
  (`configs/lotus_v5_stance_relief.yaml`). The cut side is assumed to be
  posterior, because the fibular head is posterolateral.

### Limb and tissues

* Limb geometry: the team's CAD, designed from MRI-based approximations
  (`data/limb/`). Bone labels (tibia, fibula, femur, patella) are heuristic.
* Skin and fat thickness are measured at each skin vertex from the
  segmented inner surfaces; muscle fills the space to the bones.
* Bones are rigid and fixed (patella free); the proximal cut is an axial
  roller.
* Skin: HGO, Ni Annaidh 2012. Fat: Ogden, Sengeh 2016. Bulk modulus of
  skin and fat: 100 × shear modulus (near-incompressible, approximation).
* **Muscle is incompressible: its volume does not change under
  compression, but it still deforms freely in shape.** Its shear response
  is unchanged (Ogden, Sengeh 2016, shear modulus 5.2 kPa). Volume is held
  constant with FEBio's penalty method: bulk modulus **5.2 MPa** (1000 ×
  shear, Poisson ratio 0.4995), on UT4 nodally integrated tetrahedra, which do not lock at this
  ratio. See *Muscle incompressibility check* below.
* Silicone pad: neo-Hookean, shear modulus 23.3 kPa (median of silicone-gel
  liners, Sanders 2004), near-incompressible (k = 100 × shear).
* Distal cushion: 85A TPU, 2 walls plus 8% gyroid infill, homogenised as a
  compressible neo-Hookean solid, E = 71 kPa (range 36 to 252 kPa),
  Poisson ratio 0.3.

### Contact and solver

* Contact: FEBio `sliding-elastic`, penalty method (scale 0.1 of auto
  penalty), Coulomb friction μ = 0.5 between skin and silicone (screening
  value, not measured).
* Static, quasi-static loading; adaptive time stepping, full Newton
  (`max_ups: 0`), energy tolerance 0.01, up to 10 retries per step.

### Tests

| Config | Step 1 | Step 2 | Mesh (skin / bone size) | Purpose |
|---|---|---|---|---|
| `lotus_v5_stance.yaml` | lace to 60 lbf | push the socket 10 mm onto the limb (displacement control) | 7 / 3.5 mm | grip at the design tension |
| `lotus_v5_stance_relief.yaml` | lace to 60 lbf | same 10 mm push | 9 / 6 mm | the same, with the approximate fibular relief |
| `lotus_v5_stance_25lbf.yaml` | lace to 25 lbf | same 10 mm push | 9 / 6 mm | grip per lbf of lace (scaling) |
| `lotus_v5_socket.yaml` | lace to 60 lbf | 2 mm pull-out | 9 / 6 mm | swing pull-out check |

`scripts/stance_check.py RUN_DIR` reads a stance run and reports, against
the requirements:

* whether the socket carries **1.2 kN** axially before the run ends, and
  after how much push;
* the **distal share** of the axial load at 1.2 kN (pass if below 15%);
* the petal load, slip area fraction, peak contact pressure and minimum
  volume ratio J (tissue inversion check).

Displacement control is used instead of force control so that a socket
that slips still yields a load-displacement curve and the slip load.

### Muscle incompressibility check

A 20 mm cube of the muscle material, on the same UT4 tetrahedra as the limb,
is held at its sides and base and indented 5 mm (25%) by a 10 × 10 mm
indenter, so it can keep its volume only by flowing up around the
indenter. Run `python scripts/muscle_incompressibility_check.py`:

| Bulk modulus / shear | Bulk modulus (MPa) | Poisson ratio | Volume change at full indent |
|---|---|---|---|
| 100 (previous setting, Sengeh 2016) | 0.52 | 0.4950 | −1.71% |
| **1000 (current)** | **5.2** | **0.4995** | **−0.42%** |
| 10 000 | 52 | 0.5000 | −0.36% |

Going from 1000 to 10 000 hardly changes the result, so the remaining
0.4% is mesh discretisation rather than material compressibility. A
higher ratio would only make the solver stiffer. The shape response is
unchanged: the shear parameters are the same in every row.

Under the lateral pad's top edge, the earlier 75 N run (bulk = 100 × shear)
showed muscle elements compressed to 41% of their volume. With the
incompressible setting that tissue has to move aside instead of
compacting. The stance runs are being repeated with this setting; runs
made before the change are kept with the suffix `_k100`.

## Materials (Scenario A: human tissue; Scenario B: prototype)

Both sets are in `configs/materials/`. Each lists every parameter with its
value, unit, source (DOI), convention conversion, calibration range and
uncertainty. The FEBio strain-energy forms and their mapping to the sources
and to ABAQUS are documented in
[docs/solver_selection.md](docs/solver_selection.md). In short, skin is HGO
(Ni Annaidh 2012), fat and muscle are 2-term Ogden (Sengeh et al. 2016,
in-vivo transtibial indentation), and silicone and the cushion are
uncalibrated assumptions. TPU uses the Reppel & Weinberg 2018 Mooney-Rivlin
fit, which is not inferred from Shore hardness and must be calibrated for
your infill.

## Loads and contact

* The lace is an ordered list of guide points on bodies. Free spans pull both
  ends with the span tension. Optional capstan guide friction can be applied.
  Per-petal resultant force and moment go to the rigid petals. Sum checks are
  run (the lace is an internal load, so the sums must be zero).
* `cable.strand_tension` is an **example** value unless `tension_basis` says
  otherwise. `tension_basis: design` (used by the LOTUS V5 configs) requires
  `cable.design_source`. Rated 25/50/75/100% cases require
  `cable.rated_tension: {value, unit, source, meaning: strand_tension}`
  and `tension_basis: rated`.
* Contact is `sliding-elastic`: finite sliding, separation, compression only
  and Coulomb friction. Bonded interfaces only happen with `contact.bonded: true`.
* `load_model: prescribed_traction` is a labelled experiment. It applies
  pressures from a rigid-petal moment balance and cannot establish grip.
* After tightening, at most one of the following can be run: pull-out (cup
  displaced distally, with the reaction recorded), body weight (proximal
  push), bending or torsion (prescribed cup rotation). Suction requires a
  justified pressure and a sealing assumption.

## Verification

`rlfea verify` runs 47 checks against analytical or independently computed
references. Each check exercises the same FEBio adapter that is used for the
limb. The checks cover:

* HGO, Ogden and Mooney-Rivlin stresses vs an independent strain-energy
  implementation, in tension, compression, shear, volumetric and rotated
  stretch, plus rigid rotation.
* Coulomb sliding, contact opening, the rigid hinge, cable force and moment
  balance, and the energy balance.
* Restart after SIGKILL.
* Near-incompressible locking (indentation at k = 100 mu): standard tet4 is
  expected to lock and does; ut4 and tet10 must be within 10% of the
  three-field hex8 result at the same mesh density.

`pytest` runs the unit tests, plus the FEBio benchmarks when `febio4` is
available.

## Results

The current status is recorded in [docs/RESULTS.md](docs/RESULTS.md) (the full synthetic-run report is [docs/synthetic_report.pdf](docs/synthetic_report.pdf)),
generated from the run directories.

## License

MIT for this code. FEBio is distributed under its own licence (MIT, see the
FEBio repository).
