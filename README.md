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
* `cable.strand_tension` is an **example** value. Rated 25/50/75/100% cases
  require `cable.rated_tension: {value, unit, source, meaning: strand_tension}`
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
