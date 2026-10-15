# residual-limb-fea

This toolkit runs full residual-limb finite element analysis for an adjustable
transtibial petal socket. It covers geometry audit, conforming layered
skin/fat/muscle/bone meshing, nonlinear frictional-contact solves in
**FEBio 4.13**, solver-field visualisation (static figures, VTK/ParaView and
an interactive browser viewer) and reproducible LaTeX reports.

> **Scope.** The results are *numerical*. Verification (the equations are
> solved correctly) is kept separate from validation (the model represents a
> real limb), and **no clinical validation has been done**. A converged run
> does not establish grip, comfort or tissue safety. See
> [LIMITATIONS.md](LIMITATIONS.md).

The repository contains **no patient geometry**. The bundled example is a
synthetic limb (`rlfea.geometry.synthetic`). Supplied scans or CAD stay in a
private directory referenced by a private configuration file.

## Pipeline

```
geometry (STEP B-rep | STEP AP242 tessellated | STL/OBJ/PLY)  -- immutable source, SHA-256 recorded
  -> audit (units, bbox, closure, components, orientation, self-intersection, containment, empty bodies)
  -> logged repairs/transforms (before/after audit per operation)
  -> layered conforming mesh: skin + fat wedges (penta6, three-field) offset from the skin,
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

Alternatively, run `docker build -t rlfea .` for a complete pinned environment
(FEBio, Python and TeX).

## End-to-end commands

```bash
rlfea capabilities                                  # FEBio found? every required feature exercised
rlfea verify --out verification_runs                # analytical / independent-reference benchmarks
rlfea run configs/synthetic_human.yaml --runs runs  # synthetic full limb: tightening (75 N example) + pull-out
rlfea status runs/<run_id>                          # prepared | running | converged | failed | interrupted
rlfea post runs/<run_id>                            # re-derive fields, VTK, figures, report from raw output
rlfea solve runs/<run_id> --restart runs/<run_id>/solver/model.dmp [--set etol=0.01]
                                                    # resume from the last converged step (patched FEBio, see docs)
rlfea gui --runs runs --port 8080                   # interactive viewer at http://localhost:8080
rlfea audit path/to/limb.step path/to/bone.step     # geometry audit only (no modification)
pytest -q                                           # unit tests (+ FEBio benchmarks when febio4 is available)
```

### Your own (private) geometry

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

## Results

The current status is recorded in [docs/RESULTS.md](docs/RESULTS.md),
generated from the run directories.

## License

MIT for this code. FEBio is distributed under its own licence (MIT, see the
FEBio repository).
