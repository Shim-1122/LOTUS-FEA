# Solver backend selection

## Decision
**FEBio 4.13** (febio_spec 4.0), built from the official source tag `v4.13`
(commit `32ae206ff4881dfb54f62296cd1558e58ed9fcc6`) with Intel MKL PARDISO, is the
production backend. Optional ABAQUS input generation is provided as a separate
adapter (`rlfea.solvers.abaqus`) for users with a licence; this repository has
**never executed ABAQUS**.

## Evidence used (official documentation and source of the exact version built)
Capabilities were verified against the FEBio 4.13 User's Manual shipped in the
source tree (`Documentation/FEBio_User_Manual.pdf`, "User's Manual Version 4.13")
and the FEBioMech/FEBioXML source of the same commit, then exercised by
`rlfea capabilities` and `rlfea verify`.

| Requirement | FEBio 4.13 feature | Where verified |
|---|---|---|
| HGO skin with explicit fibres and dispersion | `Holzapfel-Gasser-Ogden` (uncoupled), params `c,k1,k2,kappa,gamma,k`; per-element axes via `<ElementData type="mat_axis">` | UM §4.1.2.8; `FEHolzapfelGasserOgden.cpp`; `FEBioMeshDataSection4.cpp` |
| Fat / muscle hyperelasticity | `Ogden` (uncoupled), `Mooney-Rivlin` (uncoupled) | UM §4.1.2; source |
| Optional passive muscle anisotropy | `trans iso Mooney-Rivlin` (available; disabled until calibrated parameters exist) | UM §4.1.2 |
| Near incompressibility without locking | three-field (Simo-Taylor) element for hex8/penta6 with uncoupled materials; `ut4-solid` nodally averaged tet4 (Gee et al.) | UM §1, §3.7.1, §4.1.2 |
| Finite-sliding contact, separation, compression-only, friction | `sliding-elastic` (facet-on-facet, `fric_coeff`, `tension=0`) | UM §3.14 |
| Explicit bonded interfaces | `tied-elastic` (only when `contact.bonded: true`) | UM §3.14 |
| Rigid bones / petals, hinges | `rigid body`, `rigid revolute joint`, `rigid spring`, `rigid_force`, `rigid_moment` | UM §3.12 |
| Incremental loading with adaptive steps | auto time stepper (`dtmin`, `dtmax` must-point curve, `max_retries`, `opt_iter`) | UM §3.3.2 |
| Restart | `-dump` / `-r` command-line options | `FEBioApp.cpp` |

## Conventions that differ between solvers (do NOT assume identical HGO definitions)
* FEBio HGO ground matrix `W = c/2 (I1~-3)`: FEBio `c` = ABAQUS `2*C10` = Ni Annaidh `mu`.
* FEBio HGO volumetric term (pressure_model 3, set by the HGO class): `U = k/2((J^2-1)/2 - ln J)` = ABAQUS GOH with `D = 2/k`.
* FEBio Ogden: `W = sum c_i/m_i^2 (sum_j l~_j^m_i - 3)`; small-strain shear modulus `mu0 = sum c_i / 2`. ABAQUS Ogden uses `2 mu_i/alpha_i^2`.
* Default FEBio uncoupled volumetric term: `U = k/2 (ln J)^2`.
* Fibres: two families at `+/-gamma` from local `e1` in the `e1-e2` plane; active only when `E_a >= 0`.
These were confirmed numerically: FEBio Cauchy stresses match an independent
implementation of the published equations to ~1e-12 for tension, compression,
shear, volumetric and general deformations (`rlfea verify`).

## Build notes
* MKL 2025+ removed the legacy sparse BLAS routine `mkl_dcsrtrsv` used by FEBio 4.13's ILU0 preconditioner; MKL 2024.2.2 is pinned.
* FEBio's CMake looks only for static MKL on Linux; `scripts/build_febio.sh` applies a one-line documented patch to accept shared libraries.
* An empty `<Globals><Constants/></Globals>` block makes FEBio 4.13 drop the Material section; the writer omits Globals.
* Load-controller ids must be consecutive; the writer renumbers them.

## Not supported / not used
* No constitutive model is substituted silently: if a requested model is unavailable the configuration is rejected.
* Contact stick status is written by FEBio only to the binary `.xplt` plot file. Reports use a stick/slip classification **derived** from the logged contact tractions (`|t_t| >= 0.98 mu p` = slip) and label it as derived. The `.xplt` is preserved for inspection in FEBio Studio.
