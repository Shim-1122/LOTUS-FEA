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
  grip vs gait load at a chosen lace tension, with the FEA deformation,
  contact-pressure, shear and stick/slip fields of the stance runs replayed
  on the limb (`scripts/export_gait_fields.py` inlines them). It runs no FEA
  solve and is not clinically validated.
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
  points.
* **Force model: perfectly tangential** (`cable.force_model: tangential`,
  team 2026-10-26). At each guide point only the horizontal circumferential
  (hoop) component of the lace tension is applied; the axial and radial
  chord components are dropped. Each petal therefore receives the hoop
  resultant, which points inward. At 25 lbf: anterior 159 N, medial 298 N,
  posterior 281 N, lateral 193 N (about 930 N in total).

### Socket geometry (LOTUS V5 CAD, used tentatively)

* Source: `LOTUS_V5_Closed (4).step`, all parts in mm. The STEP file is not
  published (`data/socket/` is git-ignored).
* **Four petals:** anterior (carries the dial), medial, posterior, lateral.
  Petals are **rigid solids on frictionless hinges** (axis from each axle's
  principal axis) that **rotate at most 2° inward** (team, 2026-10-26),
  measured from the CAD closed pose (`socket.cad.rotation_limit`). The limit
  is a hinge stop: a rigid 30 mm tab below each hinge (part of the petal)
  meets a stiff bumper fixed to the cup when the petal has turned 2° inward.
  Without the stop, the earlier runs turned the medial petal 11 to 14° and
  the lateral petal 2.5 to 5° inward at only 25 to 30 lbf, driving their top
  edges into the limb.
* **Cup:** the CAD cup wall, hinge shield and puck seat form one fixed rigid
  body. (Runs before 2026-10-26 left the cup wall unconstrained; it drifted
  away and carried no load. Those runs are kept with the suffix `_freecup`.)
* **Pose:** fitted automatically (`scripts/fit_socket_pose.py`). Each petal
  is then closed about its hinge until the 20th percentile of its pad nodes
  is 0.1 mm into the skin, so contact starts snug. That closing (0.88 to
  0.99°) counts toward the 2° limit, leaving 1.0 to 1.1° for the lace.
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
* Skin: HGO, Ni Annaidh 2012. Bulk modulus of skin and fat: 100 × shear
  modulus (near-incompressible, approximation).
* **Tissue set v2 (2026-10-26): Sengeh 2016 is no longer used** (team
  request). Fat: Ogden, shear modulus **6.55 kPa**, exponent 5 (Mukhina et
  al. 2020). Muscle: Ogden, shear modulus **9.5 kPa** (middle of the
  7.4 to 11.7 kPa shear-wave-elastography range measured in transtibial
  residual limbs, Ranger et al. 2022), exponent 5 (Mukhina 2020). Muscle
  uncertainty range 7.4 to 16.2 kPa (relaxed to contracted, Steer et al.
  2020). These are 1.3× (fat) and 1.8× (muscle) stiffer than before.
* **Muscle is incompressible: its volume does not change under
  compression, but it still deforms freely in shape.** Volume is held
  constant with FEBio's penalty method: bulk modulus **9.5 MPa** (1000 ×
  shear, Poisson ratio 0.4995), on UT4 nodally integrated tetrahedra, which
  do not lock at this ratio. See *Muscle incompressibility check* below.
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

A 20 mm cube of the muscle material (computed with the earlier shear values, c = 5.2 kPa; the ratio result carries over), on the same UT4 tetrahedra as the limb,
is held at its sides and base and indented 5 mm (25%) by a 10 × 10 mm
indenter, so it can keep its volume only by flowing up around the
indenter. Run `python scripts/muscle_incompressibility_check.py`:

| Bulk modulus / shear | Bulk modulus (MPa) | Poisson ratio | Volume change at full indent |
|---|---|---|---|
| 100 (previous setting) | 0.52 | 0.4950 | −1.71% |
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

## Mathematical model

Every calculation the LOTUS V5 runs and the Gait Tester rely on, with the
worked numbers for the team patient. Units: mm, N, MPa (kPa where noted);
1 lbf = 4.4482 N; g = 9.81 m/s².

### 1. Patient mass and standing weight

Mass lost with a transtibial amputation, from Winter's segment mass
fractions (foot 1.45%, shank 4.65% of body mass) with a fraction $f$ of the
shank removed:

$$m_{lost} = M\,(0.0145 + 0.0465\,f)$$

With $M$ = 81.6 kg (180 lb) and $f$ = 0.60 (residual tibia about 150 mm):
$m_{lost}$ = 81.6 × 0.0424 = **3.46 kg**, so the body mass after amputation
is **78.1 kg**.

Standing weight with the prosthesis ($m_p$ = 1.5 kg, assumed):

$$W = (M - m_{lost} + m_p)\,g = (78.14 + 1.5) \times 9.81 = \mathbf{781\ N}$$

### 2. Loads on the socket

**Walking stance (governing load).** The vertical ground reaction force
peaks at about 1.2 W at loading response and push-off. The prosthesis's own
weight does not pass through the socket-limb interface:

$$F_{stance} = 1.2\,W - m_p g = 1.2 \times 781.3 - 14.7 = \mathbf{923\ N}\ (207\ \text{lbf})$$

**Swing pull-out.** The prosthesis hangs from the limb and swings about the
knee (peak knee angular velocity $\omega$ ≈ 6.5 rad/s, prosthesis centre of
mass $r$ ≈ 0.30 m from the knee):

$$F_{swing} = m_p\,(g + \omega^2 r) = 1.5\,(9.81 + 6.5^2 \times 0.30) = \mathbf{34\ N}$$

**Sit-to-stand.** Seat-off peak about 0.6 W = **469 N**.

**Design target.** The team asked for "a bit of headroom" over the largest
load, taken as 1.3×:

$$F_{target} = 1.3 \times 923 = \mathbf{1200\ N}$$

The distal end may carry at most 15% of it, so the petals must carry at
least 0.85 × 1200 = **1020 N**.

**Gait curve (Gait Tester).** Load as a function of cycle position $\phi$:

$$F(\phi) = \begin{cases} \max\left(0,\ a(\phi)\,W - m_p g\right) & 0 \le \phi \le 60\% \ (\text{stance}) \\ -s(\phi)\,F_{swing} & 60\% < \phi \le 100\% \ (\text{swing}) \end{cases}$$

$a(\phi)$ is a double hump peaking at 1.2 (12% and 48%) with a midstance
valley near 0.78; $s(\phi)$ rises to 1 at 88%. Both are editable arrays
(`GAIT_CURVE`, `SWING_CURVE`) and are estimates.

### 3. Lace tension needed (friction-only estimate)

The FEA tube-lace run gives the total pad normal force as close to
proportional to the line tension $T$: $N = k_N T$ with $k_N$ = 7.3, 7.8,
8.0 and 8.1 at 25, 50, 75 and 100% of the tension. Using $k_N$ = 8.1, the
petals can resist sliding up to

$$F_{petals} = \mu\,k_N\,T$$

Setting $F_{petals}$ = 1020 N:

$$T = \frac{0.85\,F_{target}}{\mu\,k_N}$$

| μ | T for 1× (784 N on the petals) | **T for 1.3× (1020 N)** | T for 2× (1569 N) |
|---|---|---|---|
| 0.5 | 194 N (44 lbf) | **252 N (57 lbf)** | 387 N (87 lbf) |
| 0.8 | 121 N (27 lbf) | 157 N (35 lbf) | 242 N (54 lbf) |

Swing needs $T = F_{swing} / (\mu k_N)$ = 5 to 8 N (1 to 2 lbf), so it never
governs. The design tension, **60 lbf (267 N)**, is 30% of the 200 lbf
RevoFit rating.

The stance capacity in the friction model is $F_{petals}/(1 - s_d)$ with
$s_d$ = 0.15, and the pull-out capacity is $F_{petals}$ alone (the distal
end cannot resist pull-off). Margin:

$$\text{margin}(\phi) = \frac{\text{capacity}}{|F(\phi)|},\qquad \text{slip if margin} < 1,\quad \text{target} \ge 1.3$$

The minimum safe tension is the smallest $T$ (searched 1 to 300 lbf in
0.5 lbf steps) with margin ≥ 1 at every point of the cycle: **44 lbf** at
μ = 0.5 for the team patient, 57 lbf with the 1.3× headroom.

### 4. Lace forces on the petals

The lace is an ordered list of guide points $\mathbf{x}_k$ on the petals
(the tube entries and exits). Each free span between consecutive points
carries the line tension $T$ and pulls both of its ends toward each other.
In the general ("span") model the force on the body that owns point $k$ is

$$\mathbf{F}_k = T\left(\hat{\mathbf{e}}_{k,k-1} + \hat{\mathbf{e}}_{k,k+1}\right),\qquad \hat{\mathbf{e}}_{k,j} = \frac{\mathbf{x}_j - \mathbf{x}_k}{|\mathbf{x}_j - \mathbf{x}_k|}$$

**Tangential model (used for LOTUS V5).** Only the hoop component is kept,
with $\hat{\mathbf{t}}_k = (-y_k, x_k, 0)/\sqrt{x_k^2 + y_k^2}$:

$$\mathbf{F}_k^{t} = (\mathbf{F}_k \cdot \hat{\mathbf{t}}_k)\,\hat{\mathbf{t}}_k$$

For a petal whose two lace points sit at angles $\pm\varphi$ from its
centre line, the two tangential pulls add to the classic band (capstan-free
hoop) result $F_r = 2T\sin\varphi$, directed inward.

Each rigid petal receives the resultant $\sum \mathbf{F}_k$ and moment. In
the span model the lace is an internal load, so the sum over all bodies is
checked to be zero. In the tangential model the dropped components make the
sum non-zero; that difference is carried by the joints.

**Petal kinematics.** Each petal is a rigid body on a frictionless revolute
joint (its hinge axle). Inward rotation $\varphi$ is limited to
$\varphi_{max}$ = 2° from the CAD pose by a hinge stop: a rigid tab of
length $L$ = 30 mm below the hinge swings outward by
$\delta = L \sin\varphi$, and a stiff bumper on the cup is placed at the
gap

$$g = L\left[\sin\varphi_{max} - \sin\varphi_0\right] \approx L \sin(\varphi_{max} - \varphi_0)$$

where $\varphi_0$ (0.88 to 0.99°) is the snug-fit closing already applied.
With $L$ = 30 mm, $g$ = 0.53 to 0.58 mm. The tab-bumper contact is
frictionless with a stiff penalty. A small
radial stabiliser spring (0.05 N/mm per petal) is reported separately.

### 5. Contact

FEBio `sliding-elastic`: penalty normal contact (scale 0.1 on the automatic
penalty) and Coulomb friction,

$$|\mathbf{t}_t| \le \mu\,p_n,\qquad \mu = 0.5$$

A contact face sticks while $|\mathbf{t}_t| < \mu p_n$ and slips at equality.
Slip area fraction = slipping contact area ÷ total contact area. Distal
share = axial contact force on the cushion ÷ axial reaction on the cup.

### 6. Materials

**Small-strain moduli and Poisson ratio.** For each nearly incompressible
material with shear modulus $\mu$ and bulk modulus $\kappa$:

$$\nu = \frac{3\kappa - 2\mu}{2\,(3\kappa + \mu)}$$

| κ / μ | ν | Used for |
|---|---|---|
| 100 | 0.4950 | skin, fat, silicone pad |
| 1000 | 0.4995 | muscle (incompressible, see §8) |

Under a mean pressure $p$, the volume change is $\Delta V/V \approx -p/\kappa$.
For muscle at 100 kPa that is 19% with κ = 0.52 MPa, but 1.9% with
κ = 5.2 MPa.

**Skin (HGO, Ni Annaidh 2012).**

$$W = \frac{c}{2}(\bar I_1 - 3) + \frac{k_1}{2k_2}\sum_{a=1,2}\left[e^{k_2 \langle E_a\rangle^2} - 1\right] + U(J),\quad E_a = \kappa_d(\bar I_1 - 3) + (1 - 3\kappa_d)(\bar I_{4a} - 1)$$

with $c$ = 0.2014 MPa, $k_1$ = 24.53 MPa, $k_2$ = 0.1327, $\kappa_d$ = 0.1404,
fibres at ±41° from the circumferential axis, and
$U = \tfrac{k}{2}\left[\tfrac{J^2-1}{2} - \ln J\right]$, $k$ = 100 c = 20.14 MPa.

**Fat and muscle (2-term symmetric Ogden; tissue set v2).**

$$W = \sum_{i=1}^{2}\frac{c_i}{m_i^2}\left(\bar\lambda_1^{m_i} + \bar\lambda_2^{m_i} + \bar\lambda_3^{m_i} - 3\right) + \frac{k}{2}(\ln J)^2,\quad c_1 = c_2 = c,\ m_2 = -m_1$$

Small-strain shear modulus $\mu = \tfrac{1}{2}(c_1 + c_2) = c$. Fat:
$c$ = 6.55 kPa (Mukhina 2020), $m$ = 5, $k$ = 0.655 MPa (100 c). Muscle:
$c$ = 9.5 kPa (Ranger 2022), $m$ = 5 (Mukhina 2020), $k$ = **9.5 MPa**
(1000 c). Sengeh 2016 (c ≈ 5.2 kPa, m ≈ 4.8) is no longer used.

**Silicone pad (neo-Hookean, Sanders 2004).** FEBio's uncoupled
Mooney-Rivlin has $\mu = 2(c_1 + c_2)$. With $c_2$ = 0 and the median
silicone-gel liner shear modulus μ = 23.3 kPa: $c_1$ = 23.3 / 2 = **11.65 kPa**,
$k$ = 100 μ = 2.33 MPa.

**Distal cushion (85A TPU, 2 walls, 8% gyroid).** Solid TPU small-strain
Young's modulus from the NinjaFlex Mooney-Rivlin fit ($c_1$ = 0.385,
$c_2$ = 1.47 MPa):

$$E_s = 6(c_1 + c_2) = 6 \times 1.855 = 11.1\ \text{MPa}$$

Gibson-Ashby scaling for a bending-dominated lattice at relative density
$\rho^*$ = 0.08:

$$E^* = C\,E_s\,(\rho^*)^n = 1 \times 11.1 \times 0.08^2 = \mathbf{71\ kPa}$$

Range: $n$ = 1.5 gives 252 kPa; $C$ = 0.5 gives 36 kPa. Poisson ratio 0.3
(compressible lattice).

**Bone.** Rigid in the LOTUS V5 runs (E = 17 GPa, ν = 0.3 when deformable).

### 7. Geometry

**Layer thickness.** At each outer skin vertex, the skin thickness is the
closest distance to the segmented inner skin surface, and skin plus fat is
the closest distance to the outer muscle surface. Where the layers would
self-intersect or reach bone, they are thinned locally to at most 90% of the
available depth.

**Petal closure.** Each petal is rotated about its hinge axis (the axle's
principal axis) until the 20th percentile of its pad nodes sits 0.1 mm
inside the skin. The rotations come out at 0.88 to 1.01°.

**Fibular relief trim.** Along the petal width, $u \in [0, 1]$ from the
posterior edge, with depth $d$ = 16 mm, flat fraction $f$ = 0.45 and ramp
fraction $r$ = 0.25:

$$w = \tfrac{1}{2} - \tfrac{1}{2}\cos\!\left(\pi\,\mathrm{clip}\!\left(\frac{u - f}{r}, 0, 1\right)\right),\qquad z_{top}(u) = z_{top,0} - d\,(1 - w)$$

On a petal about 100 mm wide at the top, that is 45 mm flat at full depth,
a 25 mm ramp, and 30 mm uncut.

**Proposed lip flare (not yet modelled).** Over the top $h$ = 15 mm, the
inner wall follows an arc of radius $R$ tangent to the wall. The outward
offset at height $y$ above the start of the flare is

$$\delta(y) = R - \sqrt{R^2 - y^2},\qquad \delta(15) = 4.0,\ 5.0,\ 6.8\ \text{mm for } R = 30,\ 25,\ 20\ \text{mm}$$

### 8. Muscle incompressibility check

`scripts/muscle_incompressibility_check.py`: a 20 mm cube of muscle on UT4
tetrahedra, with its sides and base held, is indented 5 mm by a 10 × 10 mm
indenter. The volume change is computed from the deformed tetrahedra,
$V = \sum \tfrac{1}{6}\,|(\mathbf{b} - \mathbf{a}) \times (\mathbf{c} - \mathbf{a}) \cdot (\mathbf{d} - \mathbf{a})|$:
−1.71% at κ = 100 μ, **−0.42% at 1000 μ** (used), −0.36% at 10 000 μ.

### 9. Post-processing quantities

* Volume ratio $J = \det \mathbf{F}$ (inversion when $J \le 0$; crushing
  shows as small $J$).
* Principal log strain $\ln \lambda_i$ and Green-Lagrange strain
  $\tfrac{1}{2}(\lambda_i^2 - 1)$ from $\mathbf{C} = \mathbf{F}^T\mathbf{F}$.
* Hydrostatic pressure $p = -\tfrac{1}{3}\mathrm{tr}\,\boldsymbol\sigma$.
* Force imbalance: |sum of reactions and applied loads| ÷ applied load
  (pass below 1e-3).
* **Stance check** (`scripts/stance_check.py`): the axial load $F(t)$ is
  read at each output time of the push. The state at the target load is
  linearly interpolated between the bracketing outputs $A$ and $B$:
  $w = (F_{target} - F_A)/(F_B - F_A)$, and each reported quantity is
  $q = q_A + w\,(q_B - q_A)$.

### 10. Gait Tester formulas

**Pressure estimate beyond the FEA fields** (calibrated on the 25 lbf run at
the end of tightening, $T$ = 111.2 N):

$$p_{peak} = c_{pk}\,\frac{k_{N,FEA}\,T + c_t\,F_{petals,axial}}{A}$$

* $k_{N,FEA}$ = 966.8 N / 111.2 N = **8.7**
* $A$ = 26 233 mm² (contact area)
* $c_{pk}$ = peak ÷ mean = 478 kPa ÷ (966.8 N / 26 233 mm²) = 478 / 36.9 = **13.0**
* $c_t$ = rise in normal force ÷ rise in petal axial load during the push =
  (1101 − 967) / (175 − 75) = **1.3**

**FEA fields.** For tension $T$ between saved states $i$ and $i+1$, every
nodal field (displacement, pressure, traction) is interpolated linearly in
tension:

$$q(T) = q_i + \frac{T - T_i}{T_{i+1} - T_i}\,(q_{i+1} - q_i)$$

Above the last solved state the fields are **not extrapolated**: the page shows that state and labels it NOT SOLVED.
A stance load $L$ adds the push increments of the 25 lbf run,
$q_{push,k} - q_{t=1}$, interpolated by $L$ against their axial-load
increments (40, 72 and 116 N for the 1, 2 and 3 mm push states). Above
116 N the push field is held at its last state and labelled. Nodal pressure
is the maximum of the adjacent contact faces, so the field peak matches
FEBio's face peak.

**Pressure gradient** (indicator of shear loading of the tissue):

$$|\nabla p|_{max} \approx \max_{\text{skin edges } (i,j)} \frac{|p_i - p_j|}{|\mathbf{x}_i - \mathbf{x}_j|}$$

## Materials (Scenario A: human tissue; Scenario B: prototype)

Both sets are in `configs/materials/`. Each lists every parameter with its
value, unit, source (DOI), convention conversion, calibration range and
uncertainty. The FEBio strain-energy forms and their mapping to the sources
and to ABAQUS are documented in
[docs/solver_selection.md](docs/solver_selection.md). In short, skin is HGO
(Ni Annaidh 2012), fat and muscle are 2-term Ogden (tissue set v2: Ranger
et al. 2022 and Mukhina et al. 2020), and silicone and the cushion are
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
