# Review of the supplied hand-off package (provenance only)

The hand-off archive `Claude_FEA_Handoff_No_V6_CAD.zip` was inspected before
any implementation choices were made. Its content is summarised below. No
patient geometry and no historical result files are included in this
repository. The historical work was used only to understand earlier
assumptions and failure modes. It was **not** used as a reference solution,
and the new model was **not** tuned to reproduce it.

## 01_Anatomy (supplied limb model, private)

| File | Content found | Assessment |
|---|---|---|
| `Limb w- Bone - Part 1.step` | B-rep, 106 B-spline faces, unit METRE; OCC tessellation closed after orientation repair; 143.7 x 152.1 x 289.1 mm; 3.149e6 mm^3 | skin envelope with flat proximal cut |
| `Limb w- Bone - Part 2.step` | AP242 tessellated shell (gmsh/OCC importer crashes, so it is read by `rlfea.geometry.step_tessellated`), closed, 95 x 74 x 115 mm, 3.00e5 mm^3 | shape of a **distal femur** with condyles (heuristic, unconfirmed) |
| `Limb w- Bone - Part 3.step` | tessellated, closed, 48 x 24 x 47 mm, 2.55e4 mm^3 | shape of a **patella** (heuristic, unconfirmed) |
| `Limb w- Bone - Part 4.step` | tessellated, closed, 28 x 26 x 100 mm, 1.99e4 mm^3, posterior | **unidentified** slender bone |
| `Limb w- Bone - Part 1 (1).step` | header and product data only; **no geometry** | EMPTY body |
| `Limb_with_Bones.zip` | the same five STEP files | duplicate container |

**Missing anatomy:** no body matches a tibia. No bone articulates with the
femoral condyles, and nothing has the size of a proximal tibia. A fibula
cannot be confirmed either. The bones are modelled as supplied. Missing
anatomy is reported explicitly and never created.

All bones lie inside the skin envelope with no intersections. The minimum
clearances are 1.47 mm for the femur (at the proximal cut), 3.95 mm for the
patella and 8.06 mm for part 4. No duplicate bodies were found apart from the
zip container. Handedness (left or right) cannot be determined from the files.

## 03_Previous_Tissue_Analysis
These are patch-level studies (tied planar patches, an Ogden or HGO custom
solver, and a "screening" load factor). The ABAQUS decks in
`FEA_Preliminary_Package/abaqus_decks` were **not executed**. The material
values listed there (skin mu 22 kPa with alpha 6, and others) were described
as "user's provisional ranges" with no primary source. They are **not reused**;
the parameter set in `configs/materials/human_tissue.yaml` was sourced
independently, with DOIs.

## 04_Full_Limb_Experimental
These are custom-solver experiments with prescribed dead tractions, not
socket contact. One snapshot reports 23.5 mm maximum displacement at 7.5% of
its normal load. Another reports 12.3 mm at 5.625 N strand tension on 8,193
nodes and 4,665 quadratic tetrahedra. These are unvalidated historical results
from an experimental formulation. They are listed here for provenance only.
`Full_Limb_HGO_Abaqus_Comparison_UNEXECUTED.inp` was never run in ABAQUS.

## 05_Tissue_Report
This is a previously delivered PDF and its LaTeX. It predates any full-limb
target-load solution.

## Design intent
V6 socket CAD is excluded and rejected by the user. The socket in this
repository is a fresh parametric model built from the stated concept (petals,
silicone interfaces, RevoFit lace, hamstring clearance, ribbed hinges, distal
cushion and cup, and an Anclote puck, plate and Bulldog adapter below). It is
labelled PROVISIONAL until verified component dimensions and routing are
supplied.
