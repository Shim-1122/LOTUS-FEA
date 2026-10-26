# How hard must the socket grip? (walking estimate)

These are hand calculations with rounded literature values. They size the
lace tension and are not a substitute for a body-weight FE case or a bench
test.

## Patient and prosthesis

Patient: 50 years old, 175 cm (5'9"), 81.6 kg (180 lb) before amputation.

| Quantity | Value | Basis |
|---|---|---|
| Mass lost with the limb | 3.5 kg | foot 1.45% + 60% of shank 4.65% of body mass (Winter segment fractions); the residual tibia is about 150 mm long |
| Body mass after amputation | 78.1 kg | |
| Prosthesis mass (socket, pylon, foot) | 1.5 kg | **assumed**; weigh the device |
| Standing weight with prosthesis | 781 N | |

## Loads the socket-limb interface carries in walking

| Phase | Load on the socket | Estimate |
|---|---|---|
| Stance, peak (loading response / push-off) | axial, pushes the limb INTO the socket | 1.2 × 781 N − prosthesis weight ≈ **923 N (207 lbf)** |
| Late swing | axial, pulls the socket OFF the limb | 1.5 kg × (g + ω²r) = 1.5 × (9.81 + 6.5² × 0.30) ≈ **34 N (8 lbf)** |

* **Stance.** 1.2 BW is a conservative peak vertical ground reaction force
  for level walking (the prosthetic side is usually somewhat lower).
* **Swing.** The swing estimate uses a peak knee angular velocity of about
  6.5 rad/s, with the prosthesis centre of mass about 0.30 m from the knee.
* **The distal end is not load-bearing.** The cushion is comfort only, so
  the stance load has to go into the petals through friction and the
  socket's taper. **Stance, not swing, is the case that drives the design.**
* **Design target at 2 × the requirement:** about **1.85 kN** axial in
  stance and **70 N** pull-out in swing.

## Lace tension that would provide it (friction-only estimate)

In the LOTUS V5 run with tube lacing, the total pad-to-skin normal force
grows almost in proportion to the line tension. The ratio is
7.3, 7.8, 8.0 and 8.1 N per N of line tension at 25, 50, 75 and 100% of
75 N. If the whole interface could slide at once, the friction capacity
would be μ × ΣN ≈ μ × 8.1 × T.

| Skin-silicone friction μ | T for 1 × stance (923 N) | T for 2 × stance (1.85 kN) |
|---|---|---|
| 0.5 (model default) | 228 N (51 lbf) | **456 N (102 lbf)** |
| 0.8 | 142 N (32 lbf) | 285 N (64 lbf) |

The 2 × target needs about 64 to 102 lbf of line tension. That is within
the 200 lbf RevoFit rating, and well within the roughly 300 lbf the team
reports from independent testing. Swing needs only about 2 to 4 lbf, so it
never governs.

## Caveats

* **Estimate, not proof.** Friction-only capacity ignores the socket taper,
  which helps in stance and makes this estimate conservative. It also
  assumes the whole interface slides at once.
* **Pressure check.** Scaling the 75 N result linearly, 456 N of line
  tension gives a peak pad pressure of about 1.2 MPa (196 kPa × 6.1). That
  is far above the 75 N run and likely uncomfortable, especially at the
  pad edge where the model already crushes muscle at 75 N. The real
  tissue response is nonlinear, so this must be checked in the FE model
  before choosing a tension.
* **Next steps.** Weigh the device, measure skin-silicone friction, and run
  the FE model at the proposed tension with a force-controlled 1.85 kN
  body-weight case.
