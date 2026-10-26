# Stance-grip check: run `lotus_v5_stance_25lbf`

Run state: **failed**, last converged time 1.3 (t = 1 is the end of lace tightening, t = 2 is the full 10 mm push).
Lace tension: 111 N (25 lbf). Target axial load: 1200 N; distal share limit: 15%.

| t | push (mm) | axial load (N) | petals (N) | distal (N) | distal share | cup wall (N) | slip area | peak pressure (kPa) | min J |
|---|---|---|---|---|---|---|---|---|---|
| 0.25 | 0.0 | 24 | 19 | 5 | 21% | 0 | 0.26 | 77 | 0.85 |
| 0.50 | 0.0 | 54 | 43 | 11 | 20% | 0 | 0.15 | 149 | 0.78 |
| 0.75 | 0.0 | 80 | 63 | 17 | 21% | 0 | 0.11 | 287 | 0.70 |
| 1.00 | 0.0 | 100 | 75 | 25 | 25% | 0 | 0.09 | 478 | 0.62 |
| 1.10 | 1.0 | 140 | 110 | 30 | 21% | 0 | 0.11 | 519 | 0.59 |
| 1.20 | 2.0 | 172 | 135 | 36 | 21% | 0 | 0.11 | 569 | 0.61 |
| 1.30 | 3.0 | 216 | 175 | 41 | 19% | 0 | 0.12 | 647 | 0.59 |

## Verdict

* **Target NOT reached:** the largest axial load is 216 N (18% of the target) at a 3.0 mm push.
* The run stopped before the full push, so this is a lower bound on the grip, not its capacity.

Quasi-static axial push only: no gait bending moment, torsion or dynamics. Grip here means the model's frictional and geometric resistance to the socket sliding proximally, with penalty contact (which underestimates peak pressure).
