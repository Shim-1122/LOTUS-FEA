# Stance-grip check: run `lotus_v5_stance_25lbf`

Run state: **converged**, last converged time 2.0 (t = 1 is the end of lace tightening, t = 2 is the full 10 mm push).
Lace tension: 111 N (25 lbf). Target axial load: 1200 N; distal share limit: 15%.

| t | push (mm) | axial load (N) | petals (N) | distal (N) | distal share | cup wall (N) | slip area | peak pressure (kPa) | min J |
|---|---|---|---|---|---|---|---|---|---|
| 0.25 | 0.0 | 13 | 10 | 3 | 22% | 0 | 0.22 | 55 | 0.92 |
| 0.50 | 0.0 | 24 | 21 | 3 | 14% | 0 | 0.14 | 81 | 0.88 |
| 0.75 | 0.0 | 32 | 28 | 4 | 12% | 0 | 0.13 | 91 | 0.87 |
| 1.00 | 0.0 | 38 | 34 | 4 | 11% | 0 | 0.12 | 94 | 0.86 |
| 1.10 | 1.0 | 89 | 78 | 9 | 10% | 2 | 0.28 | 100 | 0.86 |
| 1.20 | 2.0 | 142 | 121 | 14 | 10% | 7 | 0.42 | 97 | 0.86 |
| 1.30 | 3.0 | 198 | 161 | 20 | 10% | 16 | 0.50 | 109 | 0.85 |
| 1.40 | 4.0 | 256 | 200 | 27 | 11% | 31 | 0.53 | 121 | 0.83 |
| 1.50 | 5.0 | 317 | 234 | 34 | 11% | 49 | 0.55 | 124 | 0.80 |
| 1.60 | 6.0 | 382 | 268 | 42 | 11% | 73 | 0.54 | 129 | 0.76 |
| 1.70 | 7.0 | 452 | 301 | 50 | 11% | 99 | 0.51 | 137 | 0.72 |
| 1.80 | 8.0 | 525 | 337 | 60 | 11% | 129 | 0.49 | 142 | 0.67 |
| 1.90 | 9.0 | 601 | 369 | 70 | 12% | 161 | 0.49 | 157 | 0.63 |
| 2.00 | 10.0 | 680 | 406 | 81 | 12% | 195 | 0.48 | 174 | 0.58 |

## Verdict

* **Target NOT reached:** the largest axial load is 680 N (57% of the target) at a 10.0 mm push.

Quasi-static axial push only: no gait bending moment, torsion or dynamics. Grip here means the model's frictional and geometric resistance to the socket sliding proximally, with penalty contact (which underestimates peak pressure).
