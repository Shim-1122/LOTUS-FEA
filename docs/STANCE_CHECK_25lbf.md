# Stance-grip check: run `lotus_v5_stance_25lbf`

Run state: **converged**, last converged time 2.0 (t = 1 is the end of lace tightening, t = 2 is the full 10 mm push).
Lace tension: 111 N (25 lbf). Target axial load: 1200 N; distal share limit: 15%.

| t | push (mm) | axial load (N) | petals (N) | distal (N) | distal share | cup wall (N) | slip area | peak pressure (kPa) | min J |
|---|---|---|---|---|---|---|---|---|---|
| 0.25 | 0.0 | 16 | 12 | 4 | 27% | 0 | 0.21 | 57 | 0.92 |
| 0.50 | 0.0 | 31 | 23 | 6 | 21% | 1 | 0.14 | 79 | 0.89 |
| 0.75 | 0.0 | 43 | 34 | 7 | 17% | 2 | 0.16 | 111 | 0.85 |
| 1.00 | 0.0 | 55 | 44 | 8 | 15% | 3 | 0.16 | 148 | 0.82 |
| 1.10 | 1.0 | 107 | 86 | 13 | 12% | 8 | 0.25 | 129 | 0.85 |
| 1.20 | 2.0 | 162 | 126 | 18 | 11% | 17 | 0.35 | 125 | 0.85 |
| 1.30 | 3.0 | 218 | 165 | 24 | 11% | 28 | 0.40 | 114 | 0.85 |
| 1.40 | 4.0 | 276 | 202 | 31 | 11% | 43 | 0.43 | 121 | 0.82 |
| 1.50 | 5.0 | 339 | 238 | 38 | 11% | 63 | 0.44 | 124 | 0.78 |
| 1.60 | 6.0 | 405 | 273 | 46 | 11% | 87 | 0.43 | 129 | 0.74 |
| 1.70 | 7.0 | 475 | 307 | 55 | 11% | 114 | 0.42 | 137 | 0.70 |
| 1.80 | 8.0 | 534 | 333 | 66 | 12% | 134 | 0.32 | 183 | 0.65 |
| 1.90 | 9.0 | 613 | 367 | 76 | 12% | 171 | 0.35 | 200 | 0.61 |
| 2.00 | 10.0 | 694 | 400 | 87 | 12% | 207 | 0.36 | 220 | 0.57 |

## Verdict

* **Target NOT reached:** the largest axial load is 694 N (58% of the target) at a 10.0 mm push.

Quasi-static axial push only: no gait bending moment, torsion or dynamics. Grip here means the model's frictional and geometric resistance to the socket sliding proximally, with penalty contact (which underestimates peak pressure).
