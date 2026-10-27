# Stance-grip check: run `lotus_v5_stance_25lbf`

Run state: **converged**, last converged time 2.0 (t = 1 is the end of lace tightening, t = 2 is the full 10 mm push).
Lace tension: 111 N (25 lbf). Target axial load: 1200 N; distal share limit: 15%.

| t | push (mm) | axial load (N) | petals (N) | distal (N) | distal share | cup wall (N) | slip area | peak pressure (kPa) | min J |
|---|---|---|---|---|---|---|---|---|---|
| 0.25 | 0.0 | 12 | 8 | 4 | 30% | 0 | 0.25 | 80 | 0.96 |
| 0.50 | 0.0 | 29 | 22 | 7 | 24% | 0 | 0.19 | 126 | 0.95 |
| 0.75 | 0.0 | 48 | 36 | 11 | 23% | 1 | 0.16 | 181 | 0.93 |
| 1.00 | 0.0 | 69 | 50 | 15 | 22% | 4 | 0.15 | 234 | 0.92 |
| 1.10 | 1.0 | 177 | 136 | 26 | 15% | 14 | 0.31 | 240 | 0.86 |
| 1.20 | 2.0 | 290 | 216 | 37 | 13% | 35 | 0.43 | 236 | 0.81 |
| 1.30 | 3.0 | 407 | 291 | 50 | 12% | 66 | 0.50 | 256 | 0.76 |
| 1.40 | 4.0 | 533 | 359 | 64 | 12% | 109 | 0.54 | 269 | 0.70 |
| 1.50 | 5.0 | 670 | 421 | 79 | 12% | 166 | 0.58 | 289 | 0.65 |
| 1.60 | 6.0 | 815 | 480 | 97 | 12% | 235 | 0.60 | 308 | 0.60 |
| 1.70 | 7.0 | 966 | 537 | 116 | 12% | 313 | 0.61 | 329 | 0.54 |
| 1.80 | 8.0 | 1126 | 593 | 138 | 12% | 398 | 0.62 | 349 | 0.49 |
| 1.90 | 9.0 | 1226 | 623 | 170 | 14% | 433 | 0.49 | 339 | 0.43 |
| 2.00 | 10.0 | 1322 | 647 | 206 | 16% | 469 | 0.40 | 385 | 0.37 |

## Verdict

* **Target reached:** the socket carries 1200 N after a 8.7 mm push.
* **Distal share at target:** 13% (PASS against < 15%). Petals: 616 N.
* **At target:** peak pressure ~342 kPa, slip area 0.52, min J 0.44.

Quasi-static axial push only: no gait bending moment, torsion or dynamics. Grip here means the model's frictional and geometric resistance to the socket sliding proximally, with penalty contact (which underestimates peak pressure).
