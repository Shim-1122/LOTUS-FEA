# Load check: run `lotus_v5_stance_60lbf_disp`

Run state: **failed**, last converged time 1.79091 (t = 1 is the end of lace tightening, t = 2 is the full 10 mm push).
Lace tension: 267 N (60 lbf). Target axial load: 1200 N; distal share limit: 15%; petal share minimum: 85%.

| t | push (mm) | axial load (N) | petals (N) | petal share | distal (N) | distal share | cup wall (N) | slip area | peak pressure (kPa) | min J |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.25 | 0.00 | 35 | 26 | 75% | 9 | 25% | 0 | 0.21 | 143 | 0.94 |
| 0.50 | 0.00 | 83 | 57 | 69% | 18 | 22% | 8 | 0.19 | 267 | 0.90 |
| 0.75 | 0.00 | 143 | 83 | 58% | 29 | 20% | 31 | 0.22 | 384 | 0.85 |
| 1.00 | 0.00 | 203 | 105 | 52% | 40 | 20% | 57 | 0.16 | 478 | 0.79 |
| 1.10 | 1.00 | 343 | 197 | 57% | 53 | 15% | 93 | 0.21 | 487 | 0.74 |
| 1.20 | 2.00 | 484 | 282 | 58% | 66 | 14% | 136 | 0.27 | 489 | 0.70 |
| 1.30 | 3.00 | 610 | 364 | 60% | 83 | 14% | 163 | 0.23 | 505 | 0.64 |
| 1.40 | 4.00 | 764 | 437 | 57% | 96 | 13% | 230 | 0.28 | 567 | 0.60 |
| 1.50 | 5.00 | 919 | 509 | 55% | 112 | 12% | 298 | 0.28 | 653 | 0.55 |
| 1.60 | 6.00 | 1083 | 589 | 54% | 130 | 12% | 365 | 0.28 | 672 | 0.51 |
| 1.70 | 7.00 | 1252 | 662 | 53% | 150 | 12% | 444 | 0.27 | 829 | 0.46 |

## Verdict

* **Target reached:** the socket carries 1200 N after a 6.7 mm push.
* **Distal share at target:** 12% (PASS against < 15%).
* **Petal share at target:** 53% = 639 N (FAIL against >= 85%, i.e. >= 1020 N). Cup wall: 419 N.
* **At target:** peak pressure ~780 kPa, slip area 0.27, min J 0.47.

Quasi-static axial push only: no gait bending moment, torsion or dynamics. Grip here means the model's frictional and geometric resistance to the socket sliding proximally, with penalty contact (which underestimates peak pressure).
