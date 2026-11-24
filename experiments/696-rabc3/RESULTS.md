# 696 — R-ABC3 (R-ABC2 + K 3,400 on def_day, A/B/C): NOT BUILT, K not adopted (readiness [4592ee], 2026-10-18)

**Outcome: never built.** R-ABC3 was preregistered as Branch A of 10-19 s2: it would be built only if read711 recorded "ADOPT K 3,400".
The coordinator confirmed that read711 says K stays 2,600 (64d49fda), so the waiter was stopped and no package exists.

## The read that decided it
- ku38 (K 4,200, CPM 3,800 on 446bj2z): entry aiUKL3IOpD4miuD5LpEu, 10-18 s2, overall 0.2736785. `read711.py` was run by readiness
  (read-only) and by the coordinator, with the same verdict.
- Falsifiers pass: Δpds raw 0, Δmse raw −2.2e-8, |Δnmae| 0.00044 < 0.002.
- Amendment 2 dial test: board Δreach raw vs 446bj2z = **+0.00942 > +0.0019**, so the f 0.61 dial is refuted. Condition (i) at the
  preregistered f 0.87 stands, and reach701.json says it fails (H1 reach falls with K, 3/3 seeds).
- Amended (ii) step-down: Δoverall K 3,400 (ku30) +0.0011785, K 4,200 +0.0008526. The best size is 3,400, its condition fails, and the
  step-down ends at 2,600.
- Raw DE vs 446bj2z at K 4,200: reach +0.0094, fid +0.0005, jac −0.0025, nmae +0.0004.

## What it means
- The board's reach rises with the call-set size (ku30 +0.0054, ku38 +0.0094 raw).
- Neither held-out H1 dial predicts that rise: f 0.87 says reach falls, and f 0.61 predicted −0.0016 at K 4,200.
- So the size gain is board-only, and per AGENTS §11.1c it does not enter October.
- Class: scorer mechanics (call-set size). Tries: 5 board sizes (1,600 / 2,000 / 2,600 / 3,400 / 4,200), 2 dials.
- What stays open: H1's call-set regime does not match the board's under size changes. f 0.87 retrodicted neither ku30 nor ku38, and
  f 0.61 retrodicted ku30 but failed ku38. Any later October call-set lever needs an instrument that retrodicts both points first.

## Actions (10-18)
- 03:22:14Z: rabc3_696 waiter (pid 60608, armed 10-17 06:33Z, never past its GO-rabc3 wait) stopped with SIGTERM. Its log ends with the
  trap's "NOT READY unexpected exit", then a STOPPED line with the reason, then `END`. GO-rabc3 was never touched.
- 03:22:43Z: `october_recipe.py decide rc_callset_cpm off`. The value is unchanged (DECIDED off since 10-06 on 523's K pair); the
  evidence now carries 711's rule and keeps 523's. The resolved October recipe (env + args) is byte-identical before and after;
  `check` rc 0, pending [].
- The coordinator moved R-ABC2H (a43a64d0) from 10-20 s1 to 10-19 s2, the slot Branch A would have used.
