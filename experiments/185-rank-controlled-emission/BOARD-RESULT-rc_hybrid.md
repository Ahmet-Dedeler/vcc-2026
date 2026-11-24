# `rc_hybrid` board result — **+0.0093325**, and rule 13's `reach` toll is REFUTED for profile-preserving arms

Entry `nQUJZhZb3MW2uFJHBX5e`, model name `The Mitochondria Declined To Comment`, **empty
description**, submitted 2026-09-03 11:20 UTC, published 11:41 UTC. **Overall `0.1834190095`
against a15's `0.1740865087` — `+0.0093325`, rank 19 of 631.**

Preregistration (committed before any line of the builder changed, unedited since):
`PREREGISTRATION-rc_hybrid.md`. Build and falsifiers: `experiments/191-rc-hybrid-build/RESULTS.md`.

## Member by member

| member | a15 scaled | rc_hybrid scaled | Δ scaled | Δ overall | a15 raw | rc raw | Δ raw |
|---|---:|---:|---:|---:|---:|---:|---:|
| `pds` | 0.745871 | 0.7458708290005253 | **−0.000000** | −0.000000 | 0.83712 | 0.83712 | **0.00000** |
| `mse` | 0.120057 | 0.12005661519665856 | **+0.000000** | +0.000000 | 0.87367 | 0.87367 | **0.00000** |
| `nmae` | 0.094239 | 0.094290 | +0.000051 | +0.000009 | 0.94330 | 0.94327 | −0.00003 |
| `fid` | −0.030397 | −0.029758 | +0.000639 | +0.000106 | 0.50337 | 0.50356 | +0.00019 |
| **`reach`** | 0.115368 | **0.170308** | **+0.054940** | **+0.009157** | 0.18117 | **0.22982** | **+0.04865** |
| `jac` | −0.000619 | −0.000253 | +0.000366 | +0.000061 | 0.03017 | 0.03031 | +0.00013 |
| | | | **+0.055995** | **+0.009333** | | | |

## The four pre-committed falsifiers, adjudicated

1. **`pds` or `mse` moves at all → the construction is broken.** **PASS, in the strongest possible
   form.** Both reproduce a15's board value **to sixteen significant figures**
   (`0.7458708290005253` and `0.12005661519665856`). The per-gene panel totals were digest-identical
   on 900/900 blocks locally and the board confirms the pseudobulk is the same object, on an
   independent scorer, against a reference we have never seen.
2. **raw `nmae` within 0.0010 of a15's.** PASS: **−0.00003**, thirty times inside the bound.
3. **challenge `n_pred` within ±2%.** PASS at build time: +0.13% (4,480 vs 4,474), and `fid`'s
   +0.00019 raw confirms the yield term did not engage.
4. **board `jac` within ±0.002 scaled.** PASS: **+0.000366**, and it moved the *right* way.

## The two hypotheses, adjudicated

| | predicted Δ overall | predicted Δ raw `reach` | realised |
|---|---:|---:|---|
| **H_transfer** | +0.0028 | +0.0151 | — |
| **H_toll** | −0.0024 | −0.0129 | — |
| **realised** | **+0.0093** | **+0.04865** | **H_transfer, by 3.33×** |

**H_toll is refuted.** Rule 13 charges emission-law arms a **−0.0146 scaled `reach`** toll, calibrated
on a21 and a22. Both of those *also* moved the delivered profile (a21: board `nmae` −0.064, `n_pred`
6,034 → 1,913, recall 0.669 → 0.271). `rc_hybrid` moved neither, and the toll's sign is wrong for it
by a wide margin. **The toll is a property of profile-changing emission arms, not of emission arms.**

## The finding that outlasts the arm: h1-test UNDERSTATES board `reach` by 3× on this class

h1-test read Δ`reach` **+0.03771 raw** (rc_hybrid 0.21681 against a15's twin mean 0.17910). The board
delivered **+0.04865 raw**. The implied transfer factor is **1.290** — against the **0.4277** that
164 measured and that every `reach` prediction in this project has passed through, and against the
**0.390** the a8/a15 pair implied.

**So the factor is not a constant, and for a designed-order arm it is an amplification, not a
shrink.** Two candidate mechanisms, neither yet separated:
- 191 measured that the incumbent's ordering defect is **2.6× smaller** on the challenge than on
  h1-test (Spearman −0.152 against −0.391), with a delivered span 82% of h1-test's. That predicts a
  *smaller* board gain, so it works against the observation and something else must dominate.
- The challenge's `n_real` distribution differs from h1-test's (median ~89–110 against 222) and
  `R_p = k*/n_real` is largest where `n_real` is smallest. `AGENTS.md` §2.1 says the panel was
  **deliberately binned by DE-gene count**, so the challenge has a guaranteed low-`n_real` stratum
  that h1-test does not match — and `docs/BOARD-LEDGER.md` §5 measures that ~three quarters of
  `reach` is between-perturbation heterogeneity. **This is the leading explanation and it is
  testable.**

## Where it leaves us

Board raw `reach` **0.18117 → 0.22982**. The shrunk cohort (31 teams at raw `pds` > 0.80 and raw
`mse` < 1.0) averages **0.2106**, so we have gone from **fourth-worst of 41** on the conditional
residual to **above the cohort mean**. Field best is Aginglab's 0.2621 and `luka`'s 0.2908.

Overall **0.18342, rank 19 of 631**; leader `cqawesome` 0.21712. **Gap 0.0337, down from 0.0438.**

**Largest single-member move this project has recorded** (`reach` +0.054940 scaled, against a8's
previous record of +0.021398) and the **second-largest single-arm gain** (behind a11's +0.01193).
