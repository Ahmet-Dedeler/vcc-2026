# BOARD RESULT — `h1union` (g=0.28): **0.2118310, rank 3 of 656. Branch A fires.**

Entry `iUziqE1gewfcPKzYUDaB`, `The Nucleus Declined To Comment`, description empty.
Submitted 2026-09-04 03:50:39Z, published ~04:16Z. Baseline `cd4f034` = **0.2049196**.

## 1. **Δ = +0.0069114. All six members up.**

| member | `cd4f034` | `h1union` | Δ |
|---|---:|---:|---:|
| `pds` | +0.73780 | +0.74537 | **+0.00757** |
| `mse` | +0.15928 | +0.17957 | **+0.02029** |
| `nmae` | +0.10486 | +0.11074 | **+0.00588** |
| `fid` | −0.02013 | −0.01648 | **+0.00365** |
| `reach` | +0.21930 | +0.22294 | **+0.00364** |
| `jac` | +0.02841 | +0.02884 | **+0.00043** |
| | | **sum** | **+0.04147 → /6 = +0.0069114** |

raw `mse` **0.83654 → 0.81726**; raw `pds` **0.83228 → 0.83572**.
Implied cosine **0.09018 → 0.09947, +10.3045% relative.**

**Rank 3 of 656.** Gap to `cqawesome` **+0.01248**, was +0.01939 — **36% of the gap to first
closed in one slot.** Cosine ratio leader/us **1.343×**, was 1.481×.

## 2. Branch determination — **A**, per `DECISION-TREE-AMENDED.md`

Overall **+0.0069 ≥ +0.0020** and raw `pds` **+0.00344 ≥ +0.0010**. **Rule *and* sources
validated.** Next slot = **223's ranked stack as riders, re-priced against `h1union`** (208's rider
rule), each owing a leave-top-3-out retention test. **225 is already executing exactly this.**

**Branch F did not fire, but the concentration warning was directionally right.** `pds` realised
**70%** of its 3-cancer prediction (+0.00344 vs +0.004887) while cosine realised ~100%. So the
concentrated channel *did* underperform the distributed one — just not catastrophically. **Keep the
retention test; treat it as a discount, not a veto.**

## 3. 🔴 My reissued band was WRONG, and the original preregistration was right

`DECISION-TREE-AMENDED.md` reissued **+0.0010 to +0.0045, central +0.0030**. The truth is
**+0.0069** — outside my band, inside the *original* preregistered **+0.0055 to +0.0105**.

**What I got right:** `kolf2` does inflate the 4-fold Δ`pds`, and the 3-cancer `pds` number was the
honest one — indeed it was still 43% optimistic (realised 70% of it), so the 4-fold number would
have been far worse.

**What I got wrong:** I built the overall band from the **macro** cosine (+7.6803%) multiplied by
216's `k_cos` 1.384. I should have used the **aggregate** cosine (+10.3065%), which 222 supplied in
the same table and which landed at **+10.3045%** — a match to four significant figures. Applying a
transfer factor to the aggregate on top of that was double-counting.

**And `nmae` rose (+0.00588) rather than falling.** 222's exact-identity instrument predicted
−0.000108 overall. **222's own caveat was right and its instrument was wrong** — it warned the
instrument "mis-signs its own family's board point" and told me to read it as a magnitude bound, not
a sign. The magnitude bound held; the sign did not. **Trust that caveat in future.**

## 4. 🔴 The open question this creates — and it is worth more than the +0.0069

**Two arms in the same family now have two different PLOLO→board cosine factors, on the matched
`cos_agg` functional both times:**

| arm | PLOLO rel `cos_agg` | board rel cosine | factor |
|---|---:|---:|---:|
| `cd4f034` (216 §4.4, matched functional) | +17.98% | +24.89% | **1.384** |
| `h1union` (222, three cancer folds) | +10.3065% | +10.3045% | **1.000** |

**This is not a validated instrument and must not be used as one until the discrepancy is
explained.** Candidate explanations, none yet tested: the two numbers may be aggregated over
different fold sets (222 excluded `kolf2` from its headlines; 216's ladder may not have); `cos_agg`
may be computed over different gate/tag conventions — **the exact trap 216 itself named and then may
have re-entered**; or the factor may genuinely depend on where the added source sits relative to the
existing stack.

**Sizing every future arm depends on resolving this.** 221 owns local→board calibration and has been
sent both points.
