# 696 — R-ABC3: R-ABC2 + K 3,400, on def_day, on A/B/C (10-19 s2 Branch A; preregistered 2026-10-16 ~13:00Z, readiness [4592ee])

Written before ku38's read (10-18 s2), before read711, and before R-ABC2 exists. Built ONLY if the coordinator records read711's
"ADOPT K 3,400" (711 Amendment 2: ku38 validates the f 0.61 dial). Otherwise it is never built, and the coordinator fills 10-19 s2.
Branch B (a seed-only re-roll of R-ABC2) was dropped: the builder has no seed knob (SEED and RC_JITTER_SEED_BASE are module
constants), so it would need production edits in the release window. Rule 13's re-roll sd is already confirmed on 12 board pairs.

## What is built
R-ABC2's resolved recipe (695 `work/recipe695.json`) + `DEF_DAY_RC_CALLSET_CPM=3000` (K = 400 + 3000 = 3,400), the same args,
built exactly as R-ABC2 (def_day.sh on 180's def-bundle, 453's operator 66610be7 pre-placed). The driver stops unless the resolved env
is R-ABC2's env plus exactly that one export. It is the set-value rehearsal of the K switch on def_day before the 10-20 freeze.
The export is env-only (530), so def_day.sh is not edited and 690k's R-KA proof stands.

## Prediction
pred = R-ABC2's board read + ku30's board Δoverall (+0.0011785: K 3,400 vs 446bj2z's K 2,600). K changes only the DE call set,
so pds and mse must equal R-ABC2's exactly (bulk identity is checked before READY).

## Bins on D = R-ABC3 − pred (overall, board; a K change re-draws pass 1: re-roll sd ≈ 0.00025, ±0.0005 = 2 sd)
| bin | D | consequence |
|---|---|---|
| ABOVE | ≥ +0.0005 | K 3,400 ships; record the interaction |
| **COMPOSES** | (−0.0005, +0.0005) | K 3,400 ships |
| BELOW | ≤ −0.0005 | K reverts to 2,600 for R0 (the latest-adopted switch reverts first) |

## Falsifiers
- Pre-slot: the bulk is R-ABC2's (bulkid658: identical per-(context, target, gene) sums vs R-ABC2's restored package). Otherwise NOT READY.
- Pre-slot: the content-aware config diff vs R-ABC2 (scripts/def_metadiff.py) is exactly `rc_callset_cpm`, with value 3000.
- Board: pds and mse raw identical to R-ABC2's (exact members). A mismatch means the bulk moved; diagnose before reading D.

## Class and tries
Scorer mechanics (the call-set size; nothing fitted to A/B/C beyond 711's rule). K's board points: 1,600, 2,000, 2,600, 3,400 and 4,200
(711). 711's rule picks among 3 sizes with a +0.0005 step margin and an H1 mechanism gate.

## Resolved input (10-18 02:2xZ, before read711 and before any R-ABC3 build)
R-ABC2's board read = 0.2739849 (10-18 s1, nnHQfcZHsT3as9EafA3i) ⇒ **pred = 0.2751634** if built. COMPOSES (0.2746634, 0.2756634);
BELOW ≤ 0.2746634; ABOVE ≥ 0.2756634. pds and mse raw must equal R-ABC2's exactly: 0.8646413972500929 and 0.624402700065685.
