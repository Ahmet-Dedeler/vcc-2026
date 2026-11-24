# 694 — R-ABC: the October recipe as resolved, on def_day, on the real A/B/C bundle (2026-10-16, readiness [4592ee])

Preregistration: `PREREGISTRATION.md` (a37cbcae, written before any build). **Board read: 0.2743318, COMPOSES** (below). The coordinator submitted it
for 10-16 s1 under a documented override of gate v1 (`OVERRIDE-gate694.md`); see "The gate's false FAIL" below.

## What was built
- `october_recipe.py command` (plain, all switches decided) = 524dr2's env + `DEF_DAY_DE_UPRN=b_k100_q20r3p10` +
  `DEF_DAY_DE_DNMAG=d_m200`, args `1.75 0.5 noprog 4`. def_day.sh ran on 180's def-bundle with 453's operator 66610be7 pre-placed.
  The tree was builder e99bbd4a (dn_mag2 b9d6cfcf pinned) and emitter 0f56209d, after ka693b's KAU2 + DNV2 PASS (693 RESULTS).
- Package `submissions/694rabc-october-recipe-abc.vcc`, sha256 `125c1db94c3843f6d4f422a777f421180d030e065c75604c30d4a6bbd51a48d8`;
  content identity 494b7cce…; novelty rc 0.
- Run 10-15 23:39:53Z → 10-16 03:50:10Z (def_day 250 min), then relabel, prep, novelty, F1 and gate by 03:51:36Z.

## Checks (all on the finished build)
- **F1 bulk identity vs the board-scored 524dr2: IDENTICAL**: 0 of 900 blocks differ in any (gene) sum (bulkid658 on 524dr2's
  restored package). The re-level replayed 524dr2's path exactly. lvl0's rho equals 524dr2's lvl0 to all 16 digits
  (D 0.8716545108846331, E 0.8023769331241494, F 0.7724595749424453), and the factors went 1/1/1 → 1.06/1.1048/1.126 →
  1.06/1.1113/1.126, the same as 524dr2. So UPRN and DNMAG move no per-gene total, as designed, and pds / mse must read
  524dr2's to the board's precision.
- CLIFF OK (rho D 0.9735, E 0.9794, F 0.9885), converged; postbuild 42/42; zero-gene audit 0.
- The meta differs from 524dr2's in exactly de_dnmag + de_uprn (content-aware; see below), with module shas 7e6b3e52 / b9d6cfcf.
- DNMAG delivery (690h, per context): mean dn_frac 0.195 / 0.224 / 0.187, mean dn_moved 1,529 / 1,513 / 1,542, LOW_CONTEXTS none,
  0 invariant failures. That matches the DN build's A 0.193 / 1,512 etc.

## The gate's false FAIL (v1) and the fix (v2)
gate v1 ended the build "NOT READY" (03:51:36Z): its strict diff vs 524dr2 listed `bulk_load_sha256` and `sameline` besides
the two riders. Both are artifacts of a NEW work dir, not differences. The stack dump's `env` string records VCC_BULK_DUMP
(the path), and the same-line npz's `provenance` string records identity_source (the path). Every array in both npz files is
equal to 524dr2's (NaN-aware: pdown's NaN pattern is identical), and both strings are equal once each run's work dir is
replaced by one token. Science [d80a10] and the coordinator found the same by hand; the coordinator overrode and submitted.
- gate694 v2 accepts those two keys only after that content check. test694 11/11 (G6 provenance-only → PASS; G7 a changed
  array, G8 a provenance change beyond the path, G9 a changed same-line w → FAIL). `ready694r2.sh` re-gated the finished build
  (package sha unchanged, novelty re-run, F1 quoted): ALL GATES PASS, `READY rabc 125c1db9…` in `work/rabc694.r2.log`. That is the
  same package; nothing was armed on that log.
- The same false FAIL would have hit any replacement final entry diffed against R0 (each def_day run has its own work dir) and
  the rider rehearsal diffed against R-NEW. `scripts/def_metadiff.py` (content-aware diff, tests 8/8 in
  `scripts/test_def_metadiff.py`) is now the diff in OCT22-CHECKLIST §C 6a and in rider690.sh.

## Timing (for October)
def_day with UPRN + DNMAG v2: dump + lvl0 81 min, lvl1 84, lvl2 84 = 250 min, vs 524dr2's 101 (~34 per level). DNMAG v2 accounts for
~15 min per level (693), so most of the extra ~50 min per level build is UPRN's rank-neutral UP pass. R0 on the day ≈ release +
4 h 15 m (OCT22-CHECKLIST §D). Speeding UPRN up is a 693-style candidate (byte-identical rewrite), if R0's latency matters.

## Prediction and read (bins frozen in PREREGISTRATION.md)
Δ vs 524dr2 0.2730889, predicted +0.0013033 → 0.2743922 (all in nmae). COMPOSES [+0.0010, +0.0016]; SUBADDITIVE
[+0.0004, +0.0010); NULL/HARM < +0.0004; SUPERADDITIVE > +0.0016. F1 (pds / mse identical) is expected to hold on the board, since
the bulk sums are identical.

## Board read (10-16 s1, entry PooEbFClgTLMQoXiOCQP; captured `data/board-entries/PooEbFClgTLMQoXiOCQP.json`)
**0.2743318, Δ +0.0012429 vs 524dr2 (dVHLws6K, 0.2730889): COMPOSES** ([+0.0010, +0.0016]), 95 % of the additive prediction
+0.0013033 (read 0.2743318 vs predicted 0.2743922).

| member | raw R-ABC | raw 524dr2 | Δ raw | Δ scaled | Δ overall |
|---|---:|---:|---:|---:|---:|
| pds | 0.8653102935711632 | same | **0 exactly** | 0 | 0 |
| mse | 0.6240714820375378 | same | **0 exactly** | 0 | 0 |
| fid | 0.5143453637781704 | same | **0 exactly** | 0 | 0 |
| jac | 0.04081154789269244 | same | **0 exactly** | 0 | 0 |
| nmae | 0.889059778339481 | 0.8935180202928171 | −0.004458 | +0.007472 | +0.0012453 |
| reach | 0.28654278771012615 | 0.286553828222712 | −1.10e-5 | −1.42e-5 | −0.0000024 |

- **F1 held** (pds and mse identical to every digit, as the bulk identity predicted) and **F2 held** (fid and jac identical:
  0 call flips). **F3 is reported**: |Δreach raw| 1.10e-5 is just above the 1e-5 report line (DN on 446bj2z's path read +9.6e-8).
  It costs −0.0000024 overall and is the only member the two riders move besides nmae.
- nmae raw −0.004458 vs the additive −0.004696 (UPRN −0.00305 on 446bj2z + DN −0.001646 on 630arm): the riders compose on
  the def_day path with tercilec + the CLEAN4 hook at 95 % of their separate reads.
- Class: scorer mechanics (nothing fitted to A/B/C). Tries: UPRN and DN were each read once on the board, and this is the
  first read of both on the def_day path. Held-out: H1 (630 / 631), which over-read both. No sealed confirmation (R7 is spent
  on F + the hook).
- **Consequence:** the board has now scored the exact recipe R0 runs (same resolved env, same tree, same def_day path).
  October ships as decided.

## test694 G6 was flaky (found and fixed 2026-10-18 ~12:00Z, readiness)
- test694 G6 ("new work dir, arrays identical → ALL PASS") failed in some processes: 10/11 on 10-18, 11/11 on 10-16. It depended on
  PYTHONHASHSEED (seeds 0–3 PASS, 4–5 FAIL). gate694's `_norm` replaced the work-dir tokens in SET order. In a temp dir under
  macOS's `/var` → `/private/var` symlink, the unresolved `/var/...` token could replace inside the resolved `/private/var/...`
  string first, leaving `/private<WORK>` ("env differs beyond the work dir").
- Fix: tokens are replaced longest first, as `scripts/def_metadiff.py._tokens` already does. test694 11/11 for seeds 0–7.
- **Verdicts stand.** A wrong normalisation leaves the strings unequal, so this bug could only cause a false FAIL, never a false
  PASS. R-ABC's gate PASS (ready694r2) and its board read are unaffected. No live gate uses gate694 any more (695 onward use
  def_metadiff).
- Separately, the same day, def_metadiff gained an absolute float floor (`--float-atol` 1e-12; f9b6e09d, 698's false FAIL). It does
  not change this test's outcome: G6 failed identically with the old and new def_metadiff before this fix.
