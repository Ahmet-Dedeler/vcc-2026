# 458 — the reach contingency: a design sketch (written 2026-09-26 12:01Z, before any 458/458d read)

Use this only if 458d reads fid ≥ +0.0015 and reach ≤ −0.005 raw with overall ≤ 0 (`PREREGISTRATION.md`, addendum
11:31Z). Nothing here is built.

## Why a programme prior can cost reach
- `reach` walks the reference's DE genes in OUR significance order and keeps the deepest prefix that is ≥ 90 % sign-pure.
- The emitter (`src/vcc2026/rc_emit_k.py`, block branch around line 950) sets each block gene's target LFC to
  `psign · |E0|`. The design's delivered LFC `E0` survives only for the strongest `1 − q` of the block (`protect_q` 0.8).
- In a small-gate target (n_real 2–17), a reference DE gene is often a weak block gene, so it is unprotected. If the
  prior flips a correct design sign on a gene that ranks first, reach goes from 1 to ~0.
- 462 §11: that one mechanism is all of Pan's H1 reach loss (4 targets, n_real 2/2/8/17). 446bja's board read had
  reach −0.0079 raw.

## Two candidate fixes (each needs a knob, default off, a KA, an H1 price, then a board slot)
1. **Demote the flips.** Where the prior disagrees with the design sign (`sign(psign) != sign(E0)`), deliver the
   flipped gene at reduced magnitude, e.g. `0.5 · |E0|` or the MDE floor.
   - It stays called (fid keeps the sign) but sinks in the significance order, below the design-signed genes.
   - reach's top prefix is then design-signed. Cost: a few flipped genes may fall out of significance (n_pred ≫
     n_real, so fid barely moves), plus the nmae of their smaller magnitude.
2. **Per-target rank protection.** Protect the design sign on each target's top-m ranked block genes (by `|E0|` or by
   the pass-2 |z|), not a global quantile.
   - m ≈ 10 covers reach's first-error depth.
   - This is `protect_q` computed per target over its most significant genes rather than over the block's `|E0|`
     distribution.

## How to price
- 462l (`drive462.sh`, queued) gives the H1 profile of Arc's own programme through the shipped harness. If reach drops
  there too, the mechanism is quality-independent and fix 1 or 2 is needed for any prior.
- Implement the knob in `rc_emit_k` behind `lfc_target["demote_flips"]` / `["protect_top_m"]` (default absent =
  bit-identical, proven by 435's H1 known answer). Price it with `run462.py` against `arcfull` on H1 (3 seeds; fid,
  reach, nmae, jac).
- Then a board arm: 458d + the knob, vs 458d.

## Update 12:07Z 09-26: what the emitter and the record say (read before building either fix)
- **The emitter's order is the design's, not the prior's.** Called genes rank by `|key|` (`rc_hybrid_block_k`:
  `ranked = argsort(-(order_abs + jitter + block_bonus))`), and each gene's Wilcoxon z follows that rank. The prior
  changes only the fold-change sign through the A-channel mean. A flipped gene therefore keeps its rank.
- **Fix 1 has a closed cousin.** "A reach order boosted by design–prior agreement" lost reach monotonically on H1
  (390 §6–8, −0.013 at β 4; FINDINGS §2, 09-13 list). Demoting flips through `order_key` is that construction's mirror.
  Do not build it without a new argument for why it would behave differently.
- **Fix 2 does not reach the failing genes.** A small-gate target's few reference DE genes can sit hundreds of places
  down our 2,600-gene call set, so protecting our top-m ranked genes would not cover them.
- **Net reading:** a better prior raises gate-gene sign accuracy on average (441: called gate genes 0.567 → 0.598 with
  CLEAN4). Its reach cost is small-gate volatility, near coin-flip per target. With 900 target-contexts on the board, a
  reach drop within ~2 × the ordering channel's noise (sd 0.0045 overall) is not evidence of a systematic cost.
- So if 458d reads fid-positive with reach ≤ −0.005 raw, first check whether the loss concentrates in a few small-gate
  targets. The board gives no per-target data, so this is only argued from H1. If it cannot be separated from noise,
  prefer the dose-curve partner (458e) to a new emitter knob.

## Update 13:23Z 09-26: reach and fid should move TOGETHER for a programme change (406's descriptive table)
- On H1, reach comes almost entirely from weak targets. By n_real tercile: weak (mean n_real 6) reach 0.666; mid (187)
  0.059; strong (2,216) 0.016.
- In weak targets the breaking gene is **never in the head**: 60 % block, 40 % uncalled. Block genes carry the
  prior's sign (unprotected mostly).
- So a programme that is better for the recipient should RAISE reach along with fid. A worse one lowers both.
  446bja's reach −0.0079 fits Nadig's Jurkat programme being worse than CLEAN4 on A's weak-target genes.
- Reading rule for 458d: fid and reach moving the same way is the expected, interpretable case. fid up with reach
  down (the contingency) would mean the programme helps the many calls but hurts the few genes weak targets hinge on.
  That is a mixed signal, and 458e (w 0.5) is its first follow-up.

## CORRECTION 14:10Z 09-26: fid and reach do NOT co-move for a programme prior (462l supersedes the update above)
- The update above expected a better programme to raise reach along with fid. 462l refutes it. On H1, Arc's own clean
  programme (oracle-grade on that line) raises fid +0.0198 and lowers reach −0.0388.
- The loss comes from the same four small-gate targets as Pan's: CNN3, COL6A1, ZNF664, AASS. Their true DE genes run
  against the line's common programme, so any per-gene prior pays that reach.
- Reading rule for 458d, replacing the one above: expect fid up (and nmae better) with a reach cost of roughly
  446bja's size. Adopt on fid ≥ +0.0015 with overall ≥ 0, as preregistered.
  - The contingency still applies when overall ≤ 0.
  - Its follow-up stays 458e. A per-gene prior cannot know target-specific exceptions, and the emitter fixes above
    have a closed cousin or miss the genes.
