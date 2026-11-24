"""533: Draft v15 text blocks, part C (§3(k) additions and the new §3(l)). See v15_blocks_a.py for the conventions."""

S3K_ADD = """

*Six more single-channel entries (2026-09-27 … 09-29; C-83; Figure F-12b).* Raw members, context-averaged, arm minus its
base (bold: the member the change was built to move):

| entry (base) | change | pds | mse | nmae | fid | reach | jac |
|---|---|---:|---:|---:|---:|---:|---:|
| 446bj2 (446b) | per-target bulk key, leave-one-donor-out | +0.000208 | **−0.014473** | +0.000060 | −0.000022 | −0.001963 | +0.000016 |
| 453dr (446b) | the day-of refit path at CD4 ×1 (§3(h)) | −0.001345 | **+0.015033** | −0.000188 | −0.000036 | −0.003436 | −0.000004 |
| 457dr (453dr) | CD4 ×1.5 on the refit path | +0.001713 | **−0.014313** | +0.000200 | −0.000034 | +0.000711 | +0.000004 |
| 446bd2 (446b) | a fifth donor line in the bulk | **+0.002148** | **−0.000471** | +0.000121 | −0.000020 | −0.002314 | +0.000000 |
| 446br2 (446b) | zero-gene repair floor −2 | +0.000007 | +0.000070 | **−0.001413** | 0 | −0.000039 | 0 |
| 446bj2z (446bj2) | exact-rank magnitude pass on down genes | 0 (bit) | 0 (bit) | **−0.002179** | 0 (bit) | 0 (bit) | 0 (bit) |

- The four bulk-only changes moved nmae by at most 0.0002 raw, fid by at most 0.00004 and jac by at most 0.000016. Reach
  moved by up to 0.0034 raw. Reach is the ordering channel, whose re-roll sd between digest-identical entries is about 0.027
  in scaled reach (0.0045 overall), and a bulk change re-draws which genes head the |z| order.
- The exact-rank DE change left all five other members bit-identical to 16 figures, fid, reach and jac included. That is the
  strongest form of C-62.
- The repair of zeroed genes is not a pure DE change: it edits the totals of the few genes it repairs, and the bulk moved by
  0.00007 raw.

*What the profile-error member reads once the bulk is free (518; C-84, C-85, C-93).*
- Each context is re-levelled to its cliff, the level (set from control cells) at which the member's jackknife credit
  exactly offsets the prediction's spread energy. There the member's raw value is 1.032 − 2X/ΣD. X is the energy-weighted
  agreement of predicted and true responses, and ΣD the truth's summed energy. The constant carries the mean row and the
  observed controls' noise (2.2–2.6 % of ΣD).
- The identity rebuilds H1's raw member exactly from its measured terms. On the board, solving it for X from 446bj2's raw
  value gives 2X/ΣD = 0.406. Spread energy is free and only direction pays: +1 % relative X is −0.0041 raw, +0.0007
  overall.
- At the same constants, the leaders' raw values imply +13.3 % … +13.8 % (two teams) and +33.6 % (the leader) more aggregate
  cosine than ours. That assumes their mean-row and control terms equal ours.
- The magnitude pot: at our directions, a per-target magnitude oracle raises X by +23.2 / +25.1 / +27.1 % on three held-out
  cancer lines (+52 % on an iPSC line, +26.5 % on H1). That is −0.10 … −0.13 raw, +0.018 … +0.023 overall, and pds-free by
  construction. Truly strong targets are under-scaled (median ideal-to-current scale 1.31 in H1's top decile by true energy).
- One transferred step came from that pot. A leave-one-donor-out (LODO) key weights a target by how strongly each donor's row agrees
  with the prediction built without that donor. It moved mse −0.014473 raw on the board (446bj2), past its preregistered
  −0.004 … −0.012, while pds moved +0.0002 against +0.001 … +0.004.
- After a raised cap on that key and a response-strength key (both priced on local folds only; the cap's board read is
  pending), about 20–23 % remains. Half of it sits in 11–19 targets per held-out line that no donor-side, CD4, VIPerturb, KOLF2, expression or
  cross-fold feature predicts out of fold.
- The direction pot is far larger: on H1 the per-target direction oracle at our energies takes the aggregate cosine from 0.113
  to 0.818. It is the information limit of §3(h)."""

S3L = """

### 3(l) What the rank and set members reward

**F20 — The rank member (reach) is limited by signs on weak responses, not by the order of the calls. The call set's size is
a bet on the reference's distribution of response sizes, and per-target information cannot hedge it across lines.
Provenance: candidate NEW (a property of the scorer at this pipeline's design), pending audit. Confirmation: analytic on H1
with exact known-answer rows (513, 521); a board dose pair on the global call-set size is preregistered (523) and unread.**

*Where reach lives (513; C-86).* Reach scores how deep a prediction's ranking of the reference's significant genes stays at
≥ 90 % sign purity. On H1, from the production emitter's saved statistics (the four DE members reproduced to 2 × 10⁻¹⁶):
- targets with fewer than 20 significant genes (32 of 82) score 0.643, with a walk sign accuracy of 0.81;
- 20–100 significant genes (9 targets) score 0.078 (walk accuracy 0.51); 100–1,000 (19) score 0.027; ≥ 1,000 (22) score
  0.014. An estimate from the board's members puts the challenge contexts' median target at 31–45 significant genes (192;
  disputed below), in the band where this pipeline gets 0.08.
- Sign accuracy decays down the pipeline's own order, from 0.80 in its first ten ranks to 0.55 beyond rank 250 (targets with
  ≥ 100 significant genes). The order carries confidence but never reaches the 0.9 the walk needs. Ordering by the truth's
  own |fold change| with our signs is worse (0.234).
- Demoting the high-expression block, whose signs (0.58 accurate in the 20–100 band) crowd the head, loses at every block
  size and scheme: Δreach −0.0007 … −0.038. Behind the ~2 per-target genes per target that are 0.95 accurate, the rest are at
  chance.
- No identifiable subset of useful size is pure. Genes on which ≥ 3 donors agree with |consensus| ≥ 1 log2 are 0.972
  accurate but number 72 over 82 targets, almost all in pools of ≥ 100.

So reach is sign-limited, not order-limited. Rivals' reach moved by +0.03 raw with their means pinned (fixed-mean rank
edits); such an edit needs a head several genes deep at ≥ 0.9 purity, which the public sources supply only for the strongest
targets.

*Sizing the call set per target (521; C-87).* A response-strength key (a noise-normalised count of strongly responding genes
in the donors' own rows) predicts H1's number of significant genes at r +0.60 (K562 alone) and +0.55 (four sources). The
earlier design predictor reached +0.24, and the key's accuracy does not depend on the truth's power.
- On held-out cancer lines it predicts only r +0.43 (K562), +0.25 … +0.27 (HCT116) and +0.19 … +0.25 (HEK293T); R² 0.02–0.18.
- The payoff is asymmetric. Shrinking a negligible target's call set to 10 genes gains +0.0004 overall on average; shrinking
  a non-negligible one loses −0.0034 (8.5 : 1). A shrink therefore needs 89–97 % of its targets to be negligible.
- The cancer lines put 7–46 % of their weakest-ranked 5–10 % of targets in the truth's top half, and every small-call-set
  policy prices −0.0005 … −0.016 on every fold.
- Per-target sizing stays closed, now by its payoff rather than by its input.

*The global size is a regime bet (521, 523; C-88).* fid rewards covering the reference's significant genes and jac penalises
over-calling, and their balance moves with the distribution of response sizes.
- On H1 a global call set of 1,600 genes instead of 2,600 is −0.0155 [−0.0251, −0.0072] overall at the reference's own
  power, and +0.0048 [+0.0009, +0.0100] with the reference thinned to 35 % of its significant genes.
- The board's regime is disputed: a fit to fid's coverage term says ~0.5, the jac level ~0.33. So the answer is a board
  measurement: a preregistered dose pair (1,600 and 2,000 genes called per target, DE-only, pds and mse pinned) is
  scheduled from 10-06 and unread.
- By-products. The board's weak end is thinner than proportional thinning assumes (nmae cohorts put 13–30 % of targets below
  10 gate genes, against H1's 32–42 %), so small-pool levers are over-priced on thinned H1. And "an unexpressed target has a
  null response" is refuted: panel targets below 5 CPM in a line's controls are null 26 % of the time on the K562 fold (11 of
  42), about as often as expressed ones (23 %).

*Why it matters for benchmark design.* Two of the four DE members turn on quantities a model cannot learn from public data.
Reach turns on the few significant genes of weak responses, where it needs ≥ 90 % sign purity. The fid/jac balance turns on
the reference's distribution of response sizes. A benchmark that wants its DE channel to reward biology rather than
calibration could report reach by response-size band and fix the call-set size relative to the reference's own counts."""

ANCHORS_C = [
    ("  (median K 4 → 6) moved nmae by +0.0013 (446b).", "after", S3K_ADD),
    ("within the published format contract and disclosed here, not as a biological result.", "after", S3L),
]
