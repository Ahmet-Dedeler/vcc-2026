# 721 — The public board as the instrument for the DE channel (2026-10-15 02:40–04:10Z, session [020377])

Class: measurement (no build, no `vcc`, no `cell-eval2` scoring run; one in-process call of cell-eval2's metric functions on a
synthetic table). Preregistered in `PREREGISTRATION.md` before any fit. Inputs: `data/leaderboard/*.json` (252 snapshots) plus
the 10-14 22:14Z capture → 18,634 unique scored entries from 1,841 teams (`work/entries_all.json`); the installed cell-eval2
0.16.0 metric code; the spec's constants; our three known-K board entries; H1-test's real DE table
(`experiments/170-gamma-mse/work/de_real_h1.npz`, 85 targets × 10,779 genes). No FINDINGS number was an input.

## 0. Answer
**Our DE channel is dominated at equal bulk, and the lever is sign quality on the genes we call, not the call-set size.**
Among the 65 teams whose raw pds is within ±0.01 of ours, the median team scores fid 0.5232 / reach 0.302 / nmae 0.868 raw
against our 0.5144 / 0.286 / 0.898; only jac is ours (0.0408 vs 0.0337). Converted through the board's own scales, the field
median at our bulk is **+0.0124 overall** above us on the four DE members (+0.0225 at the 75th percentile), i.e. 0.275 → 0.287–
0.297, the 0.30 pack. Our DE members have not moved since 09-22 (fid sd 0.0017 over 22 entries; every gain was bulk).
A K arm is NOT the move: the K series is concave (+0.0035 then +0.0012 overall per +1,000 genes) and the forward model shows
fidelity is insensitive to K once a perturbation's yield saturates. The pack's fidelity is sign accuracy on called genes
(the ~95 % of our 2,600 calls that are not truly DE land at ≈ 0.50–0.52 against the real sub-threshold direction; the pack
lands ≈ 0.53, Illumina 0.56), the same quantity that drives reach, which is why fid and reach move together within teams.

## 1. The forward model, transcribed from the code and verified against it
From `cell_eval2/metrics/direction.py` (`_direction_frame`, `_components`, `_fidelity_family`) and `de.py` (`de_sig_jaccard`):
- called set per perturbation: `{p_adj_pred < 0.05}` ∩ {real lfc defined}, own target gene removed → `n_pred`;
- `k` = called genes whose predicted sign equals the REAL sign (any gene, significant or not);
- **scored fidelity** `de_wilcoxon_direction_fidelity_yield_raw = k / max(n_pred, n_conf)`, `n_conf` = real-significant count;
- `jac = |called ∩ real-sig| / |called ∪ real-sig|`; `reach_raw = k*/n_conf`, k* = deepest prefix of the PREDICTION's ranking
  (p_adj, p, |lfc|) over REAL-significant genes whose cumulative sign purity ≥ 0.90.
**Known answer:** a synthetic 3-target table with n_conf 10/30/50, n_pred 25, half the calls right-signed:
`de_direction_fidelity_yield_raw` returned 0.52 / 0.4333 / 0.26 and `de_sig_jaccard` 0.40 / 0.8333 / 0.50, equal to the
formulas (`work/` log; the readers warned only about the absent p_value column). Transcription verified.

Consequence (the dilution identity): with TP true-DE calls at sign accuracy a1 and K−TP non-DE calls at accuracy a0 (the real
direction of a sub-threshold gene is defined for 100 % of them on H1),
`fid_p = [a1·TP + a0·(K−TP)] / max(K, n_conf)`. At our precision (TP/K ≈ 0.04–0.05) the second term is ~95 % of the mass:
fid ≈ a0 + 0.045·(a1 − a0). Moving fid by +0.01 needs a0 up by ≈ +0.009, or a1 up by ≈ +0.2.

## 2. Scaling constants recovered from the board (all 18,634 entries; residual = clamps only)
fid: b 0.51246, r 0.80993 (spec 0.505–0.522 / 0.795–0.832 per context) · jac: b 0.03039, r 0.39766 (spec 0.021–0.037 /
0.375–0.423) · reach: b 0.07898, r 0.96522 (spec 0.047–0.097 / 0.958–0.978) · pds: b 0.50047, r 0.95127. Exchange rates:
+0.01 raw fid = +0.0336 scaled = +0.0056 overall; +0.01 raw jac = +0.0272 scaled = +0.0045; +0.01 raw reach = +0.0113 scaled =
+0.0019; −0.01 raw nmae = +0.0151 scaled = +0.0025.

## 3. Hypotheses, read against the preregistered rules
- **H1 dominance — TRUE.** 27 latest entries within ±0.02 pds / ±0.03 mse (scaled) of our best: median fid 0.0379 scaled
  (ours 0.0056), 89 % above ours; reach 78 % above; nmae 81 % above; jac 0 % above; overall median 0.2879 vs 0.2749, 74 %
  above. (`work/descriptive-frontier.txt`.) Only one entry on the whole board dominates us on all six members (Illumina).
- **H2 regime — FALSE** (the prereg's TRUE reading needed K ≥ 4,000). (i) Our own dose curve: n_pred 1,600 → 2,000 → 2,600
  gave overall +0.00348 then +0.00119 per +1,000 genes; fid +0.0059 then +0.0028; jac +0.0009 then −0.0003; reach flat
  (`work/kseries.txt`). The marginal genes' sign accuracy was 0.52 both times. (ii) The yield term is small: H1-test's
  n_conf is median 89, mean 765, p90 2,691; mean yield at K 1,600 / 2,000 / 2,600 is 0.925 / 0.947 / 0.972
  (`work/h1-nconf.txt`), and A/B/C's n_conf is ≈ 0.4× H1-test's (b_jac 0.021–0.037 vs 0.071), so beyond 2,600 the yield left
  to buy is ≤ 1–2 % of fid. (iii) The field: fid rises with jac across teams (0.0305 → 0.0345 → 0.0375 → 0.062 median jac
  by fid bin), not the inverse that a pure K-up would produce. A naive linear extrapolation of our curve (+0.0007 at K 3,200,
  UNVALIDATED, concave) is the most a K arm could claim; it is not worth a slot.
- **H3 consistency — TRUE in sign, size consistent.** The H1 n_conf tail predicts that K 1,600 and 2,000 lose yield on the
  15–18 % of targets with n_conf above K; both board reads lost fid (−0.0041, −0.0017 raw) and left jac within 0.0004. Sizes
  are below H1's yield loss × fid (≈ −0.012 / −0.006 for an H1-shaped tail), as expected for a tail ~0.4× H1's.
- **H4 our accuracy — not separable from board data alone.** The dilution identity gives one equation in (a0, a1); with
  a0 = 0.50, a1 ≈ 0.83; with a0 = 0.51, a1 ≈ 0.60. The repo's H1 instrument puts the called block's precision at 0.534 vs
  CLEAN4 0.522 and "sign accuracy decays 0.80 → 0.55 down our order" (FINDINGS), consistent with a0 ≈ 0.50–0.52. The
  measurement that settles it is on H1 with the production emitter: sign accuracy on called genes split by real significance
  (handed to the science lane, §6).

## 4. What the field says fidelity IS
- Across the 339 latest entries above 0.2: Spearman fid–reach +0.57, fid–mse +0.50, fid–nmae +0.47, fid–pds +0.37, fid–jac +0.32
  (`work/fid-vs-pds.txt`). fid_raw ≈ 0.3545 + 0.1929·pds (R² 0.11): our residual −0.0074 (−0.8 sd); with jac −0.0092 (−1.0 sd).
- Within teams (238 teams with ≥ 10 entries), consecutive changes in fid co-move with reach (median Spearman +0.42, positive in
  87 % of teams) and jac (+0.37, 82 %), weakly with nmae (+0.12) and pds (+0.16) (`work/comovement-and-opportunity.txt`).
  So fid and reach are one underlying quantity, the accuracy and ordering of predicted signs, and bulk moves are mostly
  separate. Illumina's 10-14 jump moved mse 0.305 → 0.241 raw with fid +0.013, reach +0.037, nmae −0.033 together: a model
  change, not an emission change.

## 5. Decision implied for the final recipe
1. No call-set-size arm. `rc_callset_cpm` stays off (K 2,600 total); 523's conclusion stands, now for the right reason.
2. The DE channel is worth +0.012 overall at the field median for our bulk and +0.022 at p75; the raw targets are
   fid 0.523, reach 0.302, nmae 0.868 (jac may fall to 0.034 without net loss). Every candidate must be priced on these three
   raw members, and "sign accuracy on called non-DE genes" (a0) is the quantity to move: +0.01 of a0 ≈ +0.0056 overall by itself.
3. Candidates that act on a0 and ordering, in the order the board evidence supports: a sign source for the called block that
   is right on sub-threshold genes (PIE's `lfc_pred`/`p_de` sign for the block, 700; the bulk delta's own sign where
   |delta| is confident), and an ordering of calls by sign confidence (reach). Each is a prior-only DE arm (re-roll ≈ 1e-5),
   so one board read each is enough, before Oct 21.

## 6. Validation and what would show this is right, other than us
- The transcription reproduces cell-eval2 on a synthetic truth (§1). The recovered b/r sit inside the spec's per-context
  ranges (§2). The H1 n_conf tail predicts the sign of the two held-out K reads (§3 H3).
- **Replication:** the coordinator assigned [bf003b] to re-derive §3–§4 from `work/entries_all.json` with its own code
  (the `.txt` files carry the exact numbers to match: 27 / 89 % / 0.0379; 65 teams, median fid 0.5232; +0.0124). A
  Cursor cloud agent (a different model) was attempted at 04:20Z and refused: the account is on the free plan now
  (`HTTP 403 "Cloud Agent is not available for free users"`); the prompt is kept in this session's log if Ahmet re-upgrades.
- **The measurement that decides §5.3 (science lane):** on H1 with the production emitter's real DE table, compute per
  target (a) sign accuracy on called ∩ real-sig, (b) sign accuracy on called ∩ non-sig, (c) the same for PIE's signs on the
  same called set. If (b) < 0.52 and PIE's (b) ≥ 0.54, the block sign source is the arm.
- The final check is the Nov board: nothing here is adopted on this analysis alone (rule 11.1c).

## 7. Falsifiers that did not fire
A single `n_conf` distribution could reproduce the fid–jac–reach pattern of the field only through a0/a1 differences, so the
"instrument void" clause (prereg) did not trigger. The one thing the board cannot resolve is a0 vs a1 separately (§3 H4).

## 8. Addendum 2026-10-15 ~05:30Z: the decision measurement came back ([bf003b], `experiments/700-pie-transfer/signsplit700.py`,
`work/signsplit700.json`; h1-test, production emitter 435 s0 call set, adjudicable genes)
| sign source on OUR called set | (b) called & non-sig, n = 195,858 | (a) called & sig, n = 25,061 |
|---|---:|---:|
| emitted (production) | **0.5154** | **0.6007** |
| CLEAN4 prior | 0.5171 | 0.5778 |
| design key | 0.5076 | 0.5667 |
| PIE `lfc_pred` | 0.5115 | 0.5594 |
| PIE (p_de ≥ 0.5) else CLEAN4 | 0.5140 | — |
- **Cross-check of the forward model:** the board-only dilution identity (§3 H4) gave, at a0 = 0.51, a1 ≈ 0.60. H1's direct
  measurement is a0 0.5154 / a1 0.6007. Two instruments, one decomposition.
- **Decision (prereg §6 rule):** PIE's (b) 0.5115 < 0.54, so the PIE block-sign arm is OFF. CLEAN4 remains the best public
  sign source for sub-threshold directions. The pack's a0 ≈ 0.53 is not reachable with any public consensus tested; it is
  per-target model skill on small effects (what also moves their reach and nmae), which is the transfer problem itself.
- What the board still leaves open for the DE channel at equal bulk: the pack's nmae 0.868 vs our 0.898 (magnitude on called
  genes, FINDINGS lead 631) and reach 0.302 vs 0.286 (sign-confidence ordering). Neither is a K or a sign-source question.

## 9. Replication (06:05Z): §2 and the DE-opportunity table reproduced by [bf003b] with independent code
`experiments/700-pie-transfer/replicate721.py` → `work/r721/replicate721.txt` (4aeb2aa2). Scaling constants exact to five
decimals for fid, jac, reach and pds (rms residuals 3e-4 … 4e-3); the pooled linear fits for nmae and mse are poor, as §2
says (per-context scaling plus clamps). DE opportunity at equal pds, n = 65, reference 498p2g: fid +0.0051, reach +0.0029,
nmae +0.0077, jac −0.0033 → **+0.0123 overall** (here +0.0124); the only difference is the jac median, 0.0085 vs 0.0091, a
latest-entry tie. §3–§4 hold. With §8's H1 cross-check this closes the preregistered validation list except the final board.
