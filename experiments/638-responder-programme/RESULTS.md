# 638 — Responder-weighted programmes do not beat CLEAN4 on gate signs (NEGATIVE); CLEAN4's own definition on the full gene axis adds coverage (tiny value)

> **VOID ROWS (2026-10-02T19:10Z; [21dd4c]'s 616 audit, 7dac7a98, `experiments/616-clean-orion-folds/AUDIT-orion-dependence.md`):** the thinned Orion gate tables' DE truth was ~98 % selection artefact (HEK293T: 37,748 calls on the thinned table vs 774 on clean raw cells; 13 / 132 targets scoreable). In this experiment, the hct116 / hek293t rows are VOID. What remains: h1 / k562, still negative.

2026-10-02 ~12:05Z, council session [d4cd5f], portfolio item 2 (candidate sign sources on 636's leave-one-line-out
instrument). Analysis only. `resp638.py` 17 s, peak RSS 6.3 GB (mostly mmap pages of the Orion deltas; over the
council's 4 GB announce line, and not announced, so I'm noting it here; later heavy runs go through
`experiments/602-council-gate/gate602.sh`). `ext638.py` 1 s, 0.5 GB.

## A. Responder-weighted programmes (`resp638.py`, `work/resp638.json`)
Gate genes are responders, but CLEAN4 counts every perturbation equally. Per bank, with s_p = RMS of perturbation p's
delta (no truth used), P_down was weighted by s_p (W1) or s_p² (W2), restricted to the top 10 % / 25 % of
perturbations by s_p (TOP10/TOP25), or taken as the sign of the mean delta (MEAN). Each form was pooled over (i) all
banks allowed on the surface and (ii) CLEAN4's own banks.

Δ gate sign accuracy vs CLEAN4-LOLO on shared genes, all targets (95 % target bootstrap):

| form | h1 | hct116 | hek293t | k562 |
|---|---|---|---|---|
| U, allowed banks | −0.007 | −0.032 | −0.064 | −0.032 |
| W2, CLEAN4's banks (best overall) | −0.001 | +0.007 [−0.002, +0.016] | +0.004 [−0.002, +0.011] | −0.023 |
| TOP10, allowed banks | +0.004 [−0.005, +0.012] | +0.002 | −0.028 | −0.020 |
| MEAN, CLEAN4's banks | −0.003 | −0.020 | −0.039 | −0.014 |

No form beats CLEAN4 by the portfolio's bar (≥ +0.01 on most folds). Pooling more banks again loses everywhere.
**Closed: perturbation-strength weighting of a programme** (with 634's gene-level DE-event weighting, both ways to
make the programme "responder-aware" from public banks lose or tie).

## B. CLEAN4-ext: the same definition on each bank's full gene axis (`ext638.py`, `work/ext638.json`)
371 built CLEAN4 on the 7-bank gene intersection (6,346 genes), so it covers only 53 % of h1's gate. Recomputed from
the banks' own axes (k562_gwps, k562_essential, nadig_jurkat; nadig_hepg2 has no delta on disk; on k562 only Jurkat is
allowed), it covers more gate genes, with signs better than all-down on every surface:

| surface | τ | added gate coverage | CLEAN4-ext acc on added genes | all-down on them |
|---|---|---:|---:|---:|
| h1 | 0 / 0.05 | +19.4 % / +8.1 % | 0.539 / 0.559 | 0.520 / 0.513 |
| hct116 | 0 / 0.05 | +11.5 % / +5.1 % | 0.521 / 0.524 (w 0.609) | 0.487 / 0.536 |
| hek293t | 0 / 0.05 | +12.0 % / +4.7 % | 0.658 / **0.808** | 0.313 / 0.327 |
| k562 | 0 / 0.05 | +27.0 % / +16.9 % | 0.579 / 0.596 | 0.477 / 0.479 |

**Value.** The uncalled-gene hook delivers MDE-sized magnitudes, about 4 % of the (2a−1)·|l| pot on its coverage (H1:
−0.0015 raw nmae from a 0.036 share). The extension adds about a quarter to the hook's share, so it is worth about
+0.0001 overall on the board. That is too small for an arm. It is a free coverage fix to fold into any build that
consumes CLEAN4 (the hook, the block prior), with a KA that leaves CLEAN4's existing genes unchanged.
