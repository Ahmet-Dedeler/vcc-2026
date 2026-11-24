# 636 — The council's leave-one-LINE-out instrument for DE sign sources (built, KA passes) + its reference rows

> **VOID ROWS (2026-10-02T19:10Z; [21dd4c]'s 616 audit, 7dac7a98, `experiments/616-clean-orion-folds/AUDIT-orion-dependence.md`):** the thinned Orion gate tables' DE truth was ~98 % selection artefact (HEK293T: 37,748 calls on the thinned table vs 774 on clean raw cells; 13 / 132 targets scoreable). In this experiment, the hct116 / hek293t surfaces are VOID. What remains: h1, k562 (and Broad-16 where used).

2026-10-02 ~11:45Z, council session [d4cd5f] (portfolio item 2 in `experiments/600-breakthrough-council/SYNTHESIS.md`).
Analysis only; 7.6 s, 3.2 GB peak. Script `inst636.py` (importable `Inst` class plus a CLI); reference rows in
`work/ref636.json`.

## What it is
It scores any per-(target, gene) prediction of a held-out line's DE signs (P_down) and, optionally, DE membership on
four surfaces whose gates were built with the scorer's own rules (480's `gates480.py`, KA'd in 480 §1).

| surface | file (480/work) | targets | gate pairs | notes |
|---|---|---:|---:|---|
| h1 | gate-h1_full_ka.npz | 85 (67 non-essential) | 65,044 | Arc Flex, full power (= 323's h1-test reference) |
| hct116 | gate-orion_hct116.npz | 132 | 11,309 | Orion 3′, thinned to 150 cells |
| hek293t | gate-orion_hek293t.npz | 132 | 37,748 | Orion 3′, thinned |
| k562 | gate-k562_gwps.npz | 267 | 3,179 | low power |

Gate signs are the target because 635 showed the shared-control-noise objection does not touch them, and that the
reachable programme prize lives there. Leave-one-line-out is enforced per surface (`BANK_FORBIDDEN`): never the
surface's own line, never a same-study sibling (430), never an H1 bank on h1.

Outputs per source: gate sign accuracy (plain, |lfc|-weighted, on truth-up and truth-down genes, all and non-essential
targets), coverage, the paired Δ vs a reference source on shared genes with a 2,000-draw target bootstrap, and for a
membership score, per-target AUROC (gate vs other tested genes) and precision at n_real.

## Known answer
h1 CLEAN4 gate sign accuracy 0.5734, equal to 634/635's (same banks, same gate). PASS.

## Reference rows (all targets; Δ vs CLEAN4-LOLO on shared genes, 95 % target-bootstrap CI)

| surface | all-down | CLEAN4 LOLO (cov) | best single unconditional bank | unconditional mean of allowed banks |
|---|---:|---:|---|---|
| h1 | 0.525 | **0.573** (0.53) | nadig_jurkat −0.014 / k562_essential −0.008 | −0.007 [−0.016, +0.003] |
| hct116 | 0.598 | **0.691** (0.78) | nadig_jurkat +0.006 [−0.002, +0.015] | −0.032 [−0.040, −0.023] |
| hek293t | 0.562 | **0.663** (0.74) | **nadig_jurkat +0.057 [+0.053, +0.061]** | −0.064 [−0.070, −0.059] |
| k562 | 0.456 | **0.603** (0.55; Jurkat + HepG2 only) | nadig_jurkat −0.020 [−0.040, 0.000] | −0.032 [−0.061, −0.001] |

## Reading
- CLEAN4 is a strong bar on the cancer folds (0.60–0.69), and much stronger than on H1 (0.57). Pooling more banks
  LOSES on every surface: the equal-weight unconditional mean of all allowed banks is −0.006 … −0.064. Programme sources
  must be weighted by quality, not count (consistent with 371 and 430).
- Nadig Jurkat alone is the best single source on the Orion folds, decisively on HEK293T, but not on K562 or H1. That
  is a per-surface pattern to watch in any learned weighting, not an adoption.
- The instrument is ready for [7772dc]'s supervised DE model (621). Format: `work/pred-<surface>.npz` with targets,
  genes, P_down (T × G, NaN = no call), optional member, plus a json sidecar listing the excluded banks.

## CAVEAT for every reading on this instrument (13:50Z; [4d0207]'s 609 and SCIENCE §5)
The Orion folds' gate genes are down-skewed on highly expressed genes: CPM-block gate genes are 72–74 % DOWN on HCT116 /
HEK293T, vs 0.56 on H1 and 0.51 on K562 (uncalled genes 25–28 %). A comparison between sign sources whose down-rates
differ on highly expressed genes is therefore tilted on Orion toward down-heavy priors (CLEAN4 is ~84 % down) and toward
anything correlated with expression. 639d's "block route" was this artifact (withdrawn: 639f).
**Rule from now on:**
- verdicts lead with H1 (and Broad-16 where it applies);
- Orion rows are shown only beside an expression-only control (zc + 0.5·z(log1p own CPM)) and the source with
  expression partialled out;
- report each fold's gate down-fraction by region (head / block / uncalled) next to any regional reading;
- never read an absolute magnitude off Orion (AGENTS 6.5).
