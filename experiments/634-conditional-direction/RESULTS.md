# 634 — A DE-conditional direction prior from public banks LOSES to CLEAN4's unconditional sign on H1's gate. NEGATIVE.

2026-10-02 ~11:30Z, council session [d4cd5f] (experiment 600's breakthrough council). Analysis only: no build, no `vcc`,
no `cell-eval2`. Script `cond634.py`, log `work/cond634.log`, table `work/cond634.json`. Runtime 33 s, peak RSS 4.2 GB
for a few seconds (the Orion row subsample), slightly over the council's 4 GB announce line. Memory gate before the run:
~26 GiB free+inactive, no builder running.

## Question
CLEAN4's per-gene P_down (371) is the fraction of ALL perturbations in a bank with delta < 0. That is an UNCONDITIONAL
programme lean, dominated by non-responding (target, gene) pairs. nmae and reach score only the reference's DE genes, so
the sign they need is P(down | the gene responds). A board rival (chill) posted nmae 0.835 on its first entry with
"separate DE membership/direction ranking" and "signed conditional response estimates". Does a gene-level direction
prior estimated from DE EVENTS in public banks predict gate-gene signs better than CLEAN4?

## Construction
- Banks (never H1): k562_gwps, orion_hct116, orion_hek293t, rpe1, nadig_jurkat, k562_essential, viperturb_k562.
- Event in a bank: robust z = (delta − median_g) / (1.4826·MAD_g) over the bank's perturbations (centre/scale from a
  3,000-row subsample, seed 634); |z| ≥ 3, 4 or 6; perturbations all, or non-essential only.
- P_cond(g) = (down events + 0.5) / (events + 1); family-balanced (equal weight per bank where ≥ minev events) or pooled.
- Truth: h1-test reference (323), 85 targets × 10,779 gated genes, gate = q < 0.05, finite non-zero lfc, gene ≠ target
  (65,045 gate pairs; up-fraction 0.475 = 1 − all-down's accuracy).

## Result (gate-gene sign accuracy; "acc_w" weights by |lfc|, the nmae-relevant reading)

| source | targets | coverage | acc | acc_w | acc on truth-up | acc on truth-down |
|---|---|---:|---:|---:|---:|---:|
| all-down | all 85 | 1.000 | 0.525 | 0.542 | 0 | 1 |
| **CLEAN4 (unconditional, 4 banks)** | all 85 | 0.531 | **0.573** | **0.572** | 0.483 | 0.660 |
| unconditional mean, 7 banks | all 85 | 0.997 | 0.555 | 0.557 | 0.358 | 0.734 |
| COND family-bal, all perts, z4 | all 85 | 0.971 | 0.519 | 0.494 | 0.455 | 0.578 |
| COND pooled, all perts, z3 (best COND) | all 85 | 0.971 | 0.521 | 0.494 | 0.392 | 0.638 |
| COND, non-essential perts, z4 | all 85 | 0.968 | 0.503 | 0.480 | 0.295 | 0.693 |
| **head to head, same genes: CLEAN4** | all 85 | 0.531 | **0.573** | **0.572** | | |
| **head to head, same genes: COND z4** | all 85 | 0.531 | 0.541 | 0.534 | | |
| COND where it disagrees with CLEAN4 | all 85 | 0.145 | **0.441** | 0.435 | | |
| head to head: CLEAN4 | 67 non-essential | 0.504 | 0.558 | 0.555 | | |
| head to head: COND z4 | 67 non-essential | 0.504 | 0.519 | 0.515 | | |

Every conditional variant (z 3/4/6, min events 3/10, family-balanced or pooled, all or non-essential perturbations,
confidence ≥ 0.1) loses to CLEAN4. Its |lfc|-weighted accuracy is BELOW 0.5 (0.47–0.50), so on H1's large gate effects a
gene's typical direction-when-responding in other lines points the wrong way more often than not. Where the two
sources disagree, the conditional sign is right 44 % of the time.

## Reading
- Gate-gene signs on an unseen line follow the line's broad UNCONDITIONAL lean (CLEAN4's programme) better than the
  gene's DE-event direction elsewhere. Responses that cross a significance bar in other lines are line-specific in
  direction, while the weak shared lean transfers.
- So the rivals' "conditional response estimates" are not a gene-level direction table. If chill's edge is
  conditional direction, it is per TARGET (P(up | DE, target, gene) from same-target source rows), which 516 already
  priced at ≈ 0 on cancer folds for public sources. Its edge is more likely the membership model + count generation.
- Closed at this construction: gene-level DE-conditional direction priors from public banks (z-threshold events). Not
  tested: per-target conditional signs with a membership model, and a Flex-only bank (VIPerturb alone).
- Side note: the unconditional mean over 7 banks (adding Orion, RPE1, VIPerturb) is WORSE than CLEAN4's 4 banks on H1
  (0.555 vs 0.573 at much higher coverage), in line with 371/430: more banks do not make a better programme.
