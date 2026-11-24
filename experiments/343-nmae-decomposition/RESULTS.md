# 343 — The H1 NMAE deficit is a sign problem at the design stage, not shrinkage and not emission

September 9, 2026, 22:54–22:56 UTC. **DESCRIPTIVE ONLY**: arithmetic on 323's saved cap15 outputs
(`stats-cap15.npz`, `reference_de.parquet`, `effects.npz::cap15`). No emission, no scorer call, no gate,
no new data. Scripts `decompose343.py`, `signs343.py`; outputs in `work/`. The recomputed mean NMAE over
the 58 eligible targets is **0.950836**, matching 323's recorded value exactly.

## Why this was asked
On the 2026-09-09 21:10Z board the 0.0444 gap to the leader is NMAE +0.023 and reach +0.019 (scaled
contributions). Across the top four teams we are last on NMAE (0.108 scaled vs 0.170–0.247) while leading
MSE and Jaccard. NMAE_p = Σ|l̂−l| / Σ|l| over the reference's own significant set, so it can be lost
three ways: predicting nothing on those genes, predicting too small, or predicting the wrong sign.

## Result, cap15 on the reused 85-target H1 panel (58 NMAE-eligible targets, mean gate 1,120 genes)

| quantity | value |
|---|---:|
| NMAE, observed | 0.9508 |
| magnitude ratio m = Σ|l̂| / Σ|l| on the gate | 0.205 |
| gate mass predicted exactly zero | 0.0000 |
| gate mass inside our own significant set | 0.346 |
| |l|-weighted sign accuracy, delivered lfc, all gate genes | **0.588** |
| … genes we call significant / genes we do not | 0.638 / 0.552 |
| |l|-weighted sign accuracy of the **designed** effect (pre-emission) | **0.600** (0.649 / 0.558) |
| weighted agreement designed-sign vs delivered-sign | 0.867 |

Exact per-target counterfactuals, then averaged:

| counterfactual | NMAE |
|---|---:|
| our magnitudes, every sign correct | **0.803** |
| our signs, best single scalar per target (grid to 6×; mean best 1.48×) | 0.928 |
| our signs, least-squares scalar (mean 3.1×) | 1.005 |
| every sign correct and best scalar | 0.658 |

## Reading
1. **Not a coverage problem.** Nothing on the gate is predicted as zero; the emitter's keep-top-8000 leaves
   every gate gene with a nonzero fold change.
2. **Not fixable by amplitude.** At 59% weighted sign accuracy the NMAE-optimal gain is 1.5× for −0.023
   raw, and the least-squares gain makes it worse. 317 already measured on the board that amplitude costs
   more MSE than it buys; this shows why: the wrong-sign 41% of gate mass pays 1+a.
3. **Not the emission.** Delivered signs agree with designed signs on 87% of gate mass, and the design's own
   weighted sign accuracy is only 0.600. Emission costs about one point; the other forty are upstream.
4. **It is the transfer information.** Sign accuracy on genes the reference actually moves is 60% at the
   design stage on a held-out line. That is the quantity 332–341 tried to raise via source-stage models and
   could not (none beat the per-gene diagonal on retrieval; 341 raised weighted sign +0.033 only with the
   recipient's own training labels). Reach is the same quantity read as a 90%-pure prefix.

## Limits
Reused H1 development panel, one line, 58 targets; not an A/B/C forecast. Weighted sign accuracy over
the gate is not the official fidelity (which reads our call set) or reach (which reads the reference pool
in our confidence order). No claim about which mechanism the leader uses.

## Addendum 23:05Z — the record already contains these conclusions; 343 confirms them on the 319 composition
- **227** (`experiments/227-nmae-axis/RESULTS.md:151-183`) stratified the exact NMAE numerator and denominator by
  control CPM: every stratum's numerator share tracks its denominator share within 1%, an oracle per-CPM-stratum
  scalar is worth only +0.0089 raw, and sign correction is worth +0.10 to +0.24. Same conclusion as §"Reading" above.
- **317** showed every emission knob (overdispersion, keep/tail spread, keep-top) is digest-identical on delivered
  counts and amplitude trades 28–41:1 against on the board.
- **094/155** measured the direction asymmetry: predicting up was right 10.8% on 094's panel; the reference pool is
  59.7% down and a constant all-down vector beat the incumbent's sign accuracy.
- New here, descriptive only: on genes with recipient basal CPM in the 5–16 band the emitter delivers a designed
  increase as a decrease 40% of the time (35% even for top-decile designed effects) versus 0–7% in every other band,
  with 6.7% of the reference's significant calls and 15% of NMAE mass in the sub-10-CPM band where our delivered
  sign accuracy is 0.508 and magnitude ratio 0.033. 342 (Hamilton totals) moved NMAE by −0.0006 on this composition,
  and 227's stratum ceiling bounds any fix there at under +0.009 raw. Not an arm.
