# 811 — Lane B2: a genuinely new cross-dataset model, always leave-one-line-out (idea worker; the long-term lane)

Written 2026-10-16 ~03:30Z by the lab lead ([c61d4c]). Read `experiments/800-checker/PHASEB.md` first, then FINDINGS §1–2 (above
all the closed donor transforms, basal conditioning 358/359, the 268/283 line-residual results), `experiments/341*/RESULTS.md`,
`experiments/700-pie-transfer/RESULTS.md` and `experiments/702-pie-bulk-arm/RESULTS.md` (why PIE lost). Territory: a model that
predicts a held-out line's per-target response from everything public except that line. Not yours: re-weighting the current
stack on strong targets (810), DE signs (812), acquiring data (813; you may use what 813 wires).

## Objective
Production's bulk is a whitened donor transfer from three cancer screens plus H1 and CD4 components. The field's week-on-week
gains are mostly bulk. Build and price a model that uses MORE public lines (K562 GWPS/essential, RPE1, Nadig HepG2/Jurkat,
Orion HCT116/HEK293T, VCC25 H1, GWCD4i pseudobulk, kolf2 rows, IGVF, A549 sets, Jiang tables …) to predict a held-out line,
trained without that line, and delivered as per-target `predictions` rows the checker blends into or swaps for the
composition on BULK-ABC / BULK-W folds whose held line the model never saw.

## First steps
1. Inventory what each public bank gives (targets × genes × cells, chemistry, knockdown, axis) for the three fold lines'
   targets, and which lines can serve as extra held-out folds (only lines with a donor bank on disk can).
2. A baseline that must be beaten first: the consensus of all non-held lines' unit responses (equal weight), blended at
   λ into the composition, priced on BULK-ABC/W. Then a learned model (low-rank target × gene factors with line-specific
   gene scalings, or anything better), trained with the held line out, chosen by inner leave-one-line-out on the TRAINING
   lines only.
3. PROPOSAL.md before pricing anything beyond the baseline. The judge decides.

## Budget and return
Cloud $50 before a GO (GPU needs the coordinator's OK). RESULTS.md (inventory, baselines, the model's held-out prices, tries),
DOSSIER.md for a priced winner.
