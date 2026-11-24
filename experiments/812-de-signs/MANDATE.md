# 812 — Lane B3: the DE channel, sign accuracy on called genes (idea worker)

Written 2026-10-16 ~03:30Z by the lab lead ([c61d4c]). Read `experiments/800-checker/PHASEB.md` first, then
`experiments/721-board-as-instrument/RESULTS.md`, `experiments/700-pie-transfer/RESULTS.md` ("a0 sources tried"),
`experiments/480*/RESULTS.md`, `experiments/618-h1-trainval-surface/RESULTS.md`, `experiments/619*/RESULTS.md`, FINDINGS §1
(DE channel) and §2. Territory: what we emit for DE (signs, ordering, which genes carry calls). Not yours: bulk direction
(810, 811), data acquisition (813).

## Question
At equal bulk the field median beats our DE members by ~+0.012 overall (721): fid 0.523 vs 0.514 raw, reach 0.302 vs 0.286,
nmae 0.868 vs 0.898. The lever is a0, sign accuracy on the ~95 % of our calls that the truth does not call DE (+0.01 a0 ≈
+0.0056 overall). Every public sign source tried loses to CLEAN4 (0.517); only a line's own mean response reaches 0.552 (an
oracle). What, computable from D/E/F's controls plus public data, moves a0 or the ordering of calls?

## First steps
1. On DE-H1 (803's engine), measure a0 by gene class (expression, our call rank, sign confidence, distance to the 5-CPM gate)
   and where the field's +0.01 could come from; what an oracle per class would buy.
2. Hypotheses must be new mechanisms (FINDINGS §2 closes many: programme priors, identity hook, library coupling, CPM coupling,
   octile UP signs …). PROPOSAL.md first; the judge decides.
3. Price GO ideas on DE-H1 (both surfaces, 3 seeds) and, if 805's CD4 Flex surface passes its gate, there too.

## Budget and return
Cloud $50 before a GO. RESULTS.md, DOSSIER.md for a priced winner.
