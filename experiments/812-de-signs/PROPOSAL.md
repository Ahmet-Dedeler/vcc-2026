# 812 PROPOSAL: the DE channel's signs, order and call composition. Verdict: NO-GO on all four (for the judge)

Worker 812, 2026-10-16 08:10Z. Territory: what we emit for DE (signs, ordering, which genes carry calls). The numbers
behind every line are in `RESULTS.md`. Every screen below was analytic on DE-H1's realized stats, and each carried a
known-answer row (fid / reach / nmae recomputed equal the arm's recorded raw values to ≤ 1e-12, or the scorer's lfc
recomputed from cells to 1e-5).

**Honesty about order.** The mandate asks for the proposal before compute. These four screens were run as exploratory
step-1 measurements (about $0.40 of EC2 plus 1 to 13 s gated jobs on the Mac), before this file existed. None of them asks
for a GO, so nothing is preregistered here. Their try counts are listed so that any later reuse pays its best-of-N.

| # | hypothesis (class) | mechanism | analytic H1 price, overall-eq (f 1) | tries | verdict |
|---|---|---|---|---|---|
| H-NTCB | NTC-guide control bias (scorer mechanics; day-of from D/E/F's own `ntc_id`) | a guide's own off-target effect shifts the scorer's control mean by its share; every target's truth leans against it; a Huber centre of the guide means isolates it | +0.00017 (h1test) / +0.00016 (h1tv): 5–6 flagged genes, M0 sign accuracy 0.69 / 0.66 vs CLEAN4 0.43 / 0.44 on them | 4 estimators × 3 thresholds | NO-GO: real but its footprint is 0–5 genes per context (A 5, B 1, C 1 at \|z\| ≥ 2) |
| H-CBLOCK | block composed by CLEAN4 sign confidence instead of CPM (scorer mechanics) | swap low-confidence block genes for confident lower-CPM genes, so a0 rises | fid +0.0034 … +0.0090 raw, jac −0.0037 … −0.029: net −0.0091 … +0.0004 | 9 (3 floors × 3 swap shares) | NO-GO: the CPM block's genes are more often truly DE, so jac pays for fid |
| H-CISFRONT | expressed cis neighbours (< 1 / 5 kb TSS) at the front of the call order with a DOWN sign (transfer, needs a GENCODE table) | dCas9-KRAB spreading; reach reads the front | +0.00003 (h1test) / +0.0004 (h1tv) | 2 distances × 3 seeds | NO-GO: 425 already carries them, and h1test's 12 real-significant pairs were all right-signed already |
| H-HEAD | head size (key-ranked calls) other than 400 (scorer mechanics) | a key-ranked head should buy precision × a1 | H 200 / 800 / 1,200: −0.0004 / −0.0004 … −0.0009 / −0.0018 … −0.0020 | 7 sizes | NO-GO: H 400 is the optimum, and the key-ranked head has LOWER precision than the CPM block |

**What I recommend instead (no compute asked).** The emission side of the DE channel is at its optimum for the sign
information we have. `rates812` puts the field's DE gap in sign-skill units (RESULTS §4): the field median's +0.012
overall is reproduced by about +0.04 … +0.05 of sign accuracy on the truth's DE genes (called a1 and uncalled u1) plus
+0.005 on a0. That is per-target model skill on responding genes, the same size as 635's clean-oracle edge of a line's
own responding programme (0.618 vs CLEAN4 0.573). It belongs to lane 811 (a new model). I propose that 811 adopt
`rates812`'s classes as its DE objective: per-class sign accuracy on held-out lines, priced at H1's per-cell rates of
roughly 2.5e-6 (a0), 5e-6 (a1), 9e-6 (head a1) and 2.7e-6 (u1) overall per corrected cell, before any emitter run.

Not asked: no GO, no cloud beyond what was spent, and no dossier.
