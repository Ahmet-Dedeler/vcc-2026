# 626 — The identity rider is robust to SISTER references; a cheap guard: drop a reference whose control-deviation similarity to the context exceeds 0.3

2026-10-02, council session [7772dc] (assigned by [4d0207]). Analysis only: Broad-16's 12 working lines (reserved 4 never
read) + Arc / Satija Flex control profiles. `sister626.py`, `work/sister626.json`. < 1 GB, 24 s.

## Question
On Oct 22 the identity guard drops only an exact same-call reference. If a D/E/F line is a SISTER of a reference (same
lineage: two T-ALL lines, two squamous lines), the identity contrast shrinks on the lineage genes and could flip the
rider's sign there.

## Measure
Similarity between a context and a reference = Spearman of their log1p control-CPM DEVIATIONS from the pool mean, over
genes every profile measures. It is truth-free and computable on the day among the Flex contexts.

## Broad-16 results (gate-sign of z(C4) + 0.5·z(I) vs CLEAN4; authors' fdr < .05 gates)
- **Similarities.** Only KYSE140~KYSE450 (esophageal, same lineage) is high: **0.343**. Every non-sister pair has a
  median of −0.08 and a maximum of 0.256. The cervical pair C4I~SKGII is NOT similar in its controls (0.047).
- **The rider stays positive with the sister in the pool.** All 11 references: KYSE140 +0.0050, KYSE450 +0.0149.
  Without the sister: +0.0059 / +0.0127.
- **A sister in a SMALL pool** (3 references = the sister + 2 random, against 3 random non-sisters, 20 draws):
  - KYSE140: +0.0013 vs +0.0029;
  - KYSE450: +0.0125 vs +0.0133;
  - SKGII: +0.0070 vs +0.0110;
  - C4I: −0.0020 vs −0.0047.
  It dampens the gain by ≈ 0.001 … 0.004 and never flips a positive line negative. Identity-alone accuracy drops by
  0.02–0.03 for the KYSE pair.
- **Guard "drop references with similarity > t".** At t 0.3 it drops only the KYSE sister (2/12 lines), and the mean
  rider gain changes by −0.0001. At t ≥ 0.4 it never fires.

## Cost on Arc's Flex pool
- Deviation similarities within {A, B, C, H1} are all NEGATIVE: A~B −0.21, A~C −0.41, A~H1 −0.39, B~C −0.08,
  B~H1 −0.42, C~H1 −0.26. The guard never fires on the board's A/B/C pool, so its cost there is 0.
- In a cross-lab pool, same-lab profiles cluster: Satija's HAP1~HEK293 +0.47, K562~HAP1 +0.27. The measure partly
  tracks lab. On the day every reference is Arc Flex, so lab cannot drive it, and only true line similarity remains.

## Recommendation for `def_identity_prior` (to [4d0207])
1. Keep the exact-call guard (622e).
2. **Add: drop any reference whose control-deviation similarity to the context exceeds 0.3**, computed over the day-of
   Flex pool (A, B, C, the other D/E/F contexts, H1). Log the full similarity row as 597/590 log top-5 calls.
3. Never drop below 2 references. If the guard would, keep the 2 least similar.

A CCLE-lineage guard is NOT recommended. Same lineage does not mean similar controls (the cervical pair: 0.047), and it
would drop useful references.

Class: day-of instrument guard (scorer-independent; nothing fitted to A/B/C). Small stakes either way. The point is that
the October rider cannot be flipped by a sister reference.
