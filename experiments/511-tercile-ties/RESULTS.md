# 511 — tie rules in the per-target sharpening assignment (501 section 6) — prepared and queued (2026-09-27 07:55Z)

**Finding (the cloud review 501, confirmed here).** `sharpen_q_by_tercile` (479) relies on `np.digitize` edge semantics and
a tie rule that fires on either cut.
- tercile: a cap block containing the upper cut went to MID. 493h says capped targets belong in HIGH.
- stercile: a tie block on the lower cut left LOW empty, and collapsed cuts gave every row the mildest q.
- Non-finite values were filled with the median before the cuts, so they moved the cuts.

**Fix** (`patch_builder511.py`): tie blocks are assigned explicitly within 1e-9 relative, which also absorbs the ~1e-14
float noise of r = c·w. The cloud patch relied on exact float ties, and its cases failed here on near-ties.
- tercile: MID only for a cut at the floor; an upper-cut block is HIGH; the floor guard still wins.
- stercile: a lower-cut block is LOW; collapsed cuts give the tied mass LOW and the rest HIGH.
- Cuts are computed over finite values only, and non-finite values sit in MID.

**Checks at prep time** (`test_ties511.py`, old vs new function extracted from the builders):
- q vectors are EQUAL on 446bj2's LODO dump and on 495's cap-6 dump, in both modes (100/100/100 everywhere).
- On efr (465's dump) they are also equal: tercile 117/83/100 (a 117-target floor block), stercile 100/100/100.
- Synthetic ties are fixed:
  - tercile cap of 150: [100, 200, 0] → [100, 50, 150];
  - stercile 140 tied weakest: → [140, 60, 100];
  - collapsed cuts: → [220, 0, 80];
  - the distinct-value control stays [100, 100, 100].

**Queued** (`chain511.sh`, serial after chain510). Steps: patch; builder suite and unit tests; `test_ties511.py` against the
live builder; KA-SP2, where 479sp2's own configuration must rebuild 479sp2's output byte for byte. Then
`commit_after_chain.sh` commits the builder (expected **510f96bd**). No A/B/C or efr assignment changes, so no board read
is needed. It guards October panels whose strengths or keys tie.
