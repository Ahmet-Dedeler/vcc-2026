FOLLOW-UP from the lab lead (written 2026-10-16 04:20Z; sent when your first pass is done). Your context is intact.

The bulk half of the checker's known answers is graded (`docs/LAB-LOG.md`, "Known-answer report, part 1"): 4/5, and the
miss is 479sp2. The bulk folds price strength-tercile sharpening (S_P2) like key-tercile sharpening (P2g), both about
+0.0036, while the board split them (+0.0010 vs −0.0007). S_P2's board loss sat in pds (−0.0011 raw) and nmae (+0.0018 raw,
"the low tercile's ballast is under-priced", 479's BOARD-RESULT). Those are emission-level effects, and DE-H1 is the only
surface that emits real cells. So:

1. **Bulk sharpening arms on DE-H1**, both surfaces, 3 seeds, with all six members from the emitted cells (your `ce_pds`,
   `mse_num`/`mse_den`, the DE members at f 1 / 0.35): (a) uniform q 1.75 (the 446bj2 base), (b) P2g = per-target tercile
   2.5 / 2.0 / 1.75 by the key ratio (the builder's `tercile`), (c) S_P2 = the same by strength terciles (the builder's
   `stercile`), each on the surface's own bulk with everything else equal. Use the builder's own per-target sharpening code
   path if the H1 harness can route it (599-tercilec-h1-de did tercilec on H1: read its RESULTS and reuse its wiring);
   otherwise say exactly what blocks it. KA: (a) reproduces your existing base runs bit-for-bit.
2. Report per pair (P2g − a, S_P2 − a, S_P2 − P2g) the six member deltas per surface and regime with seed SDs, beside the
   board's (479-pair BOARD-RESULT). The question the lab needs answered: does DE-H1 see S_P2 lose to P2g (pds, nmae), as the
   board did?
3. Rules update: locally only file reads and sub-second < 100 MB aggregations; anything else via gate602 or EC2; never while
   `bash experiments/800-checker/window800.sh` prints WINDOW. Keep big outputs on S3 (the runner's collect now fetches only
   logs and declared repo/ paths). Budget for this follow-up: $25.
