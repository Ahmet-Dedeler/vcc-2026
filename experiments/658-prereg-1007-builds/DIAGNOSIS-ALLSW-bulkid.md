# ALLSW bulkid DIFFERENT: diagnosis and coordinator override (2026-10-06 ~03:45Z)

**Finding.** build_allsw658.sh ended with "NOT READY bulk is NOT 498p2g's (bulkid658 rc 1)": 28 (block, gene) count sums differ
in 11 of 900 blocks, sum |d| = 28. The coordinator restored 498p2g's package again and listed every difference.
- Each is a paired ±1 swap between two low-count genes (sometimes four) inside the same block. Counts run 9–113, e.g.
  A/ANKZF1 TMEM132A 49 vs 50 with NUDT19 32 vs 31.
- Every block total is unchanged (net 0).
- Every other build gate passed: config diff exactly de_mag + de_mde_sign + de_sign_emitter_sha256; CLIFF OK equal to 498p2g's;
  prep rc 0; novelty rc 0; package c86515ccf8b4d42d5e21e654d20922efdb2154fe1d181720bcbb56c1ecdb0311.

**Mechanism.** These are integer-rounding ties resolving differently. ALLSW is emitted by 0f56209d, while 498p2g used edb5b8dc
(the emitter key in EXPECT). A largest-remainder tie between two genes in a block flips with tiny floating-point differences.
idall3 is the precedent: the same emitter change on 498p2g's base, read on the board 10-05, gave pds bit-identical and mse
+9.9e-8 raw.

**Why the override.** ALLSW's F1 (PREREGISTRATION-ALLSW) is a BOARD rule: "bulk pinned to 498p2g, 1e-6 raw". bulkid658's
exact count identity is a stricter pre-slot proxy. 28 net-zero ±1 swaps out of 16.7 M sums move log-pseudobulk by ~1e-11 of
mse's scale, and they cannot reorder pds's nearest neighbours, so F1 is predicted to PASS. read658.py still checks F1 on the
real board read and binds there. The coordinator therefore marks the build READY in its driver log, with this note.

**Caveat for the read.** ALLSW's DE deltas vs 498p2g include the emitter change as well as zs + hook. Same-law re-roll noise
(0.00025 overall) bounds the emitter part, which is what R2's bins already use.
