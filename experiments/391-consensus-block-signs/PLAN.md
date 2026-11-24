# 391 — PLAN: cross-line sign prior on the called CPM block through the unwelded emitter, six members on H1 (PI-run)

Hypothesis: 390 C5b — CLEAN4 consensus signs (371: K562-GWPS, K562-essential, Nadig Jurkat, Nadig HepG2; no H1 data) on
the 2,200 called CPM-block genes at the delivered magnitude, down-signed magnitudes shrunk by one per-target factor so the
block's arithmetic mass is unchanged — delivers analytic H1 fid +0.0116 [+0.0057, +0.0177] with pds/mse pinned by
construction. The end-to-end emitter must reproduce that within its steering/placement losses.
Inputs: 389's emitter snapshot (`k_emit391.py`, sha recorded), 323 reference/effects/anchors (read-only), 371
`consensus_prior.npz`. Harness: 389's exact chain (385's cap15 sequence) with K = 3 classes at h1 control quartiles.
Arms: cap15 (KA: must match 323 bit-for-bit), k3c5b (the arm), k3c5bshuf (prior deranged across gated genes, seed
20260913: the targeting control — must not gain fid), k3down (all-down block, same balance: 390 R2 analogue).
Readout: six raw members vs cap15, paired 10,000-target bootstrap, H1-local scaling (report389 method).
Kill: k3c5b fid gain < +0.004 raw or not above k3c5bshuf; or pds/mse not exactly pinned; or overall CI entirely < 0.
Scope: H1 development evidence (one line, 85 targets), not an A/B/C forecast. No builder, no submission.
