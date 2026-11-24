# 614: plan and checkpoints (K562 real-emitter DE harness)

Owner: [21dd4c]. Spec: `experiments/631-nmae-headroom/SPEC-614-K562-harness.md`. Gate: `TRUSTGATE.md`, fixed first.
Each stage commits its outputs and a line here, so [d4cd5f] or anyone can resume.

| stage | what | output | status |
|---|---|---|---|
| 1 | K562 truth: k562harness480 with N 250, MIN 100, OUT here | work/k562-gwps-c250.h5 | DONE: 223 targets (100–250 cells, median 188), 18,400 NTC, 8,248 genes, integer counts |
| 2 | target subset: seeded 100 of the targets with ≥ 100 cells (seed 614) | work/targets614.txt | DONE |
| 3 | reference: evaluate323.reference() logic without H1 asserts (`ref614.py`) | work/ref_k562/ | DONE. KA: the H1 build is IDENTICAL to 323's reference. K562: median n_real 5 (q90 64, max 386): LOW POWER |
| 4 | key: Orion-only consensus proxy, scale s fitted on H1 only (`key614.py`) | work/key614.npz | DONE: s 1.124; on H1 the donor mean tracks cap15 at r 0.73; 100/100 targets keyed |
| 5 | wrapper `run614.py` (swaps the reference loader, key, CLEAN4 banks, CLASSES3; NaN-safe metrics); KA vs 435 | chain614b | DONE: KA PASS (max \|Δ\| 0) |
| 6 | gate runs G1 / G2 / G3, 3 seeds (chain614e2, after amendment A1's sub-gate reservoir) | work/g1, g2*, g3cp | DONE |
| 7 | gate verdict (`gate614.py`, mechanical) | work/gate614.json | **FAIL** (G1 PASS; G2, G3 FAIL): 614 is DESCRIPTIVE (RESULTS.md) |
| 8 | (descriptive) H1@11k control | work/ref_h1t11k, g1t, g2t_* | running (chain614d) |
| 9 | levers: no PASS, so hook forms run as descriptive rows (chain614h; PRED-614-levers.md registered first) | hk_* | running |

Rules: gate602 for jobs > 4 GB; the coordinator's build windows (10-04 05:00Z – 10-05 18:00Z); own dir only;
commits with a pathspec. Target: the gate verdict by 10-04 05:00Z.
