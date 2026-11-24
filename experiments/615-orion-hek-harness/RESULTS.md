# 615 — the 614 real-emitter DE harness on Orion RAW HEK293T

Coordinator's assignment (2026-10-02 ~18:00Z), after 614's K562 gate failed on truth power.
- **Order:** (0) a preregistered technical-blob test; (1) the K562-only proxy key with all four CLEAN4 banks; (2) the
  SAME gate as 614, with bands committed before any run; (3) HCT116 raw only if HEK passes and time allows.
- **Heavy stages** go through gate602, outside the coordinator's 10-04 05:00Z – 10-05 18:00Z window.

## Stage 0 — technical-blob test: **PASS** (`blob615.py`, `work/blob615.json`; prereg `PREREGISTRATION-blob.md`, b91a9b59)
**Extraction** (`extract615.py`, 178 s, peak 7.1 GB, ungated: my slip, all later stages are gated):
- 223 batches; 66,289 cells of the 300 panel targets.
- 235 targets have ≥ 100 cells (median 200).
- 218,838 controls, of which 300 per batch are kept as the harness pool's source.
- **Depth on the scorer's 18,533-gene axis: 15.0k** (panel cells) / 15.7k (controls). Orion's quoted 20.7k is
  all-feature `total_counts`. That is 0.75× A/B/C's 20.1k and 1.33× K562's 11.3k.

**The statistic** (known-answer rows first, `blobka615.py`): this definition reproduces the published shares 44.2 /
26.5 / 3.8 % vs 44.1 / 26.4 / 3.6 %.

**Decision statistic: the split-half share** on the harness subset (100 targets, ≤ 250 cells, median 220).

| surface / control profile | split-half share | plain share | mean pair cos | top-PC share (centred) |
|---|---|---|---|---|
| **raw HEK, P′ (one pool, the harness's batch mixture) — decides** | **0.0597** | 0.0143 | +0.0050 | 0.176 |
| **raw HEK, BM (each target's own batch mixture) — decides** | **0.0598** | 0.0144 | +0.0050 | 0.168 |
| raw HEK, P (all controls) | 0.0637 | 0.0150 | +0.0057 | 0.176 |
| raw HEK, G = the old greedy builder's selection, from raw cells | **0.8158** | **0.4369** | +0.4156 | 0.165 |
| reference: orion-hek293t-thinned | 0.8460 | 0.4420 | +0.4233 | 0.133 |
| reference: orion-hct116-thinned | 0.8760 | 0.2646 | +0.2529 | 0.067 |
| reference: h1-test (Arc's chemistry) | 0.0487 | 0.0376 | +0.0257 | 0.159 |

All 235 eligible targets: P′ 0.0430, BM 0.0441.

- **Verdict: PASS** (P′ 0.060 and BM 0.060, both ≤ 0.10). Raw HEK with pooled controls carries about h1-test's blob
  (0.049), not the thinned surfaces' (0.85).
- **Mechanism, now nailed:** G reproduces the thinned HEK truth from raw cells (plain 0.437 vs 0.442; split-half 0.816
  vs 0.846).
  - G's selection: controls from the first 19 batches in file order, the first 150 cells per target.
  - So the Orion blob was the old builder's greedy control selection (a batch/depth mismatch), not the platform.
  - Per-target batch matching adds nothing over a pooled sample (BM ≈ P′ ≈ P), because each target's cells spread over
    the batches like the pool's.
- **Descriptive caveat, carried forward:** the dominant direction of BETWEEN-target variation is still the long-gene
  axis.
  - The top PC of the mean-removed effects has |cos| 0.87 with the thinned blob direction (P′ 0.871, BM 0.868), and it
    carries 17.6 % of the mean-removed energy.
  - FINDINGS H.4 saw the same in the Orion donor tables (top-PC cosine 0.872, 19.7 %).
  - It is per-target spread along that axis, not a shared offset, and batch matching does not remove it.
  - The residual mean vector is anti-aligned with the blob direction (cos −0.39), so there is no greedy-type shared
    offset.
  - For DE-channel members this is a per-target technical or biological signal that the key cannot predict. It lowers
    achievable fid/sign accuracy on the affected genes, but it does not bias every target the same way.

## Stage 1 — truth, reference, keys (chain615a, 18:13–19:02Z)
- **Truth** (`truth615.py`, `work/truth615.json`):
  - 100 targets, 220 cells median (min 102), 18,400 P′ controls from 222 batches.
  - Library median 15,036 perturbed / 15,179 control (ratio 0.99; the old greedy harness had 0.725).
- **Knockdown QC** (`kd615.py`): median 48 % of the target's own expression remains (52 % CRISPRi knockdown); 7 / 100
  keep > 70 %. The labels are sane.
- **Reference** `work/ref_hek` (ref614.py): 11,063 genes in the DE universe.
- **Proxy key** (`key615.json`): k562_gwps × s = 0.496 (fitted on H1: 629,677 cells, Pearson 0.449, sign agreement
  0.663). 86 / 100 targets have a key row.
- **Hook signs** (`idsign615.py`): c4h 2,617 signs (1,530 down); idall3 6,346; they agree on 0.992 where both sign.

## The low-power diagnostic, read BEFORE any gate run: **K562-like (worse)** (`diag615.py`, `work/diag615.json`)
| truth | median n_real | q90 | targets with n_real ≥ 10 (nmae-scored) | DE calls | DOWN share |
|---|---|---|---|---|---|
| **HEK 615 (raw, 220 cells, 15k)** | **1** | 26 | **18 / 100** | 1,472 | 0.666 (0.650 without own genes) |
| K562 614 (GWPS, ≤ 250 cells, 11k) | 5 | 64 | 34 / 100 | 2,213 | 0.435 |
| H1 (h1-test, 400 cells, 20k) | 89 | 2,691 | 58 / 85 | 65,130 | 0.526 |

- **Not a construction artefact.**
  - The target's own gene is significant in 69 / 96 targets (median lfc −1.04), so the test has power for real
    effects.
  - Downstream responses are tiny: 30 / 100 targets have 0, and the q75 is 5.
  - Four targets (SFXN4 299, MED13 230, PARP3 178, FAM120A 100) hold half of all calls.
- Orion's weak knockdown (52 %, against ~80 %+ in Replogle's CRISPRi lines) is a likely contributor. So is the panel's
  design: targets were drawn to span negligible-to-strong responses in H1, not in HEK293T.
- **The thinned HEK gate table (480) had a median of 141.5 calls per target.** This clean truth has 1, with more cells
  and more depth. So the thinned Orion DE truths were almost entirely the selection artefact (quantified in 616).
- Per TRUSTGATE the gate still runs with unchanged bands; the coordinator was told before it started (19:05Z). Expect a
  power FAIL like 614's.

## Gate stopped; wrapper KA PASS (19:40Z)
- The coordinator stopped the gate before any gate row (TRUSTGATE.md, written reason).
- The wrapper KA ran anyway: `run615 h1ka` == 435's H1 stats on both stats files (sums equal, max |Δ| p / q / lfc 0,
  same NaN pattern). run615 is validated for any later use of this harness, e.g. a high-signal subset.
- **Disk:** the raw extractions `work/hek-raw.h5` (4.2 GB) and `hct-raw.h5` were deleted after their consumers ran (615's
  truth, 616's clean truths and gate tables). `extract615.py` regenerates each in ~3 min.
