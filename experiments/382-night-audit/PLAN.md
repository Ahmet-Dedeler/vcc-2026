# 382-night-audit — PLAN

Independent re-derivation of five headline numbers from raw saved files, code written
from scratch. The audited *scripts* (`evaluate323.py`, `exact_de323.py`, `run370.py`,
`run345*.py`, 363/346 runners) will not be opened until my numbers exist; only
reports (`RESULTS.md`), metric specs (`docs/reference/vcc2026-metrics-BRIEF.md`),
and raw data files are read beforehand.

## Hypotheses (one per item — each is "recorded value reproduces")

1. **363 fid/nmae.** From `stats-{cap15,block_down,block_random}.npz` (own p/q/lfc per
   target×universe gene) + 323 `reference.npz` (reference q/lfc, genes, keep, targets)
   I can reproduce the recorded raw fid and nmae in
   `experiments/363-block-down-sixmember/work/metrics-*.json` to tol 1e-9, using only
   the BRIEF spec: fid F_p=k/max(n_pred,n_real), own target excluded, adjudicable =
   finite nonzero reference lfc, ours significant = own q<0.05, 0/0 dropped;
   nmae = Σ|ours−ref|/Σ|ref| over the reference gate (q<0.05, finite ref lfc,
   target removed), gate ≥10 genes. Cross-check gate sizes against
   `reference_de.parquet`. pds/mse are not recomputable from these files (no design
   effects, no emitted counts) — will state that, not attempt it.
2. **346 scaled Δ.** From the raw members in 346 `metrics-*.json` + 323
   `h1-anchors.json`, scaled=(raw−b)/(r−b) with mse clamped to [0,1] (clamp direction
   verified: cap15 raw mse 0.9749 → scaled 0.0, so clamp-after-scale, lower bound 0),
   overall=mean of six. Expect recomputed Δ to match `results.json` contrasts at full
   precision; the assignment's 4-decimal headlines (+0.0071/−0.0149/−0.0123) can only
   match to ~3e-5 (rounding) — will report both comparisons separately.
3. **370 even/odd prior.** Own code: context-A mean CPM via h5py over
   `X/data|indices|indptr` in row chunks (≤6 GiB, free is 80%); block = top-2,200 by
   A-CPM among Jurkat-mapped genes; per-gene P(down) fit on even-indexed targets,
   evaluated held-out on odd (and vice versa), eval set = block minus target head
   (top-400 |delta|), exact-zero deltas dropped. Expect prior ≈0.5766, all-down
   ≈0.5147, tol ±0.003. Known answers first: A shape 18400×18533, median lib 20109,
   gene-order head TSPAN6/TNMD/DPM1, Jurkat own-target negativity ≈99%.
4. **Board deltas.** Entry files identified via `docs/BOARD-LEDGER.md` notes:
   319=`Q519I2pLWRCW7PAEx5ul`, 319p=`SvzrTOQJAAaNjawzSyrC`,
   319pA=`sEa0LAsWFnXOoXSPJZ59`, 254=`PVHrJHweHD5gvOf43ftL`,
   cd4svd5=`Xv0ayvuoHpsKYWv2iu42`, 319d3=`5LOjb2VLiIBh0JG6aAgg`.
   Recompute 319p−319, 319pA−319p, 254−cd4svd5, 319d3−319 from `score_avg`, and check
   each `score_avg` == mean of its six `score_*` (tol 1e-12).
5. **345 withdrawal.** K562 all-row mean of `data/donors/k562_gwps/delta.npy`
   (mmap), mapped by gene name to the H1 universe; |lfc|-weighted sign accuracy per
   target on the parquet gate (58 eligible targets), macro mean. Expect ≈0.527 macro
   / ≈0.547 pooled (PI), not 0.275 (345). Sign convention verified by own-target
   negativity first.

## Kill/falsify criteria

- Any item whose recomputation misses tolerance = MISMATCH: re-check my code once
  (spec re-read, off-by-one, axis alignment), and if it still misses, report
  MISMATCH with both values — do not "fix" by reaching for the audited scripts'
  choices beforehand. Audited scripts may be read only after numbers are final,
  solely to explain a residual discrepancy.
- If a known answer (item 3 shape/median, item 1 cap15 as known-answer anchor, 345's
  gate rebuild of 58 targets) fails, stop that item and report.

## Inputs (all read-only)

- `experiments/323-current-family-validation/work/{reference.npz,reference_de.parquet,metrics-cap15.json,h1-anchors.json}`
- `experiments/363-block-down-sixmember/work/{metrics,stats}-{cap15,block_down,block_random}.json/npz`
- `experiments/346-same-line-calibration-h1/work/{metrics-{cap15,same_ridge_b05,same_ridge_b1,diag_cross_b1}.json,results.json}`
- `data/donors/nadig_jurkat/{delta.npy,genes.txt,perts.txt}`,
  `data/donors/k562_gwps/{delta.npy,genes.txt,perts.txt}`,
  `data/challenge/context_A.h5ad`, `data/challenge/gene_names.csv`
- `data/board-entries/{Q519I2pLWRCW7PAEx5ul,SvzrTOQJAAaNjawzSyrC,sEa0LAsWFnXOoXSPJZ59,PVHrJHweHD5gvOf43ftL,Xv0ayvuoHpsKYWv2iu42,5LOjb2VLiIBh0JG6aAgg}.json`

## Outputs (all in `experiments/382-night-audit/`)

`audit{1..5}.py`, `audit{1..5}.json`, `RESULTS.md`. No downloads. Peak RSS expected
< 2 GiB (largest single array: K562 delta 323 MB via mmap; H1 stats trivial).
`NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4`. Run serially, one process at a time.
