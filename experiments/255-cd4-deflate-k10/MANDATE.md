# 255 — BUILD `cd4svd10` (k = 10 deflation) overnight so branch A of `E249/DECISION-TREE.md` has its package ready. DO NOT SUBMIT.

You are a build subagent. Work in `experiments/255-cd4-deflate-k10/`. **Edit nothing outside it**: the builder
already has `VCC_CD4_DEFLATE_K` (249, `be029c8`) and its flag-off byte-identity was re-proven an hour ago by 254
(`8df4dadd…` through the current HEAD builder) — do not rebuild the incumbent, cite 254 §KA 1. You may build into
`submissions/` under the `255-` prefix only. **NEVER `vcc submit`.** `uv run python`. Builds ~20 GiB and `vcc prep`
~25 GiB: strictly serial, never alongside another build/prep (none is running; check anyway), gated on
`memory_pressure` ≥ 40% free AND `vm_stat` free+inactive ≥ 30 GiB; record peak RSS.

## Why
If tomorrow's `cd4svd5` board point lands in branch A (`experiments/249-cd4-deflate/DECISION-TREE.md`), the next
one-variable step is the k ladder at **k = 10**: 246 §3.3/§4 measured `B2_cd4_svd10_resolved` at `pds` +0.01281
(3/3), macro −1.67%, aggregate −1.54%, price +0.0015 / +0.0029 / +0.0075 at 0.7/1.0/2.0× `pds` transfer; 250 found
the `pds` plateau off A/B/C sits at k 10–15 with the price flat past 5–10. A board point at k = 10 maps the
k-response for October. The package must exist before the reading, not after.

## READ FIRST
`docs/READ-FIRST.md`; 249's `PREREGISTRATION.md`, `RESULTS.md` §1–§5, `build249.sh`, `stage249.py`, `repro249.py`,
`metadiff249.py`, `profile249.py`, `probe249.py` (reuse them with k = 10); 246 §3.3 and §4 (the k = 10 rows);
243 §5.2 (measure the board cosine off the package); `AGENTS.md` §11.1 rule 10, §14.

## Known-answer rows (abort on failure)
1. The staged arm at `VCC_CD4_DEFLATE_K=10` reproduces 246's `B2_cd4_svd10_resolved` on the three cancer folds
   through 162's instrument on the builder's own staged effects — macro −1.67%, aggregate −1.54%, Δ`pds` +0.01281
   (read 246's exact digits from `price246.json`/`analyse246.json`) — to the digit; `energy_removed_frac` ≈ 0.286.
2. `VCC_CD4_DEFLATE_K=5` still stages to 249's `a28d5be2…` (the k = 5 path untouched; nothing changed, so this is
   a cheap re-hash of 249's staged file if it exists, else a re-stage).

## Build
`submissions/255-cd4svd10-f034.{h5ad,vcc}` with `VCC_CD4_DEFLATE_K=10` and everything else the shipped
configuration; 249's full gate set with `submissions/249-cd4svd5-f034.h5ad.meta.json` as the incumbent reference
(rule-10 diff = exactly `cd4_deflate_k` 5 → 10 plus its CD4-frame consequences; profile split by carried row,
scale-corrected; non-vacuity; `n_pred` ±2%; delivered amplitude via 243's stream, report the predicted raw `mse`;
`vcc prep --dry-run`). Record sha256 and bytes. Write `PREREGISTRATION-DRAFT.md` for the lead: six-member
prediction vs `cd4svd5` from 246 §4 at three `pds` transfer factors, and five falsifiers in 249's form.

## Deliverable
`RESULTS.md`: the two rows, every gate's number, the package, the draft preregistration, and one line
**BUILDABLE / NOT BUILDABLE**. No model name, no description, no submission. Then stop.
