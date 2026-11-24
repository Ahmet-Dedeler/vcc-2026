# 553 — files pruned to save disk (2026-10-01 07:50Z; the disk had 36 GiB free)

Regenerable intermediates whose jobs are finished; each one's json record (with its sha256) stays:

- `h1t7/h1t7.h5ad` (1.5 GB): `uv run python experiments/553-r7b-harness/prep553.py h1t7 build` rebuilds it bit for bit
  (h1r7's target rows + all h1-test controls; check against `h1t7/h1t7.json` sha256 3a2a27d2…). Its access output
  `access553_h1t7/` is kept.
- `r7bd7000/context_R_7000.h5ad`, `r7bd7000/controlref553.npz` (0.41 GB): `prep553.py subset r7bd7000` then
  `prep553.py controlref r7bd7000` rebuild them (shas in the json files beside them). Only the floor sweep used them; its
  diag records `diag553_r7bd7000_*.json` are kept.
- `h1t7s2/h1t7s2.h5ad` (1.5 GB, 11:05Z): `prep553.py h1t7 h1r7s2` rebuilds it (sha in `h1t7s2/h1t7s2.json`); its access
  output `access553_h1t7s2/` is kept. The later ablation files `h1t7s3..5` are pruned the same way after their accesses.
- `h1t7s3/h1t7s3.h5ad` (1.51 GB, 11:44Z, prune553.py): regenerable by its build step; sha256 in the json beside it.
- `h1r7s2/h1r7s2.h5ad` (0.67 GB, 12:26Z, prune553.py): regenerable by its build step; sha256 in the json beside it.
- `h1r7s3/h1r7s3.h5ad` (0.67 GB, 12:30Z, prune553.py): regenerable by its build step; sha256 in the json beside it.
- `h1t7s4/h1t7s4.h5ad` (1.50 GB, 13:27Z, prune553.py): regenerable by its build step; sha256 in the json beside it.
- `h1r7s4/h1r7s4.h5ad` (0.67 GB, 13:28Z, prune553.py): regenerable by its build step; sha256 in the json beside it.
- `h1t7s5/h1t7s5.h5ad` (1.51 GB, 15:11Z, prune553.py): regenerable by its build step; sha256 in the json beside it.
- `h1r7s5/h1r7s5.h5ad` (0.67 GB, 15:11Z, prune553.py): regenerable by its build step; sha256 in the json beside it.
- `h1r7/h1r7.h5ad` (0.67 GB, 15:11Z, prune553.py): regenerable by its build step; sha256 in the json beside it.
- `h1r7c0A/h1r7c0A.h5ad` (0.61 GB, 16:02Z, prune553.py): regenerable by its build step; sha256 in the json beside it.
- `h1r7c0B/h1r7c0B.h5ad` (0.61 GB, 16:35Z, prune553.py): regenerable by its build step; sha256 in the json beside it.
- `h1r7c1A/h1r7c1A.h5ad` (0.61 GB, 17:09Z, prune553.py): regenerable by its build step; sha256 in the json beside it.
- `h1r7c1B/h1r7c1B.h5ad` (0.61 GB, 17:40Z, prune553.py): regenerable by its build step; sha256 in the json beside it.
- `h1r7c2A/h1r7c2A.h5ad` (0.61 GB, 18:12Z, prune553.py): regenerable by its build step; sha256 in the json beside it.
- `h1r7c2B/h1r7c2B.h5ad` (0.61 GB, 18:44Z, prune553.py): regenerable by its build step; sha256 in the json beside it.
- `h1r7/h1r7.h5ad` (0.67 GB, 00:59Z, prune553.py): regenerable by its build step; sha256 in the json beside it.
- `h1r7c0A/h1r7c0A.h5ad` (0.61 GB, 00:59Z, prune553.py): regenerable by its build step; sha256 in the json beside it.
- `h1r7c0B/h1r7c0B.h5ad` (0.61 GB, 00:59Z, prune553.py): regenerable by its build step; sha256 in the json beside it.
- `h1r7c1A/h1r7c1A.h5ad` (0.61 GB, 00:59Z, prune553.py): regenerable by its build step; sha256 in the json beside it.
- `h1r7c1B/h1r7c1B.h5ad` (0.61 GB, 00:59Z, prune553.py): regenerable by its build step; sha256 in the json beside it.
- `h1r7c2A/h1r7c2A.h5ad` (0.61 GB, 00:59Z, prune553.py): regenerable by its build step; sha256 in the json beside it.
- `h1r7c2B/h1r7c2B.h5ad` (0.61 GB, 00:59Z, prune553.py): regenerable by its build step; sha256 in the json beside it.

## Second prune, 2026-10-02 01:21Z (`prune553b.py --delete`; approved by the coordinator for 10-04/10-05 scratch)

676 arrays (16.53 GiB), each hashed before deletion: `work/PRUNED-1002.tsv` (category, path, bytes, sha256). Per category: access 391 files 5.81 GiB; draw-arm 174 files 5.55 GiB; chalf-surf 6 files 1.02 GiB; chalf-arm 72 files 2.30 GiB; r7b-arm 5 files 0.19 GiB; sweep 28 files 1.66 GiB.
Kept: R7b's access inputs (`r7b/`, `r7bd5000*`), the h1r7 re-pricing inputs (the five draws' controlref and callblock, `h1r7[sN]_cpm`, `h1/bulk553.npz`), the known answers' files (`h1_*`, `access553_h1`, `access553_rehearsal5`, `access553_r7bd5000-xBST2`, `h1cb*`, `h1ka`, `test553`) and every json / md / csv / log.
Regenerate: arms with `run553.py <surface> <MODE> base <seed> ...` (emission is deterministic: KA-FINAL); control halves with `prep553.py chalf h1r7 <k>`; the floor sweep with `prep553.py subset r7bd<floor>` + `prep553.py controlref r7bd<floor>`; an access's caches by re-running the access (its results json is kept).
