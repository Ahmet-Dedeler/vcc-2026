# 801 — the lab's fresh held-out folds and the bulk fold engine (running record)

Worker: Claude Opus 5.5 under `MANDATE.md` (lab lead [c61d4c]). Writes only in `experiments/801-lab-folds/` (plus the lab
runner's own logs, `experiments/800-checker/cloud/runs*` and `docs/LAB-LOG.md`'s cost table, which `labcloud.sh` appends).

## 0. Answer (2026-10-16 02:50Z)
- **The draws are done.** Eight Arc-like panels were drawn by rule arc800 (`draw801.py`, seeds 800101–800108), `work/draws.json`
  sha256 `8c8e67d2cecd26ef16a34ba8c3fa33996c61975bed5ea714bd900f21fd3cd499`.
  - L1–L6 are WORKING: bundles, CD4 banks, 18 held-out operators and 18 fold dumps.
  - L7–L8 are SEALED: only their target lists, seeds and sha256 exist. No bundle, bank, operator, dump or statistic was made.
  - Every panel matches A/B/C on donor coverage (K562 269, CD4 283) and on the coverage × strength strata. Two Q4 slots moved
    to Q5 and one O1 slot to O2, by the rule frozen before the first draw.
- **The engine reproduces its sources.**
  - `engine801.py` on 522's own efp dumps (Mac), and on fresh EC2 dumps of efp made by `dump801.py`: every 546o, 551 and 522
    row reproduces, worst |d| **0.0 (Mac) and 1.6e-15 (EC2)**, 48 rows. This includes the incumbent (551 `h2c1|V1pct+ptP2g`).
  - dump801's own known-answer block (dump522's) passes on EC2 with worst |d| 1.6e-15 on 18 rows.
  - On A/B/C's four PLOLO folds the engine matches 540c / arms540 / 522 at |d| 0.0 (48 rows).
- **Registry:** `work/registry.json` (sha256 `8630992e…`) holds 25 folds: L1–L6 × 3 working, ABC × 4, efp × 3 as known-answer
  folds. L7–L8 are listed as sealed-drawn. Operators are kept on S3 (`s3://…/output/lab/ops801/L<k>/`) under their sha256;
  only the dumps came to the Mac.
- **Smoke:** the incumbent vs `cap3` and vs `pct+S_P2` on all 22 folds (`work/smoke/`), fold-level r and Δpds only, no board
  units. On L1–L6 the incumbent beats cap3 by r +13 … +35 % on every fold; on A/B/C by +2.1 … +7.0 %. Against pct+S_P2 it is
  within ±0.6 % everywhere.
- **Cloud: $10.55 for the first pass (budget $40); +$0.83 for the 710 designs follow-up (budget $15, §11).** That is eight launches, five of which failed on path and environment
  defects (§4). No instance of mine is still running.
- **Three findings for the lab lead:**
  - **dump708.py crashes after saving its npz.** Its json record uses `default=float`, which fails on the Path-valued
    `fold_ops` that its own `--fold-ops held` CLI builds. The 10-22 day-of run would exit with rc 1. Fixed only in dump801,
    not in 708.
  - **Operator fits are not bit-reproducible across machines or thread counts.** The efp held-K562 fit gave three different
    sha256s: Mac `16ea9d34…`, EC2 8xlarge / 5 threads `2a0247d2…`, EC2 4xlarge / 6 threads `f2fd179f…`. The same instance
    type and thread count does reproduce: L1 held-K562 fitted twice gave `a97587a7…` both times. The prices do not care: the
    refit moves cos_agg by ≤ 1.4e-11 and moves pds and K by 0 on seven arms (`work/refit801-k562_gwps.json`). Record
    operators by sha256; a refit is a different file but the same price.
  - **A/B/C-matched panels still have much heavier percentile-cap tails than A/B/C itself.** The fold cap `hi` on L1–L6 is
    7.8–19.9; on A/B/C it is 5.2–7.0. So the incumbent's key gain on fresh folds is 3–10× its A/B/C value, as 522 §4b saw on
    its Arc-like re-draws (hi 9.5–10.9). Matching coverage and the CD4 strength mix does not reproduce A/B/C's tail. Either
    A/B/C's panel is unusual, or Arc's selection used something the GWCD4i strength bins do not capture. The checker should
    not read L-fold key gains as board-sized without the over-read ratio (DESIGN §7).

## Status log
- 2026-10-15 22:54Z: `draw801.py` (sha256 `c4b8507e…`, `work/draw801.sha256`) dry run → `work/dryrun801.txt`; draw at
  22:54Z → `work/draws.json` sha256 `8c8e67d2cecd26ef16a34ba8c3fa33996c61975bed5ea714bd900f21fd3cd499`.
- 22:55Z: `diag801.py` → `work/diag801.json` sha256 `44c7e0e1524f6010088e0cb5b4f4091f453981bfb8dea9dba75cf60cd4fa08f9`
  (L1–L6, A/B/C, efp, efr only; written before any score).

## 1. The draws (rule "arc800", `draw801.py`)
- Exclusion at draw time: 24 `pert_counts.csv` files (19 distinct panels), 4,373 targets; K562 essential list 2,021 genes.
- Basal threshold: A/B/C panel 5th percentile of pseudobulk control CPM (mean over A/B/C's context files) = 13.4122 CPM.
- Strength edges (A/B/C quintiles of GWCD4i Rest `n_downstream`, 283 A/B/C targets carried): 0.4 / 2 / 7 / 31.6.
- Pool at the first draw: 5,125 targets. Strata = coverage class (C3 / K2 / O2 / O1) × strength bin (Q1–Q5, NA), quota =
  A/B/C's count, capacity = floor(pool cell / 8). Three cells were short and moved by the frozen rule (dry run, before the
  first draw): C3|Q4 → C3|Q5, K2|Q4 → K2|Q5, O1|Q2 → O2|Q2. K562 coverage (269) and CD4 coverage (283) are exact.

| panel | seed | status | sha256 (sorted targets, `\n`-joined) |
|---|---|---|---|
| L1 | 800101 | working | `c09d6f1db42d96cd…` |
| L2 | 800102 | working | `c8a597c4d001a352…` |
| L3 | 800103 | working | `e587ba0fb6e5db48…` |
| L4 | 800104 | working | `cbdb34a7781dc945…` |
| L5 | 800105 | working | `deaaa8293a148f51…` |
| L6 | 800106 | working | `100b0f49231f314f…` |
| L7 | 800107 | **sealed** | `97c30e7d11dddc01…` |
| L8 | 800108 | **sealed** | `32fd27d28425b12e…` |

Full hashes in `work/draws.json`. Every panel: 300 unique, 0 essential, 0 overlap with any panel file or earlier lab draw,
all on the official axis, strata counts == quotas (asserted). Bundles for L1–L6 at `work/L<k>/bundle/` (control files are
APFS clones of 199's relabelled `def-bundle` files, sha256 re-checked). L7–L8: only `work/sealed/L7|L8/pert_counts.csv`
(target lists, so any later "every pert_counts.csv on disk" exclusion covers the sealed reserve).

### Diagnostics (`work/diag801.txt`; L1–L6 only)
| panel | ess | K562 | HCT | HEK | CD4 | Q1 | Q2 | Q3 | Q4 | Q5 | NA | C3 | K2 | O2 | O1 | n_down median | CPM q5 | CPM q50 | < thr |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| L1 | 0 | 269 | 291 | 293 | 283 | 57 | 28 | 84 | 55 | 59 | 17 | 259 | 10 | 31 | 0 | 4.0 | 16.98 | 64.68 | 0 |
| L2 | 0 | 269 | 290 | 295 | 283 | 57 | 28 | 84 | 55 | 59 | 17 | 259 | 10 | 31 | 0 | 4.0 | 18.74 | 61.80 | 0 |
| L3 | 0 | 269 | 290 | 294 | 283 | 57 | 28 | 84 | 55 | 59 | 17 | 259 | 10 | 31 | 0 | 4.0 | 15.67 | 58.81 | 0 |
| L4 | 0 | 269 | 290 | 294 | 283 | 57 | 28 | 84 | 55 | 59 | 17 | 259 | 10 | 31 | 0 | 3.0 | 17.07 | 62.82 | 0 |
| L5 | 0 | 269 | 290 | 294 | 283 | 57 | 28 | 84 | 55 | 59 | 17 | 259 | 10 | 31 | 0 | 4.0 | 16.38 | 63.82 | 0 |
| L6 | 0 | 269 | 291 | 293 | 283 | 57 | 28 | 84 | 55 | 59 | 17 | 259 | 10 | 31 | 0 | 3.0 | 16.26 | 64.38 | 0 |
| A/B/C | 0 | 269 | 290 | 294 | 283 | 57 | 28 | 84 | 57 | 57 | 17 | 259 | 10 | 30 | 1 | 4.0 | 13.41 | 56.99 | 15 |
| efp (534's efp rule) | 0 | 146 | 293 | 296 | 183 | 52 | 40 | 46 | 19 | 26 | 117 | 140 | 6 | 153 | 1 | 1.0 | 0.00 | 11.51 | 156 |
| efr (534's efr rule) | 0 | 300 | 300 | 300 | 295 | 80 | 59 | 78 | 39 | 39 | 5 | 300 | 0 | 0 | 0 | 2.0 | 0.03 | 49.88 | 68 |

Reading: L1–L6 match A/B/C on coverage and strength by construction (two targets Q4 → Q5, one O1 → O2). The two older
rules miss A/B/C on both axes: efp has half the K562 coverage, 117 atlas-less targets and half its targets below A/B/C's
basal floor; efr is weaker (Q4+Q5 78 vs 114) and has no K562-lacking targets.

## 2. Engine and dump step (code; every copied file names its source and git blob in its header)
| file | what | source (git blob) |
|---|---|---|
| `guard801.py` | write guard (any write under the repo outside 801 raises) + data-closure recorder (audit hook + h5py) | guard534.py `545769bb…` |
| `lib801.py` | lib438 verbatim, caches under `work/lib/` | lib438.py `6995209e…` |
| `dump801.py` | dump708 verbatim through lib801 (+ `--registry`, `--ka522` = dump522's own KA block verbatim) | dump708.py `175666ad…`, dump522.py `42ff0d2f…` |
| `engine801.py` | the bulk fold engine (below) | fold522.py `070ffac1…`, sp2_522.py `1c89e5c3…`, arms522c.py `f7a42e63…` |
| `cd4_801.py` | CD4 banks for one WORKING panel (refuses L7/L8) | cd4_534.py `24d6206e…` |
| `discover801.py` | the closure finder: dump801 / fold522.Abc with an identity stub operator (numbers meaningless; files only) | new |
| `test_local801.py` | engine801's KA on dumps already on disk (522's dump522-efp-*) | new |
| `smoke801.py`, `cloud801.sh` | the smoke; the EC2 side (verify / panel / refit-efp / ka-efp / smoke) | new |

### Engine KA on the Mac, before any cloud run (2026-10-15 23:13–23:15Z; `work/ka/local-*.json`)
engine801 on 522's own local dumps (dump522-efp-*.npz) and on A/B/C's 518 dumps (fold522.Abc). **Every row |d| = 0.0.**
- efp ×3 folds, 16 rows each: PT0 (cap3, pct), 546o `base q1.75` / `fixed 6` / `p96.6` / `p96.6 + P2g` / `p96.6 + S_P2`,
  551 `h2c1|V1_5_3+sh1.75` / `V1pct+sh1.75` / `V1pct+ptP2g` (the incumbent) / `V1pct+ptS_P2` / `V1_5_3+ptP2g` / `V1_5_3+ptS_P2`,
  522 arms522 `cap3` / `cap6` / `pctPROD` (cos_agg, pds, K_mean). Peak 0.61–1.33 GiB, 6–16 s.
- ABC ×4 folds (k562, hct116, hek293t, kolf2), 12 rows each: PT0 ×2, fold522.Abc's own (arms540 V1_0.5_3 / V1_0.5_6),
  540c `V1_0.5_6 + P2g` / `+ S_P2` / `V1_0.5_3 + P2g` / `+ S_P2` / `V1_0.5_3` / `V1_0.5_6`, 522 arms522 ×3. Peak 1.0–1.3 GiB.
- 708's rehearsal row: **not present** (708's rehearsal waiter is still armed behind the science chain; no stored row exists).

### Data closure (`work/closure/*.txt`, discover801, 23:05Z)
dump801 on efp (3 folds) and fold522.Abc (4 folds) opened 88–90 repo files each (code excluded below). Data: data/challenge
(gene axis, panel, context files via the frame caches), data/donors/{k562_gwps, orion_hct116, orion_hek293t, h1_2025,
cd4_deflation_basis}, data/external/perturbmatch, 115's compact H1 donor, the panel's CD4 banks, the fold operators, 410's
build_b meta + check410, 414's arms414 / beta414 / decomp414, 420's comp420 + plolo frames, 438/546o/540o/551/522 stored
rows (KA only), 518's dumps (ABC only), and lib caches (bankstats cloned from 438's, sha-identical; abc_frames). Peak RSS of
the stub runs 1.05–2.4 GiB (the two Orion folds read 2.25 / 2.41 GiB: slightly over the 2 GB local guideline, ~6 s each).

### Staging (S3 mirror `input/lab/repo/`; the Mac uplink measured ~0.9–2.4 MB/s, so big inputs come server-side)
- Server-side S3→S3 copies from `backup-20261002/` (no Mac traffic): 46 files (data/challenge, the five donor dirs,
  perturbmatch, 199's GWCD4i h5ad 16.8 GB); every one is re-hashed on EC2 against the Mac's sha256
  (`work/manifest/local_sha256_backup.txt`) before any job uses it (`cloud801.sh verify`).
- Uploaded from the Mac (< 1 GB in all, finished before 23:30Z): 199's three relabelled control files (632 MB; L1–L6 and
  efp's bundles all carry these bytes), 115's compact donor (209 MB), ~170 MB of small records (list in
  `work/manifest/small_stage.txt`, `stage2.txt`), bankstats, abc_frames.
- NOT uploaded: the L1–L6 CD4 banks (re-extracted on EC2 from the mirrored h5ad; every non-json file must hash equal to the
  Mac's extraction, `work/manifest/local_sha256_cd4.txt`); operators (fitted on EC2 and kept on S3).

## 3. The engine's interface (`engine801.py`)
```python
import engine801 as E
F = E.load("L3/orion_hct116")          # or "ABC/kolf2", "efp/k562_gwps"; check=True runs the fold's known answers first
F.targets, F.n, F.kv, F.carried, F.strength, F.frames, F.cpm, F.scored, F.nonpanel, F.fin   # fold objects (dump fields)
arm = E.named(F, "incumbent")          # cap3, cap6, pct, pct+P2g (= incumbent), pct+S_P2, cap3+S_P2, cap3+P2g, cap6+P2g, cap6+S_P2
arm = E.Arm(w=..., q=1.75 or vector, fin=None or rows)   # any bulk candidate; E.blend_rows(F.fin, ext, lam) builds `fin`
summary, arrays = E.price(F, arm, base=E.named(F, "incumbent"))
```
- `w`: per-target key weights (n). `q`: a scalar (fold522's whole-bulk sharpen) or an n-vector (rows sharpened per q group,
  sp2_522 / comp546o); KA-PT0 proves the two paths are bit-identical at a constant q. `fin`: replacement h2c1 composition
  rows (n × 18,533; panel columns zero), e.g. external rows blended or swapped per target.
- `arrays` (frames × targets unless noted): `x` ⟨P,T⟩, `a` ⟨P,P⟩ over the held line's scored genes, `d` ⟨T,T⟩ (n),
  `pds` literal per-target PDS (the scorer's rank rule, ties ½), `pdsl` board-emulated PDS (414's β, λ 0.5), `K` ballast K
  (KMAX 16); with `base`: `base_x/base_a/base_pds/base_pdsl/base_K`, and the change energy |P_arm − P_base|² split
  `dE_scored`, `dE_unscored` (non-panel genes off the held line's axis), `dE_panel`, `dE_cpm0..3` (non-panel genes by the
  frame's A/B/C control CPM bin [0,5), [5,20), [20,100), [100,∞)), `dE_total`.
- `summary`: cos_agg (= mean over frames of Σx / √(Σa Σd)), pds, pds_lam0.5, pds_lam1, K_mean, AB, per_frame; with a base
  `vs_base` = r (cos_agg / cos_agg_base − 1), dpds, dpds_lam05, dK_mean, dK_rel, share_unscored, share_cpm.
- `E.cos_agg_arr(x, a, d, idx)` recomputes cos_agg on a target resample (the bootstrap identity, boot522).

## 4. Cloud runs (all through `experiments/800-checker/cloud/labcloud.sh`; on-demand; costs from `labcloud.sh wait`)
| run id | type | jobs | outcome | actual |
|---|---|---|---|---|
| `20261015T231934Z-lab801-folds-a` | r7i.8xlarge | verify, panel ×6 (CD4 + 3 fits + 3 dumps), refit-efp k562 | verify PASS (47/47 inputs sha-identical to the Mac); CD4 re-extraction byte-identical ×6; 18 fits done (460–620 s at 5 threads) but **every dump died on import** (guard801 looked for AGENTS.md, absent from the cloud tarball) and the operators died with the instance (they were uploaded only at the end of the job). efp refit: **sha `2a0247d2…` ≠ the Mac's `16ea9d34…`** | $1.31 |
| `20261015T235614Z-lab801-folds-b` | r7i.12xlarge | verify, prep ×6, fold ×18 | prep died at once (guard801 refused plolo's import-time `mkdir` of experiments/162-…/work, which exists on the Mac only); the 18 fold jobs waited 30 min for a prep marker that the failure path never wrote | $2.12 |
| `20261016T003830Z-lab801-folds-c` | r7i.16xlarge | same | prep died (nested import-time mkdir, 218's work/plolo); fail markers now written by an ERR trap, so the run ended in 6 min | $0.76 |
| `20261016T005128Z-lab801-folds-d` | r7i.16xlarge | same | prep died in `bind`: the synthetic bundle's frames computed on EC2 ≠ the A/B/C frame cache computed on the Mac (staged). **Frames (log1p of pseudobulk) are not bit-identical across Mac and EC2** | $0.66 |
| `20261016T010308Z-lab801-folds-e` | r7i.16xlarge | same, A/B/C frames now computed on EC2 (serialised by flock) | frames bound; the stub dump path reached its last line and died in `json.dumps` of the spec: **dump708.py's own record crashes on Path-valued fold_ops** (default=float), i.e. on its own `--fold-ops held` CLI path, after the npz is saved | $0.69 |
| `20261016T011500Z-lab801-folds-f` | r7i.16xlarge | same | prep ×6 PASS (incl. the stub dump path); 18 fits done, then **the instance role may not write `s3://…/lab-ops/`** (AccessDenied), so every fold job died before its dump; operators lost again | $2.71 |
| `20261016T015504Z-lab801-folds-g` | r7i.16xlarge | same, operators kept under `output/lab/ops801/` (a write probe in prep) | **DONE 0**: 18 fits (908–1169 s at 3 threads, peak 11.99–12.17 GiB) kept on S3, 18 dumps (18–52 s, peak 3.64–4.01 GiB), dump internal KA (per-donor mean == LODO key) passed on each | $1.98 |
| `20261016T022732Z-lab801-ka-h` | r7i.4xlarge | verify, ka-efp ×3, refit-study, smoke ×3 (L folds) | **DONE 0**: efp dumps + engine KA PASS (worst 1.6e-15); refit study; the L-fold smoke | $0.31 |
| **total (this mandate, `labcloud.sh wait` actuals)** | | | | **$10.55** |

Fixes, each in the file it concerns: guard801's repo check uses pyproject.toml; `GUARD801_ALLOW_MKDIR=1` on EC2 lets a NEW
empty directory through (file writes outside 801 still raise); prep writes a fail marker on every error path and runs the
whole dump path with an identity stub operator before any fit; fold jobs upload their operator to
`s3://…/output/lab/ops801/<L>/` the moment it exists and later runs fetch instead of re-fitting; prep writes an S3 probe
there first; dump801's record serialises Paths.

## 5. Known-answer table (abort unless |d| < 1e-12; every row passed)
| surface | instrument | rows | worst \|d\| | record |
|---|---|---|---|---|
| efp ×3 (522's dump522 files, Mac) | engine801 | 16 per fold: PT0 ×2, 546o ×5, 551 ×6, arms522 ×3 | **0.0** | `work/ka/local-efp-*.json` |
| ABC ×4 (518's dumps, Mac) | engine801 via fold522.Abc | 12 per fold: PT0 ×2, fold522.Abc own, 540c ×6, arms522 ×3 | **0.0** | `work/ka/local-abc-*.json` |
| efp ×3, fresh EC2 dump from 238's operators (run H) | dump801 `--registry --ka522` | 6 per fold: per-donor mean == LODO key (rel), 438 `h2c1\|V1_5_3+sh1.5` (realised + design), 546o `base q1.75` / `fixed 6` / `p96.6`, 540o `V1_0.5_6` | **1.6e-15** | `work/cloud/runH/efp-dump-*.log` |
| efp ×3, the same EC2 dumps | engine801 | the 16 rows above | **1.6e-15** | `work/ka/ka801-efp-*.json` |
| L1–L6 ×3 (EC2) | dump801 internal | per-donor mean == LODO key | < 1e-12 (asserted; the npz exists only if it passed) | `work/cloud/runG/L*-dump-*.log` |
| L1–L6 ×3 (EC2) | engine801 on load | PT0 ×2 (per-target path at constant q == scalar path, bit for bit) | **0** | `work/smoke/smoke801-L-*.json` (`ka_worst`) |
| 708's rehearsal row | — | **not present**: 708's rehearsal waiter is still armed behind the science chain and no row is stored | — | `experiments/708-def-panel-repricing/work/rehearse708.log` |

EC2 reproduces the Mac's stored rows to 1.6e-15 even though its frames differ from the Mac's in the last bits (run D) and
its operators would too (§0). So a fold's prices are platform-stable to ~1e-15. Bitwise equality holds only within one
platform.

## 6. Registry (`work/registry.json`, sha256 `8630992e1914853d02c11d9505d7f17c2e89d19dd2c206910aa068ff1dfc68de`)
Each entry: dump path + sha256 (checked against the EC2 record), operator (S3 URI or repo path) + sha256 + axis sha256 +
donors, target count, carried count, coverage, status. "carried" = targets that both of the fold's source donors carry
(the LODO key's carried set).
| fold | status | targets | carried | dump sha256 | operator sha256 |
|---|---|---|---|---|---|
| L1/k562_gwps | working | 269 | 259 | `11a35da67abec27d…` | `a97587a76a0920bb…` |
| L1/orion_hct116 | working | 291 | 259 | `9c943d16d22cbcd0…` | `9fca9af379af950c…` |
| L1/orion_hek293t | working | 293 | 259 | `02ea64e804034a15…` | `26bff5fa85a8e813…` |
| L2/k562_gwps | working | 269 | 259 | `b819a616d0eb20c9…` | `896c9aad5377f59c…` |
| L2/orion_hct116 | working | 290 | 259 | `73ae7edf2437515d…` | `fe452a39cf2e7d74…` |
| L2/orion_hek293t | working | 295 | 259 | `0e7d891474fb8acf…` | `2bc13b0f85887705…` |
| L3/k562_gwps | working | 269 | 259 | `bad4bd95d3166fae…` | `76ff2df001ee8c03…` |
| L3/orion_hct116 | working | 290 | 259 | `e734cd7de7ffb0dd…` | `5de5998446f876b5…` |
| L3/orion_hek293t | working | 294 | 259 | `97167b12ac137aef…` | `bf35cc9005bb0f9b…` |
| L4/k562_gwps | working | 269 | 259 | `c4d84778c0a005d7…` | `f345c46db8c641ae…` |
| L4/orion_hct116 | working | 290 | 259 | `6e389f8daa8983d4…` | `785068f375107acc…` |
| L4/orion_hek293t | working | 294 | 259 | `526b139cdd19843a…` | `1138a9c941589f44…` |
| L5/k562_gwps | working | 269 | 259 | `dfd49502a5ccc8c3…` | `b1875eb29c7f5d30…` |
| L5/orion_hct116 | working | 290 | 259 | `e0d4fbc91cc3c5c5…` | `6f523ee01821a237…` |
| L5/orion_hek293t | working | 294 | 259 | `9593b04268ebcbf1…` | `46bc1a32cbaebc36…` |
| L6/k562_gwps | working | 269 | 259 | `3f8b96fcc3422543…` | `074a7a8e05d21be5…` |
| L6/orion_hct116 | working | 291 | 259 | `49c0d4aad0eeea90…` | `2eec98fd1bf4833d…` |
| L6/orion_hek293t | working | 293 | 259 | `e72bd5db918fe033…` | `fdca7d400c02209c…` |
| ABC/k562_gwps | abc | — | — | `7fefe5640cd240ef…` | `13d988fcc3b129d2…` |
| ABC/orion_hct116 | abc | — | — | `e5f0581d79303769…` | `7cecbf7f36e40d22…` |
| ABC/orion_hek293t | abc | — | — | `c169e600dab1d0ec…` | `97e486e0d1a805e1…` |
| ABC/kolf2 | abc | — | — | `5e464e6dcc2c76b4…` | `e2f65ecf69b7e0b0…` |
| efp/k562_gwps | ka | 146 | 140 | `9d2cca5d8cfeca7b…` | `None…` |
| efp/orion_hct116 | ka | 293 | 140 | `46c2e722b3397b8c…` | `None…` |
| efp/orion_hek293t | ka | 296 | 140 | `b34ab6200dd7361d…` | `None…` |
| L7 | sealed-drawn | 300 | — | target list `97c30e7d11dddc01…` (`work/sealed/L7/pert_counts.csv`) | — |
| L8 | sealed-drawn | 300 | — | target list `32fd27d28425b12e…` (`work/sealed/L8/pert_counts.csv`) | — |

The registry above was regenerated after the efp dumps arrived (25 folds). The ABC folds' operator sha256 come from
plolo.OPERATOR_REGISTRY (062's held fits; kolf2 = production `e2f65ecf…`).

## 7. Pricing smoke (`work/smoke/<fold>-<arm>.npz`, `smoke801-{abc,L-k562_gwps,L-orion_hct116,L-orion_hek293t}.json`)
The arm is priced with the incumbent as its base, so the npz change-energy arrays are arm − incumbent. r = cos_agg(incumbent)
/ cos_agg(arm) − 1; Δpds = literal pds(incumbent) − pds(arm). "out" = the share of the change energy off the held line's
scored axis. A/B/C folds ran on the Mac (light); L folds ran on EC2 (run H). Fold-level only: no board units, no bins.

| fold | n | hi (incumbent's cap) | incumbent K | vs cap3: r | Δpds | out | vs pct+S_P2: r | Δpds |
|---|---|---|---|---|---|---|---|---|
| L1/k562_gwps | 269 | 13.23 | 5.60 | +19.86 % | +0.0040 | 0.23 | -0.14 % | +0.0007 |
| L2/k562_gwps | 269 | 8.84 | 6.40 | +17.40 % | +0.0149 | 0.22 | -0.41 % | +0.0038 |
| L3/k562_gwps | 269 | 7.82 | 6.01 | +13.25 % | +0.0098 | 0.23 | -0.35 % | +0.0024 |
| L4/k562_gwps | 269 | 15.01 | 5.95 | +34.70 % | +0.0189 | 0.23 | +0.37 % | +0.0023 |
| L5/k562_gwps | 269 | 11.70 | 5.65 | +19.72 % | +0.0105 | 0.24 | -0.54 % | -0.0005 |
| L6/k562_gwps | 269 | 9.17 | 5.99 | +23.14 % | +0.0132 | 0.24 | -0.62 % | +0.0003 |
| L1/orion_hct116 | 291 | 15.47 | 5.83 | +18.01 % | +0.0080 | 0.01 | +0.36 % | -0.0038 |
| L2/orion_hct116 | 290 | 14.09 | 5.98 | +14.94 % | +0.0136 | 0.01 | -0.29 % | +0.0052 |
| L3/orion_hct116 | 290 | 19.87 | 6.41 | +21.92 % | +0.0074 | 0.01 | +0.03 % | -0.0035 |
| L4/orion_hct116 | 290 | 14.07 | 5.58 | +17.29 % | +0.0134 | 0.01 | +0.37 % | +0.0000 |
| L5/orion_hct116 | 290 | 12.38 | 6.39 | +15.42 % | +0.0088 | 0.01 | -0.07 % | -0.0014 |
| L6/orion_hct116 | 291 | 17.18 | 6.15 | +21.39 % | +0.0091 | 0.01 | +0.35 % | -0.0050 |
| L1/orion_hek293t | 293 | 11.57 | 6.71 | +20.24 % | +0.0076 | 0.01 | +0.15 % | +0.0045 |
| L2/orion_hek293t | 295 | 12.13 | 6.65 | +15.64 % | +0.0069 | 0.01 | +0.56 % | -0.0067 |
| L3/orion_hek293t | 294 | 8.36 | 6.86 | +19.13 % | +0.0071 | 0.01 | +0.62 % | -0.0067 |
| L4/orion_hek293t | 294 | 11.62 | 6.87 | +20.96 % | +0.0161 | 0.01 | +0.03 % | -0.0050 |
| L5/orion_hek293t | 294 | 10.43 | 6.99 | +14.17 % | +0.0064 | 0.01 | +0.50 % | -0.0055 |
| L6/orion_hek293t | 293 | 12.89 | 6.94 | +26.65 % | +0.0109 | 0.01 | +0.06 % | -0.0008 |
| ABC/k562_gwps | 269 | 6.15 | 7.74 | +2.12 % | +0.0076 | 0.20 | +0.64 % | +0.0026 |
| ABC/orion_hct116 | 290 | 6.52 | 7.93 | +5.00 % | +0.0094 | 0.01 | +0.08 % | +0.0008 |
| ABC/orion_hek293t | 294 | 7.03 | 7.99 | +3.75 % | +0.0111 | 0.01 | +0.17 % | -0.0004 |
| ABC/kolf2 | 273 | 5.20 | 8.06 | +7.02 % | +0.0031 | 0.01 | -0.82 % | -0.0013 |

Reading (no aggregation, fold level only):
- vs cap3, the incumbent is ahead on all 22 folds in r and in Δpds.
- vs pct+S_P2 it is within ±0.6 % in r on every fold. Δpds is mixed: mostly + on the K562 folds, mostly − on the Orion folds.
- On the K562 folds, 22–26 % of the change energy lies off the scored axis (K562's 8.2k-gene axis). Under 702 rule 3 that is
  "partial-price" territory, so the checker should report it.

## 8. Deviations and notes
- **Local jobs while production waiters ran.** The draw, the diagnostics, the CD4 extraction, the closure discovery, the engine's
  Mac KA and the A/B/C smoke ran on the Mac. They were light (0.33–1.33 GiB, seconds each), but they ran while the production
  release watcher, the R0 waiter and an arm build were alive. Amendment 1's machine rule forbids any lab job then. The lab
  lead's resume note (22:58Z) allowed light local jobs outside the reserved window, and every one of them ran outside it.
  - The closure discovery's two Orion folds peaked at 2.25 / 2.41 GiB, over the 2 GB guideline, for about 6 s.
  - Inside the 23:30–01:45Z window the only local processes were file utilities: cp -c, shasum and aws s3 calls. There was
    also one 0.1 s `python -c ast.parse` syntax check of guard801.py at ~00:42Z.
- **Uploads.** < 1 GB from the Mac before 23:30Z. The efp operators (3.6 GB) went at 01:46–02:05Z, with vcc
  `pending_uploads` empty. Everything else came server-side from `backup-20261002/` and was sha256-verified on EC2 before use.
- **Nothing under `data/` was written.** The guard enforces this. Outside 801 the only writes were the lab runner's own:
  `experiments/800-checker/cloud/{runs.jsonl, runs/<rid>/, staged.tsv}` and `docs/LAB-LOG.md`'s cost rows, all appended by
  `labcloud.sh`.
- **Five failed launches (B–F) cost ~$7.** The cause was environment differences that the Mac tests could not show: a
  missing AGENTS.md, import-time mkdir, frame ulps, a Path in the json record, and S3 IAM. Each was fixed where it occurs, and
  prep now runs the full dump path with a stub operator plus an S3 write probe before any fit.
- **The A/B/C PLOLO folds load 518's dumps, which are Mac-made.** They are not re-dumped on EC2: their operators (062's held
  fits) are not on S3. The engine reads them at |d| 0 against 540c, arms540 and 522.
- **Sealed reserve.** L7/L8 exist only as target lists and sha256. No line of code computed a statistic on them; diag801
  hard-codes L1–L6. Their `pert_counts.csv` files make any later "every panel file on disk" exclusion exclude them.
  `docs/CONFIRMATION-LEDGER.md` was not touched; the lab lead records "Lab reserve LZ".

## 9. What is left
- The checker (`check800.py`, the lab lead) can read the folds through `engine801.load` + `price` (§3).
- Two jobs are left for the A/B/C panel:
  - An EC2 re-dump of A/B/C's own folds with dump801 under the D/E/F rules, so ABC is on the same footing as the L folds.
    It needs 062's held operators (3.9 GB) staged.
  - A Mac-fit check of an L-panel operator.
- dump708.py's json crash (§0) should go to the production coordinator before 10-22: a one-line fix, `default=str` for Path.
- The heavy-tail finding (§0) deserves its own look before L-fold key gains are read: compare A/B/C's tail anatomy (522 §4a)
  with L1–L6's.

## 10. Files
- Code: `draw801.py` (sha256 `c4b8507e…`), `diag801.py`, `cd4_801.py`, `guard801.py`, `lib801.py`, `dump801.py`, `frames801.py`,
  `discover801.py`, `engine801.py`, `test_local801.py`, `smoke801.py`, `refit801.py`, `registry801.py`, `cloud801.sh`, `.gitignore`.
- Records: `work/draws.json`, `work/dryrun801.txt`, `work/draw801.log`, `work/diag801.{json,txt}`, `work/registry.json`,
  `work/ka/*.json`, `work/smoke/*.json` (+ npz), `work/refit801-k562_gwps.json`, `work/closure/*.txt`, `work/manifest/*`
  (the Mac sha256 manifests the EC2 verify reads), `work/cloud/jobs-*.txt`, `sync-*.txt`, and `work/cloud/runG/`, `runH/`
  (every EC2 job log, fit log, closure and sha record).
- Data (gitignored): `work/dumps/dump708-{L1..L6,efp,efprefit}-*.npz` (1.9 GB), `work/L*/bundle/context_*.h5ad` (APFS clones),
  `work/L*/cd4/` (1.4 GB), `work/lib/`.

## 11. Follow-up (lab lead, 10-16 02:55Z): composition rows for row 710 (H1 ×2.5), per fold
`designs801.py` (new), run on EC2 through `cloud801.sh designs-L` / `designs-ABC`. Output `work/designs/<fold with / as __>.npz`
holds {`targets`, `fin_<design>` (float64, official 18,533-gene axis, panel columns zero; asserted per file)} plus a
`.json` per fold (npz sha256, per-design row sha256, KA, provenance). `work/designs/SHA256SUMS.txt` lists all 22 npz.
Everything is staged to the mirror: `s3://…/input/lab/repo/experiments/801-lab-folds/work/designs/` (44 files, copied
server-side from the run outputs).

### L1–L6 × 3 cancer folds, D/E/F rules: `fin_h2c1` (the incumbent) and `fin_h2.5c1`
- **Route:** `dump801.main(…, designs=["h2.5c1"])`, which is dump708's own extra-design path (`FC.final("h2.5c1")`). Each
  fold is re-run in full with its kept operator (mirrored from `output/lab/ops801/`; sha256 asserted equal to the
  registry's before the run).
- **KA: the re-run dump equals the registry dump on all 25 arrays** (fin, fPP, kv, V, truth, betas, frames, scored, …), worst
  |d| **0.0** on all 18 folds, **bit-identical**. That held on the same instance type (r7i.16xlarge) and the same 3 threads
  as run G. So `fin_h2c1` is byte for byte the registry dump's `fin`.
- `fin_h2.5c1` differs from `fin_h2c1` on every fold (row sha256 differ).

### A/B/C × 4 folds, A/B/C rules: `fin_h2c1.5` (ALLSW's stack) and `fin_h2.5c1.5`
- **518's `fST` is h2c1.5.** dump518 builds it as `compose_arm(…, m_cd4=1.5, g_h1=2 × UNION_G)`, with `COMPS["stack"]` =
  (1.5, 0.56) in lib420, i.e. H1 ×2 (VCC_H1_UNION_G 0.56) and CD4 ×1.5. So I went ahead.
- **Route:** lib801 bound to `"abc"` (438's abc binding). That means data/challenge, 062's held operators (production
  `e2f65ecf…` for kolf2), A/B/C's CD4 banks and SVD-5 deflation. Then build_fold438, Donors438, residuals, pre_normrel,
  op.apply, Finals, and `FC.final("h2c1.5")` / `FC.final("h2.5c1.5")`. FC.meta confirms (m_h1, m_cd4) = (2, 1.5) / (2.5, 1.5);
  2.5 × 0.28 = 0.70 = `VCC_H1_UNION_G`.

| fold | n | operator | h2c1.5 vs 518 fST, max \|d\| | h1c1 vs 518 fPP | targets |
|---|---|---|---|---|---|
| ABC/k562_gwps | 269 | 062 held_k562_gwps | 8.8e-15 | 7.0e-15 | equal |
| ABC/orion_hct116 | 290 | 062 held_orion_hct116 | 7.6e-15 | 8.4e-15 | equal |
| ABC/orion_hek293t | 294 | 062 held_orion_hek293t | 3.4e-14 | 3.8e-14 | equal |
| ABC/kolf2 | 273 | production union_pruned | 8.7e-15 | 9.7e-15 | equal |

- These rows are not bitwise equal: 518's dumps were made on the Mac, these on EC2. All are inside the 1e-12 tolerance.
- Every input the run read was sha256-verified against the Mac first (109 files). The newly uploaded ones are 062's three
  held operators + jsons (3.6 GB, 03:09–03:35Z, window CLEAR, vcc `pending_uploads` empty) and the kolf2 bank (30 MB). A/B/C's
  CD4 banks and the production operator were copied server-side from `backup-20261002/` and verified.

### Runs and cost (this follow-up; `labcloud.sh wait` actuals)
| run id | type | jobs | outcome | actual |
|---|---|---|---|---|
| `20261016T031543Z-lab801-designs-i` | r7i.16xlarge | — | terminated by me 1 min in: its sync list lacked two verified inputs, so verify would have failed | $0.02 |
| `20261016T031833Z-lab801-designs-j` | r7i.16xlarge | verify, prep ×6, designs-L ×18 | DONE 0, 10 min | $0.72 |
| `20261016T034108Z-lab801-designs-k` | r7i.2xlarge | verify, designs-ABC ×4 (serialised on the frame-cache lock) | DONE 0, 10 min | $0.09 |
| **total** | | | | **$0.83** (budget $15); this mandate overall **$11.38** |

### Notes
- **Local work in this follow-up:** file reads and shasum only. I queued one local stub-closure run through gate602 at
  02:58Z; it never started (the gate held it on swap free 0 MB) and I killed it at 03:10Z. Its closure was found on EC2
  instead.
- `.gitignore` already covers work/ npz/h5ad (written 10-15).
- No instance of mine is running (`labcloud.sh status`).

## 12. Follow-up 2 (lab lead, 10-16 08:05Z): is the A/B/C-vs-L cap gap the panel or the operator recipe? (`ABCp`)
**Answer: about one third operator recipe, two thirds panel.**
- Refitting A/B/C's held-out operators by the DAY recipe raises the fold's percentile-cap hi:
  - k562 6.15 → 6.18;
  - hct116 6.52 → 7.84;
  - hek293t 7.03 → 8.39.
- The L folds stay well above even the day-recipe A/B/C values:
  - k562 7.8–15.0 (median 10.4);
  - hct116 12.4–19.9 (median 14.8);
  - hek293t 8.4–12.9 (median 11.6).
- So 062's recipe (all donor rows, essentials kept, no cap) does compress A/B/C's tail on the two Orion folds. The panel
  still accounts for most of the A/B/C-vs-L gap. The day-recipe A/B/C folds are now on the same instrument as the L folds and
  as the 10-22 day path.

### What was built (all on EC2; code `dump518p.py`, `cloud801.sh fit-abc / abcp`, `engine801.py` `ABCp` loader)
1. **Day-recipe operators.** `scripts/build_union_whitener.py --bundle data/challenge --donors <the other two>
   --basis-prune essential --max-rss-gib 14` (12k per-donor cap; the L folds' flags), r7i.8xlarge, 8 threads, 322–338 s
   each. Kept at `s3://…/output/lab/ops801/ABCp/` (and mirrored under `input/lab/repo/experiments/801-lab-folds/work/ABCp/ops/`).
   Their axes equal 062's registered axes (`6eaba848…` k562 fold, `edf6ba88…` Orion folds); only the basis differs.

   | fold | day operator sha256 | basis rows kept (per source) | 062 operator sha256 |
   |---|---|---|---|
   | k562_gwps | `48cb88b68764576479516a5215b4f710dd79a1815d3536611da88c1adaa646e1` | hct116 10,807 + hek293t 10,803 | `13d988fc…` |
   | orion_hct116 | `a0e27947cfd72fcd6de9ec29a69759e32c946f7bfbf0e698dd472921c0583318` | k562 7,784 + hek293t 10,776 | `7cecbf7f…` |
   | orion_hek293t | `5d982462bf45380e43af7664e7b6c4249c0218930ed753c61ced339cfc183db7` | k562 7,784 + hct116 10,807 | `97e486e0…` |
2. **The folds, under the board's A/B/C rules.** `dump518p.py` is a copy of dump518.py (blob `d5096113…`): lib420's path,
   CD4 ×1.5 all-rows sizing, SVD-5 re-fitted, H1 ×2. The operator is swapped in-process (`plolo.FOLDS[held]`).
   - **KA first, with 062's operators.** The rebuild reproduces 518's stored dump on EVERY array (fPP, fST, kv1, carried,
     truth, betas, …): worst |d| 8.3e-14 / 8.4e-15 / 6.4e-14. dump518's own KAs also pass: KA-PP; V1 vs 418's `V1__PP` at
     3.1e-14 / 1.5e-14 / 1.6e-14 (not bitwise, Mac → EC2); arms540 V1_0.5_3 and V1_0.5_6 at |d| 0.0. The fold's hi
     reproduces too (6.1494 / 6.5208 / 7.0282).
   - **Then with the day operators** → `work/dumps/abcp/dump518-<held>.npz`:
     - k562_gwps `169279402c589a38…`;
     - orion_hct116 `10f5860ce7191307…`;
     - orion_hek293t `dc7f4df69cd26fc9…`.
   - Full sha256 and KA in `work/dumps/{abc062,abcp}/dump518p-*.json`.
   - The day ops go through plolo's `production` operator branch (self-described json sha256, axis re-derived and checked).
     plolo's `held` branch wants 062's metadata format, which a build_union_whitener json does not carry.
3. **Engine and registry.** `engine801.load("ABCp/<held>")` uses fold522.Abc's route (score518.Fold) on the new dump.
   `F.panel` stays `"ABC"`, so ka804 applies its A/B/C arm set: the incumbent is cap6+P2g(abc), and the K562 fold gets
   PIE_ABC. The stored-row KAs are at the dump stage; on load only PT0 runs. `work/registry.json` (sha256 `ce2ced4c…`) now
   lists `ABCp/<held>` (status `abc-dayops`) with dump and operator sha256, both hi values and the recipe.
4. **ka804's arm set on ABCp** (`experiments/804-known-answers/ka804.py`, unmodified, `KA804_OUT` outside the repo).
   Standard-format arrays are in `work/ka-abcp/{bulk,meta}/` (1.8 MB) and on S3 at
   `input/lab/repo/experiments/801-lab-folds/work/ka-abcp/`. For the lab lead's regrade, a look at r vs each arm's base:
   - cap6 vs cap3: +1.52 / +3.35 / +3.66 %;
   - 707: +2.03 / +4.02 / +3.45 %;
   - 705: +0.41 / +0.13 / −0.17 %;
   - 702 (K562 only): −1.65 %.

### Percentile-cap hi per fold (P 96.64697522769474 of key/median over carried)
| fold | A/B/C, 062 operator | A/B/C, day operator (ABCp) | L1 | L2 | L3 | L4 | L5 | L6 |
|---|---|---|---|---|---|---|---|---|
| k562_gwps | 6.15 | 6.18 | 13.23 | 8.84 | 7.82 | 15.01 | 11.70 | 9.17 |
| orion_hct116 | 6.52 | 7.84 | 15.47 | 14.09 | 19.87 | 14.07 | 12.38 | 17.18 |
| orion_hek293t | 7.03 | 8.39 | 11.57 | 12.13 | 8.36 | 11.62 | 10.43 | 12.89 |

### Runs (this follow-up; actuals)
| run | outcome | $ |
|---|---|---|
| `20261016T080732Z-lab801-abcp-a` | 3 day fits DONE (kept on S3); dumps died: sens415's `full-*.npz` caches absent on EC2 | 0.70 |
| `…083058Z-…-b` | fits reused; 415's `prep-*.npz` absent (staged 415/418/420 work dirs after this) | 0.18 |
| `…084540Z-…-c` | dump518's KA V1 is a bitwise assert; it differs in the last bits Mac → EC2, so it is now |d| < 1e-12 (recorded) | 0.13 |
| `…085444Z-…-d` | 062 KA PASS; the day ops were refused by plolo's `held` branch (062 metadata format), so they go through the `production` branch | 0.15 |
| `20261016T090515Z-lab801-abcp-e` | **DONE 0**: 062 KA PASS ×3, ABCp dumps ×3, ka804 ×3 | 0.16 |
| **total** | | **$1.32** (budget $15) |

Uploaded from the Mac (window CLEAR, no pending vcc uploads): 415/418/420's work dirs (671 MB, the lib420 path's caches),
418's arrays, 207/169/410 small records. No instance of mine is running.

## 13. Follow-up 3 (lab lead, 10-16 09:31Z): 710's compositions on the ABCp folds
The ask: 062's `fin_h2c1.5` differed from the ABCp fold's `fin` (the compositions live in the operator's whitened space), so
row 710 needed designs built with the DAY operator.

**Delivered:** `work/designs/ABCp__<held>.npz` with keys `targets`, `fin_h2c1.5`, `fin_h2.5c1.5` (same keys and targets as
ABC__; panel columns zero, asserted), plus `.json` per fold.
- Code: designs801.py's ABC path (lib801 bound to "abc", A/B/C rules, `Finals.final`), new modes `ABC062check` and `ABCp`.
  `ABCp` swaps the fold's held-out operator for the day operator (`work/ABCp/ops/held_<held>.npz`), the same file
  dump518p.py's ABCp fold used. That is done through the bind spec's `fold_ops` under its own cache tag `abcp`.
- Mirrored to `s3://…/input/lab/repo/experiments/801-lab-folds/work/designs/` (6 files).

| fold | npz sha256 | row sha256 h2c1.5 / h2.5c1.5 | operator sha256 (day) |
|---|---|---|---|
| ABCp/k562_gwps | `dc0ced1283a75768518f4445725849bbab9e0dd95b7083b61cf474189cbff848` | `e39443d6…` / `7eab4d7b…` | `48cb88b6…` |
| ABCp/orion_hct116 | `4e9680c49d8a4625409821e74112aebe9da95fbb273213ff25be41c291cb4e65` | `44349419…` / `1a35ad82…` | `a0e27947…` |
| ABCp/orion_hek293t | `78bd43a28c8c7266d42576e5f0595d07ad446455bc82f2d635e2a118caf56892` | `c474c9be…` / `8e091231…` | `5d982462…` |

**Known answers (all pass):**
- `h2c1.5` == the ABCp fold's `fin` (dump518p's `fST` in `work/dumps/abcp/dump518-<held>.npz`, which is what
  `engine801.load("ABCp/<held>").fin` returns): |d| = **0.0, bit-identical**, ×3. `h1c1` == its `fPP`: 0.0. Targets equal.
- With 062's operator, the same code reproduces the ABC__ files' row sha256 for `h2c1.5` (and `h2.5c1.5`): **match ×3**
  (`work/designs/check/ABC062__<held>.json`; same instance type r7i.2xlarge and 4 threads as the ABC__ run).
- FC.meta: (m_h1, m_cd4) = (2, 1.5) / (2.5, 1.5).

**Runs:**
- `20261016T092938Z-lab801-designs-abcp`: terminated by me at once. Its sync list lacked the 415 caches that the verify
  manifest now lists. $0.02.
- `20261016T093113Z-lab801-designs-abcp2`: DONE 0 in 9 min. $0.08.
- Total **$0.10** (budget $2). No instance of mine is running.

**Local work under the hold:** cp -c, file copies, json reads. Plus one 0.1 s `python -c ast.parse` syntax check of
designs801.py at ~09:29Z, which strictly counts as a local job during the R7a hold.

I did not run price710 / price710p.

## 14. Follow-up 4 (lab lead, 10-16 10:40Z): A/B/C's panel through the L-fold recipe (`LABC`); is it the targets or the rules?
**Answer: the panel's targets.**
- On the key arm, LABC sits with ABCp/ABC062 at the bottom of each line, far below all six L folds. 495 (cap6 vs cap3), r:

  | held line | LABC | ABCp | ABC062 | L1–L6 range |
  |---|---|---|---|---|
  | k562 | +2.03 % | +1.52 % | +0.82 % | +8.4 … +14.0 % |
  | hct116 | +3.71 % | +3.35 % | +3.84 % | +5.0 … +9.2 % |
  | hek293t | +4.19 % | +3.66 % | +3.20 % | +7.4 … +13.6 % |

- A/B/C's fold rules (CD4 ×1.5 all-rows, SVD-5, H1 ×2 stack) move it by only ≈ +0.4 … +0.5 pt (ABCp → LABC). Pct vs cap3 is
  +2.0 / +3.6 / +4.0 %, about the same as cap6, so the fold cap hi on LABC stays near 6, as on A/B/C.
- **The truth-strength profile does not name the property.** The three A/B/C fold variants share one truth (same targets,
  same held line), so their profile is identical. That profile sits inside the L spread on every statistic:
  - median d: k562 46.1 vs L 39.7–44.2; hct116 50.3 vs L 49.4–53.5; hek293t 37.8 vs L 37.1–40.6;
  - top-30 share: 0.292 / 0.310 / 0.309 vs L 0.304–0.363 / 0.293–0.341 / 0.281–0.340;
  - n > 4 × median: 7 / 11 / 16 vs L 7–18 / 7–15 / 6–18.
- Spearman with the 495 gain, pooled over 27 folds: median d −0.38 (p 0.05), mean d −0.33, top-30 share +0.30, n > 4 × median
  +0.20, mean/median +0.21. Within a line the signs disagree (k562: median d −0.86, mean/median +0.76; hct116: mean d −0.80;
  hek293t: all |ρ| ≤ 0.20).
- So the truth's strength distribution does not separate A/B/C from the L folds. What differs is the source-side key: A/B/C's
  LODO key has the light tail, with cap hi ≈ 6 vs 8–20 on L. The key gain needs a heavy-tailed key whose tail the held truth
  supports. That points at how A/B/C's targets look in the two SOURCE donors, not at how strongly the held line responds.
  This reading is mine; the profile itself is in `work/profile801.{txt,json}`.

### What was built (EC2 run `20261016T104243Z-lab801-labc`, r7i.4xlarge, DONE 0, 11 min, **$0.17**; budget $6)
- **LABC/<held>:**
  - `dump801.py --tag LABC --bundle data/challenge` (D/E/F rules, the L-fold recipe), with `--cd4-root work/LABC/cd4/cd4_compact`
    and `--fold-ops held --held-ops work/ABCp/ops`.
  - The panel is A/B/C's 300 targets; the held truths are the same.
  - The CD4 banks come from `cd4_801.py LABC`, the L folds' route, and are byte-identical on every non-json file to A/B/C's
    production banks (`data/external/cd4t/cloud-compact/cd4_compact`).
- **The operator is not fresh: the day-recipe fit for this panel was reused.** These are ABCp's
  `build_union_whitener.py --bundle data/challenge --basis-prune essential` (12k cap) fits from 08:13Z, the same command and
  panel. It keeps LABC vs ABCp a pure rules contrast. A refit would differ only in its last bits (see §0: a refit moves
  prices ≤ 1.4e-11). Shas `48cb88b6…` / `a0e27947…` / `5d982462…`.
- **Dumps:**

  | fold | dump sha256 |
  |---|---|
  | LABC/k562_gwps | `eba095c118837ae9a881f65d8dca1906d491733ac7d200b0427757963e0d20d3` |
  | LABC/orion_hct116 | `e97f1fbd6bd86626533bf733c5d69aef3de2f5bf199943c03c0bc7e90dfb8337` |
  | LABC/orion_hek293t | `34ecd80e26bff4dc147de8d23326181898d3d0bb17381758134ed4638e60fb95` |

- **KAs as for the L folds:** the dump's per-donor mean == LODO key, max rel |d| 3.4e-15 / 1.7e-15 / 6.8e-16. engine801 PT0
  ×2 = 0 on load.
- **Registry:** `work/registry.json` (sha256 `d4ee481c…`, 31 folds) has `LABC/<held>` with status `labc`.
- **ka804, unmodified, L arm set (inc = pct+P2g)** → `work/ka-labc/{bulk,meta}/`, mirrored to
  `input/lab/repo/experiments/801-lab-folds/work/ka-labc/` (34 objects).
  - **Caveat:** on LABC/k562 the 702 arm is a no-op. ka804 picks PIE_L (L panels) when `panel != "ABC"`, and PIE_L carries no
    A/B/C targets (n_seen 0, r 0.00 %). The other arms are unaffected.
  - Arms vs their bases, r on k562 / hct116 / hek293t:
    - 705: +0.55 / +0.04 / −0.38 %;
    - 707: +1.49 / +3.62 / +2.74 %;
    - pct+P2g vs pct: +1.64 / +0.85 / +0.13 %.
- **Truth-strength profile:** `profile801.py`, run on EC2 from 804's standard arrays (L, ABC062), `ka-abcp` and `ka-labc`.

No local job ran under the hold: only cp, aws s3 calls, file reads, and one sub-second json edit of the registry. No
instance of mine is running.

## 15. Follow-up 5 (lab lead, 10-16 11:49Z): the "A-only" draws A1–A6 and their 18 folds
### Draws (EC2 run `20261016T115240Z-lab801-draw-a`, written before any fold existed)
- **Rule (`draw801a.py`):** draw801's Rule class, imported unchanged, plus the ONE pool change: also remove
  `data/donors/rpe1/perts.txt` ∪ `data/donors/nadig_jurkat/perts.txt`.
  - The two lists' sha256 are `4377d2dc…` and `f5a5c1bb…`; their union is 2,351 genes.
  - **HepG2 list: absent.** There is none under data/donors, and 814 names none. LOG l.200–212 lists Nadig HepG2 at 2,394,
    the same DepMap set and count as Jurkat's.
- **Independent draws from one pool:** seeds 801301–801306. The exclusion is every `pert_counts.csv` on disk except 801's
  own L1–L6 bundles: 29 files, 4,973 targets.
  - The file list was made on the Mac and is sha-checked on EC2: `work/manifest/panel_files_801a.sha256`.
  - It includes L7/L8 (sealed), LABC/A/B/C, R6, R7a and the rest.
- **Pool and quotas:** pool 4,355 (4,525 before the screen cut). Capacity is each draw's own need (no disjointness), so no
  cell was short and **no moves** were needed. Every panel has A/B/C's exact strata: Q1 57, Q2 28, Q3 84, Q4 57, Q5 57, NA 17;
  C3 259, K2 10, O2 30, O1 1. K562 269, CD4 283.
- **Asserted at draw time:** 0 K562-essential, 0 essential-screen, 0 overlap with the exclusion (L7/L8 included), every quota
  met.
- **Record:** `work/draws801a.json` sha256 `d4565ad28185f16792e12bce4627b761e6c28e8fd2d283a547d9eef50314c844` (per draw: sha256,
  overlaps, `describe` = strength bins, coverage, basal CPM quantiles).

| panel | seed | sha256 | overlap L1..L6 | A/B/C | earlier A | HCT / HEK coverage | basal CPM q5 |
|---|---|---|---|---|---|---|---|
| A1 | 801301 | `10ee7d1845ec57e8…` | 29 22 31 36 31 37 | 0 | — | 289 / 293 | 17.27 |
| A2 | 801302 | `406ed30b29b92a36…` | 31 33 31 29 32 32 | 0 | A1 43 | 290 / 294 | 16.91 |
| A3 | 801303 | `2781c64279935a76…` | 37 30 30 26 22 35 | 0 | 41 44 | 291 / 295 | 17.64 |
| A4 | 801304 | `56c9ad5f91b50bf3…` | 33 32 26 31 29 34 | 0 | 33 32 34 | 290 / 295 | 16.81 |
| A5 | 801305 | `ae4b0f5717c3508f…` | 30 31 44 26 32 25 | 0 | 35 37 30 25 | 290 / 294 | 18.53 |
| A6 | 801306 | `f1188d7b8a9ff89a…` | 32 29 22 32 29 33 | 0 | 25 30 36 28 25 | 291 / 292 | 17.94 |

### Folds (EC2 run `20261016T120234Z-lab801-aonly`, r7i.16xlarge, DONE 0, 37 min)
- **The L-fold recipe, unchanged.** Per panel: CD4 banks via `cd4_801.py`, frames, the stub dump-path check, then an S3 write
  probe. Per fold: a fresh day-recipe held-out fit (build_union_whitener `--basis-prune essential`, 12k cap, 3 threads,
  1,022–1,222 s, kept at `s3://…/output/lab/ops801/A<k>/` the moment it existed), then `dump801`.
- **CD4 banks have no Mac reference.** Under the local hold no Mac extraction was possible. Their per-file sha256 are recorded
  (`work/cloud/runAonly/A<k>-cd4.sha256`); the L folds had already shown EC2 == Mac on 6 panels.
- **KAs:**
  - dump per-donor mean == LODO key, max rel |d| ≤ 2.9e-14 on all 18;
  - engine801 PT0 ×2 = 0 on load (ka804's meta).
- **Registry:** `work/registry.json` (sha256 `08807731…`, 49 folds) has `A<k>/<held>` with status `a-only`, plus the dump and
  operator sha256, axis sha256, donors and target counts.
- **ka804, unmodified, L arm set (inc = pct+P2g):** `KA804_OUT=/work/ka804-<fold>` → `work/ka-aonly/{bulk,meta}/` (11 MB),
  mirrored to `input/lab/repo/experiments/801-lab-folds/work/ka-aonly/` (204 objects). The dumps are mirrored too.
  - 702 note: PIE_L carries 67–82 of each A panel's K562-fold targets, through overlap with L1–L6, so 702 is a PARTIAL arm
    on A-folds.
- **A first look, for orientation only (the lead grades with panelrank800 / overread800).** 495 cap6 vs cap3, r on A folds:

  | held line | A folds | ABCp | LABC |
  |---|---|---|---|
  | k562 | +3.8 … +9.9 % | +1.5 % | +2.0 % |
  | hct116 | +3.0 … +7.1 % | +3.4 % | +3.7 % |
  | hek293t | +9.9 … +13.5 % | +3.7 % | +4.2 % |

### Cost
| run | $ |
|---|---|
| draw-a (r7i.xlarge, 7 min) | 0.03 |
| aonly (r7i.16xlarge, 37 min; launch cap $6.04, actual below the $5 budget) | 2.65 |
| **total** | **$2.68** |

No local job ran: only cp, aws, find/shasum file reads and a sub-second registry json edit. No instance of mine is
running.

## 16. Follow-up 6 (lab lead, 10-16 12:01Z; science 713 part 2): LABC with H1 ×2.5 / ×3 designs (`LABCd`)
**Done.** The three LABC dumps were rebuilt under the new tag `LABCd`: `dump801.py --tag LABCd --designs h2.5c1,h3c1`, with
follow-up 4's spec unchanged (bundle data/challenge, D/E/F rules, ABCp's day operators, the LABC CD4 root re-extracted). The
re-extraction is byte-identical to A/B/C's production banks again. The originals are untouched.

| fold | new dump sha256 (`work/dumps/dump708-LABCd-<held>.npz`) | json sha256 | KA vs current LABC dump |
|---|---|---|---|
| k562_gwps | `ba5433a259efcc0b006110197f60cfe15a2b1d660e1d4d0ffac911856ee0440f` | `1479b0b4…` | every array **bit-identical** (vs `eba095c1…`) |
| orion_hct116 | `f48b21de135abda5714e9d66947c98f952028667a4dc4bae1a7e35572b202da4` | `7b6f5ed7…` | every array **bit-identical** (vs `e97f1fbd…`) |
| orion_hek293t | `bc524dc5b88d49a9563e4c37eafc3f1e7f8b6a402aae6f50e715d2de88c5634d` | `e49c5283…` | every array **bit-identical** (vs `34ecd80e…`) |

- **KA:** every non-design array of the new dump equals the current LABC dump byte for byte; fin == the original fin
  (h2c1), |d| = 0. The only extra arrays are `fin_h2.5c1` and `fin_h3c1`, and each differs from fin.
- **Design row sha256** (float64 bytes):

  | fold | fin_h2.5c1 | fin_h3c1 |
  |---|---|---|
  | k562 | `e696e2e4…` | `aa9dc70c…` |
  | hct116 | `83cf21a2…` | `93edbbcb…` |
  | hek293t | `97538520…` | `abb71b17…` |

  KA records: `work/cloud/runLABCd/LABCd-ka-<held>.json`.
- **Format and mirror:** the dump json records `designs: ["h2.5c1", "h3c1"]`. Mirrored to
  `s3://…/input/lab/repo/experiments/801-lab-folds/work/dumps/` (6 files).
- **Run:** `20261016T124527Z-lab801-labcd` (r7i.4xlarge, the same type and 3 threads as follow-up 4), DONE 0, 7 min,
  **$0.13** (budget $1).
- I did not price the design arms. No local job ran, and no instance of mine is running.

## 17. Follow-up 7 (lab lead, 10-16 13:22Z): the panel kit for 10-22, and its known-answer dry run
### The command
```
bash experiments/801-lab-folds/kit/kit801.sh <TAG> <pert_counts.csv> <gene_names.csv> [--ref REF] [--no-wait]
bash experiments/801-lab-folds/kit/kit801-collect.sh <TAG> <RUN_ID>        # only after --no-wait
# 10-22:  bash experiments/801-lab-folds/kit/kit801.sh DEF data/challenge-def/pert_counts.csv data/challenge-def/gene_names.csv
```
**Mac side (`kit801.sh`):** shell only (cp, cmp, aws, the launcher), so it runs under the local-job hold. It reads the panel
files and writes nothing beside them; `data/challenge-def/` is untouched.
- Checks: the tag is new and not reserved (L1–L8, A1–A6, LABC/LABCd, ABC, ABCp, efp, efr, real, sealed are refused, and an
  existing `work/<TAG>` stops it); `pert_counts.csv` has a `target_gene` column and 300 targets; `gene_names.csv` is
  byte-identical to `data/challenge/gene_names.csv` (the engine's official axis), otherwise it stops and asks.
- Then: copies the two csv files to `work/<TAG>/bundle/`, stages them, writes the job list and launches one r7i.8xlarge (cap
  60 min).

**EC2 side (`cloud801.sh` subcommands, each already exercised on L/A/LABC):**
- `verify`.
- `prep <TAG>`: manifest (A/B/C's, relabelled D/E/F), the relabelled A/B/C controls (frames as on every lab fold), CD4 banks
  via `cd4_801.py`, the frame cache, the stub dump path, an S3 probe.
- `fold <TAG> <held>` ×3: a FRESH day-recipe held-out operator (`build_union_whitener.py --basis-prune essential`, 12k cap,
  `FIT_THREADS=8`), kept at `s3://…/output/lab/ops801/<TAG>/`, then `dump801` (D/E/F rules, `T801_THREADS=3`).
- `ka804 <TAG> <held>` ×3 with `KA_DIR=<tag>`: ka804 unmodified, the L arm set (inc = pct+P2g) → `work/ka-<tag>/{bulk,meta}`.
- `regfrag <TAG>`: `work/registry-<TAG>.json`, status `kit`.
- With `--ref`: `kacmp <TAG> <REF>`.

**Collect (`kit801-collect.sh`):** `labcloud wait`; outputs copied into work/ (APFS clones); `ka-<tag>`, the dumps, the
fragment and the operators mirrored server-side to `s3://…/input/lab/repo/experiments/801-lab-folds/work/`; the fragment
merged into `work/registry.json` when window800 prints CLEAR (kept as a fragment otherwise).

**Files** (`kit/`; sha256 prefixes at writing):
- `kit801.sh` (`df4f01f0`), `kit801-collect.sh` (`fcd52585`);
- `kit810.sh` (`e812d71d`), `kit811ext.py` (`2a545547`), `kit811.py` (`58eae10f`), `run-kit811.sh` (`9ccd12f1`);
- `make-jobs.sh` (`8824b2cf`), `sync-base.txt`;
- `cloud801.sh` gained `tag_ok`, `regfrag`, `kacmp`, the kit manifest step, `FIT_THREADS` and `KA_DIR`. `cd4_801.py` and
  `frames801.py` now accept any non-sealed tag.

### Dry run: TAG = LABCk on `data/challenge`'s files, `--ref LABC` (run `20261016T132753Z-lab801-kit-labck`, DONE 0)
| held | fresh operator == ABCp's (LABC's) | dump vs `dump708-LABC-<held>.npz` | ka804 arms vs `ka-labc` |
|---|---|---|---|
| k562_gwps | **yes**, `48cb88b6…` bit-identical | **bit-identical**, every array (vs `eba095c1…`) | 11 arms, worst \|d\| 2.3e-13 |
| orion_hct116 | **yes**, `a0e27947…` | **bit-identical** (vs `e97f1fbd…`) | 10 arms, worst 2.3e-13 |
| orion_hek293t | **yes**, `5d982462…` | **bit-identical** (vs `34ecd80e…`) | 10 arms, worst 2.8e-13 |

- **The fresh fits reproduce LABC's operators bit for bit.** The kit fits on the same instance type (r7i.8xlarge) with the
  same 8 threads as ABCp's fits, which shows a day-recipe fit is deterministic given (instance type, threads). The dumps are
  therefore byte-identical too; dump KA per-donor mean == LODO key at 3.4e-15 / 1.7e-15 / 6.8e-16.
- **The ka804 arrays differ by ≤ 2.8e-13 (< 1e-12, not bitwise).** Cause: the pricing threads. The kit's ka804 runs at 2
  threads, matching the A-only and 804's L runs; LABC's ran at 3. Every number agrees to ~1e-15 relative.
- The registry has `LABCk/<held>` (52 folds). Arrays are in `work/ka-labck/` and mirrored. The KA record is
  `work/cloud/kit-LABCk/kacmp-LABCk-vs-LABC.json`.

### Runtime and cost per panel
16 min from launch to DONE:
- ~6 min instance and venv setup, plus the 113-file input verify;
- ~5 min for the 3 fits in parallel (295–314 s);
- ~0.5 min of dumps;
- ~1.5 min of ka804.

**$0.51 per panel** (r7i.8xlarge, actual). On 10-22, add ~1 min for `kit801-collect.sh`.

### The 810 / 811 job files the kit feeds (for the lead to launch through 810 / 811)
`bash experiments/801-lab-folds/kit/make-jobs.sh <TAG>` writes `kit/jobs810-<TAG>.txt`, `kit/sync810-<TAG>.txt`,
`kit/jobs811-<TAG>.txt` and `kit/sync811-<TAG>.txt`. The script's header carries the two launch commands; the sync files
are one path per line, fed as `$(sed 's/^/--sync /' file)`. They are written for LABCk as the example.
- **810 (sgu1, 707): `kit810.sh`.** A copy of `job810aonly.sh` (blob `c9586229…`) with the tag as a parameter. Operators come
  from 801's mirror of the kit's fits (`work/<TAG>/ops`) and are checked against `work/registry-<TAG>.json`. Outputs go to
  `810/work/kit-<TAG>/`. price810 / decomp810 are called unchanged; they key only on `F.panel == "ABC"`.
- **811 (c0_l25, c1_l25, src_l25): `run-kit811.sh <TAG>`.**
  - `kit811ext.py` is a copy of ext811a (blob `55fbda14…`) with `KIT811_PANELS`, keeping its L1 CD4 known answer.
  - `kit811.py` is a copy of aonly811 (blob `873db6c1…`); only the fold assert is opened to the kit tags, and it keeps the
    L1/k562 known-answer fold.
- **Not run here:** these job files were not launched (budget and remit). Their own KAs gate them on the day: 810's decomp
  `parts sum == fin`; 811's L1 CD4 and the L1/k562 arm reproduction.

**Notes for the day:**
- The kit's frames are A/B/C's relabelled controls, as on every lab fold. A read on D/E/F's own control frames would need a
  bundle with those contexts (708's route); the kit does not do that.
- 702 is a no-op on a fresh panel (PIE_L holds none of its targets), as on LABC.

Cost of this follow-up: **$0.51** (budget $1). No local job ran (shell file operations, aws, one sub-second registry json
merge with window800 CLEAR). No instance of mine is running.

## 18. Follow-up 8 (lab lead, 10-17 10:18Z): `--frames <bundle dir>` for the kit (D/E/F's own control frames on 10-22)
### The option
```
bash experiments/801-lab-folds/kit/kit801.sh <TAG> <pert_counts.csv> <gene_names.csv> --frames <bundle dir> [--ref REF] [--no-wait]
# 10-22:  bash experiments/801-lab-folds/kit/kit801.sh DEFo data/challenge-def/pert_counts.csv data/challenge-def/gene_names.csv --frames data/challenge-def
```
- **Mac side (`kit801.sh`, shell + jq only).** `<bundle dir>` must be repo-relative, with a `manifest.json` that lists 3
  contexts.
  - The kit writes `work/<TAG>/bundle/frames_src.txt` (the dir) and `frames.sha256` (manifest + the three
    `context_<c>.h5ad`).
  - It stages each of those files to the S3 mirror only where the mirror lacks a same-size copy (on 10-22 that is ~0.66 GB,
    under the 1 GB Mac-transfer rule) and adds them to the run's sync list.
  - Nothing is written beside the source files.
- **EC2 side.**
  - `cloud801.sh prep` copies the dir's manifest and context files into the fold bundle (no relabel) and aborts unless every
    sha256 equals the Mac's.
  - `fold` runs dump801 without `--synthetic`.
  - `frames801.py` and `discover801.py`'s stub path bind with `synthetic = not frames_src.txt`. lib801.bind then reads the
    frames from those controls: per-gene pseudobulk b_ctrl and the MIN_CPM 10 gate, the same arithmetic as for the default
    frames.
- **Default (no `--frames`)** is unchanged: A/B/C's controls relabelled D/E/F (199's def-bundle).
- **Basal-expression floor:** the kit derives none (the panel is given). Only draw801/draw801a do, from A/B/C's context files,
  and they are not part of the kit.

### KA: `--frames data/challenge` (TAG LABCf, `--ref LABC`; run `20261017T102000Z-lab801-kit-labcf`, DONE 0, 16 min, $0.50)
| held | fresh operator == LABC's | dump file sha256 | ka804 arms vs ka-labc |
|---|---|---|---|
| k562_gwps | yes (`48cb88b6…`) | **`eba095c1…` = dump708-LABC's, byte-identical file** | 11 arms, ≤ 2.3e-13 |
| orion_hct116 | yes (`a0e27947…`) | **`e97f1fbd…` = LABC's** | 10 arms, ≤ 2.3e-13 |
| orion_hek293t | yes (`5d982462…`) | **`34ecd80e…` = LABC's** | 10 arms, ≤ 2.8e-13 |

- The ka804 arrays differ by the ≤ 2.8e-13 pricing-thread difference already seen in §17 (2 threads vs LABC's 3).
- Prep recorded "frames = data/challenge's own controls (A,B,C), sha256 == the Mac's".
- The registry has `LABCf/<held>` (55 folds); the arrays are in `work/ka-labcf/`.
- `data/challenge-def` was not touched (it does not exist).

### What the frames change (for the 10-22 read of the two fold variants)
**In the dump (`dump708-<TAG>-<held>.npz`):**
- **Changed:**
  - `b_A|B|C`: the frames' per-gene control pseudobulk, log1p(TS · counts / total);
  - `gate_A|B|C`: genes with control CPM > 10 and b > 0 in that frame;
  - `beta_A|B|C`: 414's board-emulated β pattern (A/B/C's artifact, zeroed on the panel), RESCALED by the fold's realised
    PP|K3 energy, which is rendered on the frames. So **on D/E/F frames the λ 0.5 board emulation still uses A/B/C's β
    pattern**, only re-scaled.
- **Unchanged:** every composition and key array, because they live in effect space and are built before any render: `fin`,
  `fPP`, `fin_<design>`, `kv`, `V`, `dnorm`, `pnorm`, `cosPd`, `cos_dd`, `Fnorm`, `carried`, `avail`. The held-line truth
  (`truth`, `scored`) is unchanged too. CD4 banks and operators do not see the frames.

**In engine801 / ka804 arrays:** everything rendered changes:
- `x`, `a`, `pds`, `pdsl` and `K`, because the carriers' profiles P_DE, the row render (minnorm on the frame's b and gate),
  the per-frame sharpening gate and ballast realisation all read the frames;
- the change energies `dE_*` and their control-CPM bins.

Unchanged: `d` (truth energy) and the arm definitions (w, q, fin). The percentile cap hi is also unchanged, since it is a
function of `kv`.

**So the two variants differ only in how a fixed predicted effect is written into each context's counts and read back.**
The key, the composition and the truth are shared. Δ(variant) on an arm therefore measures frame dependence, not panel or
model dependence.

Cost of this follow-up: **$0.50** (budget $0.6). No local job ran (shell, jq, aws, shasum, and a sub-second registry json
merge by `kit801-collect.sh` with window800 CLEAR). No instance of mine is running. File sha256 prefixes:
`kit/kit801.sh` `1ad7545c`, `cloud801.sh` `d57f8ef9`, `frames801.py` `1e5b9bf7`, `discover801.py` `607b3553`.
