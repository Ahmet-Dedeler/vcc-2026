# 805 — a fresh Flex DE surface from GWCD4i primary CD4 T cells: RESULTS (written as it goes)

Worker for the lab lead [c61d4c]; mandate `MANDATE.md`. Stage 0 is preregistered in `PREREGISTRATION-stage0.md`
(written before any cell-level expression value was read).

## Status
- 2026-10-16 03:23Z: **stage 0 FAILS in both conditions on the blob bar (Rest 0.120, Stim48hr 0.123 > 0.10) while
  power passes in both (n_real median 34.5 / 36; 94 % / 88 % ≥ 10).** By the prereg, 805 stops at stage 0. The lab lead was
  messaged at each result and confirmed "no post-hoc relaxation". The blob-direction diagnostic (descriptive) is done: the
  CD4 blob is not Orion's long-gene axis but is partly a library-depth axis (16–21 % of the mean effect); with that axis removed
  the share is still 0.107–0.109. Total 805 cloud cost $1.90.

## Stage 0 — known answers (run on the instance before any CD4 cell is read)
| row | required | measured (runs s0c and the clean s0e, identical; `work/s0_clean/early/s0_ka/ka805.json`) | verdict |
|---|---|---|---|
| mock test (`test805_mock.py`) | builder counts, labels, axis map, integer counts | MOCK OK | PASS |
| h1-test sha256 on the instance | 06c6baea… | equal | PASS |
| KA-1 (Amendment 2): ref614 on h1-test vs 614/323 | same rows, n_real identical, calls identical, lfc ≤ 1e-9, p rel ≤ 1e-3 (≤ 1e-9 for ≥ 99.9 %) | rows same; **n_real identical for 85/85 (median 89, 58 ≥ 10)**; calls identical (916,215 tests); lfc 2.2e-16; p rel max 1.18e-4, 99.94 % ≤ 1e-9 | PASS |
| KA-2: blob615 split-half / plain share on h1-test | 0.0487 / 0.0376 ± 0.002 | 0.048693 / 0.037613 | PASS |

KA history, in full (the prereg's amendments explain each): launch 1 (`20261016T021043Z-s0-805`) stopped at KA-1 because
the x86 DE parquet sha (796a7e29…) differs from the arm64 Mac's 1bb8c76d… (n_real identical); launch 2
(`20261016T021703Z-s0b-805`) stopped at Amendment 1's 1e-9 p-value bound (521 of 916,215 p-values differ by up to 1.2e-4
relative, tie structure; zero calls move). Both stops happened before any CD4 file was opened. The x86 frame is
deterministic: launches 2 and 3 produced the same sha 796a7e29….

## Stage 0 — what the files are (from obs/var only; launch 4's `select.log`, `work/s0_diagD/early/s0_Rest/`)
- `D{1..4}_Rest.assigned_guide.h5ad`: 3,074,496 / 2,940,194 / 2,778,524 / 2,693,903 cells; X nnz 11.85e9 / 11.67e9 /
  11.47e9 / 9.84e9 (≈ 3,850 nonzeros per cell). Download of all four (504 GiB) on i4i.4xlarge NVMe: 6 min.
- `guide_type`: targeting / non-targeting / '' (unassigned: 794k / 469k / 288k / 401k cells). Usable non-targeting
  (single guide, not low quality): 76,543 / 83,987 / 76,823 / 76,636 → π = 0.2438 / 0.2675 / 0.2447 / 0.2441; control
  quotas 4,485 / 4,922 / 4,502 / 4,491, per-target quotas 97 / 107 / 98 / 98.
- Gene axis: each file has 18,130 genes (no duplicate symbols), the same set in a different ORDER in D3; **17,772 of the
  official 18,533 genes are covered (761 missing → zero columns)**.
- Targeting guide ids per donor differ (22,305 / 21,500 / 20,693 / 23,733); every one is found in the sgRNA and kd tables.
  Many targets are carried by only some donors: under the original shortfall rule ~170 of ~260 selected targets could
  not fill their donor quotas (up to 244 of 400 cells moved). Hence Amendment 4 (every donor fills its quota), so a
  target's donor mix always equals the controls'.

## Stage 0 — the power gate

### Rest (clean run `20261016T024934Z-s0e-805`, Amendment 4; `work/s0_clean/early/s0_Rest/{select,gate805}.json`)
**Verdict: FAIL by the preregistered rule. Power PASSES, the blob test FAILS (0.120 > 0.10).**

| | P (decides; 100 random) | S (A/B/C targets that qualify; reported) | h1-test (reference) |
|---|---|---|---|
| targets | 100 (P_sha 7bda4936…; 256 pool targets skipped for donor quota) | 40 (of 160 with ≥ 400 cells in any donor mix) | 85 |
| n_real median | **34.5** | 63 | 89 |
| share n_real ≥ 10 / ≥ 20 | **0.94** / 0.79 | 0.975 / 0.825 | 0.68 / — |
| n_real q10 / q25 / q75 / q90 / max | 11.9 / 21.75 / 87.25 / 200.7 / 2,426 | 15.9 / 30.25 / 183.5 / 535.6 / 3,070 | |
| DE calls total, DOWN share | 11,513, 0.393 | 10,141, 0.496 | 65,130, 0.526 |
| Spearman vs published DESeq2 n_downstream | 0.598 (published median 2) | 0.764 (3.5) | |
| own-gene pseudobulk remaining q10/q50/q90 | 0.061 / 0.114 / 0.376 | 0.055 / 0.136 / 0.324 | |
| own-gene lfc q50; significant | −3.30; 98 / 99 in universe | −3.05; 39 / 39 | |
| DE universe (> 5 CPM in controls) | 9,244 genes | 9,244 | 10,779 |
| median control library (official axis) | **7,578** (no thinning: f = 1) | same controls | 20,190 |
| **blob split-half share** (plain; mean pair cos; top-PC) | **0.1201** (0.0464; 0.0410; 0.092) | — | 0.0487 (0.0376; 0.0257) |

Truths (on S3 under `output/lab/20261016T024934Z-s0e-805/big/Rest/`): `truth_P.h5ad` 58,400 × 18,533
sha256 8fc61ddd910b3302f6506f60abe371d3e1972332e5bb2b4cbdad3155d7d26494; `truth_S.h5ad` 34,400 × 18,533
sha256 11786f40f6aafed3954b6aee19190edb8a9619f08688ae0445db499ade198abe.

Reading (no rule attached): resting CD4 T cells carry real DE power at 400 cells even at 7.5k UMI: Wilcoxon on single
cells calls ~17× the published pseudobulk DESeq2 count for the same targets, ranked consistently (Spearman 0.60–0.76).
The blob share is 2.5× h1-test's, at 0.38× its depth. The split-half form divides by Σ⟨d_A, d_B⟩, which shrinks when
per-target effects are weak, so part of the 0.12 may be low signal rather than a shared technical axis. The rule says FAIL,
and the record keeps it.

### Diagnostic (launch 4, `20261016T023623Z-s0d-805`, the pre-Amendment-4 shortfall rule; never the verdict)
P read n_real median 97, but split by donor-quota shortfall: **no-shortfall targets 35 (n 26), shortfall targets 121
(n 74)**; S 58 vs 179; Spearman(shortfall size, n_real) +0.32 (P), +0.58 (S). Pooling cells from donors in proportions
that differ from the controls' manufactures DE calls (donor genes). Every CD4 truth must hold each target's donor mix
equal to the controls'. Its blob share was 0.0946. Files: `work/s0_diagD/`.

### Stim48hr (the prereg's one fallback; same instance, same rule; `work/s0_clean/early/s0_Stim48hr/`)
**Verdict: FAIL. Power PASSES, the blob test FAILS (0.123 > 0.10). Stage 0 stops here (no third condition, no re-draw).**

| | P (decides; 100 random) | S (A/B/C targets that qualify) |
|---|---|---|
| targets | 100 (P_sha ba0b018b…; 319 pool targets skipped for donor quota; pool 4,769) | 41 (of 160) |
| n_real median | **36** | 51 |
| share n_real ≥ 10 / ≥ 20 | **0.88** / 0.71 | 0.90 / 0.83 |
| n_real q10 / q25 / q75 / q90 / max | 7.9 / 17 / 101.75 / 260 / 1,600 | 11 / 25 / 117 / 347 / 1,831 |
| DE calls total, DOWN share | 13,141, 0.529 | 6,950, 0.475 |
| Spearman vs published n_downstream | 0.576 (published median 1) | 0.619 (3) |
| own-gene remaining q10/q50/q90; lfc q50; significant | 0.057 / 0.156 / 0.360; −2.81; 100 / 100 | 0.071 / 0.159 / 0.438; −2.82; 40 / 41 |
| DE universe; median control library | 9,472 genes; **10,910** (f = 1) | same |
| **blob split-half share** (plain; mean pair cos; top-PC) | **0.1227** (0.0478; 0.0379; 0.065) | — |

Stim48hr facts: cells per donor 2,805,354 / 2,863,571 / 2,607,532 / 2,815,784; usable controls 73,534 / 87,747 /
77,334 / 83,474 (quotas 4,201 / 5,013 / 4,418 / 4,768; per target 91 / 109 / 96 / 104); the same 17,772-gene axis.
Truths on S3 `output/lab/20261016T024934Z-s0e-805/big/Stim48hr/`: `truth_P.h5ad` sha256
cb6cba78638afa6e6d35c52d7d3e31fdd5a398266939e3971f6ceed2102bb8c1 (58,400 × 18,533); `truth_S.h5ad`
7d8101868f2fbf2ca160d011b73bd5101d0afe6ec341835356fc71adaa6edfd8 (34,800 × 18,533).

### Stage 0 conclusion
GWCD4i has the DE POWER the mandate asked for in both conditions: at 400 cells, Wilcoxon against 18,400 controls gives a
median of 34.5 (Rest) / 36 (Stim48hr) calls per random knocked-down target, 88–94 % with ≥ 10, against H1's 89 / 68 %,
K562's 5 and clean Orion's 1. It does so at 0.38× / 0.55× Arc's depth, and against a published pseudobulk median of 1–2.
It FAILS the preregistered technical-blob bar in both (split-half share 0.120 / 0.123 vs ≤ 0.10; h1-test 0.049), so
under the mandate 805 does not build the surface. The blob-direction diagnostic below asks whether the shared axis is
technical or biological. It is descriptive and does not reopen the verdict.

### Blob-direction diagnostic (descriptive; the lab lead's request; run `20261016T033002Z-blob805b`, `work/blobdir/*.json`)
`blobdir805.py`, in 615's space (log1p(5e4·CPM) pseudobulk, target genes excluded). KA rows on the same instance:
h1-test split-half 0.048693 (= 615 exactly); thinned HEK self-cosine 1.0000 and its share 0.8460 (= 615's 0.8460).

| surface | split-half share | … depth axis removed | cos(μ, HEK long-gene blob) | \|cos(PC1, HEK blob)\| | cos(μ, own depth axis) | \|cos(PC1, depth)\| | μ energy on depth | Spearman(lib ratio, depth proj.) | target/control library q10/q50/q90 |
|---|---|---|---|---|---|---|---|---|---|
| h1-test (Arc H1) | 0.0487 | 0.0480 | +0.037 | 0.089 | −0.120 | 0.058 | 1.4 % | 0.46 | 0.937 / 0.990 / 1.030 |
| Orion HEK thinned (known technical blob) | 0.8460 | 0.8436 | +1.000 | 0.898 | −0.199 | 0.003 | 3.9 % | 0.55 | 0.921 / 1.002 / 1.082 |
| **CD4 Rest P** | **0.1201** | **0.1086** | +0.050 | 0.093 | **−0.399** | **0.624** | **15.9 %** | 0.57 | 0.919 / 0.982 / 1.090 |
| **CD4 Stim48hr P** | **0.1227** | **0.1075** | +0.015 | 0.039 | **−0.456** | **0.567** | **20.8 %** | 0.59 | 0.911 / 0.971 / 1.041 |

(The depth axis = b(top-quartile-library controls) − b(bottom-quartile controls) of the same truth.)

Reading:
- The CD4 blob is NOT Orion's technical blob (cosine with the long-gene axis 0.02–0.05, as low as h1-test's).
- It is partly a library-depth axis: the shared mean effect points down the controls' own depth direction (perturbed cells
  median 2–3 % shallower than controls; h1-test 1 %). At 7.5–11k UMI that direction carries 16–21 % of the mean
  effect's energy (h1-test 1.4 %), and the leading between-target direction follows it (|cos| 0.57–0.62).
- Projecting the depth axis out lowers the share to 0.107–0.109, still just above 0.10. The rest is a weak shared
  "knocked-down vs NTC" direction (mean pair cos 0.038–0.041 vs h1-test's 0.026). That may be biology (an RNA-content or
  stress response common to CRISPRi knockdowns in primary T cells) or residual technical structure. This instrument cannot
  separate the two.
- Mechanism of the depth part: CRISPRi cells carry slightly less RNA, and at low depth log-CPM of low-count genes is
  steep in library size. Arc's 20k-UMI lines are flatter there. A per-target library-matched control pool (or depth
  thinning of controls to each target's library distribution) would be the obvious repair. It is NOT tried here: the
  prereg forbids a post-hoc fix, and the scorer itself uses one shared control pool.

## Files and code
| file | role |
|---|---|
| `PREREGISTRATION-stage0.md` | the rule, the KA rows, the gate, written first |
| `probe805.py` | HDF5 layout probe by S3 range reads (structure, obs, var, indptr only) |
| `s0build805.py` | select / extract (×4 donors, parallel) / merge (thinning, the truths P and S as h5ad) |
| `s0stats805.py` | `ka` (KA-1 DE sha on h1-test, KA-2 blob statistic) and `gate` (n_real, knockdown, blob, verdict) |
| `s0job805.sh` | the EC2 job for one condition (mock test → KA → download → select → extract → merge → ref614 ×2 → gate) |
| `test805_mock.py` | mock test of the builder on four synthetic donor files (runs first on the instance) |
| `blobdir805.py`, `blobjob805.sh` | the blob-direction diagnostic (HEK long-gene axis, own depth axis) with its KA rows |
| `work/s0_clean/` | the clean run's small outputs (select, gate, reference.json, logs, SHA256SUMS.big.txt) |
| `work/s0_diagD/`, `work/s0_ka_fail{1,2}/` | the diagnostic run (pre-Amendment-4) and the two KA stops |

Big files (S3 only, sha256 in `work/s0_clean/collected/s0_<cond>/SHA256SUMS.big.txt`):
`s3://<your-s3-bucket>/output/lab/20261016T024934Z-s0e-805/big/{Rest,Stim48hr}/` (truth_P/S.h5ad,
ref_P/S/{reference.npz, reference_de.parquet, reference.json}, colmap_D*.npy, select.json); copies of truth_P and
select.json at `input/lab/repo/experiments/805-cd4-flex-surface/work/truths/<cond>/`. Nothing > 2 GB on the Mac; no
download to the Mac beyond small json/logs and one 17 MB parquet (`work/s0_ka_fail2_reference_de_x86.parquet`).

Inputs staged to the lab's S3 mirror (`s3://<your-s3-bucket>/input/lab/repo/…`):
`data/harness/h1-test.h5ad` (S3→S3 copy from `input/harness/h1-test.h5ad`, same size 676,544,876 B; sha256 asserted on
the instance = 06c6baea8c62ef3d38cacdbf58e104908eb93f0ef3b25aca56b86d28240f9692, the Mac's copy),
`experiments/614-k562-real-emitter/work/targets_h1.txt` (f49f8343…d6bd1f), the three GWCD4i suppl tables
(`data/external/cd4_cell2026/suppl_tables_2026-08-24/{DE_stats,guide_kd_efficiency,sgrna_library_metadata}.suppl_table.csv`)
and `data/donors/k562_essential/perts.txt`.

## Cloud cost (805 only)
| run | type | cap | actual |
|---|---|---|---|
| 20261016T014750Z-probe805 | r7i.xlarge | 30 min / $0.18 | $0.136 (S3 range reads too slow; terminated at 29 min, no output) |
| 20261016T021043Z-s0-805 | i4i.4xlarge | 240 min / $5.75 | $0.098 (KA-1 sha stop) |
| 20261016T021703Z-s0b-805 | i4i.4xlarge | 240 min / $5.75 | $0.104 (KA-1 1e-9 stop) |
| 20261016T022334Z-s0c-805 | i4i.4xlarge | 240 min / $5.75 | $0.242 (var-order stop) |
| 20261016T023623Z-s0d-805 | i4i.4xlarge | 240 min / $5.75 | $0.449 (diagnostic, complete) |
| 20261016T024934Z-s0e-805 | i4i.4xlarge | 240 min / $5.75 | $0.845 (the gate: Rest + Stim48hr, DONE 0) |
| 20261016T032411Z-blob805 | r7i.xlarge | 45 min / $0.25 | $0.010 (KA rows passed; instance role cannot read output/ → no CD4 rows) |
| 20261016T033002Z-blob805b | r7i.xlarge | 45 min / $0.25 | $0.015 (blob-direction diagnostic, DONE 0) |
| **total 805** | | | **$1.90** of the $80 budget (stage 0 ≤ $10) |

## What stage 1–3 would have needed, and what is left
- Stages 1–3 (the 300-target surface, the real-emitter harness, the trust gate) were NOT built: stage 0 failed its
  preregistered blob bar in both conditions. There is no `run805e.py` driver and no trust-gate table.
- What a later attempt should take from here (each needs its own preregistration):
  1. **Donor mix:** keep Amendment 4 (every donor fills its quota). Only ~28 % (Rest) / ~24 % (Stim48hr) of the eligible
     pool qualifies, and a 300-target panel with strength matching may run short.
  2. **The blob:** test a depth-matched control pool (controls thinned or drawn to match the targets' library
     distribution) as a preregistered variant, judged by the same 0.10 bar with h1-test's KA row. The depth axis explains
     part of the 0.12, not all of it (0.107–0.109 after removal).
  3. **Depth:** the surface is shallower than Arc's (Rest 0.38×, Stim48hr 0.55×), so depth-dependent DE levers (MDE scale
     ∝ 1/√count) would read on a different scale than on A/B/C. Stim48hr is the closer of the two.
  4. Power is not the obstacle: 88–94 % of random knocked-down targets have ≥ 10 calls, and the A/B/C-panel subset has
     a median of 51–63.
- One honest note: the S set (A/B/C targets) was read as a secondary. Its targets are A/B/C's panel, so no lab surface
  may use them; reading it chose nothing.
