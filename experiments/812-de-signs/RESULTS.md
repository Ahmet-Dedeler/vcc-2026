# 812: Lane B3, the DE channel's sign accuracy on called genes. RESULT: a well-supported negative, plus the exchange rates

Worker 812, 2026-10-16 07:06–08:15Z. Mandate: `MANDATE.md`; preamble `experiments/800-checker/PHASEB.md`. Analysis only:
- no emitter run, no build, no `vcc`, no cell-eval2 CLI, no git;
- nothing was written outside this directory, and nothing under `data/` was written.

Cloud: 7 launches, all r7i.xlarge on-demand, **$0.42** actual (b7 = skill812, $0.03) (runs b1–b6 in `cloud/wait-*.log`; `labcloud.sh status`
was clean at the end). Local work: three jobs through `gate602.sh` (12.8 s / 0.53 GB peak, 4.4 s / 0.29 GB, ~13 s / ~0.5 GB; a fourth,
`ntcvar812`, was withdrawn from the queue and run on EC2) and one compact npz extraction (0.1 GB). No local job ran in a reserved window.

## 0. Answer
1. **The a0 map.** a0 is the sign accuracy on the ~88 % of our calls that the truth does not call. It sits at **0.5153
   (h1-test) / 0.5174 (h1-trainval)** and is flat across seeds.
   - **It is 0.50 for every source on weak targets** (n_real ≤ 10 / 22): 0.502 / 0.504, against 0.530 / 0.533 on strong ones.
   - **It is 0.50 on the 42 % of M0 whose truth p > 0.5.** Those signs are the truth's sampling noise, so no predictor can
     move them.
   - What is movable lives in the strong two thirds of targets and in genes near significance. On truth-p < 0.05,
     emitted 0.545–0.557 vs CLEAN4 0.551–0.581.
2. **No emission-side lever in this territory clears noise.** Four new mechanisms were screened analytically on both
   surfaces. All are NO-GO (`PROPOSAL.md`):
   - NTC-guide control bias: real, +0.0002 overall, with a footprint of 0–5 genes per context;
   - a sign-confidence block: net ≈ 0, because fid gains and jac pays;
   - cis neighbours at the front of the order: ≤ +0.0004;
   - head size: 400 is the optimum.
3. **Where the field's DE edge comes from, in units we can act on (§4).** The field median's +0.012 overall at equal bulk
   is reproduced by about **+0.04–0.05 of sign accuracy on the truth's DE genes** (called a1 and uncalled u1) **plus
   ~+0.005 on a0**. That is per-target model skill on responding genes, the size of 635's clean-oracle edge of a line's
   own responding programme. It is lane 811's objective, not an emission trick. The per-class exchange rates below price
   any future sign source before an emitter run.

## 1. Inputs (path, sha256 prefix) and known answers
- Truth and stats:
  - h1-test: `experiments/323-current-family-validation/work/reference.npz` (304aa511…), `reference_de.parquet` (1bb8c76d…).
    Production-emitter zs stats: `435-zs-on-coarse/work/h1/stats-base_s{0,1,2}_bps_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.npz`
    (9cee6074…, 00ea20db…, 7accfe37…).
  - h1-trainval: `618-h1-trainval-surface/work/ref_tv/reference.npz` (5a567f09…), `reference_de.parquet` (ccddb720…).
    Stats: `g1/stats-base_s{0,1,2}_bps_zr-5_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.npz` (d622cf79…, c11303c0…, 4e56346a…).
  - These are the hook-free zs arms, i.e. the incumbent's calls and called-gene signs. The hook moves only uncalled genes,
    and UPRN / DN are rank-neutral.
- Keys: `323/work/effects.npz` (1ab617c1…), `618/work/tv/effects.npz` (8119e81b…, = run618e's recorded key_sha256).
- CLEAN4: `391/work/clean4_prior_challenge.npz` (7e532f84…).
- Cells and annotation: `data/harness/h1-test.h5ad` (06c6baea…); `data/challenge/context_{A,B,C}.h5ad` (f22e7196…,
  41557555…, e090e7c7…); GENCODE v44 TSS `425/work/tss_v44.tsv` (8b4948a1…).
- H1 NTC guide ids: `work/stage/h1ctrl_guides.tsv` (51102f53…). It maps all 18,400 h1-test control barcodes into the VCC2025
  `adata_{Test,Training,Validation}.h5ad` obs: 31 guides, 243–1,331 cells each, 48 batches; only obs was read.
- EC2 copies of the small reference arrays: `work/stage/ref_{h1test,h1tv}.npz` (64f393a3…, 4278bdf1…). They hold the same
  8 arrays as the reference.npz files and are gitignored.
- **Known answers, all PASS:**
  - fid recomputed from the stats equals the arm's recorded raw fid to 1e-12 on all 6 (surface, seed) rows (h1test s0
    0.5089291277653674, h1tv s0 0.5052857989716864);
  - reach and nmae recomputed equal the recorded values to 1e-12 (`cis812`, `rates812`);
  - the scorer's lfc recomputed from h1-test's cells against the control mean equals the reference table to 1.04e-5
    (3 targets; 491 got 1.1e-5).

## 2. Step 1: a0 by gene class (`a0class812.py` → `work/a0class812.json`, 9b10977…; seeds 0–2 agree to ±0.0001)
| | h1-test s0 | h1-trainval s0 |
|---|---|---|
| a0 (called, truth non-sig) / n | 0.5153 / 193,554 | 0.5174 / 282,584 |
| a1 (called, truth sig) / n | 0.6005 / 24,954 | 0.6216 / 41,184 |
| share of calls not truly DE | 0.886 | 0.873 |
| a0 head / block | 0.520 / 0.514 | 0.521 / 0.517 |
| a0: emitted vs CLEAN4 on CLEAN4-covered genes | 0.5163 vs 0.5171 | 0.5212 vs 0.5225 |
| a0 on CLEAN4-uncovered genes (18 % of M0) | 0.5106 | 0.5004 |
| emitted DOWN share on M0 vs the truth's | 0.644 vs 0.511 | 0.643 vs 0.508 |
| block genes delivered with the CLEAN4 sign | 0.891 | 0.891 |
| a0 by CLEAN4 confidence \|P−0.5\| < 0.05 / 0.05–0.1 / 0.1–0.2 / ≥ 0.2 (CLEAN4's own sign) | 0.506 / 0.527 / 0.531 / 0.585 | 0.510 / 0.532 / 0.544 / 0.564 |
| a0 by target strength weak / mid / strong | 0.503 / 0.519 / 0.530 | 0.504 / 0.520 / 0.533 |
| a0 by truth p (share): < 0.01 (2.5 %), 0.01–0.05 (8.4 %), 0.05–0.2 (19 %), 0.2–0.5 (28 %), > 0.5 (42 %) | 0.557 / 0.545 / 0.530 / 0.517 / 0.499 | 0.551 / 0.545 / 0.529 / 0.514 / 0.507 |
| a0 by control CPM | flat (0.511–0.526 in every bin with n > 1,000) | flat (0.507–0.524) |

**What an oracle per class would buy (fid at f 1, raw → overall-eq).** Each row sets the truth's sign on that class's
called non-sig genes:

| class | h1-test | h1-trainval |
|---|---|---|
| all M0 | +0.424 → +0.235 | +0.415 → +0.230 |
| block M0 | +0.357 → +0.198 | +0.348 → +0.193 |
| head M0 | +0.066 → +0.037 | +0.065 → +0.036 |
| M0 with truth p < 0.05 | +0.043 → +0.024 | +0.045 → +0.025 |
| M0 with truth p < 0.2 | +0.120 → +0.066 | +0.121 → +0.067 |
| the line's own leave-one-target-out mean response (a0 0.5405 / 0.5406) | +0.023 → +0.012 | +0.020 → +0.011 |

The last row is the "0.552 oracle" of 700 / 721. 635 showed that it is mostly the shared control-sample noise. That noise
cannot be predicted from the controls themselves; §3.1 tests the one part of it that could be.

Real-significant genes by CPM (`sigcpm812.py` → `work/sigcpm812.json`): the up share is 0.39–0.55 in every bin, so
there is no Wilcoxon-sparsity UP lean worth a prior. Uncalled real-significant genes carry 74–76 % of the nmae
denominator. Their CLEAN4 sign accuracy is 0.55–0.59 against 0.47–0.56 emitted pre-hook; that is the incumbent hook's
territory.

## 3. The four screens (details and try counts in `PROPOSAL.md`)
### 3.1 H-NTCB: NTC-guide-specific effects bias the scorer's control mean (`ntcbias812.py`, `ntcvar812.py`, `ntcflag812.py`)
- **Overdispersion between guides is at the null in A/B/C and H1.**
  - The ratio R (n-weighted between-guide mean square over cell variance) is ~1 under pure sampling, with q10/q90
    ≈ 0.73/1.27 for 46 guides.
  - B: median 1.005, q90 1.30. C: 1.001, 1.30. A: 1.056, 1.42. H1: 0.989, 1.37.
  - H1's 48 BATCHES do deviate (median 1.215, q99 7.4), but batch composition is shared by targets and controls.
    Batch-robust estimators predict nothing: accuracy 0.47–0.51 in every \|z\| bin.
- **The mechanism is real where it acts.**
  - On H1 the guide-Huber estimator flags 5–6 genes.
  - On their called non-sig cells, −sign(e) is right 0.69 (h1test, n 293) / 0.66 (h1tv, n 426). CLEAN4 gets 0.43 / 0.44
    on the same cells.
  - Overriding the emitted sign there gives fid +0.00030 / +0.00029 raw = **+0.00017 / +0.00016 overall**.
- **Footprint on A/B/C: 5 / 1 / 1 genes at \|z\| ≥ 2 and 0 / 0 / 1 at \|z\| ≥ 3.** It is worth ~0 on the board. NO-GO.
- **By-product (science).** NTC guide `non-targeting-30` represses **CST3** in all three contexts (to 0.07 / 0.36 / 0.44 of
  the pool mean; 30–59 SE). Other single-guide repressions in A: NUP37 (NTC-30), SMARCD1 (NTC-28), SLC12A7 (NTC-39),
  RAD50 (NTC-68). These are reproducible guide-specific (seed-like) off-target effects inside the NTC set, consistent
  with F16 / 643. The scorer's control mean of CST3 is about 1–2 % low on every context, so every target's real CST3
  lfc leans UP (`work/ntcflag812.log`, 3f6831f9…).

### 3.2 H-CBLOCK: compose the block by CLEAN4 confidence (`pool812.py` → `work/pool812.json`, d73ea935…)
- CLEAN4 accuracy rises with confidence equally at every expression level. On all non-sig genes it is 0.507 at
  \|P−0.5\| < 0.05 and 0.52–0.54 at 0.1–0.2.
- Re-choosing 25 / 50 / 100 % of the block by confidence (CPM floors 10 / 20 / 50, same K, emulated with CLEAN4 signs)
  raises a0 by +0.003 … +0.010 and fid by +0.0034 … +0.0090 raw.
- But jac falls −0.0037 … −0.029: the CPM block's high-expression genes are more often truly DE.
- **Net −0.0091 … +0.0004.** NO-GO. It is the mirror of 528 (call propensity), where jac gains and fid pays.

### 3.3 H-CISFRONT: cis neighbours at the front of the order (`cis812.py` → `work/cis812.json`, 17ebefd4…)
- TSS ≤ 1 kb pairs: 14 (h1test) / 18 (h1tv); 12 / 15 of them are real-significant.
- Before the move, our signs were right on 12 / 12 and 14 / 15; the median order rank was 142–152 (h1test) and 5–8 (h1tv).
- Moving them to the front with lfc ≤ −1: reach +0.00014 / +0.0021 raw; **+0.00003 / +0.0004 overall** (3 seeds
  identical). At ≤ 5 kb: +0.00003 / +0.00034.
- NO-GO. It confirms 425 for the order too.

### 3.4 H-HEAD: head size at fixed K 2,600 (`head812.py` → `work/head812.json`, c0e24c10…)
| H | 0 | 200 | 400 (today) | 800 | 1,200 | 1,600 | 2,600 |
|---|---|---|---|---|---|---|---|
| precision (TP / calls), h1test | 0.113 | 0.114 | 0.113 | 0.111 | 0.108 | 0.102 | 0.078 |
| Δoverall-eq (fid + jac), h1test | −0.0017 | −0.0004 | 0 | −0.0004 | −0.0018 | −0.0040 | −0.0148 |
| Δoverall-eq, h1tv | −0.0015 | −0.0005 | 0 | −0.0009 | −0.0020 | −0.0043 | −0.0166 |

- 400 is the optimum.
- **A key-ranked head has LOWER precision than the CPM block.** So "the pack calls genes that are truly DE" is not
  reachable by ranking our key further down.

## 4. Exchange rates: what +x of sign accuracy buys, per class (`rates812.py` → `work/rates812.json`, b1bff0f4…)
Method: H1, f 1, analytic. Calls, order and magnitudes are fixed. A random share of the class's wrong-signed cells is
flipped to the truth's sign (3 draws). Overall-eq uses the board's INV (fid 3.32482, reach 1.12927, nmae 1.64787) / 6.

| class (accuracy now: h1test / h1tv) | +0.01 → overall-eq h1test / h1tv | +0.05 → h1test / h1tv | per corrected cell |
|---|---|---|---|
| a0 called non-sig (0.515 / 0.517) | +0.0048 / +0.0048 (fid only) | +0.0242 / +0.0238 | 2.5e-6 |
| a1 called sig (0.600 / 0.622) | +0.0012 / +0.0021 (fid, reach, nmae) | +0.0073 / +0.0088 | ~5e-6 |
| a1head (0.722 / 0.738) | +0.0003 / +0.0005 | +0.0023 / +0.0018 | ~9e-6 |
| u1 uncalled sig (0.537 / 0.544) | +0.0011 / +0.0012 (nmae, reach) | +0.0061 / +0.0061 | 2.7e-6 |

- **The field's DE gap at equal bulk, read through these rates.** 721 measured it as fid +0.0051, reach +0.0029 and
  nmae +0.0077 overall.
- a0 alone would need +0.011 for the fid part. It cannot produce the nmae or reach part.
- The nmae part alone (−0.028 raw) needs u1 +0.07–0.08, or a1 +0.12, or both about +0.05.
- **+0.05 on a1 and u1 plus +0.005 on a0 gives +0.0158 … +0.0173 overall.** That sits between the field median
  (+0.012) and its p75 (+0.022).
- So the pack's edge is consistent with gate-wide per-target sign skill on responding genes, about the +0.045 that 635's
  clean oracle of a line's own responding programme holds over CLEAN4 (0.618 vs 0.573). This is not an emission or
  call-composition effect.
- Caveat: these are H1 rates. Board transfer differs by class (`docs/H1-BOARD-CALIBRATION.md`: sign priors 0.07–0.67×
  on fid), so the decomposition gives orders of magnitude, not a fit.

## 5. Tries, wear and what was not done
- **Tries (all analytic, on the worn DE-H1 surfaces):**
  - H-NTCB: 4 estimators × 5 bins, and 3 override thresholds;
  - H-CBLOCK: 9 configurations;
  - H-CISFRONT: 2 distances × 3 seeds;
  - H-HEAD: 7 sizes.
  - None was selected as a winner, so no best-of-N correction applies. Any later reuse of a configuration counts these
    tries.
- **Not done:**
  - No emitter run: no hypothesis reached a size worth one (803's engine cost ~1 h per surface-seed).
  - No CD4 Flex surface (805).
  - No checker `price` call: nothing was priced as a candidate, so there is no wear-ledger entry.
- **Rule notes:**
  - `PROPOSAL.md` was written after the analytic screens, not before. They were step-1 measurements and ask for no GO.
  - `surf812.py` gained the compact-ref path and the reach helpers after `a0class812` / `pool812` / `sigcpm812` ran. The
    loaded values are unchanged: the 8 arrays are identical, and the KA re-passed in `cis812` and `rates812`.
  - `ntcflag812.py` wrote no json on EC2 (no `work/` dir there). Its full output is the job log `work/ntcflag812.log`.

## 6. What is left (for the lead)
- **Close this lane's emission-side search for October.**
  - The DE channel's emission (call set, order, head size, block composition, CPM-coupled and control-derived sign
    sources) is at its optimum for the sign information we hold.
  - Today's closures: §3. Earlier closures: 344, 390, 406, 440, 491, 513, 528, 601, 610, 613, 635, 651, 657, 700.
- **Hand §4 to lane 811 as its DE objective.** Score a new model's per-target signs on held-out lines by class (a0 / a1 /
  head / u1) and price them with §4's per-cell rates before any emitter run. A model that gains ≥ +0.03 on truth-DE genes'
  signs over CLEAN4 leave-one-line-out is worth an emitter run on DE-H1. The block / hook prior and the head signs are
  the delivery routes; they already exist as builder knobs (`VCC_DE_SIGN_PRIOR`, the hook).
- The CST3 / NTC-30 off-target is a small science note for FINDINGS (production's call; the lab does not edit shared docs).

## 7. skill812: the DE sign-skill instrument (follow-up from the lab lead, 2026-10-16 08:15Z; run b7 08:25–08:31Z, $0.03)
**What it is.** `skill812.py` (29b5f134…) scores ANY source's per-target signs on H1's 211 targets (h1-test 85 +
h1-trainval 126). The input is one npz per surface (`targets`, `genes`, `pred`; 0 / NaN = no sign).
- Classes come from the incumbent's realized call set (the zs arm, s0): a0, a1, a1head, u1, each also by control CPM bin.
- Every class is reported against CLEAN4 (paired cells) and against the incumbent's delivered sign: the hook arm's realized
  lfc (h1-test `612/work/h1_c4kgall` s0, h1-trainval `618/work/hk_c4h` s0).
- Statistics: paired target bootstrap (4,000 draws), and a pooled equal-weight mean of the two surfaces.
- It prices each class at rates812's "+0.01" H1 rates, in overall-eq units.
- **Caveat (printed on every table): the rates are H1's.** Board transfer differs by class (sign priors 0.07–0.67× on fid).
- Usage: `uv run python experiments/812-de-signs/skill812.py --name N --pred h1test=A.npz --pred h1tv=B.npz --h1-free yes`.

**Known answers, both PASS** (`work/skill/ka_clean4.json` 6b376edb…, `ka_emitted_zs.json` ac160add…):
- CLEAN4's a0 is **0.5171 / 0.5225**. Source rows: 700's `signsplit700.json` row `callset_head|block|all|b_nonsig`
  clean4 0.5171, and a0class812's `a0_c4_on_c4cov` 0.5171 / 0.5225.
- The zs arm's own realized signs reproduce rates812's class accuracies to 5e-5: a0 0.5153 / 0.5174, a1 0.6005 / 0.6216,
  a1head 0.7222 / 0.7383, u1 0.5366 / 0.5444.
- Also on record: the incumbent's hook moves u1 by +0.007 / +0.009 over the pre-hook emitted sign.

**Step-3 sources (information only; tries 4; inputs in `work/sources/`, gitignored).**
- (a) `bank812.py` (01d10eac…) extracts the H1 targets from 16 non-H1 banks of 12 lines; manifest
  `work/h1banks/manifest.json` (eaf6d4df…).
  - K562: gwps 196, essential 43, VIPerturb 201, Jiang 20 targets. HCT116 203, HEK293T 204, Jurkat 56, RPE1 57, DLD1 64,
    A549 15 + 20, BXPC3 / HAP1 / HT29 / MCF7 20 each, TeloHAEC pilot 41.
  - VCC25 H1 and Pan's H1 are never read. CD4 (GWCD4i) and kolf2's full rows are not included: they are not in donor
    format on disk.
  - `sources812.py` (e324956e…) builds the per-target mean over lines (banks averaged within a line) in three forms:
    raw, bank-RMS-scaled, and per-line centred (each bank's all-perturbation mean row subtracted).
  - The mean covers 99.7 % of cells, with 4.1–4.2 lines per signed cell.
- (b) The production bulk delta's own sign = the sign of the ungated cap15 composition. ps_bulk is cap15 times a positive
  per-target scalar, and sharpening keeps signs.

Pooled overall-eq [95 %], source replacing the comparator on the cells it signs:

| source | a0 vs CLEAN4 | a0 vs incumbent | a1 vs incumbent | a1head vs incumbent | u1 vs CLEAN4 | u1 vs incumbent |
|---|---|---|---|---|---|---|
| pubmean_raw (`ce10e51a…`) | −0.0025 [−0.0043, −0.0007] | −0.0020 [−0.0034, −0.0005] | −0.0011 [−0.0035, +0.0014] | −0.0004 [−0.0009, +0.0001] | **+0.0013 [+0.00007, +0.0025]** | **+0.0022 [+0.0010, +0.0034]** |
| pubmean_rms (`4e91d6c6…`) | −0.0023 [−0.0041, −0.0005] | −0.0016 [−0.0030, −0.0002] | −0.0014 [−0.0037, +0.0009] | −0.0001 [−0.0004, +0.0002] | +0.0012 [−0.0001, +0.0025] | +0.0022 [+0.0009, +0.0034] |
| pubmean_ctr (`c38f9a8b…`) | −0.0036 [−0.0055, −0.0016] | −0.0029 [−0.0045, −0.0013] | −0.0041 [−0.0063, −0.0018] | −0.0005 [−0.0007, −0.0002] | +0.0007 [−0.0008, +0.0019] | +0.0013 [−0.00003, +0.0025] |
| bulk_sign (`2017fb34…`) | −0.0041 [−0.0059, −0.0024] | −0.0040 [−0.0052, −0.0027] | −0.0052 [−0.0066, −0.0037] | 0 (it IS the head's sign) | −0.0003 [−0.0013, +0.0005] | −0.0001 [−0.0005, +0.0004] |

(a1 / a1head "vs CLEAN4" is not a useful comparison: CLEAN4 is never the head's sign, and every key-like source beats it there
by +0.15–0.19 accuracy.)

**Reading.**
- **The bulk sign beats CLEAN4 nowhere it could be used.**
- **The public per-target mean beats CLEAN4 in exactly one paying class, u1 (the hook's genes).**
  - raw: +0.0013 [+0.00007, +0.0025]. The interval clears 0 only before the best-of-3 shift (≈ −0.0005).
  - Against the incumbent's delivered hook sign it gains +0.0022 [+0.0010, +0.0034]. That survives the shift.
  - Most of it sits at control CPM 5–20, where the key is gated to 0 and the hook has no confident CLEAN4 sign.
- It loses on a0, so it is never a block prior.
- Per-line centring removes about half of the u1 gain and deepens the a0 loss. The line programme is part of what helps on
  H1's uncalled genes. That is the warning sign: 516 showed this family's H1 gain does not carry to held-out cancer lines.
- Hence `PROPOSAL-u1hook.md`. The worker's prior is NO-GO for October. It asks the judge for a cheap held-out-line stage 0
  (< $2) before any DE-H1 price (~$5).

**For 811.** skill812 is the DE objective. Hand it a model's H1-free per-target predictions and it returns the class table
and the price.
- Bar for an emitter run: pooled u1 + a1 price vs the incumbent ≥ +0.002 with lower bound > 0, a0 not below CLEAN4 by more
  than its CI, and confirmation on a held-out line.
