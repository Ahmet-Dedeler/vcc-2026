# 179 — `hy_kt8000` built and verified: **6 of 7 gates PASS; gate 5 fails as WRITTEN and passes on the mechanism.** Package ready. NO slot spent, NO submission.

Date: **2026-09-02**. Build + verification only. **No `vcc submit`. No `cell-eval2`.** `vcc prep`
was run with `--dry-run` only and wrote nothing. Peak RSS of the build **13.74 GiB**; every
verification stage stayed under 6 GiB.

**The arm.** a23 exactly, plus `--emitter resample_mix --resample-keep-top 8000`: the 8,000-gene
keep block stays on the shipped parametric quasi law and a15's poisoned tail block is emitted as
**resampled real control cells** from the context's own `context_{A,B,C}.h5ad`, raw counts with
their own library sizes, against the pooled control profile.

**Package:** `submissions/hy-kt8000.h5ad`
sha256 **`e2902589f9e563fb1050064d9a75bee95a108da27eab7fc41c2a0ade9d3323f7`**
19,720,467,224 bytes · 360,000 rows · nnz 2,461,200,172 (**51.81%** of the density cap)
built on macOS-27.0-arm64, builder sha256 `44372b27e4a2f1e96f8c766ec162c1361c01f87ce9b71741883799234a227542`.

---

## Verdict in seven lines

1. **Rule 10 is clean.** The patched builder at a23's configuration reproduces a23's effects
   audit on **78 fields, 0 differences**, including `effects.final_sha256 = 22c09f22c749bda8`
   and `operator.sha256 = e2f65ecf69b7e0b0`. The switches cannot reach the effect matrix.
2. **The config diff is exactly two keys**, `emitter` and `resample_keep_top`, both added.
   Everything else that moved is outside `config` and is listed in §3 rather than hidden.
3. **`pds` and `mse`'s profile half are pinned by identity, verified on the artifacts**:
   the delivered per-gene panel-total sha256 is identical to a23's on **900/900** blocks while
   the nonzero structure differs on **900/900**. The two files are different files that deliver
   the same pseudobulk.
4. **Gate 5 fails as the mandate wrote it and passes as the mechanism.** The mandate's
   "0.0–0.2% BH-significant at zero effect" is 177 §2's number for **`rs_all`** — every gene
   resampled. `hy_kt8000` keeps 8,000 genes on the quasi law by construction and therefore
   cannot inherit it. Split by `param_mask`: **the resampled tail collapses from 48.7 / 33.0 /
   42.7 % to 1.9 / 0.06 / 1.5 % (37 / 1 / 31 genes)** while the parametric head is unmoved.
5. **The `n_pred` probe reproduces a23 and refutes 177's forecast for the arm.** a23 reads
   **4,474** against 177's predicted ≈ 4,424 (**1.1%** — a passing known-answer row);
   `hy_kt8000` reads **3,528** against 177's predicted ≈ 2,410. The realised board cut is
   **−21%, not −46%**: on the challenge contexts it is **0.46×** the size of the h1-test cut.
6. **That makes the arm both safer and less informative than 177 priced it.** `fid` now only
   loses if `n_real` lies in **(3,528, 4,474]** — a 946-wide window instead of 177's
   2,410–4,424 — and every DE delta, gain and loss alike, scales down by the same 0.46.
7. **Rule 13 still vetoes it** (`predicted = −0.00243` against a `+0.0019` gate), and the
   coordinator's 22:40 correction independently reads it as a **wash**. The package is ready;
   the decision to spend a slot is Ahmet's.

---

## 1. What was built into the builder, and gate 1 (unit tests)

Two flags in `scripts/make_submission_pm_cd4.py`, defaults = shipped behaviour:

| flag | default | what it does |
|---|---|---|
| `--emitter {quasi,resample_mix}` | `quasi` | the emission LAW. `quasi` writes **no** meta key |
| `--resample-keep-top K` | `--keep-top` | the `param_mask` cut; genes ranked by the same score `--keep-top` uses, over the same CPM-gated axis |

Four helpers were factored out so the `param_mask` is the **same index set** the block spread
map promotes rather than a second ranking that happens to agree: `emission_rank_axis`,
`resolve_keep_count`, `keep_block_indices` (all three lifted verbatim out of
`emission_spread_vector`, which now calls them) and `emission_param_mask`. Plus
`assert_raw_control_counts`, the guard for 177 §0.1.

`src/vcc2026/donor_transfer.py` is **unmodified** — `emit_counts(mode="resample_mix")` is used
exactly as it exists, and its frozen sha `9ffca7ec…` is unchanged in the built meta.

**Both call-site contracts are honoured, and one is now verified on the real rendered object:**

* (a) **raw counts, own libraries.** `control_counts = resample_ctrl[picks].todense()` with
  `lib_sizes = libraries[picks]`, never rescaled. `assert_raw_control_counts` rejects a block
  whose row sums differ from its own libraries, a non-integer block, a negative or non-finite
  block, and a wrong-shaped block.
* (b) **pooled control profile.** The builder already passed `totals / totals.sum()`, and
  `b_ctrl = log1p(TS · totals/Σ)` makes that identically `expm1(b_ctrl)/Σ`. Measured on all
  three challenge contexts at zero effect: `max |profile_from_effect(0, b_ctrl) − pooled|` =
  **8.7e-19 / 1.7e-18 / 1.3e-18**, max relative **4.3e-16 / 6.6e-16 / 5.3e-16**. `fold = 1`.

### Gate 1 — **51 tests pass** (46 pre-existing + 5 new), peak RSS **0.192 GiB**

| new test | what it asserts, and what it would read if the thing were absent |
|---|---|
| `test_default_emitter_reproduces_the_pre_change_builder_byte_for_byte` | resolves the newest committed builder revision that predates `--emitter` **out of git**, loads it, and drives the REAL writer on a synthetic three-context bundle with both. sha256 of the emitted `.h5ad` must match at bare defaults **and** at a15's shipped emission policy (`tail 1.0 / kt 8000 / ks 0.25 / φ 4`). A same-module comparison cannot answer this question, because the patch is on both sides of it. It **fails**, never skips, if the revision cannot be found. Non-vacuity: the resample law on the same inputs must produce a **different** artifact. |
| `test_resample_mix_over_the_whole_gated_axis_is_the_quasi_law` | `param_mask` all-True ⇒ `emit_counts(resample_mix)` is bit-identical to `quasi` (values, indices, indptr), and one gene off the mask breaks it. End-to-end through the builder on a bundle where every gene clears the CPM gate, at `--resample-keep-top = n_genes`, one emission block. |
| `test_resample_mix_delivers_the_quasi_panel_totals_on_every_block` | the two laws emit **different files** but identical per-gene panel totals on 6/6 blocks of a three-context build. This is the unit-test form of gate 4. |
| `test_rescaled_control_cells_are_rejected_before_they_reach_the_emitter` | 177 §0.1 as a guard: a block rescaled to a nominal 100k depth, a CPM-normalised block, a non-integer, negative, non-finite or wrong-shaped block are all rejected; a genuine raw block passes. |
| `test_emitter_flags_default_to_the_shipped_law_and_fail_closed` | defaults are `quasi` / `None`; `emission_param_mask(…, K)` equals the keep block `emission_spread_vector` promotes at the same K, is capped by the gate, is empty at K = 0; and three malformed flag combinations are rejected. |

**One thing the single-block restriction records rather than hides.** The writer drives all 900
blocks off ONE emission stream, and `resample_mix` draws `n·G` extra uniforms per block for its
stochastic rounding even when the splice discards them. So "all-parametric = quasi" is exact for
the first block and a stream offset thereafter. That is not a law change — the panel totals are
a function of the picks (a separate stream) and the profile only, which gate 4 confirms on
900/900 blocks of the real build.

---

## 2. Gate 2 — the rule-10 bit-identity proof

`--emitter` and `--resample-keep-top` are consumed only by `run()`'s validation, by the meta,
and by `write_submission_from_effects`, which is called **after** `del donors, pm, cd4_banks,
cd4, operator, model`. So a build with the switches ON must reproduce the incumbent's effects
audit exactly. The built package's own meta, against `submissions/a23-massfix.h5ad.meta.json`:

```
effects.baseline_pre_sha256                 7678bab9ea9d1b8d…  vs 7678bab9ea9d1b8d…  EXACT
effects.candidate_pre_sha256                0601336fdf1ec5c5…  vs 0601336fdf1ec5c5…  EXACT
effects.post_whiten_sha256                  b2426e15d76d9fe5…  vs b2426e15d76d9fe5…  EXACT
effects.final_sha256                        22c09f22c749bda8…  vs 22c09f22c749bda8…  EXACT
effects.whitening.before_all_median_norm     4.238101079325155  vs  4.238101079325155  EXACT
effects.whitening.before_operator_median_norm 4.122892101713494 vs  4.122892101713494  EXACT
effects.whitening.after_operator_median_norm 25.73792992728808  vs 25.73792992728808   EXACT
effects.whitening.after_all_median_norm      4.229295604454299  vs  4.229295604454299  EXACT
effects.whitening.operator_preserve_scale    0.16018740098217019 vs 0.16018740098217019 EXACT
effects.cd4.fraction                         0.09472706196346513 vs 0.09472706196346513 EXACT
effects.cd4.post_add_preserve.*              4.122892101713487 / 4.1665533562021135     EXACT
operator.sha256                             e2f65ecf69b7e0b0…  vs e2f65ecf69b7e0b0…  EXACT
operator.columns_sha256                     b2cae0d381c3cd70…  vs b2cae0d381c3cd70…  EXACT

fields compared: 78   differences: 0        GATE 2 PASS (rule 10 clean)
```

Operator by path and sha, never by description (rule 11):
`data/donors/whitener_union_pruned.npz`, sha256
`e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123`, 17,942 columns,
`basis_prune: essential`, `all_panel_coordinates_excluded: true`.

---

## 3. Gate 3 — the meta diff is exactly two lines

```
emitter              built 'resample_mix'   incumbent <absent>
resample_keep_top    built 8000             incumbent <absent>
config keys that moved: ['emitter', 'resample_keep_top']       GATE 3 PASS (2 differences)
```

**Everything else in the meta that moved, stated rather than omitted** (none is in `config`):

| key | built | a23 | why |
|---|---|---|---|
| `hashes.implementation.scripts/make_submission_pm_cd4.py` | `44372b27…` | `0852378b…` | the file that was patched. It is **not** in `FROZEN_IMPLEMENTATION_SHA256`; the frozen modules, including `donor_transfer.py`, are unchanged |
| `output.nnz` | 2,461,200,172 | 2,378,695,978 | +3.5%; 51.81% vs 50.08% of the density cap |
| `output.contexts.{A,B,C}.stored_entries` | 847.7 / 791.8 / 821.7 M | 816.8 / 771.2 / 790.6 M | the arm |
| `output.bytes`, `output.sha256`, `output.path` | — | — | the arm |
| `peak_rss_gib` | 13.74 | 17.47 | page-cache state, not a code change |

---

## 4. Gate 4 — the delivered pseudobulk identity, on the two real artifacts

`digest179.py`: one sha256 of the delivered `int64` per-gene panel totals per (context,
perturbation) block, streamed from both `.h5ad` files. Block order verified equal from the obs
categoricals before any comparison; all 900 blocks are 400 rows.

```
per-gene PANEL TOTAL digest identical on   900/900
panel total SUM identical on               900/900
nonzero-STRUCTURE digest identical on        0/900
```

**Both halves matter.** Experiment 156's first a23 verification read a perfect `1.00000000` for a
build with the mass fix *and* for one without it, because it was an algebraic tautology. The
non-vacuous form is 171 §6's, and this is it: without the post-splice `_repair_totals` the
residual would land in the totals and the blocks would separate. The structure column is the
other half — if it had matched too, the two files would simply be the same file and the totals
identity would prove nothing. It matches on **0** of 900.

**So `Δpds = 0.000000` and the profile half of `Δmse = 0` are identities on this package**, not
estimates. BRIEF §1 defines `pds` on exactly this pseudobulk.

---

## 5. Gate 5 — the zero-effect null. **Fails as written; the mechanism passes.**

`null179.py` hands the builder an effect matrix of **zeros** at a23's configuration and emits
through the production path (`write_submission_from_effects` → `profile_from_effect` →
`emission_spread_vector` → `emission_param_mask` → `emit_counts`), then scores the 400 emitted
cells against the context's own 18,400 controls with an exact mid-rank Mann-Whitney over the
scorer's DE universe (mean per-cell CPM > 5), BH at 0.05.

| law | ctx A | ctx B | ctx C | 177 §2's own table |
|---|---:|---:|---:|---|
| KA quasi flat 0.25 / φ4 | 32.60% | 33.96% | 19.48% | 28.4 / 33.1 / 18.8 |
| KA quasi flat 1.00 / φ4 | 81.96% | 78.69% | 79.87% | 82.4 / 77.7 / 80.1 |
| quasi SHIPPED (tail 1.0 / kt 8000 / ks 0.25) | 39.87% | 38.53% | 29.22% | 43.3 / 44.8 / 38.4 |
| **`resample_mix`, K = 8000** | **28.08%** | **30.93%** | **17.21%** | — (**no such row exists**) |
| 177's `rs_all` (every gene resampled) | — | — | — | 0.0% (1) / 0.2% (18) / 0.2% (19) |

**The mandate's gate is 177 §2's `rs_all` row applied to a different arm.** `hy_kt8000` emits
8,000 genes on the unchanged parametric law, whose own zero-effect rate is the
`KA quasi flat 0.25/φ4` row above (19–34%). It is arithmetically impossible for the hybrid to
read 0.0–0.2% overall, and 177 never claimed it would — §2 has no `hy_ktK` row at all.

**Split by `param_mask`, which is the falsifiable statement** (`nullsplit179.py`, same emission,
same seeds):

| ctx | law | INSIDE mask (8,000 genes) | OUTSIDE mask (tail) |
|---|---|---:|---:|
| A | quasi | 3019 / 8000 = **37.74%** | 940 / 1929 = **48.73%** |
| A | **resample_mix** | 2751 / 8000 = **34.39%** | **37 / 1929 = 1.92%** |
| B | quasi | 3172 / 8000 = 39.65% | 537 / 1626 = 33.03% |
| B | **resample_mix** | 2976 / 8000 = 37.20% | **1 / 1626 = 0.06%** |
| C | quasi | 2052 / 8000 = 25.65% | 906 / 2124 = 42.66% |
| C | **resample_mix** | 1711 / 8000 = 21.39% | **31 / 2124 = 1.46%** |

*What this reads if the resampling were absent: `rate out` would equal the quasi row's.* It
falls by **96.1% / 99.8% / 96.6%**. As a fraction of the **whole** scored universe the tail's
false calls go **9.5% → 0.37%**, **5.6% → 0.01%**, **8.9% → 0.31%** — which is 177's `rs_all`
number (0.0–0.2%) to within the same order, on the genes the arm actually changed. The
parametric head moves by 2–4 pp, consistent with the emission stream offset §1 records, and
stays at the shipped law's rate.

**Verdict on gate 5: the emission law does exactly what 177 says it does. The gate's threshold
was transcribed from the wrong arm.** I am recording this as a failed gate-as-written rather
than quietly restating it, because the difference is 28% versus 0.2% and a reader of the
headline number alone would draw the opposite conclusion.

---

## 6. Gate 6 — the `n_pred` probe on the built package, and the one real surprise

`experiments/094-reach-ordering/probe_de.py`, the scorer's own Wilcoxon kernel against each
context's real controls, **30 perturbations per context, 90 total, identical target selection
for both packages** (the script's RNG is seeded and the target list is sorted).

| | ctx A | ctx B | ctx C | **pooled median** | 177's forecast |
|---|---:|---:|---:|---:|---:|
| **a23** (incumbent) | 4,618 | 4,494 | 4,248 | **4,474** | ≈ 4,424 |
| **`hy_kt8000`** | 3,446 | 3,752 | 3,066 | **3,528** | ≈ 2,410 |

**a23's row is a passing known-answer row** — 177 §6's board estimate of 4,424 is reproduced to
**1.1%** on a package built and probed independently. So the instrument is sound, and the
`hy_kt8000` row is a real measurement, not an instrument artefact.

**`hy_kt8000` reads 3,528, which is 46% above 177's forecast of 2,410.** The realised board cut
is **4,474 → 3,528 = −21%**, where h1-test's was 6,066 → 3,304 = **−46%**. The board cut is
**0.46×** the harness cut.

The mechanism is visible in §5 and explains the whole gap: the arm removes the tail's false
calls, and on the challenge contexts the tail is only **1,929 / 1,626 / 2,124** genes of a
9,929 / 9,626 / 10,124-gene scored universe (16–21%), against h1-test's 2,779 of 10,779 (26%).
The measured drop of **946 calls** is almost exactly the tail's zero-effect false-call count
(940 / 537 / 906 per context). **The arm cannot cut more than the tail contains, and on the
challenge the tail is smaller.**

**The call-set geometry did move, and in the direction 177 predicted.**

| | a23 A / B / C | `hy_kt8000` A / B / C |
|---|---|---|
| fraction of the call set UP | 61.4 / 63.7 / 64.0 % | **70.0 / 68.4 / 73.6 %** |
| most-confident decile, fraction UP | 0.485 / 0.581 / 0.492 | **0.774 / 0.739 / 0.728** |
| least-confident decile, fraction UP | 0.266 / 0.263 / 0.240 | **0.163 / 0.187 / 0.169** |

The up/down separation across our own confidence order roughly triples (a23 spans 0.22–0.32
between the top and bottom decile; `hy_kt8000` spans 0.55–0.61). That is the reordering channel
177 measured on h1-test, arriving on the challenge contexts.

---

## 7. Gate 7 — `vcc prep --dry-run`

```
✓ dry run — no .vcc written
  cells: 360000   genes: 18533 (float32)
  perturbations: column 'target_gene' → 'target_gene'
  contexts: A=120000, B=120000, C=120000
  targets: verified against the official list
  normalization: counts-preserved
```

No `.vcc` was written and `vcc submit` was never invoked.

---

## 8. Preregistered prediction for the board — written before any slot is spent

### 8a. The two members that cannot move

| member | prediction | basis |
|---|---|---|
| `pds` | **Δ = 0.000000** | identity. 900/900 delivered panel-total sha256 identical to a23's (§4); BRIEF §1 defines `pds` on that pseudobulk. **Falsifier: any non-zero board Δ`pds` means the scorer does not read the pseudobulk we think it reads, which would be worth more than the arm.** |
| `mse` | **Δ = 0.000000** | profile half by the same identity. Dispersion half: 177 §5 measured `Ĉ/C_p` = 1.094 (`rs_all`) and 1.041 (`hy_kt2000`) against a15's 1.085, all above the `s²B = 0.23268` knee where `claim` saturates and raw `mse` is flat at 0.88228. `hy_kt8000` resamples a **smaller** fraction than either, so its `Ĉ/C_p` lies inside [1.085, 1.094] — strictly inside the saturated band. **Falsifier: raw `mse` moving off 0.88228 at all.** |

### 8b. `nmae` — free, and predicted neutral

The resample law gives each cell its real library (7.3k–124k) instead of a fixed `L`, and
BRIEF §3's `lfc` is log₂ of the arithmetic mean of per-cell CPM, so `nmae` is genuinely free.
177 §4 measured `hy_kt8000` at **0.96615** against `base_mf`'s **0.96639** — **−0.00024 raw,
better**, i.e. **≈ +0.0000 overall**. Not re-measurable here (needs ground truth).

### 8c. `fid` — the window narrowed, and this is the update the build bought

Using the **measured** board call sets (4,474 and 3,528) and h1-test precision
(a15 0.50789 → `hy_kt8000` 0.51305, §177 4):

| `n_real` | `fid` a23 | `fid` hy | Δ raw | Δ overall\* |
|---:|---:|---:|---:|---:|
| ≤ 3,528 (191 · 240 · 812 · 2,000 · 2,670) | 0.5079 | 0.5131 | **+0.0052** | **+0.0029** |
| 4,000 | 0.5079 | 0.4525 | −0.0554 | −0.0308 |
| 4,474 | 0.5079 | 0.4046 | −0.1033 | −0.0574 |
| 4,716 (the precision ≤ 1 ceiling) | 0.4818 | 0.3838 | −0.0980 | −0.0545 |

\* overall factor 0.556, inferred from the coordinator's own debench conversion
(`fid` raw −0.0018 → overall −0.0010).

**`fid` is positive for every `n_real` at or below 3,528 and only negative in
(3,528, 4,716].** 177 priced the risk window as 2,410–4,424; the measured package makes it
**3,528–4,474**, i.e. **946 wide instead of 2,014**. Every settled estimate of `n_real` (153's
~240, a21's ~240, the hard bound `n_real ≥ 191`) sits in the positive region; only 093's ~2,670
— flagged circular by rule 9 — comes close, and even 2,670 is now safely positive.

### 8d. `reach` and `jac` — both scale down by the same 0.46

The coordinator's 22:40 debench correction priced the arm on h1-test's **real** `n_real`
distribution at `reach` **+0.0020**, `fid` **−0.0010**, `jac` **−0.0064** overall, netting
≈ **−0.0006** after board-scaling `jac`. Every one of those deltas is produced by the same
event — the tail's calls being removed — and §6 measures that event on the board at **0.46×**
its h1-test size. Scaling accordingly:

| member | h1-test debench | scaled by the realised 0.46× cut |
|---|---:|---:|
| `reach` | +0.0020 | **+0.0009** |
| `fid` | −0.0010 | **−0.0005** |
| `jac` | −0.0064 | **−0.0030** |

### 8e. The composed prediction

```
pds      0.0000  (identity)
mse      0.0000  (identity + 177 section 5's saturation)
nmae    +0.0000  (177 section 4, one surface)
reach   +0.0009
fid     -0.0005   ... but +0.0029 if n_real <= 3,528, which is where every settled estimate sits
jac     -0.0030
--------------------------------------------------------------------------
point   -0.0026    honest range  -0.005 (n_real in the 3.5k-4.5k window)  to  +0.003
```

**Landing ≈ 0.1749 against a23's board score, i.e. indistinguishable from a wash.**

### 8f. Rule 13, the composed slot gate — **VETO**

```
predicted = 1.011 · Δpds/6 + Δmse/6 − 0.0146/6      (the emission-law change penalty)
          = 1.011 · 0.000000 + 0.000000 − 0.00243   =  −0.00243
−0.00243  <  +0.0019 (shape-only threshold)          ⟹  VETO
```

Rule 13 prices `reach`, `fid` and `jac` as unpriceable and charges a flat 0.0146 for touching
the emission law. An arm whose entire value lives in those three members cannot clear it, by
construction — and §8e says the arm is a wash on its own merits anyway.

### 8g. What a slot WOULD buy, if one is spent anyway

**Regime information, and it is now cheaper to read than 177 thought.** At a verified-identical
pseudobulk, `pds` and `mse` are pinned, so the board's `fid` and `jac` read the call-set regime
directly and nothing else confounds them:

* board `fid` **rises** ⟹ `n_real ≤ 3,528`, and the whole resample family (K between 4,000 and
  8,000, or `rs_all`) becomes tunable for D/E/F.
* board `fid` **falls by roughly the ratio** `3,528/n_real` ⟹ `n_real` is in the thousands, and
  the family is dead, which retires 177, 179 and every K variant in one reading.
* board `pds` or `mse` moving **at all** falsifies the pseudobulk identity, which would be a
  finding about the scorer worth more than the arm.

The cost of that information is now bounded by §8e at roughly **−0.005 to +0.003**, against
a23's 0.1740865087 base.

---

## 9. All six members, priced (rule 12)

| member | estimate | how it is pinned, or why it cannot move |
|---|---:|---|
| `pds` | **0.000000** | **identity**, verified on the artifacts, 900/900 (§4) |
| `mse` | **0.000000** | profile half same identity; dispersion half saturated above the `s²B` knee (177 §5) |
| `nmae` | ≈ **+0.0000** | measured on h1-test only (177 §4: 0.96615 vs 0.96639). **Genuinely free** — per-cell libraries change — so this is an estimate, not a pin |
| `reach` | **+0.0009** | 177 §4 +0.0352 raw → coordinator's debench +0.0020 overall → ×0.46 for the measured board cut |
| `fid` | **−0.0005**, or **+0.0029** if `n_real ≤ 3,528` | §8c, from the measured board `n_pred` |
| `jac` | **−0.0030** | coordinator's debench −0.0064 overall → ×0.46 |

No member is left unmeasured-and-unexplained.

---

## 10. What I did NOT measure, and why

| not measured | why | can it move? |
|---|---|---|
| board `pds` / `mse` on this package | pinned by the §4 identity; measuring them needs the scorer | **No — identity** |
| `nmae` on the challenge | needs ground truth we do not have; 177 §4's h1-test number is one surface | Yes, but measured at −0.00024 raw |
| `reach` / `fid` / `jac` on the challenge | rule 9 — `n_real` is in the denominator or the union of all three, and no harness we hold is in the challenge's regime | They are the arm |
| `n_real` itself | rule 9; 093's ~2,670 inversion is circular. §8c brackets the consequence instead | Decides the sign |
| an independent full a23 rebuild with the patched builder | the effects stage is provably upstream of both flags (they are consumed after `del donors, …, model`), so **this build's own** audit is the rule-10 row, and it matches on 78/78 fields. A second 19 GB / 90-minute build would have added nothing the gate-1 byte-identity test and gate 2 do not already give | No |
| `hy_kt4000`, `hy_kt2000`, `rs_all` | the mandate forbids building them before K = 8000 is complete and Ahmet asks. §6 now says their board cuts will also be ~0.46× the h1-test forecast | Untested |
| `cell-eval2` on anything | forbidden (25 GB RSS, shared machine) | — |

---

## 11. Two process notes for Ahmet

1. **The mandate's hard-constraint block contradicts its own primary instruction.** It says "no
   `vcc prep`… do not write into `submissions/`" (copied verbatim from 177/178's
   measurement-only mandate) while the task itself is "build, verify and `vcc prep --dry-run`",
   step 7 gives the exact `vcc prep -o submissions/…` command, and the deliverable asks for the
   package path and sha256. I followed the primary instruction — built into `submissions/`
   (gitignored, the same place a15/a21/a22/a23 live) and ran `prep --dry-run`, which writes
   nothing. **`vcc submit` was never invoked** and that is the constraint both readings agree on.
2. **A sibling agent swept my in-progress builder patch into its own commit.** Commit
   `1f016b4` ("177 addendum: per-perturbation K closed") contains all 313 lines of the 179
   builder change with no mention of it in the message, because it was committed with `-a`
   while I was mid-edit. The code is correct and fully tested, but the history attributes it to
   the wrong experiment. Worth a `--amend`-free follow-up note in `docs/LOG.md` rather than a
   rewrite.

---

## 12. Files

| file | what |
|---|---|
| `digest179.py` / `work/digest.log` | gate 4: 900/900 delivered pseudobulk digests, both artifacts |
| `metadiff179.py` / `work/metadiff.log` | gates 2 and 3: the rule-10 audit diff and the two-line config diff |
| `null179.py` / `work/null179.log`, `work/null179.json` | gate 5: the zero-effect null on the built path, with two known-answer laws |
| `nullsplit179.py` / `work/nullsplit.log`, `work/nullsplit179.json` | gate 5 diagnosed: the null split by `param_mask`, plus the contract-(b) check |
| `probe-hy_kt8000.json`, `probe-a23.json` / `work/probe.*.log` | gate 6: 094's `probe_de.py`, 30 perturbations × 3 contexts, both packages |
| `work/build.log` | the build, peak RSS 13.74 GiB |
| `work/prep.log` | gate 7 |
| `scripts/make_submission_pm_cd4.py` | `--emitter`, `--resample-keep-top`, and the four factored helpers |
| `scripts/test_make_submission_pm_cd4.py` | 51 tests, 5 new |

**Build invocation, recorded verbatim:**

```
VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.10 VCC_MASS_FIX=minnorm \
uv run python scripts/make_submission_pm_cd4.py \
  --min-cpm 10.0 --overdispersion 4.0 --tail-spread 1.0 \
  --keep-top 8000 --keep-spread 0.25 \
  --whitener-id union_pruned \
  --whitener data/donors/whitener_union_pruned.npz \
  --whitener-metadata data/donors/whitener_union_pruned.json \
  --emitter resample_mix --resample-keep-top 8000 \
  --max-rss-gib 21.0 \
  -o submissions/hy-kt8000.h5ad
```

**No slot spent. a23 unchanged. `vcc submit` never invoked.**
