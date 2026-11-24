# 210 — `cd4f034` BUILD RESULT: **BUILDABLE**

207's dose arm, built to the preregistration at `PREREGISTRATION.md`, one variable against the
incumbent `rc_callset`: **`VCC_CD4_TARGET_F` 0.10 → 0.34**. The known-answer row reproduces the
incumbent **to the byte**, and every falsifier checkable before the board **PASSES** — with one
number that came out off its predicted value and that the coordinator must read before submitting
(§5). The package is built, prepped-clean and **NOT submitted**; no model name and no description
were chosen (`AGENTS.md` §14 — the coordinator's call).

```
VCC_CD4_TARGET_F = 0.34      (was 0.10; the builder solves it to fraction 0.32207201067578145)
VCC_ALPHA        = 0.982     unchanged, NOT re-solved (207 §6.2: alpha* = 0.982 exactly)
emitter          = rc_callset            unchanged
everything else  = rc_callset, unchanged
```

| | |
| --- | --- |
| package | `submissions/210-cd4f034.h5ad` |
| **sha256** | **`7036e89f8a03c82d69e79a5c6690b275adae52dbb3ff72bb6ab23e7ca622bf7a`** |
| bytes | 12,648,807,800 |
| meta | `submissions/210-cd4f034.h5ad.meta.json` |
| peak RSS | **19.833 GiB** (21,295,775,744 B; mandate ≤ 22 GB), 1,288.45 s wall, built serially |
| known-answer rebuild | `submissions/210-ka-rc_callset-f010.h5ad`, 12,648,807,800 B, 19.833 GiB, 1,408.60 s |
| `.vcc` | dry-run only, **not written** |

**No code changed.** `src/` and `scripts/` are untouched; the arm is one environment variable.
`experiments/210-cd4-build/build210.sh` is `build202.sh` with three executable lines different —
the dose (the variable) and the output-path plumbing the known-answer row needs:

```
3a4
> OUT="${2:-submissions/210-${ARM}.h5ad}"
12c13
< export VCC_CD4_TARGET_F=0.1
---
> export VCC_CD4_TARGET_F="${VCC_CD4_TARGET_F:-0.34}"
24c25
<   -o "submissions/202-${ARM}.h5ad"
---
>   -o "$OUT"
```

The output path does not reach the emitted bytes — 191 and 202 wrote `rc_hybrid` to two different
paths and got one sha256, and the known-answer row below re-proves it here.

---

## 1. Known-answer rows

### Row 1 — the incumbent must rebuild through this script, bit for bit

`VCC_CD4_TARGET_F=0.1 build210.sh rc_callset` re-runs 202's exact configuration through this
script. Verified by hashing the files on disk, not by reading a log:

```
2a894aa543db7f574ebc024dbb3ebd5e36f923827c4e9a01f263092fe49692e0  submissions/210-ka-rc_callset-f010.h5ad
2a894aa543db7f574ebc024dbb3ebd5e36f923827c4e9a01f263092fe49692e0  submissions/202-rc_callset.h5ad
```

Both **12,648,807,800 bytes**; the digest is the mandate's, character for character.
**PASS.** This environment is 202's environment and everything below is readable.

*(Worth recording: the builder's own `[rss]` instrumentation reported **19.83 GiB** on both of
today's runs against 202's recorded **16.95 GiB** for the same configuration — a 2.9 GiB offset
that is environmental and not a code path, which is exactly what a byte-identical artifact
proves. Do not treat 202's `peak_rss_gib` as a reproducible quantity.)*

### Row 2 — the profile instrument, against 202's published answer

`profile210.py` run on 202's own pair (`202-rc_callset` vs `202-rc_hybrid`), whose answer 202
published, before it was pointed at the arm (`work/profile210-ka.log`):

| quantity | 202 published | measured here | |
| --- | --- | --- | --- |
| per-gene panel-total digests identical | 900/900 | **900/900** | PASS |
| nonzero-structure digests identical | 0/900 | **0/900** | PASS |
| digest value `(A, ABCD1)` | `bb6befc01caf7493` | **`bb6befc01caf7493`** | PASS |
| digest value `(C, ZRANB3)` | `3c9d8f92f3fa4d93` | **`3c9d8f92f3fa4d93`** | PASS |
| per-cell library totals (the weld) | min 17,695 max 22,323 sd 0.394 | **17,695 / 22,323 / 0.394** | PASS |
| \|Δ log2FC\| over the DE universe | max 8.2e−05, mean 2.0e−06 | **max 8.2e−05, mean 2.0e−06** | PASS |
| ‖ΔE‖/‖E‖, delivered effect | (identity) | **0.000000**, 900/900 blocks exactly 0 | PASS |

The instrument reads the profile-**preserving** case as exactly zero displacement, so a non-zero
reading on the arm is the arm and not the tool.

### Row 3 — the DE probe, against 202's published answer

`probe210.py` is `202/probe191.py` verbatim except for one guard (marked `210 CHANGE`): 202's
copy aborts when a package's `effects.final_sha256` differs from the cached `eff.a15.npy`, which
is correct for a profile-preserving arm and would abort on *this* arm working. In the same pass
that reads the arm, the probe re-reads the incumbent:

| | 202 published | measured here |
| --- | ---: | ---: |
| `rc_callset` `n_pred` A / B / C | 2,599 / 2,599 / 2,599 | **2,599 / 2,599 / 2,599** |
| `rc_callset` ρ(−log p, \|H\|) mean over contexts | +0.9993 | **+0.9993** |

**PASS.** The instrument reading the arm is the instrument that produced the baseline.

---

## 2. The six falsifiers

| # | falsifier | measured | verdict |
| --- | --- | --- | --- |
| 1 | rule-10 config diff is exactly the pre-declared one | **1 config key, 7 effects keys, all pre-declared** | **PASS** |
| 2 | every meta value in 207 §8.2's table | **11 of 11 exact, plus the dose algebra** | **PASS** |
| 3a | panel-total digests must **differ** on 900/900 | **900/900 differ**, structure 0/900 identical | **PASS** |
| 3b | delivered median row norm ≈ 1.0009×, `A` within 1.0% | **0.99367×** and **0.97970× (−2.03%)** | **FAIL — see §5** |
| 4 | non-vacuity: the delivered effect actually moved | **‖ΔE‖/‖E‖ = 0.235141**, per-block min **0.0706**, zero blocks at 0 | **PASS** |
| 5 | challenge `n_pred` within ±1% of 2,599 | **2,599** in every context, **+0.00%** | **PASS** |
| 6 | `vcc prep --dry-run` clean; nnz against the cap | clean; **33.2154%** | **PASS** |

### Falsifier 1 — the rule-10 config diff, verbatim

`metadiff210.py`, both expected sets pre-declared in source before the diff was read. 202's arm
changed the emitter, so its effects audit had to be **zero** differences; this arm changes the
**profile**, so the audit must differ in exactly the deterministic consequences of the dose and
in nothing else.

```
=== A: effects + operator audit, cd4f034 vs rc_callset (rule 10) ===
  fields compared: 78   differences: 7
    expected      effects.cd4.fraction                                      0.32207201067578145          vs 0.09472706196346513
    expected      effects.cd4.post_add_preserve.candidate_pre_preserve_median_norm4.479655759355644            vs 4.1665533562021135
    expected      effects.cd4.post_add_preserve.fraction                    0.32207201067578145          vs 0.09472706196346513
    expected      effects.cd4.post_add_preserve.post_add_global_preserve_scale0.9203591354319881           vs 0.989521013951823
    expected      effects.cd4.residual_frame.solved_fraction                0.32207201067578145          vs 0.09472706196346513
    expected      effects.cd4.residual_frame.target_f_requested             0.34                         vs 0.1
    expected      effects.final_sha256                                      87a4b6d181298f63826dc9caeb   vs 22c09f22c749bda82d12a0862b
  A PASS

  PINNED audit fields (must be EXACT):
    effects.baseline_pre_sha256                                   7678bab9ea9d1b8d951b0b   vs 7678bab9ea9d1b8d951b0b   EXACT
    effects.candidate_pre_sha256                                  0601336fdf1ec5c533dcad   vs 0601336fdf1ec5c533dcad   EXACT
    effects.post_whiten_sha256                                    b2426e15d76d9fe59f794a   vs b2426e15d76d9fe59f794a   EXACT
    effects.whitening.after_all_median_norm                       4.229295604454299        vs 4.229295604454299        EXACT
    effects.whitening.before_all_median_norm                      4.238101079325155        vs 4.238101079325155        EXACT
    effects.whitening.after_operator_median_norm                  25.73792992728808        vs 25.73792992728808        EXACT
    effects.whitening.before_operator_median_norm                 4.122892101713494        vs 4.122892101713494        EXACT
    effects.whitening.operator_preserve_scale                     0.16018740098217019      vs 0.16018740098217019      EXACT
    effects.whitening.gamma                                       0.65                     vs 0.65                     EXACT
    effects.whitening.lambda                                      0.0                      vs 0.0                      EXACT
    effects.cd4.fully_supported_targets                           276                      vs 276                      EXACT
    effects.cd4.supported_entries                                 2712528                  vs 2712528                  EXACT
    effects.cd4.mode                                              gene_norm                vs gene_norm                EXACT
    effects.cd4.panel_residual_nonzero_entries                    0                        vs 0                        EXACT
    effects.cd4.residual_frame.cd4_fraction_constant              0.05                     vs 0.05                     EXACT
    effects.cd4.residual_frame.effective_f_at_cd4_fraction        0.05278322684523277      vs 0.05278322684523277      EXACT
    effects.cd4.residual_frame.f_per_unit_fraction                1.0556645369046553       vs 1.0556645369046553       EXACT
    effects.cd4.residual_frame.fraction_source                    VCC_CD4_TARGET_F         vs VCC_CD4_TARGET_F         EXACT
    effects.cd4.residual_frame.whitened_base_median_row_norm      4.229295604454299        vs 4.229295604454299        EXACT
    effects.cd4.residual_frame.residual_median_row_norm_at_unit_fraction4.464717385709141        vs 4.464717385709141        EXACT
    effects.cd4.post_add_preserve.baseline_output_median_norm     4.122892101713487        vs 4.122892101713487        EXACT
    operator.sha256                                               e2f65ecf69b7e0b0b817d0   vs e2f65ecf69b7e0b0b817d0   EXACT
    operator.columns_sha256                                       b2cae0d381c3cd7065a3de   vs b2cae0d381c3cd7065a3de   EXACT
    operator.metadata.sha256                                      e2f65ecf69b7e0b0b817d0   vs e2f65ecf69b7e0b0b817d0   EXACT
    operator.metadata_sha256                                      838f263c6a4e1bed069d5f   vs 838f263c6a4e1bed069d5f   EXACT
  PINNED PASS

=== B: config diff (rule 10's intended-difference count) ===
  cd4_target_f               cd4f034 0.34                   rc_callset 0.1
  config keys that moved: ['cd4_target_f']
  pre-declared expected : ['cd4_target_f']
  B PASS (1 differences)
```

**One config key moved.** `alpha` 0.982, `emitter` `rc_callset`, `whitener_id` `union_pruned`,
`min_cpm` 10, `tail_spread` 1.0, `keep_top` 8000, `keep_spread` 0.25, `overdispersion` 4.0,
`cd4_fraction` 0.05 (the legacy knob, untouched), and all four `rc_*` keys and all three
`rc_callset_*` keys are the incumbent's. **`operator.sha256` is
`e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123`, identical.
`effects.whitening.after_all_median_norm` is `4.229295604454299` on both sides — unmoved.
`effects.cd4.fully_supported_targets` is 276.** So preregistration falsifier 6's four named
conditions all hold.

The seven audit differences are the dose and nothing else. `effects.post_whiten_sha256` — the
whitened base, hashed **before** the CD4 add — is **identical**, which is the sharpest available
proof that the arm touched nothing upstream of the residual. `effects.final_sha256` **must**
differ and does; a `final_sha256` that had not moved would mean the residual never applied.

Preregistration falsifier 6 says "exactly 2 fields (`cd4_target_f` and its deterministic
consequence `effects.cd4.fraction`/`solved_fraction`)". Read as *knobs*, that is what happened:
**one** environment variable moved and **one** derived number, the solved fraction, which the
meta records under three names (`cd4.fraction`, `cd4.post_add_preserve.fraction`,
`residual_frame.solved_fraction`) and which drags two arithmetic consequences with it (the
pre-preserve norm and the preserve scalar it implies) plus the final hash. Every one of those
seven is entailed by the single knob; none is a second variable. **Adjudicated PASS**, with the
count stated as 1 + 7 rather than 2 so nothing is hidden.

### Falsifier 2 — 207 §8.2's required meta values, on the built meta

```
  config.cd4_target_f                                    0.34                    want 0.34                  OK
  config.alpha                                           0.982                   want 0.982                 OK
  config.emitter                                         rc_callset              want rc_callset            OK
  config.whitener_id                                     union_pruned            want union_pruned          OK
  effects.cd4.fully_supported_targets                    276                     want 276                   OK
  effects.cd4.residual_frame.target_f_requested          0.34                    want 0.34                  OK
  effects.cd4.residual_frame.fraction_source             VCC_CD4_TARGET_F        want VCC_CD4_TARGET_F      OK
  effects.cd4.residual_frame.effective_f_at_cd4_fraction 0.05278322684523277     want 0.05278322684523277   OK
  effects.cd4.residual_frame.f_per_unit_fraction         1.0556645369046553      want 1.0556645369046553    OK
  effects.whitening.after_all_median_norm                4.229295604454299       want 4.229295604454299     OK
  operator.sha256                                        e2f65ecf69b7e0b0b817…   want e2f65ecf69b7e0b0…     OK
  solved_fraction == target_f / f_per_unit: 0.32207201067578145 == 0.34/1.0556645369046553   OK
  dose ratio delivered: 0.322072011 / 0.094727062 = 3.400000x   (requested 0.34/0.10 = 3.4x)
```

`solved_fraction` is **0.32207201067578145**, i.e. **0.322** — 207 §8.2 predicted "≈ 0.32" on the
challenge panel and the builder landed there from `f_per_unit = 1.0556645369046553`, unchanged.
`effective_f_at_cd4_fraction` is unmoved at 0.05278322684523277, confirming it is a diagnostic of
the legacy `VCC_CD4_FRACTION` and not of this arm, exactly as the mandate says.

**Two further algebra checks the mandate did not ask for, because "solved" is not "delivered":**

```
                     residual median row norm    f delivered    achieved post-preserve median row norm
arm  f = 0.34         1.4379605055144615          0.34           4.122892101713487
inc  f = 0.10         0.42292956044542984         0.09999999999999999   4.122892101713487
```

The **delivered** `f` in 099's own definition (residual median row norm over whitened-base median
row norm) is **0.34 exactly**, and the achieved post-add median row norm in effect space is
**bit-identical** to the incumbent's, 4.122892101713487 on both. 207 §6.2's claim 1 — "the
preserve scalar pins the effect-space median row norm at every dose" — holds here to sixteen
digits. This matters for §5: what moves is downstream of the renderer, not the effect.

### Falsifier 3a — the delivered profile moved, on 900/900 blocks

`profile210.py`, streaming both artifacts (`work/profile210.log`):

```
per-gene PANEL TOTAL digest identical on   0/900   (MUST be 0 -- this arm moves the profile)
panel total SUM identical on               12/900
nonzero-STRUCTURE digest identical on      0/900
```

**Zero blocks match.** This is the exact inverse of 202's row, on the same instrument that read
202's 900/900 half an hour earlier. The 12/900 agreeing panel-total *sums* are the scalar total
count per block colliding, which is expected — the weld pins each cell's library to `T = S//400`,
so the block sums live in a narrow integer range.

The weld survives: per-cell library totals min **17,695**, max **22,323**, sd median **0.398**
(incumbent 0.394). The 1,000,000 format cap is **45× away**.

### Falsifier 4 — non-vacuity, the number most likely to catch a silent failure

200 and 207 both found the CD4 residual can be identically zero while every other check passes,
because the sizing rule takes a median. Measured on the delivered effect
`E = log1p(TS·p̂) − b_ctrl`, `TS = 5e4`:

```
NON-VACUITY  ||E_arm - E_inc||_F / ||E_inc||_F = 0.235141
  per-block relative displacement:  min 0.070633   p01 0.071448   median 0.230833   max 0.794375
  blocks with displacement exactly 0: 0/900        blocks below 1e-6 relative: 0/900
  ctx A 0.238384   ctx B 0.233086   ctx C 0.233828
  restricted to the builder's 10-CPM gate: 0.235270  (per-block min 0.070602)
```

**Every one of the 900 blocks moved, the smallest by 7.06%.** Not a single silent zero, and the
minimum is what proves it — a median or a mean would hide one.

**And it landed on 207's own dose-response line.** 207 §"non-vacuity" measured
`‖delivered(f) − delivered(0)‖/‖delivered(0)‖` at **0.1023** for `f = 0.10` and **0.3312** for
`f = 0.34`. The residual is linear in the fraction, so the arm-vs-incumbent displacement should be
`0.3312 − 0.1023 = 0.2289`. **Measured 0.2351**, +2.7%, on a different frame (the shipped
challenge artifacts, post-emission) with an independently derived integer-quantisation inflation
of ~0.2% in the squared norm. The dose that was requested is the dose that was delivered.

**Call-set exposure.** `rc_callset` prescribes `top 400 by |H| ∪ top 2,200 by control CPM`. The
CPM half is a function of the control cells alone and cannot move with the dose, so the arm's
entire call-set exposure is the `|H|` head. Measured head overlap with the incumbent: **min 184,
median 314, max 399 of 400**, implying a call-set jaccard **≥ 0.9360 at the median block** and
≥ 0.8466 at the worst. The set moves, but by ≤ 216 genes of 2,600 in the worst block.

### Falsifier 5 — challenge `n_pred`

`probe210.py`, 40 targets × 3 contexts, Wilcoxon + BH against **the real control cells the scorer
will use**, target row excluded:

| | A | B | C | pooled |
| --- | ---: | ---: | ---: | ---: |
| `cd4f034` | **2,599** | **2,599** | **2,599** | **2,599** |
| `rc_callset` | 2,599 | 2,599 | 2,599 | 2,599 |

**+0.00%**, against a ±1% bound. The call set does not move in size. The builder's own diagnostics
over all 900 blocks agree and are *tighter* than the incumbent's: `call_set_prescribed`
2600/2600/2600 and `call_set_missed` 0/0/0 are **identical** to the incumbent's (they do not
appear in the meta diff at all), while `call_set_delivered.max` fell **2,752 → 2,689**,
`call_set_unsuppressible.max` **152 → 89** and `call_set_jaccard.min` rose **0.9448 → 0.9669**.
Nothing prescribed is ever missed and the unsuppressible floor got smaller.

Fraction of the call set called UP: **49.8 / 49.9 / 48.0** against the incumbent's
49.7 / 49.8 / 47.8 — unchanged.

**On the delivered order.** ρ(−log p, |H|) inside the call set reads **+0.9401 / +0.9382 /
+0.9353, mean +0.9379** for the arm against the incumbent's +0.9993. **This is a cross-key
number, not an order regression**, and it must not be read as one: `eff.a15.npy` is the
`f = 0.10` key, so the arm is being scored against the *incumbent's* ordering key while
delivering its own. The drop is a measurement of how far the key moved with the dose — the same
thing the 314/400 head overlap measures — and the arm's own order is |H|-monotone by exactly the
construction 202 verified (`unsuppressible` median 0, `missed` 0, so the delivered set *is* the
prescribed set). The arm's own effect matrix is not written to disk by the builder, so no
same-key ρ is available without a code change, which the mandate forbids and which is not needed:
ρ against the universe also stays high (+0.7072 / +0.7188 / +0.7209 against +0.7536 / +0.7612 /
+0.7648) and ρ against expression is unchanged (−0.268 / −0.259 / −0.276 against −0.255 / −0.255
/ −0.273).

### Falsifier 6 — format contract and density

```
vcc prep submissions/210-cd4f034.h5ad -g data/challenge/gene_names.csv \
  --perts data/challenge/pert_counts.csv -o submissions/210-cd4f034.vcc --dry-run

✓ dry run — no .vcc written
  cells:         360000
  genes:         18533  (encoding: float32)
  contexts:      A=120000, B=120000, C=120000
  targets:       verified against the official list
  normalization: counts-preserved
```

| | `cd4f034` | `rc_callset` |
| --- | ---: | ---: |
| stored entries A / B / C | 538,693,542 / 512,429,392 / 526,606,559 | 539,067,895 / 512,298,146 / 526,660,526 |
| **nnz** | **1,577,729,493** | 1,578,026,567 |
| **% of the 4.75e9 cap** | **33.2154%** | 33.2216% |
| nnz per cell | 4,382.6 | 4,383.4 |
| max cell total | 22,323 | 22,323 |
| bytes | 12,648,807,800 | 12,648,807,800 |

**The two files are exactly the same size and are not the same file.** `X/data` holds
1,577,729,493 float32 against 1,578,026,567 — a 297,074-element difference, which is 0.28 of the
dataset's 1,048,576-element HDF5 chunk and is absorbed by the tail chunk's allocation. The
digests, the nnz and the sha256 all differ; only the on-disk footprint collides. Recorded because
a byte count is the first thing anyone checks and here it is uninformative.

---

## 5. ⚠️ THE ONE NUMBER THAT CAME OUT OFF-PREDICTION — the delivered amplitude

The mandate's falsifier 3 asks for **delivered median row norm ≈ 1.0009×** and **delivered `A`
within 1.0%**, both from 207 §6.2. On the shipped artifacts:

| | 207 predicted | **measured on the shipped pair** |
| --- | ---: | ---: |
| delivered median row norm | 1.0009× (0.99987× on amp207's instrument) | **0.99367×** (−0.63%) |
| delivered `A = mean ‖row‖²` | 0.99036× (−0.96%) | **0.97970×** (**−2.03%**) |

**The `A` bound is missed by 2×, and the drop is real, not an emission artefact.** `amp210.py`
measures both statistics over five nested gene sets, from the noisiest to the most
quantisation-immune, with the integer-emission noise energy computed analytically per set:

```
gene set                            n (A/B/C)   median|row| inc     arm   ratio    A inc     arm   ratio  q-noise/A
all genes                   18533/18533/18533          3.454774 3.432898 0.99367 14.708266 14.409672 0.97970    0.22%
builder 10-CPM gate            9135/8668/9127          3.447739 3.425875 0.99366 14.656333 14.357508 0.97961    0.02%
top 2,000 by control CPM       2000/2000/2000          1.714140 1.687206 0.98429  3.634070  3.518322 0.96815    0.00%
panel total >= 1,000 counts    2343/2186/2209          1.812379 1.786859 0.98592  4.051548  3.927502 0.96938    0.00%
panel total >= 10,000 counts         26/36/39          0.220573 0.217898 0.98787  0.075873  0.071160 0.93788    0.00%
```

Quantisation adds a fixed energy to **both** rows and therefore biases the ratio **towards 1**; it
cannot manufacture a drop, and it is 0.22% of `A` on the full axis (analytic bound, matching the
0.2% derived from `|dt| ≤ 1` count and `TS/T = 6.25e-3`) and ~0% on the high-count sets. The ratio
**survives and deepens** on the quantisation-immune sets. The drop is in the delivered profile.

**The likely mechanism, stated as a hypothesis and not as a result:** `amp207.py` renders with
`deliver(..., "minnorm")` — experiment 156's **mass-fix** renderer — while production ships
`VCC_MASS_FIX` **unset** and uses the implicit renormalisation. On 156's own measurement that
implicit correction is 4.0–4.5% of the row norm in a nearly-shared direction, and it is
mass-ratio dependent, so it moves with the dose. 207 also measured on PLOLO fold estimates rather
than the production panel. Two different rendered objects; production's is the one that ships.
**This is a correction to 207 §6.2, not to the arm.** Note what did *not* move: the effect-space
post-preserve median row norm is bit-identical between the two builds (§2, falsifier 2), so 207's
claim about the preserve scalar is intact — the divergence is entirely downstream of it.

### Two consequences the coordinator owns

**1. Preregistration falsifier 2 may fire on the LOW side, for a favourable, identified reason.**
Applying **207's own conversion** (`amp207.py`: `d_raw = ΔA/A × 0.2331583/0.229141`) to the
measured ratio:

```
  A ratio 0.97970  ->  implied d raw mse from the A channel alone   -0.020657
  207's  0.99036   ->  implied                                      -0.009809
  extra, uncredited, favourable                                     -0.010848 raw
```

207 deliberately did **not** credit the `A` drift, calling it conservative. If the extra −0.0108
transfers, raw `mse` lands near **0.8333** against the point prediction of **0.84413** — **below
preregistration falsifier 2's window of [0.836, 0.852]**. That falsifier would then fire on a
mechanism that is measured, understood and *in our favour*. This is mixing two instruments
(207's ΣD normalisation with a ratio measured on the shipped artifact), so treat it as a flag,
not a re-price — but read it before reading the board, not after.

**2. The "shape-only" classification is no longer measured, it is argued.** 207 §7 put the arm in
rule 13's shape-only class (+0.0019 gate) on the strength of "delivered median row norm 1.0009×,
delivered `A` within 1.0%". Those numbers are 0.99367× and 0.97970×. The arm still changes **no**
amplitude parameter — `alpha` 0.982, `gamma` 0.65, `lambda` 0.0, `effect_scale` 1.0 all identical,
effect-space median row norm pinned to sixteen digits — and 2.03% remains 7× smaller than the 14%
that made 170 §2 raise the trap in the first place. But the literal criterion is not met. The
preregistered +0.00437 clears the shape-only gate by 2.3× and **fails** the amplitude gate
(+0.0089) by 2×, so the class assignment is the whole decision. **I am not authorised to
re-classify a binding preregistration and I have not.** The measurement is here; the call is the
coordinator's.

---

## 6. Preregistration §5's six falsifiers — what is checkable now

| # | pre-committed falsifier | checkable pre-board | status |
| --- | --- | --- | --- |
| 1 | raw `mse` not below `rc_callset`'s 0.8738045 | no — board only | pre-board evidence: the profile moved on 900/900 blocks by 0.235 relative, at the intended dose; the `A` channel moved favourably (§5) |
| 2 | raw `mse` outside [0.836, 0.852] | no — board only | **at risk on the LOW side**, §5 |
| 3 | `pds` falls by more than 0.004 raw | no — board only | — |
| 4 | `nmae` raw rises by more than 0.002 | no — board only | the delivered log2FC moved: \|Δlfc\| mean **0.012889** over the DE universe, median-of-block-means 0.012267, signed block mean median **+0.000166** (range −0.0058 to +0.0071). Against 202's profile-preserving 2.0e−06 this is 6,400× larger, i.e. `nmae` is genuinely in play, as 207 priced it (−0.000152 overall) |
| 5 | `fid`, `reach` or `jac` moves by more than 0.004 scaled | no — board only | `n_pred` **unchanged at 2,599** (+0.00%), call-set jaccard ≥ 0.936 at the median block, up-fraction unchanged, `unsuppressible` floor smaller. The call set did **not** collapse |
| 6 | config diff is not exactly 2 fields, with `operator.sha256`, `fully_supported_targets` = 276, `alpha` = 0.982, `after_all_median_norm` unmoved | **yes, fully** | **PASS** — §2 falsifier 1, all four named conditions exact |

*(One bookkeeping note, immaterial to any threshold: falsifier 1's stated incumbent raw `mse`,
**0.8738045**, is `rc_hybrid`'s board value; `BOARD-RESULT.md`'s table prints **0.87367** for both
`rc_hybrid` and `rc_callset` in the same column. The two disagree in the 4th decimal. Both sit far
above the 0.84413 point prediction, so falsifier 1 is unaffected either way, but the ledger
should be reconciled before anyone inverts a cosine from it.)*

---

## 7. Verdict

**BUILDABLE.** The known-answer row reproduces the incumbent to the byte on 12,648,807,800 bytes;
two measurement instruments reproduce 202's published answers before being pointed at the arm; the
config diff is one environment variable and its seven entailed consequences with the operator sha,
the whitening norm, the whitened-base hash and `fully_supported_targets` all exact; the delivered
profile moved on 900/900 blocks with **no** silent zero and by **0.2351** relative — within 2.7% of
207's own dose-response line; `n_pred` is unchanged at 2,599; and `vcc prep --dry-run` is clean at
33.2154% of the density cap.

**One reservation, flagged not buried:** the delivered `A` fell **2.03%**, twice 207's predicted
0.96% and outside the mandate's 1.0% bound. The drop is real (it survives on quantisation-immune
gene sets), it is favourable for `mse`, and its likely cause is that 207 measured through the
mass-fix renderer that production does not ship. It puts preregistration falsifier 2 at risk on
the low side and it removes the *measured* basis for the shape-only classification that the rule-13
gate turns on. §5 has the numbers; the classification call is the coordinator's.

**Not submitted.** `vcc submit` was never run. No model name and no description were chosen. The
`.vcc` was dry-run only and not written. No file outside `experiments/210-cd4-build/` and
`submissions/` was created or modified; `src/` and `scripts/` are untouched.

### Files

| path | what |
| --- | --- |
| `build210.sh` | the arm's build script — `build202.sh`, one variable |
| `work/build-cd4f034.log` | the arm's build, meta and `/usr/bin/time -l` |
| `work/build-ka-f010.log` | known-answer row 1's rebuild of the incumbent |
| `work/ka-sha256.txt`, `work/arm-sha256.txt` | `shasum` on the files, not read from a log |
| `metadiff210.py`, `work/metadiff210.log` | falsifiers 1 and 2, all five sections |
| `profile210.py`, `work/profile210.log`, `work/profile210.json` | falsifiers 3a and 4 |
| `work/profile210-ka.log`, `work/profile210-ka.json` | known-answer row 2 — the instrument against 202's pair |
| `amp210.py`, `work/amp210.log`, `work/amp210.json` | §5, the delivered amplitude over five gene sets |
| `probe210.py`, `work/probe210.log`, `work/probe210.json` | falsifier 5 and known-answer row 3 |
| `work/prep210.log` | falsifier 6 |
| `submissions/210-cd4f034.h5ad` | **the package** |
| `submissions/210-ka-rc_callset-f010.h5ad` | the known-answer rebuild, byte-identical to `202-rc_callset.h5ad` |
