# 191 — `rc_hybrid` BUILT and all six falsifiers PASS. **BUILDABLE.** No slot spent, no submission.

Date: **2026-09-03**. Build + verification only. **No `vcc submit`. No `cell-eval2` process.**
`vcc prep` was run `--dry-run` only and wrote nothing. Nothing written outside
`experiments/191-rc-hybrid-build/`, `scripts/`, `src/` and `submissions/`.

**Verdict in one line: BUILDABLE. Every pre-committed falsifier passes, the ported emitter
reproduces all fourteen of 185's h1-test known-answer numbers exactly, and the ordering
delivers on the real challenge contexts at Spearman +0.974 against a +0.95 bar.**

Three things the coordinator did not have before and needs before the slot is spent:

1. **`a15_local` is not merely a local reference — it is a15.** The rebuilt package is
   **byte-identical** to the shipped a15 artifact: sha256
   `26fb38341c007661a780d275f00123ca3c911f3c7ea407fc625a34bea1bc85e3` on all
   19,049,370,200 bytes, verified with an independent `shasum`. The mandate's premise
   ("a15's recorded artifact cannot be byte-reproduced here") is **false**: the arm64/x86_64
   drift is real in the effect matrix (`effects.final_sha256` differs) but it is entirely
   absorbed by the integer quantisation `floor(n·L·profile)`. So falsifier 1 is not a
   comparison against a proxy — it is a comparison against a15 itself.
2. **The incumbent's ordering defect is ~2.6× SMALLER on the challenge than on h1-test.**
   Spearman(−log p, `|H|`) inside a15's own call set reads **−0.152** on A/B/C against
   **−0.391** on h1-test, and the expression correlation **+0.279** against **+0.486**. The
   mandate expected ≈ −0.39 and it is not there. `rc_hybrid` still lands at **+0.974**, so
   the delivered *span* is 1.126 Spearman units on the challenge against h1-test's 1.375 —
   **82%**. Whatever the h1-test → board `reach` shrink is made of, part of it is this.
3. **The artifact is 31% sparser** (nnz 2.378 G → 1.647 G, 50.1% → 34.7% of the density cap,
   19.0 GB → 13.2 GB). That is the exact-marginal emitter concentrating each gene's counts
   into fewer cells, it breaks no contract, and `vcc prep --dry-run` passes — but it is a
   visible property of the package and it is stated rather than buried.

---

## 0. The two packages

| | path | sha256 | bytes | nnz | % of density cap | peak RSS |
|---|---|---|---:|---:|---:|---:|
| `a15_local` | `submissions/191-a15_local.h5ad` | `26fb38341c007661a780d275f00123ca3c911f3c7ea407fc625a34bea1bc85e3` | 19,049,370,200 | 2,377,814,190 | 50.06% | 14.81 GiB |
| `rc_hybrid` | `submissions/191-rc_hybrid.h5ad` | `e0b0677f984a94a1c4a7af2add196908cdab90d7acb463eb1d6e0397c2d870a8` | 13,202,460,120 | 1,647,149,337 | 34.68% | 19.83 GiB |

Metas at `<path>.meta.json`. Both built on macOS-27.0-arm64, builder sha256
`90c4f49eebc1c706ab0a3673f5e79688d79d5b73251a06580f76b1a1e529342a`,
`src/vcc2026/rc_emit.py` sha256
`3203e44670eb6ba8d165e7c1594d9adeab185174770aad99bf517d6602c92024`.
`src/vcc2026/donor_transfer.py` is **unmodified** (`9ffca7ec…`, its frozen value).

**`a15_local` reproduces a15's recorded artifact byte for byte.** Every `output` field of the
meta matches a15's recorded meta exactly: `sha256`, `bytes`, `nnz` (2,377,814,190),
`nnz_per_cell` (6605.0394166666665), `pct_of_density_cap`, and all three contexts'
`stored_entries` (816,513,520 / 770,941,878 / 790,358,792).

**Build command** (`experiments/191-rc-hybrid-build/build191.sh`), a15's configuration read
off `submissions/a15-cd4f010.vcc.meta.json` field by field:

```
VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.1 uv run python scripts/make_submission_pm_cd4.py \
  --whitener data/donors/whitener_union_pruned.npz \
  --whitener-metadata data/donors/whitener_union_pruned.json --whitener-id union_pruned \
  --min-cpm 10 --tail-spread 1.0 --keep-top 8000 --keep-spread 0.25 --overdispersion 4.0 \
  [--emitter rc_hybrid] -o submissions/191-<arm>.h5ad
```

The bracketed flag is the **only** difference between the two invocations, and §5 proves it is
the only difference in the two written metas outside `output`.

---

## 1. The falsifier table — all six PASS

| # | falsifier | bound | measured | verdict |
|---|---|---|---|:-:|
| **1** | delivered per-gene panel totals identical, per target per context, by digest | 900/900 | **900/900 identical**; nonzero-structure digest identical on **0/900** | **PASS** |
| **2** | `nmae` within 0.0010 raw of `a15_local` | ±0.0010 | **+0.00011 raw** (h1-test, direct, with truth); challenge-side driver **0.564×** h1-test's ⇒ ≈ **+0.00006** | **PASS** |
| **3** | challenge `n_pred` within ±2% | ±2% | **4,480 vs 4,474 = +0.13%** (40 targets × 3 contexts); full-panel pass-1 median **4,436** | **PASS** |
| **4** | Spearman(−log p, `\|H\|`) in the call set, challenge contexts, ≥ +0.95 | ≥ +0.95 | **+0.9761 / +0.9713 / +0.9755, mean +0.9743** | **PASS** |
| **5** | rule-10 config diff: exactly the emitter-related differences; operator sha256 matches; `effects.whitening.after_all_median_norm` has not drifted | exactly 5 pre-declared keys, 0 elsewhere | **78 effects/operator fields, 0 differences**; config diff **exactly the 5 pre-declared keys**; operator sha `e2f65ecf…` EXACT; `after_all_median_norm` 4.229295604454299 EXACT | **PASS** |
| **6** | `vcc prep --dry-run` passes | clean | **clean on `rc_hybrid`** (and on `a15_local`): 360,000 cells, 18,533 genes, A/B/C = 120,000 each, targets verified, counts preserved | **PASS** |

---

## 2. Gate 0 — the known-answer row. All fourteen numbers reproduce EXACTLY

The mandate's corollary of rule 10: the ported emitter must reproduce 185's h1-test row
before it is trusted on the challenge. `ka191.py` drives the **production** module
(`vcc2026.rc_emit`, the object the challenge build calls) over h1-test with 185's own control
block, a15's own per-gene totals, 169's `|H|`, 185's seeds and 185's spread policy. Scored
with 185's own instruments (`debench.Surface`, `six185.walk/members/spearmans`).

| arm | quantity | measured | 185 | verdict |
|---|---|---:|---:|:-:|
| a15 | `reach` | **0.18283** | 0.18283 | PASS |
| a15 | `n_pred` (debench convention) | **6034** | 6034 | PASS |
| a15 | `jac` | **0.06824** | 0.06824 | PASS |
| a15 | `nmae` | **0.96659** | 0.96659 | PASS |
| a15 | `P(1)` | **0.73171** | 0.7317 | PASS |
| a15 | Spearman(`\|H\|`) in call set | **−0.3906** | −0.391 | PASS |
| a15 | Spearman(expression) in call set | **+0.4856** | +0.486 | PASS |
| **rc_hybrid** | **`reach`** | **0.21681** | 0.21681 | PASS |
| **rc_hybrid** | **`n_pred`** | **6043** | 6043 | PASS |
| **rc_hybrid** | **`jac`** | **0.06827** | 0.06827 | PASS |
| **rc_hybrid** | **`nmae`** | **0.96670** | 0.96670 | PASS |
| **rc_hybrid** | **`P(1)`** | **0.87805** | 0.8780 | PASS |
| rc_hybrid | Spearman(`\|H\|`) in call set | **+0.9837** | +0.984 | PASS |
| rc_hybrid | Spearman(expression) in call set | **+0.0073** | +0.007 | PASS |

Plus a15's raw-BH median `n_pred` **6035** (185 §1's value) and the delivered per-gene
panel-total digest identical to a15's on **85/85** targets for every arm
(`b53d02aebfd0816a…`, 171 §6's digest).

**Seeding.** 185's `realize185.py` opened a fresh `default_rng(90000+i)` for its incumbent arm
and again for its rc arm; the production writer drives all 900 blocks off one stream. Both
were run:

| | 185's two-stream bookkeeping | production's single stream | \|diff\| |
|---|---:|---:|---:|
| `reach` | 0.21681 | 0.21681 | **0.00000** |
| `jac` | 0.06827 | 0.06827 | 0.00000 |
| `nmae` | 0.96670 | 0.96670 | 0.00000 |
| `fid` | 0.50693 | 0.50694 | 0.00001 |

Against 185 §3e's own twin spread for this family (7e−5) and a15's seed spread (0.00746).
**The seeding is not a variable.**

The per-block diagnostics reproduce 185 §2d row for row on h1-test — forced placements **0**,
maximum pre-repair residual **0**, repair moves **0**, unresolved **0**, sub-gate reservoir
**111 counts/cell** (185: 113), designed-vs-realised `z` max **0.0048** (185: 0.0044), per-cell
total sd **0.475** (185: 0.47), `T` = 20,614 on AASS (185: 20,614).

---

## 3. Falsifier 1 — the pseudobulk identity, on the two real artifacts

`digest191.py`, one sha256 of the delivered `int64` per-gene panel totals per (context,
perturbation) block, streamed from both files. Block order verified equal from the obs
categoricals before any comparison; all 900 blocks are 400 rows.

```
per-gene PANEL TOTAL digest identical on   900/900
panel total SUM identical on               900/900
nonzero-STRUCTURE digest identical on        0/900
```

**Both halves matter.** Without the exact-total placement the residual would land in the
totals and the blocks would separate; and if the structure column had matched too, the two
files would simply be the same file and the totals identity would prove nothing. It matches on
**0** of 900. Sample:

```
(A, ABCD1)   totals bb6befc01caf7493 vs bb6befc01caf7493   struct e84ccddbfb81fafa vs 7dcfab48b00cd46d
(A, ACLY)    totals b17f1e11a0126af9 vs b17f1e11a0126af9   struct 6d4656237d7aeef7 vs 8217fd951c710082
(C, ZRANB3)  totals 3c9d8f92f3fa4d93 vs 3c9d8f92f3fa4d93   struct 15365cb91c79d2fe vs 4cacd84fcc550a15
```

**So Δ`pds` = 0.000000 and the profile half of Δ`mse` = 0 are identities on this package**, and
because `a15_local` is byte-identical to a15, they are identities against **a15's own board
row**, not against a local proxy.

**The weld** (same pass): every `rc_hybrid` cell carries one library total, the `S mod 400`
remainder cells one count above it.

| arm | min cell total | max cell total | median per-block sd | max per-block sd |
|---|---:|---:|---:|---:|
| `a15_local` | 17,076 | 25,061 | 342.9 | 543.3 |
| `rc_hybrid` | 17,695 | 22,323 | **0.394** | **0.500** |

Max cell total 22,323 against the 1,000,000 format cap.

**Placement, all 900 challenge blocks** (from the built meta's `output.rc_hybrid`):

| diagnostic | median | max |
|---|---:|---:|
| forced placements | 0 | **0** |
| max pre-repair row residual | 0 | **0** |
| sub-gate repair moves | 0 | **0** |
| unresolved rows | 0 | **0** |
| designed-vs-realised \|z\| | 0.00898 | 0.01966 |
| `clip_frac` (\|z\| missed by > 0.5) | 0.0035 | 0.0098 |
| sub-gate reservoir, counts/cell | 85.0 | 95.8 |
| sub-gate spread (a15's tail law) | 1.0 | 1.0 |
| per-cell total `T` | 20,001 | 22,322 |

The greedy placement balanced **every row of every one of the 900 blocks on its own**; the
sub-gate slack was never spent. `sub_gate_spread_max = 1.0` on 900/900 confirms every sub-gate
gene ran a15's own tail law.

---

## 4. Falsifiers 2, 3, 4 — measured on the challenge contexts

### 4a. Falsifier 3 — `n_pred`

`probe191.py` recomputes the scorer's own Wilcoxon against each context's real control cells
from `data/challenge/context_{A,B,C}.h5ad` — which, per INBOX-3, *is* the comparison group on
both sides, so this is exact and not a proxy. 40 targets per context, the same 120 blocks for
both packages, self-gene excluded (BRIEF §4–7).

| | ctx A | ctx B | ctx C | **pooled median** |
|---|---:|---:|---:|---:|
| `a15_local` | 4,642 | 4,470 | 4,354 | **4,474** |
| `rc_hybrid` | 4,651 | 4,480 | 4,354 | **4,480** |
| Δ | +9 | +10 | 0 | **+0.13%** |

**Two independent known-answer rows fall out of this.** `a15_local`'s pooled 40-target median
is **4,474**, which is exactly 179's independent probe of a23 (4,474). And the `rc_hybrid`
build's own pass-1 diagnostic — a15's call count on **all 900** blocks, raw BH — reads median
**4,436**, which is exactly 185 §6a / 096's recorded a15 challenge `n_pred` of **4,436**.
185 predicted 4,439 vs 4,436 for this arm; the realised move is +6 on the probe's subsample.

### 4b. Falsifier 4 — the delivered order, and the finding the mandate did not anticipate

Median over probed targets of Spearman(−log₁₀ p, x) **inside our own call set**:

| arm | ctx | `\|H\|` in call set | `\|H\|` over the universe | expression in call set |
|---|---|---:|---:|---:|
| `a15_local` | A | −0.2244 | −0.0235 | +0.2673 |
| `a15_local` | B | −0.0184 | +0.0960 | +0.2441 |
| `a15_local` | C | −0.2139 | −0.0193 | +0.3268 |
| `a15_local` | **mean** | **−0.1522** | — | **+0.2794** |
| **`rc_hybrid`** | A | **+0.9761** | +0.5260 | +0.0483 |
| **`rc_hybrid`** | B | **+0.9713** | +0.5966 | +0.0342 |
| **`rc_hybrid`** | C | **+0.9755** | +0.5589 | +0.0475 |
| **`rc_hybrid`** | **mean** | **+0.9743** | — | **+0.0433** |

**PASS at +0.974 against the +0.95 bar. The order delivers on the real contexts, and the
expression defect is removed there** (+0.279 → +0.043, the same collapse 185 measured on
h1-test at +0.486 → +0.007). The universe value of ≈ +0.55 is by design: `|H|` orders inside
the call block and inside its complement, with the block boundary set by a15's own p.

**But `a15_local` reads −0.152, not the ≈ −0.39 the mandate expected.** That number was
h1-test's; no challenge measurement of it existed. This is a new measurement and it is a
partial falsification of the assumption that h1-test's incumbent ordering defect is
representative of the challenge's. Consequences, stated and not resolved here:

* the **span** the arm delivers is 1.126 Spearman units on the challenge against 1.375 on
  h1-test — **82%**;
* part of the h1-test → board `reach` shrink (164's 0.4277, the a8/a15 pair's 0.390) may
  therefore be *this* — a smaller defect to fix — rather than a transfer failure;
* it does not change the arm's sign, and it does not weaken falsifier 4: the ordering is
  delivered at 97% of a perfect rank correlation.

**Corroboration that the reordering arrived** — fraction of called genes UP by our own
confidence decile, averaged over A/B/C (1 = most confident):

| arm | d1 | d2 | d3 | d4 | d5 | d6 | d7 | d8 | d9 | d10 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `a15_local` | 0.515 | 0.718 | 0.733 | 0.614 | 0.480 | 0.388 | 0.330 | 0.280 | 0.255 | 0.246 |
| `rc_hybrid` | **0.757** | 0.754 | 0.638 | 0.482 | 0.312 | 0.259 | 0.322 | 0.355 | 0.367 | 0.302 |

a15's most-confident decile sits at chance (0.515); `rc_hybrid`'s carries the signal (0.757),
and the top-to-bottom span goes 0.269 → 0.455. This is the same shape 179 measured for the
resample arm (0.774 / 0.739 / 0.728), reached here **without changing the call set**: the
fraction of the call set that is UP is unmoved (60.6/63.4/63.7% → 60.4/63.2/63.7%).

### 4c. Falsifier 2 — `nmae`

`nmae` needs ground truth, which the challenge does not have; the reference cancels out of the
*difference* between two predictions, which it does have. Both were measured.

**Direct, on h1-test, with truth:** a15 **0.96659** → `rc_hybrid` **0.96670** =
**+0.00011 raw**. 185 predicted +0.0002 against a15's twin mean. **Inside the 0.0010 bound by
9×.**

**On the challenge, reference-free**, over all 900 blocks and every gene past the scorer's
5-CPM gate:

```
|d lfc| : max 0.060072   mean 0.001341   median of block means 0.001307
signed  : median -0.000296   range [-0.000891, -0.000128]
```

The same statistic on h1-test is **0.002375**, so the challenge driver is **0.564×** h1-test's,
and the transferred prediction is **+0.00006 raw**.

*Stated because it is the weaker half:* the triangle-inequality bound
`Σ_{S_p}|Δlfĉ| / Σ_{S_p}|lfc|` reads 0.0043 on h1-test (max 0.0096) where the realised Δ`nmae`
is +0.00011 — the bound overstates by 40×, because the two arms' lfc errors move together
relative to the truth. Transferred at 0.564× it reads 0.0024 on the challenge. **The bound is
not the instrument; the realised measurement is.** Both are reported.

### 4d. `mse`'s dispersion half — measured on the shipped artifacts, not on a reconstruction

The digest identity pins `mse`'s profile half and says nothing about `C_hat`, the jackknife of
our own emitted noise, which this emitter changes on purpose. `claim191.py` measures it with
the library's own `jackknife_correction` kernel (the one the builder's knee path already
imports — no `cell-eval2` process), against 096's `C_p` proxy from each context's own 46 NTC
constructs. 16 targets per context.

| ctx | `C_p` | cap `K·C_p` | a15_local `Ĉ/C_p` | rows at cap | rc_hybrid `Ĉ/C_p` | rows at cap |
|---|---:|---:|---:|---:|---:|---:|
| A | 14.4398 | 14.4398 | 1.192 | 16/16 | **2.648** | 16/16 |
| B | 14.8057 | 14.8057 | 0.999 | 6/16 | **2.526** | 16/16 |
| C | 16.4397 | 16.4397 | 1.093 | 16/16 | **3.184** | 16/16 |
| **mean** | | | **1.094** | 38/48 | **2.786** | **48/48** |

Known-answer rows: 185 §6a measured a15 at `Ĉ/C_p` **1.085** on a reconstruction; the shipped
artifact reads **1.094** (0.8% apart). 185 bracketed the rc family at **1.9–3.3**; `rc_hybrid`
reads **2.786**, inside it. 096 quoted `C_p` = 14.544 / 14.772 / 16.372 from 16 constructs;
all 46 give 14.440 / 14.806 / 16.440.

**Reading.** `rc_hybrid` is at the cap on 48/48 rows; `a15_local` is at it on 38/48 and a hair
below on ten context-B rows, so `claim` rises by **+0.15%** overall. That is the same
configuration 185 §6a priced (a15 `claim` 0.24581, rc arms at the 0.24653 cap) and scored at
Δ raw `mse` = **0.88228 → 0.88228, i.e. 0.000000**. So `mse` is pinned by 185's measured route
rather than by an identity, and the pin is now measured on the artifact that would be
submitted.

---

## 5. Falsifier 5 — the rule-10 diff, VERBATIM

### 5a. Effects + operator audit, `rc_hybrid` vs `a15_local`

```
effects.baseline_pre_sha256                             7678bab9ea9d1b8d951b0b24edebde4939  vs 7678bab9ea9d1b8d951b0b24edebde4939  EXACT
effects.candidate_pre_sha256                            0601336fdf1ec5c533dcad69f6cb6c09d6  vs 0601336fdf1ec5c533dcad69f6cb6c09d6  EXACT
effects.cd4.fraction                                    0.09472706196346513                 vs 0.09472706196346513                 EXACT
effects.cd4.post_add_preserve.baseline_output_median_norm  4.122892101713487                vs 4.122892101713487                   EXACT
effects.cd4.post_add_preserve.candidate_pre_preserve_median_norm 4.1665533562021135         vs 4.1665533562021135                  EXACT
effects.final_sha256                                    22c09f22c749bda82d12a0862bacf25b55  vs 22c09f22c749bda82d12a0862bacf25b55  EXACT
effects.post_whiten_sha256                              b2426e15d76d9fe59f794a59bf5198b535  vs b2426e15d76d9fe59f794a59bf5198b535  EXACT
effects.whitening.after_all_median_norm                 4.229295604454299                   vs 4.229295604454299                   EXACT
effects.whitening.after_operator_median_norm            25.73792992728808                   vs 25.73792992728808                   EXACT
effects.whitening.before_all_median_norm                4.238101079325155                   vs 4.238101079325155                   EXACT
effects.whitening.before_operator_median_norm           4.122892101713494                   vs 4.122892101713494                   EXACT
effects.whitening.operator_preserve_scale               0.16018740098217019                 vs 0.16018740098217019                 EXACT
operator.sha256                                         e2f65ecf69b7e0b0                    vs e2f65ecf69b7e0b0                    EXACT
operator.columns_sha256                                 b2cae0d381c3cd70                    vs b2cae0d381c3cd70                    EXACT

fields compared: 78   differences: 0        PASS (rule 10 clean)
```

Operator by path and sha, never by description (rule 11):
`data/donors/whitener_union_pruned.npz`, sha256
`e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123`, 17,942 columns,
`basis_prune: essential`, `all_panel_coordinates_excluded: true`. Identical in both packages
and identical to a15's recorded value.

### 5b. The config diff — exactly the five pre-declared keys

The expected set was written into `metadiff191.py` **before** the diff was read.

```
emitter                    rc_hybrid 'rc_hybrid'            a15_local <absent>
rc_jitter_seed_base        rc_hybrid 1000                   a15_local <absent>
rc_key                     rc_hybrid 'gated_effect'         a15_local <absent>
rc_scorer_min_cpm          rc_hybrid 5.0                    a15_local <absent>
rc_sign_policy             rc_hybrid 'abs'                  a15_local <absent>

config keys that moved: ['emitter', 'rc_jitter_seed_base', 'rc_key', 'rc_scorer_min_cpm', 'rc_sign_policy']
pre-declared expected : ['emitter', 'rc_jitter_seed_base', 'rc_key', 'rc_scorer_min_cpm', 'rc_sign_policy']
PASS (5 differences)
```

`alpha` 0.982, `min_cpm` 10.0, `keep_top` 8000, `keep_spread` 0.25, `tail_spread` 1.0,
`overdispersion` 4.0, `whitener_id` union_pruned, `cd4_target_f` 0.1, `spread_map` block,
`effect_scale` 1.0, all `shrink_*` null, `mass_fix` absent, donor weights all 1.0, seeds
0 / 0 / 1000003 — **all unchanged**. `resample_keep_top` is absent from both: the builder
refuses it under `rc_hybrid`, so two emission channels cannot be live at once.
`hashes.inputs` and `hashes.implementation` are **identical dictionaries** between the two
packages.

### 5c. Everything else in the meta that moved, stated rather than omitted

| key | `rc_hybrid` | `a15_local` | why |
|---|---|---|---|
| `output.sha256` / `output.path` | `e0b0677f…` | `26fb3834…` | the arm |
| `output.bytes` | 13,202,460,120 | 19,049,370,200 | −30.7%; the emitter concentrates counts |
| `output.nnz` | 1,647,149,337 | 2,377,814,190 | 34.68% vs 50.06% of the density cap |
| `output.contexts.{A,B,C}.stored_entries` | 567.9 / 529.0 / 550.2 M | 816.5 / 771.2 / 790.4 M | the arm |
| `output.rc_hybrid.*` (42 keys) | present | absent | the new per-block diagnostics block, deliberately outside `config` |
| `peak_rss_gib` | 19.834 | 14.808 | allocator / page-cache state, not a code change — see §7 |

Nothing else differs anywhere in either meta.

---

## 6. `a15_local` against a15's RECORDED audit — the arm64/x86_64 drift row

| scalar | `a15_local` (arm64) | a15 recorded (x86_64) | rel drift |
|---|---:|---:|---:|
| `effects.whitening.after_all_median_norm` | 4.2292956044543 | 4.2292956044543 | 4.2e−16 |
| `effects.whitening.before_all_median_norm` | 4.23810107932515 | 4.23810107932515 | 0.0 |
| `effects.whitening.after_operator_median_norm` | 25.7379299272881 | 25.7379299272881 | 2.8e−16 |
| `effects.whitening.before_operator_median_norm` | 4.12289210171349 | 4.12289210171349 | 0.0 |
| `effects.whitening.operator_preserve_scale` | 0.16018740098217 | 0.16018740098217 | 3.5e−16 |
| `effects.cd4.fraction` | 0.0947270619634651 | 0.0947270619634652 | 4.4e−16 |
| `effects.cd4.residual_frame.effective_f_at_cd4_fraction` | 0.0527832268452328 | 0.0527832268452327 | 3.9e−16 |
| `effects.cd4.residual_frame.f_per_unit_fraction` | 1.05566453690466 | 1.05566453690465 | 4.2e−16 |
| `effects.cd4.residual_frame.solved_fraction` | 0.0947270619634651 | 0.0947270619634652 | 4.4e−16 |
| `effects.cd4.residual_frame.whitened_base_median_row_norm` | 4.2292956044543 | 4.2292956044543 | 4.2e−16 |
| `effects.cd4.residual_frame.residual_median_row_norm_at_unit_fraction` | 4.46471738570914 | 4.46471738570914 | 0.0 |
| `effects.cd4.post_add_preserve.baseline_output_median_norm` | 4.12289210171349 | 4.12289210171349 | 6.5e−16 |
| `effects.cd4.post_add_preserve.candidate_pre_preserve_median_norm` | 4.16655335620211 | 4.16655335620212 | 4.3e−16 |
| `effects.cd4.post_add_preserve.post_add_global_preserve_scale` | 0.989521013951823 | 0.989521013951823 | 2.2e−16 |
| **worst relative drift, 14 scalars** | | | **6.5e−16** |

| digest | `a15_local` | a15 recorded | |
|---|---|---|---|
| `effects.baseline_pre_sha256` | `7678bab9ea9d1b8d` | `7678bab9ea9d1b8d` | **EXACT** |
| `effects.candidate_pre_sha256` | `0601336fdf1ec5c5` | `0601336fdf1ec5c5` | **EXACT** |
| `effects.post_whiten_sha256` | `b2426e15d76d9fe5` | `51a640a84c2eec4b` | DIFFERS |
| `effects.final_sha256` | `22c09f22c749bda8` | `2980d70655f2a9d4` | DIFFERS |
| `operator.sha256` | `e2f65ecf69b7e0b0` | `e2f65ecf69b7e0b0` | **EXACT** |
| **`output.sha256`** | **`26fb38341c007661…`** | **`26fb38341c007661…`** | **EXACT (19,049,370,200 bytes)** |

**Two things this row establishes that 170 §1 could not.** First, the drift is *localised*: the
two PRE-whitening digests are exact and the two POST-whitening digests differ, so the whole
difference enters at `whitener.apply` — exactly the pair of 17,942² GEMMs 170 identified, and
nothing else. Second, and more usefully, **the drift does not reach the artifact at all**: the
delivered counts are `floor(n·L·profile)` and a 6.5e−16 relative perturbation of the profile
cannot move an integer whose scale is 10⁴, so the emitted `.h5ad` is byte-identical across the
two architectures. Any future build of a15's configuration on this box can be checked against
a15's recorded `output.sha256` directly.

---

## 7. What was built into the builder, and the byte-identity proof

**One new production module and one new flag. `src/vcc2026/donor_transfer.py` is untouched**
(its frozen sha `9ffca7ec…` is unchanged in both metas), which is the same containment 179
chose for `--emitter resample_mix`.

* **`src/vcc2026/rc_emit.py`** (873 lines) is **generated**, not retyped, by
  `experiments/191-rc-hybrid-build/assemble191.py`, which slices the two implementations that
  produced 185's numbers out of their source files verbatim and appends the driver:
  `experiments/092-call-set-truncation/de_kernel.py` (sha256 `48b0fbe5…`) →
  `ctrl_tiesums` / `wilcoxon_z`; `experiments/185-rank-controlled-emission/eme185.py`
  (sha256 `decb6917…`) → the whole EME; `rc_driver.py.in` (sha256 `36a2a1a1…`) → the only new
  code. The assembler asserts that the two files' binary searches agree behaviourally before
  taking one copy of them. The mandate said "port it, do not re-derive it"; this is the form of
  that instruction which a reader can check.
* **`--emitter rc_hybrid`** in `scripts/make_submission_pm_cd4.py`. It fails closed on every
  knob that would make the rendered effect something other than `where(gate, effect, 0)` —
  `--effect-scale ≠ 1`, either `--shrink-*` key, `VCC_MASS_FIX`, `VCC_RESAMPLE_X`,
  `--resample-keep-top` — because an ordering emitter pointed at the wrong key vector would
  look exactly as reasonable as one pointed at the right one.

**The byte-identity is proven, not asserted.** `test_rc_hybrid_is_off_by_default_and_the_default_path_is_byte_identical`
resolves the newest committed builder revision that predates `rc_hybrid` **out of git**, loads
it, drives the **real writer** on a synthetic three-context bundle with both modules, and
compares sha256 of the emitted `.h5ad` — at the bare default **and** at a15's shipped emission
policy (`tail 1.0 / kt 8000 / ks 0.25 / φ 4`). A same-module comparison cannot answer this
question, because the patch is on both sides of it. The test **fails**, never skips, if the
revision cannot be found, and it checks non-vacuity by requiring `resample_mix` on the same
inputs to produce a *different* artifact.

Two more tests were added:
`test_rc_hybrid_delivers_the_quasi_pseudobulk_and_welds_every_cell` (6/6 blocks of a real
three-context build: identical panel totals, different file, every cell welded, all four
placement diagnostics zero) and
`test_rc_hybrid_fails_closed_on_every_knob_that_would_move_its_key`.

**`uv run --with pytest python -m pytest scripts/test_make_submission_pm_cd4.py -q` →
60 passed** (57 pre-existing + 3 new).

---

## 8. What went over budget, and what is NOT measured

**Peak RSS exceeded the mandate's 14 GB: 14.81 GiB (`a15_local`) and 19.83 GiB (`rc_hybrid`).**
Reported rather than smoothed over. The attribution is unambiguous and it is not the arm:

* the `rc_hybrid` run's high-water was already **19.742 GiB at `before challenge emission`**,
  i.e. it was set in the frozen effects stage before a single line of `rc_emit` had run;
* the rc control block added **0.092 GiB** on top of it (19.742 → 19.834), and nothing after;
* the *same* effects stage, unchanged code, measured **14.81 GiB** in the `a15_local` run
  minutes earlier, **17.47 GiB** for a23 and **13.74 GiB** for `hy_kt8000` (179) on this
  machine. It is allocator / page-cache state on a 5.8 GB donor-table load, not a code change.

Lowering it means changing frozen model code, which this mandate forbids. The two builds were
run **serially**, never concurrently, and the box has 64 GB.

| not measured | why | can it move? |
|---|---|---|
| board `reach` / `fid` / `jac` for either package | rule 9: no harness we hold reads them; and the board is the coordinator's slot, not mine | this is the whole point of the slot |
| challenge-side `nmae`, `reach`, `fid`, `jac` absolutely | no ground truth for A/B/C exists | — |
| `cell-eval2` on either package | forbidden by the mandate; not needed for any falsifier | — |
| a re-measured `claim` for a15 at 096's exact 16-construct `C_p` convention | all 46 constructs were used instead, which is the builder's own convention; the two agree to 0.7% | no |
| the `signed` sign policy variant | the preregistration binds this arm to `abs`; 185 §7's design note records `signed` as the immune variant at −0.007 raw `reach` | not this arm |

**Two things worth flagging to the coordinator, neither a falsifier:**

1. `sign_flip_frac` — the share of designed genes whose Wilcoxon `z` sign opposes `sign(H)`
   under the `abs` policy — is **0.439 median on the challenge** against 0.370 on h1-test.
   `cell-eval2` reads direction from the lfc sign and never from the `z` sign (185 §2c,
   `cell_eval2/metrics/direction.py`), so it costs nothing under the current scorer; it is the
   exposure 185 §7's design note names.
2. The 31% nnz drop is the visible signature of the arm. It is inside every format constraint
   and `vcc prep` passes, but a package one third sparser than every previous one is a thing a
   reviewer should be told about before, not after.

---

## 9. Verdict

# BUILDABLE.

`submissions/191-rc_hybrid.h5ad`, sha256
`e0b0677f984a94a1c4a7af2add196908cdab90d7acb463eb1d6e0397c2d870a8`. All six pre-committed
falsifiers PASS, all fourteen of 185's h1-test known-answer numbers reproduce exactly through
the production emitter, and the reference package `submissions/191-a15_local.h5ad` is
byte-identical to a15 itself. The package is prepped-clean and **not submitted**; no model
name and no description were chosen — `AGENTS.md` §14 and the slot are the coordinator's.

---

## 10. Files

| file | what it produces |
|---|---|
| `assemble191.py` | generates `src/vcc2026/rc_emit.py` from `de_kernel.py` + `eme185.py` + `rc_driver.py.in`, verbatim |
| `rc_driver.py.in` | 191's only new emitter code: the control block, the two-pass driver, the CSR conversion |
| `build191.sh` | a15's configuration, one variable; `build191.sh {a15_local,rc_hybrid}` |
| `ka191.py` / `work/ka191.log`, `work/de191.*.npz` | gate 0 — the production emitter over h1-test, 3 arms × 85 perturbations |
| `score191.py` / `work/score191.json` | gate 0 — 185's row read off those tables with 185's own instruments |
| `digest191.py` / `work/digest191.{log,json}` | falsifiers 1 and 2 — 900/900 digests, the weld, `Δlfc` |
| `probe191.py` / `work/probe191.{log,json}` | falsifiers 3 and 4 — `n_pred` and the delivered order on A/B/C |
| `metadiff191.py` / `work/metadiff191.log` | falsifier 5 and the arm64/x86_64 drift row |
| `claim191.py` / `work/claim191.{log,json}` | `mse`'s dispersion half on the shipped artifacts |
| `work/vccprep191.log` | falsifier 6 |
| `work/build-{a15_local,rc_hybrid}.log` | the two builds |
