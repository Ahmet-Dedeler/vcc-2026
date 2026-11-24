# 202 — `rc_callset` BUILD RESULT: **BUILDABLE**

192 §5.3's rider on `rc_hybrid`, built to the preregistration at `PREREGISTRATION.md`, one
variable against the incumbent: the **call set**. Two known-answer rows and all six falsifiers
**PASS**. The package is built, prepped-clean and **NOT submitted**; no model name and no
description were chosen (AGENTS.md §14 — that is the coordinator's call).

```
call set per (context, perturbation) = top 400 by |H|  U  top 2,200 by control CPM   (2,600)
order                                = 191's rc_hybrid law, unchanged
everything else                      = rc_hybrid, unchanged
```

| | |
| --- | --- |
| package | `submissions/202-rc_callset.h5ad` |
| sha256 | `2a894aa543db7f574ebc024dbb3ebd5e36f923827c4e9a01f263092fe49692e0` |
| bytes | 12,648,807,800 |
| meta | `submissions/202-rc_callset.h5ad.meta.json` |
| peak RSS | **16.95 GiB** (mandate ≤ 20 GB), 1,399 s wall, built serially |
| incumbent rebuilt alongside | `submissions/202-rc_hybrid.h5ad`, 13,202,460,120 bytes, peak 17.43 GiB |

---

## 1. Known-answer rows (rule 10's corollary)

### Row 1 — the incumbent must survive the patch, bit for bit

`build202.sh rc_hybrid` reruns 191's exact invocation through the **patched** builder and
emitter. Verified by hashing the two files on disk, not by reading the build log:

```
e0b0677f984a94a1c4a7af2add196908cdab90d7acb463eb1d6e0397c2d870a8  submissions/202-rc_hybrid.h5ad
e0b0677f984a94a1c4a7af2add196908cdab90d7acb463eb1d6e0397c2d870a8  submissions/191-rc_hybrid.h5ad
```

Both 13,202,460,120 bytes; the required digest is the mandate's, character for character.
**PASS.** The change cannot have perturbed the incumbent at any scale on any of 360,000 cells.

**Row 1b, the same claim on h1-test at table precision** (`work/ka202_emit.log`): `rc_hybrid`
driven through the patched emitter reproduces 191's DE table **bit-identically on all
85 × 10,779 entries** of `pred_p`, `pred_padj`, `pred_lfc` and `pred_z` (max |diff| **0** on each),
and its per-gene panel-total digest matches 191's on **85/85**.

**Row 1c, proven against git rather than asserted** — the mandate's "prove it against the current
builder resolved out of git". Three new tests were added to `scripts/test_make_submission_pm_cd4.py`;
the first,
`test_rc_callset_is_off_by_default_and_every_existing_path_is_byte_identical`, resolves the
newest committed revision of **both** `scripts/make_submission_pm_cd4.py` and
`src/vcc2026/rc_emit.py` that predates `rc_callset` — both land on **`3da402d`** — loads them,
drives the **real writer** on a synthetic three-context bundle, and compares sha256 of the emitted
`.h5ad` at the bare default **and** at a15's shipped policy (`tail 1.0 / kt 8000 / ks 0.25 / φ 4`),
on the default, `rc_hybrid` and `rc_hist_a` paths. It fails, never skips, if the revision cannot be
found, and checks non-vacuity by requiring `rc_callset` on the same inputs to produce a *different*
artifact.

```
uv run --with pytest python -m pytest scripts/test_make_submission_pm_cd4.py -q
65 passed in 62.86s          (62 pre-existing + 3 new)
```

**Provenance.** `src/vcc2026/rc_emit.py` is generated, not hand-edited. Re-running
`assemble202.py` reproduces the shipped module byte for byte —
sha256 `6f5e360b729558aa7caeb702c79218948997e4c6e584c457a2602f0df0cb4456` — from the recorded
sources `de_kernel.py 48b0fbe5…`, `eme185.py decb6917…`, `rc_driver.py.in 975fff7b…`.

### Row 2 — 192's own h1-test row must reproduce, or the call set is not 192's

`work/ka202_prescribed.log`, the prescribed set applied to 191's frozen h1-test tables at
`K` = 2,600, `M` = 400, `f` = 0.50, marginal to `rc_hybrid` at `n_pred` 4,474:

| member | base | arm | Δ measured here | 192 §5.3 | |
| --- | ---: | ---: | ---: | ---: | --- |
| `fid` | 0.50771 | 0.50868 | **+0.00097** | +0.00097 | PASS |
| `jac` | 0.04058 | 0.05346 | **+0.01287** | +0.01288 | PASS |
| `reach` | 0.26126 | 0.27249 | **+0.01123** | +0.01123 | PASS |

Every base and every arm value reproduces 192's to all five printed decimals. The mandate quotes
the row as `fid` +0.0010 / `jac` +0.0129 / `reach` +0.0112; these are those numbers.
**Row 2: PASS** — the call set built here is 192's call set and not a near neighbour of it.

### Row 3 (unplanned, and the strongest single check on the challenge pipeline)

The falsifier probe re-measures **`rc_hybrid` on the real contexts** in the same pass as the arm.
It reproduces 191's recorded numbers exactly: `n_pred` A **4,651** / B **4,480** / C **4,354**,
pooled **4,480**; ρ(−log p, |H|) A **+0.9761** / B **+0.9713** / C **+0.9755**, mean **+0.9743** —
the mandate's quoted baseline, to four decimals. So the instrument reading the arm is the
instrument that produced the baseline.

*(191 ran its probe at 40 targets per context. A first run here at 30 targets drew a different
subset and read `rc_hybrid` at +0.9754 / 4,466; the 40-target run below is the comparable one and
both logs are kept, `work/probe202.log` and `work/probe202-perts30.log`.)*

---

## 2. The six falsifiers

| # | falsifier | measured | verdict |
| --- | --- | --- | --- |
| 1 | panel-total digests identical on 900/900; structure on 0/900 | **900/900** and **0/900** | **PASS** |
| 2 | challenge `n_pred` within ±3% of 2,600 | **2,599** pooled (**−0.04%**) | **PASS** |
| 3 | the two halves of the union do different jobs | **117 of 400** shared | **PASS** |
| 4 | delivered order preserved vs `rc_hybrid`'s +0.9743 | **+0.9993** (rose +0.0250) | **PASS** |
| 5 | rule-10 config diff is exactly the pre-declared one | **4 keys, 4 expected** | **PASS** |
| 6 | `vcc prep --dry-run` clean; nnz against the cap | clean; **33.22%** | **PASS** |

### Falsifier 1 — the profile is pinned by identity

`digest191.py`, reused verbatim (sha `f346c36c…`, byte-identical to 191's copy), streaming over
`202-rc_callset.h5ad` against `191-rc_hybrid.h5ad`. *The tool prints its two positional arguments
as `rc_hybrid` and `a15_local`; here those slots hold `rc_callset` and `rc_hybrid`.*

```
per-gene PANEL TOTAL digest identical on   900/900
panel total SUM identical on               900/900
nonzero-STRUCTURE digest identical on      0/900
```

The individual digests are also the same *values* 191 recorded against a15
(`(A, ABCD1) bb6befc01caf7493`, `(C, ZRANB3) 3c9d8f92f3fa4d93`, …), so the chain
**a15 = `rc_hybrid` = `rc_callset`** holds on the delivered pseudobulk, and `rc_hybrid`'s board
point already proved the first link reproduces `pds` and `mse` to sixteen significant figures.
0/900 on the nonzero structure makes the identity non-vacuous: these are genuinely different files.

Two riders from the same pass:

* **the weld survives** — per-cell library totals min 17,695, max 22,323, sd 0.394, identical to
  `rc_hybrid`'s row; every cell welded, and the 1,000,000 format cap is 45× away.
* **`nmae` is pinned harder than `rc_hybrid`'s was** — reference-free |Δ log2FC| over all 900
  blocks and the whole DE universe: **max 8.2e−05, mean 2.0e−06**. The a15 → `rc_hybrid` pair
  measured max 0.060072 / mean 0.001341 on the same instrument and the board returned
  Δ`nmae` = −0.00003 raw. This arm is **670× smaller in the mean**.

### Falsifier 2 — challenge `n_pred`

`probe191.py`, reused verbatim (sha `81535b8e…`), 40 targets × 3 contexts, Wilcoxon + BH computed
against **the real control cells the scorer will use**, target row excluded as every member does:

| | A | B | C | pooled |
| --- | ---: | ---: | ---: | ---: |
| `rc_callset` | **2,599** | **2,599** | **2,599** | **2,599** |
| `rc_hybrid` | 4,651 | 4,480 | 4,354 | 4,480 |

**2,599 against a prescribed 2,600 is −0.04%**, inside the ±3% bound by two orders of magnitude,
and −42.0% against `rc_hybrid`. The single missing gene is the **target's own row**, which the
prescribed set contains (the panel was selected for high basal expression, §2.1) and which every
metric excludes; the builder's own diagnostic, which does not exclude it, reports the set delivered
at exactly 2,600.

The builder's diagnostics over **all 900 blocks** (not a 120-block sample) close the question:

```
call_set_prescribed     min 2600   median 2600   max 2600
call_set_delivered      min 2600   median 2600   max 2752
call_set_missed         min    0   median    0   max    0
call_set_unsuppressible min    0   median    0   max  152
call_set_jaccard        min 0.9448 median 1.0000 max 1.0000
call_set_z_scale        min 0.353  median 0.575  max 0.660
```

Nothing prescribed is ever missed; the only excess is the **unsuppressible floor** the emitter
documents — genes whose delivered panel total is far enough from the control that no allocation
over 400 cells can make them look like it. On the challenge that floor is **median 0**, against
h1-test's median 35, because the challenge's delivered effect is gentler. This was the arm's
single real exposure and the challenge did not realise it.

### Falsifier 3 — call-set composition, and whether the halves do different jobs

`compose202.py` (the one new tool: `probe191.py`'s I/O loop with the call set recorded as a mask
instead of a count), 40 targets × 3 contexts = 120 blocks, medians:

```
prescribed  |top 400 by |H|  U  top 2,200 by CPM|      2,600
delivered by BH                                        2,599
  jaccard(delivered, prescribed)                       0.9996
  prescribed but MISSED                                1  (max 1 — the excluded target row)
  delivered but UNSUPPRESSIBLE                         0  (p90 0, max 112)

THE TWO HALVES
  |top 400 by |H||                                       400
  |top 2,200 by CPM, head excluded|                    2,200
  |halves SHARED| (head n unrestricted top-2,200 CPM)    117
  the |H| head's median rank in the CPM order          3,992 of 9,929
  the CPM half's median rank in the |H| order          4,587 of 9,929

AGAINST `rc_hybrid`'s OWN DELIVERED CALL SET
  `rc_hybrid` n_pred                                   4,480
  |delivered n rc_hybrid|                              1,194
  jaccard(delivered, rc_hybrid)                        0.2072
  of the 400-gene |H| head, in rc_hybrid's set           272
  of the 2,200-gene CPM half, in rc_hybrid's set         928
```

**The halves are different objects.** Only **117 of the 400** `|H|`-head genes would have been
picked by control CPM at all — **283 of them, 71%, are reachable only by the head** — and the
typical head gene sits at rank **3,992 of 9,929** in the expression order, i.e. the 40th
percentile, nowhere near the top 2,200. Symmetrically the typical CPM-half gene sits at rank
**4,587 of 9,929** in the `|H|` order. So the union is not a re-spelling of either key, which is
exactly the premise 192 §5.3's `+0.011 reach` half rests on. *(The companion statistic "2,200 of
2,200 CPM-half genes lie beyond the head's `|H|` rank" is reported in the log but is a tautology of
the construction — the head takes the top 400 and the CPM half excludes it — and carries no
information; the median-rank numbers above are the real ones.)*

Against the incumbent the arm is a **genuine set change, not a truncation**: jaccard 0.2072, and
at the median block **1,405 of the 2,599 called genes are ones `rc_hybrid` never called** while
**3,286 of `rc_hybrid`'s calls are dropped**. Both halves contribute new genes — 272 of the 400
head genes and 928 of the 2,200 CPM genes are already in `rc_hybrid`'s set, leaving ~128 and
~1,272 new. *(These are medians of separate per-block quantities, so they sum to 1,400 rather than
exactly to 1,405; the per-block values are in `work/compose202.json`.)*

### Falsifier 4 — the delivered order

ρ(−log p, `|H|`) inside the call set, on the challenge contexts, same pass as the baseline:

| | A | B | C | mean |
| --- | ---: | ---: | ---: | ---: |
| **`rc_callset`** | **+0.9994** | **+0.9996** | **+0.9990** | **+0.9993** |
| `rc_hybrid` | +0.9761 | +0.9713 | +0.9755 | +0.9743 |

**The order did not drop; it rose by +0.0250**, essentially to 1. Preregistration falsifier 5 —
"a real `reach` loss would mean the order was disturbed by the set change" — is not triggered.
The mechanism is falsifier 2's last row: with `unsuppressible` at median 0 the delivered set *is*
the prescribed set, whose interior is `|H|`-monotone by construction, so nothing breaks the
monotone map.

**One place this differs from h1-test, stated rather than buried.** On h1-test the same ρ *fell*,
+0.9837 → +0.9212 (`work/rho202.log`), because that surface carries a median of 35 unsuppressible
genes per perturbation that enter the call set out of order; over the prescribed set alone it read
+0.9719. The challenge does not have that tail. Two further order statistics from the same pass:
over the whole gene universe ρ rises **+0.5605 → +0.7599**, and ρ(−log p, expression) moves
**+0.0433 → −0.2611** — the sign of the intended trade, the CPM half buying membership while the
`|H|` head keeps the front of the walk. The call set's up-fraction moves 62.4% → 49.1%.

### Falsifier 5 — the rule-10 config diff, verbatim

`metadiff202.py`, `EXPECTED_CONFIG_DIFF` pre-declared in source before the diff was read:

```
=== A: effects + operator audit, rc_callset vs rc_hybrid (rule 10) ===
  effects.final_sha256                      22c09f22c749bda8… vs 22c09f22c749bda8…  EXACT
  effects.post_whiten_sha256                b2426e15d76d9fe5… vs b2426e15d76d9fe5…  EXACT
  effects.whitening.after_all_median_norm   4.229295604454299 vs 4.229295604454299  EXACT
  operator.sha256                           e2f65ecf69b7e0b0  vs e2f65ecf69b7e0b0    EXACT
  operator.columns_sha256                   b2cae0d381c3cd70  vs b2cae0d381c3cd70    EXACT
  fields compared: 78   differences: 0
  A PASS (rule 10 clean)

=== B: config diff (rule 10's intended-difference count) ===
  emitter                    rc_callset 'rc_callset'  rc_hybrid 'rc_hybrid'
  rc_callset_cpm             rc_callset 2200          rc_hybrid '<absent>'
  rc_callset_head            rc_callset 400           rc_hybrid '<absent>'
  rc_callset_z_cap           rc_callset 1.9           rc_hybrid '<absent>'
  config keys that moved: ['emitter', 'rc_callset_cpm', 'rc_callset_head', 'rc_callset_z_cap']
  pre-declared expected : ['emitter', 'rc_callset_cpm', 'rc_callset_head', 'rc_callset_z_cap']
  B PASS (4 differences)
```

Exactly the pre-declared call-set difference and nothing else. `meta['operator']['sha256']` is
**identical**. `effects.whitening.after_all_median_norm` has **not drifted** — 4.229295604454299
on both sides, and against a15's *recorded* meta the worst relative drift over the 14 audit
scalars is **6.5e−16**, last-bit arm64/x86_64 noise (170 §1). The four shared `rc_*` keys
(`rc_key`, `rc_sign_policy`, `rc_scorer_min_cpm`, `rc_jitter_seed_base`) did **not** move, which is
what makes this one variable rather than two: the ordering key, the sign policy, the scorer gate
and the tie-breaking seed are the incumbent's.

Section C of the log states every remaining meta difference rather than omitting it; all of them
are the output's own byte/nnz/diagnostic fields, listed in falsifier 6.

### Falsifier 6 — format contract and density

```
vcc prep submissions/202-rc_callset.h5ad -g data/challenge/gene_names.csv \
  --perts data/challenge/pert_counts.csv -o submissions/202-rc_callset.vcc --dry-run

✓ dry run — no .vcc written
  cells:         360000
  genes:         18533  (encoding: float32)
  contexts:      A=120000, B=120000, C=120000
  targets:       verified against the official list
  normalization: counts-preserved
```

| | `rc_callset` | `rc_hybrid` |
| --- | ---: | ---: |
| stored entries A / B / C | 539,067,895 / 512,298,146 / 526,660,526 | 567,940,597 / 528,966,536 / 550,242,204 |
| **nnz** | **1,578,026,567** | 1,647,149,337 |
| **% of the 4.75e9 cap** | **33.22%** | 34.68% |
| nnz per cell | 4,383.4 | 4,575.4 |
| bytes | 12,648,807,800 | 13,202,460,120 |

**Correction to the mandate's parenthetical.** "A smaller call set does not change the emitted
counts" is right about the *per-gene panel totals* — falsifier 1 proves those to the digest — but
the **nnz does move**, −4.2%, because the same totals are *placed across cells* differently once
the target `|z|` of the complement is scaled. It moves the arm further from the cap, so it is
harmless here, but it is a real difference and it is why the file is 554 MB smaller.

---

## 3. What the local surface says the arm will do

Not a falsifier, and pinned to the mandate's instruments rather than offered as a board
prediction — `fid`, `jac` and `reach` are unreadable in absolute terms on every local harness
(rule 9). The **emitted** h1-test arm (`work/score202.log`, real counts through the patched
emitter, not the prescribed-set retrodiction of row 2):

| `f` | `n_pred` | Δ`fid` | Δ`jac` | Δ`reach` | Δ`nmae` | overall @1.290 | @0.4277 |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.35 | 2,634 | +0.00167 | +0.01411 | +0.02105 | +0.00000 | +0.01243 | +0.00902 |
| 0.40 | 2,634 | +0.00167 | +0.01435 | +0.01955 | +0.00000 | +0.01218 | +0.00900 |
| 0.50 | 2,634 | +0.00143 | +0.01420 | +0.00810 | +0.00000 | +0.00920 | +0.00789 |
| 0.65 | 2,634 | −0.00152 | +0.01227 | +0.00125 | +0.00000 | +0.00503 | +0.00482 |

Positive at every `n_real` level and at both `reach` transfer factors, as 192 reported. The
emitted arm reads `fid` **1.48×** and `jac` **1.10×** its prescribed prediction and `reach`
**0.72×** — i.e. the emission realises rather more of the `jac` half and rather less of the
`reach` half than the frozen-table retrodiction implied, which is worth knowing because ~88% of
the predicted board gain is `jac`. Δ`nmae` is +0.00000 at every level, consistent with falsifier 1.

**A mechanism suite** (`work/mech202.log`, 6 groups, all PASS) checks on a synthetic block that
the two laws draw the same RNG stream in the same order (88,636 calls), that the totals identity
is not vacuous, that the prescribed set is delivered when recomputed independently of the emitter,
that the complement's order is `rc_hybrid`'s to within 0.02, that `z_cap` is the only place the
p-value multiset moves, and that a prescribed set with `call_set_first=False` is refused.

---

## 4. Verdict

**BUILDABLE.** Two known-answer rows and all six falsifiers PASS, on the measured numbers above.
The incumbent is reproduced to the byte through the patched builder; the delivered profile is
digest-identical to it on 900/900 blocks while the files differ on 0/900 structure digests; the
challenge call set lands at 2,599 against a prescribed 2,600; the delivered order rose rather than
fell; and the config diff is exactly the four pre-declared keys with the operator sha and the
whitening norm untouched.

**Not submitted.** `vcc submit` was never run. No model name and no description were chosen.
Package: `submissions/202-rc_callset.h5ad`, sha256
`2a894aa543db7f574ebc024dbb3ebd5e36f923827c4e9a01f263092fe49692e0`. The `.vcc` was dry-run only
and not written.

### Files

| path | what |
| --- | --- |
| `work/build-rc_callset.log` | the arm's build, meta and `/usr/bin/time -l` |
| `work/build-rc_hybrid.log` | known-answer row 1's rebuild of the incumbent |
| `work/digest202.json` | falsifier 1's 900-block digests |
| `work/probe202.log`, `work/probe202.json` | falsifiers 2 and 4, 40 targets × 3 contexts |
| `work/probe202-perts30.log`, `.json` | the earlier 30-target draw, kept |
| `work/compose202.log`, `work/compose202.json` | falsifier 3 |
| `work/metadiff202.log` | falsifier 5, all four sections |
| `work/pytest202.log` | 65 passed |
| `work/ka202_prescribed.log` | known-answer row 2 and 192's full (K, M) grid |
| `work/ka202_emit.log` | known-answer row 1b and the emitted h1-test tables |
| `work/score202.log`, `.json` | §3's emitted h1-test arm |
| `work/mech202.log`, `work/rho202.log`, `work/lfc202.log` | mechanism, order, fold-change |
| `probe191.py`, `digest191.py` | verbatim copies of 191's tools, so their output lands here |
| `compose202.py` | the one new tool (falsifier 3) |
