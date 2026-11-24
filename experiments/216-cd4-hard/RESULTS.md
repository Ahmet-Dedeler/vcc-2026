# 216 — The CD4 dose cannot be what the leaders did. **The published atlas is byte-for-byte the object we already ship, nothing landed in the doubling window, and the dose lever's ceiling is +32.4% cosine on PLOLO and +46.0% projected onto the board, against the leaders' +106%.** Re-priced on `cd4f034`'s realised board rates, **`f = 0.34` is the argmax and no dose clears the gate** — and the per-context spread is 3.2%, so the D/E/F row stays a constant.

**Date:** 2026-09-04 · four PLOLO folds (162's instrument, imported verbatim) · **no `vcc`, no
`cell-eval2`, no build, no submission, no production change, nothing posted anywhere, and no
query carried anything about our method** (`AGENTS.md` §14).

**Peak RSS: 11.70 GB** (`op216.py`, the eight 17,942² eigendecompositions), then **10.71 GB**
(`gate216.py`, which held all four fold operators at once). Both are inside the mandate's
12 GB ceiling and both are **above the standing ~10 GB house guidance, flagged rather than
buried**; everything else peaked at **3.3 / 4.99 / 6.18 / 6.32 / 8.55 GB**. **Three memory
excursions happened and all three are on the record rather than buried:** a first `op216.py`
reached **18.5 GB** and was killed within minutes, because `np.linalg.eigh` dispatches to
LAPACK `dsyevd`, whose workspace is ~2n² doubles on top of the input and the eigenvector
output (§5.4 records the three fixes — in-place spectral reconstruction, an F-ordered input,
and `scipy`'s `overwrite_a`); the discarded first `opsix216.py` run reached **16.2 GB** for the
same reason and is kept as `work/opsix216-BUGGY.log`; and a duplicated background chain briefly
ran two copies of one job, killed before either wrote a result. **The kept runs all sit at or
below 11.70 GB.**

**Downloaded to disk: 32.2 MB** (8 files) plus **30.8 MB of HTTP range reads never written to
disk**, against a 60 GB budget — the inventory in §3.1 was read before a byte of bulk data was
fetched, and none of it was needed.

Everything written is inside `experiments/216-cd4-hard/`, except those 32.2 MB of public
supplementary tables in `data/external/cd4_cell2026/`, which the mandate names as the one
exception.

---

## 0. The answer in fourteen lines

1. **Phase 1 settles the date coincidence: it is a coincidence.** The object our CD4 banks
   were built from — `GWCD4i.DE_stats.h5ad`, 16,786,240,107 B — is **byte-identical on the
   live bucket at 12 probed ranges spanning offset 0 to 16.79 GB**, `LastModified`
   **2026-05-28T23:43:31Z**. **Zero objects on that bucket have a modification date inside
   2026-08-26 → 09-03.** §3.
2. **The one artefact that did change (2026-08-24, four days before publication) added no
   data.** `DE_stats.suppl_table.csv` gained 14 columns — all 14 were **already `.obs`
   columns of the h5ad we hold**, same 33,983 rows, same perturbation set. §3.2.
3. **Panel coverage is already at the published ceiling.** Ours: **291 union / 276
   all-three-conditions**. The published DE table: **291 / 276**. Identical. The nine
   never-covered panel genes have **zero rows** in the published table. §3.3.
4. **The one materially different processed object is the pseudobulk** —
   `GWCD4i.pseudobulk_merged.h5ad`, **17,771 official genes against our 10,113**, a strict
   superset, public since **2026-05-28**. Its extra 7,658 genes hold **11.94%** of the gated
   truth energy we cannot reach. **Measured, not assumed: cutting that tail off the residual
   makes the cosine BETTER, not worse** — the elasticity is **−0.55 to −0.97**, so the
   pseudobulk route is worth **+3.2% cosine under the most generous law available and
   negative under the measured one**. §3.4–3.5.
5. **The dose curve turns over at `f* = 1.018`** (interior maximum, quadratic in log f, all
   four folds agree in sign), at **+32.4% cos_agg over the shipped `f = 0.10`** and **+47.0%
   over CD4 off**. Past it the curve falls monotonically to +1.4% at `f = 40`. §4.1.
6. **So the lever's ABSOLUTE ceiling is +32.4% on PLOLO and +46.0% projected onto the board,
   against the leaders' measured +106%.** The dose can account for at most **30.6%** of their
   move on PLOLO rates, and only at **−17.2σ** on `pds`. §4.2.
6b. 🔴 **`cd4f034`'s board result landed mid-run: +0.0130557, rank 7, 3.0× its prediction — and
   re-pricing the whole ladder on its realised rates gives the SAME verdict from a different
   direction.** The measured factors are **k_cos = 1.384** (not INBOX-2's 1.80, which compares
   212's *macro per-target* cosine with 197's *energy-weighted aggregate* — the same functional
   trap as §2.1's four gates) and **k_pds = 1.987** (INBOX-2's reading, confirmed). The `mse`
   gain grows 1.38× and the `pds` cost grows 1.99×, **and they cancel**: on board rates
   `f = 0.34` is again the argmax and **no dose on the 31-point grid clears the +0.0019 gate**.
   INBOX-2's `f ≈ 2.6` target is falsified by the turnover — at 2.6 the projected board cosine
   is **lower** than at 1.0. §4.4.
6c. **The board-rate ceiling is board cosine 0.1056 at `f ≈ 1.0` = ×1.460 — 0.799× the leaders'
   0.13211 — and it costs −0.01195 overall to stand there.** §4.4.
7. **The incumbent is not merely at the constrained optimum — it is the MAXIMUM, under both
   calibrations.** On the
   exact six-member price (188's `nmae` identity, not a proxy), `f = 0.34` is the argmax over
   the entire 31-point grid: the neighbours are **−0.00024** (`f = 0.28`) and **−0.00029**
   (`f = 0.40`), and the cosine argmax is **−0.00984**. 207's three published selection
   values reproduce exactly. §6.0.
8. **The library form re-opens and loses again — but 207's stated mechanism was the wrong
   one and is corrected here.** 207 re-solved α on the *pre-whitening* norm; the delivered
   norm stayed at 0.45–0.70×. Re-solving on the **delivered** norm hits **1.0000×** and
   moves the family a long way (`lib1`: cos_agg 0.014980 → 0.016227, `pds` 0.7276 → 0.7352).
   It still loses: on the `pds`+`mse` route the best amplitude-matched library arm is
   **+0.00021** against the residual's **+0.00452**, and every pure-library arm is negative —
   and once the exact `nmae` is added (§6.0) the library arms fall to **−0.05 and worse**. §5.
9. **The 1/k route is the worst of the three.** CD4 as one more equal line (`lib1@adel`) is
   **−31.5% cosine against the incumbent** and needs **α = 4.17** to restore delivered
   amplitude. The arithmetic that predicted this: one more equal line is worth **×1.133**
   under 159's law (√(4/3) = ×1.155), while the dose alone already delivers **×1.180**. §5.3.
10. 🔴 **The one genuinely new mechanism is the OPERATOR REFIT, and it is the best arm here —
    and it still fails three ways.** Putting CD4 into the whitener's own basis at **0.87% of
    the trace** rescues the library form completely: α falls **4.17 → 1.27** and `pds` goes to
    **0.754993, +9.85σ above the incumbent**, the largest `pds` move in this experiment. It
    fails because (a) priced on all six members it is **+0.00394 against the +0.0089 amplitude
    gate — 0.44×**;
    (b) **its permutation control is the arm itself** — `B₄ᵀB₄` is invariant to row
    permutation, verified at 8.5e−16 on all four folds, so the gain is provably not
    target-specific and is GEOMETRY, the channel the board priced at **−0.01352** on a19; and
    (c) it moves three variables at once. §5.4.
10b. **INBOX-3 answered, and it retires an October risk.** The per-context cosine maximisers
    are **A 1.004 / B 1.015 / C 1.036 — a ×1.032 spread** — the curve is flat to six decimals
    across that window (transfer loss **0.000%** in every direction), and **A, the T-lineage
    context, has the LOWEST optimum, not the highest.** `docs/DEF-RUNBOOK.md` §3's
    `VCC_CD4_TARGET_F` row **stays a single global constant; no conditional is needed.** §4.5.
11. **VERDICT: NO-GO. Nothing clears its gate. No slot.** One arm out of everything measured
    is positive against the incumbent at all — the operator refit, at **+0.00394 against a
    +0.0089 gate** — and it is the one that cannot be falsified. Everything else is negative.
    The incumbent is `cd4f034`, already submitted. **The CD4 territory is closed at
    `f = 0.34`.** The leaders' ×2.06 is not a CD4 dose, not a CD4 form, not a CD4
    re-extraction, and not anything that arrived on 2026-08-28. §8.

---

## 1. Rule 11 — every operator this experiment touched, by path and sha256

| role | path | sha256 | n | max eig |
|---|---|---|---:|---:|
| **production / `kolf2` fold** | `data/donors/whitener_union_pruned.npz` | **`e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123`** | 17,942 | 5.894570 |
| PLOLO fold `k562_gwps` | `experiments/062-cd4-post-whiten/cloud-output/whiteners/held_k562_gwps.npz` | `13d988fcc3b129d2…` | 17,806 | 8.425882 |
| PLOLO fold `orion_hct116` | `…/held_orion_hct116.npz` | `7cecbf7f36e40d22…` | 17,942 | 8.036340 |
| PLOLO fold `orion_hek293t` | `…/held_orion_hek293t.npz` | `97e486e0d1a805e1…` | 17,942 | 3.695842 |

Production's basis is **29,394 × 17,942** read from `whitener_union_pruned.json`
(`basis: [29394, 17942]`, `basis_prune: essential`, rows kept 7,784 / 10,807 / 10,803) —
asserted, not inferred, before the file was loaded. §6's refit operators are derived from
these four and are identified by the parent sha plus the stated mixing weight; none is
written to disk.

The CD4 source is pinned separately:
`experiments/199-def-runbook/work/GWCD4i.DE_stats.h5ad`,
sha256 `c355f535ff32cf7ba1edc49cf9c6039fe84f2c9ebe4d005515cba75790cfbb62`, 16,786,240,107 B.

---

## 2. Phase 0 — known-answer rows. **All pass: 20 instrument gates, the byte-for-byte bank re-extraction, and 215's two windows.**

`gate216.py` → `work/gate216.json`; `regate216.py`; `verify216.py`; the bank re-extraction.
Two further gate sets are carried inside the analyses they protect: **5/5** on the 1/k law
(§5.3) and **2/2** on the library form (§5.2).

### 2.1 The instrument — 20 gates

| # | gate | expected | measured | verdict |
|---|---|---:|---:|---|
| KA-0 | production operator sha256 + basis | `e2f65ecf` / [29394, 17942] | identical | **PASS** |
| KA-1 | all four PLOLO folds **rebuilt from scratch** vs 207's caches | max\|diff\| = 0 | **0.000e+00** on 4/4 | **PASS** |
| KA-2 | PLOLO a15 `pds`, `k562_gwps` fold | 0.79578409 | **0.79578409** | **PASS** |
| KA-3 | PLOLO a23 `pds`, same fold | 0.78912131 | **0.78912131** | **PASS** |
| KA-3b | 204's `q`, same fold | 0.50164311 | **0.50164311** | **PASS** |
| KA-3c | 204's `q_w`, same fold | 0.50487994 | **0.50487994** | **PASS** |
| KA-4 | a15 → a23 mass-fix price, scaled, 4-fold | **+0.004014** | **+0.004014** | **PASS** |
| KA-5 | 207's sweep at `f` = 0, 0.10, 0.34 × {`cos_agg`, `pds`} | 6 published values | all 6 exact to 6 d.p. | **PASS ×6** |
| KA-6 | 212's macro `cos` at `f = 0.10` | 0.0218418 | **0.0218418** | **PASS** |
| KA-6b | 212's macro `cos` at `f = 0.34` | 0.0248622 | **0.0248622** | **PASS** |
| KA-7 | 212's forward dose gain | **+13.83%** | **+13.83%** | **PASS** |
| KA-8 | 197's one-sided floor at a15's raw `mse` 0.8738045 | 0.07221 | **0.07221** | **PASS** |

> **A trap worth recording, because it made four correct gates look like failures.**
> `score_effects` produces every endpoint twice — `a23_*` (production's mass fix) and `a15_*`
> (without). **207 publishes the a23 tag; 204 and 212 publish the a15 tag.** Compared against
> the wrong tag, 204's `q_w` reads 0.46995 against 0.50488 and 212's cosine misses by 0.196%.
> Against the right tag all four are exact. `regate216.py` records the re-scoring; no
> measurement changed.

### 2.2 The CD4 bank re-extraction, byte-for-byte

`experiments/047-cd4-direct/extract_cd4.py` re-run against the pinned 16.8 GB source, output
compared file by file against the shipped banks by sha256:

```
42 files:  38 IDENTICAL,  4 differ
the 4:     cd4_{rest,stim8hr,stim48hr}/meta.json  and  manifest.json
the diffs: source_local_path, panel_csv, genes_csv (absolute paths on this machine)
           and peak_rss_gib.  No data field differs.
```

**All 36 data files plus `gene_map.csv`, `target_row_map.csv` and both `*.txt` axes reproduce
to the byte.** 199's "36 of 36" is confirmed, and the four JSON diffs are enumerated above so
the next agent does not have to re-derive which they are.

### 2.3 215's Phase-2 windows, verified independently

The mandate requires these to be owned by this experiment before anything is built on them.
Re-derived from the raw snapshots in `data/leaderboard/` plus 215's `work/lb-live.json` with
197's identity `cos_floor = (1.023611 − raw_mse) / 2.074495` and **none of 215's code**:

| team | window | floor from → to | ×  | Δ raw `pds` | leak bound | × bound | ρ |
|---|---|---|---:|---:|---:|---:|---:|
| **`cqawesome`** | `m0830b` → `m0903c` | **0.06427 → 0.13211** | **2.056** | +0.04116 | 0.00332 | **12.4** | **0.94** |
| **`jw_13288`** | `willow_34` → `juniper_46` | **0.06333 → 0.13117** | **2.071** | +0.04947 | 0.00332 | **14.9** | **1.13** |
| `Aginglab.com` | `GeroAI_v11` → `GeroAI_v19` | 0.03096 → 0.09887 | 3.193 | +0.00314 | 0.00332 | 0.9 | 0.07 |
| us | first informative → latest | 0.02873 → 0.07291 | 2.538 | +0.00240 | — | 1.1 | 0.08 |

**Every digit of 215's two headline windows reproduces**, including the discriminants. One
divergence, in 215's favour: taking each team's *global* minimum floor rather than its
chronologically first informative row gives `jw_13288` **×2.243** (`cedar_90` → `juniper_46`),
so 215's ×2.07 is the conservative reading. **The number this experiment is measured against
is +106%, and it is now this experiment's own.**

### 2.4 Territory check — is the sweep run in the right configuration?

The mandate warns that if PerturbMatch and the CD4 dose are substitutes, the sweep must be run
with PerturbMatch off. **They are not substitutes.** 214's `work/report214.json` (its
`RESULTS.md` was not yet written when this was read, so the numbers were taken from its JSON,
not from a claim) measures, **at the `cd4f034` dose**:

| arm | Δ raw `pds` | Δ overall | Δcos | per-fold signs |
|---|---:|---:|---:|---|
| **drop PerturbMatch** | **−0.002842** | **−0.001059** | **−2.05%** | `− − 0 −` |

PerturbMatch is still *positive* at the tripled CD4 dose. The sweep below is therefore run
**with PerturbMatch on**, which is production, and no re-run is required.

---

## 3. Phase 1 — the published atlas against our banks. **It is the same object.**

`atlas216.py`, `remote216.py`, `energy216.py`, `coverage216.py`. Two public read-only network
surfaces (Crossref; the S3 bucket the paper's own data-sharing readme names). Nothing posted.

### 3.1 The release, and its dates

`10.1016/j.cell.2026.08.002` — *"Genome-scale perturb-seq in primary human CD4+ T cells maps
context-specific regulators of T cell programs and human immune traits"*, Cell, Crossref
record **created 2026-08-28T15:06:03Z**, deposited 2026-08-29, VOR licence **CC BY-NC 4.0**
effective 2026-08-05. First author Zhu Ronghui.

The full bucket inventory, read **before** any bulk download (mandate: report before
exceeding 10 GB — this experiment never approached it):

| object | size | LastModified |
|---|---:|---|
| `D{1..4}_{Rest,Stim8hr,Stim48hr}.assigned_guide.h5ad` × 12 | 110–161 GiB each, **1,684 GiB total** | **2025-12-19** |
| `GWCD4i.DE_stats.h5ad` ← **our source** | **15.63 GiB** | **2026-05-28** |
| `GWCD4i.DE_stats.by_donors.h5mu` | 15.71 GiB | **2026-05-28** |
| `GWCD4i.DE_stats.by_guide.h5mu` | 27.40 GiB | **2026-05-28** |
| `GWCD4i.pseudobulk_merged.h5ad` | 41.51 GiB | **2026-05-28** |
| `data_sharing_readme.md` | 26 KB | 2026-05-29 |
| `suppl_tables/*` (28 files) | 170 MB | **2026-08-24** |
| `metadata/*.jsonld` (12) | — | 2025-12-20 / 12-23 |

**Objects by modification date: 2025-12-19 (12), 12-20 (8), 12-23 (4), 2026-05-28 (4),
05-29 (1), 08-24 (28). Objects modified inside 2026-08-26 → 09-03: ZERO.**

### 3.2 Is it the same bytes? Yes — 12 ranges, offset 0 to 16.79 GB

| | |
|---|---|
| local pinned file | 16,786,240,107 B, sha `c355f535…` |
| remote object | 16,786,240,107 B, ETag `c9ff52fc…-2002` (multipart), LastModified 2026-05-28T23:43:31Z |
| byte ranges probed | **12** × 64 KiB at offsets 0, 4,096, 5.60 G, 5.90 G, 8.39 G, 10.00 G, 11.19 G, 11.25 G, 12.26 G, 15.44 G, 16.27 G, 16.79 G |
| **result** | **12 / 12 MATCH** |

**The `.h5ad` the publication points at is the file our banks were built from, unchanged since
three months before publication.**

The one artefact that *did* move — `DE_stats.suppl_table.csv`, replaced 2026-08-24 — is a
tabular copy of that file's `.obs`. It went 16 → 28 columns:

- **added (14):** `neighboring_gene_KD`, `distal_offtarget_flag`, `low_target_gex`,
  `n_guides`, `single_guide_estimate`, `guide_correlation_signif(+_pval)`,
  `guide_correlation_all(+_pval)`, `guide_n_signif_ontarget`,
  `donor_correlation_all_{mean,min}`, `donor_correlation_hits_{mean,min}`
- **dropped (2):** `offtarget_flag`, `ontarget_effect_category`
- **rows:** 33,983 → 33,983, **identical (`target_contrast`, `culture_condition`) set**
- **`guide_kd_efficiency` and `K562_comparison` supplementary tables: sha256-identical to the
  copies we downloaded on 2026-08-22.**

**All 14 "new" columns are already `.obs` columns of the pinned h5ad** (verified directly
against the local file). **Zero columns are new to us.** The 08-24 refresh brought a CSV into
sync with an h5ad that had carried those fields since May.

### 3.3 Target coverage against the 300 panel — we are at the ceiling

| | union | all three conditions | Rest | Stim8hr | Stim48hr |
|---|---:|---:|---:|---:|---:|
| **our banks** | **291** | **276** | 283 | 287 | 284 |
| **the published DE table** | **291** | **276** | 283 | 287 | 284 |
| the pseudobulk object, `keep_for_DE` | 291 | — | 283 | 287 | 286 |
| the pseudobulk object, all rows | 297 | — | — | — | — |
| the sgRNA library (designed / validated) | 297 | — | — | — | — |

The nine panel genes we never cover — `ABCD1`, `C5orf22`, `EEF1A2`, `EPHB2`, `FZD2`,
`KIF21B`, `LENG1`, `NICN1`, `REV1` — have **zero rows** in the published per-perturbation
table. Six of the nine are in the guide library; they produced no DE-eligible result. The
three that are not in the library (`EEF1A2`, `EPHB2`, `FZD2`) were never targeted.

**There is no target coverage to gain.** The only difference anywhere is Stim48hr 284 → 286 in
the pseudobulk's own eligibility flags, i.e. **two targets**, which would move
all-three-conditions coverage by at most 276 → 278.

### 3.4 The gene axis — the one real difference, and it is not new

| | genes | on the official 18,533 axis |
|---|---:|---:|
| `GWCD4i.DE_stats.h5ad` (our source) | 10,282 | **10,113** |
| `GWCD4i.DE_stats.by_guide.h5mu` (per-guide) | 10,282 | **10,113** — *same axis, no gain* |
| **`GWCD4i.pseudobulk_merged.h5ad`** | **18,129** | **17,771** |
| extra official genes in the pseudobulk | | **7,658** (strict superset: 0 DE genes missing) |

The pseudobulk carries 278,684 samples (guide × donor × condition, 4 donors, 12,731 unique
perturbed genes, median 2 guides per panel target) at raw UMI counts, so it also allows a
different normalisation and a different aggregation. **It has been public since 2026-05-28.**
Reading its `.obs` and `.var` cost **16.9 MB over HTTP range requests**; the per-guide
modality cost 14.0 MB.

**How much scored energy is behind that difference** (`energy216.py`, four PLOLO folds, truth
side, under production's CPM > 10 gate, mean over the three challenge frames):

| | share of gated truth energy |
|---|---:|
| reachable by our CD4 gene axis today | **0.8713** |
| the pseudobulk's extra 7,658 genes | **0.1194** |
| unreachable by either | 0.0093 |

### 3.5 …and the extra genes are worth **nothing**. Measured, not assumed.

An energy share is not a cosine. `coverage216.py` measures the elasticity by cutting the CD4
residual's gene support **from the low-expression end** — which is exactly what the extra
genes are, the ones the atlas's own DESeq2 expression filter dropped — at fixed `f`, so the
delivered dose is constant across rows.

| axis energy dropped | gated energy kept | macro `cos_agg` @ `f = 0.34` | rel |
|---:|---:|---:|---:|
| 0% | 0.8714 | 0.023698 | 0.00% |
| 2% | 0.8540 | 0.023787 | **+0.37%** |
| **5%** | 0.8278 | **0.023863** | **+0.70%** |
| 10% | 0.7843 | 0.023853 | +0.65% |
| 15% | 0.7407 | 0.023843 | +0.61% |
| 20% | 0.6972 | 0.023739 | +0.17% |
| 30% | 0.6100 | 0.023509 | −0.80% |
| 45% | 0.4793 | 0.023129 | −2.40% |

**Removing the low-expression tail makes the cosine BETTER** up to a 15% energy cut, on 3 of
4 folds. The measured local elasticity `d log(gain) / d log(reachable energy)` is **−0.554**
(at `f = 0.34`) and **−0.965** (at `f = 0.10`), against a √-law's +0.500.

Extrapolating the pseudobulk's ×1.1370 of reachable energy:

| law | gain at `f = 0.34` | cos_agg | vs today |
|---|---:|---:|---:|
| **measured (local slope −0.554)** | +0.005222 | 0.023313 | **−1.62%** |
| all-points slope +0.205 | +0.005756 | 0.023847 | +0.63% |
| √ law (+0.500) | +0.005979 | 0.024070 | +1.57% |
| **linear (+1.000) — the most generous bound available** | +0.006375 | 0.024466 | **+3.24%** |

**At the most generous bound the pseudobulk re-extraction is worth +3.24% of PLOLO cosine =
+0.0030 of board cosine = +0.0011 overall — under the +0.0019 shape gate — and at the measured
elasticity it is negative.** A 41.5 GiB download and a new bank build buy nothing.

*(By-product, recorded and not pursued: **pruning the low-expression tail off the CD4 residual
is itself worth +0.70% macro cosine** at a 5% energy cut, positive on 3 of 4 folds. That is
+0.00065 of board cosine ≈ +0.00024 overall — 8× under the shape gate. Not an arm.)*

### 3.6 Phase 1 verdict, stated against the mandate's two branches

> *"If the published version is materially different and better, that is the finding… If it is
> the same object, say so — then the date coincidence is a coincidence and Phase 2 is the
> whole experiment."*

**It is the same object.** 12 byte ranges, identical size, unchanged since 2026-05-28. Panel
coverage identical at 291/276. The 08-24 supplementary refresh added zero data. The one
processed object we do not use (the pseudobulk) has been public since May and is measured to
be worth nothing. **The date coincidence is a coincidence, and Phase 2 is the whole
experiment.**

### 3.7 Licence position — RECORDED, NOT RESOLVED (Ahmet's call, as 182 and 205 did)

- The **article** carries **CC BY-NC 4.0** (Crossref, VOR, effective 2026-08-05).
- The **S3 bucket carries no `LICENSE` object** and `data_sharing_readme.md` states no
  licence terms. The 12 `.jsonld` metadata sidecars accompany only the cell-level files.
- Challenge rule §7 permits public lab/literature data to **train or fine-tune** a model and
  forbids using it to **revise** model output. We ship this atlas as a fixed learned residual
  added post-whitening — a trained/conditioned component, not a per-prediction lookup.
- The live tension is **NC vs a prize with a cash component**, which is unchanged by anything
  in this experiment and is not mine to settle. **No action taken.**

---

## 4. Phase 2.1 — the dose curve, to its limit

`sweep216.py` (176 arms: 31 doses + 13 shuffled controls × 4 folds), `analyse216.py`,
`price216.py`. Peak RSS 4.99 GB. Primary endpoint `cos_agg = √(1 − sse)`, the
energy-weighted aggregate cosine 197 inverts off the board; 212's macro `cos` carried
alongside on the a15 tag. Every one of 207's published rows is present and reproduces.

### 4.1 The curve and its turnover

| `f` | fraction | **cos_agg** | rel vs shipped | macro `cos` (a15) | rel | `pds` | Δ`pds`/σ | `q_w` |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0.0000 | 0.018091 | −9.9% | 0.019776 | −9.5% | 0.735107 | −4.69 | 0.488256 |
| 0.05 | 0.0528 | 0.019121 | −4.8% | 0.020871 | −4.4% | 0.739579 | −1.70 | 0.488636 |
| **0.10 ← shipped a15** | 0.1057 | **0.020086** | 0.0% | **0.021842** | 0.0% | **0.742112** | 0.00 | 0.488706 |
| 0.20 ← `pds` argmax | 0.2113 | 0.021802 | +8.5% | 0.023408 | +7.2% | **0.743382** | +0.85 | 0.488788 |
| **0.34 ← `cd4f034`, INCUMBENT** | 0.3593 | **0.023698** | **+18.0%** | **0.024862** | **+13.8%** | 0.740275 | **−1.23** | 0.488959 |
| 0.40 | 0.4226 | 0.024330 | +21.1% | 0.025269 | +15.7% | 0.738481 | −2.43 | 0.488986 |
| 0.55 | 0.5811 | 0.025490 | +26.9% | 0.025868 | +18.4% | 0.733126 | −6.01 | 0.489247 |
| 0.70 | 0.7396 | 0.026166 | +30.3% | **0.026068 ← macro-cos argmax** | +19.3% | 0.727455 | −9.81 | 0.489512 |
| 0.80 ← 207 stopped here | 0.8453 | 0.026416 | +31.5% | 0.026059 | +19.3% | 0.723664 | −12.35 | **0.489590 ← `q_w` argmax** |
| 0.90 | 0.9510 | 0.026550 | +32.2% | 0.025982 | +19.0% | 0.719814 | −14.92 | 0.489573 |
| **1.00 ← cos_agg argmax** | 1.0566 | **0.026597** | **+32.4%** | 0.025860 | +18.4% | 0.716369 | **−17.23** | 0.489395 |
| 1.30 | 1.3736 | 0.026420 | +31.5% | 0.025361 | +16.1% | 0.707224 | −23.35 | 0.489279 |
| 2.00 | 2.1132 | 0.025436 | +26.6% | 0.024149 | +10.6% | 0.692686 | −33.08 | 0.489167 |
| 3.00 | 3.1699 | 0.024203 | +20.5% | 0.022915 | +4.9% | 0.681086 | −40.85 | 0.488713 |
| 7.00 | 7.3964 | 0.022071 | +9.9% | 0.020943 | −4.1% | 0.664866 | −51.70 | 0.488267 |
| 15.0 | 15.849 | 0.020998 | +4.5% | 0.019977 | −8.5% | 0.657353 | −56.73 | 0.487956 |
| **40.0** | 42.265 | 0.020363 | **+1.4%** | 0.019413 | −11.1% | 0.653192 | −59.52 | 0.487548 |

*(the full 31-dose grid with paired-target bootstrap CIs is in `work/analyse216.json`;
`d_pds/σ` uses 162's LOO σ = 0.001494 raw, measured against the shipped `f = 0.10`.)*

**The four maximisers, located by a quadratic in log f through the argmax and its neighbours:**

| endpoint | `f*` | value | rel vs shipped | rel vs CD4 off |
|---|---:|---:|---:|---:|
| **`cos_agg` (197's functional) — PRIMARY** | **1.018** | 0.026598 | **+32.42%** | **+47.02%** |
| macro `cos` (212's functional) | 0.739 | 0.026075 | +19.35% | +31.81% |
| PLOLO `pds` | 0.175 | 0.743514 | +0.17% | +1.13% |
| `q_w` | 0.827 | 0.489593 | +0.18% | +0.27% |

**The curve does turn over, and it turns over well inside the range.** There is a genuine
interior maximum at `f ≈ 1.0` — the dose at which the residual's median row norm equals the
whitened base's — and beyond it the cosine decays monotonically for two more decades. The
answer to the mandate's "beyond if the curve is still rising" is: **it is not still rising.**

### 4.2 🔴 The arithmetic that decides this experiment

| quantity | value |
|---|---:|
| the leaders' measured move (§2.3, this experiment's own numbers) | **×2.056 / ×2.071 = +106%** |
| the CD4 dose at its absolute cosine maximum, vs shipped | **+32.4%** |
| **fraction of the leaders' move the dose can explain, at ANY dose** | **30.6%** |
| the dose inside 162's Gate-1, vs shipped | +18.0% |
| **fraction of the leaders' move the dose can explain inside the guard** | **17.0%** |
| implied board cosine at the cosine argmax (`k = 1`, conservative) | **0.10494** |
| the same, at 212/207's measured PLOLO→board factor `k_cos` = **1.2598** | 0.11161 |
| the leaders' floor | **0.13211** |
| **share of the direction gap closed at the cosine argmax** | **48.6%** |
| **share of the direction gap closed inside Gate-1 (`f = 0.34`)** | **27.0%** |

`k_cos` = **1.2598**, re-derived here from this experiment's own a14 (fraction 0.05) and a15
(f 0.10) rows — 207 published **1.2598**.

> **Superseded in §4.4, and by a board point rather than an argument.** `cd4f034` landed while
> this experiment was running, so the PLOLO→board cosine factor is now *measured* at **1.384**
> for this arm family and the `pds` factor at **1.987**. §4.4 re-prices the whole ladder on
> those and reaches the same operating point. The `k = 1.0` and `k = 1.2598` columns above are
> kept because they are what 207's shipped decision was made on, and because a conclusion that
> survives three independent calibrations is worth more than one stated under the best of them.

**Not one dose anywhere on a 31-point grid spanning `f` = 0 to 40 doubles our cosine.** The
lever's ceiling is a third of what the board measured the leaders doing, and the leaders were
*below* us five days before they did it.

### 4.3 The priced grid, and the incumbent's position on it

096/197's identity through `X` only (`claim`, `Kc`, `ΣD` and `A` all fixed — the emission law,
spread, `keep_top` and `min_cpm` are untouched and the residual's preserve scalar restores the
base's median row norm), `pds` through PLOLO ×1.011, `nmae` through 187 §5's gene-wise
identity at fixed α:

| `f` | ovr `pds` | ovr `mse` | ovr `nmae` | **TOTAL (k=1)** | TOTAL (cal) | inside Gate-1? |
|---:|---:|---:|---:|---:|---:|:---:|
| 0.20 | +0.00047 | +0.00247 | +0.000001 | **+0.00295** | +0.00359 | **Y** |
| 0.28 | −0.00009 | +0.00414 | +0.000001 | **+0.00405** | +0.00513 | **Y** |
| **0.34 ← INCUMBENT** | −0.00068 | +0.00521 | +0.000002 | **+0.00452** | **+0.00588** | **Y** |
| 0.40 | −0.00135 | +0.00612 | +0.000002 | **+0.00477** | +0.00636 | N (−2.43σ) |
| 0.55 | −0.00335 | +0.00779 | +0.000005 | +0.00444 | +0.00647 | N |
| 0.80 | −0.00688 | +0.00912 | +0.000008 | +0.00225 | +0.00462 | N |
| 1.00 | −0.00960 | +0.00938 | +0.000006 | **−0.00021** | +0.00223 | N (−17.23σ) |
| 2.00 | −0.01843 | +0.00771 | +0.000004 | −0.01071 | −0.00871 | N |

207 published its in-sample argmax at `f = 0.40`, **+0.00477** — reproduced here exactly, and
its Gate-1 choice `f = 0.34` at +0.00452 against 207's +0.00437 (207 used 188's exact `nmae`
identity where this table uses 187's proxy; §6 runs the exact one).

**The critical re-measurement the mandate asks for: the guard, re-referenced to the CURRENT
incumbent.** 207's guard was measured against the shipped `f = 0.10`. `cd4f034` is now on the
board, so the guard must be re-taken from `f = 0.34`:

| step from `f = 0.34` | Δ raw `pds` | σ | inside Gate-1? | Δ overall, proxy `nmae` | **Δ overall, EXACT `nmae` (§6.0)** |
|---|---:|---:|:---:|---:|---:|
| → `f = 0.40` | −0.001794 | **−1.20** | **Y** | +0.00025 | **−0.00029** |
| → `f = 0.47` | −0.004301 | −2.88 | N | +0.00019 | — |
| → `f = 0.55` | −0.007149 | −4.79 | N | −0.00008 | **−0.00221** |
| → `f = 1.00` | −0.023906 | −16.00 | N | −0.00473 | **−0.00984** |

**The only further step the guard permits is `f = 0.34 → 0.40`, and on the exact `nmae`
identity it is NEGATIVE (−0.00029); on the proxy it is +0.00025, 7.6× under the +0.0019 shape
gate. Either way it is vetoed.** The dose is exhausted and `cd4f034` is the operating point —
§6.0 shows it is in fact the argmax of the whole grid.

---

### 4.4 🔴 INBOX-2 — re-priced on the REALISED BOARD RATES. **Same verdict, from a completely different calibration.**

`cd4f034`'s board result landed during this run (`experiments/210-cd4-build/BOARD-RESULT.md`,
entry `QSe3JlNUYZ2BzvlcsL14`): **+0.0130557, rank 7 — 3.0× its own preregistered +0.00437.**
Both transfer factors this experiment's grid was built on are now measurable against a board
point, so `board216.py` re-prices the whole ladder against them.

**A correction to INBOX-2 first, because it changes the extrapolation by 30%.** INBOX-2 reads
the PLOLO→board factor as **1.80** from *"212 measured +13.83%, the board delivered +24.9%"*.
Those are **two different cosine functionals**: 212's +13.83% is the **macro per-target**
cosine on the a15 tag (§2.1 KA-6), while 197's board floor inverts the **energy-weighted
aggregate** `X/(√ΣD · mean_c √A_c)`, which is `cos_agg`. On the matched functional:

| | `f = 0.10` | `f = 0.34` | relative |
|---|---:|---:|---:|
| PLOLO `cos_agg` | 0.020086 | 0.023698 | **+17.98%** |
| board floor, `(1.023611 − raw mse)/2.074495` | 0.072278 | 0.090268 | **+24.89%** |
| | | **k_cos** | **1.3840** |

*(This is the same tag/functional trap that made four of my Phase-0 gates look like failures.
The `mse` leg needs no model at all here: the identity reproduces **both** board points —
0.873670 vs 0.87367 and 0.836350 vs 0.83635 — so board raw `mse` is exactly affine in the
projected cosine.)*

**INBOX-2's second factor is confirmed and is the one that matters.**

| | PLOLO | board | ratio |
|---|---:|---:|---:|
| Δ raw `pds`, `f` 0.10 → 0.34 | −0.001837 | **−0.00365** | **k_pds = 1.987** |

207 and 162 priced `pds` at ×1.011. **The realised cost is 1.97× that**, exactly as INBOX-2
says. So the correct re-pricing boosts the `mse` gain by 1.384× **and** the `pds` cost by
1.987× — and they very nearly cancel:

| `f` | PLOLO `cos_agg` | **board cosine** | board raw `mse` | Δ raw `pds` (board) | ovr `pds` | ovr `mse` | **TOTAL vs `cd4f034`** | × leader |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 (CD4 off) | 0.018091 | 0.060710 | 0.89767 | −0.01027 | −0.00379 | −0.01076 | **−0.01454** | 0.460 |
| 0.10 (a15) | 0.020086 | 0.071225 | 0.87585 | +0.00365 | +0.00135 | −0.00693 | −0.00558 | 0.539 |
| 0.20 | 0.021802 | 0.080273 | 0.85709 | +0.00617 | +0.00228 | −0.00364 | −0.00136 | 0.608 |
| 0.28 | 0.022958 | 0.086368 | 0.84444 | +0.00318 | +0.00117 | −0.00142 | −0.00025 | 0.654 |
| **0.34 ← `cd4f034`** | **0.023698** | **0.090268** | **0.83635** | 0 | 0 | 0 | **0** | **0.683** |
| 0.40 | 0.024330 | 0.093602 | 0.82943 | −0.00356 | −0.00131 | +0.00121 | **−0.00010** | 0.709 |
| 0.55 | 0.025490 | 0.099713 | 0.81676 | −0.01420 | −0.00524 | +0.00344 | −0.00180 | 0.755 |
| 0.80 | 0.026416 | 0.104597 | 0.80662 | −0.03300 | −0.01217 | +0.00521 | −0.00696 | 0.792 |
| **1.00 ← cosine argmax** | 0.026597 | **0.105551** | 0.80465 | −0.04750 | −0.01751 | +0.00556 | **−0.01195** | **0.799** |
| 2.00 | 0.025436 | 0.099433 | 0.81734 | −0.09455 | −0.03486 | +0.00333 | −0.03153 | 0.753 |
| 2.60 (INBOX-2's target) | ≈0.0249 | ≈0.0958 | ≈0.8248 | ≈−0.109 | ≈−0.040 | ≈+0.0021 | **≈−0.038** | ≈0.725 |
| 40 | 0.020363 | 0.072685 | 0.87283 | −0.17302 | −0.06380 | −0.00640 | −0.07020 | 0.550 |

**Argmax of the board-rate price: `f = 0.34`. Doses clearing the +0.0019 shape gate: NONE.**

> **INBOX-2's corrected extrapolation is falsified by measurement, not by a rate argument.**
> It projects that +106% needs `f ≈ 2.6` at the board slope. But the *curve turns over at
> `f ≈ 1.0`* (§4.1), so the board slope stops applying long before 2.6 — **at `f = 2.6` the
> board-projected cosine is 0.0958, LOWER than the 0.1056 at `f = 1.0`**, and the arm prices at
> ≈ −0.038. **The dose's board-rate ceiling is ×1.4603 on cosine (0.072278 → 0.105551), which
> is 0.799× the leaders' 0.13211 and 71% of their ×2.056 move — and it costs −0.01195 overall
> to stand there.** INBOX-1's "excluded by two orders of magnitude" was too strong; INBOX-2's
> "back in play at f ≈ 2.6" is too weak; **the measured answer is that the lever reaches 80% of
> where they already stand and cannot be afforded past `f = 0.34`.**

**Two things the board says that this grid deliberately does not extrapolate**, recorded
because they are the most favourable numbers available and I am not crediting them:

| member | realised, `f` 0.10 → 0.34 | what PLOLO predicted |
|---|---:|---|
| `nmae` | **+0.001762** overall | 188's exact identity said **−0.000152** — *sign inverted* |
| `fid` + `reach` + `jac` | **+0.006103** overall | unpriced (rule 9); 207's falsifier 5 fired **favourably**, 12× the a14→a15 pair |

Those two together are **+0.0079 of the realised +0.0131**, i.e. **60% of `cd4f034`'s board
gain came from members no local instrument could price.** That is the honest reason the
prediction missed by 3.0× — and it is a reason to distrust *any* six-member projection in this
family, in **both** directions. It does not rescue `f = 0.40`: that step's priced legs are
−0.00010 and its `pds` cost is now known to be 2× what PLOLO says.

### 4.5 INBOX-3 — the dose curve per context. **The spread is 3.2% and the risk is retired.**

`ctx216.py` re-runs the identical grid keeping the three challenge frames separate (known-answer
row: the pooled column reproduces `sweep216.json` at **max\|diff\| = 0.000e+00**), `ctxanalyse216.py`.

| context | lineage | cosine maximiser `f*` | value |
|---|---|---:|---:|
| **A** | **JURKAT — T-lymphoblastic, the donor's own lineage** | **1.004** | 0.026345 |
| B | HELA — cervix | 1.015 | 0.026578 |
| C | CAL33 — head-and-neck squamous | 1.036 | 0.026868 |
| **pooled (what we would ship)** | | **1.018** | 0.026597 |

1. **The spread between the three maximisers is `f` 1.004 → 1.036 = ×1.032.** Three percent.
2. **A is the LOWEST, not the highest.** The dose optimum is not pulled up by lineage match, so
   **lineage is irrelevant to the dose as well as to the existence of the gain** — INBOX-3 names
   this as the better of the two outcomes, and it is the one that happened. (Note the *level*
   ordering runs the other way: C > B > A at every dose, which is truth quality, not lineage.)
3. **The transfer loss from shipping the pooled optimum into any single context is 0.000%**, and
   the worst cross-context loss (evaluating each context at another's maximiser) is likewise
   **0.000%** — the curve is flat to six decimals across the whole 1.004–1.036 window.
4. The per-context `pds` cost curves are near-identical: at `f = 0.40`, board-projected
   Δ raw `pds` is **−0.00391 / −0.00223 / −0.00455** for A / B / C; at `f = 1.00`,
   −0.04860 / −0.04498 / −0.04891. **B is consistently the cheapest and A and C track each
   other**, and the ordering never changes with dose.

> **Consequence for October, stated as INBOX-3 asked.** `docs/DEF-RUNBOOK.md` §3's row for
> `VCC_CD4_TARGET_F` **stays a single global constant. No conditional is needed.** A ×1.032
> spread with a curve that is flat to six decimals across it cannot produce an `abc-knee`-shaped
> failure on D/E/F, whatever their lineages turn out to be. **This is a risk retired, and it was
> retired for the cost of a `groupby`.**

---

## 5. Phase 2.2/2.3 — the FORM, re-opened with the α re-solve corrected

`form216.py` (252 arms: 21 configurations × 3 α policies × 4 folds), `formanalyse216.py`.
Peak RSS 6.18 GB.

### 5.1 207's stated mechanism was the wrong α, and this experiment fixes it

207 §5.4 point 2 blamed the library form's failure on delivered amplitude collapsing to
**0.523×** even at its re-solved α — and its re-solve (`library207.py:103`) matched the
**pre-whitening** median row norm. That is not the norm that collapsed. Three policies:

| policy | what α matches | delivered ratio at `lib1`, macro |
|---|---|---:|
| `@a0` | nothing (α = 0.982) | **0.3195×** |
| `@apre` | 207's pre-whitening median row norm | **0.3884×** |
| **`@adel`** | **the DELIVERED median row norm, the object the board reads** | **1.0000×** |

**The correction is large inside the library family.** At `lib1`, going `@a0` → `@adel` moves
cos_agg **0.014980 → 0.016227 (+8.3%)** and `pds` **0.727623 → 0.735226 (+0.0076 raw)**. So
the mandate's suspicion was right: 207 did not re-solve the α that mattered.

### 5.2 It does not rescue the form

Known-answer rows first: **`lib0@a0` reproduces `f0` (CD4 off) at max\|diff\| = 0.000e+00**,
and **`lib0.5@a0` on `orion_hct116` reads cos_agg 0.024300 against 207 §5.5's published
0.024304** — the one fold on which 207 found the library form *not* bad.

Amplitude-matched arms only (`@adel`, delivered ratio 1.0000×), macro over four folds:

| arm | α | cos_agg | vs shipped | vs `cd4f034` | `pds` | Δ`pds`/σ vs inc. | **ovr (k=1)** |
|---|---:|---:|---:|---:|---:|---:|---:|
| `lib0` ≡ CD4 off | 0.984 | 0.018091 | −9.9% | −23.7% | 0.735095 | −3.47 | −0.00549 |
| `lib0.125` | 1.252 | 0.018290 | −8.9% | −22.8% | 0.739048 | −0.82 | −0.00373 |
| `lib0.25` | 1.694 | 0.017991 | −10.4% | −24.1% | 0.739849 | −0.29 | −0.00386 |
| `lib0.5` | 2.626 | 0.017379 | −13.5% | −26.7% | 0.739209 | −0.71 | −0.00498 |
| **`lib1` — the 1/k route** | **4.170** | **0.016227** | **−19.2%** | **−31.5%** | 0.735226 | −3.38 | **−0.00813** |
| `lib0.25 + res f=0.20` | 1.590 | 0.021853 | +8.8% | −7.8% | 0.735841 | −2.97 | **+0.00021 ← best** |
| `lib0.25 + res f=0.55` | 1.393 | 0.023960 | +19.3% | **+1.1%** | 0.713225 | −18.11 | −0.00518 |
| `lib0.5 + res f=0.34` | 1.962 | 0.021804 | +8.6% | −8.0% | 0.709046 | −20.90 | −0.00985 |
| **residual alone, `f = 0.34` (INCUMBENT)** | **0.982** | **0.023698** | **+18.0%** | 0.0% | **0.740275** | 0.00 | **+0.00452** |

**Every pure-library arm is below CD4-off-with-a-corrected-α except `lib0.125`, and the best
composed arm is +0.00021 against the residual's +0.00452 — 21× worse.** The single
configuration that beats the incumbent on cosine (`lib0.25 + res f=0.55`, +1.1%) pays
**−18.1σ** on `pds` and prices at −0.00518.

**The mechanism, now stated at matched amplitude, which 207 could not do:** the library
form's cosine is **monotonically decreasing in `w` above `w ≈ 0.125`**, and at `w = 1` it is
**below CD4 off**. Mixing CD4 into the pre-whitening consensus is worse than not having CD4,
because the operator then whitens it — 207 §5.2's finding, confirmed with the amplitude
confound removed rather than merely asserted.

### 5.3 The 1/k route, and the arithmetic that predicted it

`keff216.py`, 159's law as 207 §7.2 calibrated it (`1/cos²_surface = 66.5 + 1920/k`, board→
surface ÷2.078). **Known-answer gates: 5/5 pass** — it returns 207's published `k_eff` of
2.46 (CD4 off), 3.09 (shipped), 4.49 (`f = 0.34`), 3.84 (best 4-line equal consensus) and
10.41 (`jw_13288`'s floor).

| configuration | board cosine | `k_eff` |
|---|---:|---:|
| CD4 off | 0.07138 | **2.46** |
| shipped `f = 0.10` | 0.079249 | **3.09** |
| `f = 0.34` — INCUMBENT | 0.093502 | **4.49** |
| `f = 1.00` — cosine argmax, at −17.2σ `pds` | 0.104940 | **5.90** |
| `f = 1.00`, `k_cos`-calibrated | 0.111607 | 6.85 |
| best equal consensus of every panel-covering public line (197 §4.2) | 0.08730 | 3.84 |
| **`cqawesome` floor** | **0.13211** | **10.61** |

| move | × on board cosine |
|---|---:|
| ONE more equal line, from our `k_eff` = 3.09 → 4.09 | **×1.133** |
| 159's √(4/3) rule of thumb for 3 → 4 lines | ×1.155 |
| **the DOSE alone at `f = 0.34`** | **×1.180** |
| the DOSE at its cosine argmax | ×1.324 |
| **the leaders' measured move** | **×2.056** |

**The 1/k route cannot explain ×2.06, and the measurement confirms it cannot even deliver its
own ×1.13** — `lib1@adel` lands at ×0.685 of the incumbent, because the fourth line goes in
*before* the whitener. The dose we already ship is worth more than a fourth equal line
(×1.180 vs ×1.133), which is 207's reading — CD4 is a better-than-one-line estimator of the
shared component — confirmed forward here.

**To stand where the leaders stand needs `k_eff` = 10.61 against our 3.09: 7.5 more effective
lines.** 166 proved there is no fifth panel-covering public CRISPRi screen. Whatever they did,
it is not more donor lines and it is not this atlas at a different dose or in a different form.

### 5.4 The operator RE-FITTED with CD4 in its basis — **the best arm in this experiment, and it still fails, three ways**

`op216.py` (40 arms; 8 eigendecompositions of 17,806²/17,942², ~6–7 min each), `perminv216.py`,
`opanalyse216.py`. Peak RSS **11.70 GB**.

**Construction, stated because rule 11 forbids inferring an operator.** The stored `.npz`
holds the eigendecomposition of `C = Bᵀ B / n`, so `C` is recoverable exactly. CD4's own Gram
`C₄ = B₄ᵀ B₄ / m` (m = 251–271 carrying rows, on the fold's `wcols` axis) is mixed in at a
stated **energy** share, which avoids depending on donor row-count conventions the frozen
files do not record:

```
C' = trace(C) · [ (1−β)·C/trace(C)  +  β·C₄/trace(C₄) ]
```

Two β per fold: **`natural`** = m/(n_basis+m), the weight CD4's rows would get if simply
appended (0.0070–0.0095); and **`equal`** = 1/(k+1), CD4 as one more equal line (0.25–0.333).

**Known-answer gates: 2 per fold, 8/8 pass.** (a) reconstruct-then-re-eigendecompose
preserves `(C)^(−γ/2)` at this operator's real conditioning — tested on a synthetic matrix
carrying the fold's own quantile-matched spectrum at n = 4,000: rel max\|diff\| **1.4–2.6e−07**;
(b) trace conservation on the full problem: **8.9e−15 to 8.7e−11**. Every refit operator's
condition number is **1.000e+12**, identical to its parent — the floor is relative, so the
refit does not change *how much* the operator amplifies, only **which directions** it
amplifies.

Macro over four folds, against the incumbent (shipped operator, residual `f = 0.34`,
`cos_agg` 0.023698, `pds` 0.740275):

| arm | β | α | delivered | `cos_agg` | vs inc | `pds` | Δ`pds`/σ | ovr `pds` | ovr `mse` | **TOTAL vs inc** | class | × gate |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|
| **`refit-natural / lib1`** | 0.0087 | **1.270** | 1.0000 | 0.022777 | −3.9% | **0.754993** | **+9.85** | **+0.00480** | +0.00388 | **+0.00416** | AMPLITUDE | **0.47** |
| `refit-natural / lib0.5` | 0.0087 | 1.150 | 1.0000 | 0.021006 | −11.4% | 0.748802 | +5.71 | +0.00249 | +0.00133 | −0.00070 | AMPLITUDE | −0.08 |
| `refit-natural / res f=0.34` | 0.0087 | 0.982 | 1.0006 | 0.023565 | −0.6% | 0.740140 | −0.09 | −0.00074 | +0.00501 | −0.00024 | shape | −0.13 |
| `refit-natural / res f=0.10` | 0.0087 | 0.982 | 1.0005 | 0.019925 | −15.9% | 0.742022 | +1.17 | −0.00003 | −0.00023 | −0.00479 | shape | −2.52 |
| `refit-natural / CD4 off` | 0.0087 | 0.982 | 0.9978 | 0.017922 | −24.4% | 0.734787 | −3.67 | −0.00273 | −0.00312 | −0.01037 | shape | −5.46 |
| `refit-equal / lib1` | 0.3125 | 1.242 | 1.0000 | 0.019741 | −16.7% | 0.745657 | +3.60 | +0.00132 | −0.00050 | −0.00370 | AMPLITUDE | −0.42 |
| `refit-equal / res f=0.34` | 0.3125 | 0.982 | 1.0020 | 0.022936 | −3.2% | 0.739118 | −0.77 | −0.00112 | +0.00411 | −0.00153 | shape | −0.81 |
| `refit-equal / CD4 off` | 0.3125 | 0.982 | 0.9958 | 0.017155 | −27.6% | 0.732974 | −4.89 | −0.00341 | −0.00422 | −0.01215 | shape | −6.40 |

**The mandate's distinction is real and it is the largest single effect measured here.** An
operator that has *never seen* CD4 destroys the library form — it amplifies CD4's near-null
directions by up to 7,943× (`AGENTS.md` §8.2), which is why `lib1` through the shipped
operator needs **α = 4.17** and still reads `cos_agg` 0.016227 and `pds` 0.735226. The same
form through an operator that HAS seen CD4, at only **0.87% of the trace**, needs **α = 1.27**
and reads **`pds` 0.754993 — +9.85σ above the incumbent**, the largest `pds` move anywhere in
this experiment. 207 §5.2's "CD4 must enter un-whitened" is therefore **incomplete**: CD4 must
not be pushed through an operator blind to it, which is a different and fixable statement.

**And it still fails, three independent ways.**

1. **It is 0.44× its own gate.** α moves 0.982 → 1.270, so rule 13's amplitude class applies
   (+0.0089) and **211's standing rule forbids any amplitude credit** — `claim` is at
   100.0000% of its cap in the binding context. Priced on `pds` + `mse`: +0.00416; with 188's
   exact `nmae` added (§5.5): **+0.00394 against +0.0089.**
2. 🔴 **It cannot be falsified, because its permutation control is the arm itself.**
   `perminv216.py`: the refit depends on CD4 only through `B₄ᵀ B₄`, and
   `(P B₄)ᵀ(P B₄) = B₄ᵀ B₄` for any permutation `P`. Measured on all four folds — **every one
   of the 251/267/271/253 CD4 rows moves under 207's derangement, and the Gram is identical to
   8.5e−16 – 1.2e−15.** So the refit operator built from **target-shuffled** CD4 is the *same
   operator, bit for bit*. **Whatever this arm buys is provably NOT target-specific biology.**
   In 215's channel language it is **GEOMETRY** — and the board has already priced a geometry
   arm in exactly this family: **a19, the whitener ridge, realised −0.01352 overall**, and
   215's instrument classified it GEOMETRY at ρ = 0.02.
3. **It is three variables at once** — a new operator, a new form, and α ×1.29 — so it
   violates rule 5 and could not be shipped as one arm even if it cleared.

Two supporting negatives, both measured: the `equal` weight is worse than `natural` on both
endpoints on every fold (β = 1/3 costs −16.7% cosine at `lib1`), and the refit is
**neutral-to-negative for the residual form** (`res f=0.34`: −0.00024; CD4 off: −0.01037).
**The refit does not add information — it removes an amplification pathology that only the
library form suffers from.**

### 5.5 The refit arm's `nmae`, so rule 12 is closed on the one arm that got close

`opsix216.py`: 188's exact identity and the top-89 ordering precision, re-fitting the
`natural` operator on each fold. **Known-answer rows: the two shipped-operator reference arms
reproduce `six216.py` to five decimals** — a15 `nmae` **0.98476** / precision **0.58995**, and
the incumbent **0.98530** / **0.59752**.

| arm | α | delivered | `nmae` wald89 | Δ vs incumbent | **ovr `nmae`** | top-89 precision | Δ prec | `nmae` lfc89 | `nmae` h89 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| a15 `f = 0.10` shipped | 0.982 | 1.0000 | 0.98476 | −0.00055 | +0.000152 | 0.58995 | −0.00757 | 0.99219 | 0.97873 |
| **`res f=0.34` INCUMBENT** | 0.982 | 1.0009 | **0.98530** | 0 | 0 | **0.59752** | 0 | 0.99239 | 0.97838 |
| `refit-nat / res f=0.34` | 0.982 | 1.0007 | 0.98520 | −0.00010 | **+0.000028** | 0.59782 | +0.00031 | 0.99235 | 0.97807 |
| **`refit-nat / lib1 @adel`** | 1.270 | 1.0063 | 0.98610 | +0.00080 | **−0.000222** | **0.59872** | **+0.00120** | 0.99292 | 0.98146 |
| `refit-nat / lib0.5 @adel` | 1.150 | 1.0012 | 0.98431 | −0.00099 | +0.000275 | 0.59378 | −0.00374 | 0.99195 | 0.97732 |

`nmae` is **not** the killer here — at α = 1.27 the damage is −0.000222, three orders of
magnitude smaller than the −0.1249 the same form suffers at α = 3.85 through the shipped
operator (§6.0b), which is itself the cleanest confirmation that §6.0b's mechanism is right.
All three of 188's pools agree in sign. **Final six-member price of the best arm in this
experiment:**

| member | value vs `cd4f034` | route |
|---|---:|---|
| `pds` | **+0.00548** | PLOLO 4 folds, +0.01472 raw = **+9.85σ** |
| `mse` | **−0.00133** | 096/197 `X`-only identity at matched delivered amplitude, Δcos −3.89% |
| `nmae` | **−0.00022** | 188 exact identity, wald89 pool; lfc89 and h89 agree in sign |
| `reach` / `fid` | unpriced; ordering half **+0.00120** (favourable) | rule 9; h1-test vacuous for the residual and 4/85 for the library form (§6.1) |
| `jac` | unpriced, and not exempt | the call set moves with the profile |
| **TOTAL** | **+0.00394** | against the **+0.0089** amplitude gate → **0.44×** |

**VETO at 0.44× its gate, before the two unpriced DE members, and it could not survive a
permutation control that its own algebra makes vacuous.**

> **A process note worth recording, because a known-answer row is what caught it.** The first
> run of `opsix216.py` built the refit operator *before* scoring the two shipped-operator
> reference arms — and the spectral reconstruction scales the eigenvector matrix **in place**
> (that is what keeps the peak under 8 GB), so those two arms were scored through a wrecked
> operator. Its a15 row read `nmae` **1.05241** against `six216.py`'s **0.98476** and
> precision **0.54319** against **0.58995**. Nothing about the output looked malformed; the
> refit rows were fine and only the baseline was wrong, so **every delta in that table was
> wrong by a plausible-looking amount.** It was caught solely because the a15 row is a known
> answer. The buggy output is kept at `work/opsix216-BUGGY.{json,log}` and the fixed script
> carries the reason in a comment. This is `AGENTS.md` §11.1 rule 10's corollary, exactly.

---

## 6. Phase 3 — six members and the gate

`six216.py` (188's exact `nmae` identity, 4 folds, challenge frame A, production key `|H|`),
`final216.py`. Peak RSS 6.32 GB.

### 6.0 The six-member price, and **`cd4f034` is its argmax over the whole grid**

`price216.py`'s wide grid uses 187 §5's `nmae` proxy, which is ≈0 by construction. The table
below replaces it with **188's exact identity** at each candidate point. Known-answer rows:
207 §8.1's own selection table gives +0.00413 at `f = 0.28`, **+0.00437 at `f = 0.34`** and
+0.00408 at `f = 0.40` — **3/3 reproduce exactly.**

| operating point | ovr `pds` | ovr `mse` | ovr `nmae` | **TOTAL (k=1)** | TOTAL (cal) | Δ top-89 precision | α | class |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| `f = 0.10` (a15, shipped) | 0 | 0 | 0 | 0 | 0 | 0 | 0.982 | shape |
| `f = 0.20` | +0.00047 | +0.00247 | +0.000191 | +0.00314 | +0.00378 | +0.00346 | 0.982 | shape |
| `f = 0.28` | −0.00009 | +0.00414 | +0.000078 | +0.00413 | +0.00520 | +0.00687 | 0.982 | shape |
| **`f = 0.34` ← `cd4f034`** | **−0.00068** | **+0.00521** | **−0.000152** | **+0.00437** | **+0.00572** | **+0.00757** | 0.982 | shape |
| `f = 0.40` | −0.00135 | +0.00612 | −0.000683 | +0.00408 | +0.00567 | +0.00783 | 0.982 | shape |
| `f = 0.55` | −0.00335 | +0.00779 | −0.002279 | +0.00216 | +0.00418 | **+0.00875 ← peak** | 0.982 | shape |
| `f = 0.80` | −0.00688 | +0.00912 | −0.004298 | −0.00205 | +0.00032 | +0.00718 | 0.982 | shape |
| `f = 1.00` (cos argmax) | −0.00960 | +0.00938 | −0.005256 | **−0.00547** | −0.00303 | +0.00567 | 0.982 | shape |
| `lib w=1` (1/k) `@adel` | −0.00257 | −0.00556 | **−0.124934** | **−0.13306** | | +0.01077 | **3.849** | AMPLITUDE |
| `lib w=0.5 @adel` | −0.00108 | −0.00390 | −0.049776 | −0.05476 | | +0.01066 | 2.626 | AMPLITUDE |
| `lib w=0.5 + res f=0.34 @adel` | −0.01233 | +0.00248 | −0.036762 | −0.04661 | | +0.00650 | 1.962 | AMPLITUDE |

**Against the incumbent, every candidate on the board is negative:**

| candidate | Δ overall vs `cd4f034` | gate | × gate | verdict |
|---|---:|---:|---:|---|
| `f = 0.28` | **−0.00024** | +0.0019 | −0.13 | **VETO** |
| `f = 0.40` | **−0.00029** | +0.0019 | −0.15 | **VETO** |
| `f = 0.20` | −0.00123 | +0.0019 | −0.65 | VETO |
| `f = 0.55` | −0.00221 | +0.0019 | −1.16 | VETO |
| `f = 0.80` | −0.00642 | +0.0019 | −3.38 | VETO |
| `f = 1.00` | −0.00984 | +0.0019 | −5.18 | VETO |
| `lib w=0.5 + res f=0.34 @adel` | −0.05098 | +0.0089 | −5.73 | VETO |
| `lib w=0.5 @adel` | −0.05913 | +0.0089 | −6.64 | VETO |
| **`lib w=1 @adel` — the 1/k route** | **−0.13743** | +0.0089 | −15.44 | **VETO** |

**`cd4f034` is the argmax of the exact six-member price over the entire 31-point dose grid
and all three forms.** With the exact `nmae` in place, even the one step the `pds` guard
permits (`f = 0.40`) is **−0.00029**, i.e. **negative**, not the +0.00025 the proxy `nmae`
suggested in §4.3. The dose is not merely exhausted; the incumbent is the maximum.

### 6.0b 🔴 A correction to 207's stated mechanism for the library form's `nmae`

207 §5.4 point 3 attributed the library form's `nmae` explosion to the delivered amplitude
collapsing to 0.523×. **That attribution is wrong, and the corrected α proves it by making
`nmae` far worse:**

| `lib w=1` | delivered / a15 | α | `nmae` (wald89) | ovr `nmae` |
|---|---:|---:|---:|---:|
| `@a0` (amplitude collapsed) | 0.330 | 0.982 | 1.05642 | −0.019907 |
| **`@adel` (amplitude restored)** | **1.048** | **3.849** | **1.43452** | **−0.124934** |

Restoring the delivered amplitude multiplies the `nmae` damage **6.3×**. 187's marginal
refund rate `e` collapses to **0.013** (169 §5's threshold is `e > ½`). The mechanism is not
an amplitude deficit — **the library direction is bad, and bringing it up to full amplitude
amplifies the error rather than the signal.** That is the same conclusion 207 reached, for a
different and now-measured reason.

*(Two footnotes so §5.2 and §6.0 do not look contradictory, and because the second one makes
the library form look BETTER than it is.*

*(a) `form216.py`'s `delivered_ratio` is measured with production's mass fix over all three
frames and reads 1.0000×; `six216.py`'s `dlv/a15` column is measured in frame A only, without
the mass fix, over the operator columns. Different normalisation, same object.*

*(b) The α secant is clipped at 6× its start point. In `form216.py` it starts from the
pre-whitening α and converges on **4 of 4 folds** — per-fold α 1.960 / 4.423 / 7.173 / 3.122,
mean **4.170**, delivered ratio 0.9999–1.0000 everywhere. In `six216.py` it starts from
α = 0.982, so the clip binds on `orion_hek293t` at 5.892 against the 7.173 that fold needs,
and its mean α reads **3.849**. **That fold's `lib w=1 @adel` row is therefore
under-amplified, and since amplifying this direction is what destroys `nmae`, the −0.1249
reported above UNDERSTATES the damage.** Stated rather than fixed, because the arm loses by
15.4× the gate either way.)*

### 6.0c `nmae`'s turn-over, which is what closes the dose past `f = 0.40`

| `f` | 0.20 | 0.28 | **0.34** | 0.40 | 0.55 | 0.80 | 1.00 |
|---|---:|---:|---:|---:|---:|---:|---:|
| `nmae` wald89 | 0.98407 | 0.98448 | **0.98530** | 0.98722 | 0.99296 | 1.00023 | 1.00368 |
| 187's `e` (improves iff `e > ½`) | 0.717 | 0.549 | **0.435** | 0.306 | 0.153 | 0.087 | 0.040 |
| ovr `nmae` vs a15 | +0.000191 | +0.000078 | **−0.000152** | −0.000683 | −0.002279 | −0.004298 | −0.005256 |

**`e` crosses ½ between `f` = 0.28 and 0.34 and is 0.04 by `f` = 1.00.** All three of 188's
pool constructions (`wald89`, `lfc89`, `h89`) agree in sign for every residual arm, so 188's
sign-inversion hazard does not bite. **The `nmae` leg alone costs −0.0053 at the cosine
argmax, which is 55% of the `mse` gain there.**

### 6.1 Non-vacuity, asserted before any member is read

`vacuity216.py` → `work/vacuity216.json`.

**(a) h1-test, residual form — EXACTLY VACUOUS, re-measured not inherited.**

```
h1-test: 85 targets, 5 of them in the 300 panel, CD4 carries all three conditions on 4
         nonzero residual row norms: 4 of 85
         MEDIAN row norm of the residual = 0.0   (mean 0.2539, max 6.7707)
         -> f_per_unit = 0.0 -> solved fraction = 0.0 -> the residual is IDENTICALLY ZERO
```

Mechanism: production's `solve_cd4_fraction` sizes the residual by a **median** row norm over
every row; with 81 of 85 rows empty the median is 0. Production **raises** on this input;
`plolo.solve_fraction` returns 0.0 silently (207 handed that back, unchanged here — no
production file was touched). **So `nmae`, `reach`, `fid` and `jac` cannot be read on h1-test
for any residual-form arm in this experiment. They are reported as UNMEASURED with that
mechanism, exactly as rule 12 requires.**

**(b) h1-test, library form — non-vacuous, but on 4 of 85 rows (4.7%).** The library form is
not sized by a median, so it does change the surface — on 4.7% of it. A member read there is
a 4-row measurement and is reported as unusable rather than as a number.

**(c) the four PLOLO folds — non-vacuous at every dose this experiment reports.**

| | `k562_gwps` | `orion_hct116` | `orion_hek293t` | `kolf2` |
|---|---:|---:|---:|---:|
| targets | 269 | 290 | 294 | 273 |
| rows CD4 carries (all three conditions) | **251** | **267** | **271** | **253** |
| median residual row norm | **4.4677** | **4.4934** | **4.4718** | **4.4491** |
| median whitened base row norm | 4.9114 | 4.9002 | 4.9460 | 4.1403 |

Delivered displacement from the `f = 0` artifact, 4-fold mean Frobenius relative norm:

| `f` | 0.05 | 0.10 | 0.20 | **0.34** | 0.40 | 0.55 | 0.80 | 1.00 | 2.00 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| solved `fraction` | 0.0528 | 0.1057 | 0.2113 | **0.3593** | 0.4227 | 0.5811 | 0.8453 | 1.0566 | 2.1133 |
| relative displacement | 0.0514 | 0.1023 | 0.2017 | **0.3312** | 0.3828 | 0.5001 | 0.6624 | 0.7628 | 1.0495 |

**Every row in this document moved the object by a stated, large amount. None is a silent
zero.** All four of 207's published rows reproduce to five decimals.

---

## 7. Phase 4 — the four falsifiers

### 7.1 Falsifier 1 — target-shuffled CD4 at every dose. **PASSES, and now over 13 doses.**

207's construction imported by value (whole rows permuted by a derangement across targets,
seed 20726; geometry, per-row norms, support pattern and carrying-row count all preserved).
Δ vs CD4 off, macro over four folds; CIs are paired-target bootstraps (6,000 replicates):

| `f` | real Δcos_agg | **shuffled Δcos_agg** | shuffled 95% CI | real Δ`pds` | **shuffled Δ`pds`** | gap |
|---:|---:|---:|---|---:|---:|---:|
| 0.05 | +0.001030 | **−0.000049** | [−0.00021, +0.00001] | +0.004472 | **−0.001992** | +0.001079 |
| 0.10 | +0.001995 | **−0.000149** | [−0.00048, −0.00005] | +0.007005 | **−0.004641** | +0.002144 |
| 0.20 | +0.003711 | **−0.000495** | [−0.00115, −0.00034] | +0.008275 | **−0.012463** | +0.004206 |
| **0.34** | +0.005607 | **−0.001248** | [−0.00228, −0.00101] | +0.005168 | **−0.027585** | +0.006855 |
| 0.40 | +0.006239 | **−0.001643** | [−0.00279, −0.00136] | +0.003375 | **−0.034064** | +0.007882 |
| 0.55 | +0.007399 | **−0.002740** | [−0.00411, −0.00231] | −0.001981 | **−0.049150** | +0.010139 |
| 0.80 | +0.008325 | **−0.004660** | [−0.00619, −0.00395] | −0.011443 | **−0.071111** | +0.012985 |
| **1.00** | +0.008506 | **−0.006090** | [−0.00765, −0.00514] | −0.018737 | **−0.085850** | +0.014596 |
| 1.30 | +0.008329 | **−0.007920** | [−0.00951, −0.00669] | −0.027883 | **−0.103887** | +0.016249 |
| 2.00 | +0.007345 | **−0.010860** | [−0.01253, −0.00926] | −0.042421 | **−0.133905** | +0.018205 |
| 3.00 | +0.006112 | **−0.013071** | [−0.01500, −0.01138] | −0.054020 | **−0.159667** | +0.019183 |
| 7.00 | +0.003980 | **−0.015543** | [−0.01836, −0.01429] | −0.070240 | **−0.199017** | +0.019523 |
| 15.0 | +0.002907 | **−0.016567** | [−0.01982, −0.01557] | −0.077753 | **−0.217373** | +0.019474 |

**The shuffled donor is negative on both endpoints at all 13 doses, with no permutation tie
anywhere**, and its CI excludes zero from `f = 0.10` upward. 207's twelve rows reproduce
exactly and the extension to `f = 15` behaves the same way. The real-vs-shuffled gap widens
monotonically to `f = 7` and then plateaus (0.019523 → 0.019474) as both collapse. **The
cosine gain is target-specific biology at every dose, including the new tail.**

### 7.2 Falsifier 2 — the `pds` guard, at every dose. **Both operating points given.**

σ = 162's LOO **0.001494** raw; Gate-1 threshold **1.645σ**. Against the shipped `f = 0.10`
the largest dose inside the guard is **`f = 0.34`** (−1.23σ; `f = 0.40` fails at −2.43σ) —
**207's operating point, selected independently here.** Against the current incumbent
`f = 0.34` the largest permitted step is **`f = 0.40`** (−1.20σ), worth **+0.00025**.

The cosine optimum (`f ≈ 1.0`) sits **far** past the guard at **−17.23σ**. As the mandate
requires, both are priced and neither is chosen here:

All three priced on the exact six-member route (§6.0), which is the one the coordinator should
read:

| operating point | Δ`pds`/σ vs incumbent | ovr `pds` | ovr `mse` | ovr `nmae` (exact) | **TOTAL vs shipped** | **vs incumbent** | class / gate |
|---|---:|---:|---:|---:|---:|---:|---|
| **A — `f = 0.34` (`cd4f034`, ON THE BOARD)** | 0.00 | −0.00068 | +0.00521 | −0.000152 | **+0.00437** | **0** | shape / +0.0019 — **already shipped, and the grid argmax** |
| **B — `f = 0.40`, the largest the guard allows** | −1.20 | −0.00135 | +0.00612 | −0.000683 | +0.00408 | **−0.00029** | shape / +0.0019 → **VETO, negative** |
| **C — `f = 1.00`, the cosine argmax** | −16.00 | −0.00960 | +0.00938 | −0.005256 | −0.00547 | **−0.00984** | shape / +0.0019 → **VETO, −5.2× the gate** |

**There is no operating point for the coordinator to choose between.** The mandate's
contingency — *"if the cosine optimum is far past the `pds` guard, give both operating points
and let the coordinator choose"* — resolves to a single point: the cosine optimum IS far past
the guard (−16σ), and priced on all six members it is **worse than the incumbent by 2.3× the
gate**, because `pds` (−0.0096) and `nmae` (−0.0053) together exceed the `mse` gain (+0.0094).

### 7.3 Falsifier 3 — leak-free fold, and all four folds

`kolf2` is not a production donor and contributed no row to any production operator, so it is
the leak-free fold. Per-fold cosine curves:

| held fold | lineage | `cos_agg` at shipped | cosine argmax `f` | cosine at argmax | rel |
|---|---|---:|---:|---:|---:|
| `k562_gwps` | myeloid | 0.028762 | 1.10 | 0.039716 | **+38.1%** |
| `orion_hct116` | epithelial | 0.023709 | 1.00 | 0.031432 | **+32.6%** |
| `orion_hek293t` | epithelial | 0.020058 | 0.70 | 0.023574 | **+17.5%** |
| **`kolf2` (leak-free)** | pluripotent | 0.007814 | 1.45 | 0.012135 | **+55.3%** |

**All four folds are positive at the shipped dose and all four have an interior cosine
maximum in `f` ∈ [0.70, 1.45].** The direction of the finding does not depend on any one
fold, and **not one fold — including the most favourable, `kolf2` at +55.3% — reaches the
+106% the board measured.**

### 7.4 Falsifier 4 — 215's own caveat, respected

`cqawesome` had 5 submissions inside a 2.5-day snapshot gap, so its doubling may have been one
arm or five. **No argument in this document requires it to be one.** Everything here is a
statement about the *ceiling* of the CD4 lever — +32.4% against +106% — which is unchanged
whether the leaders got there in one step or in five, and the same conclusion follows from
`jw_13288`'s independent ×2.071 window measured on a different account.

---

## 8. Verdict, in overall units

All priced on the exact six-member route of §6.0 (188's `nmae` identity, PLOLO rates), against
the incumbent `cd4f034`. The dose rows are also re-priced on the **realised board rates** in
§4.4 (`k_cos` 1.384, `k_pds` 1.987); **both calibrations pick `f = 0.34` as the argmax and
neither lets any dose clear the gate**, which is the strongest form this conclusion takes:

| candidate | what it is | overall vs `cd4f034` | gate | × gate | verdict |
|---|---|---:|---:|---:|---|
| **`refit-natural / lib1`** | **the operator refitted with CD4 in its basis, + the 1/k form — the best arm here, priced on all six** | **+0.00394** | +0.0089 (α ×1.29) | **0.44** | **VETO — and its permutation control is the arm itself (§5.4)** |
| `refit-natural / res f=0.34` | the refit operator alone, incumbent form | −0.00021 | +0.0019 | −0.11 | VETO |
| `refit-equal / lib1` | the refit at equal donor weight | −0.00370 | +0.0089 | −0.42 | VETO |
| **`f = 0.28`** | one dose step down | **−0.00024** | +0.0019 | −0.13 | **VETO** |
| **`f = 0.40`** | one dose step up — the largest the `pds` guard allows | **−0.00029** | +0.0019 | −0.15 | **VETO** |
| `f = 0.55` | past the guard | −0.00221 | +0.0019 | −1.16 | VETO |
| `f = 1.00` | the cosine argmax | **−0.00984** | +0.0019 | −5.18 | **VETO** |
| *the same, on realised board rates (§4.4)* | *`f` = 0.28 / 0.40 / 0.55 / 1.00 / 2.60* | *−0.00025 / −0.00010 / −0.00180 / −0.01195 / ≈−0.038* | *+0.0019* | *≤ −0.05* | ***VETO — argmax is again `f = 0.34`*** |
| `lib w=0.25 + res f=0.20 @adel` | best library arm on the `pds`+`mse` route | ≤ −0.00431 | +0.0089 (α ×1.62) | ≤ −0.48 | **VETO** — its exact `nmae` was not run; every measured library `nmae` is negative, so the bound holds a fortiori |
| `lib w=0.5 + res f=0.34 @adel` | best composed library arm priced on all six | −0.05098 | +0.0089 (α ×2.00) | −5.73 | VETO |
| `lib w=0.5 @adel` | library alone | −0.05913 | +0.0089 (α ×2.67) | −6.64 | VETO |
| **`lib w=1 @adel`** | **the 1/k route, CD4 as a fourth equal line** | **−0.13743** | +0.0089 (α ×3.9–4.2) | −15.44 | **VETO** |
| pseudobulk re-extraction | +7,658 genes, most generous law | ≤ +0.0011 | +0.0019 | ≤ 0.58 | **VETO — and negative at the measured elasticity** |
| low-expression pruning of the residual | by-product of §3.5 | ≈ +0.00024 | +0.0019 | 0.13 | **VETO** |

**Exactly one candidate is positive against the incumbent — the operator refit, at 0.44× its
gate and with a vacuous permutation control. Every other candidate is negative, and `cd4f034`
is the argmax of the dose grid and of all three forms.**

**NO-GO. Nothing clears. No slot. No builder invocation and no preregistration are written,
because nothing in this experiment should be built.**

**What the experiment establishes, and it is a well-supported negative:**

1. **The 2026-08-28 publication released no data.** The object is unchanged since 2026-05-28
   at 12 byte ranges; nothing on the bucket moved inside the window; the one refreshed
   artefact duplicates columns we already held. **The date coincidence in 215 §3.3 is a
   coincidence, and it can now be closed rather than carried.**
2. **The CD4 lever is exhausted at `f = 0.34`.** Its ceiling over 31 doses spanning `f` = 0 to
   40 is +32.4% cosine; the leaders moved +106%. The incumbent already sits at the largest
   dose the `pds` guard permits, and on the exact six-member price it is the argmax of the
   whole grid — every neighbouring dose is negative against it.
3. **The form question is settled with the confound removed.** The α re-solve 207 ran was on
   the wrong norm; corrected, the library family moves a long way and still loses by 21×, and
   at `w = 1` it is worse than having no CD4 at all.
4. **The one new mechanism found is about the OPERATOR, not about CD4.** The shipped whitener
   amplifies directions its own basis never measured by up to 7,943×, and that — not the
   residual form — is what forces every new source to enter after the whitener. Fixing it
   moves `pds` by **+9.85σ**, more than anything else here, but the fix is target-blind by
   algebra and therefore geometry, not direction. §5.4, §9 item 3b.
5. **Therefore the leaders' direction did not come from this atlas.** Not from its dose, not
   from its form, not from a better-processed version of it, and not from anything published
   on 2026-08-28. **The remaining direction hypothesis is elsewhere**, and §9 says where I
   would look.

---

## 9. What I could not resolve, and what I would hand on

1. **What the leaders actually did.** This experiment closes one hypothesis by measurement; it
   proposes no replacement. The `k_eff` arithmetic still says they hold ~10.6 effective lines
   against our 3.09 and 166 says there is no fifth public panel-covering CRISPRi screen — so
   the remaining explanations are an *estimator* change (194's territory) or a source outside
   the CRISPRi census (Tahoe-100M, 184; the observational atlases). **Neither is mine.**
2. **The pseudobulk object was never downloaded.** Its `.obs`/`.var` were read over 16.9 MB of
   range requests and its value was bounded from below by a measured elasticity rather than by
   building a bank from it. If a future agent wants the bank anyway, the bound to beat is
   **+0.0011 overall at the most generous law**, and the download is 41.5 GiB.
3. **`GWCD4i.DE_stats.by_donors.h5mu` (15.7 GiB) is untouched.** It gives per-donor-pair DE on
   the same 10,282-gene axis. It cannot add genes or targets, but it would allow a
   reliability-weighted or shrunk CD4 estimator — the `.obs` reliability columns
   (`donor_correlation_hits_mean`, `donor_correlation_all_min`, `guide_correlation_signif`,
   `n_guides`, `single_guide_estimate`) have been in the file since 2026-05-28 and **we have
   never used any of them.** That is one live thread this pass leaves, and it is an
   *estimator* question (194's territory), not a dose question.
3b. 🔴 **A fact about the OPERATOR that its owner should have, and that this pass found by
   accident.** The shipped whitener amplifies its own near-null directions by up to **7,943×**,
   and CD4 has real energy there because the three donor lines never measured those
   directions. **That is what destroys anything added BEFORE the whitener**: the library form
   needs α = 4.17 through the shipped operator and α = 1.27 through one that has seen CD4 at
   0.87% of the trace, and its `pds` moves +9.85σ between the two (§5.4). The arm itself is
   vetoed — its permutation control is the arm — but the *diagnosis* generalises to any future
   pre-whitening component, and it says the operator's blindness, not the residual form, is
   what forces every new source to enter post-hoc. **That belongs to 193/198, not to me**, and
   it is the sharpest thing in this document that is not a negative.
4. **The licence tension is untouched** (§3.7). Recorded for Ahmet.
5. **Why `cd4f034` beat its prediction by 3.0×, and why I would not trust the next one
   either.** Of its realised **+0.0130557**, **+0.0079 (60%) came from `nmae`, `fid`, `reach`
   and `jac`** — one member whose sign PLOLO's exact 188 identity got *wrong* (+0.001762
   realised against −0.000152 predicted) and three that no local harness can read at all
   (rule 9). The two members every local instrument *can* price, `pds` and `mse`, together
   realised **+0.00519** against a predicted +0.00453 — i.e. **the priceable half was accurate
   to 15% and the whole miss lives in the unpriceable half.** That is worth knowing in both
   directions: it is why this document's vetoes are stated on the priceable legs and why none
   of them is closer than 4× to its gate. **The open number for whoever prices the next arm in
   this family is why the DE members move favourably with dose** — 216 measured the ordering
   half rising monotonically to `f ≈ 0.55` (§6.0), which is consistent, but the yield halves
   are unmeasured and they are where the board's surplus came from.
6. **The `refit-natural` operators were not saved.** Eight 17,942² eigendecompositions were
   computed and discarded (2.5 GB each). If anyone revisits §5.4, budget ~7 min per operator
   and note that `op216.py` rebuilds them deterministically from the frozen parent plus β.

---

## 10. Files

| script | what it produces |
|---|---|
| `gate216.py` | Phase 0: folds rebuilt from scratch, the operator census, 16 known-answer gates → `work/gate216.{json,log}` |
| `regate216.py` | the four gates re-scored against the tag their sources publish (a15 vs a23) |
| `verify216.py` | 215's Phase-2 windows re-derived from the raw snapshots with 197's identity → `work/verify216.json` |
| `atlas216.py` | Phase 1: the 12 byte-range identity probes, the release inventory by date, the suppl-table diff, panel coverage → `work/atlas216.json` |
| `remote216.py` | `.obs`/`.var` of the pseudobulk and per-guide objects over HTTP range requests, 30.8 MB fetched → `work/remote216.json` |
| `energy216.py` | how much gated truth energy the CD4 gene axis reaches, and what the pseudobulk would add → `work/energy216.json` |
| `coverage216.py` | the gene-coverage elasticity, by cutting the residual's support from the low-expression end → `work/coverage216.json`, digested to `work/covanalyse216.json` |
| `sweep216.py` | the 31-dose extended sweep + 13 shuffled controls, 176 arms → `work/sweep216.json` |
| `analyse216.py` | the curve, its turnover, the bootstrap CIs, falsifiers 1–3 → `work/analyse216.json` |
| `price216.py` | the priced grid through 096/197's identity, `k_cos` re-derived → `work/price216.json` |
| `form216.py` | the library form at three α policies including the delivered-norm solve, 252 arms → `work/form216.json` |
| `formanalyse216.py` | the form comparison, priced, with 207's two known-answer rows → `work/formanalyse216.json` |
| `op216.py` | the operator re-fitted with CD4 in its basis, 2 β × 4 folds × 5 arms, 8 eigendecompositions → `work/op216.json` |
| `opanalyse216.py` | the refit priced on `pds` + `mse` against a same-fold shipped-operator baseline → `work/opanalyse216.json` |
| `opsix216.py` | 188's exact `nmae` and the ordering half for the surviving refit arm → `work/opsix216.json` (and `work/opsix216-BUGGY.*`, kept as the record of §5.5's process note) |
| `perminv216.py` | falsifier 1 applied to the refit, where the Gram makes it vacuous by algebra → `work/perminv216.json` |
| `six216.py` | 188's exact `nmae` identity and the ordering half of `reach`/`fid` at every candidate dose and form → `work/six216.json` |
| `final216.py` | the joined six-member price and the gate verdicts, 3 known-answer gates vs 207 §8.1 → `work/final216.json` |
| `vacuity216.py` | non-vacuity on h1-test (both forms) and on all four PLOLO folds at every dose → `work/vacuity216.json` |
| `keff216.py` | the 1/k arithmetic under 159's law, 5 known-answer gates → `work/keff216.json` |
| `board216.py` | **INBOX-2**: the ladder re-priced on `cd4f034`'s realised board rates; `k_cos` and `k_pds` measured, not assumed → `work/board216.json` |
| `ctx216.py` | **INBOX-3**: the same grid with the three challenge frames kept separate → `work/ctx216.json` |
| `ctxanalyse216.py` | the per-context maximisers, their spread, the transfer loss, and the per-context `pds` cost → `work/ctxanalyse216.json` |

**Downloaded (the mandate's one exception, `data/external/cd4_cell2026/`): 36.6 MB** —
`data_sharing_readme.md` and seven supplementary tables, plus 30.8 MB of HTTP range reads that
were never written to disk. **Against a 60 GB budget. No bulk object was fetched.**

No `vcc`, no `cell-eval2`, no build, no submission, no production file touched, nothing posted
anywhere, and no query carried anything about our method.
