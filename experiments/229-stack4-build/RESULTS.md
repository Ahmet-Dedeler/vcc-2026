# 229 — The registry is IN, COMMITTED and byte-inert. **The arm is NOT.** 225's construction is not expressible in a single-union build, and every construction that IS expressible dies on 216's exact `nmae` identity — the six-member composite is **−0.00219**, not the rate card's +0.00268.

**Date:** 2026-09-04 · deterministic · **no `vcc submit`, no `cell-eval2`, no model name, no
description.** One file changed outside this experiment directory,
`scripts/make_submission_pm_cd4.py`, **now committed as `0204168`** with its known-answer hashes
in the message (INBOX-1). Nothing under `data/`, `src/`, or any other experiment was written —
223's, 225's, 120's and 166's source tables were **read** and hash-verified, never written.

---

## 0. Verdict

```
REGISTRY:  LANDED AND COMMITTED (0204168).  All 3 mandate known-answer rows PASS, plus
           INBOX-2's inertness row and two new ones (KA-4, KA-5) at 0.0e+00.

ARM:       NOT BUILDABLE — do not spend the slot.
           Nothing a single-union builder can express clears the +0.0019 shape-only gate
           once `nmae` is MEASURED instead of rate-carded.

           six-member composite vs the shipped `h1union`, 3 cancer folds, g = 0.35:
             ARM6 (all six sources)          -0.00219   (rate card said +0.00268)
             ARM4 (225's four names)         -0.00589   (rate card said -0.00102)
             ARM3 (H1 pair + dld1_lowmoi)    +0.00067   (0.35x gate — under it)

           Two packages were built and fully verified before this was measured.  They are
           on disk, correct, and NOT RECOMMENDED:
             submissions/229-stack6-g035.vcc  7a3c2d3e414f8d9dfba0bd1c42f2b93994dedfda73cf207cea9b7aabdc78291c
             submissions/229-stack4-g035.vcc  efa5cf7b92d33ff0be23ed77ae6478e8ee82d9b4597e9149b334d5df8e0e1674
```

**Two independent findings, either of which alone would have stopped the arm.**

**FINDING 1 — 225's arm is not expressible.** 225 priced its stack as a **second residual riding
on top of `h1union`**: `base222 = cd4@0.34 + h1pair@0.28`, then `+ stack4@0.35`. That is **two
residuals at two doses.** The builder carries **one** extra residual at **one** dose, and the
authorised edit (a registry plus a selector) does not change that. Naming 225's four sources
**replaces** the H1 pair instead of joining it. Measured on 225's own instrument, reproducing
225's own rows at `0.0e+00`:

| construction | expressible? | PLOLO price vs `h1union`, `g = 0.35` (rate card) |
|---|---|---:|
| **`CUM4`** — `h1union` + stack4 as a **second** residual (225's actual arm) | **NO** | **+0.00367** (1.93× gate) |
| **`CUM1`** — `h1union` + `dld1_lowmoi` (223's arm, INBOX-3 §4's fallback) | **NO** | **+0.00228** (1.20× gate) |
| `ARM6` — all six in one union at one dose | yes | +0.00268 (1.41×) |
| `ARM5` — H1 pair + dld1 + jiang + igvf | yes | +0.00257 (1.35×) |
| `ARM3` — H1 pair + dld1 | yes | +0.00148 (0.78×) |
| `ARM4` — 225's four names, H1 pair dropped | yes | **−0.00102** (−0.53×) |
| `ARM1` — `dld1_lowmoi` alone | yes | **−0.00277** (−1.46×) |

**INBOX-3 §4's fallback fails for the same reason the four-source arm does.** "A single source, so
a union of one expresses it exactly" is not right: 223's and 225's `+0.00228` for `dld1_lowmoi` is
**marginal over `h1union`**, and a union of one *replaces* the H1 pair. Measured, `ARM1` is the
**worst arm in the table at −0.00277**, and it costs **−27.9 % of cosine on `kolf2`.**

**FINDING 2 — and it is the one that ends it.** `six229.py` runs 216's exact `nmae` identity
(188's split), the instrument 228 used yesterday to kill 224's arm. **Nine known-answer rows
against `216/work/six216.json` and `222/work/six222.json` reproduce at exactly `0.00e+00`.** On
that instrument:

| arm | `nmae wald89` | Δ vs incumbent | overall `nmae` | 188's `e = ΔM/ΔU` |
|---|---:|---:|---:|---:|
| `h1union g=0.28` (**incumbent**) | 0.98569 | — | — | — |
| **`stack6` g=0.35** | **1.00291** | **+0.01722** | **−0.004783** | **−0.033** |
| `stack4` g=0.35 | 1.00357 | +0.01788 | −0.004966 | −0.125 |
| `stack3` g=0.35 | 0.98842 | +0.00273 | −0.000758 | +0.160 |

Break-even efficiency is **0.5**; 224's arm died at **0.072**; **`ARM6` is at −0.033** — its
matched mass `M` does not merely fail to grow, it **falls** (ΔM = −31.7) while ΔU = +959.7. **The
arm moves MAGNITUDES, not signs**, which is the answer to INBOX-4 §2 and the opposite of what the
precedent predicted. All three `nmae` pools agree and all three cross **1.0** (`wald89` 1.00291,
`lfc89` 1.00099, `h89` 1.00997) — the arm becomes worse than the reference on the identity.

**"3/3 cosine, 3/3 `pds`" was true of 228's arm too.** `ARM6` also raises top-89 ordering
precision (+0.00115) and `pds` (+0.001755 on a second, independent scoring path). It does not
matter: `nmae` alone gives back **−0.004783** against a whole-arm price of **+0.00268**.

**Recommendation: leave the slot unspent, or spend it on something else.** `AGENTS.md` rule 7 —
a slot must buy information — and INBOX-3 §5 both point the same way.

---

## 1. The change — one file, one hardcoded tuple replaced by a registry

`scripts/make_submission_pm_cd4.py`: `71cf433d…` (222's) → **`30fd9ab34d1e540e3329919bced4d4cfcfc19c81d34213d00b67289fd64532d0`**,
committed as **`0204168`**, **168 insertions, 7 deletions**, all inside the `VCC_H1_UNION_G` block.

| added | what it does |
|---|---|
| `UNION_SOURCE_REGISTRY` | name → `{path, sha256: {delta.npy, b_ctrl.npy, genes.txt, perts.txt}}`. **Ten entries, 40 hashes.** Rule 11: a path and a sha, never a description. |
| `VCC_UNION_SOURCES` | ordered comma-separated list of active entries; unknown or repeated names raise at import. **Unset = `("gse295214_h1", "h1_2025")`**, so `VCC_H1_UNION_G=0.28` alone still reproduces the shipped package bit-for-bit. |
| `verify_union_source_hashes()` | called at the top of `prepare_h1_union`. Every selected bank's four files are re-hashed and compared **before a byte of a 25-minute build is spent**, and the build **raises** on a mismatch. |

Four of the ten banks live under `experiments/*/work/`, a tree other agents write to. Recording a
sha catches a changed file only if something compares it — this compares it.

**Order is part of the arm and is recorded.** `union_residual` accumulates a float64 sum per
coordinate, so a permutation of the same members can differ in the last bit. The env var's order
is the union's order and it is written into `config.h1_union_sources`.

**`solve_cd4_fraction` was not modified.** Not one line. §4 proves the consequence.

```bash
experiments/229-stack4-build/build229.sh submissions/229-ka-cd4f034-f034.h5ad       # KA row 1
VCC_H1_UNION_G=0.28 experiments/229-stack4-build/build229.sh \
  submissions/229-ka-h1union-g028.h5ad                                              # KA row 2
VCC_H1_UNION_G=0.35 \
VCC_UNION_SOURCES=gse295214_h1,h1_2025,dld1_lowmoi,jiang_5line_noK562,igvf_cardiomyocyte,myllia_a549 \
  experiments/229-stack4-build/build229.sh submissions/229-stack6-g035.h5ad
```

### 1.1 The registry — ten entries, every path and every `delta.npy` sha

| name | path | `delta.npy` sha256 | rows | panel now | pess. Oct |
|---|---|---|---:|---:|---:|
| `gse295214_h1` | `experiments/115-gse295214-h1/work/compact-donor` | `d0fa6555f9a62ce503d2ad55c8bfda954322657b52dcd48a24e29c330e718c7a` | 2,354 | 64 | 34 |
| `h1_2025` | `data/donors/h1_2025` | `269628f6eae749857e619bad4556ce9ac8b8810882e7af8a9925e3178eec0cc8` | 255 | 20 | 3 |
| `dld1_lowmoi` | `experiments/223-silent-zero-audit/work/tables/dld1_lowmoi` | `fd7bd0637cf9949971040f6c184b28f288391191f2e956c7cb16bb92d266aaee` | 4,715 | 67 | **73** |
| `jiang_5line_noK562` | `experiments/225-source-sweep/work/tables/jiang_5line_noK562` | `11c4e9ec149419ec6e49d3ea8334900da51b9ceb370042fcb9fdd286b52eed97` | 218 | 9 | 4 |
| `igvf_cardiomyocyte` | `experiments/225-source-sweep/work/tables/igvf_cardiomyocyte` | `d30b4930671c5f4b18c11ac88bb40b3c27a4d9c0a03875f48ba4ed2258736935` | 1,888 | 87 | 36 |
| `myllia_a549` | `experiments/120-myllia-a549/work/compact-donor` | `e81436089579914637bfcb54df454e5bf4b3759fba44b8d8e48e5204ed1c9ed6` | 175 | 14 | 4 |
| **`nadig_hepg2`** *(October, never selected)* | `experiments/166-donor-acquisition/donors/nadig_hepg2` | `ca4565e84045b91dbfca35280c01ae73f9a67112b3de2c8ba799cfda3fbce567` | — | **0** | 33 |
| **`nadig_jurkat`** *(October)* | `data/donors/nadig_jurkat` | `1fa16024c2d7115fc4b7c06f0ab7d9afe8ba6ebd0ae6dcef62c4955636e29129` | — | **0** | 33 |
| **`rpe1`** *(October)* | `data/donors/rpe1` | `253cc16e51f4635f23e8e597d91bd38a84162a1102486d4a5764c48c379c2b4a` | — | **0** | 32 |
| **`k562_essential`** *(October)* | `data/donors/k562_essential` | `4daf73cb932fdb089857c8629a70e83d8d3090570aa1f7d678c831b620ea281b` | — | **0** | 30 |

**INBOX-2 answered: all four October banks are already in the standard compact-donor layout**
(`delta.npy`, `b_ctrl.npy`, `genes.txt`, `perts.txt` all present, all four hashed). **None needs a
build step on Oct 22** — selecting one is an env var. The two H1 entries' four shas each reproduce
222's recorded values exactly; the full 40-hash table is in the source and in every arm build's
`effects.h1_union.sources[]`.

---

## 2. Known-answer row 1 — the default path is byte-identical to `cd4f034`

Rebuilt with the **changed** builder, both flags unset:

| | value |
|---|---|
| `shasum -a 256 submissions/229-ka-cd4f034-f034.h5ad` | **`7036e89f8a03c82d69e79a5c6690b275adae52dbb3ff72bb6ab23e7ca622bf7a`** |
| required | `7036e89f8a03c82d69e79a5c6690b275adae52dbb3ff72bb6ab23e7ca622bf7a` |
| bytes | **12,648,807,800** (required 12,648,807,800) |
| **verdict** | **PASS** |

Its `.meta.json` differs from `210-cd4f034.h5ad.meta.json` in **exactly three keys**, all
non-substantive: the builder's own implementation hash, `output.path`, `peak_rss_gib`
(19.833 → 19.598). `output.sha256`, `output.nnz`, every `effects.*`, every `operator.*`, every
input hash: **identical**.

## 3. Known-answer row 2 — the H1 pair at `g = 0.28` reproduces the shipped `h1union`

`VCC_H1_UNION_G=0.28` with `VCC_UNION_SOURCES` **unset**:

| | measured | required | verdict |
|---|---|---|---|
| `229-ka-h1union-g028.h5ad` sha256 | **`8df4dadd5a63b76c7ca3e3d57a511e7e52000b79d016d140897c55dd24c3ecf6`** | `8df4dadd…` | **PASS** |
| `.h5ad` bytes | **12,648,807,800** | 12,648,807,800 | **PASS** |
| **`229-ka-h1union-g028.vcc` sha256** | **`f13ce40bca908af7cb85b8e84887af473aae4588410c238529139b7d1f245b1e`** | **`f13ce40b…`** | **PASS** |
| **`.vcc` bytes** | **2,824,785,920** | **2,824,785,920** | **PASS** |

Its `.meta.json` differs from the shipped `222-h1union-g028.h5ad.meta.json` in **exactly three
keys** — implementation hash, `output.path`, `peak_rss_gib`. Every `config.*` including
`h1_union_sources = ['gse295214_h1', 'h1_2025']`, the whole `effects.h1_union.*` block and its
four-file source shas, every `operator.*`, `output.sha256`, `output.nnz`: **identical.**

### 3.1 The effect matrix, proved seven more times, in separate processes

`stage229.py` (222's `stage222.py`, copied) imports a builder **by path** and stages the
production effect matrix without emitting, so the pre-change and post-change builders run in
different processes on identical inputs. The pre-change control is the byte copy taken **before
the first edit** (`work/make_submission_pm_cd4.PRE229.py`, `71cf433d…`, verified equal to
`git show HEAD:` at git `d1fc9859` at that moment — the trap `63400de` records).

| tag | builder | `VCC_H1_UNION_G` | `VCC_UNION_SOURCES` | `effects.final_sha256` |
|---|---|---|---|---|
| `pre` | `PRE229.py` (`71cf433d…`) | unset | — | `87a4b6d181298f63826dc9caeb9674510a5dbc0d3c388e9e4fbcb359f3cd2de2` |
| `off` | 6-entry registry (`a8785cab…`) | unset | unset | **same** |
| `off_reg10` | **10**-entry registry (`30fd9ab3…`) | unset | unset | **same** |
| `pre_on028` | `PRE229.py` | 0.28 | — | `abeb0749f8ecc9410fcec0d9d0cebd2569c97ca456d73f9e6b3bfd359164c315` |
| `off_on028` | 6-entry | 0.28 | unset | **same** |
| `on028_reg10` | 10-entry | 0.28 | unset | **same** |
| `stack6g035` | 6-entry | 0.35 | the six | `c5721825c5cc1f3b3e3369b471a6209259acb744971c8ad9dbd1721214233c6a` |
| `stack6_reg10` | 10-entry | 0.35 | the six | **same** |
| `stack4g035` | 6-entry | 0.35 | the four | `0f45331f7c01eeeeb7878220a958ddc43e940ab868f90e774812921823bcbee1` |

`87a4b6d1…` is `210-cd4f034.h5ad.meta.json`'s own `effects.final_sha256`; `abeb0749…` is
`222-h1union-g028`'s; `c5721825…` and `0f45331f…` equal the two arm builds' own.

## 4. Known-answer row 3 — `solve_cd4_fraction` did not move by one bit

| package | `effects.cd4.fraction` | `target_f_requested` | `f_per_unit_fraction` | `fully_supported_targets` |
|---|---|---|---|---|
| `210-cd4f034` | 0.32207201067578145 | 0.34 | 1.0556645369046553 | 276 |
| `222-h1union-g028` (shipped) | 0.32207201067578145 | 0.34 | 1.0556645369046553 | 276 |
| **all six 229 packages** | **0.32207201067578145** | **0.34** | **1.0556645369046553** | **276** |

**The stop condition did not fire: the source-list change moves CD4's solved fraction by zero
bits.** Every arm here is one variable.

---

## 5. Three more known-answer rows

### 5.1 KA-4 — the production union IS 225's union, and KA-5 — 225's ladder reproduces

`price229.py`, 4 PLOLO folds, 327 s, peak **6.21 GiB**, **0 failures**.

For each of six selections the union `scripts/make_submission_pm_cd4.py` builds is compared to the
union 225 priced (`resid218.source_residual` + `compose218.union_rows` — an independently written
implementation):

| selection | folds | `max abs Δ mean_effect` | `max abs Δ votes` | carrying (k562 / hct116 / hek293t / kolf2) |
|---|---:|---:|---:|---|
| `h1pair` | 4 | **0.000e+00** | **0** | 74 / 76 / 76 / 71 |
| `stack1` (dld1) | 4 | **0.000e+00** | **0** | 60 / 66 / 65 / 62 |
| `stack3` | 4 | **0.000e+00** | **0** | 115 / 123 / 122 / 114 |
| `stack4` | 4 | **0.000e+00** | **0** | 130 / 135 / 135 / 126 |
| `stack5` | 4 | **0.000e+00** | **0** | 161 / 168 / 168 / 155 |
| `stack6` | 4 | **0.000e+00** | **0** | 165 / 172 / 172 / 159 |

**KA-5.** 225's `CUM1` and `CUM4` rows — the second-residual frame — reproduced through the
production union on 3 cancer folds × 6 doses × 2 arms against `225/work/stack225d.json`:
**Δ(Δcos) = 0.0e+00 and Δ(Δ`pds`) = 0.0e+00 on all 36 cells.** Baselines also match 225 §1 to the
last printed digit (`h1union` cos/pds 0.039436/0.805397, 0.032080/0.792213, 0.026157/0.763623).

### 5.2 INBOX-2's row — REGISTRY MEMBERSHIP CANNOT AFFECT A BUILD

The six-source package was rebuilt end-to-end with the four October banks present in the registry
and none of them selected:

| | 6-entry registry | 10-entry registry | verdict |
|---|---|---|---|
| `.h5ad` sha256 | `e2b9238723fea123d1130c88dbfcf7518d5a5fa44b9a27988abfc459de1936cd` | **identical** | **PASS** |
| bytes | 12,657,196,408 | **identical** | **PASS** |
| `.meta.json` differences | — | **3**: implementation hash, `output.path`, `peak_rss_gib` | **PASS** |

**Pre-registering a source we cannot yet price is safe.** Nothing reads an entry
`VCC_UNION_SOURCES` does not name, and `verify_union_source_hashes` only hashes the selected
banks. (The rebuilt `.h5ad` was deleted afterwards as a proven byte-duplicate; its `.meta.json` is
kept at `work/229-stack6-g035-reg10.h5ad.meta.json`.)

### 5.3 `six229.py`'s own nine known-answer rows

`nmae_wald89`, `top89_precision` and `delivered_norm_ratio` for three arms, against
`216/work/six216.json` and `222/work/six222.json`: **9/9 PASS at `|Δ| = 0.00e+00`.**
`a15 f=0.10` 0.98475707 / 0.58994688 / 1.00000000 · `cd4 f=0.34` 0.98530336 / 0.59751696 /
1.00094558 · `h1union g=0.28` 0.98569263 / 0.60312409 / 1.00333918.

---

## 6. The price — and the correction that ends the arm

Frame, functional, gate and tag, named: **162's PLOLO folds; macro per-fold cosine (`a23`); the
0.0019 shape-only gate; fold set `k562_gwps` + `orion_hct116` + `orion_hek293t`; `kolf2` reported
separately in §7.** **PLOLO units, never converted to board points** — the factor is not
identified (1.384 on `cd4f034`, 1.000 on `h1union`, matched `cos_agg`).

### 6.1 The full single-union lattice, marginal over the shipped `h1union` (rate card)

| arm | union | `g` | rel Δcos | Δ`pds` | Δ`q_w` | cos+ | pds+ | price | × gate | `kolf2` rel Δcos |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `ARM1` | dld1 | 0.35 | −5.504 % | −0.00177 | −0.00036 | 0/3 | 0/3 | −0.00277 | −1.46× | **−27.87 %** |
| `ARM3` | H1 pair + dld1 | 0.28 | +0.356 % | +0.00140 | +0.00002 | 3/3 | 3/3 | +0.00114 | 0.60× | −3.24 % |
| `ARM3` | " | 0.35 | +1.406 % | +0.00148 | −0.00001 | 3/3 | 3/3 | +0.00148 | 0.78× | +1.66 % |
| `ARM4` | 225's four | 0.35 | −2.593 % | −0.00032 | +0.00004 | 1/3 | 2/3 | −0.00102 | −0.53× | −18.45 % |
| `ARM5` | + jiang + igvf | 0.35 | +1.891 % | +0.00285 | +0.00028 | 3/3 | 3/3 | +0.00257 | 1.35× | +3.97 % |
| `ARM6` | all six | 0.28 | +1.790 % | +0.00259 | +0.00036 | 2/3 | 3/3 | +0.00233 | 1.23× | −1.25 % |
| **`ARM6`** | all six | **0.35** | **+2.865 %** | **+0.00269** | +0.00025 | **3/3** | **3/3** | **+0.00268** | **1.41×** | +3.81 % |
| `ARM6` | all six | 0.45 | +3.891 % | +0.00182 | +0.00019 | 3/3 | 2/3 | +0.00229 | 1.21× | +10.14 % |
| *`CUM1`* | *2nd residual* | *0.35* | *+1.405 %* | *+0.00256* | *+0.00016* | *3/3* | *3/3* | *+0.00228* | *1.20×* | *−0.23 %* |
| *`CUM4`* | *2nd residual* | *0.35* | *+4.037 %* | *+0.00361* | *+0.00055* | *3/3* | *3/3* | *+0.00367* | *1.93×* | *+8.02 %* |

`ARM6` at `g = 0.35` is the best expressible construction at every dose ≥ 0.28 and the only one
that beats 223's `dld1_lowmoi` recommendation. **`g = 0.35` was not tuned** — it is 225's own
pre-registered dose, and the whole ladder was priced before selection (`PREREGISTRATION.md`).

### 6.2 🔴 The six members with `nmae` MEASURED — INBOX-4's requirement, and the verdict

The rate card's `nmae` leg is not a measurement: it is `0.000980 × (Δcos / 0.01)`, so it is
**positive by construction whenever cosine rises.** 216's exact identity measures it directly.

| member | `ARM6` rate card | **`ARM6` measured** | `ARM4` measured | `ARM3` measured |
|---|---:|---:|---:|---:|
| `mse` | +0.000306 | +0.000306 | −0.000353 | +0.000171 |
| `reach` | +0.000206 | +0.000206 | −0.000237 | +0.000115 |
| `fid` | +0.000071 | +0.000071 | −0.000082 | +0.000040 |
| `jac` | +0.000008 | +0.000008 | −0.000010 | +0.000005 |
| `pds` | +0.002004 | +0.002004 | −0.000239 | +0.001100 |
| **`nmae`** | **+0.000083** | **🔴 −0.004783** | **−0.004966** | **−0.000758** |
| **total** | **+0.002678** | **−0.002188** | **−0.005887** | **+0.000673** |
| **× gate** | 1.41× | **−1.15×** | −3.10× | 0.35× |

**The measured `nmae` leg is the opposite sign and 58× the magnitude of the rate card's.** Nothing
expressible clears the gate.

### 6.3 Signs or magnitudes — 188's split, the statistic that killed 224's arm

| arm | ΔU | ΔM | **`e` = ΔM/ΔU** | verdict |
|---|---:|---:|---:|---|
| break-even | — | — | **0.50** | — |
| *224's arm (228 measured)* | — | — | *0.072* | *died* |
| **`stack6` g=0.35** | **+959.66** | **−31.70** | **−0.033** | **MAGNITUDES** |
| `stack4` g=0.35 | +827.56 | −103.78 | −0.125 | MAGNITUDES |
| `stack3` g=0.35 | +257.84 | +41.18 | +0.160 | MAGNITUDES |

**`ARM6` is worse than 224's arm on the statistic that killed 224's arm.** Its matched mass does
not merely fail to grow — it **falls**. And this is not an amplitude artefact: the delivered norm
ratio *drops* (−0.0028 vs the incumbent), so the +960 of unmatched magnitude is a redistribution
onto coordinates whose signs do not match truth, not a global inflation.

All three `nmae` pools agree and all three cross **1.0** — the arm is worse than the reference on
the identity: `wald89` 0.98569 → **1.00291**, `lfc89` 0.99244 → **1.00099**, `h89` 0.98203 →
**1.00997**.

**One correction to INBOX-4's premise, because it changes the reading.** INBOX-4 says *"`h1union`
… its `nmae` went UP +0.00588"*. On this instrument `h1union`'s **`nmae` went slightly DOWN**
(raw **+0.00039** worse than `cd4f034`, overall −0.000108 — 222 §8 reports exactly this). What
went up by ~+0.0056 was **top-89 ordering precision** (0.59752 → 0.60312). So the precedent is not
"a donor addition raises `nmae`"; it is "a donor addition can raise ordering precision while
leaving `nmae` flat". **`ARM6` reproduces the precision half (+0.00115) and breaks the `nmae`
half by 44×.**

### 6.4 How far outside the instrument's calibrated range this is

222 §8 records that this instrument **fails its own family's board retrodiction** — for the CD4
dose step it gives −0.000152 overall where the board delivered +0.001762, wrong in sign and 11.6×
in magnitude. **So read `nmae` here as a magnitude bound, not a sign.** As a magnitude bound it is
still decisive:

| step | raw Δ`nmae` | this instrument's overall | board realised |
|---|---:|---:|---|
| `a15` → `cd4f034` | +0.00055 | −0.000152 | **+0.001762** (sign inverted) |
| `cd4f034` → `h1union` | +0.00039 | −0.000108 | part of **+0.0069114** total |
| **`h1union` → `stack6`** | **+0.01722** | **−0.004783** | **unknown — 31–44× outside the calibrated range** |

Both calibration points sit at ~+0.0004–0.0006 raw. `ARM6` is at **+0.01722**. Extrapolating a
sign-unreliable instrument 31–44× past every point at which it has ever been checked, in the
direction it says is bad, with an independently damning mechanism (`e` = −0.033), is not a basis
for spending a slot.

### 6.5 Per fold, and where each arm sits relative to 223's silent zero

| fold | `ARM6` rel Δcos | Δ`pds` | fraction | carrying |
|---|---:|---:|---:|---|
| `k562_gwps` | +0.611 % | +0.002395 | 0.26506 | 165/269 |
| `orion_hct116` | +3.321 % | +0.001491 | 0.26524 | 172/290 |
| `orion_hek293t` | +4.662 % | +0.004175 | 0.26881 | 172/294 |

`k562_gwps` is `ARM6`'s weak fold and `dld1_lowmoi`'s strongest (225: +1.646 % for DLD-1 alone) —
one shared dose dilutes DLD-1 exactly where it is best, which is the mechanism behind the
`CUM4` → `ARM6` gap.

| arm | carrying /300 | coverage | median carrying | median all-rows | fraction | **OLD all-rows rule would solve** |
|---|---:|---:|---:|---:|---:|---:|
| `h1union` (shipped) | 77 | 25.7 % | 6.1654 | 0.0000 | 0.19207 | **0.00000** |
| `ARM4` @0.35 | 139 | 46.3 % | 7.1468 | **0.0000** | 0.20712 | **0.00000 — identically absent** |
| `ARM6` @0.35 | **176** | **58.7 %** | 6.4400 | 4.5282 | 0.22986 | 0.32690 (**1.42×** over-dose) |

`ARM6` is the first arm this project has built that **crosses 50 % coverage**, so both sizing rules
return something and they disagree by 42 %. 222's coverage-aware rule is what makes the number
0.22986 rather than 0.32690.

---

## 7. `kolf2`, reported separately

| arm | `kolf2` rel Δcos | `kolf2` Δ`pds` | 3-cancer-fold rel Δcos | ratio |
|---|---:|---:|---:|---:|
| `ARM1` (dld1 alone) @0.35 | **−27.87 %** | −0.02347 | −5.504 % | 5.1× worse |
| `ARM4` @0.35 | **−18.45 %** | −0.01370 | −2.593 % | 7.1× worse |
| `ARM6` @0.35 | **+3.81 %** | +0.00400 | +2.865 % | **1.33×** |
| *`CUM4` @0.35* | *+8.02 %* | *+0.00491* | *+4.037 %* | *1.99×* |

**`ARM6` does not inflate on `kolf2`** — 1.33× against the 3–5× warned about and `CUM4`'s 2.0×.
The mechanism is the cleanest evidence in this experiment that **the H1 pair is load-bearing**:
`kolf2` is pluripotent, the H1 pair is the only pluripotent thing in the union, and it sits in
**both** the baseline and `ARM6`, so its advantage cancels. Remove it — `ARM4`, or `ARM1` — and
`kolf2` collapses by **−18 % to −28 % of cosine.** Any future arm that drops the H1 pair must
carry this number.

---

## 8. The five falsifiers — all pass on the built packages

The packages are correct; it is the arm that is not worth a slot. Every falsifier was run and
passed *before* §6.2 was measured.

### 8.1 Falsifier 1 — the config diff is exactly the source-list change and its consequences

`metadiff229.py` (222's, plus an `--against h1union` mode). **Known-answer row for the tool
itself:** re-run on 222's own packages it reproduces 222's falsifier-1 PASS unchanged.

**vs `cd4f034`:** config ADDED exactly `['h1_union_g', 'h1_union_sizing_rule',
'h1_union_sources']`, 0 changed, 0 removed; audit CHANGED exactly
`['effects.cd4.post_add_preserve.candidate_pre_preserve_median_norm',
'effects.cd4.post_add_preserve.post_add_global_preserve_scale', 'effects.final_sha256']`; 109
keys added, all under `effects.h1_union.`; 0 unexpected; **0 frozen-input differences**; 1
implementation hash; 0 unclassified.

**vs `h1union`:** config CHANGED exactly `['h1_union_g', 'h1_union_sources']`, 0 added, 0 removed;
**0 stray audit changes** outside `effects.h1_union.` and the two shared preserve scalars; 0
frozen-input differences; 1 implementation hash; 0 unclassified.

**54 of 54 pinned scalars byte-identical on both arms against both incumbents** (zero `FAIL` lines
across all four diffs), including `operator.sha256`
**`e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123`**, `whiten_gamma` **0.65**,
`whiten_lambda` **0.0**, `alpha` **0.982**, `cd4_target_f` **0.34**, `emitter` **`rc_callset`**,
`whitener_id` `union_pruned`, `effects.post_whiten_sha256`,
`effects.whitening.after_all_median_norm` 4.229295604454299, `frozen_bundle.sha256`,
`output.n_rows` 360,000. **PASS.**

### 8.2 Falsifier 2 — the digest test, and the toll it implies

Per-gene panel-total digest on all 900 (context, target) blocks: **0/900 identical** for
`stack6` vs `cd4f034`, `stack6` vs `h1union`, and `stack4` vs `cd4f034`. Nonzero-structure
digests: 0/900.

**The delivered profile changes on every block; the emission LAW does not.** Rule 13 is
single-valued and the classification sets the toll:

| rule 13's class | this arm? |
|---|---|
| emission-law, profile-CHANGING → **−0.0146 scaled** | **no** — `config.emitter` is `rc_callset` in both, and every `rc_callset_*` / `rc_key` / `emission_seed` constant is byte-identical |
| emission-law, profile-PRESERVING → 0 | n/a |
| **profile-class** (operator, donors, magnitude, CD4 dose) → **0** — 169 falsified the toll at +0.0006 | **yes.** A donor/residual arm. |

**Toll = 0.** Supporting evidence the emission law is untouched: per-cell library totals identical
in range (17,695–22,323, `sd median` 0.398 → 0.400) and `n_pred` unchanged (§8.5).

### 8.3 Falsifier 3 — non-vacuity. There is no zero.

| comparison | ‖ΔE‖/‖E‖ | **per-block minimum** | blocks exactly 0 | blocks < 1e−6 |
|---|---:|---:|---:|---:|
| `stack6` vs `cd4f034` | 0.259663 | **0.055972** | **0 / 900** | 0 / 900 |
| `stack6` vs `h1union` | 0.234744 | **0.044751** | **0 / 900** | 0 / 900 |
| `stack4` vs `cd4f034` | 0.232363 | **0.043691** | **0 / 900** | 0 / 900 |

Per-block (`stack6` vs `cd4f034`): p01 0.056624, median 0.185608, max 1.623789; per context A
0.260395 / B 0.260714 / C 0.257914. **Every delivered block moved. PASS.**

### 8.4 Falsifier 4 — carrying vs non-carrying separation

**At the effect-matrix level the statement is exact.** `effects-stack6g035.npy` vs
`effects-off.npy`, 300 rows, 176 carried:

| row class | n | rel. displacement | after dividing out the shared scalar **0.944814213** | cosine(arm, inc) |
|---|---:|---|---|---|
| **carrying** | 176 | min 0.056105 · median 0.313764 · max 1.378387 | min **1.02e−02** · median 3.11e−01 | 0.6475 – 0.9999 |
| **NOT carrying** | 124 | **exactly 0.055186 for all 124** | max **9.29e−17** | **1.0000000000** exactly |

`0.055186 = |1 − 0.944814|`, and the scalar's spread across all 124 non-carrying rows is **zero to
9 decimal places**. **Separation 1.1 × 10¹⁴.** Against `h1union`: scalar 0.956245406, separation
9.9 × 10¹³. For `stack4` vs `cd4f034`: 0.957131074 and 8.9 × 10¹³.

**At the delivered-count level the emitter blurs it, more than in 222:**

| block class | n | cosine min / median | `‖ΔE − s·E‖/‖E‖` min / median / max |
|---|---:|---|---|
| **carrying** | 528 | 0.62528515 / **0.95263490** | 0.020448 / **0.311104** / 1.640143 |
| **NOT carrying** | 372 | **0.99926525** / 0.99971083 | 0.009212 / **0.022732** / 0.036241 |

Median separation **13.7×**, cosines separate cleanly, but the **min-carrying / max-non-carrying
test OVERLAPS at 0.6×** (222 got 3.3×). The reason is measured: 222's shared preserve scalar moved
**1.2 %**, this one moves **5.5 %**, so `rc_callset` re-draws integers under 4× the amplitude
change while the thinnest carrying blocks move only 2.0 %. **PASS at the level where the claim is
exact (1e−16); REPORTED NOT PASSED in the delivered min/max form.**

### 8.5 Falsifier 5 — `vcc prep --dry-run` clean, and the `.vcc` pre-prepped

```
✓ dry run — no .vcc written
  cells: 360000   genes: 18533 (float32)
  perturbations: column 'target_gene' → 'target_gene'
  contexts: A=120000, B=120000, C=120000
  targets: verified against the official list
  normalization: counts-preserved
```

| cap | limit | `stack6` | margin |
|---|---:|---:|---|
| stored entries in `.X` | 4,750,000,000 | **1,578,196,231** | 33.2 % of cap |
| counts per cell | 1,000,000 | **22,323** (min 17,695) | 44.8× under |
| rows | 360,000 | 360,000 | exact |
| cells per (context, target) block | 400 | 400 on all 900 | exact |

Both arms' `.vcc` are written and on disk. Bonus (222's falsifier 4): challenge `n_pred` pooled
median **2,599** for both `stack6` and `h1union` — **+0.00 %** against ±2 % — and delivered order
ρ(−log p, |H|) inside the call set **+0.9257** against the incumbent's +0.9306.

### 8.6 F1 and F3 — the source controls, on what is BUILT

`falsif229.py`, 48 arms, 87 s, marginal over `h1union`, three cancer folds. F1 is 218's derangement
(geometry, per-row norms, support pattern and carrying count preserved; only which target each row
belongs to moves), two seeds. F3 is CD4's own information on exactly the union's rows at the same
dose under the same rule, scored as a **difference** (225's correction to 223's ratio convention).

| arm | `g` | real | SHUF `20726` | SHUF `771103` | row-matched CD4 | **F1** | **F3** |
|---|---:|---:|---:|---:|---:|---|---|
| `ARM6` | 0.35 | **+2.865 %** | −10.267 % | −10.017 % | −3.923 % | **PASS**, 12.9 pp | **PASS** |
| `ARM6` | 0.28 | +1.790 % | −9.216 % | −8.944 % | −4.114 % | **PASS**, 11.0 pp | **PASS** |
| `ARM4` | 0.35 | −2.593 % | −10.744 % | −8.669 % | −5.220 % | beats controls, arm negative | — |

`ARM6`'s margin over its own shuffle is **65× 204's ±0.20 % resolution.** **The cosine gain is
real.** It is simply not worth what `nmae` charges for it — which is the whole lesson.

---

## 9. Concentration — leave-top-3-out on cosine and `pds` separately

Three cancer folds, `g = 0.35`, against **both** baselines.

| leave-top-k-out | `ARM6` vs `h1union`: Δcos | Δ`pds` | `ARM6` vs `cd4f034`: Δcos | Δ`pds` |
|---|---:|---:|---:|---:|
| full | 100 % (+0.000842) | 100 % (+0.002687) | 100 % (+0.003199) | 100 % (+0.007574) |
| k = 1 | 70.2 % | 26.7 % | 90.7 % | 66.8 % |
| **k = 3** | **29.1 %** | **−40.9 %** | **76.4 %** | **27.7 %** |
| k = 5 | 7.0 % | −89.3 % | 67.1 % | 0.0 % |
| k = 10 | −24.9 % | −171.9 % | 49.2 % | −46.2 % |

**The MARGINAL gain over the incumbent is more concentrated than 222's was** (222: 72.7 % cos /
14.6 % `pds`). Against `cd4f034` the picture is almost exactly 222's — **the six-source arm's
TOTAL gain is as broad as the shipped one's; it is the last increment that is narrow.** The board
realised ~**100 %** of `h1union`'s cosine prediction and **70 %** of its `pds` one, and `ARM6`'s
rate-card price was **75 % `pds`-carried**, so the retention discount lands on the larger half.
This was already a caution before §6.2 turned it moot.

**Spillover:** 99.3 % of Δcos and 99.8 % of Δ`pds` sit on the rows the union carries.

---

## 10. October — each source's projected coverage under `DEF-RUNBOOK` §5's two draws

`octo229.py`. On 2026-10-22 the panel is replaced with 300 targets **disjoint** from A/B/C's
(measured: A/B/C ∩ either draw = **0**) and only the last entry after that counts.

| source | bank rows | A/B/C (now) | realistic draw | pessimistic draw | retained |
|---|---:|---:|---:|---:|---:|
| `gse295214_h1` | 2,354 | 64 | 57 | 34 | 53 % |
| `h1_2025` | 255 | 20 | 7 | 3 | **15 %** |
| `dld1_lowmoi` | 4,715 | 67 | 67 | **73** | **109 %** |
| `jiang_5line_noK562` | 218 | 9 | 3 | 4 | 44 % |
| `igvf_cardiomyocyte` | 1,888 | 87 | 60 | 36 | 41 % |
| `myllia_a549` | 175 | 14 | 5 | 4 | 29 % |
| `nadig_hepg2` *(registered, unpriceable today)* | — | **0** | — | **33** | ∞ |
| `nadig_jurkat` *(registered)* | — | **0** | — | **33** | ∞ |
| `rpe1` *(registered)* | — | **0** | — | **32** | ∞ |
| `k562_essential` *(registered)* | — | **0** | — | **30** | ∞ |
| **UNION `h1pair`** (incumbent) | | **77** | **59** | **36** | **47 %** |
| **UNION `stack4`** | | **139** | **120** | **101** | **73 %** |
| **UNION `stack6`** | | **176** | **162** | **120** | **68 %** |

The incumbent's union loses **53 %** of its coverage at the pessimistic draw; `stack6` loses
**32 %** and still carries **120/300 — more than the incumbent carries today.** `dld1_lowmoi` is
the only registered source that **gains** (67 → 73), because it is a 4,715-row genome-scale bank;
`h1_2025` (255 rows, 20 → 3) is nearly worthless in October and is in the union only because it is
free and it is what makes known-answer row 2 exact.

**This does NOT rescue the arm.** October coverage is a statement about how many rows a source can
carry, not about whether carrying them helps; §6.2 says on today's panel it does not. But it is
the reason the registry is worth landing even with the arm dead, and the four pre-registered banks
are the reason Oct 22 is an env var rather than four bank builds under time pressure.

Two October items, sized and not retired:
- At **120/300 = 40 %** coverage a six-source union falls **below** the 50 % threshold again, so
  production's all-rows median would return exactly 0 and the residual would be *identically
  absent* under the pre-222 rule. **`DEF-RUNBOOK` §5.-1's stale-rehearsal item now has to cover a
  ten-entry registry**, not the H1 pair it was written for.
- **CD4 still over-doses ~1.24×** at the pessimistic draw. Re-measured here on all four folds and
  unchanged: 1.0120–1.0171× today, **1.2369–1.2453×** at the `q0.14` draw.

---

## 11. What is left on the table, named precisely

`CUM4` prices **+0.00367** on the rate card against `ARM6`'s **+0.00268**; the difference is that
`CUM4` sizes the H1 pair at `g = 0.28` on **its own** carrying rows and the four new sources at
`g = 0.35` on **theirs**, as two residuals. The minimal edit that would express it:

- `VCC_H1_UNION_G` + `VCC_UNION_SOURCES` become **group-wise** — e.g.
  `VCC_UNION_GROUPS="gse295214_h1,h1_2025@0.28; dld1_lowmoi,jiang_5line_noK562,igvf_cardiomyocyte,myllia_a549@0.35"`
  — each group unioned and sized independently by `solve_carrying_row_fraction`, then all groups
  **and** CD4 merged into the **same single add** with production's one global preserve rescale.
  `compose218.merged()` already chains, so this is an accumulation over groups and **not** a second
  add — the property falsifier 2 rests on survives.
- Known-answer rows survive unchanged: one group holding the H1 pair at 0.28 must still reproduce
  `f13ce40b…`, and the flag-off path must still reproduce `7036e89f…`.

🔴 **But do NOT build it on the strength of +0.00367.** That number is a rate card, and this
experiment's central result is that the rate card's `nmae` leg is **positive by construction** and
was wrong by **58×** on the first arm it was checked against. **Anyone who opens the group-wise
dose must run `six229.py` on `CUM4` FIRST** — it is 6 minutes and it is the difference between an
arm and a regression. `CUM4` carries the same four magnitude-heavy sources at a *higher* effective
dose than `ARM6` does, so the prior should be that it is worse on `nmae`, not better.

---

## 12. The packages

| artifact | bytes | sha256 |
|---|---:|---|
| `submissions/229-stack6-g035.vcc` — **NOT RECOMMENDED** | 2,825,646,080 | `7a3c2d3e414f8d9dfba0bd1c42f2b93994dedfda73cf207cea9b7aabdc78291c` |
| `submissions/229-stack6-g035.h5ad` | 12,657,196,408 | `e2b9238723fea123d1130c88dbfcf7518d5a5fa44b9a27988abfc459de1936cd` |
| `submissions/229-stack4-g035.vcc` — **NOT RECOMMENDED** | 2,825,809,920 | `efa5cf7b92d33ff0be23ed77ae6478e8ee82d9b4597e9149b334d5df8e0e1674` |
| `submissions/229-stack4-g035.h5ad` | 12,657,196,408 | `cf9951a52beeb6744492ab10326a57f17e20d6ff98e0c1c3cdffe7b86f583805` |
| `submissions/229-ka-h1union-g028.vcc` (KA row 2; proven duplicate of `222-h1union-g028.vcc`, **deletable**) | 2,824,785,920 | `f13ce40bca908af7cb85b8e84887af473aae4588410c238529139b7d1f245b1e` |
| `submissions/229-ka-h1union-g028.h5ad` (**deletable**) | 12,648,807,800 | `8df4dadd5a63b76c7ca3e3d57a511e7e52000b79d016d140897c55dd24c3ecf6` |
| `submissions/229-ka-cd4f034-f034.h5ad` (KA row 1; proven duplicate of `210-cd4f034.h5ad`, **deletable**) | 12,648,807,800 | `7036e89f8a03c82d69e79a5c6690b275adae52dbb3ff72bb6ab23e7ca622bf7a` |

Deleting the three known-answer artifacts frees **40.9 GB**; nothing downstream reads them. The
two arm packages are **53.9 GB** and can also go once the coordinator accepts the verdict.
**Nothing was deleted here** except the reg10 rebuild (a proven byte-duplicate, §5.2).

**No model name chosen, no description written, nothing submitted.**

## 13. Files

| script | what it does | runtime / peak RSS |
|---|---|---|
| `stage229.py` | 222's `stage222.py` — stages the production effect matrix from a builder given **by path**, without emitting | 9–15 s / 19.9 GiB |
| `price229.py` | KA-4, KA-5, the full single-union lattice, `kolf2`, the concentration arrays | 327 s / 6.21 GiB |
| **`six229.py`** | **216's exact `nmae` identity + 188's `e = ΔM/ΔU` split on six arms — the result** | ~3 min |
| `falsif229.py` | F1 (two seeds) and F3 on the built constructions | 87 s / 7.21 GiB |
| `octo229.py` | the October coverage table against 199's two synthetic draws | 2 s |
| `metadiff229.py` | 222's falsifier-1 diff plus an `--against h1union` mode | 1 s |
| `profile229.py` | 222's `profile222.py` with the union as an argument — falsifiers 2, 3, 4 | 60 s / 3.6 GiB |
| `probe229.py` | 222's `probe222.py`, unmodified — challenge `n_pred` and delivered order | 60 s / 13.2 GiB |
| `build229.sh` · `runbuilds229.sh` · `runarm6.sh` · `runreg10.sh` · `runpreps229.sh` | five builds and three preps | 23–30 min per build |
| `PREREGISTRATION.md` | the arm and its falsifier, with an honest note on when it was written | — |
| `work/make_submission_pm_cd4.PRE229.py` | **the git control**, `71cf433d…`, taken before the first edit | — |

## 14. Honest limits — and two operational corrections the runbook needs

1. 🔴 **Two builds ran CONCURRENTLY for 27 minutes.** `runarm6.sh` guards with
   `pgrep -f "make_submission_pm_cd4.py --whitener"` on a 60 s poll and slipped into the 9-second
   gap between one build's python exiting and the next starting. Peaks 12.79 + 14.73 GiB ≈
   27.5 GiB combined on a 64 GB machine also running four other agents. **No result was
   corrupted** — the overlapping build reproduced `8df4dadd…` exactly — and no single process
   exceeded 19.9 GiB, but the intent was strict serialisation. **Use a PID file, not a `pgrep`
   poll.**
2. 🔴 **`vcc prep` peaks at 24.9–26.8 GiB, not 22.33 GiB.** Measured on three preps here
   (24.92 / **26.77** / 24.84 GiB) against 222's recorded 22.33. `DEF-RUNBOOK` §6 budgets Oct 22
   around that number and it is **20 % low**. Nothing else was large during these preps.
3. **The `nmae` instrument's SIGN is not trustworthy** (222 §8: it mis-signs its own family's
   board point by 11.6×). §6.4 argues the magnitude is decisive anyway, 31–44× outside its
   calibrated range with an independent mechanism. **A cloud `cell-eval2` run (221's work) would
   settle this properly and nothing here substitutes for it.**
4. **`ARM5` and `ARM1` were not put through `six229.py`** — only `stack3`, `stack4`, `stack6`.
   `ARM5` is 0.94× of `ARM6`'s carrying set and `ARM1` is already negative on the rate card, so
   neither could clear once `ARM6` and `ARM3` bracket the range, but this is an inference, not a
   measurement.
5. **The delivered amplitude moves −2.6 %** (`A = mean|row|²` 14.453 → 14.082 vs `h1union`;
   222's arm moved +0.3 %). 187 found α = 0.982 sits on context A's ρ cliff. It was not
   re-calibrated — that would have been a second variable.
6. **No `cell-eval2` run and no board-side member measured.** `fid`, `jac` and `reach` carry
   `n_real` in a denominator or a union and are unreadable on every harness we hold (rule 9).
   `mse` is not measured here at all.
7. **`PREREGISTRATION.md` post-dates the builds.** INBOX-3 asked for it beforehand; it arrived
   after. The guard's substance holds — the selection instrument is 225's PLOLO harness, which
   reproduces 225 at 0.0e+00, and the built packages were **never scored** — but the file was not
   written first and a reader should discount accordingly.
