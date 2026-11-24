# 248 — The second union residual is IN the builder, committed and byte-inert (`2632b92`). The rider it expresses — `dld1_lowmoi @0.35` on top of the shipped stack — reproduces 225/229's `CUM1` row through the builder's own code to 0.0e+00 and **does NOT clear the mandate's gate**: +0.00216 under 243's class rule (1.14× the shape-only gate, below the required 1.5×) and a `pds` gain that keeps 25.2% leave-top-3-out (bar 50%). **No arm package was built. RIDER NOT BUILDABLE (verdict line at the bottom, §11).**

**Date 2026-09-04.** Deterministic. **No `vcc submit`, no `cell-eval2`, no model name, no
description.** One file changed outside this directory, `scripts/make_submission_pm_cd4.py`,
committed as **`2632b92`** (+239 / −40, file sha256 `f7a12a70…`) with its known-answer hashes in the
message (`work/commit-message.txt`); the package-level row landed after that commit and is in this file and the evidence commit. Nothing under `submissions/` carries a `248-` arm package; the
only `248-` artifact is the known-answer rebuild of the incumbent in §2.2 (deleted on match, meta
kept). Every ≥ 15 GiB stage ran alone under the mandate's gate (`memory_pressure` ≥ 40% free AND
`vm_stat` free+inactive ≥ 30 GiB AND zero other builders AND zero `vcc prep`, the working-tree
builder pinned to this experiment's bytes); peak RSS per stage in §9.

| row | required | measured | verdict |
|---|---|---|---|
| **KA 1** flags unset, staging | `effects_array_sha256` `abeb0749…` (222 `stage-on028.json`) | **`abeb0749f8ecc9410fcec0d9d0cebd2569c97ca456d73f9e6b3bfd359164c315`** through the changed builder; **and** through the pre-248 byte copy with `VCC_UNION2_*` set (the flag does not exist there) | **PASS** |
| **KA 1** flags unset, package | `.h5ad` `8df4dadd…` on 12,648,807,800 B | **`8df4dadd5a63b76c7ca3e3d57a511e7e52000b79d016d140897c55dd24c3ecf6` on 12,648,807,800 B** through `build248.sh` with both flags unset (1477.9 s / 16.72 GiB; `def_postbuild` 32/32, `moved=[]`; no `union2` key anywhere in the meta) | **PASS** |
| **KA 2** order invariance | union 1 UNSET, union 2 = the H1 pair @0.28 → the shipped staged effects bit for bit | **`abeb0749…`**, `union2.solved_fraction` 0.19207161964444322 = the shipped `h1_union` fraction, CD4 0.32207201067578145; at unit level union 2 alone == union 1 alone bit for bit (array + audits) | **PASS** |
| **KA 3** the rider reproduces 225/229's `CUM1` | rel macro cos +1.405%, Δ`pds` +0.00256 at 0.35, 3 cancer folds, through 162's instrument on the builder's own staged effects | **+1.405% / +0.002557**, 3/3 folds; 225's `stack225d` rows (real, SHUF, row-matched CD4) reproduced at every one of its six doses to |Δ| < 1e−12 (**120/120** rows over 4 folds); base == 218 at 0.0e+00; builder output == `compose218.merged()` chained, bit for bit | **PASS** |
| price (i) 243's class rule | ≥ 1.5 × 0.0019 = **+0.00285** | **+0.00216** (point class; best class +0.00309) | **FAIL** |
| `pds` ≥ 0 on 3/3 cancer folds | 3/3 | 3/3 (+0.004582 / +0.001555 / +0.001532) | PASS |
| leave-top-3-out retention on `pds` | ≥ 50% | **25.2%** (top-3 share 81.9%); cosine 54.5% | **FAIL** |
| price (ii) with 216's exact `nmae` | reported separately | **+0.00107** (`nmae` identity −0.000856, `e` = 0.18: magnitudes) | — |

---

## 1. The change — `scripts/make_submission_pm_cd4.py`, commit `2632b92`

`patch248.py` applies 17 exact-anchor replacements to git HEAD `f59b653b` (builder `2d63b979…`,
the byte copy `work/make_submission_pm_cd4.PRE.py`) and must reproduce the working tree byte for
byte (`work/make_submission_pm_cd4.248.py`, `work/builder.intended.sha256` = `f7a12a70…`); every
gate in this experiment re-hashed the working tree against that value before running.

- **`VCC_UNION2_SOURCES`** (ordered registry names) + **`VCC_UNION2_G`** (dose), both unset by
  default. Import-time rules, each observed firing (`envcheck248.py`, **11/11**): both or neither
  (a dose without sources, or sources without a dose, is a knob that is set and cannot act — 242
  #17 — and is refused); names must be in `UNION_SOURCE_REGISTRY`, may not repeat, may not be an
  empty list, and may not also be in union 1 while `VCC_H1_UNION_G` is set (a source in both
  unions would be dosed twice); `0 < g < 5`. The H1 pair as union 2 with union 1 unset is
  allowed (that is KA 2).
- **`prepare_h1_union(..., sources=None)`**, **`verify_union_source_hashes(sources=None)`**,
  **`h1_union_support_audit(union, details, *, sources=None, target_g=None)`**: `None` means
  union 1's globals, read at call time (so 229/243's harnesses that rebind `H1_UNION_SOURCES`
  still drive them); union 2 passes its own. Union 2's banks are hash-verified before the build
  exactly as union 1's are.
- **`apply_frozen_cd4_residual(..., extra2=None, extra2_audit=None, audits_out=None)`**: the extra
  residuals are accumulated in order — union 1 (`H1_UNION_G`), then union 2 (`UNION2_G`) — each
  sized by `solve_carrying_row_fraction` at its own dose against the same whitened base and added
  to the merged residual with the same operations the union-1 path always used (`merged_effect =
  cd4·frac; merged_effect += u1·f1; merged_effect += u2·f2`; support OR'd; counts `np.maximum`),
  then ONE `add_post_whiten_residual` and ONE global preserve. The return arity is unchanged; union
  2's audit block comes back through `audits_out["union2"]` and `build_production_effects` writes
  it as **`effects.union2`** with the same fields as `effects.h1_union` (sources + four-file shas,
  `carrying_rows`, `g_per_unit_fraction`, `solved_fraction`, `shipped_median_rule_would_solve`, the
  merged preserve) plus `targets_shared_with_h1_union` / `entries_shared_with_h1_union`. CD4's own
  solve (`solve_cd4_dose`, 247's router; `solve_cd4_fraction` underneath) is **not modified**.
- **`run()`** prepares union 2 after union 1 and threads it through both effect builds (the main
  one and `VCC_ALPHA=knee`'s per-context rebuild); the config gains `union2_g`, `union2_sources`,
  `union2_sizing_rule` **only when on** (the `h1_union_g` convention), so a default build's config
  block and audit are byte-identical to every prior build's.

**Unit level** (`unit248.py`, two builders in one process, the working tree vs the byte copy,
synthetic residuals with complete / partial / untouched rows through `apply_frozen_cd4_residual`):
**27/27**. OFF: output array, CD4 audit and union audit bit-identical to pre-248, with and without
union 1. SWAP: union 2 alone == union 1 alone bit for bit (array, CD4 audit, frame key for key,
preserve). ON: union 1 + union 2 == `compose218.merged()` chained twice + one add + one preserve
bit for bit; CD4's fraction and union 1's frame untouched by union 2; union 2's fraction ==
`solve_carrying_row_fraction(whitened, u2, g2)`; rows union 2 does not carry == the union-1-only
output × one scalar (≤ 1e−12); rows it carries move in direction; the one-union-one-dose form
(229's `ARM3` shape) is a different object; `extra2` without `audits_out`, without a dose, or with a
panel coordinate in its support is refused.

**Regression suites** under the changed builder: the builder's own 57-test suite **56/57 under
this builder and 56/57 under the byte copy, 0 regressions, 0 fixed** (`suite248.py`; the one
failure, `test_panel_exclusion_for_operator_and_cd4_residual`, predates 222 — 247 §5.1);
235's hardening tests **9/9** (`work/tests235.log`), including "the second residual is sized by the
carrying-row rule and CD4 is not" and the shipped H1-pair frame to the last bit.

---

## 2. Known-answer row 1 — flags unset is byte-inert

### 2.1 Staging (`stage248.py` = 222's `stage222.py` + union 2 + saved residuals; `stage248.sh`, each row alone and gated; `work/stage/stage248.log`)

| row | builder | env | `effects_array_sha256` | `cd4.fraction` | union fractions | peak RSS |
|---|---|---|---|---|---|---:|
| `off` | working tree `f7a12a70…` | `VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.34 VCC_H1_UNION_G=0.28` | **`abeb0749…` MATCH** | 0.32207201067578145 | h1 0.19207161964444322 | 17.81 GiB |
| `pre_u2` | pre-248 byte copy `2d63b979…` | + `VCC_UNION2_SOURCES=dld1_lowmoi VCC_UNION2_G=0.35` (unknown to it) | **`abeb0749…` MATCH** | same | same | 19.91 GiB |
| `swap` | working tree | `VCC_H1_UNION_G` **unset**, `VCC_UNION2_SOURCES=gse295214_h1,h1_2025 VCC_UNION2_G=0.28` | **`abeb0749…` MATCH** (KA 2) | same | **union2 0.19207161964444322**, 77 carrying | 19.95 GiB |
| `cum1` | working tree | `VCC_H1_UNION_G=0.28 VCC_UNION2_SOURCES=dld1_lowmoi VCC_UNION2_G=0.35` | `6225a58217e72ae04753f9a1d83458cb7c1d0ccd33d17a1a9bc49cec8703d086` (the rider) | same | h1 0.19207161964444322, **union2 0.19847054819075655**, 67 carrying | 19.99 GiB |
| `cum4` | working tree | `VCC_UNION2_SOURCES=dld1_lowmoi,jiang_5line_noK562,igvf_cardiomyocyte,myllia_a549 VCC_UNION2_G=0.35` | `2f9eba51d7c7316dd81e9b541200ae848abaee2e8578aceabdb2f8e36242d7e7` | same | h1 0.19207161964444322, union2 0.20712104380637059, 139 carrying | 19.99 GiB |

The only thing that moves the bytes is the new code responding to the two variables; CD4's and
union 1's solved fractions are byte-identical in all five rows (the single-variable constraint).
`cum4`'s union-2 fraction, 0.20712104380637059, is exactly 229's `ARM4 @0.35` fraction (229 §6.5):
the same four-source object sized by the same rule — what 248 changes is that it now rides on top
of the H1 pair instead of replacing it.

### 2.2 Package level (`run248ka.sh` → `build248.sh` = 222's `build222.sh` + the two passthroughs, both unset, gated, alone)

Both flags unset through the changed builder (commit `2632b92`, builder `f7a12a70…`), 222's build
flags verbatim, strictly serial and alone (the gate opened at free+inactive 30.3 GiB, `memory_pressure`
77% free, zero other builders, zero preps; `work/build/run248ka.log`):

| what | measured | required |
|---|---|---|
| `submissions/248-ka-h1union-g028.h5ad` sha256 / bytes | **`8df4dadd5a63b76c7ca3e3d57a511e7e52000b79d016d140897c55dd24c3ecf6` / 12,648,807,800 — MATCH** | `8df4dadd…` / 12,648,807,800 (the shipped `h1union`; 240 row 1, 247 §3.4, 249 §2.2) |
| `output.nnz` / `effects.final_sha256` | 1,577,411,397 / `abeb0749…` | the shipped meta's |
| `effects.cd4.fraction` / `effects.h1_union.residual_frame.solved_fraction` | 0.32207201067578145 / 0.19207161964444322 | equal |
| `config.union2_*` / `effects.union2` | **absent / absent** | absent (the `h1_union_g` convention) |
| `def_postbuild.py` against the 222 incumbent meta | **32/32 checks pass, `moved=[]`** | — |
| `hashes.implementation` for the builder | `f7a12a70…` (this commit's bytes) | changed, as 235's / 247's / 249's did; the h5ad did not |
| cost | **1477.9 s**, peak RSS **16.72 GiB** (`/usr/bin/time -l`; the meta's own `peak_rss_gib` 16.72) | 249: 1345 s / 20.00 GiB |

The h5ad was re-hashed and deleted as a verified byte duplicate of `submissions/222-h1union-g028.h5ad`;
its meta is kept at `submissions/248-ka-h1union-g028.h5ad.meta.json` and `work/build/`, and is the
same-builder reference a future rule-10 diff of a `248-` arm would be run against. The `.vcc` was not
re-prepped: it is a deterministic transform of the h5ad (240 row 1, 235 row 2, 249 §2.2 all took
`8df4dadd…` to `f13ce40b…`), so the byte identity rests on the h5ad row above.

---

## 3. Known-answer row 3 — the rider through the builder reproduces 225/229's `CUM1` (`repro248.py`, one fold per process; `work/repro248-*.json`, `repro248-summary.log`)

162's PLOLO with the production objects substituted in: both unions from `prepare_h1_union(sources=…)`
and the whole of `apply_frozen_cd4_residual` (CD4 by its rule at 0.34, union 1 = the H1 pair at 0.28,
union 2 at `g`, one merged add, one preserve) for 225's `merged()` chain. Fold caches are byte copies
of 229's (sha-identical). **120 known-answer rows over four folds, all pass:**

- **KA-A** the incumbent (`cd4@0.34 + h1pair@0.28`) through the builder == 218 `compose218`'s
  `cd40.34+h1PAIRCARRY0.28` on `cos_agg`, `a23_cos`, `pds`: **0.00e+00** on all four folds; CD4
  fraction 0.373761 / 0.370782 / 0.376056 / 0.316408 = 218's.
- **KA-B** `CUM1` real, `SHUF` (seed 20726) and row-matched `CD4ROWS` at 225's six doses (0.10 …
  0.60) == 225's `stack225d.json` rows on `d_cos` and `d_pds` to |Δ| < 1e−12, three cancer folds:
  **54/54**; the same for `CUM4`: **54/54**.
- **KA-D** the builder's output at 0.35 == `compose218.merged(merged(cd4, f, h1, f1), 1, u2, f2)` +
  one add + one preserve, **bit for bit**, both unions, all four folds.
- **KA-E** leave-top-3-out on `dld1_lowmoi @0.35` == 225's `conc225.json`: `pds` retained 25.2% /
  cosine 54.5%, same top-3 targets.

The mandate's row, three cancer folds at `g` = 0.35, per-fold mean of ratios:

| fold | union-2 fraction | carrying | rel macro cos | rel `cos_agg` | Δ`pds` | Δ`q_w` |
|---|---:|---:|---:|---:|---:|---:|
| `k562_gwps` | 0.22528 | 60/269 | +1.646% | +2.349% | +0.004582 | +0.00030 |
| `orion_hct116` | 0.23063 | 66/290 | +1.565% | +2.192% | +0.001555 | +0.00009 |
| `orion_hek293t` | 0.23347 | 65/294 | +1.003% | +2.257% | +0.001532 | +0.00008 |
| **mean (3 cancer)** | | | **+1.405%** | **+2.266%** | **+0.002557** (3/3 > 0) | +0.00016 |
| `kolf2` (apart) | | 62/273 | −0.234% | +0.733% | −0.002536 | |

The mandate quotes the row as "rel `cos_agg` +1.405%"; the +1.405% is 225/229's **macro** (per-target
mean) cosine, the functional those tables print. On 243's functional, the energy-weighted aggregate,
the same construction reads **+2.266%**. Both are reported everywhere below; the class rule is priced
on `cos_agg`, as 243 §4.4 states it.

---

## 4. The price — both instruments, all six members, every falsifier (`price248.py` → `work/price248.{json,log}`)

### 4.1 The dose ladder and the falsifiers, through the builder (3 cancer folds; `kolf2` apart)

`CUM1` = `dld1_lowmoi` as union 2 on top of `cd4@0.34 + h1pair@0.28`:

| variant | `g` | rel macro | rel `cos_agg` | Δ`pds` | +`pds` | Δ`q_w` | vs 225 | `kolf2` macro | `kolf2` Δ`pds` |
|---|---:|---:|---:|---:|---:|---:|---|---:|---:|
| real | 0.10 | +0.785% | +0.956% | +0.001266 | 3/3 | +0.00005 | ✓ | +0.226% | −0.000557 |
| real | 0.17 | +1.141% | +1.478% | +0.001834 | 3/3 | +0.00017 | ✓ | +0.229% | −0.001046 |
| real | 0.28 | +1.403% | +2.054% | +0.002495 | 3/3 | +0.00021 | ✓ | +0.015% | −0.001975 |
| **real** | **0.35** | **+1.405%** | **+2.266%** | **+0.002557** | **3/3** | +0.00016 | ✓ | −0.234% | −0.002536 |
| real | 0.45 | +1.238% | +2.376% | +0.002504 | 3/3 | +0.00016 | ✓ | −0.704% | −0.003156 |
| real | 0.60 | +0.715% | +2.133% | +0.002164 | 2/3 | +0.00007 | ✓ | −1.575% | −0.004089 |
| SHUF 20726 | 0.35 | −1.216% | −1.159% | −0.000472 | 2/3 | +0.00018 | ✓ | −0.600% | +0.001230 |
| SHUF 771103 | 0.35 | −1.954% | −1.717% | −0.004065 | 0/3 | +0.00013 | — | −2.558% | −0.004404 |
| row-matched CD4 | 0.35 | +0.401% | **+3.066%** | **−0.005434** | 0/3 | +0.00021 | ✓ | −0.437% | −0.002141 |

Target-shuffle (218's derangement, two seeds) is negative on cosine at **every** dose (seed 20726:
−0.02% … −3.30% macro over 0.10 … 0.60; seed 771103: −0.25% … −4.43%) and on `pds` at every dose
for seed 771103 (seed 20726's `pds` is +0.0003 … −0.0027, inside the instrument's ±0.0015 LOO
sigma at the low doses). **The row-matched CD4 control** (CD4's own information on exactly the
rows DLD-1 carries, same rule, same dose — 225's F3) loses to the real source on macro cosine
(+0.40% vs +1.41%) and on `pds` (−0.0054 vs +0.0026) at every dose, **but wins on `cos_agg`**
(+3.07% vs +2.27% at 0.35, at every dose). On the functional the class rule prices, DLD-1 adds less
aggregate cosine than the free alternative; what it adds that CD4 cannot is direction on the
`pds` channel — and that channel is the concentrated one (§4.2).

Non-carrying rows (218's bar 0.018 pp): carrying rows rel cos **+5.126%**, NOT-carrying **+0.01003%**
(largest excursion 0.01254 pp); Δ`pds` carrying +0.011337, NOT-carrying +0.000028. **Bounded.**

`CUM4` = 225's four as union 2 (`dld1_lowmoi, jiang_5line_noK562, igvf_cardiomyocyte, myllia_a549`):

| variant | `g` | rel macro | rel `cos_agg` | Δ`pds` | +`pds` | Δ`q_w` | vs 225 | `kolf2` macro | `kolf2` Δ`pds` |
|---|---:|---:|---:|---:|---:|---:|---|---:|---:|
| real | 0.28 | +3.897% | +5.172% | +0.003759 | 3/3 | +0.00052 | ✓ | +7.304% | +0.004547 |
| **real** | **0.35** | **+4.037%** | **+5.779%** | **+0.003609** | **3/3** | +0.00055 | ✓ | +8.018% | +0.004906 |
| SHUF 20726 | 0.35 | −3.719% | −2.993% | −0.010058 | 0/3 | −0.00005 | ✓ | −0.701% | +0.001791 |
| SHUF 771103 | 0.35 | −1.608% | −2.423% | −0.004083 | 1/3 | +0.00029 | — | −1.917% | +0.003183 |
| row-matched CD4 | 0.35 | +1.316% | +4.750% | −0.008944 | 0/3 | +0.00024 | ✓ | −0.134% | −0.004139 |

(the full six-dose ladders for both, all variants, are in `work/repro248-summary.log`). `CUM4`'s
non-carrying rows move **+0.04082%** (largest excursion **0.04688 pp — above 218's 0.018 pp bar**):
its merged preserve scalar moves 4.4% (§6) against `CUM1`'s 1.6%, and the renderer's nonlinearity
turns a larger shared rescale into a larger off-row excursion. `CUM4`'s `kolf2` gain (+8.0%) is 2.0×
its cancer-fold gain — 229 §7's inflation warning, reproduced.

### 4.2 Concentration — leave-top-k-out on cosine and `pds` separately (222/225's definition; the gate reads k = 3)

| arm @0.35 | k | Δ`pds` full | Δ`pds` LOO | **`pds` retained** | top-k share | Δcos full | Δcos LOO | **cos retained** |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `CUM1` | 1 | +0.00256 | +0.00170 | 66.6% | 34.7% | +0.000471 | +0.000370 | 78.6% |
| **`CUM1`** | **3** | +0.00256 | +0.00064 | **25.2%** | **81.9%** | +0.000471 | +0.000257 | **54.5%** |
| `CUM1` | 5 | +0.00256 | −0.00019 | −7.3% | 118.2% | +0.000471 | +0.000169 | 35.9% |
| `CUM1` | 10 | +0.00256 | −0.00168 | −65.6% | 185.3% | +0.000471 | +0.000009 | 2.0% |
| `CUM4` | 1 | +0.00361 | +0.00174 | 48.2% | 63.4% | +0.001222 | +0.000987 | 80.8% |
| **`CUM4`** | **3** | +0.00361 | −0.00007 | **−1.8%** | **115.4%** | +0.001222 | +0.000641 | 52.4% |
| `CUM4` | 5 | +0.00361 | −0.00131 | −36.3% | 151.4% | +0.001222 | +0.000438 | 35.9% |

Top-3 Δ`pds` targets, `CUM1`: `ZNF324, LENG1, TRIM35` (k562) · `TARBP2, ZNF324, TTBK2` (hct116) ·
`UROS, PBLD, CYTH2` (hek293t). Three targets of ~290 carry 82% of DLD-1's `pds` gain; the cosine
gain is distributed (54.5% retained). That is 225 §8.2's finding, reproduced through the builder to
the digit, and it is the second of the two conditions the rider fails.

### 4.3 (i) 243's class rule — direction legs with the arm's delivered amplitude, plus the DE-four propagation

Δcos = k_cos × rel `cos_agg` × 0.107332 (h1union's cosine **measured** off its package, 243 — never
the floor); `mse` leg = 197's exact identity at the arm's own terms (h1union's `A_c`/`B_c` scaled by
the arm's delivered amplitude `r`, all three contexts on the `B` leg) minus h1union's raw 0.8172630
(the identity reproduces that board value at h1union's terms to 0.0e+00); `pds` leg = k_pds ×
Δ`pds` / 0.452 / 6; DE-four = 0.239 × Δcos (the h1union step's realised sign-propagation rate;
0.173–0.423 across the three board steps). **No package exists, so `r` is the effect-level delivered
median row norm on frame A (216's instrument): `CUM1` r = 1.00423 (+0.42%), `CUM4` r = 0.99440
(−0.56%)** — the package-level `A_c` stream (`totals243`) is the post-build falsifier and was not
reached. The amplitude leg alone is worth +0.000147 (`CUM1`) / −0.000194 (`CUM4`) overall.

| arm | class (k_cos / k_pds) | Δcos | pred. raw `mse` | `mse` | `pds` | DE-four | **TOTAL (i)** | × gate | `nmae` identity | **TOTAL (ii)** |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **`CUM1`** | **new-source, point (0.94 / 0.70)** | +0.002286 | 0.81183 | +0.00095 | +0.00066 | +0.00055 | **+0.00216** | **1.14×** | −0.00086 | **+0.00107** |
| `CUM1` | new-source, a15-tag k (1.00 / 0.70) | +0.002432 | 0.81153 | +0.00100 | +0.00066 | +0.00058 | +0.00225 | 1.18× | −0.00086 | +0.00114 |
| `CUM1` | new-source, HIGH band (1.00 / 1.15) | +0.002432 | 0.81153 | +0.00100 | +0.00108 | +0.00058 | +0.00267 | 1.41× | −0.00086 | +0.00156 |
| `CUM1` | re-dose class (1.36 / 1.05), sensitivity | +0.003308 | 0.80977 | +0.00131 | +0.00099 | +0.00079 | +0.00309 | 1.63× | −0.00086 | +0.00190 |
| **`CUM4`** | **new-source, point** | +0.005831 | 0.80676 | +0.00184 | +0.00093 | +0.00139 | **+0.00417** | **2.19×** | **−0.00462** | **−0.00106** |
| `CUM4` | new-source, HIGH band | +0.006203 | 0.80602 | +0.00197 | +0.00153 | +0.00148 | +0.00499 | 2.62× | −0.00462 | −0.00028 |
| `CUM4` | re-dose class, sensitivity | +0.008436 | 0.80157 | +0.00275 | +0.00140 | +0.00202 | +0.00617 | 3.25× | −0.00462 | +0.00067 |

Which class DLD-1 belongs to is not clean: of its 67 rows, 48 carry no union residual today (the
"first source into empty rows" case, k_cos 0.94–1.00) and 19 are already H1-pair rows (a re-dose,
k_cos 1.36); all but a handful are CD4-occupied, as h1union's were when 243 measured its 0.94. The
point class is the headline; the re-dose reading is the sensitivity, and even it leaves the second
gate condition (§4.2) failed. DE-four sensitivity at the point class: +0.00040 (low rate) … +0.00097
(high rate) against the +0.00055 carried.

**The controls priced under the same rule** (point class): `CUM1` row-matched CD4 **+0.00058**
(`cos_agg` +3.07% but `pds` −0.0054), SHUF 20726 −0.00067, SHUF 771103 −0.00193; `CUM4` row-matched
CD4 +0.00032, SHUF −0.00457 / −0.00269. The real source beats its row-matched control by +0.00158
(`CUM1`) / +0.00385 (`CUM4`) on the composite — the gain is real; it is the *size* and the *shape*
that fail.

### 4.4 (ii) 216's exact `nmae` identity (`six248.py` → `work/six248.{json,log}`; every arm through the builder's `apply_frozen_cd4_residual`; **34/34 known-answer rows at 0.00e+00** against 216's `six216.json`, 222's `six222.json` and 229's `six229.json`, including 229's `stack6/4/3` single-union arms re-expressed through the new merge path)

| arm | dlv/a15 | `nmae` wald89 | Δ vs incumbent | **overall `nmae`** | Δ top-89 precision | `pds` 3f | ΔU | ΔM | **`e` = ΔM/ΔU** |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `h1union g=0.28` (incumbent) | 1.0033 | 0.98569 | — | — | — | 0.78655 | — | — | — |
| **`CUM1` dld1@0.35** | 1.0076 | 0.98877 | **+0.00308** | **−0.000856** | +0.00186 | 0.78905 | +373.4 | +67.0 | **0.179** — magnitudes |
| `CUM1` dld1@0.28 | 1.0106 | 0.98832 | +0.00262 | −0.000729 | +0.00158 | 0.78895 | +308.4 | +57.4 | 0.186 |
| **`CUM4` four@0.35** | 0.9977 | 1.00233 | **+0.01664** | **−0.004622** | +0.00240 | 0.78942 | +930.1 | −17.8 | **−0.019** — magnitudes |
| `stack6 g=0.35` (229 ARM6) | 1.0006 | 1.00291 | +0.01722 | −0.004783 | +0.00115 | 0.78830 | +959.7 | −31.7 | −0.033 |
| `stack3 g=0.35` (229 ARM3) | 1.0069 | 0.98842 | +0.00273 | −0.000758 | +0.00182 | 0.78779 | +257.8 | +41.2 | 0.160 |

**Two-for-two mis-sign record, as required.** This instrument gave −0.000152 for the CD4 dose step
(board realised +0.001762, R216:555) and −0.000108 for the h1union step (board realised +0.00588
scaled = +0.00098 overall, BR222:49) — wrong in sign both times. So `CUM1`'s −0.000856 is a
**magnitude bound**, not a sign: this arm cannot move `nmae` by more than ~0.0009 overall either
way, and `CUM4`'s −0.0046 is 17× outside the instrument's calibrated range (both calibration points
sit at ~+0.0004–0.0006 raw; `CUM4` is at +0.0166) with 188's mechanism pointing the same way as
229's `ARM6` (`e` < 0: matched mass falls while unmatched mass rises by 930). Composite (ii) replaces
the `nmae` share (43.2%) of the DE-four propagation with the identity's number: **`CUM1` +0.00107,
`CUM4` −0.00106.**

---

## 5. The gate (mandate C), evaluated

| condition | `CUM1` (dld1 @0.35) | `CUM4` (225's four @0.35) |
|---|---|---|
| (i) ≥ 1.5 × 0.0019 = +0.00285 | **FAIL** — +0.00216 at the point class (1.14×); +0.00309 only under the re-dose sensitivity | PASS — +0.00417 (2.19×) |
| `pds` ≥ 0 on 3/3 cancer folds | PASS — 3/3, mean +0.002557 | PASS — 3/3, mean +0.003609 |
| leave-top-3-out retention ≥ 50% on `pds` | **FAIL** — 25.2% | **FAIL** — −1.8% |
| **build?** | **NO** | **NO** |

`CUM1` fails two of three conditions; `CUM4` clears the price and fails the concentration
condition, and the exact `nmae` identity charges it −0.0046 with the same sign-independent
mechanism that killed 229's `ARM6` (`e` = −0.02 vs −0.03), so its composite (ii) is negative.
**No package was built**, per the mandate's "if it does not clear, do not build". The lead's
decision after tonight's `ARM6` reading is unaffected by anything here: there is no third
candidate to weigh.

Two readings that would have changed the verdict, stated so they are not re-derived: (a) pricing
on the macro cosine instead of `cos_agg` lowers (i) further (+1.405% vs +2.266%); (b) the re-dose
class raises `CUM1` to +0.00309, which clears the price condition and still fails retention.
Nothing on this table reaches the gate on both conditions.

---

## 6. The staged rider's geometry at the effect level (`effdiff248.py` → `work/effdiff248.{json,log}`; the split the mandate asks for on the package, stated where the arm exists)

| | `CUM1` | `CUM4` |
|---|---|---|
| staged `effects_array_sha256` | `6225a582…` | `2f9eba51…` |
| union 2 carries / both unions / union 2 only | 67 / 19 / 48 of 300 | 139 / 40 / 99 |
| merged preserve scalar (off 0.9093569351340488 →) | 0.8944861646401419 (s = 0.98365) | 0.8696966167072415 (s = 0.95639) |
| union-2-carried rows: rel displacement min / median / max | 0.139 / 0.304 / 0.544 | 0.045 / 0.313 / 1.259 |
| union-2-carried rows: cosine(on, off) min / max | 0.8824 / 0.9903 | 0.6697 / 0.99995 |
| NOT carried: rel displacement | **exactly 0.016353 on all 233** = |1 − s| | exactly 0.043614 on all 161 |
| NOT carried: cosine(on, off); max ‖on − s·off‖/‖off‖ | 1.0000000000; **9.7e−17** | 1.0000000000; 9.0e−17 |
| scale-corrected separation (carried min / non-carried max) | **1.4e+15×** | 1.0e+14× |
| carried rows' share of ‖ΔE‖²; ‖ΔE‖/‖E‖; rows with zero displacement | 98.975%; 0.1404; **0/300** | 97.998%; 0.2293; 0/300 |
| most displaced targets | LENG1 0.544, SH2D3A 0.506, CHD3 0.490, PRNP 0.457, MAP2K2 0.454 | ZNF507 1.259, PLAGL2 1.087, MYRF 1.027, GLI3 1.016, DNTTIP1 0.768 |

The preserve scalar is one number written three times (`effects.cd4.post_add_preserve`,
`effects.h1_union.merged_preserve`, `effects.union2.merged_preserve`), asserted equal. The rider
applied on every row it carries, nowhere else, and it is not vacuous: the shipped all-rows median
rule would have solved **0.0** for it (`shipped_median_rule_would_solve`, 22.3% coverage), which is
why the carrying-row rule and the second slot exist.

---

## 7. What the rule-10 diff, the profile split, `n_pred` and `vcc prep` would have been

Not run: they are the post-build gate set and no arm package was built. What stands in for them:
§2.2's same-builder rebuild of the incumbent diffs `moved=[]` against the 222 meta through
`def_postbuild.py` (the default path's config block is byte-identical), §6 is the profile split at
the effect level, and §2.1's five staged rows pin every scalar the metadiff would have pinned (CD4's
fraction and frame, union 1's frame, the operator `e2f65ecf…`, α 0.982, `cd4_target_f` 0.34,
`h1_union_g` 0.28) byte-identical across off / cum1 / cum4. The pre-declared diff for a future
`248-` build is: config ADDED exactly `union2_g`, `union2_sources`, `union2_sizing_rule`; audit
CHANGED exactly `effects.final_sha256` and the two merged-preserve scalars under `effects.cd4` and
`effects.h1_union`; audit ADDED only under `effects.union2.`; frozen inputs 0.

---

## 8. Honest limits

- **No board member was measured and no `cell-eval2` was run.** Every price is PLOLO through 243's
  class rule and 216's identity; the rule's factors rest on three board points and are class
  constants, not scalars, and DLD-1's class is mixed (§4.3).
- **The delivered amplitude is effect-level**, not streamed off a package (`totals243`): 249
  measured that instrument reading −0.31% where the package said −0.41%. The amplitude leg here is
  ±0.0002 overall, and the gate verdict does not turn on it.
- **The row-matched CD4 control beats DLD-1 on `cos_agg` at every dose** (§4.1). Under 225's F3
  convention (macro cosine, and `pds`) the source wins; under the functional the class rule prices,
  the free alternative adds more aggregate cosine. This was not a mandate falsifier but it is a
  reason not to read (i)'s `mse` leg as DLD-1-specific.
- **`kolf2` reads negative for `CUM1`** (macro −0.23%, `pds` −0.0025) and inflated for `CUM4`
  (+8.0%, 2.0× the cancer folds); both are reported apart, as every experiment since 218 has done.
- **The instrument's `pds` is concentrated on the same handful of targets whatever the source**
  (225 §8.2 item 3; `ZNF324`, `TARBP2`, `TRIM35` recur here). The retention condition is the
  mandate's; the caveat that it may be an instrument property rather than a source property
  stands, and it cuts the other way for the arm — a gain the instrument cannot resolve is not a
  gain to spend a slot on.
- **`g` = 0.35 was not re-tuned** (225's pre-registered dose); the ladder shows 0.28–0.45 within
  ±0.0001 of `pds` and ±0.3% of cosine of it.

---

## 9. Cost and memory (every ≥ 15 GiB stage alone; the gate logged before each)

| stage | wall | peak RSS |
|---|---:|---:|
| `envcheck248.py` (11 imports) / `unit248.py` | ~20 s / 1.0 s | — / 0.21 GiB |
| staging `off` / `pre_u2` / `swap` / `cum1` / `cum4` (each alone) | 9–10 s after inputs | 17.81 / 19.91 / 19.95 / 19.99 / 19.99 GiB |
| `repro248.py` × 4 folds (one per process) | 42 / 58 / 57 / 52 s | 4.69 / 4.80 / 4.81 / 4.78 GiB |
| `six248.py` (4 folds, 9 arms) | 70.6 s | 7.00 GiB |
| `suite248.py` (57 tests × 2 builders) / 235's tests | 58 s / ~2 min | 0.32 / 1.06 GiB |
| `effdiff248.py`, `price248.py` | seconds | < 0.5 GiB |
| KA build (`run248ka.sh`, both flags unset; the 8df4dadd… rebuild) | **1477.9 s** | **16.72 GiB** |
| arm build / `vcc prep` | **not run** (the gate did not clear) | — |

Nothing from 249 or the armed submissions was touched; `submissions/249-cd4svd5-f034.vcc` (armed for
the second 2026-09-05 slot) and `229-stack6-g035.vcc` (armed for the first) are untouched.

## 10. Files

```
patch248.py        the ONE builder change as 17 exact-anchor replacements; --check reproduces work/make_submission_pm_cd4.248.py
envcheck248.py     the import-time rules, every env combination in subprocesses     -> work/envcheck248.json
unit248.py         27 unit rows, working tree vs the pre-248 byte copy               -> work/unit248.json
stage248.py/.sh    KA rows 1-2 at staging level + the two arms' effect matrices       -> work/stage/
repro248.py        KA row 3 + the PLOLO price, every falsifier, one fold per process  -> work/repro248-<fold>.json, repro248-summary.{json,log}
six248.py          216's exact nmae identity, 9 arms through the builder              -> work/six248.{json,log}
price248.py        both instruments, all six members, the controls, the gate          -> work/price248.{json,log}
effdiff248.py      the staged rider's effect-level split                              -> work/effdiff248.{json,log}
suite248.py        the builder's 57 tests under both builders                         -> work/suite248.{json,log}; 235's tests -> work/tests235.log
build248.sh        222's build flags + the two passthroughs
run248ka.sh        the gated package-level KA rebuild                                 -> work/build/
work/make_submission_pm_cd4.PRE.py   the git control, HEAD 2d63b979..., taken before the first edit
work/make_submission_pm_cd4.248.py   the intended builder (= the working tree = commit 2632b92)
work/builder248.diff                 the change, as git saw it
work/commit-message.txt              the builder commit's message with the known-answer hashes
```

```bash
cd ~/Code/virtual-cell-challenge-2026
git show f59b653b:scripts/make_submission_pm_cd4.py > experiments/248-second-residual/work/make_submission_pm_cd4.PRE.py
uv run python experiments/248-second-residual/patch248.py --check experiments/248-second-residual/work/make_submission_pm_cd4.PRE.py
uv run python experiments/248-second-residual/envcheck248.py
uv run python experiments/248-second-residual/unit248.py
bash experiments/248-second-residual/stage248.sh                       # off / pre_u2 / swap / cum1 / cum4, gated
cp experiments/229-stack4-build/work/plolo/*.np* experiments/248-second-residual/work/plolo/
for f in k562_gwps orion_hct116 orion_hek293t kolf2; do uv run python experiments/248-second-residual/repro248.py --fold $f; done
uv run python experiments/248-second-residual/repro248.py --summary
uv run python experiments/248-second-residual/six248.py
uv run python experiments/248-second-residual/price248.py
uv run python experiments/248-second-residual/effdiff248.py
uv run python experiments/248-second-residual/suite248.py; uv run python experiments/235-october-hardening/tests/test_october_hardening.py
bash experiments/248-second-residual/run248ka.sh                       # the package-level KA row, gated
```

## 11. Verdict

**RIDER NOT BUILDABLE.** The builder can now express 225/229's `CUM` construction (`2632b92`, inert
by default at unit, staging and package level), and expressing it shows the rider is not worth a
slot: `dld1_lowmoi @0.35` on top of the shipped stack prices **+0.00216** under 243's class rule
(1.14× the shape-only gate, below the mandated 1.5×) with **25.2%** of its `pds` gain surviving
leave-top-3-out (bar 50%); `CUM4` prices +0.00417 but keeps −1.8% of its `pds` gain and reads
−0.0046 on the exact `nmae` identity. **Nothing was built and nothing was submitted.**
