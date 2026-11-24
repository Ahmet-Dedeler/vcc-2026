# 319 — the suppression-only variant `normrel15-cap` (`VCC_NORMREL_MODE=cap`): BUILT AND VERIFIED

> **Stage update, September 9:** this is the original build handoff. The
> coordinator subsequently submitted 319 on September 8; it scored
> **0.21865973944346428**. Read `BOARD-RESULT.md` for the published record and
> `COORDINATOR-REVIEW-2026-09-09.md` for interpretation limits. Source `.h5ad`
> files listed here were removed in the documented September 8 storage sweep;
> the `.vcc` and metadata remain.

**319 (`normrel15-cap`): BUILT AND VERIFIED.**

All three known-answer rows pass — the first two **bitwise, 0.000000e+00** — the cap is a **real
arm and not a near-incumbent** (it touches **49.94 %** of the real panel's (target, donor) cells and
keeps **52.5 %** of the shipped arm's median weight move away from uniform), the single full build
landed with `config.normrel_exp = -1.5` **and** `config.normrel_mode = "cap"`, postbuild is **33/34
with exactly one FAIL — the expected config-drift line, now naming the new key**, the config differs
from `submissions/312-normrel15.h5ad.meta.json` in **exactly one key**, prep succeeded, and the
novelty check returns `OK: novel package`.

**Nothing was submitted. `vcc submit` was never run.** The coordinator preregisters the reading.

| artifact | value |
|---|---|
| h5ad | `experiments/319-cap-variant/work/build/the-parrot-approved-my-expenses.h5ad`, 12,648,807,800 B, sha256 **`f216eea700614b396603496332410543851359d90bc36534941a53929a6515c5`** |
| package | `submissions/319-normrel15-cap.vcc`, 2,824,294,400 B, sha256 **`57585d558f9777f066c02f5b299cacccca0dd4954b4201b8fa653bd0f6816dd2`** |
| meta beside it | `submissions/319-normrel15-cap.h5ad.meta.json`, sha256 `3389633241a0ac63…` (byte copy of the build's meta) |
| builder, patched | `scripts/make_submission_pm_cd4.py` **`1318cd3c7dbe8ed1…`** (was `d87396e5c9467e38…`, HEAD `6095137`) |
| driver, patched | `scripts/def_pipeline.py` **`3fb0c4bcc1e29c2c…`** (was `b356cae45fff0eef…`) |
| revert copies | `work/builder.PRE319.py`, `work/def_pipeline.PRE319.py` (byte copies taken before the patch) |
| the patch | `apply_cap_patch319.py` (idempotent), diff at `work/patch319.diff` (52 lines) |
| `nnz` | 1,577,383,789 = 33.2 % of the density cap (312: 1,577,401,690 · 316: 1,577,896,070 · `cd4svd5`: 1,577,615,064) |

---

## 1. The diff

Two files, both authorised; `git diff scripts/` is saved verbatim at `work/patch319.diff`. The
builder change is **312's staged `apply_cap_patch.py` applied unmodified**; `apply_cap_patch319.py`
is that script with two changes, neither touching the arm: the pre-patch byte copy goes into
**this** experiment's `work/` instead of 312's, and step 4 registers the knob in the driver (the
mandate asks for it; 312's staged script did not do it).

### (a) `scripts/make_submission_pm_cd4.py:574-576` — the knob

```python
NORMREL_MODE = (os.environ.get("VCC_NORMREL_MODE", "plain").strip() or "plain").lower()
if NORMREL_MODE not in ("plain", "cap"):
    raise SystemExit(f"VCC_NORMREL_MODE must be 'plain' or 'cap'; got {NORMREL_MODE!r}")
```

Read once at import; default **`plain`** = today's shipped behaviour. Measured at import, all four:

| env | measured |
|---|---|
| unset | `NORMREL_MODE 'plain'` (and `NORMREL_EXP None`) |
| `VCC_NORMREL_MODE=` (empty) | `'plain'` — the driver's "this arm is off" convention |
| `VCC_NORMREL_MODE=CAP` | `'cap'` — the value is lower-cased |
| `VCC_NORMREL_MODE=squash` | `SystemExit: VCC_NORMREL_MODE must be 'plain' or 'cap'; got 'squash'` |

*Cosmetic note for whoever reads this next:* the block lands **between** the `VCC_NORMREL_EXP`
comment and the `NORMREL_EXP_ENV` line it documents. That is where 312's staged patch put it and I
left it there rather than deviate; it is a comment-adjacency wart, not a behaviour.

### (b) `scripts/make_submission_pm_cd4.py:1531-1541` — the arm, four lines inside `_normrel_pre_whiten_effects`

```python
             if med > 0.0 and nv > 0.0:
-                w[name] = float((nv / med) ** exponent) * float(
+                ratio = nv / med
+                if NORMREL_MODE == "cap":
+                    # suppress only the OVERSIZED rows; never up-weight a small one
+                    ratio = max(1.0, ratio)
+                w[name] = float(ratio ** exponent) * float(
                     active_weights.get(name, 1.0)
                 )
```

`nv` and `med` are Python floats, so `max(1.0, ratio)` is exact. Nothing else in the function moves:
the fallback for `med <= 0` or `nv <= 0` (uniform weight), the renormalisation inside
`model.effect`, and the norm match are untouched, and the clamp sits in the one expression both
renderings (`VCC_NORMREL_PM` unset and `=blend`) share.

### (c) `scripts/make_submission_pm_cd4.py:5521` — the meta config

```python
            **({} if NORMREL_EXP is None else {"normrel_mode": NORMREL_MODE}),
```

The `cd4_sizing` / `normrel_exp` convention: written only when the per-target weighting is ON. A
build with `VCC_NORMREL_EXP` unset — `cd4svd5` and every arm before 312 — writes no new key at all.
Note the consequence, which is intended: a **future `plain` build at `exp = -1.5`** now writes
`normrel_mode: "plain"` and will therefore show one more key than `312-normrel15` did.

### (d) `scripts/def_pipeline.py:126-131` — `ARM_BOARD_DELTA`

`VCC_NORMREL_MODE` registered in exactly the shape 312 used for `VCC_NORMREL_EXP` and 316 for
`VCC_NORMREL_PM`, so the driver announces it as a declared arm rather than an unexplained key.

---

## 2. The three known-answer rows

`test319.py`, on **313's saved fixture** (`work/fixture313.npz`, `723ce13c0c3b2e63…`, seed 3130907,
20 targets × 60 genes × 3 donors), through the **real production functions** on a real
`MultiDonorModel`, against a **byte copy of the builder as it stood before the patch**
(`work/builder.PRE319.py`, `d87396e5c9467e38…` — the 316-patched builder that produced
`316-normrel15-pm`) loaded as an independent module, 316's pattern. Output `work/test319.json`,
`work/test319.log`. Peak RSS **0.22 GiB**.

| row | claim | measured | tol | verdict |
|---|---|---|---|---|
| **KA-1** | `VCC_NORMREL_MODE` **plain**, exp −1.5 == the PRE-319 builder's arm, **bitwise** | **0.000000e+00** | 0 | **PASS** |
| KA-1′ | the `pm_audit` dict is unchanged by the patch | identical | — | PASS |
| KA-1″ | the module's default `VCC_NORMREL_MODE` is `plain` | `'plain'` | — | PASS |
| **KA-2** | `VCC_NORMREL_MODE=cap` with **exp 0** == the PRE-319 builder's **incumbent**, **bitwise** | **0.000000e+00** | 0 | **PASS** |
| KA-2′ | ... and its `pm_audit` is the incumbent's | identical | — | PASS |
| **KA-3a1** | every donor with **r ≤ 1** keeps **exactly** its uniform raw weight under `cap` | **0.000000e+00** (27 cells) | 0 | **PASS** |
| **KA-3a2** | every donor with **r > 1** is **strictly down-weighted** (exp −1.5) | 0 exceptions (26 cells) | 0 | **PASS** |
| **KA-3b** | the shares the **patched builder actually used**, recovered from its own output by least squares (313's KA-2b recipe), == the capped rule | **7.71e-09** | 1e-5 | **PASS** |
| KA-3b′ | ... and, on the identifiable subspace, for the fixture's case-(d) row `T005` (one donor row is numerically zero, so least squares cannot see its weight) | 5.38e-09 | 1e-5 | PASS |
| KA-3c | the capped shares are **not** the plain rule's — the arm is not a relabelling | max \|dW\| **0.4232** | >1e-4 | PASS |
| KA-4 | 316's recorded shipped-arm delivered-norm medians on this fixture reproduce | 0.000000e+00 (both) | 1e-6 | PASS |

**12/12.** KA-1 and KA-2 are the mandate's rows 1 and 2 discharged at **fixture scale, bitwise**:
each production-scale byte-for-byte row would need its own ~27-minute build and the mandate budgets
one. What replaced them at production scale is §4 — one-key config diff, an identical PerturbMatch
audit block, and the stage digests.

**KA-3a is the mandate's paragraph, checked**: "a donor at or below its own panel median keeps
weight 1 and only oversized rows are down-weighted" is true of the **raw** `w` dict, exactly (0.0e+00
on all 27 at-or-below cells). It is *not* true of the **share** each donor ends up with: `model.effect`
renormalises `w` over the carrying donors, so suppressing one donor raises the others' shares. Both
statements matter and they are different; the arm acts through the shares.

---

## 3. How much the cap touches — measured on the REAL 300-target panel, not on the fixture

`capshare319.py`. Production builds the three-donor model with `gene_norm=False, eb=False`
(`MB:4972-4977`), so `model.s[name]` **is** `Donor.load(path).decompose()[1]` and the ratio
`r(t, j) = ‖s_j[row(t)][dcols_j]‖ / median over the panel rows j carries` can be read off the donor
banks without the 20 GiB build. Known-answer rows: the same function reproduces `test319.json`'s
fixture block from the real `MultiDonorModel` (**0.0e+00**) and the memory-light row-wise path
reproduces the full `decompose()` path on a real bank (k562_gwps, 269 panel rows, **0.0e+00**).
Peak RSS **5.85 GiB**, 2/2 KA. Output `work/capshare319.json`, `work/capshare319.log`.

| donor | panel targets carried | panel median row norm | r median | r range | **r > 1** | q05 / q25 / q75 / q95 |
|---|---:|---:|---:|---|---:|---|
| `k562_gwps` | 269/300 | 6.938586 | 1.000000 | 0.3463 – 2.7970 | **134/269 = 49.8 %** | 0.7204 / 0.8769 / 1.1920 / 1.7218 |
| `orion_hct116` | 290/300 | 7.194266 | 1.000000 | 0.4766 – 2.6649 | **145/290 = 50.0 %** | 0.6612 / 0.8269 / 1.3302 / 1.9201 |
| `orion_hek293t` | 294/300 | 6.251894 | 1.000000 | 0.4935 – 2.9446 | **147/294 = 50.0 %** | 0.6524 / 0.8351 / 1.3333 / 1.9996 |

> **The cap touches 426 of 853 (target, donor) cells on the real panel = 49.94 %.** It could hardly
> be otherwise: `r` is normalised by *that donor's own panel median*, so **half of every donor's
> cells sit above 1 by construction**. The mandate's escape hatch — "if the cap touches almost
> nothing, say so" — **does not fire**. This is a half-panel intervention.

What that means per target, and this is the part that is specific to `cap`:

| | count |
|---|---:|
| targets where **no** carried donor is capped → the weights collapse to uniform and the row is the **incumbent's**, bit for bit | **77/300** |
| targets where **some but not all** are capped | **149/300** |
| targets where **every** carried donor is capped | **74/300** |
| capped/carried histogram | `0/1:4, 0/2:13, 0/3:60, 1/1:2, 1/2:8, 1/3:68, 2/2:14, 2/3:73, 3/3:58` |

**Weight move away from the incumbent's uniform shares**, per target, max over donors:

| rule | median | mean | max |
|---|---:|---:|---:|
| `cap` (this arm) | **0.04336** | 0.06317 | 0.30902 |
| `plain` (the shipped arm, `Q519I2pLWRCW7PAEx5ul`) | **0.08261** | 0.10024 | 0.38560 |
| `cap` vs `plain` | 0.03171 | 0.03980 | 0.29908 |

**The cap keeps 52.5 % of the shipped arm's median weight move.** Read it as: the same exponent, but
only the suppression half of the rule fires, so the arm is roughly *half the shipped dose,
one-sided* — plus a quarter of the panel handed back to the incumbent untouched. That is the object
the board will price; it is **not** a smaller version of `exp = -0.75`, because the targets it
reverts and the targets it moves are chosen by the sign of `log r`, not by dose.

---

## 4. The build

One invocation of the gated two-invocation pattern (`run319.sh`, modelled on 316's, which is
312's), on the real A/B/C bundle. `work/build/run319.log`, `work/build/driver.log`,
`work/build/logs/`, `work/build/timings.json`.

```
2026-09-08T04:30:36Z ARMED  VCC_NORMREL_EXP=-1.5  VCC_NORMREL_MODE=cap  VCC_NORMREL_PM=<unset>
  builder 1318cd3c7dbe8ed1…  driver 3fb0c4bcc1e29c2c…  git 4b8bb13
  gate: free 81%  free+inactive 31.1 GiB  other builders 0  vcc preps 0  disk 23.8 GiB
2026-09-08T04:57:53Z invocation 1 rc=1   (the expected postbuild FAIL; artifact written)
  restored scripts/make_submission_pm_cd4.py from pin-artifacts.builder-before.py
    (2b9f5065019449d6… -> 1318cd3c7dbe8ed1…, identical)
```

**The postbuild line, verbatim:**

```
33/34 checks pass       FAIL lines: 1
  FAIL  config differences vs 249-cd4svd5-f034.h5ad.meta.json are only
        ['contexts', 'expect_contexts', 'targets', 'whitener_id']:
        moved=['expect_contexts', 'normrel_exp', 'normrel_mode', 'whitener_id']
```

The same single FAIL 312 got (`33/34`, `moved=['expect_contexts', 'normrel_exp', 'whitener_id']`)
with `normrel_mode` added: the guard noticing a deliberate arm. `run319.sh` refuses to prep unless
the FAIL count is exactly 1 **and** that FAIL is the config-drift line. `[arms]` passed all seven
value assertions (alpha 0.982, `cd4_target_f` 0.34, `h1_union_g` 0.28, the two union sources,
`cd4_sizing` absent, `cd4_deflate_k` 5, `cd4_deflate_basis` absent).

**The arm landed:** `config.normrel_exp = -1.5`, `config.normrel_mode = "cap"`, `normrel_pm`
**absent** (the shipped weight-first rendering, i.e. 312's, not 316's).

**The one-key config diff** (`configdiff319.py`, run as a hard gate before the prep; it was
smoke-tested first on 316's meta, where it correctly refuses):

| vs | keys differing | |
|---|---:|---|
| **`312-normrel15`** (entry `Q519I2pLWRCW7PAEx5ul`) | **1** | `normrel_mode: <absent> → 'cap'` |
| `249-cd4svd5-f034` | 4 | `normrel_exp`, `normrel_mode`, plus `expect_contexts` and `whitener_id union_pruned → union_pruned_def` — the same two driver keys 312 and 316 moved |

**Where the change entered, from the meta's own stage digests:**

| digest | **319-cap** | 312-normrel15 | 316-parity | 249-cd4svd5 |
|---|---|---|---|---|
| `baseline_pre_sha256` | `b2522cc4337389be` | `9e63576df8436257` | `92555af13c543f2d` | `7678bab9ea9d1b8d` |
| `candidate_pre_sha256` | `571be04401b79f87` | `9b19eff922fc2278` | `92555af13c543f2d` | `0601336fdf1ec5c5` |
| `post_whiten_sha256` | `09c0de412da6bb2c` | `f3df0f6426ff1bdc` | `218e399b4d83df6e` | `b2426e15d76d9fe5` |
| `final_sha256` | `db1da6e5362deec7` | `ab8997931a48a191` | `6c25cdc875b92878` | `a28d5be2e9cb3ed1` |

For 319 `baseline_pre ≠ candidate_pre`, as for 312 and `cd4svd5` and unlike 316 — the machine-checkable
signature that PerturbMatch was applied **outside** the norm match, i.e. that this is one variable
away from the shipped arm and not two. The PerturbMatch audit block is **identical across all four
builds** (9,373 observed genes, scale 4.562123, 21 / 3 preserved), so the cap did not touch that
channel. `effects.whitening.after_all_median_norm` 4.292440 (312: 4.272970, `cd4svd5`: 4.229296,
316: 4.221715) — rule 10's cheap alpha falsifier, no drift of a size that reads.

**Resources**, every heavy step serial behind the ≥ 40 % free, ≥ 23 GiB free+inactive,
no-other-builder, no-other-prep gate:

| stage | s | peak GiB |
|---|---:|---:|
| whitener | 328.5 | 12.32 |
| repin-panel-support | 10.7 | 19.44 |
| **build** | **1274.2** | **20.14** |
| postbuild | 0.8 | 0.12 |
| (total, 10 stages) | 1630.4 | — |

My own analysis jobs peaked at **0.22 / 5.85 / 1.75 GiB**; the build and the prep are the mandate's
explicitly authorised heavy steps.

---

## 5. Prep and novelty

```
✓ prep complete
  input:         experiments/319-cap-variant/work/build/the-parrot-approved-my-expenses.h5ad
  output:        submissions/319-normrel15-cap.vcc
  cells:         360000
  genes:         18533  (encoding: float32)
  perturbations: column 'target_gene' → 'target_gene'
  contexts:      column 'context' → A=120000, B=120000, C=120000
  targets:       verified against the official list
  normalization: counts-preserved
       18.44 real   maximum resident set size 11,968,364,544 B (11.15 GiB)
```

The driver's `prep` stage is `vcc prep --dry-run` (validation only, writes nothing), so the package
was written by a real `vcc prep`, gated, in its own step.

```
$ uv run python scripts/check_submission_novelty.py check submissions/319-normrel15-cap.vcc
hashing submissions/319-normrel15-cap.vcc (2.82 GB)...
prediction identity: f216eea700614b39... (output.sha256 of 319-normrel15-cap.h5ad.meta.json (container shape agrees))
content identity: 359b446f3165b3e1... (content of pred.h5ad.zst (X, obs index, target_gene, context, var index))
OK: novel package, sha256 57585d558f9777f0...
                                                                              rc=0
```

**STOP. `vcc submit` was never run and no slot was spent.**

---

## 6. Addendum (not asked for, and it is the number the mandate's Correction-5 rule exists to protect)

Every per-target norm this family has quoted as "delivered" came from 313's synthetic fixture until
318 measured it on the real thing. So the cap arm's delivered profile is measured here **where it
acts** — on the emitted counts of the 12.6 GB predictions — with 318's own recipe
(`deliver319.py`, `p(t,c)` summed over the 400 emitted cells, `effect = log1p(50000·p/Σp) − b_ctrl(c)`,
900 blocks). **Known-answer row: 318's published shipped/parity summary reproduces at
0.000e+00 on all three of median, sd(log) and whole-matrix.** Peak RSS 1.75 GiB, 21.6 s.

| pair (real panel, 900 blocks) | m median | m range | sd(log m) | whole matrix | row cosine median / min |
|---|---:|---|---:|---:|---|
| **`cap` / incumbent** | **0.994774** | 0.9368 – 1.3368 | **0.038050** | **1.0103362** | **0.996614** / 0.829256 |
| shipped / incumbent | 0.990479 | 0.9132 – 1.3328 | 0.044363 | 1.0050947 | 0.987622 / 0.765240 |
| `cap` / shipped | 1.009427 | 0.8989 – 1.0353 | 0.015394 | 1.0052149 | 0.998002 / 0.881484 |
| shipped / parity (318's pair, the KA) | 1.011353 | 0.9321 – 1.2349 | 0.032895 | 1.0085497 | 0.998378 / 0.948492 |

Two readings, both on the real panel:

1. **The cap rotates the delivered rows about half as far as the shipped arm does.** Median
   `1 − cos` against the incumbent is 0.003386 for `cap` against 0.012378 for the shipped arm; as an
   angle that is **0.52×**, which agrees to the second digit with the **52.5 %** median weight move
   measured independently in §3 from the donor banks alone. Two different instruments, one number.
2. **The cap is not the parity fix.** Its delivered amplitude spread against the incumbent
   (sd(log m) 0.0381) is of the same order as the shipped arm's (0.0444), and its whole-matrix ratio
   is *further* from 1 (1.01034 vs 1.00509). `VCC_NORMREL_MODE` is a **weighting** hypothesis; the
   amplitude question is 316's and it is separately answered (and lost, −0.00075 on the board).

Nothing here is a price. The arm is unpriced by design — the mandate reserves the preregistration
and the reading for the coordinator.

---

## 7. State left behind

* **`scripts/make_submission_pm_cd4.py` and `scripts/def_pipeline.py` are PATCHED**
  (`1318cd3c7dbe8ed1…`, `3fb0c4bcc1e29c2c…`), and **committed** with the known-answer hashes in the
  message, staging those two paths and this experiment directory explicitly — AGENTS §11.1b ("a
  build agent that modifies production code MUST commit it as part of *done*", explicit paths, never
  `git add -A` with agents running). Byte copies of the pre-patch files are at
  `work/builder.PRE319.py` / `work/def_pipeline.PRE319.py`. `apply_cap_patch319.py` is idempotent
  and re-derives the patch from the pre-patch text.
* **Under `submissions/`: only** `319-normrel15-cap.vcc` and its sibling meta. The system preamble's
  blanket "never write into `submissions/`" is overridden for those two files by this mandate's
  explicit instruction ("`vcc prep` to `submissions/319-normrel15-cap.vcc`, copy the meta beside
  it"), the same way 312 and 316 wrote theirs. **`vcc submit` was never run.**
* **Disk:** the build directory holds 13.2 GiB (h5ad 11.8 + whitener 1.2 + `cd4_compact` 0.2). Free
  space on the volume went 23.8 → 16.8 GiB during the build and reads **24.4 GiB** now; nothing was
  deleted. One operational note: `run319.sh`'s gate carries a **20 GiB free-disk** condition sized
  for the whole build, and it therefore blocked the *prep*, which needs 2.6 GiB, once the build had
  spent the space. I killed the waiter rather than editing a running script (the standing rule) and
  ran the prep as its own gated step; the log records that.

## 8. Limitations

* Mandate rows 1 and 2 are **fixture-scale bitwise**, not production-scale byte-for-byte: two more
  full builds would be needed and the mandate budgets one. §4 lists what replaced them at
  production scale.
* 313's fixture is synthetic (313 §9). Structural statements from it — that `plain` reproduces the
  shipped arm, that exponent 0 is the identity, which donors are capped, that r ≤ 1 keeps weight 1 —
  are exact and transfer. **Its magnitudes do not**, and none is quoted here as delivered: every
  delivered number in §3 and §6 is measured on the real 300-target panel and labelled with its
  stage (donor-bank row norms in §3, emitted-count log-effects in §6).
* The arm is **unpriced**. No PLOLO retrodiction, no `nmae` identity leg, no six-member forecast,
  no board reading — the mandate stops at the package.
* `VCC_NORMREL_MODE=cap` was built at **one** exponent, −1.5, deliberately, so the cap is the only
  variable against `Q519I2pLWRCW7PAEx5ul`. The hypothesis that motivated it — that suppression alone
  lets the *dose go further* — is **not** tested by this build; testing it needs a second point at a
  larger |exp| under `cap`, and §3's numbers say what to expect there (at −1.5 the cap already
  delivers half the shipped weight move, so a `cap` dose curve is shifted, not truncated).
* `BETA = 0.0` throughout, so the `beta * m_full` offset inside the norm match never fires.

## 9. Files

| file | note |
|---|---|
| `apply_cap_patch319.py` | the patch, idempotent; saves the pre-patch byte copies; registers the arm in the driver |
| `lib319.py` | the independent statement of the weight rule (both ratio paths, raw and renormalised weights) |
| `test319.py` | the three known-answer rows on 313's fixture, 12/12 |
| `capshare319.py` | the cap's reach on the real 300-target panel, from the donor banks |
| `configdiff319.py` | the one-key config gate, run before the prep |
| `deliver319.py` | §6's delivered-magnitude profile, with 318's summary as its known-answer row |
| `run319.sh` | the gated build (+ the prep step, which was re-armed by hand — §7) |
| `work/patch319.diff` | `git diff scripts/`, 52 lines |
| `work/test319.json`, `work/test319.log` | KA rows, cap reach on the fixture, fixture delivered norms |
| `work/capshare319.json`, `work/capshare319.log` | per-donor `r`, per-target capped counts, weight moves |
| `work/deliver319.json` | the four delivered-magnitude pairs, per-block rows for `cap`/incumbent |
| `work/build/run319.log`, `driver.log`, `logs/postbuild.log`, `timings.json`, `prep.log` | the build and the prep |
| `work/build/the-parrot-approved-my-expenses.h5ad{,.meta.json}` | the artifact, `f216eea700614b39…` |
| `work/builder.PRE319.py`, `work/def_pipeline.PRE319.py` | revert copies |

Sources read: `docs/READ-FIRST.md`, `AGENTS.md` (§6, §7, §11, §14, §6.2b),
`experiments/312-normrel-arm/{apply_cap_patch.py,apply_patch.py,test312.py,run312.sh}`,
`experiments/313-normrel-parity/{parity313.py,RESULTS.md}`,
`experiments/316-parity-arm/{apply_patch316.py,test316.py,run316.sh,RESULTS.md,work/test316.json}`,
`experiments/318-per-target-magnitude/{RESULTS.md,profile318.py,work/profile318.json}`,
`scripts/{make_submission_pm_cd4.py,def_pipeline.py,check_submission_novelty.py}`,
`src/vcc2026/{multi_donor.py,donor_transfer.py}`,
`submissions/{312-normrel15,316-normrel15-pm,249-cd4svd5-f034}.h5ad.meta.json`.
