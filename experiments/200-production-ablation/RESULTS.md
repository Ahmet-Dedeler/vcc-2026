# 200 — Leave-one-component-out on the shipped stack

**Every component of a15 has now been removed and re-measured on all six members. Six are
load-bearing, four are inside their gate, and none is net negative outside it. The stack is
coherent and every component earns its place or costs nothing.** That is the headline, and it is
the first time anyone here has been able to say it.

Two things were found on the way, and both are worth more than the headline:

1. **`--overdispersion 4.0`'s `reach` value is not a property of the knob.** In a8's flat-0.25
   context it is **+0.00393 raw** on h1-test and retrodicts the board's +0.021398 scaled in
   sign. In a15's tiered context the *same knob* is **−0.02610 raw** (t = −7.49). The board
   point that justified `--overdispersion 4.0` was measured in a context we no longer ship, and
   in the context we do ship it points the other way.
2. **One shipped *setting* is a live candidate worth up to +0.0056 overall, and it is exactly
   the one the mandate flagged.** Replacing the tiered spread policy with a **flat 1.00** spread
   holds `mse` at the identical cap (claim 0.24653, branch `s²B`, raw `mse` 0.88228 — bit-equal
   to a15's) and holds `pds` by the per-gene-total identity, so it is a free four-member trade.
   It scores **−0.005580** for the tier under the incumbent `jac` convention and **+0.002585**
   for the tier under the call-set-scaled one. **The verdict flips on a `jac` convention this
   experiment cannot adjudicate**, so it is handed to 192 rather than recommended.

`a23`'s mass fix — the other component the mandate singled out — **has never cleared, on any
instrument**: the board measured it at **−0.00073 overall**, this experiment at **+0.00063**,
and the shape gate is ±0.0019. It is a wash, not a `reach` problem.

---

## Phase 0 — known-answer rows

All gates PASS. Every one is a *reproduction*, not a re-derivation: where a prior experiment's
artifact exists on disk it is loaded and scored, not rebuilt.

### `pds` — PLOLO, 162's pre-committed 4-fold unweighted mean, challenge frame

| row | value | expected | |
|---|---:|---:|---|
| SHIP absolute, `orion_hct116` / `orion_hek293t` / `k562_gwps` / `kolf2` | 0.791823 / 0.777375 / 0.789121 / 0.610128 | — | 4-fold mean **0.742112** |
| SHIP vs SHIP-dup (field-identical configs, state-carry control) | worst \|Δ\| **0.00e+00** | 0 | **PASS** |
| **a15 → a23 mass fix** (`noMassFix` LOO) | **−0.004014 scaled** | −0.004014 | **PASS, exact** |
| 169's `q = 1.00` identity vs 162's own `work/arms.json` | **24/24 cells, worst \|Δpds\| 0.000e+00** | 0 | **PASS** |

The mandate specified 16 cells for the `q = 1.00` identity; the four folds × three challenge
frames × two arms available here give 24, and all 24 are bit-identical.

### `mse` — 096's `claim` route on the real challenge panel

| gate | got | expected | |
|---|---:|---:|---|
| **`claim(a15 as shipped)`** | **0.24581** | 0.24581 | **PASS (+0.00%)** |
| **a15 raw `mse` / `score_mse`** | **0.88228 / 0.11112** | 0.88228 / 0.11112 | **PASS** |
| flat 0.25 / φ0 | 0.11860 | 0.11860 | PASS (−0.00%) |
| flat 0.25 / φ4 (a8's law) | 0.16652 | 0.16478 | PASS (+1.06%) |
| flat 1.00 / φ0 | 0.19438 | 0.19216 | PASS (+1.16%) |
| A/B round trip at s = 1 | A 0.283614, B 0.281236 | board_solve | PASS (+0.000%) |

### The four DE members — h1-test, emitted

The mandate's expected row was `n_pred` **6035** / `reach` **0.18283**. Only the second is
right, and the first is right in a different counter:

| gate | got | expected | source |
|---|---:|---:|---|
| a15 `reach` f = 1 | **0.18283** | 0.18283 | 157 `power_h1test.log`, 164 `decomp.log`:6 |
| a15 median `n_pred` | **6034** | **6034** | 164 `decomp.log`:29 |
| 169's five-seed ladder, d`reach` | **+0.00818** | +0.00818 | 169, on 169's own convention |
| 169's five-seed ladder, d`nmae` | **+0.04640** | +0.04640 | 169 |

Two conventions had to be separated before these could pass, and both are now recorded in
`six200.py`:

* **`n_pred` 6035 vs 6034 is the self-gene.** `realize164.log` and `realize200.log` print the
  median call count *before* the BRIEF 4-7 self-exclusion; `debench.Surface.counts()` prints it
  after. 164's own `decomp.log` line 29 reads **6034** on the same file. Same table, two
  counters, one gene apart.
* **`reach` 0.18283 vs 0.17638 is the empty-pool rule.** At f = 1, three of h1-test's 85
  perturbations have a native BH set consisting of their own target gene alone, which the
  self-exclusion then empties. `debench.Surface.reach()` **drops** those from the mean; 169's
  `per_pert_reach` **counts them as zero**. The two therefore differ by the exact constant
  **82/85 = 0.964706** — 0.18283 × 0.964706 = 0.176378 — and at f ≥ 2 the pools are no longer
  empty and the two agree to every printed digit. It is a property of the *reference*, identical
  for every arm, so it rescales every delta and reverses no sign. This table uses
  `debench.Surface.reach()`, the convention every board-linked `reach` number is on, and carries
  169's zero-filled variant alongside so 169's published number stays checkable.

**This was the gate that failed in this experiment's first pass.** It was not an emitter fault
and not a seed mismatch: `de200.KA_a15_raw.npz` is **byte-identical** to `de164.a15.npz`
(sha256 `79c81c2338…`), and scoring 164's own artifact reproduced 0.17638 as well. The
disagreement was between two scorers, and it is now resolved in favour of the production one.

### Rule 11 — path and sha256 for every operator and frozen input loaded

`work/inputs200.json` hashes **13 operator files** under near-identical names (12 distinct
sha256; `095/out/whitener_union.npz` and `096/out-alpha/whitener_union.npz` are the same file).
These are the five this experiment actually loaded, and which fold loaded which:

| path | sha256 (16) | max eig | n | used for |
|---|---|---:|---:|---|
| `data/donors/whitener_union_pruned.npz` | `e2f65ecf69b7e0b0` | 5.894570 | 17942 | `kolf2` fold; **matches a15's and a23's `operator_sha256`** |
| `experiments/062-…/whiteners/held_k562_gwps.npz` | `13d988fcc3b129d2` | 8.425882 | 17806 | `k562_gwps` fold |
| `experiments/062-…/whiteners/held_orion_hct116.npz` | `7cecbf7f36e40d22` | 8.036340 | 17942 | `orion_hct116` fold |
| `experiments/062-…/whiteners/held_orion_hek293t.npz` | `97e486e0d1a805e1` | 3.695842 | 17942 | `orion_hek293t` fold |
| `experiments/095-clean-harness/out-h1test/whitener_union.npz` | `f80fe8d6ba59138b` | 8.080842 | 17745 | every h1-test rebuild |

The trap the mandate warned about is real and is visible in this table:
`experiments/095-clean-harness/out/whitener_union.npz` has max eigenvalue **8.036353** against
`held_orion_hct116.npz`'s **8.036340** — the same number to five significant figures, and a
completely different sha (`5cf3899db9272ff2` vs `7cecbf7f36e40d22`). Only the sha separates them.

Frozen inputs (full record in `work/inputs200.json`, 22 entries):
`scripts/make_submission_pm_cd4.py` `90c4f49eebc1c706`, `submissions/a15-cd4f010.vcc.meta.json`
`34223442e4e8816b`, `submissions/a23-massfix.h5ad.meta.json` `1c721002e7071598`,
`data/harness/h1-test.h5ad` `06c6baea8c62ef3d`, 095's h1-test `pred.h5ad` `26df35c8a694f6e1`,
169's `h1test-effect.npz` `e4d8af8b892cb22d`, 162's `arms.json` `8401e4dc82a8a3a0`,
`de164.a15.npz` ≡ `de200.KA_a15_raw.npz`
`79c81c2338d4d489c6503caecb0a49f8f96c229f7c5411241afb4e0fdf1581e9`.

---

## Phase 3 first — the four board retrodictions

Audit 181 rule 1. Full output in `work/retro200.json`.

### R1 · a15 → a23, the mass fix — **passes, and settles the question**

| | board | here |
|---|---:|---:|
| overall | **−0.000730** | **+0.000628** |
| `reach`, scaled | −0.007933 | −0.000295 (overall units) |
| `reach`, raw | −0.007025 | −0.003989 (h1-test, sem 0.00266, t −1.50) |

Both readings are **inside the ±0.0019 shape gate**, from opposite sides. The implied h1-test →
board `reach` transfer for this profile-changing arm is **1.76** — the same order as the
`rc_hybrid` broadcast's 1.290 and nothing like the spread-policy 0.393, which is one more
demonstration that a single transfer factor does not exist.

### R2 · a14 → a15, the CD4 residual — **passes at 0.92×, the strongest retrodiction here**

Board: `pds` **+0.007257 scaled**. This experiment reconstructs a14's exact config
(`mass=None, cd4_frac=0.05`) and reads PLOLO's four folds:

| route | scaled `pds` | ratio to board |
|---|---:|---:|
| exact a14 arm vs exact a15 arm | **+0.006655** (raw +0.003009, fold sd 0.001083, n = 4) | **0.92** |
| the same move as a leave-one-out of the shipped stack (`cd4_a14`) | +0.006108 | 0.84 |
| removing CD4 outright (`noCD4`) | +0.015496 | — (larger, as it must be: a14 already carried 0.05) |

PLOLO recovers a board `pds` move to within 8% on a single-variable pair. That is the instrument
working as 162 advertised, and it is why every `pds` number below is quoted from it.

### R3 · a15 → a19, the whitener ridge — **right sign, overstated 1.37×, and it calibrates the whole table**

Board: overall **−0.013518**, recorded against `mse`. The 096 claim route gives raw `mse`
1.01377 → score 0.00000 → **−0.018520**. a15's *entire* `mse` score is +0.020009 overall, so the
claim model says a19 lost **all** of `mse`; the board says it lost **68%** of it.

**Consequence, and it applies to six rows of the main table.** The claim route has a hard cliff
where the board has a soft one. Every `+0.018520` in the leave-one-out table is a **censored
lower bound in magnitude, overstated ≈1.37×** by the one board point that tests it; read them as
≈ **+0.01352**. Applying that correction changes no verdict — all six stay far above the gate —
but it changes the sizes, and it means the six `+0.018520` rows are not distinguishable from each
other on `mse`.

### R4 · a2-prune → a8, `--overdispersion 4.0` — **passes only in a8's own context**

Board: `reach` **+0.021398 scaled = +0.018900 raw**. a8 carried a **flat 0.25** spread and no
tier, so the retrodiction has to be run there:

| context | h1-test d`reach` raw | t | sign vs board |
|---|---:|---:|---|
| **a8's flat 0.25** | **+0.003932** | +1.53 | **recovered** |
| a15's tiered 8000 | **−0.026098** | −7.49 | **reversed** |

Implied h1 → board transfer in the flat context: **4.81** (h1-test understates by ~5×). In the
tiered context no transfer factor of any sign fits.

**`--overdispersion`'s `reach` value is a property of the (spread, φ) pair, not of φ.** The
board evidence that put φ = 4 into the build was collected in a context we stopped shipping at
a12. This does not make φ = 4 wrong — it holds the whole `mse` score, see below — but it does
mean its `reach` credential expired and nobody noticed.

---

## Phase 1 — leave-one-out, all six members, overall board units

Sign convention throughout: **positive = what removing the component would cost us.**
Instruments per member exactly as the mandate specified: `pds` on PLOLO (4 folds, challenge
frame, pre-committed unweighted mean); `mse` through 096's `claim` route on the real challenge
panel; `reach`/`fid`/`jac`/`nmae` on h1-test, emitted, matched seeds, significance-selected pool,
donors `k562_gwps` + `orion_hek293t` only. Full table `work/rank200.json`, log `work/rank200.log`.

| component | `pds` | `mse` | `reach` | `fid` | `jac` | `nmae` | **SUM** | SUM, `jac(np)` | n | class |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| union whitening operator (γ 0.65) | +0.025593 | +0.018520 | +0.000777 | −0.002407 | −0.000355 | +0.001835 | **+0.043963** | +0.043805 | 6 | shape |
| equal donor weights (1,1,1) | +0.001312 | +0.018520 | −0.000408 | +0.000417 | +0.000401 | +0.009895 | **+0.030137** | +0.029655 | 6 | shape |
| `--overdispersion 4.0` | 0 | +0.018520 | −0.001931 | +0.000603 | +0.005824 | +0.000043 | **+0.023058** | +0.017087 | 6 | shape |
| tiered spread policy (8000 / 1.0) | 0 | +0.018520 | −0.002119 | −0.002743 | +0.004863 | +0.000036 | **+0.018557** | +0.011270 | 6 | shape |
| `VCC_ALPHA` 0.982 (vs v11's 1.0796) | +0.000022 | +0.018520 | — | — | — | — | **+0.018542** | +0.018542 | 2 | amp |
| CPM gate (`--min-cpm 10`) | +0.001169 | +0.018520 | +0.000562 | +0.000240 | −0.001523 | −0.000571 | **+0.018397** | +0.020258 | 6 | shape |
| post-whitening CD4 residual (f 0.10) | +0.002611 | +0.000024 | — | — | — | — | **+0.002635** | +0.002635 | 2 | amp |
| a23 mass fix (`minnorm`) | +0.000676 | +0.000073 | −0.000295 | +0.000107 | +0.000045 | +0.000022 | **+0.000628** | +0.000607 | 6 | shape |
| pruned basis (`union_pruned`) | +0.000392 | −0.000018 | — | — | — | — | **+0.000373** | +0.000373 | 2 | shape |
| PerturbMatch K562 replacement (0.5) | +0.001885 | −0.000256 | −0.001104 | +0.000075 | +0.000299 | −0.002514 | **−0.001615** | −0.001639 | 6 | shape |

### The knob neighbourhood — the shipped *setting* against its nearest measured alternative

| knob | `pds` | `mse` | `reach` | `fid` | `jac` | `nmae` | **SUM** | SUM, `jac(np)` |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| donor weights 1:1:1 vs 2:1:1 | −0.000126 | +0.015878 | +0.000261 | −0.000241 | −0.000143 | +0.004448 | **+0.020077** | +0.019761 |
| tier 8000 vs tier 14000 | 0 | +0.018520 | — | — | — | — | **+0.018520** | +0.018520 |
| γ 0.65 vs 0.50 | +0.002283 | +0.005252 | −0.000151 | −0.000361 | +0.000015 | −0.000611 | **+0.006427** | +0.006469 |
| CPM gate 10 vs 20 | +0.000684 | −0.000179 | +0.000815 | −0.000124 | +0.000632 | +0.002515 | **+0.004344** | +0.003631 |
| CPM gate 10 vs 5 (088's move) | −0.000139 | +0.004807 | +0.001183 | +0.000219 | −0.001636 | −0.000592 | **+0.003842** | +0.005709 |
| φ 4 vs φ 8 | 0 | 0 | — | — | — | — | **0** | 0 |
| **γ 0.65 vs 0.80** | −0.000826 | −0.000137 | −0.000697 | +0.000390 | +0.000136 | +0.000628 | **−0.000506** | −0.000691 |
| **tier 8000 vs tier 4000** | 0 | 0 | −0.000250 | +0.001080 | −0.004618 | +0.000024 | **−0.003765** | **+0.002196** |
| **tier 8000 vs FLAT 1.00** | 0 | 0 | −0.000344 | +0.001387 | −0.006590 | −0.000034 | **−0.005580** | **+0.002585** |

### Add-one-in from a bare consensus — the interaction sequential addition could not see

`pds` + `mse` only; there is no bare-consensus DE family, and the h1-test `reach` sem
(0.003–0.007 raw) is larger than every solo value in the table. `work/aoi200.json`.

| component | LOO `pds` | LOO `mse` | **LOO sum** | AOI `pds` | AOI `mse` | **AOI sum** | **interaction** |
|---|---:|---:|---:|---:|---:|---:|---:|
| union whitening (γ 0.65) | +0.025593 | +0.018520 | +0.044113 | +0.023557 | +0.018717 | +0.042274 | +0.001839 |
| PerturbMatch blend 0.5 | +0.001885 | −0.000256 | +0.001630 | +0.000988 | 0 | +0.000988 | +0.000642 |
| post-whitening CD4 (f 0.10) | +0.002611 | +0.000024 | +0.002635 | +0.002828 | 0 | +0.002828 | −0.000193 |
| **a23 mass fix** | +0.000676 | +0.000073 | +0.000749 | **−0.001241** | 0 | **−0.001241** | **+0.001990** |

The mass fix is the interesting row: **harmful on its own (−0.001241) and mildly helpful in the
full stack (+0.000749)**, an interaction of +0.001990 that is larger than either reading. That is
exactly the object sequential addition is blind to, and it is why "each component was validated
when it was added" was never the same statement as "each component earns its place".

---

## Phase 2 — the ranked verdict

Gate: rule 13's composed slot gate, ±0.0019 shape-only, ±0.0089 amplitude-changing.

| rank | component | six-member sum | verdict |
|---:|---|---:|---|
| 1 | union whitening operator (γ 0.65) | +0.043963 | **LOAD-BEARING** |
| 2 | equal donor weights (1,1,1) | +0.030137 | **LOAD-BEARING** |
| 3 | `--overdispersion 4.0` | +0.023058 | **LOAD-BEARING** |
| 4 | tiered spread policy | +0.018557 | **LOAD-BEARING** |
| 5 | `VCC_ALPHA` 0.982 | +0.018542 | **LOAD-BEARING** |
| 6 | CPM gate `--min-cpm 10` | +0.018397 | **LOAD-BEARING** |
| 7 | post-whitening CD4 residual | +0.002635 | noise (gate ±0.0089) |
| 8 | a23 mass fix | +0.000628 | noise |
| 9 | pruned basis | +0.000373 | noise |
| 10 | PerturbMatch K562 replacement | **−0.001615** | noise (negative, inside the gate) |

**No component is net negative outside its gate. Nothing should be dropped.**

**All six load-bearing rows carry the same `mse` term, and it is censored.** Removing any one of
them drops raw `mse` over the scoring cliff and zeroes the whole `mse` score: +0.018520, which
R3 measures as overstated ≈1.37× (read ≈ +0.01352). Ranks 1 and 2 are separated from the rest by
what they hold on *other* members — whitening by +0.025593 on `pds` (the largest single number in
this experiment, PLOLO signs `----` across all four folds), equal donor weights by +0.009895 on
`nmae` (`reb_w41` `nmae` 1.00571 against `reb_ship`'s 0.97006). Ranks 3–6 are tied inside the
`mse` instrument alone and this experiment does not separate them; what it does establish is that
each of the six is individually sufficient to lose ≈0.0135–0.0185 overall.

### The two the mandate flagged

**a23's mass fix — it never cleared, and `reach` is not why.** Board −0.00073 overall, this
experiment +0.00063, gate ±0.0019. Its `pds` credential is real and exactly reproduced
(+0.004014 scaled), its board `reach` cost is real (−0.007933 scaled), and the two very nearly
cancel. The h1-test `reach` reading is **−0.003989 raw at t = −1.50** — not resolvable — so this
experiment cannot confirm the board's `reach` cost, only that the *total* is a wash on both
instruments. **Removing it is not free score; it is a coin flip worth ±0.0007.** Leave it, and
stop citing it as a `reach` culprit.

**The tiered spread policy — the `mse` half is real, the `reach` half is a bad deal, and there is
a strictly better setting available.**

* *`mse` half.* Real and large. Dropping the tier for a flat **0.25** puts raw `mse` at 1.17174,
  score 0 — the whole `mse` score, ≈ −0.0135 after R3's correction. Flat 0.25 is not an option.
* *`reach` half.* The tier costs **−0.02863 raw h1-test `reach` at t = −7.54** against flat 0.25,
  and 171 §12c's mechanism (signal block parked at minimum-|z|, noise tail leading the walk) is
  confirmed by the sign. It also costs `fid` (−0.002743).
* *But the two halves are separable.* `mse` does not need the *tier*; it needs the spread to be
  off 0.25. **Flat 1.00 at φ 4 holds `claim` at 0.24653 on the `s²B` branch with raw `mse`
  0.88228 — identical to a15's — and holds `pds` by the per-gene-total identity.** So
  tier-vs-flat-1.00 is a pure four-member trade at zero `mse` and zero `pds` cost:

  | | d`reach` raw | `reach` | `fid` | `jac` (fixed NP) | `jac` (scaled NP) | `nmae` | SUM |
  |---|---:|---:|---:|---:|---:|---:|---:|
  | tier 8000 − flat 1.00 | −0.00465 (t −1.33) | −0.000344 | +0.001387 | −0.006590 | +0.001576 | −0.000034 | **−0.005580 / +0.002585** |

  **The whole verdict sits in `jac`, and `jac` is the member this harness models rather than
  measures.** Flat 1.00 raises median h1-test `n_pred` 6059 → 9466 (+56%) and recall 0.665 →
  0.926. 169's convention holds the board's `n_pred` fixed at 4,424 while recall moves, which
  makes the call-set expansion look free; scaling the board's `n_pred` by the measured h1-test
  ratio makes it look expensive. The only board point that discriminates is **a21**
  (`resample_x 0.14`, h1-test `n_pred` 6,034 → 1,913 at recall 0.669 → 0.271, board `jac`
  **−0.017 scaled**): fixed-NP predicts −0.0586 (right sign, 3.4× large), scaled-NP predicts
  **+0.0150 (wrong sign)**. That favours the incumbent fixed-NP convention — and therefore
  favours **dropping the tier for flat 1.00, worth ≈ +0.0056 overall** — but a21 also collapsed
  its call set rather than expanding it, so it is the wrong-direction test.
  **Unresolved, and handed to 192**, which owns the call set. Do not ship this on the strength of
  one reversed-direction retrodiction.

  For the record, `tier 4000` reads the same way and slightly smaller (−0.003765 / +0.002196),
  and it is the other `mse`-safe option: `claim` 0.24653, same branch, same raw `mse`.

### The four that are noise, and why that is not a criticism

CD4 (+0.002635), the mass fix (+0.000628), the pruned basis (+0.000373) and PerturbMatch
(−0.001615) are all inside their gates. Each of them was a correct call when it was made — CD4
retrodicts its board `pds` point at 0.92× (R2), PerturbMatch's `pds` LOO is +0.001885 with a
per-fold signs `--0-` (the `0` is the fold with no PerturbMatch coverage). They are simply small
now that the stack around them is bigger.
**Removing any of them buys nothing measurable and risks the `pds` they do hold**, so the correct
action is to leave them and stop spending slots re-litigating them.

---

## Phase 3 — falsifiers

`work/falsify200.log`, `work/falsify200.json`.

**1. Untouched members stay where they belong.** Per-gene panel totals were hashed for nine arm
configurations. Seven distinct digests over nine configurations, and the grouping is exactly what
the mechanism predicts:

| digest | configurations |
|---|---|
| `3599f15fe0cc40cf` | artifact totals (164's path) |
| `9a0b10f195f449f8` | artifact, mass fix ON |
| `20b289f756c25698` | artifact, mass fix OFF |
| `5b4c83c005ecb942` | rebuild SHIP, **rebuild noCD4** |
| `21c4d7f3128239d6` | rebuild 095 geometry, **rebuild noPM** |
| `491fc47a7a4246d6` | rebuild noWhiten |
| `1292aedb100cc654` | rebuild SHIP, gate 10 |

Every ablation that actually moves the profile breaks the digest — the mass fix, whitening, the
gate — which is the falsifier passing. **The two collisions are both real and both explained**:
`noCD4` = SHIP and `noPM` = 095 geometry because the CD4 residual is a structural no-op on this
harness (see *What could not be measured*), verified independently at the sha256 level on all
three seeds. The check therefore reports a true negative, not a blind spot.

The emission-law ablations (tier, φ) do not appear in the digest table because they cannot move
it: `realize200.py`:263 asserts `B.sum(0) == totals` on **every emitted block — 85 targets ×
68 arms = 5,780 assertions**, so the per-gene-total identity is enforced rather than sampled.
That is also why their `pds` and `mse` entries are exactly 0 by construction rather than by
measurement, and it is the same identity the `rc_hybrid` board point confirmed to sixteen
significant figures.

**2. Rule 10 config diff.** a15 → a23 differ in **exactly one** field, `mass_fix: null →
"minnorm"`; with that field removed the two configs are identical (`a15_a23_ablated_diff = {}`).
Operator sha256 identical (`e2f65ecf69b7e0b0`). `effects.whitening.after_all_median_norm`
4.2292956044543 vs 4.229295604454299, relative drift **4.4e-16** — no drift.

*Caveat recorded rather than hidden:* the two artifacts' `post_whiten_sha256` **differ**
(`51a640a84c2e` vs `b2426e15d76d`) at a config diff of one render-time flag. The eigen-multiply
is not bit-reproducible across runs while the audit numbers agree to 1e-15, so **any sha-level
comparison of two independently built effect matrices will fail for this reason alone** and must
not be used as a rule-10 check.

**3. Retrodictions.** All four above; R1, R2 and R4 pass on sign, R2 to 0.92× on magnitude, R3
right-signed and 1.37× overstated with the correction propagated into the main table.

---

## What could not be measured, and who should take it

**The CD4 residual has no measurable DE members on h1-test, and the reason is a real property of
the shipped builder.** `de200.reb_noCD4.npz` is **byte-identical** to `de200.reb_ship.npz`
(sha256 `aa28e8495ab66101…`, all three seeds), and `de200.reb_noPM.npz` to `de200.reb_e095.npz`
(`258ed003d3e3668a…`). The cause:

> production's `solve_cd4_fraction` sizes the residual by the **median** row norm of the added
> block (`plolo.py`:363–368: `rn = median(‖added‖); per_unit = rn/bn; frac = target_f/per_unit`).
> h1-test has CD4 coverage on **4 of 85** targets, so that median is exactly **0**, `per_unit` is
> 0, the guarded branch returns `frac = 0.0`, and the residual is added at strength zero.

On the challenge panel coverage is **276 of 300**, the median is non-zero, and the residual is
real — `cd4_fraction_solved` 0.094 in a15's own meta, and R2 retrodicts its board `pds` point at
0.92×. So this is an instrument limitation, not a build defect. **But it is a silent one:** any
future panel or harness with under half its targets covered gets a zero-strength CD4 residual and
no warning. Worth a one-line guard in the builder; flagged to whoever next touches
`solve_cd4_fraction`. Even with a fraction forced, h1-test could not resolve it — 4 movable rows
of 85 bound |Δ`reach`| below the 5-seed sem — so CD4's DE members are unmeasurable on this
harness by construction, and its value here is `pds` + `mse` only.

**Two consequences of that identity for the rebuild family.** `reb_e095` was intended as the
rebuild family's independent known-answer row and it is not independent: it collapses onto
`reb_noPM`. And the rebuild's median per-row cosine to the artifact's recovered `H` is **0.813**
(095 geometry 0.885), so the geometry components' DE members are measured on a *reconstruction*,
not on the shipped artifact. Their `pds` and `mse` come from the real challenge-frame path and
are not affected. Geometry DE deltas are also measured at the builder's default `--min-cpm 5.0`
(`realize200.py`:112), not the shipped 10.0; the gate is held equal within each pair, so the
deltas stand, but the shipped gate's own leave-one-out is the `gate10 vs gate0` pair used above.

**Handoffs.**

* **192 (call set)** — owns the decision above. The `tier 8000 → flat 1.00` swap is worth
  **+0.0056 overall if `jac`'s board `n_pred` is fixed and −0.0026 if it scales**, at zero `mse`
  and zero `pds` cost, with three matched seeds already emitted
  (`de200.art_flat100_mf{,_s2,_s3}.npz`). Settling which way board `jac` responds to a +56%
  call-set expansion resolves it. This is the cheapest live arm this experiment found.
* **196 (emitted depth)** — `φ 4 vs φ 8` is exactly zero on `pds` and `mse` (`claim` 0.24653,
  same branch) and unmeasured on the DE members; R4 says φ's `reach` value is context-dependent,
  so the (spread, φ) pair should be swept jointly rather than one knob at a time.
* **193 (operator eigenvalue map)** — γ 0.65 vs 0.80 is **−0.000506**, i.e. 0.80 is very slightly
  ahead and the difference is a quarter of the gate. 170's "0.65 is a six-member local optimum"
  survives as *not distinguishable from 0.80*, not as a measured optimum. γ 0.50 is clearly worse
  (+0.006427 for 0.65).
* **194 (donor tables)** — equal weights are the single largest `nmae` component in the table
  (+0.009895; `reb_w41` `nmae` 1.00571 vs `reb_ship` 0.97006) and 2:1 K562-heavy is worth
  −0.020077, so 175's proposed K562-heavy rider is priced and it is negative.

---

## Reproduction

```
uv run python inputs200.py        # rule-11 sha record          -> work/inputs200.json
uv run python folds200.py         # PLOLO fold prep             -> work/parts-*.npz
uv run python pds200.py           # 26 arms x 4 folds x 3 frames-> work/pds200.json
uv run python report_pds.py       # LOO / AOI on pds            -> work/report_pds.json
uv run python chalparts200.py     # challenge-panel parts       -> work/chalparts200.json
uv run python mse200.py           # 096 claim route, 34 arms    -> work/mse200.json
uv run python h1parts200.py       # h1-test rebuild parts       -> work/parts-h1test.npz
uv run python realize200.py       # 68 emitted arms, 3602+4164s -> de200.*.npz
uv run python six200.py           # 4 DE members, 78 tables     -> work/six200-tablev2.json
uv run python aoi200.py           # add-one-in vs leave-one-out -> work/aoi200.json
uv run python retro200.py         # the four board retrodictions-> work/retro200.json
uv run python rank200.py          # the six-member ranking      -> work/rank200.json
uv run python falsify200.py       # the three falsifiers        -> work/falsify200.json
```

No `vcc`, no `cell-eval2`, no submission build, no writes outside this directory. The heaviest
step is the one that loads the h1-test operator's 17,745² eigenbasis; `/usr/bin/time -l` measures
its maximum resident set size at **4.28 GB**, against the 12 GB budget.
