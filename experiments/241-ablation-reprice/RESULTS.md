# 241 — 200's ten-component ablation, re-priced on the renderer we actually ship

**The ordering survives. The two verdicts 219 flagged do not.**

200's table was measured with `VCC_MASS_FIX=minnorm`, a renderer only a23 ever shipped. Re-run
with it unset — the shipped configuration, and the only one the shipped `rc_callset` emitter
will even accept — **no component changes rank by more than one place, and no rank change is
larger than its own error bar.** 200's "six load-bearing, four inside their gate" ordering is
therefore **safe to quote again**, which it has not been since 219 §4 said so a month ago.

Three things did change, and all three are decision-relevant:

1. **PerturbMatch crosses out of its gate.** −0.001615 → **−0.002178** on the shipped renderer
   at the shipped dose, and **−0.003540** once the shipped `rc_callset` emitter is included
   too — **2.19×** 200's number, **1.86×** outside the ±0.0019 shape gate. It is no longer
   "negative, inside the gate"; it is **negative, outside the gate, and still not shippable** —
   because split by channel (§8.1b rule 3) **neither half clears its own bar**: DIRECTION
   +0.001403 (0.74× the ±0.0019 gate) and ORDERING +0.002137 (**0.47 sd** of the 0.004502
   ordering noise).
2. **The CD4 residual's understatement is confirmed and sized.** Its four DE members are a
   *structural* zero on h1-test (coverage 4/85 → solved fraction exactly 0 →
   `de241.reb_noCD4_nomf.npz` **bit-identical** to `de241.reb_ship_nomf.npz`). Filling them
   from the row's own PLOLO cosine move gives **+0.0055 … +0.0065** at the shipped dose against
   200's published **+0.002635** — **2.07× … 2.46× understated**, direction **UP**, exactly
   §8.4's Error 1.
3. **The censored `mse` leg is 1.616× bigger than 200 could see.** 200 read it off a15's
   `score_mse` 0.11112 → +0.018520. The stack we ship reads **0.17957** on the board, so a row
   whose removal crosses the `mse` cliff is worth **+0.029925**, not +0.018520.

One row also flips a member's sign: **the CPM gate's `pds` leg goes +0.001169 → −0.000804.**
It stays load-bearing because its value was never in `pds`; it was always the `mse` cliff.

---

## Phase 0 — retrodiction. Every published row of 200 reproduces, and two more besides

**Nothing below was re-derived; 200's own instruments were imported and re-run, and where an
artifact exists on disk it was loaded rather than rebuilt.** Rule 11 record:
`work/inputs241.json` — **7 of 7** operator/frozen shas match 200's record byte for byte
(`e2f65ecf…` production, `13d988fc…` / `7cecbf7f…` / `97e486e0…` held, `f80fe8d6…` h1-test,
harness `06c6baea8c62ef3d`, 169's effect `e4d8af8b892cb22d`). The one input that HAS drifted
is `scripts/make_submission_pm_cd4.py` (`24a30055…` against 200's `90c4f49e…`; 240 owns it) —
241 never calls its mass-fix path (`mass_fix` is `False` on every 241 emission), and the one
place it *is* used, `mse241.py`'s `MF` base, reproduced 200's `A`/`B`/`claim` exactly, which is
the evidence the drift is not in `mass_consistent_effect`. The module's `MASS_FIX` default is
asserted to be `None` at import (`emit241.py`), which is the shipped configuration.

### P0.1 · The ten-row table itself — **10 of 10 exact**

| component | 241 (200's renderer) | 200 published | \|Δ\| |
|---|---:|---:|---:|
| union whitening operator (γ 0.65) | **+0.043963** | +0.043963 | < 1e−6 |
| equal donor weights (1,1,1) | **+0.030137** | +0.030137 | < 1e−6 |
| `--overdispersion 4.0` | **+0.023059** | +0.023058 | 1e−6 |
| tiered spread policy (8000 / 1.0) | **+0.018557** | +0.018557 | < 1e−6 |
| `VCC_ALPHA` 0.982 | **+0.018542** | +0.018542 | < 1e−6 |
| CPM gate (`--min-cpm 10`) | **+0.018397** | +0.018397 | < 1e−6 |
| post-whitening CD4 residual | **+0.002635** | +0.002635 | < 1e−6 |
| a23 mass fix (`minnorm`) | **+0.000628** | +0.000628 | < 1e−6 |
| pruned basis (`union_pruned`) | **+0.000373** | +0.000373 | < 1e−6 |
| **PerturbMatch K562 replacement (0.5)** | **−0.001615** | **−0.001615** | < 1e−6 |

The four rows the mandate named are all in it. **PASS.**

### P0.2 · Per-instrument gates

| gate | got | expected | |
|---|---:|---:|---|
| `pds` — `MF:BASE` absolute, 4 folds (hct116 / hek293t / k562 / kolf2) | 0.791823 / 0.777375 / 0.789121 / 0.610128 | 200 §Phase 0 | **PASS** (worst 3.4e−07) |
| `pds` — `BASE` vs `BASE-dup`, all three bases | worst \|Δ\| **0.00e+00** | 0 | **PASS** |
| `pds` — a15 → a23 mass fix, scaled | **−0.004014** | −0.004014 | **PASS** (2.2e−07) |
| `pds` — 200's whole LOO column, 8 rows | worst \|Δ\| **4.9e−07** | 0 | **PASS** |
| DE — 164/157's a15, `reach` / median `n_pred` | **0.18283 / 6034** | 0.18283 / 6034 | **PASS** |
| DE — 200's four legs × 7 rows, re-derived independently of `six200.py` | worst \|Δ\| **4.6e−07** | 0 | **PASS** |
| `mse` — 164 §1's four `claim` gates | 0.11860 / 0.19438 / 0.16652 / 0.24581 | ≤ 3% | **PASS** |
| `mse` — A/B round trip at s = 1 | +0.000% / +0.000% | board_solve | **PASS** |
| `mse` — a15 raw / `score_mse` | **0.88228 / 0.11112** | 0.88228 / 0.11112 | **PASS** |
| `mse` — 200's stored `A`/`B`/`claim` for 7 MF arms | identical to 5 d.p. on all 21 numbers | 200 | **PASS** |
| emission — `de241.art_nomf_KA.npz` vs `de200.art_nomf.npz` | **BIT-IDENTICAL** on all 7 arrays | — | **PASS** |
| rebuild — 5 geometry row-cosines to the artifact's `H` | +0.81271 / +0.79974 / +0.88515 / +0.70662 | 200's log | **PASS** |

### P0.3 · Two *unplanned* known-answer rows, both passed

**(a) 214 measured two of these rows on the shipped renderer first.** 214 §1.1/§1.2 swept
PerturbMatch and the CD4 dose on exactly this experiment's `NOMF` and `SHIPPED` bases, in
independently written code. All five rows reproduce to < 5e−07 (`work/ka214_241.json`):

| row | 241 raw Δ`pds` | 214 | 241 overall | 214 |
|---|---:|---:|---:|---:|
| PerturbMatch, a23 renderer, f = 0.10 | +0.005058 | +0.005058 | +0.001885 | +0.001885 |
| PerturbMatch, **shipped** renderer, f = 0.10 | +0.002856 | +0.002856 | +0.001065 | +0.001065 |
| PerturbMatch, **shipped** renderer, f = 0.34 | +0.003797 | +0.003797 | +0.001415 | +0.001415 |
| CD4 residual at pm 0.50, f = 0.10 | +0.007302 | +0.007302 | +0.002722 | — |
| CD4 residual at pm 0.50, f = 0.34 | +0.007260 | +0.007260 | +0.002706 | — |

219's *"PerturbMatch reads 1.77× smaller on the shipped renderer"* is **exactly 1.77×**.

**(b) §4.2's whitening/cosine inversion.** On the shipped renderer, removing the operator
*raises* the 3-cancer-fold transfer cosine from **0.026516 to 0.036453 (+37.5%)** while costing
**+0.023685 overall on `pds`**. §4.2 measured −32% on a different fold set. Same sign, same
order — and it is the reason the rate-card fill below is **not** applied to the whitening row.

**(c) Internal consistency.** `NOMF:addMassFix` reproduces `MF:BASE` to the last printed digit
on both instruments (`pds` per fold, and `A` 0.22246 / `B` 0.22133 / `claim` 0.24588). Adding
the mass fix back to the shipped base lands on 200's base, as it must.

---

## Phase 1 — the ten-row table, before and after

Sign convention, 200's: **positive = what removing the component would cost us**, overall board
units. `n` = members priced. `Δcos %` is the component's own contribution to the **3-cancer**
PLOLO transfer cosine, in PLOLO units (falsifier 4 — see §Falsifiers). `e = ΔM/ΔU` is 188's
`nmae` identity split.

### BEFORE — 200's renderer (`VCC_MASS_FIX=minnorm`), CD4 `f = 0.10`

| # | component | `pds` | `mse` | `reach` | `fid` | `jac` | `nmae` | **SUM** | n | `pds` 3-cancer | `pds` kolf2 | Δcos % | `e` |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | union whitening operator (γ 0.65) | +0.025593 | +0.018520 | +0.000777 | −0.002407 | −0.000355 | +0.001835 | **+0.043963** | 6 | +0.030357 | +0.011300 | −39.25 | 0.293 |
| 2 | equal donor weights (1,1,1) | +0.001312 | +0.018520 | −0.000408 | +0.000417 | +0.000401 | +0.009895 | **+0.030137** | 6 | +0.001372 | +0.001133 | +6.56 | −0.743 |
| 3 | `--overdispersion 4.0` | 0 | +0.018520 | −0.001931 | +0.000603 | +0.005824 | +0.000043 | **+0.023059** | 6 | 0 | 0 | — | 0.827 |
| 4 | tiered spread policy (8000/1.0) | 0 | +0.018520 | −0.002119 | −0.002743 | +0.004863 | +0.000036 | **+0.018557** | 6 | 0 | 0 | — | 0.774 |
| 5 | `VCC_ALPHA` 0.982 | +0.000022 | +0.018520 | — | — | — | — | **+0.018542** | 2 | +0.000005 | +0.000074 | +0.10 | — |
| 6 | CPM gate (`--min-cpm 10`) | +0.001169 | +0.018520 | +0.000562 | +0.000240 | −0.001523 | −0.000571 | **+0.018397** | 6 | +0.002284 | −0.002178 | −0.84 | 0.535 |
| 7 | post-whitening CD4 residual | +0.002611 | +0.000024 | 0\* | 0\* | 0\* | 0\* | **+0.002635** | 2 | +0.002651 | +0.002491 | +9.15 | — |
| 8 | a23 mass fix (`minnorm`) | +0.000676 | +0.000073 | −0.000295 | +0.000107 | +0.000045 | +0.000022 | **+0.000628** | 6 | +0.001035 | −0.000400 | +0.54 | −0.144 |
| 9 | pruned basis (`union_pruned`) | +0.000392 | −0.000018 | — | — | — | — | **+0.000373** | 2 | — | +0.000392 | — | — |
| 10 | **PerturbMatch K562 replacement (0.5)** | +0.001885 | −0.000256 | −0.001104 | +0.000075 | +0.000299 | −0.002514 | **−0.001615** | 6 | +0.001487 | +0.003082 | +2.00 | 0.166 |

### AFTER — the shipped renderer (mass fix unset), CD4 `f = 0.10` — the single-variable flip

| # | component | `pds` | `mse` | `reach` | `fid` | `jac` | `nmae` | **SUM** | n | `pds` 3-cancer | `pds` kolf2 | Δcos % | `e` |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | union whitening operator (γ 0.65) | +0.023685 | +0.018520 | +0.000243 | −0.002217 | −0.000244 | +0.001563 | **+0.041551** | 6 | +0.028764 | +0.008448 | −37.48 | 0.316 |
| 2 | equal donor weights (1,1,1) | +0.001627 | +0.018520 | −0.000246 | +0.000117 | +0.000559 | +0.009978 | **+0.030555** | 6 | +0.001381 | +0.002364 | +6.56 | −0.742 |
| 3 | `--overdispersion 4.0` | 0 | +0.018520 | −0.001633 | +0.000549 | +0.005781 | +0.000051 | **+0.023267** | 6 | 0 | 0 | — | 0.710 |
| 4 | `VCC_ALPHA` 0.982 | +0.000123 | +0.018520 | — | — | — | — | **+0.018643** | 2 | +0.000096 | +0.000202 | +0.10 | — |
| 5 | tiered spread policy (8000/1.0) | 0 | +0.018520 | −0.002053 | −0.002745 | +0.004788 | +0.000043 | **+0.018554** | 6 | 0 | 0 | — | 0.638 |
| 6 | CPM gate (`--min-cpm 10`) | **−0.000804** | +0.018520 | +0.000196 | −0.000051 | −0.001480 | −0.000651 | **+0.015730** | 6 | −0.000087 | −0.002955 | −0.98 | 0.541 |
| 7 | post-whitening CD4 residual | +0.002722 | +0.000025 | 0\* | 0\* | 0\* | 0\* | **+0.002747** | 2 | +0.002504 | +0.003376 | +9.19 | — |
| 8 | a23 mass fix (`minnorm`) †| +0.000676 | +0.000073 | −0.000295 | +0.000107 | +0.000045 | +0.000022 | **+0.000628** | 6 | +0.001035 | −0.000400 | +0.54 | −0.144 |
| 9 | pruned basis (`union_pruned`) | +0.000355 | −0.000028 | — | — | — | — | **+0.000327** | 2 | — | +0.000355 | — | — |
| 10 | **PerturbMatch K562 replacement (0.5)** | +0.001065 | −0.000280 | −0.001195 | +0.000083 | +0.000338 | −0.002589 | **−0.002579** | 6 | +0.000342 | +0.003233 | +2.00 | 0.162 |

† the mass fix is **OFF** in the shipped stack, so this row is an **ADD-IN**: +0.000628 is what
adding it back would be worth, not what removing it would cost. (At the shipped CD4 dose it
collapses to **−0.000046**.)

### AFTER — shipped renderer **at the shipped dose `f = 0.34`**, `mse` on the board's own `score_mse` 0.17957

| # | component | `pds` | `mse` | `reach` | `fid` | `jac` | `nmae` | **SUM** | n | `pds` 3-cancer | `pds` kolf2 |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | union whitening operator (γ 0.65) | +0.019788 | +0.024793 | +0.000243 | −0.002217 | −0.000244 | +0.001563 | **+0.043928** | 6 | +0.024170 | +0.006642 |
| 2 | equal donor weights (1,1,1) | +0.002342 | +0.029925 | −0.000246 | +0.000117 | +0.000559 | +0.009978 | **+0.042675** | 6 | +0.002653 | +0.001407 |
| 3 | `--overdispersion 4.0` | 0 | +0.029925 | −0.001633 | +0.000549 | +0.005781 | +0.000051 | **+0.034672** | 6 | 0 | 0 |
| 4 | tiered spread policy (8000/1.0) | 0 | +0.029925 | −0.002053 | −0.002745 | +0.004788 | +0.000043 | **+0.029959** | 6 | 0 | 0 |
| 5 | `VCC_ALPHA` 0.982 | +0.000103 | +0.027457 | — | — | — | — | **+0.027559** | 2 | +0.000086 | +0.000154 |
| 6 | CPM gate (`--min-cpm 10`) | −0.000732 | +0.029925 | +0.000196 | −0.000051 | −0.001480 | −0.000651 | **+0.027206** | 6 | +0.000332 | −0.003927 |
| 7 | post-whitening CD4 residual | +0.002706 | +0.000143 | 0\* | 0\* | 0\* | 0\* | **+0.002849** | 2 | +0.001477 | +0.006391 |
| 8 | pruned basis (`union_pruned`) | +0.000206 | −0.000045 | — | — | — | — | **+0.000161** | 2 | — | +0.000206 |
| 9 | a23 mass fix (`minnorm`) † | +0.000007 | +0.000067 | −0.000295 | +0.000107 | +0.000045 | +0.000022 | **−0.000046** | 6 | +0.000602 | −0.001775 |
| 10 | **PerturbMatch K562 replacement (0.5)** | +0.001415 | −0.000229 | −0.001195 | +0.000083 | +0.000338 | −0.002589 | **−0.002178** | 6 | +0.000880 | +0.003022 |

`0*` = a **structural** zero, not a measurement (see §"What could not be measured"). `—` = the
member is not expressible on this instrument at all; it is left blank rather than carried at
zero, which is §8.4's rule.

### Which rows change rank

| component | 200 | shipped renderer, f = 0.10 | shipped renderer, f = 0.34 |
|---|---:|---:|---:|
| union whitening operator (γ 0.65) | 1 | 1 | 1 |
| equal donor weights (1,1,1) | 2 | 2 | 2 |
| `--overdispersion 4.0` | 3 | 3 | 3 |
| tiered spread policy | 4 | **5** | 4 |
| `VCC_ALPHA` 0.982 | 5 | **4** | 5 |
| CPM gate (`--min-cpm 10`) | 6 | 6 | 6 |
| post-whitening CD4 residual | 7 | 7 | 7 |
| a23 mass fix | 8 | 8 | **9** |
| pruned basis | 9 | 9 | **8** |
| PerturbMatch | 10 | 10 | 10 |

**Two swaps, and neither is a reordering.** `tier` ↔ `VCC_ALPHA` are **8.9e−05 apart**
(+0.018554 vs +0.018643) — both are the same censored `mse` leg plus a `pds` term smaller than
162's own 0.001494 σ. `pruned basis` ↔ `mass fix` are +0.000161 vs −0.000046, i.e. both **11×
inside** the ±0.0019 gate, and the pruned-basis row is a **single-fold** (kolf2-only) number.
**Report: the ordering is unchanged; 200's rank list is safe to quote.**

**What *did* change is magnitudes, and one of them is large.** On the stack we ship a row that
crosses the `mse` cliff is worth **+0.029925** rather than +0.018520 — **1.616×** — because our
`score_mse` went 0.11112 → **0.17957** between a15 and `h1union`. (200's R3 deflator 0.730 puts
the honest read at **+0.021845**.) Four of the six load-bearing rows are still fully censored at
the shipped dose (`equal weights`, `--overdispersion`, `tier`, `CPM gate`); the other two are
**not**, because the cliff is further away from raw 0.81726 than it was from 0.88228 —
`whitening` reads **+0.024793** and `VCC_ALPHA` **+0.027457**, both partial rather than total
losses of `mse`. Either way the six load-bearing components are individually worth **1.3–1.6×
more** than 200 could see.

---

## Phase 2 — the two rows `DEF-RUNBOOK` §3.1 says must be re-decided

### 2.1 PerturbMatch — **still net-negative, now outside its gate, and still not shippable**

| stack | DIRECTION (`pds`+`mse`+`nmae`) | ORDERING (`reach`+`fid`+`jac`) | six-member | six-member on 3 cancer folds | `pds` kolf2 leg |
|---|---:|---:|---:|---:|---:|
| 200's renderer, `f = 0.10`, `quasi` | −0.000885 | −0.000730 | **−0.001615** | −0.002013 | +0.003082 |
| shipped renderer, `f = 0.10`, `quasi` | −0.001804 | −0.000774 | **−0.002579** | −0.003301 | +0.003233 |
| shipped renderer, `f = 0.34`, `quasi` | −0.001403 | −0.000774 | **−0.002178** | −0.002713 | +0.003022 |
| shipped renderer, `f = 0.10`, **`rc_callset`** | −0.001804 | −0.002137 | **−0.003941** | −0.004664 | +0.003233 |
| **shipped renderer, `f = 0.34`, `rc_callset` — the full shipped configuration** | **−0.001403** | **−0.002137** | **−0.003540** | **−0.004076** | +0.003022 |

**Answer to the question as asked: yes, it is still net-negative, and on the complete shipped
configuration it is 1.86× (2.15× on cancer folds) OUTSIDE the ±0.0019 shape gate rather than
inside it.** 219 predicted "roughly −0.00244" from 214's 1.77× renderer factor; the renderer
alone gives **−0.002579**, and the emitter takes it to **−0.003540**.

**Reconciliation with 231, which refuted 214's Arm A.** They do **not** disagree, and the
reason is arithmetic rather than preference:

* 231 refuted Arm A's **+0.0088** because it implies a board raw `reach` of **0.31223** against a
  world best of **0.290835**. 241 measures the same arm at **+0.003540 … +0.004076** on the full
  shipped configuration — **2.2× smaller than the refuted headline**. Running 231's own
  admissibility test on 241's number: implied board raw `reach` = **0.28738 ≤ 0.290835** →
  **ADMISSIBLE, at 98.8% of the world record**. So 241's estimate is not the one 231 killed,
  and 231's refutation of the *magnitude* stands untouched — but 241's own number sits within
  0.0035 raw of the best `reach` on Earth, which is a fragility, not a comfort.
* **But 231's rule still binds, and it is why this is not a recommendation.** §8.1b rule 3 says
  price the halves separately. Doing so on the full shipped configuration: the **DIRECTION**
  half of dropping PerturbMatch is **+0.001403** — **0.74× the ±0.0019 gate**, inside it. The
  **ORDERING** half is **+0.002137** — **0.47 sd** of §8.1b's 0.004502 ordering noise, i.e.
  not distinguishable from a re-roll. **Neither half clears its own bar.** The composite only
  clears because two sub-gate halves add, and §8.1b rule 1 says price an arm in the units of
  the channel it lives in, not against a composite floor.
* **And falsifier 3 fires on this row harder than on any other.** PerturbMatch's entire `pds`
  credential is `kolf2`: **+0.003233 on the pluripotent fold against +0.000342 on the three
  cancer folds** (a 9.5× split). Read on cancer folds only — the folds that model D/E/F —
  dropping PerturbMatch is worth **+0.002713**, its largest value anywhere. The cancer/kolf2
  disagreement is the dominant uncertainty in this row, not the renderer.

**Mechanism, and it is §8.1d's law hitting a shipped component for the first time.** 188's
identity splits PerturbMatch's `nmae` cost exactly: it adds **ΔU = +0.01379** of magnitude to
the reference's DE pool and recovers only **ΔM = +0.00223** of it as sign-matched,
**`e = 0.162`** against a break-even of 0.5 — **84% of the magnitude it moves onto the scored
pool is unmatched.** That is why Δ`nmae` raw is **+0.00933** (worse) and why the row is
negative at all. **PerturbMatch is a magnitude arm.** §8.1d predicts exactly this outcome, and
this is the first time the rule has been applied to something already in the build.

**Verdict: do not drop it on this evidence, and do not defend it on 200's number either.**
It is a sub-gate negative in **both** channels even though the composite is 1.86× outside the
composite gate. Ship it, and stop citing "−0.001615, inside the floor" as the reason — the
reason is that the arm is too small to *measure*, channel by channel, and that its one
positive member is 9.5× concentrated in the fold that does not model D/E/F.

### 2.2 The CD4 residual — the correction's direction is **UP**, confirmed on the shipped renderer

| stack | `pds` | `mse` | 2-member (200's method) | Δcos, 3-cancer, rel % | DE-4 fill, `k = 1.000` | `k = 1.384` | **6-member range** |
|---|---:|---:|---:|---:|---:|---:|---:|
| 200's renderer, `f = 0.10` | +0.002611 | +0.000024 | **+0.002635** | +9.15 | +0.001067 | +0.001477 | +0.003702 … +0.004112 |
| shipped renderer, `f = 0.10` | +0.002722 | +0.000025 | **+0.002747** | +9.19 | +0.001066 | +0.001475 | +0.003813 … +0.004222 |
| shipped renderer, **`f = 0.34`** | +0.002706 | +0.000143 | **+0.002849** | **+19.90** | +0.002616 | +0.003621 | **+0.005466 … +0.006470** |

**The two-member price barely moves with the renderer (+0.002635 → +0.002747, +4%). The
understatement is not a renderer artefact — it is the four members that cannot be measured
here at all**, and they are 44–56% of the row.

* h1-test reads the four DE members at **exactly 0** for this row, and that is **structural**:
  CD4 covers **4 of 85** h1-test targets, production's `solve_cd4_fraction` sizes the residual
  by the **median** row norm of the added block, that median is 0, `per_unit` is 0, the guarded
  branch returns `frac = 0.0`, and the residual is added at strength zero.
  **`de241.reb_noCD4_nomf.npz` is bit-identical to `de241.reb_ship_nomf.npz`** (sha of
  `pred_p`‖`pred_padj`‖`pred_lfc`: `b6e5359b80c98bc8` both). 200 documented this; 241
  re-confirms it on the shipped renderer.
* Filling those four from the row's **own** PLOLO cosine move through 218's `cd4f034` rate card
  gives **+0.0055 … +0.0065** at the shipped dose = **2.07× … 2.46× 200's +0.002635**.
* **Independent corroboration from the board itself.** The same component's dose step
  `f 0.10 → 0.34` scored **+0.0130557**, and §8.4's decomposition shows the four members 200
  carried at zero delivered **+0.006103 of it, 47%**. 219's lower bound was "≥ 1.39×"; 241
  measures **2.07–2.46×** by a completely different route.

**So: 200's CD4 row is understated, the direction is UP, and its class changes.** At the shipped
dose the six-member estimate **+0.0055 … +0.0065** is still below the ±0.0089 amplitude gate a
CD4 dose change is judged against — so this does **not** promote it to "shippable"; it promotes
it from *"noise (gate ±0.0089)"* to **a real component whose price 200 could only see half of**.

---

## Phase 1b — the shipped EMITTER, not just the shipped renderer

200 emits with the `quasi` law; the shipped package emits with **`rc_callset`**. 214 §1.1
established that a23 / `rc_hybrid` / `rc_callset` share **one** effect matrix
(`effects.final_sha256` `22c09f22c749bda8`), so the emitter cannot move any component's `pds`
or `mse` — but it re-draws the **call set**, so it moves `fid`, `jac` and `reach`. Since
PerturbMatch's whole verdict lives in those plus `nmae`, that row was re-emitted through
`vcc2026.rc_emit.rc_hybrid_block` with 202's shipped `CallSet(400, 2200, 1.90)`.

**Falsifier, asserted per block rather than argued:** each of the 6 × 85 `rc_callset` blocks
asserts that its delivered per-gene panel totals equal `floor(n·L·p)` **exactly** — the same
totals the `quasi` arms delivered — so the `quasi` → `rc_callset` comparison is single-variable
by construction. Digests are identical across all three seeds within each config
(`ce142bf9402db9a8` for `ship`, `893764f3011959f6` for `noPM`) and differ between them, as
they must. Median `n_pred` lands at **2,618 / 2,620** against the prescribed call-set size of
**2,600** — the `rc_callset` signature. And **`nmae` comes back bit-identical to the `quasi`
arm** (−0.002589 both), which is the pseudobulk identity showing up in the one member that
reads only the mean profile.

**The emitter makes PerturbMatch's four DE members 1.40× worse again.**

| PerturbMatch, four DE members | `quasi` (200's law) | **`rc_callset` (shipped)** |
|---|---:|---:|
| `reach` | −0.001195 | **−0.002040** |
| `fid` | +0.000083 | +0.000059 |
| `jac` | +0.000338 | **−0.000157** |
| `nmae` | −0.002589 | −0.002589 (bit-identical) |
| **DE four** | **−0.003364** | **−0.004726** |

| six-member, PerturbMatch | renderer only | **renderer + shipped emitter** |
|---|---:|---:|
| shipped renderer, `f = 0.10` | −0.002579 | **−0.003941** (3-cancer −0.004664) |
| shipped renderer, **`f = 0.34`** — the full shipped configuration | −0.002178 | **−0.003540** (3-cancer −0.004076) |

**On the complete shipped configuration PerturbMatch reads −0.003540, which is 2.19× 200's
−0.001615 and 1.86× outside the ±0.0019 shape gate.** The channel split still does not clear:
DIRECTION **+0.001403** (0.74× the gate — and unchanged, because the emitter cannot move
`pds`, `mse` or `nmae`), ORDERING **+0.002137** (**0.47 sd** of the 0.004502 ordering noise).
231's admissibility test now reads **implied board raw `reach` 0.28738 against a world best of
0.290835 — admissible, but at 98.8% of the world record.** An arm that needs us to land within
0.0035 raw of the best `reach` team on Earth is fragile on its own terms, which is a second,
independent reason not to act on it.

---

## Falsifiers

**1 · Reproduce 200's own numbers first.** Done, and it is the largest block of this report:
**10 of 10** table rows, **8 of 8** `pds` LOO rows, **7 of 7** DE rows, all four `mse` gates,
a **bit-identical** emission table, and two unplanned cross-experiment retrodictions (214's
five shipped-renderer rows to < 5e−07, §4.2's whitening/cosine inversion in sign and order).

**2 · All six members on every row, with `nmae` on 188's identity and `e = ΔM/ΔU`.**
Every row that *can* carry six carries six. 188's identity `nmae = 1 + U/D − 2M/D` was verified
against `debench.Surface.nmae()` at **1.1e−16**. Per row, on the shipped renderer:

| component | Δ`nmae` raw | Δ`U/D` | Δ`M/D` | `e = ΔM/ΔU` | reads as |
|---|---:|---:|---:|---:|---|
| union whitening operator | −0.00562 | −0.01525 | −0.00481 | 0.316 | shrinks the pool magnitude; `nmae` improves because ΔU < 0 |
| equal donor weights | −0.03591 | −0.01445 | +0.01073 | −0.742 | the one row that buys `M` while shedding `U` — a **sign** gain |
| `--overdispersion 4.0` | −0.00006 | +0.00014 | +0.00010 | 0.710 | above break-even, and negligible either way |
| tiered spread policy | −0.00003 | +0.00012 | +0.00008 | 0.638 | above break-even, negligible |
| CPM gate (`--min-cpm 10`) | +0.00249 | −0.03042 | −0.01646 | 0.541 | at break-even |
| **PerturbMatch** | **+0.00933** | **+0.01379** | **+0.00223** | **0.162** | **magnitude — §8.1d's signature** |
| a23 mass fix | −0.00011 | −0.00009 | +0.00001 | −0.144 | too small to classify |
| post-whitening CD4 residual | 0 | 0 | 0 | — | **structural zero** |

*Reading `e` correctly:* the break-even is 0.5 **when ΔU > 0**. `nmae` improves iff
ΔU < 2ΔM, so with ΔU < 0 (whitening, weights, gate) a *small* `e` is good and with ΔU > 0
(PerturbMatch, φ, tier) a *large* `e` is good. PerturbMatch is the row where the sign of ΔU
and the size of `e` are both against it: **`e` = 0.162 at ΔU = +0.0138**, which is §8.1d's
"contribution is magnitude" case (its worked example reads 0.072).

**3 · `kolf2` separately.** Every `pds` and cosine number is reported as a 3-cancer-fold mean
**and** a kolf2 value, never merged. Three findings come only from that split:

* **PerturbMatch's `pds` is a kolf2 artefact**: +0.003233 kolf2 vs **+0.000342** cancer (9.5×).
* **The CPM gate's `pds` inverts across the split** on every base: cancer −0.000087, kolf2
  −0.002955 at `f = 0.10`; cancer **+0.000332**, kolf2 **−0.003927** at `f = 0.34`.
* **The pruned-basis row is kolf2-ONLY.** `union → union_pruned` is expressible only on the
  fold that carries the production operator, which is `kolf2`. 200's **+0.000373** — and this
  experiment's +0.000327 / +0.000161 — are **single-fold numbers from the pluripotent fold**,
  with **no cancer-fold measurement of that component in existence.** That is a bigger caveat
  on the row than the renderer ever was.

The four DE members have **no** kolf2 analogue: h1-test is one surface, and the h1 donor is
oracle leakage on it. Stated, not papered over.

**4 · No scalar PLOLO → board conversion.** Every cosine number in this report is in **PLOLO
units** and is tagged: functional `plolo.endpoints`'s macro (per-row then mean) cosine between
the **delivered** log-pseudobulk effect and the held library's own effect on that fold's scored
gene set; gate = the arm's own `--min-cpm` (10.0 except in the gate rows); frames = challenge
A/B/C, mean of three; fold set = 3 cancer folds headline, `kolf2` separate; tag = 241, base
MF / NOMF / SHIPPED as labelled. The only place a conversion appears is the CD4 DE-4 fill,
where it is given as a **range over `k_cos ∈ {1.000, 1.384}`** and flagged as unidentified
(231: 7.84× spread over four board pairs). It is never a headline and it is applied to exactly
one row.

---

## What could not be measured, and why each blank is a blank rather than a zero

| row × member | why | what would fix it |
|---|---|---|
| **CD4 residual × `reach`/`fid`/`jac`/`nmae`** | CD4 covers **4 of 85** h1-test targets → `solve_cd4_fraction`'s median row norm is 0 → residual added at strength 0. Verified by **bit-identity** of the two DE tables on both renderers. | a harness with ≥ ~50% CD4 target coverage, or a builder guard that warns instead of silently returning `frac = 0` (200 already flagged this; still unfixed) |
| **pruned basis × all four DE members, and × the three cancer folds on `pds`** | the `union → union_pruned` swap is expressible only on the fold carrying the production operator (`kolf2`); h1-test has one operator | held-out pruned operators per fold — a build, not a re-price |
| **`VCC_ALPHA` × all four DE members** | the h1-test rebuild renormalises every arm to the artifact's median row norm by construction, so amplitude is not an arm on that surface | 187 closed α on the two members that can see it; 211 closed it again |
| **every row × the `rc_callset` emitter, outside PerturbMatch** | cost. The emitter cannot move `pds`/`mse`/`nmae` (one shared effect matrix, 214 §1.1; and `nmae` came back **bit-identical** here, confirming it) but it does move `reach`/`fid`/`jac` — on PerturbMatch by **−0.001362**, a 40% worsening of that row | re-emit the rebuild family through `emitrc241.py` with `GEOMS=noWhiten,w41,…`; ~12 min per pair of configs |
| **every row × `VCC_H1_UNION_G = 0.28`** | `folds200`'s source vocabulary is the three donors; the H1 union residual is a fourth, added post-whitening. It shifts every row's base (and 218 measured +7.68% macro cosine for it), so the *levels* here are one component short of the shipped stack | port 218's `sources218`/`compose218` into the fold builder |

---

## Reproduction

```
uv run python inputs241.py        # rule-11 sha record, 7/7 match 200   -> work/inputs241.json
uv run python phase0_de241.py     # 200's DE legs, re-derived           -> work/phase0_de241.json
uv run python pds241.py           # 30 arms x 4 folds x 3 frames, 125 s -> work/pds241.json
uv run python report_pds241.py    # LOO on pds, 3 bases                 -> work/report_pds241.json
uv run python ka214_241.py        # 214's five shipped-renderer rows    -> work/ka214_241.json
uv run python cos241.py           # PLOLO transfer cosine per component -> work/cos241.json
ONLY=art_nomf_KA uv run python emit241.py     # rule-10 bit-identity row
uv run python emit241.py          # 30 arms, mass fix OFF, 690+1746 s   -> de241.*.npz
uv run python emitrc241.py        # PerturbMatch on `rc_callset`, 6 arms-> de241.rc_*.npz
uv run python six241.py           # the four DE members, before/after   -> work/six241.json
uv run python mse241.py           # 096 claim route, 3 bases, 1735 s    -> work/mse241.json
uv run python rank241.py          # the ten-row table x 3               -> work/rank241.json
uv run python phase2_241.py       # the two runbook rows                -> work/phase2_241.json
```

No `vcc`, no `cell-eval2`, no submission build, no production change, nothing written outside
this directory. Peak RSS by stage (`/usr/bin/time -l`): `pds241` **6.12 GB**, `emit241` art
**6.77 GB**, `emit241` geometry **8.10 GB**, `emitrc241` **6.53 GB**, `mse241` **5.65 GB**,
`six241`/`phase0_de241` **0.09 GB** — all under the 10 GB budget, run strictly serially, with
the memory gate re-checked before each heavy stage.

---

## Exactly which lines of `docs/DEF-RUNBOOK.md` §3.1 must change

**I have not edited the runbook.** The mandate quotes line 176; the file has moved twice since,
and while this experiment was running **someone already wrote a 241 update block into §3.1
(lines 194–206) off an intermediate version of this file.** So there are now *two* places to
fix, and the first is the new block itself.

**0 · Lines 195–200 — the `PerturbMatch` bullet of the new 241 block — carry `quasi`-emitter
numbers and must be updated to the full shipped configuration.** Three figures are superseded:

| in the runbook now | must read | why |
|---|---|---|
| six-member **−0.002178** (−0.002579 at f 0.10) | **−0.003540** (−0.003941 at f 0.10); cancer folds **−0.004076** | the row was priced with the `quasi` emission law; the shipped emitter is `rc_callset`, which moves `reach` −0.001195 → **−0.002040** and flips `jac` +0.000338 → **−0.000157** |
| ORDERING half **+0.000774 (0.17 sd)** | **+0.002137 (0.47 sd)** | same |
| implied board `reach` **0.28289** | **0.28738** — still ≤ 0.290835, but **98.8% of the world record** | same |

Everything else in that block — "SHIP IT", "retire the −0.001615 citation", "neither half clears
its own bar", "231 and 241 agree", the `kolf2` concentration, `e` = 0.162, the CD4 bullet, the
"ordering unchanged / +0.029925 censored leg" bullet — **stands unchanged and is confirmed by
the final numbers.** The direction half is emitter-invariant by construction (`pds`, `mse` and
`nmae` cannot move: one shared effect matrix, and `nmae` came back bit-identical), so the
conclusion does not change; the ordering half roughly triples and the composite worsens 1.6×.

**1 · Line 217, the `PerturbMatch` row of the §3.1 table itself, is still the OLD text** and
contradicts the block six lines above it. Currently:

> ⚠️ **200's ablation makes it the only NET-NEGATIVE component of the ten: −0.001615 six-member.**
> Inside the ±0.0019 floor as a standalone arm, which is why 200 declined to drop it.
> **Experiment 214 is testing whether it turned negative because 207 tripled the CD4 dose** —
> 064 measured the two as *synergistic* at the old dose. **If 214 finds them substitutes, this
> row and the CD4 row must both be re-decided before Oct 22.** Until then: ship it, as the
> table says. | 200, 214

Three things in it are now wrong or stale:
* **"−0.001615"** was measured on `VCC_MASS_FIX=minnorm` **with the `quasi` emitter** — a
  renderer the shipped `rc_callset` emitter literally refuses
  (`make_submission_pm_cd4.py`, the `emitter in RC_EMITTERS` guard — line 3481 at sha `24a30055…`). On the shipped renderer at the shipped dose it is
  **−0.002178**, and with the shipped emitter as well **−0.003540** (−0.004076 on cancer folds).
* **"Inside the ±0.0019 floor"** is **no longer true** — it is 1.86× outside as a composite.
  The correct reason not to drop it is the *channel split*: DIRECTION +0.001403 (0.74× the
  ±0.0019 gate) and ORDERING +0.002137 (0.47 sd of the 0.004502 ordering noise), i.e.
  **neither half clears its own bar** (§8.1b rule 3).
* **"Experiment 214 is testing whether…"** — 214 answered it (§1.2: the two are weak
  **complements** at every dose, interaction positive and growing, so the substitution
  hypothesis is **refuted**), and **231 then refuted 214's Arm A on arithmetic**. The sentence
  reads as an open question and is not one.

  Suggested replacement, same row: *"Ship it. 241 re-priced it on the shipped renderer AND the
  shipped emitter: **−0.003540** six-member at `f = 0.34` (−0.004076 on cancer folds), i.e.
  **1.86× outside** the ±0.0019 floor, but the halves are DIRECTION +0.001403 (0.74× gate) and
  ORDERING +0.002137 (0.47 sd) — **neither is shippable** (§8.1b rules 1 and 3), and 231's
  arithmetic bound admits 241's number (implied board raw `reach` 0.28738 ≤ world best
  0.290835, i.e. 98.8% of it) while still refuting 214's +0.0088. Its `pds` credential is 9.5×
  concentrated in `kolf2`; its `nmae` cost is a **magnitude** effect (`e` = 0.162,
  ΔU = +0.0138, §8.1d). Mechanical on the day; count pinned in step 8. | 200, 214, 231, 241"*

**2 · The sentence "this row and the CD4 row must both be re-decided before Oct 22" — the CD4
half is now decided.** 241 confirms the CD4 residual's price rises, not falls, on the shipped
renderer at the shipped dose (**+0.002849** two-member; **+0.0055 … +0.0065** six-member), and
214 §1.2 measured the CD4/PerturbMatch interaction **positive at every dose**. The
`VCC_CD4_TARGET_F` row (line 214) needs no change to what it ships; the clause **on line 217**
that makes it conditional on 214 should be deleted.

**3 · Add one line to §3.1, because it moves every future price on that page.** Our
`score_mse` is **0.17957** on `h1union`, not a15's 0.11112. **Any component whose removal
crosses the `mse` cliff is worth +0.029925 overall (+0.021845 after 200 R3's measured 0.730
deflator), not +0.018520.** Six of the ten are in that class. Anything on this page still
quoting +0.018520 as "the whole `mse` score" is 1.616× low. (The existing 241 block's third
bullet already says this; keep it.)

**Line numbers are as of this writing and the file is being edited concurrently — match on the
quoted text, not on the number.**

---

## Verdict

**200's ordering is confirmed and safe to quote again — no component moves more than one rank
and no rank change exceeds its own error bar — but two of its verdicts do not survive: on the
configuration we actually ship, PerturbMatch is net-negative 1.86× OUTSIDE the ±0.0019 gate
(−0.003540, 2.19× 200's number, though neither channel half clears its own bar so it still
ships), and the CD4 residual is understated 2.07–2.46× exactly as §8.4 predicted, because four
of its six members are structurally unmeasurable on h1-test rather than zero.**
