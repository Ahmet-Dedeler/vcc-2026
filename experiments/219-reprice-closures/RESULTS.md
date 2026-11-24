# 219 — Re-pricing every closure under §8.4's corrected rule. **§8.4's 4.5× is a double-count; the honest multiplier is 2.0–2.6× and it applies to a much narrower class than the section implies. Three closures re-open, none of them for the reason §8.4 gives.**

**Date:** 2026-09-04. Audit, in the tradition of 181 and 198. **No `vcc`, no `cell-eval2`, no
build, no submission, nothing written to `submissions/`, nothing written outside
`experiments/219-reprice-closures/`.** Peak RSS **4.85 GB** (`/usr/bin/time -l`, the one-off
PLOLO fold rebuild) against a 10 GB cap; every other script peaked at **0.07 GB**. All scripts
run under `uv run python`.

**40 known-answer gates, 0 fail** (`phase0_219.py`).

---

## Headline, in six lines

1. **§8.4's "≈ 4.5×" is not reproducible, and the reason is a double-count.** Its 2.51× is
   `realised total ÷ PREDICTED mse`, and 207's predicted `mse` had *already* passed the PLOLO
   cosine through a board-calibrated transfer factor `k_cos = 1.2598`. Multiplying by Error 2's
   1.80× charges the transfer error twice. **The end-to-end multiplier for an arm priced 207's
   way is 2.51×, point; range 2.00–2.56.**
2. **Error 2's 1.80× is also inflated by a functional mismatch.** It divides the board's
   *aggregate*-cosine ratio by 212's *macro-per-fold* cosine ratio. On the functional 197
   actually inverts (207's `cos_agg`) the factor is **1.38×**, and it agrees with the
   independent `a14 → a15` board calibration of **1.26×**. Range **1.17–1.39**.
3. **The correction does not rest on one board point.** `a14 → a15` is a second board
   measurement of the same knob (a 1.89× CD4 dose step) and it gives the propagation factor
   independently: **2.33×** against `cd4f034`'s **1.997×**. The two bracket it.
4. **And it does NOT apply to operator-geometry arms.** 212 §3.4 measured `d(pds)/d(cos) =
   −0.3517` against `d(mse)/d(cos) = +0.3655` along the gamma axis. Carrying the DE members
   through at the two board points' own rates, an operator-geometry arm's multiplier is
   **0.59–1.23**, not 2.51. **That kills the re-opening of 170, 193, 195, 212 and 209 on
   measured rates rather than on "the negative is large".**
5. **The one general number worth more than any individual re-price: §8.2's "a direction arm
   must deliver +31% cosine" is an `mse`-only construction and is 2.0–2.3× too demanding.**
   Corrected against the amplitude floor the bar is **+11.6% to +13.6%** of today's cosine; the
   shape-only bar is **+2.5% to +2.9%**. **A fifth independent donor line re-prices from 205's
   +0.0062 to +0.0140…+0.0163 — from 1.4× below the bar to 1.6–1.8× above it.** → 218.
6. **Nothing in the eleven direction closures reverses.** Three closures re-open, and all three
   turn on **the renderer** (214's 1.84×) or on **rule 13's floor being the RMS of channels the
   arm has since priced** — not on §8.4's DE propagation. The strongest-looking direction
   candidate this audit found (158) was killed by a measurement already on disk (165 §7).

---

## 0. Phase 0 — the correction, re-derived and given a range

`phase0_219.py` → `work/phase0_219.json`. **40 gates, 0 FAIL.** Every constant below is
re-derived from the artifacts; none is inherited from §8.4's prose.

### 0.1 The known-answer rows (rule 10, and 181's "reproduce one number from each audit")

| block | row | expected | measured |
|---|---|---:|---:|
| **181** | six scaled members average to `score_avg` on every board entry, `2026-09-02T0630Z.json` | 0 failures / 581 | **0 / 581** |
| **198** | 150's six-member average reconstructed from its own table | −0.005045 | **−0.005045** |
| **198** | 151's six-member average | −0.004344 | **−0.004344** |
| **198** | 150 ex-`reach` | +0.002476 | **+0.002476** |
| **PLOLO** | a15 `pds`, `k562_gwps` fold, equal weights (204 row 2) — **fold rebuilt from scratch** | 0.79578409 | **0.7957840906989219** |
| **PLOLO** | a23 `pds`, same fold (204 row 2b) — **rebuilt** | 0.78912131 | **0.7891213079583496** |
| **PLOLO** | a15 `q` / `q_w`, same fold (204 row 2b) — **rebuilt** | 0.50164311 / 0.50487994 | **0.5016431147 / 0.5048799396** |
| **PLOLO** | a15 → a23 mass-fix price, 4 folds, scaled (162 / 204 row 3) | +0.004014 | **+0.0040145** |
| **197** | the inversion's two constants `1.023611` / `2.074495` from `Kc`, `ΣD`, `K` | — | **exact** |
| **197** | the floor reproduces on all six of its own calibration rows | — | **6/6** |
| **207** | `cos_agg` at `f` = 0.10 / 0.34 / a14 | 0.020086 / 0.023698 / 0.019053 | **identical** |
| **212** | macro cosine `f` 0.10 → 0.34 = +13.83% | 1.1383 | **1.13828** |
| **8.2/212** | `d(overall mse)/d(cos)` = +0.3655, on **both** board points | 0.3655 | **0.3639 / 0.3633** |
| **210** | the six-member reconstruction of `cd4f034`'s +0.0130557 | +0.0130557 | **+0.0130559** |

The `mse` exchange-rate row is the strongest of these: **+0.3655 per unit of board cosine,
derived by 212 from the gamma axis, reproduces to 0.6% on two independent CD4-dose board points
that were not used to fit it.** Everything below rides on that rate, so it is worth saying that
it is now a three-way agreement, not a construction.

> **Provenance note.** The `k562_gwps` leg is **rebuilt from scratch** by `ka_plolo219.py`
> (7.6 s, peak RSS 4.85 GB): the donor sources are re-streamed and re-mixed, the CD4 banks
> re-prepared, the held-out operator `13d988fc…` re-applied, and the result re-rendered and
> re-scored. `max_abs_dev_from_162_pre = 4.13e−08` confirms the rebuilt baseline is 162's. Only
> the **other three folds** in the 4-fold mass-fix mean are read from 212's cache, whose gate 10
> reproduced 204's eleven published scalars at `rel 0.00e+00`. The fold caches were **copied**
> into `work/plolo/` before `plolo.WORK` was redirected, exactly as 212 did, so a cache miss
> could not write into a sibling experiment.

### 0.2 The 2.51×, decomposed — and where §8.4 double-counts

`cd4f034` realised **+0.0130559** against 207/210's preregistered `mse`-only price of
**+0.005206**. That ratio is **2.5078**, which is §8.4's 2.51 — and it factors **exactly**:

```
2.5078  =  P x E  =  1.9974  x  1.2556
                     ^^^^^^     ^^^^^^
                     realised   realised mse
                     total      -----------
                     ---------  predicted mse
                     realised
                     mse
```

* **`P` — the propagation.** *"A direction gain lands on five members, not one."* This is
  §8.4's Error 1 stated honestly, and it is **1.997**, not 2.51.
* **`E` — the transfer residual.** *"207's price under-read the `mse` leg."* This is Error 2
  showing up in the `mse` leg, and it is **1.256**.

**§8.4 then multiplies 2.51 by Error 2's 1.80 to get 4.5.** But 2.51 already contains `E`, and
`E` *is* the transfer error that survived 207's own `k_cos = 1.2598`. The implied total transfer
factor is `1.2598 × 1.256 = 1.582`, against 207's `cos_agg` reading of **1.384** — 14% apart,
one number, not two. **4.5× charges the same phenomenon twice.**

### 0.3 Error 2's 1.80× is measured on the wrong cosine

The board's cosine gain is read off 197's inversion of the **aggregate** cosine
`X / (√ΣD · mean_c √A_c)`. §8.4 divides it by **212's macro-per-fold cosine** ratio. Those are
different functionals of the same folds, and 207's `cos_agg` is the matched one:

| route | PLOLO Δcos | board Δcos | `k_cos` |
|---|---:|---:|---:|
| `a14 → a15`, `cos_agg`, 197's *measured* cosines | +5.42% | +6.83% | **1.261** (= 207's published 1.2598) |
| `a14 → a15`, `cos_agg`, floor-ratio | +5.42% | +6.37% | **1.175** |
| `rc_callset → cd4f034`, **`cos_agg`** (matched) | +17.98% | +24.89% | **1.384** |
| `rc_callset → cd4f034`, 212's macro ← **§8.4's 1.80** | +13.83% | +24.89% | 1.800 |

**`k_cos` = 1.17–1.39, two board points, not 1.80 on one.**

The floor-ratio and measured-cosine estimators of the board's Δcos agree to **0.47 pp** on the
`a14 → a15` pair, which is the check that licenses using the floor ratio for `cd4f034` (where
`claim` was never measured on the artifact). 197's own tightness table puts `floor/measured` at
**0.911–0.965** on the four near-cliff rows, so the ratio-of-floors is a ratio of two numbers
each 91–97% of the truth — a ±5% differential at worst.

### 0.4 The multiplier does not rest on one point, and it is dose-dependent

| board pair | Δcos (floor) | `mse` | `pds` | DE-four | total | **P = total/mse** |
|---|---:|---:|---:|---:|---:|---:|
| `a14 → a15` (`f` 0.053 → 0.10) | +0.00432 | +0.001574 | **+0.001209** | +0.000879 | +0.003662 | **2.33** |
| `rc_callset → cd4f034` (`f` 0.10 → 0.34) | +0.01799 | +0.006537 | **−0.001345** | +0.007864 | +0.013056 | **2.00** |

Per unit of board cosine:

| | `mse` | `pds` | DE-four | **total** |
|---|---:|---:|---:|---:|
| `a14 → a15` | +0.3639 | **+0.2797** | +0.2033 | **+0.8470** |
| `rc_callset → cd4f034` | +0.3633 | **−0.0748** | +0.4372 | **+0.7257** |

**Three things this says that §8.4 does not.**

1. **The `mse` rate is a constant** (0.3633 / 0.3639, against 212's 0.3655 from a third axis).
2. **The `pds` leg changes sign with dose** — 207 §4.1 measured the `pds` maximiser at
   `f = 0.20` and the cosine maximiser at `f = 0.80`, and the board straddles it. So `P` is not
   a constant of the *class*; it is a function of where on the dose curve the arm sits, and it
   falls as the dose rises. **At 217's `f = 0.80` the PLOLO `pds` cost is −0.0114 raw against
   +0.007 at `f = 0.10`, so `P` there is materially below 2.0.** That is 216/217's problem, and
   it is the one thing in this audit they should read.
3. **The DE-four rate more than doubled between the two points (0.203 → 0.437).** The likely
   mechanism is not the dose: `a14 → a15` was measured on the `quasi` emitter at `n_pred` 4,424,
   `cd4f034` on `rc_callset`'s prescribed 2,599-gene call set that took `jac` to rank 1/646. **A
   direction gain converts to DE members ~2× more efficiently now that the ordering channel is
   fixed.** Any re-price of a pre-`rc_callset` arm should use the higher rate.

### 0.5 The multiplier's range, and what it is a multiplier OF

**`M` = realised overall ÷ (096 `mse`-only price), for a cosine arm that does not move the
operator's exponent, priced 207's way (PLOLO cosine → `k_cos` 1.26 → 096 `claim` route):**

| | value |
|---|---:|
| point (the one arm where a preregistered price and a board row both exist) | **2.51** |
| range | **2.00 – 2.56** |
| ← if the old price used raw PLOLO cosine with no `k_cos` | × 1.26 → **2.5 – 3.2** |
| **operator-geometry arm** (the exponent, the eigenvalue map, a learned metric) | **0.59 – 1.23** |
| amplitude-only arm (`α`, dispersion) | **1.00** — neither `pds` nor `cos` sees a rescale |

**The uncertainty is honestly one-sided in one place and two-sided in another.** `P` has two
board points and they bracket it. `E` has **one**, and it is the residual of a transfer factor
that was itself fitted on one point — so I carry `E ∈ [1.00, 1.39]` rather than a point, and the
upper end of every re-price below is flagged. An audit that turned one board row into a fleet of
confident reversals would be worse than no audit; this one turns it into **three re-openings,
none of which are direction arms.**

### 0.6 The consequence that matters more than any single re-price

§8.2 states *"a single direction arm must deliver +31% cosine to clear rule 13's gate alone."*
That is `0.0089 / 0.00365 per 0.01 cosine`, i.e. **the `mse`-only rate**. On the measured total
rate:

| gate | required Δcos, §8.2's construction | **required Δcos, corrected** |
|---|---:|---:|
| shape-only **+0.0019** | +0.00520 abs = **5.8%** of today's cosine | **+0.00224…+0.00262 = 2.5%…2.9%** |
| amplitude **+0.0089** | +0.02435 abs = **27.0%** (published as 31% at the pre-`cd4f034` cosine) | **+0.01051…+0.01226 = 11.6%…13.6%** |

**Every direction NO-GO in this project was measured against a bar that is 2.0–2.3× too high.**
The re-prices in §2 are what happens when the bar is corrected: mostly nothing, because most of
the arms are 5–20× under it, not 2×.

---

## 1. Phase 1 — the full classification

`enumerate219.py` → `work/enumeration219.json` (196 rows, machine-readable);
`classify219.py` → `work/classification219.json` (151 rows).

| | count |
|---|---:|
| `experiments/*/RESULTS.md` on disk | **196** |
| terminal verdict (NO-GO / VETO / CLOSED / NEGATIVE) | **151** |
| verdict tokens | NO-GO 79 · CLOSED 38 · NEGATIVE 26 · VETO 8 · GO/SHIP 7 · none 38 |
| **terminal verdicts that quote a rule-13 floor (+0.0019 / +0.0089)** | **44** |

### 1.1 The mechanical partition that does most of the work

**§8.4's Error 1 is an error in a pricing *convention* — "price the arm on `mse` through 096's
`claim` route and carry the DE members at zero". That convention did not exist before experiment
162** (rule 13 was fitted 2026-09-02 against 12 board points; AGENTS.md §11.1 rule 13). An
experiment that never computed a rule-13 composite has **no `mse`-only price for §8.4 to
correct** — it was closed on `pds`, on cosine, or on a structural fact.

So the correction can only reach a closure that quotes a rule-13 floor **and** is a
direction/profile arm. Applying that test:

| class | n | what it is |
|---|---:|---|
| **PRE** | **86** | closed before the convention existed; hand-checked on the seven files the keyword screen flagged (010, 012, 014, 063, 098, 116, 121) — **zero real rule-13/096 hits in any of them** |
| **B — unaffected** | **39** | priced on all six already, or negative by more than the multiplier, or closed by a ratio/structural argument |
| **ORD — ordering, → 220** | **13** | the order half of `reach`/`fid`/`jac`; a different channel with a different instrument (§8.1 rule 1) |
| **C — closed structurally** | **7** | closed by a fact, not a number |
| **DEAD — falsifier 3** | **3** | a permutation control tied or beat the real arm |
| **A — re-priceable** | **3** | 200, 208, 214 |

### 1.2 Class A — re-priceable (3)

| exp | verdict | what makes it re-priceable |
|---|---|---|
| **200-production-ablation** | NEGATIVE | 3 of its 10 rows are priced on `pds` + `mse` **only** (`n = 2` in its own table): `VCC_ALPHA`, the **post-whitening CD4 residual**, and the pruned basis. And 214 §4 showed **every** row was measured on a renderer (`VCC_MASS_FIX=minnorm`) that has not shipped since a23 |
| **208-free-riders** | VETO | rider A (Hamilton rounding) priced **+0.001197** against a **+0.0019** floor. 214 re-measured its PLOLO `pds` on the shipped renderer at **1.84×** |
| **214-cheap-positives** | NEGATIVE | Arm A prices **all six** at **+0.0054…+0.0088** and is then vetoed by a floor that is *defined* as the realised RMS of the channels the arm has just priced |

**None of the three is a direction arm.** All three turn on the renderer or on rule 13's floor
construction. That is the audit's central negative result.

### 1.3 Class B — unaffected (39), with the reason for the ones the mandate named

| exp | family | why the multiplier does not move it |
|---|---|---|
| **193** whitener shrinkage | operator | Priced on **all six** (§5's composed table). Falsifier 1 fires. And it is an operator-geometry arm, so its multiplier is **0.59–1.23**: +0.000512 × 1.23 = **+0.00063** against +0.0019. Best of 55 maps. |
| **212 / 170** `WHITEN_GAMMA` | operator | The **ratio** argument: the operator trades cosine for retrieval at **96.2%** of the board's own rate. Instrument-independent, as the mandate says. My measured operator multiplier (0.59–1.23) is the same statement arriving from the board side. |
| **195** learned retrieval transform | operator | −0.00554 `pds`, 95% CI excluding zero, plus a measured mechanism (dies at the cell-line boundary, −0.00618 across an unseen line). Large negative under any multiplier. |
| **213** transport law | direction | **A transport law is a per-gene, sign-PRESERVING gain — 0.0000% sign flips over 7.0 M donor entries.** §8.4's propagation runs entirely through sign accuracy, so the DE-four leg is structurally ~0 here. Best arm +0.00103 against +0.0089 with `α` falling 0.982 → 0.7438. |
| **203 / 125** EB panel posterior | direction | Cosine **+1.6%**, CI including zero, against a corrected bar of +2.5–2.9%. And its own §6 measures the marginal refund rate at **`e` = −0.193**: the arm's extra L1 lands on **wrong-signed** genes. The DE-four propagation cannot be positive for an arm whose own measurement says its increment is anti-signed. |
| **158** per-target confidence scale | direction | The strongest-looking candidate in the sweep, and it dies on evidence already on disk — see §2.3. |
| **194** donor re-estimation | direction | The recoverable part is a **magnitude** correction and `pds_cosine` is exactly scale-invariant per row. +0.0003 against +0.0019; at 2.51× it is +0.00075. |
| **173 / 178** sign repair and the riders | direction | Priced on all six. 178 re-read `pds` on the board-validated PLOLO frame at **−0.0245 raw**, putting sign repair at **−0.0095** overall. Large negative. |
| **169 / 168 / 188** magnitude sharpening | profile | Priced on all six; 188 closed the `nmae` surface question **against** the arm on a third instrument. |
| **211 / 187 / 180** amplitude | amplitude | Amplitude moves neither `pds` nor the cosine (both exactly scale-invariant), so **no direction correction applies at all**. 211's +0.00245 against +0.0089 is a *floor* question, not a pricing one. |
| **097** the original direction ranking | direction | Its own recommended arm **shipped as a13 and the board scored it −0.0203 vs a12.** A board point outranks any re-price. |
| **183** unique-component boost | direction | It **lowers** the transfer cosine on 3 of 4 surfaces with truth. The multiplier applies to cosine *gains*. |
| **165** head-protected mass fix | profile | The mandated arm is measurably worse than a23 and its premise is refuted deterministically. Its `eps_tail` sibling (+0.0013 vs +0.0019) is a floor question and 165 records its value swinging **7.8×** with an unpinned calibration. |
| **150 / 151** dual moment | profile | **Already re-priced by audit 198** (R1/R2, corrected to `[−0.00042, +0.00077]` and `[+0.00013, +0.00127]`). Not re-done here — 198's method is my method and its conclusion stands. |

### 1.4 Class C — closed structurally (7)

| exp | the fact |
|---|---|
| **205** X-Cell / Pisces | the artifacts do not exist (3 files, `used_storage: 0`, unchanged since 2026-03-25) |
| **166** donor acquisition | no unexploited panel-covering CRISPRi line exists in public data (HF, 258 GEO series, 229 papers, IGVF's 6,902 measurement sets) |
| **184** Tahoe-100M | the winning law is sign-preserving and the operator is absent in CRISPRi |
| **176** LINCS/CMap | coverage cap; 181 measured the cap **generous by 40%** (landmarks are 9.5% of the reference DE pool, not 13%) |
| **182** pretrained State | structurally undeployable |
| **167** target expression | range restriction by panel construction (AGENTS.md §2.1) |
| **099** more donors | 3 of 7 prepared donors carry **0/300** panel targets by construction |

**205 is Class C on feasibility and re-prices 2.3× on value — see §2.1.**

### 1.5 DEAD — falsifier 3 (3)

Closed by a permutation control that tied or beat the real arm. **Not reopened, regardless of
scale**, exactly as the mandate requires.

| exp | the control |
|---|---|
| **209** Orion subspace deflation | the maximiser over all 36 arms **on the cosine channel** IS a permutation control |
| **190** regulon ordering | a shuffled regulon selected the same way reads **0.23699** against the real **0.23670** |
| **204** repression-weighted mixing | equal weighting is the measured **maximum** of the transfer cosine; the model's own predicted +8.5…+16.8% gain is excluded by the measurement |

### 1.6 ORD — ordering, handed to 220 (13)

153, 154, 157, 164, 171, 172, 177, 179, 185, 189, 191, 192, 202. **I do not re-price these** —
they are 220's territory and §8.1 rule 1 says never price an ordering arm on `pds`. Two facts to
carry over:

* **189 and 190's instrument tightened 77× under `rc_callset`.** 214 §5 measured the three-seed
  paired `reach` sem at **0.00009** raw against 200's **0.00695** on the `quasi` emitter. 189's
  headline was `reach` **+0.00408 ± 0.00424** — a 0.96σ result whose error bar was set by an
  emitter we no longer ship. On the current instrument that same effect would be **~45σ**.
  **That is a stronger reason to re-open 189 than anything in §8.4**, and it is 220's.
* **190 stays DEAD anyway** — its shuffled control beat it, and a 77× tighter instrument
  resolves a *negative* just as well as a positive.

---

## 2. Phase 2 — the Class-A re-prices, with the DE members ON

`reprice219.py` → `work/reprice219.json`.

### 2.1 The largest number in this audit is not a Class-A arm — it is 205's, and it is 218's

205 priced one extra effective donor line at **+21.3% cosine → +0.0062 overall** using §8.2's
`mse`-only exchange rate, and concluded it was *"a fifth of the gap, not the gap"* and **below**
the +31% single-arm bar.

| | 205's price | **corrected** |
|---|---:|---:|
| one extra effective line (+21.3% cosine) | +0.0062 | **+0.0140 … +0.0163** |
| two fully independent lines (+31.0% cosine) | +0.0102 | **+0.0203 … +0.0237** |

**A fifth panel-covering line goes from 1.4× below the amplitude floor to 1.6–1.8× above it.**
This does not reverse 205 or 166 — the artifacts still do not exist, and 166's search was
exhaustive. It reverses the **value** of continuing to look, by 2.3×. **→ 218**, which owns
donor independence and acquisition.

### 2.2 208's Hamilton rider, on the renderer that actually ships

208 priced rider A's rule-13 composite at **+0.001197** against **+0.0019** → VETO. 214 §4 then
found that 208 (and 200) measured on `VCC_MASS_FIX=minnorm`, which **only a23 ever shipped** —
and a23 scored **−0.00073** and was abandoned. On the shipped renderer:

| | 208 | 214, shipped renderer | factor |
|---|---:|---:|---:|
| Hamilton PLOLO Δ`pds`, raw | 0.001124 | **0.002069** | **1.84×** |

| | value | vs floor |
|---|---:|---|
| rule-13 priced composite, corrected | **+0.002203** | **CLEARS +0.0019** (1.16×) |
| six-member, corrected (`nmae` leg unchanged) | **+0.003909** | |
| …**with** the emission-law toll −0.00243 charged | −0.000227 | fails |

**The verdict now turns entirely on one unresolved question: does rule 13's −0.0146 `reach` toll
apply to a rounding rule whose profile change is a 0.03% tilt correction?** 208 proved the arm is
not pseudobulk-preserving (900/900 blocks differ), so the letter of the rule-13 table charges it.
But the toll's two calibration points (a21, a22) moved the profile *enormously* — a21 took
`n_pred` 6,034 → 1,913 and `nmae` −0.064 scaled — while Hamilton moves the challenge tilt from
**−0.01253 to +0.00030 log2**, i.e. it *removes* a systematic rather than adding one. Charging a
class constant fitted on two catastrophic profile moves to a 0.03% one is the exact error §8.4's
closing paragraph warns about in the other direction.

**This is a re-price, i.e. a hypothesis.** It is not a slot.

### 2.3 The candidate that looked like the audit's headline, and why it is not — 158

I want this on the record because it is the closest this audit came to a reversal and because
the way it died is the method working.

158 closed a per-target confidence rescale that delivers **+16.2% / +29.6% aggregate transfer
cosine at exactly zero `pds` cost** (verified numerically: `Δpds = +0.00e+00` to the last printed
digit for any strictly positive scale, both surfaces, every arm), with `reach` **exactly
invariant** on 132/132 targets and — under `rc_callset`'s *prescribed* call set (208 §2) —
`fid`/`jac` structurally pinned. On the corrected bar (+2.5–2.9% shape-only) a +16% cosine gain
is 6× **over** it.

**And 158's `mse` price is wrong by an order of magnitude, in the arm's favour.** Its §4 uses
`d_raw_mse = −κ·c²·((1+f)² − 1)` — the **uncorrected** `1 − cos²` form. AGENTS.md §4.1b names
exactly this: *"the bound is real and it only applies to the uncorrected squared error."* The
scored metric is `claim`-refunded and **linear** in cosine. Tested on the one board point that
can test it: 158's formula predicts `cd4f034`'s `mse` gain at **+0.000511** against a realised
**+0.006537** — **understated 12.8×**.

So the arm reads: `pds` 0 (measured), `reach` 0 (measured), `fid`/`jac` pinned, `mse` **+0.0053
… +0.0097** (corrected), `nmae` **−0.0029** (measured, re-levelled at its own optimum) → net
**+0.0024 … +0.0068** against a +0.0019 shape-only floor.

**Then falsifier 1 fires, on evidence already on disk.** 158's `f` is measured on the two Orion
surfaces, and **165 §7 retrodicted exactly this instrument against the board's own a23 point**:

> *"156 measured aggregate cosine gains of +13.6% (HCT116) and +56.0% (HEK293T); `h1-test` reads
> +0.89%; the board realised a +1.16% relative gain. The Orion pair **overstates** that by 12×
> and 48×."*

Deflated at 165's kinder factor (12×), 158's `f` is **+1.4%**, its `mse` leg is +0.00045, and the
arm is **−0.00245** — negative. Its `nmae` cost does not deflate; it was measured directly.

**158 stays closed.** The reason is not §8.4's convention and never was; it is that 158's cosine
instrument has a board point in its own family and gets it wrong by 12–48×. 165's own standing
recommendation — *"price any future profile-geometry arm on `h1-test`, and use the Orion pair
only for sign"* — is the operative rule and it was already written down.

### 2.4 200's ablation: the three rows priced on two members, and the renderer

200's LOO table prices six of ten components on all six members and **three on `pds` + `mse`
only**:

| component | 200's sum | members priced |
|---|---:|---|
| `VCC_ALPHA` 0.982 | +0.018542 | 2 — and amplitude cannot move the cosine, so no correction applies |
| **post-whitening CD4 residual (`f` 0.10)** | **+0.002635** | **2 — this is exactly §8.4's error** |
| pruned basis (`union_pruned`) | +0.000373 | 2 |

**The board has since measured half of that CD4 row directly.** The `a14 → a15` step (`f` 0.053
→ 0.10) realised **+0.003662 overall** — more than 200's price for the *entire* residual
(`f` 0 → 0.10). **200's CD4 row understates by ≥ 1.39×, and the missing part is precisely the
four DE members it carried at zero.** It is the cleanest confirmation of §8.4's Error 1 available,
and it changes 200's own verdict for that row from *"noise (gate ±0.0089)"* to a load-bearing
component.

Separately, **every** row of the ablation inherits 214's renderer finding: Hamilton reads 1.84×
larger and PerturbMatch 1.77× smaller on the shipped renderer. PerturbMatch's row moves from
−0.001615 to roughly **−0.00244** (still noise, still negative); the six-member rows' `pds` legs
all move by an unmeasured amount in an unmeasured direction. **200's "six load-bearing, four
inside their gate" ordering is not safe to quote until it is re-run with `VCC_MASS_FIX` unset.**

### 2.5 214's Arm A, and the floor that eats its own tail

214's Arm A (drop PerturbMatch) prices **all six** at **+0.0088** native / **+0.0054** after its
own measured ×1.533 h1-test deflator / **+0.0046** if `reach` transfers at 164's 0.4277 — and is
then vetoed by comparing rule 13's **two-member** composite (−0.000830) against the **+0.0089**
floor.

**Rule 13's floors are, in AGENTS.md's own words, "the realised RMS of the *unpriced* channels
for each arm class."** An arm that has priced all six has no unpriced channels. Comparing a
six-member estimate against the RMS of the channels it just priced double-charges the same
uncertainty — and comparing a *two*-member composite against it while a six-member estimate
exists throws away the measurement. **This is a defect in the gate, not in 214**, which reported
both numbers and chose the conservative one.

The honest statement for Arm A: **+0.0046 … +0.0088, with every point of it in `reach` and
`nmae`, and the one channel rule 13 does price (`pds`) negative at −0.001059.** That is a real
arm with a bad risk profile, not a vetoed one. 214 owns it.

---

## 3. Phase 3 — ranked shortlist

Ranked by **(corrected estimate × confidence) ÷ cost**. **Every row is a hypothesis. Nothing
here goes near a slot without a fresh measurement and a preregistration.**

| # | arm | corrected estimate | conf. | what must be re-measured | cost |
|---:|---|---:|:---:|---|---|
| **1** | **Re-run 200's ten-component ablation with `VCC_MASS_FIX` unset and all six members on every row** | not an arm — it re-prices ten | high | 200's own LOO harness at the shipped renderer; the CD4 row's DE-four legs; the PerturbMatch row at 214's 1.77× | **LOW** — the harness exists, no build |
| **2** | **208 rider A (Hamilton), shipped renderer** | **+0.0022** priced / **+0.0039** six-member | med | (a) 214's +0.002069 `pds` reproduced against `cd4f034` not `rc_callset`; (b) whether the −0.0146 emission-law toll is chargeable to a 0.03% tilt correction — the honest test is a per-gene panel-total digest *plus* a magnitude, not the digest alone | **LOW** |
| **3** | **The corrected +11.6–13.6% direction bar, written into AGENTS.md §8.2** | re-prices every future direction arm | high | nothing — it is Phase 0 | **ZERO** |
| **4** | **Donor acquisition, re-valued for 218** | **+0.0140…+0.0163** per independent line | med | nothing new; 166's search verdict is unchanged. The re-price only changes how much a *new* source is worth if 218 finds one | **ZERO** (→ 218) |
| **5** | **214 Arm A (drop PerturbMatch), against a gate that does not double-charge** | **+0.0046…+0.0088** | low | its `reach` transfer factor, which is the whole arm; and a `pds` re-read on `cd4f034` | **MED** (→ 214) |
| **6** | **189's per-gene reliability prior on the 77×-tighter instrument** | reach +0.00408 → resolvable | med | 220's job, not mine | **LOW** (→ 220) |
| — | 158, 193, 203, 212, 213, 194, 173/178, 097, 211 | all **negative or ≤ ⅓ of their floor** after correction | — | nothing | — |

**Rows 3 and 4 cost nothing and are worth more than rows 1, 2 and 5 combined.** That is the
audit's recommendation: the correction's value is that it re-prices the *bar*, not that it
resurrects arms.

---

## 4. Falsifiers, as required

**F1 — the correction must not reverse a verdict whose own evidence contradicts it.**
Checked on every candidate. It fires for:

| exp | measured the DE members directly? | consequence |
|---|---|---|
| **193** | yes — §5's composed six-member table | 2.51× does not apply; the operator multiplier 0.59–1.23 does |
| **169 / 168** | yes — all six, and 188 re-measured `nmae` on a third instrument | does not apply |
| **211** | yes — all six | does not apply, and amplitude cannot move the cosine anyway |
| **214 Arm A** | yes — all six | does not apply; the issue is the floor |
| **158** | yes — `pds` exactly 0, `reach` exactly invariant, `nmae` measured and surviving the oracle | does not apply; **and its own `nmae` measurement is the binding leg** |
| **203** | yes — `nmae` measured, and `e` = −0.193 says the increment is anti-signed | does not apply, and the propagation's *sign* is contradicted |
| **213** | structurally — 0.0000% sign flips over 7.0 M entries | the propagation channel does not exist for this arm |
| **150 / 151** | yes — all six (198 R1/R2) | 198 already did this work |

**F2 — a reversal that depends on the multiplier's upper end is not a reversal.** Flagged
separately: **none of the three Class-A re-opens depends on the multiplier at all.** 208's turns
on 214's measured 1.84× renderer factor; 214's and 200's turn on the floor construction and the
renderer. The only numbers in this audit that *do* depend on the multiplier are §2.1's donor
re-price and §0.6's bar — and both are quoted as ranges throughout (+0.0140…+0.0163;
+11.6%…+13.6%).

**F3 — nothing closed by a permutation control that tied or beat the real arm is re-opened.**
204, 209, 190 are classified **DEAD** and stay there. I checked their controls directly rather
than taking the label: 209's maximiser over 36 arms on the cosine channel *is* the permutation;
190's shuffled regulon reads 0.23699 against the real 0.23670; 204's equal weighting is the
measured maximum of the transfer cosine.

---

## 5. What this corrects in AGENTS.md

1. **§8.4's "Combined ≈ 4.5×" should read ≈ 2.5× (range 2.0–2.6).** The 4.5 multiplies 2.51 —
   which is `realised total ÷ predicted mse` and therefore already contains the transfer error —
   by 1.80, which *is* the transfer error. Double count. The two lines of §8.4's pricing rule are
   **alternative routes, not composable**.
2. **§8.4's "PLOLO understates the board's cosine response by 1.80×" should read 1.17–1.39×**,
   and the number should be quoted on 207's `cos_agg`, which is the functional 197 inverts.
   212's macro-per-fold cosine is a different statistic.
3. **§8.4's "validated on ONE point" is wrong — there are two.** `a14 → a15` is a board-measured
   CD4 dose step on the same knob and gives the propagation factor at **2.33×** against
   `cd4f034`'s **1.997×**.
4. **§8.4's rule should carry the class boundary.** *"Price all six members explicitly"* is right
   for a cosine arm that leaves the operator alone. For an **operator-geometry** arm the
   measured multiplier is **0.59–1.23** — the `pds` leg cancels 96% of the `mse` leg (212 §3.4)
   — and for an **amplitude** arm it is exactly 1, because both `pds` and the cosine are
   scale-invariant.
5. **§8.2's "+31% cosine" single-arm bar is an `mse`-only construction.** Corrected: **+11.6% to
   +13.6%** of today's cosine against the amplitude floor, **+2.5% to +2.9%** against the
   shape-only one. This is the single most reusable output of this audit.
6. **Rule 13's floors are the RMS of the *unpriced* channels, and the rule does not say what to
   do when an arm prices them.** Comparing a six-member estimate against them double-charges.
   214 and 211 both hit this; 214 said so.
7. **§8.4's `pds` cost rule ("2.0× the PLOLO prediction") is one-signed and one-point.** On the
   other board point PLOLO **understated a favourable** `pds` move by 1.19× (207's KA-1:
   +0.006109 predicted, +0.007257 realised). The honest form is *"PLOLO under-reads |Δpds| by
   1.2–2.0×; it does not systematically over-read the cost."*

---

## 6. Corrections to sibling experiments (flagged, not asserted)

* **158 §4's `mse` pricing route** uses the uncorrected `1 − cos²` form and understates the
  board's realised `mse` gain by **12.8×** on the one board point that tests it. Its verdict is
  unaffected (§2.3), but the formula should not be reused. AGENTS.md §4.1b already names the
  error class.
* **205's +0.0062** and every other price computed through §8.2's `+0.00365 per +0.01 cosine`
  is an `mse`-only price and should be multiplied by 2.0–2.3.
* **210's BOARD-RESULT table transcribes `rc_callset`'s raw `mse` as 0.87367**; the scaling map
  solved from a14/a15 (and the ledger's own statement that `mse` is identical to sixteen
  significant figures across a15 / `rc_hybrid` / `rc_callset`) puts it at **0.873805**. The
  cosine ratio moves 1.2489 → 1.2487, i.e. nothing depends on it — recorded only because a
  transcribed number that is quietly 0.0001 off later gets treated as a real discrepancy (198's
  own note on the ledger's 0.12423).

---

## 7. Integrity, and what I did NOT measure

No `vcc`, no `cell-eval2`, no `cell-eval`, no build, no writes to `submissions/`, no production
code touched, no file outside `experiments/219-reprice-closures/` created or modified. Peak RSS
**4.85 GB** (`/usr/bin/time -l`, `ka_plolo219.py`, the one-off PLOLO fold rebuild), against a
10 GB cap; every other script peaked at **0.07 GB**.

| not measured | why |
|---|---|
| a from-scratch rebuild of the other three PLOLO folds | The `k562_gwps` leg **was** rebuilt (§0.1). The remaining three enter only the 4-fold mass-fix mean and are read from 212's gate-10-verified cache. Rebuilding all four would have cost ~4× the memory for a known-answer row that already passes on the fold 204 published. (An earlier attempt was killed without a traceback while three sibling agents held 4.2 / 2.8 / 2.2 GB; the retry succeeded in 7.6 s once the machine had headroom, which is worth recording — a silent kill on this box is memory pressure, not a bug.) |
| any re-measurement of a Class-A arm | out of mandate — "your job is to say what is worth re-measuring" |
| 150/151 | 198 already re-priced them; repeating it would be duplicating a prior audit |
| the ordering keys (189, 190, 172, 185…) | 220's territory; handed over in §1.6 |
| the CD4 dose ladder beyond `f` = 0.34 | 216/217's territory. §0.4 point 2 is the one thing they should read: `P` falls as the dose rises because the `pds` leg turns over at `f` ≈ 0.20 |
| donor independence / acquisition | 218's; §2.1 is handed to them |
| `cd4f034`'s `claim` on its own artifact | would sharpen `E` from one point to a measurement, and it needs a build stage this mandate forbids. It is why `E` is carried as [1.00, 1.39] |

---

## 8. Files

| file | what it does |
|---|---|
| `phase0_219.py` | the 40 known-answer gates, 197's inversion re-derived, both PLOLO cosine functionals, the propagation on two board points, the 2.51/1.80/4.5 decomposition, the operator-geometry counter-rate → `work/phase0_219.json` |
| `enumerate219.py` | sweeps all 196 `RESULTS.md`, extracts verdict / priced members / zeroed members / instruments / rule-13 citations → `work/enumeration219.json`, `work/candidates219.json` |
| `classify219.py` | the A/B/C/PRE/ORD/DEAD classification of all 151 terminal verdicts, with the reason recorded per row → `work/classification219.json` |
| `reprice219.py` | the corrected §8.2 bar, the Class-A re-prices, 158's adjudication, 200's two-member rows → `work/reprice219.json` |
| `ka_plolo219.py` | the from-scratch PLOLO `k562_gwps` fold rebuild — 204's a15/a23 `pds` and `q`/`q_w` reproduced exactly, `max_abs_dev_from_162_pre` 4.13e−08 → `work/ka_plolo219.{json,log}` |

**Inputs read, never written:** `AGENTS.md`; `docs/BOARD-LEDGER.md`; `data/leaderboard/2026-09-02T0630Z.json`; every `experiments/*/RESULTS.md`; `experiments/212-gamma-cosine/work/{plolo212,cd4f212_summary}.json`; `experiments/207-cd4-dose/work/sweep207.json`; `experiments/210-cd4-build/BOARD-RESULT.md`. The fold caches under `experiments/{204,212}/work/plolo/` were **copied** into `work/plolo/`, never opened for writing; the 131 MB of clones was **deleted after the run** so it does not enter the repo (`experiments/**/work/` is not gitignored). Re-running `ka_plolo219.py` re-clones them. No `__pycache__` was written into 162's or 204's directories — the loader uses `spec_from_file_location`, verified by `find -newermt`.
