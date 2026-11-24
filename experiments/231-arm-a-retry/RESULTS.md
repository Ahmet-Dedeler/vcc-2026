# 231 — 214's Arm A is **REFUTED**, twice, on measurements the 2026-09-04 board makes available

Date **2026-09-04**. Measurement only. No `vcc`, no `cell-eval2`, no build, no `vcc prep`, no
submission, no `submissions/` write, no production change, and no file outside
`experiments/231-arm-a-retry/` created or modified. Peak RSS **5.23 GiB** (`class231.py`; the
other two heavy scripts are 5.05 and 5.18 GiB). **No slot is proposed. Nothing was built.**

---

## Verdict in ten lines

1. **The headline is refuted on arithmetic.** Arm A's `+0.0088` requires a board raw `reach` of
   **0.31223** against a world best of **0.290835** — **250% of the entire remaining gap to the
   best `reach` team on earth**, from removing one retrieval component. Even 214's own
   ×1.533-deflated `+0.0054` implies **0.29983**, still above the field best.
2. **220's bound has shrunk by 65% since it was written, and that is what kills the arm.** 220
   priced the whole ordering channel off our **0.24988** `rc_callset` row. The shipped `h1union`
   row reads **0.276545** — reach rank **7 of 656**. The remaining prize for becoming the best
   `reach` team in the world is now **+0.002689 overall**, not +0.0077.
3. **So the arm class is closed by arithmetic, permanently.** Arm A is **amplitude-changing**
   (delivered median ‖row‖ ratio **1.02819**, re-measured on the shipped stack), so its floor is
   **+0.0089**. The entire available `reach` prize is **+0.002689**. **The floor is 3.31× the
   whole prize.** No amplitude-changing arm whose gain lives in `reach` can clear its own gate at
   any transfer factor, *even at oracle*.
4. **The transfer factor — "which is the whole arm" — is NOT IDENTIFIED for Arm A's class.** Four
   consecutive single-variable board pairs give **0.390, 1.2927, 1.4346, 3.0563** (7.84× spread,
   geometric sd **×/2.343**). The only pair in Arm A's *own* class, `cd4f034 − rc_callset`, is
   **0/0**: h1-test moved by **exactly zero** (bit-identical effect matrices) while the board moved
   **+0.023416 raw**. **n = 0 usable pairs.** 214's headline used 1.290, borrowed from a different
   arm class.
5. **The `pds` re-read shrinks the penalty and leaves it negative.** On the shipped `h1union`
   stack Arm A reads `pds` **−0.003017** raw (**−0.001125** overall, signs `--0-`) against
   **−0.003797** on `cd4f034` — the new source dilutes PerturbMatch by **20.5%** — and transfer
   cosine **−1.7895%** relative, **8.9× outside** 204's ±0.20% resolution floor.
6. **Both of Arm A's costs are DISTRIBUTED, which makes the negative more credible, not less.**
   Leave-top-3-out retention: cosine **90.1%**, `pds` **112.0%** (>100% — dropping the three
   largest movers makes the loss *worse*).
7. **Priced only on shipped-stack measurements, the arm is negative.** 218/223/225's rate card —
   §8.4's rule that no member is carried at zero — gives **−0.001517** overall.
8. **Priced 214's way, 0 of 12 (pool × transfer factor) cells clear +0.0089**, and 8 of the 12 are
   refuted on arithmetic. The best unrefuted cell is **+0.002553**, **3.5× below the floor**.
9. **🔴 And 230's noise floor closes the whole channel, not just this arm.** The board's own
   unexplained `reach` spread among digest-identical arms is **sd 0.004502 overall**. The **entire
   remaining `reach` prize — 7th of 656 to world first — is +0.002689 = 0.60 sd.** *Becoming the
   best `reach` team on earth is worth less than one standard deviation of this board's noise on
   that member.* 214's published range is **1.02–1.95 sd**; the best unrefuted cell here is
   **0.57 sd**.
10. **The half of Arm A this board can measure is NEGATIVE; the half that would be positive is
   0.81 sd.** `pds` and `mse` are deterministic (`rc_hybrid` reproduced both to sixteen
   significant figures) and they say **−0.001517**. The `reach`/`nmae` half is unmeasurable at
   this precision. **Per 230's rule: not cleared, no slot.**

**219 §2.5 was right that rule 13 double-charged 214, and it does not matter.** Repairing the gate
does not save the arm: it dies on an outside bound the gate never used, on its own priced channels
re-read against the stack that actually ships, and on a noise floor larger than the entire prize.

---

## 0. Phase 0 — the reconciliation against 220's bound

`phase0_231.py`, `work/phase0_231.json`. **Seven gates, all PASS.**

| gate | got | expected | |
|---|---:|---:|---|
| KA0-c `reach` raw→scaled slope, fitted on all 656 rows | **1.1291872** | 220's 1.12927 | PASS (max resid 3.1e−03) |
| KA0-e the same slope on **our own rows only**, 45 pairs | **1.1285988** | 220's 1.12927 | PASS (range 1.1073–1.1731) |
| KA0-a1 220's `ours` raw `reach` | 0.2498832 | 0.24983 | PASS (220 quoted the live API) |
| KA0-a2 220's field best (`luka`) | 0.2908353 | 0.2908 | PASS |
| KA0-b 220's own conversion 0.0410 × 1.12927 / 6 | 0.0077110 | 0.0077 | PASS |
| KA0-d AGENTS §4.1c's board table (cqawesome 0.22431 / ours 0.20492) | both | both | PASS |

### 0.1 The board has moved under 220, and the move is the whole finding

| snapshot | our model string | our raw `reach` | scaled | field best |
|---|---|---:|---:|---:|
| `2026-09-04T0013Z` | `Peer Review Rejected This` (**`rc_callset`**) | **0.24988** | +0.19290 | `luka` 0.29084 |
| `2026-09-04T0313Z` | `The Ribosome Has Filed A Grievance` (**`cd4f034`**) | 0.27330 | +0.21930 | 0.29084 |
| **`2026-09-04T0524Z`** | **`The Nucleus Declined To Comment` (`h1union`)** | **0.276545** | **+0.222939** | 0.29084 |

**220 wrote its bound on the 0013Z row.** Two submissions later our `reach` is **+0.026662 raw**
higher and we are **7th of 656** on that member. The remaining prize:

| | raw gap | overall |
|---|---:|---:|
| 220's bound, as written (on `rc_callset`) | +0.04097 | **+0.007711** |
| **the same bound today (on `h1union`)** | **+0.014290** | **+0.002689** |

**65% of the ordering prize 220 bounded has already been collected**, and — this is the part that
decides Arm A — it was collected by two *direction* arms (`cd4f034`, `h1union`), not by an
ordering arm. `rc_callset`'s call-set arm and the `cd4f034` dose between them account for
`0.24988 → 0.27330`; `h1union` added `+0.003246` more.

### 0.2 What Arm A's claim implies for our board `reach`

Chain, explicit: `implied board raw reach = ours + (h1-test raw Δ) × (transfer factor) ÷ (deflator)`.

| variant | Δ overall | Δ board raw `reach` | **implied** | % of remaining gap | |
|---|---:|---:|---:|---:|---|
| `wald89` × 1.290, no deflator — **214's headline `+0.0088`** | +0.006717 | +0.03569 | **0.31223** | **250%** | **REFUTED** |
| `wald89` × 1.290 ÷ 1.533 — **214's `+0.0054`** | +0.004381 | +0.02328 | **0.29983** | 163% | **REFUTED** |
| `native` × 1.290, no deflator | +0.005601 | +0.02976 | 0.30631 | 208% | **REFUTED** |
| `native` × 1.290 ÷ 1.533 | +0.003653 | +0.01941 | 0.29596 | 136% | **REFUTED** |
| `wald89` × 1.000, no deflator | +0.005207 | +0.02767 | 0.30421 | 194% | **REFUTED** |
| `native` × 1.000 ÷ 1.533 | +0.002832 | +0.01505 | 0.29159 | 105% | **REFUTED** |
| `wald89` × 0.4277, no deflator | +0.002227 | +0.01183 | 0.28838 | 83% | ok |
| **`wald89` × 0.4277 ÷ 1.533 — 214's `+0.0046`** | +0.001453 | +0.00772 | 0.28426 | 54% | ok |
| `native` × 0.4277, no deflator | +0.001857 | +0.00987 | 0.28641 | 69% | ok |
| `native` × 0.4277 ÷ 1.533 | +0.001211 | +0.00644 | 0.28298 | 45% | ok |

**8 of 12 variants are refuted on arithmetic**, including **every variant at the ×1.290 factor 214
used for its headline and every variant at unity.** At the headline the arm would put us **1st in
the world on `reach`** — from 7th, by a margin (0.31223 vs 0.290835) **1.5× larger than the entire
current 1st-to-7th spread**.

Only 164's **0.4277** — a *spread-policy* factor, measured on a different arm class — leaves the
claim inside the world's measured range. §1 shows that factor is one of four mutually inconsistent
values and that **none of them belongs to Arm A's class.**

> This is not a claim that raw `reach` 0.312 is impossible — the replicate anchor is ≈0.97, so
> there is headroom in principle. It is the outside view, stated numerically: **one component
> removal cannot plausibly be worth 2.5× the entire distance between us and the best `reach` team
> among 656.** The mandate pre-committed to reading it that way.

---

## 1. The `reach` transfer factor, re-derived — **it is not identified for this arm class**

`transfer231.py`, `work/transfer231.json`. 219 called this *"the whole arm"*. **Fifteen gates, all
PASS**, including reproductions of both published values.

| gate | got | expected | |
|---|---:|---:|---|
| KA3 all eight board entries matched by id | 8/8 overall to 5e−06 | ledger | PASS |
| KA1c **185's 1.290 reproduces from its own stated inputs** | **1.2927** | 1.290 | PASS |
| KA2c **the ledger §4's a8/a15 0.390 reproduces** | **0.389962** | 0.390 | PASS |
| KA4 the ledger's own scaling check (a8 raw→scaled) | 0.124216 | 0.124349 | PASS |

### 1.1 The trajectory, with the h1-test value paired to each board row

| arm | date | overall | board raw `reach` | h1-test `reach` | class |
|---|---|---:|---:|---:|---|
| a8 | 08-27 | 0.1543607 | 0.189003 | 0.20291 | emission (spread policy) |
| a15 | 08-30 | 0.1740865 | 0.181173 | 0.18283 | profile (CD4 residual) |
| `rc_hybrid` | 09-03 | 0.1834190 | 0.229922 | 0.21681 | emission-law, profile-PRESERVING |
| `rc_callset` | 09-03 | 0.1918638 | 0.249883 | 0.223341 | emission-law, profile-PRESERVING |
| **`cd4f034`** | 09-04 | 0.2049196 | **0.273299** | **0.223341** | **PROFILE (CD4 dose)** |
| **`h1union`** | 09-04 | **0.2118310** | **0.276545** | — | **PROFILE (new source)** |

### 1.2 The factor, on every pair where both sides are measured

| pair | h1-test Δ | board raw Δ | **factor** | arm class |
|---|---:|---:|---:|---|
| a8 − a15 | +0.02008 | +0.007830 | **0.3900** | emission / spread policy |
| `rc_hybrid` − a15 twin mean | +0.03771 | +0.048749 | **1.2927** | emission-law preserving — **185's 1.290** |
| `rc_hybrid` − a15 single run | +0.03398 | +0.048749 | **1.4346** | *same pair, other baseline* |
| `rc_callset` − `rc_hybrid` | +0.00653 | +0.019961 | **3.0563** | emission-law preserving |
| **`cd4f034` − `rc_callset`** | **+0.00000000** | **+0.023416** | **UNDEFINED** | **PROFILE — Arm A's class** |
| `h1union` − `cd4f034` | *unmeasurable* | +0.003246 | — | PROFILE |

**The error bar:** geometric mean **1.2193**, geometric sd **×/2.343**, 1σ band **[0.5205,
2.8564]**, range **7.84×** — plus one undefined point. 164's 0.4277 and 185's 1.290 are both
inside this set and **3.02× apart**. *Choosing between them changes Arm A's headline by 3×.*

**Two facts that close the question:**

1. **The `cd4f034` pair divides by zero, and that is a finding, not a nuisance.** 214 measured that
   `cd4f034`'s and `rc_callset`'s h1-test effect matrices are **identical, max |diff| 0.000e+00**
   (production's median-based sizing returns fraction 0 on a surface with 4 of 85 targets covered).
   h1-test therefore predicted **exactly 0.000000** for an arm that delivered **+0.023416 board raw
   `reach` = +0.004407 overall** — **1.6× the whole ordering movement `rc_callset` produced.** The
   map from h1-test to board `reach` is not a scalar; on this arm class it is not even injective.
2. **`h1union` cannot be scored on h1-test at all.** `data/donors/h1_2025` covers **85/85** of
   h1-test's targets — AGENTS §6.4 trap 1, *"pure oracle leakage"*. The currently shipped stack
   contains a source the only `reach` instrument we own cannot carry.

**Answer to 219's question 1: for a profile-class / donor-blend arm there are ZERO usable
calibration pairs, and the surface that would provide one cannot represent the shipped stack.**

---

## 2. The `pds` re-read on the shipped stack — 219's question 2

`pds231.py`, `work/pds231.json` (51 s, peak RSS 5.05 GiB). 219 asked for a re-read *"on
`cd4f034`"*; **208's rider rule** (*"re-price every rider against the current incumbent, never
against the incumbent it was measured on"*) makes the correct baseline **`h1union`**.

**Instrument:** 162's PLOLO folds, decomposed by 200's `folds200`, assembled exactly as 225's
`harness225.Fold` does it (which is 222's shipped construction), scored by 218's `score_effects`.
Both the CD4 fraction and the union fraction are **re-solved under each `pm_blend`**, because
production solves them from the whitened consensus.

### 2.0 Seventeen known-answer rows, all PASS

| gate | got | expected | |
|---|---:|---:|---|
| KA1a 4-fold mean, incumbent `cd4f034` | **0.740255** | 214 §1.2 f=0.34 | PASS |
| KA1b 4-fold mean, Arm A | **0.736458** | 214 §1.2 | PASS |
| KA1c 214 §1.2's PerturbMatch `pds` value at f=0.34 | **0.003797** | 0.003797 | PASS |
| KA1d/e 214 §1.2's cosines | 0.024862 / 0.024350 | identical | PASS |
| KA1f 214 §2.1's Arm A Δcos % | −2.0599 | −2.0512 | PASS |
| **KA2a 222/218's `h1union` rel MACRO cosine, 3 cancer** | **0.076803** | **0.0768028** | **PASS** |
| **KA2b 222/218's `h1union` Δ`pds`, 3 cancer** | **0.004887** | **0.0048867** | **PASS** |
| KA2c the same, 4 folds | 0.009072 | 0.0090720 | PASS |
| KA2d **all four per-fold Δ`pds`** (`diag218.json`) | 0.009516 / 0.002072 / 0.003072 / 0.021628 | identical | **PASS** |
| KA3 non-vacuity, all four folds | `frac_cd4` 0.316–0.376, `frac_union` 0.190–0.225 | > 0 | PASS |

> `score_effects` renders at **m = ∞** (no integer rounding), which is 214 §2.1's *"(reference) no
> rounding at all"* line (0.740255), not its `floor` line (0.735996). All comparisons below are
> like-for-like at m = ∞. In 214's `floor` convention the same measurement reads
> **−0.002842 → −0.002258**.

### 2.1 Arm A, both baselines, shipped renderer

| baseline | folds | `pds` | Δ raw | fold sd | scaled | **overall** | signs | cos | Δcos rel | cos_agg rel |
|---|---|---:|---:|---:|---:|---:|---|---:|---:|---:|
| `cd4f034` | 4 | 0.740255→0.736458 | **−0.003797** | 0.003516 | −0.008399 | **−0.001415** | `--0-` | 0.024862→0.024350 | **−2.0599%** | −2.2921% |
| **`h1union`** | **4** | **0.747859→0.744842** | **−0.003017** | 0.002958 | −0.006674 | **−0.001125** | `--0-` | 0.027408→0.026917 | **−1.7895%** | −1.8931% |
| `cd4f034` | 3-cancer | 0.780577→0.778216 | −0.002360 | 0.002482 | −0.005221 | −0.000880 | `0--` | | −1.6250% | −1.8197% |
| `h1union` | 3-cancer | 0.784329→0.782646 | −0.001683 | 0.001561 | −0.003723 | −0.000627 | `0--` | | −1.4664% | −1.5781% |

Per fold (4-fold, `h1union`): `orion_hct116` −0.001965 · `orion_hek293t` −0.003084 ·
`k562_gwps` **+0.000000** · `kolf2` −0.007021. `k562_gwps` is the `0` — with K562 held out
PerturbMatch has no slot to occupy, so that fold is a **structural tie, not a null**, exactly as
214 found. **`kolf2`, the leak-free fold, is again the largest single mover.**

**The dilution is real and it is 20.5%.** Adding the H1-union source shrinks Arm A's `pds` penalty
from −0.003797 to −0.003017 and its cosine loss from −2.06% to −1.79%. **It does not change the
sign of either.** The cosine loss is still **8.9× outside 204's ±0.20% resolution floor** — Arm A
is a *measurable* direction loss on the stack that ships.

### 2.2 Arm A is not strictly single-variable in the delivered stack

Dropping PerturbMatch shrinks the whitened consensus, so both post-whitening residuals size **up**:

| fold | CD4 fraction | union fraction |
|---|---|---|
| `orion_hct116` | 0.370782 → 0.383442 (**+3.414%**) | 0.222539 → 0.230137 (+3.414%) |
| `orion_hek293t` | 0.376056 → 0.388512 (+3.312%) | 0.225247 → 0.232708 (+3.312%) |
| `k562_gwps` | 0.373761 → 0.373761 (+0.000%) | 0.223670 → 0.223670 (+0.000%) |
| `kolf2` | 0.316408 → 0.325394 (+2.840%) | 0.189583 → 0.194968 (+2.840%) |

**Recorded because rule 5 cares:** the "one flag" moves the CD4 dose and the H1-union dose by ~3%
as a side effect. This is inside 216's dose-turnover band (`f* = 1.018`) but it is not zero.

---

## 3. The h1-test deflator, re-measured on the shipped stack

`disp231.py`, `work/disp231.json` (12 s, 5.18 GiB). 214's instrument verbatim. **Nine gates, all
PASS**, including the reproduction of 214's own deflator to **1.53319 against 1.533**.

| surface / stack | ‖ΔE‖/‖E‖ | row median | row-cos median | rows moved |
|---|---:|---:|---:|---:|
| h1-test (2 donors, 79/85 PM-covered) | **0.41486** | 0.36445 | 0.93314 | 85/85 |
| challenge / `cd4f034` (3 donors, 251/300 PM-covered) | **0.27059** | 0.24704 | 0.97081 | 300/300 |
| **challenge / `h1union` — SHIPPED** | **0.26651** | 0.24326 | 0.97142 | 300/300 |

| | value |
|---|---:|
| 214's deflator (h1-test / `cd4f034`) | **1.5332** |
| **231's deflator (h1-test / `h1union`)** | **1.5566** |
| extra dilution `h1union` adds | **×1.0153** |

**The deflator barely moves, and the reason is informative:** the H1 union carries only **77 of 300**
panel targets and sizes to fraction **0.187038** on the challenge panel, so it displaces little of
the challenge-side effect matrix — while on PLOLO's 2–3-donor held-out folds it carries much more
weight, which is why `pds` moved 20% and the displacement moved 1.5%. **Two surfaces, two honest
answers, and they are not in conflict.**

---

## 4. Rule 13's class test, on the shipped stack

`class231.py`, `work/class231.json` (13 s, 5.23 GiB). Two gates PASS (KA1 the displacement
reproduces §3; KA2 the delivered amplitude reproduces 214's `head214` at **1.02988 against
1.03099**). No package is built, so the digest is computed on the **reconstructed** delivered
totals `floor(N·L·profile)` — the reconstruction 214 validated at `floor(x) == T` on 1.0000 of
`cd4f034`'s own blocks.

| stack | digest identical | delivered median ‖row‖ | ratio | effect ‖row‖ ratio | panel counts |
|---|---:|---|---:|---:|---:|
| `cd4f034` | **0 / 900** | 3.321334 → 3.420578 | **1.02988** | 1.02780 | ±0.0000% |
| **`h1union`** | **0 / 900** | 3.340181 → 3.434345 | **1.02819** | 1.02780 | ±0.0000% |

**Class:** Arm A is a **donor-blend** change. The digest differs — the delivered profile moves —
but the **emission law is untouched**, which is rule 13's **third** row (*profile-class: operator,
donors, magnitude, CD4 dose*), where **169 falsified the toll at +0.0006**. **Toll = 0.**

**Amplitude:** **+2.819%** on the shipped stack, against the **1.0%** bound 210's mandate used and
the **2.03%** that made 210 flag the classification. **AMPLITUDE-CHANGING → floor +0.0089**, not
+0.0019. This is the single number that makes Arm A unshippable regardless of everything else
(§0.3, §6).

**230 has since landed and the question is settled: charge 0.** Rule 13's `−0.0146` is
**mis-keyed** — its own calibration points a21 and a22 are digest-identical to the incumbent on
900/900 blocks, so by the table's own test they belong in the zero-toll row, and *the row the toll
is charged to contains zero board points.* Arm A was already toll 0 by rule 13's third row; 230
removes the ambiguity entirely. §5.3 keeps the "both ways" table only to show the question is
**immaterial to this arm** — it fails by 3.5× either way.

---

## 5. The six-member price, on the shipped baseline

`price231.py`, `work/price231.json`. **Four gates PASS**, including a reproduction of 214 §2.4's
own `wald89` SUM to **7e−06** and the board-cosine inversion under §4.1c's validity guard
(`score_mse` = 0.17957 **> 0**, so 197's inversion is licensed; our board cosine is **0.099469**,
up from `cd4f034`'s 0.09018).

> *Recorded, not tuned away:* 214's published **`native`** SUM (+0.008791) reconstructs from its
> own per-member inputs only to **+0.008692**, a 1.1e−04 (1.1%) shortfall. The `wald89` row
> reconstructs exactly. It moves no verdict here.

### 5.1 ROUTE B — priced on shipped-stack measurements only

218/223/225's rate card, which is §8.4's rule that **no member is carried at zero**: five members
priced off the measured transfer-cosine move at 0.008007 overall per 0.01 of absolute cosine.

| | value |
|---|---:|
| five members off cosine: 0.008007 × (−0.000490 / 0.01) | **−0.000393** |
| `pds`, measured on the shipped stack | **−0.001125** |
| **ROUTE B TOTAL** | **−0.001517** |
| rule 13's amplitude floor | +0.0089 |
| | **FAILS — and of the wrong sign** |

**Every input here was measured on the stack that ships.** No transfer factor, no h1-test, no
superseded baseline.

### 5.2 ROUTE A — 214's decomposition, swept over every measured transfer factor

Deflator **1.5566** (§3), `pds` **−0.001125** (§2), `mse` = 096's claim **+0.000229** plus the
re-read cosine channel at 212's rate (**−0.000179**).

| pool | factor | `reach` | `nmae` | `fid` | `jac` | `pds` | `mse` | **SUM** | implied raw `reach` | % of gap | verdict |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| native | 0.3900 | +0.001232 | +0.002032 | −0.000040 | +0.000170 | −0.001125 | +0.000050 | **+0.002319** | 0.28309 | 46% | FAILS gate |
| native | **0.4277** | +0.001351 | +0.002032 | −0.000040 | +0.000170 | −0.001125 | +0.000050 | **+0.002438** | 0.28372 | 50% | FAILS gate |
| native | 1.2193 | +0.003850 | | | | | | +0.004938 | 0.29701 | 143% | **REFUTED** |
| native | **1.2900** | +0.004074 | +0.002032 | −0.000040 | +0.000170 | −0.001125 | +0.000050 | +0.005161 | 0.29819 | 151% | **REFUTED** |
| native | 1.4346 | +0.004530 | | | | | | +0.005618 | 0.30062 | 168% | **REFUTED** |
| native | 3.0563 | +0.009652 | | | | | | +0.010739 | 0.32783 | 359% | **REFUTED** |
| wald89 | 0.3900 | +0.001480 | +0.001874 | −0.000040 | +0.000170 | −0.001125 | +0.000050 | **+0.002410** | 0.28441 | 55% | FAILS gate |
| wald89 | **0.4277** | +0.001623 | +0.001874 | −0.000040 | +0.000170 | −0.001125 | +0.000050 | **+0.002553** | 0.28517 | 60% | FAILS gate |
| wald89 | 1.2193 | +0.004627 | | | | | | +0.005558 | 0.30113 | 172% | **REFUTED** |
| wald89 | **1.2900** | +0.004896 | +0.001874 | −0.000040 | +0.000170 | −0.001125 | +0.000050 | +0.005826 | 0.30256 | 182% | **REFUTED** |
| wald89 | 1.4346 | +0.005444 | | | | | | +0.006375 | 0.30547 | 202% | **REFUTED** |
| wald89 | 3.0563 | +0.011599 | | | | | | +0.012529 | 0.33818 | 431% | **REFUTED** |

**0 of 12 cells clear +0.0089. 8 of 12 are refuted on arithmetic.** The best unrefuted cell is
**+0.002553**, **3.49× below the floor**.

### 5.3 The toll — charge 0, following 230

Arm A's class carries **toll 0**, and 230 has now shown the `−0.0146` row is mis-keyed for
*everyone*. Charging it anyway, purely to show the question is immaterial here, the four unrefuted
rows become:

| pool | factor | toll 0 | toll −0.0146 |
|---|---:|---:|---:|
| native | 0.3900 | +0.002319 | **−0.000114** |
| native | 0.4277 | +0.002438 | **+0.000005** |
| wald89 | 0.3900 | +0.002410 | **−0.000023** |
| wald89 | 0.4277 | +0.002553 | **+0.000120** |

**Both ways, every row fails the +0.0089 floor by 3.5× or more. The toll question is immaterial to
this arm** — which is itself worth telling 230.

### 5.4 🔴 The same numbers in units of 230's board `reach` NOISE — and this is the verdict

`noise231.py`, `work/noise231.json`. **Three gates PASS** (230's scaled sd ÷ 6 = its overall sd;
its "2.4× the entire floor"; Phase 0's gap reloads).

230 measured the **unexplained spread of board Δ`reach` among DIGEST-IDENTICAL arms** at
**sd 0.027011 scaled = 0.004502 overall**. 219 says every point of Arm A lives in `reach`. So:

| quantity | overall | **in sd of the board's `reach` noise** |
|---|---:|---:|
| **the ENTIRE remaining `reach` prize (ours, 7th of 656 → world first)** | **+0.002689** | **0.60 sd** |
| rule 13's shape-only floor | +0.0019 | 0.42 sd |
| rule 13's amplitude floor — **Arm A's** | +0.0089 | 1.98 sd |
| 214's headline `+0.0088` | +0.0088 | **1.95 sd** — *and refuted on arithmetic* |
| 214's deflated `+0.0054` | +0.0054 | **1.20 sd** — *and refuted on arithmetic* |
| 214's low end `+0.0046` | +0.0046 | **1.02 sd** |
| **231's best UNREFUTED cell** (`wald89` × 0.4277) | **+0.002553** | **0.57 sd** |
| 231 **route B**, shipped-stack measurements only | **−0.001517** | **−0.34 sd** |

**Read the first row.** *Becoming the best `reach` team in the world — from 7th of 656 — is worth
**0.60 standard deviations** of this board's own noise on that member.* The whole channel now sits
inside its own instrument's noise. That is a statement about the channel, not about Arm A.

The `reach` **member alone**, per unrefuted cell: **0.27 / 0.30 / 0.33 / 0.36 sd**.

### 5.5 The split 230 asks for: the measurable half against the unmeasurable half

`pds` and `mse` are **deterministic** on this board — `rc_hybrid` reproduced both to **sixteen
significant figures**. The noise is specific to the ordering channel. So Arm A splits cleanly:

| half | content | value | precision |
|---|---|---:|---|
| **MEASURABLE** | `pds` on the shipped stack | **−0.001125** | deterministic |
| | transfer cosine −1.7895% rel (−0.000490 abs), through the five-member rate card | **−0.000393** | 8.9× the resolution floor |
| | **measurable half TOTAL** | **−0.001517** | |
| **UNMEASURABLE** | the `reach`/`nmae` half, best unrefuted estimate | +0.003628 | **0.81 sd** |

**The half of Arm A that this board can measure is NEGATIVE. The half that would be positive is
0.81 sd of the noise on its own channel.** 230's rule — *"if the arm's gain is entirely in `reach`
and lands under ~2 sd, say it does not clear, even if it beats the floor arithmetically"* — applies
here in its strongest form, because Arm A does not beat the floor either, and its measurable half
has the wrong sign. **Not cleared. No slot.** This is the same call 230's coordinator made against
208's Hamilton, applied consistently.


---

## 6. Falsifiers

**F1 — price all six on the shipped baseline, in the challenge frame, never carrying a member at
zero.** §5.1 does exactly this and is negative. §5.2 prices all six the other way and is negative
against the gate at every factor. **`pds`, `mse` and the transfer cosine are measured directly on
the shipped stack; `reach`/`fid`/`jac`/`nmae` are 214's h1-test numbers re-deflated by 1.5566** —
and §1.2 records why no better is possible: **the shipped stack cannot be represented on h1-test at
all** (`h1_2025` is 85/85 oracle leakage there).

**F2 — the digest, and 230 both ways.** §4: **0/900 identical** on both stacks; class = profile
(donor blend), toll 0 by rule 13 row 3; reported with and without the −0.0146 charge in §5.3.

**F3 — leave-top-3-out retention, `reach` and cosine separately.**

| statistic | baseline | folds | full | leave-top-3-out | **retention** |
|---|---|---|---:|---:|---:|
| **cosine** | `h1union` | 4 | −0.0004905 | −0.0004419 | **90.1%** |
| cosine | `h1union` | 3-cancer | −0.0004735 | −0.0004441 | 93.8% |
| cosine | `cd4f034` | 4 | −0.0005121 | −0.0004625 | 90.3% |
| **`pds`** | `h1union` | 4 | −0.0030174 | −0.0033803 | **112.0%** |
| `pds` | `h1union` | 3-cancer | −0.0016829 | −0.0024721 | **146.9%** |
| `pds` | `cd4f034` | 4 | −0.0037970 | −0.0040784 | 107.4% |

**Both of Arm A's costs are distributed, and the `pds` cost is anti-concentrated.** Retention above
100% means the three largest-|Δ| targets contain ones where dropping PerturbMatch *helps*; remove
them and the penalty grows. This is the mirror image of 222's concentration caveat, and it points
the other way: **the negative is not three unlucky perturbations.**

*(`reach` itself cannot be leave-top-3-out'd on the shipped stack — see F1 — so the falsifier is
answered on the two channels the shipped stack can carry, and 220 §2.3 already reported the same
test for the `reach` side of the ordering family, where 60–111% of every published effect was one
perturbation.)*

**F4 — what does dropping PerturbMatch COST?** It is a retrieval component we ship deliberately,
and on the shipped stack it is buying, per §2.1:

| channel | what PerturbMatch buys | instrument |
|---|---:|---|
| `pds` | **+0.003017 raw** (+0.001125 overall), signs `--0-`, 3 of 4 folds | PLOLO, 4 folds × frames A/B/C |
| transfer cosine | **+1.7895% relative** (+0.000490 absolute), 8.9× the resolution floor | same |
| aggregate cosine | +1.8931% relative | same |
| `reach` (h1-test) | **−0.026 to −0.031 raw** — it *costs* `reach` | 214, `rc_callset` emitter |
| `nmae` (h1-test) | −0.0106 to −0.0115 raw — it *costs* `nmae` | 214 |
| `fid` / `jac` | ≈ 0 (−0.00011 / −0.00065) | 214 |

**So Arm A trades DIRECTION for ORDERING, and both sides of that trade are now sized.** §8.2
measures the direction channel at **110% of the entire remaining gap to first**; §0.1 measures the
whole ordering channel's remaining prize at **+0.002689**. AGENTS §4.1c's warning — *"Do not trade
ordering for direction. MT Zeng is the experiment showing what that costs"* — applies with the
labels swapped and the same conclusion: **the channel Arm A sells is the one that is worth more to
us than to anyone on the board, and the channel it buys is nearly exhausted.**

**F5 — the instrument is not selected on the surface it is scored on, and every gate is
independent of the arm.** Fifty-seven known-answer rows across seven scripts: 214's `pds` table,
218's four per-fold `d_pds` rows, 222's preregistration, 214's deflator, 214's amplitude, 185's
1.290, the ledger's 0.390, the ledger's scaling map, and 220's own bound. **All PASS.** Every
script aborts on failure.

---

## 7. What this changes for other experiments

* **To 219 (the closure audit), which handed this arm forward: §2.5's diagnosis is CORRECT and its
  conclusion does not follow.** Rule 13 *did* double-charge 214 — comparing a two-member composite
  against the RMS of channels a six-member estimate had already priced. **Repairing the gate does
  not resurrect the arm.** The arm dies on (a) an outside bound the gate never used, (b) a transfer
  factor that is unidentified for its class, and (c) both of its priced channels being negative
  when re-read against the stack that ships. **Row 5 of 219's shortlist is closed. The audit's own
  recommendation — that rows 3 and 4 were worth more than rows 1, 2 and 5 combined — is confirmed.**

* **220's bound needs restating, and the restatement is worth more than the bound was.** *"Matching
  the best team in the world on `reach` is worth +0.0077 overall"* was true of the `rc_callset`
  row. **On the shipped `h1union` row it is +0.002689** — we are **7th of 656** on that member, not
  29th. **And the prize was collected by direction arms, not ordering arms**, which strengthens
  220's own closing paragraph rather than weakening it.

* **A general gate, and it is arithmetic.** Rule 13's amplitude floor is **+0.0089**; the entire
  remaining `reach` prize is **+0.002689**. **The floor is 3.31× the whole prize.** *No
  amplitude-changing arm whose gain lives in `reach` can clear its own gate at any transfer factor,
  even at oracle.* Any future `reach` arm must either be **profile-preserving** (floor +0.0019,
  which the +0.002689 ceiling does still admit — barely) or carry a direction gain of its own.

* **To 185 and to AGENTS §11.1's rule-13 amendment: "the h1-test → board `reach` factor is not a
  constant" is now much stronger than when it was written.** It is **0.390 / 1.2927 / 1.4346 /
  3.0563 / undefined** across five consecutive board pairs, a 7.84× spread with a geometric sd of
  ×/2.343 — and the undefined one is the **profile class**, where h1-test read exactly zero for a
  +0.023416 board move. §11.1's *"do not apply a single transfer factor to `reach` across arm
  classes"* should be strengthened to **"for profile-class arms there is no factor at all; h1-test
  cannot see their `reach`."**

* **The shipped stack has outrun its own DE instrument.** `h1union` contains `h1_2025`, which is
  85/85 oracle leakage on h1-test. **No `reach`, `fid`, `jac` or `nmae` arm can be priced against
  the current incumbent on any surface we own.** This is a live blocker for 225/229's source
  programme, not just for Arm A, and it is the strongest argument yet for 221's cloud `cell-eval2`
  worker.

* **To 214, whose work stands.** Every number 214 reported reproduces here — its §1.2 `pds` table
  to 2e−05, its ×1.533 deflator to 1.53319, its delivered amplitude to 1.02988 vs 1.03099, its
  `wald89` six-member SUM to 7e−06. **214 chose the conservative number and was right to.** Its
  only correctable step is the borrowed 1.290.

---

## 8. What I did NOT measure, and why

| not measured | why | can it move the verdict? |
|---|---|---|
| Arm A's four DE members **on the shipped stack** | structurally impossible: `h1union` contains `h1_2025`, which covers 85/85 of h1-test's targets (AGENTS §6.4 trap 1). The only `reach` instrument we own cannot carry the incumbent | **No — and this is itself the §1.2 finding.** Route B (§5.1) prices all six without it, and is negative |
| a **challenge-side** `reach` for Arm A | needs a truth we do not have | no |
| a board point | rule 7 / the standing no-slot instruction; the arm fails its gate by 3.5× on every route | no |
| Arm A **composed with** an ordering arm | 220 closed the ordering-key family permanently; nothing is left to compose with | no |
| re-solving α for the +2.82% amplitude move | 211's standing rule — amplitude is closed, and re-solving would make this a two-variable arm (rule 5) | no; §5 states the cost rather than recovering it |
| 230's step-vs-scales toll | had not landed at 05:40Z | **no — §5.3 reports both ways and every row fails either way** |

---

## 9. Reproduction

```
uv run python phase0_231.py      # Phase 0: the reconciliation vs 220's bound   (3 s)
uv run python transfer231.py     # the reach transfer factor + error bar        (2 s)
uv run python pds231.py          # PLOLO pds/cosine re-read, both baselines     (51 s, 5.05 GiB)
uv run python disp231.py         # the h1-test deflator on the shipped stack    (12 s, 5.18 GiB)
uv run python class231.py        # rule 13's digest + amplitude class test      (13 s, 5.23 GiB)
uv run python price231.py        # the six-member price, both routes            (1 s)
uv run python noise231.py        # the price in units of 230's board reach noise (1 s)
```

`/usr/bin/time -l` maximum resident set size: `class231` **5.23 GiB**, `disp231` 5.18 GiB,
`pds231` 5.05 GiB — all far under the 10 GB budget. Total wall clock **~90 s**.

**INBOX-MEMORY compliance:** every heavy job was run **strictly one at a time**, never backgrounded
in pairs, and no sample size or falsifier was reduced to fit memory. Swap was checked before the
last launch (free 1,308 MB); the three ≥5 GiB jobs had already completed by then and the remainder
are sub-second arithmetic. No other agent's process was touched, and 229's production builds
(2.7 GB × 2 at the time) ran undisturbed throughout.

### Artefacts

| file | what it produces |
|---|---|
| `phase0_231.py` | §0 — the reconciliation, 7 gates → `work/phase0_231.json` |
| `transfer231.py` | §1 — the transfer factor on 5 board pairs, 15 gates → `work/transfer231.json` |
| `pds231.py` | §2 — Arm A's `pds`/cosine on both baselines, F3, 17 gates → `work/pds231.json` |
| `disp231.py` | §3 — the deflator on three surfaces, 9 gates → `work/disp231.json` |
| `class231.py` | §4 — the 900-block digest and the delivered amplitude, 2 gates → `work/class231.json` |
| `price231.py` | §5 — both pricing routes and rule 13's gate, 4 gates → `work/price231.json` |
| `noise231.py` | §5.4/5.5 — the price in units of 230's board `reach` noise, 3 gates → `work/noise231.json` |
| `work/make_submission_pm_cd4.SNAPSHOT.py` | byte copy of the production builder taken **2026-09-04T05:34:24Z** while 229 was editing it, sha256 `a8785cab8a6511e9809c9cd5716c2911c3b986baac8ed2c43af6e07fa4181ec8` (5,026 lines). **Read only; the live file was never modified.** |

### Integrity

Every script carries known-answer rows and aborts on failure; **57 gates, all PASS**. The PLOLO
fold caches and 200's `parts-*` caches were **copied into `work/` and every module's cache
directory re-pointed there** (`plolo.WORK`, `folds200.WORK`), so no shared cache was written:
`experiments/200-production-ablation/work/parts-challenge.npz` and `parts-h1test.npz` still carry
their 2026-09-03 mtimes, and `experiments/162-retrodicting-validator/work/` is untouched. No sealed
data opened, no `cell-eval2` run (none was running at start — checked), no build, no `vcc`
invocation, no `submissions/` write, and no file outside `experiments/231-arm-a-retry/` created or
modified. `src/` and `scripts/` were read only.

---

## VERDICT

**214's Arm A is refuted.** Its headline requires 250% of the entire remaining world gap on
`reach`; its transfer factor is unidentified for its own arm class; its amplitude floor is 3.31×
larger than the whole ordering prize that remains; and the half of it this board can actually
measure — `pds` and the transfer cosine, both deterministic to sixteen significant figures — is
**negative at −0.001517**. The half that would be positive is **0.81 sd** of the board's own
`reach` noise.

**Do not drop PerturbMatch. Do not spend a slot.**

**And the larger result, which outlives the arm:** with our shipped `reach` at 7th of 656, the
entire remaining prize on that member is **+0.002689 = 0.60 sd** of the board's noise on it.
**The ordering channel is now smaller than its own measurement error.** 220 closed the
ordering-*key* family; this closes the ordering *channel* — for amplitude-changing arms by
arithmetic against rule 13's floor, and for every arm by the noise. **219's row 5 is closed, and
the direction programme (§8.2, §4.1c) is the only place points remain.**
