# 215 — The leaders' board histories, decomposed into channels. **Their cosine was NOT there from the start: the two independent top pipelines each DOUBLED it between 2026-08-29 and 2026-09-03, and the board measured it.**

**Date:** 2026-09-04 · 17 board snapshots (2,208 distinct published observations, 644 live rows
pulled 2026-09-04) · **read-only intel.** No `vcc`, no build, no `cell-eval2`, no submission.
Nothing posted, nobody contacted, no account created, nothing joined, and **no query or public
field carried anything about our method** (`AGENTS.md` §14). Peak RSS < 0.6 GB. Everything
written is inside `experiments/215-leader-decomposition/`.

Scripts, each gated: `panel215.py`, `chan215.py`, `phase0215.py`, `phase1215.py`, `phase2215.py`,
`phase3215.py`, `valid215.py`, `sens215.py`. Logs and JSON in `work/`.

---

## 0. The answer in nine lines

1. **The instrument works, and it is calibrated on five arms whose channel we know for certain.**
   5 of 5 correctly classified, budget partitions the observed Δoverall exactly, and the three
   zero-magnitude known answers reproduce at 0.00σ. §1.
2. **It then validates out-of-sample on strangers.** On 13 DIRECTION-labelled arms belonging to
   other teams the two *independent* readouts agree at **ρ median +0.944** (theory: 1.000); on 7
   AMPLITUDE-labelled arms they disagree at **ρ median +0.010** (theory: 0.000). No free
   parameter was fitted to any of those rows. §1.4.
3. **`cqawesome`'s direction is not dark any more.** Their first informative row is `m0830b`
   (2026-08-30). From there to `m0903c` (09-03) their cosine floor goes **0.06427 → 0.13211,
   ×2.06**, and `pds` — which is exactly amplitude-free — corroborates it at **ρ = 0.94**. The
   amplitude explanation is excluded at **12.4× the conservative leak bound**. §3.
4. **206 §8's "the board carries no information about how Prima Mente got it" is half wrong, and
   the leak is their own ablation account.** `jw_13288` shrank two days before the rest of the
   cluster. Its informative window records the cluster's direction going **0.06333
   (`willow_34`, 08-29) → 0.13117 (`juniper_46`, 09-03)**, ρ = 1.13, amplitude excluded at 14.9×.
   **Both of the board's two independent top pipelines doubled their cosine in the same five days.**
5. **Nobody has ever appeared on this board already holding a leading direction.** Over 424 teams
   with a usable observation, the highest cosine anyone has *ever* carried on the earliest row we
   hold for them is **0.08616**, and among submissions #1–#6 it is **0.08614** — 65% of
   cqawesome's 0.13211. §4.1.
6. **The climbs are incremental, not acquisitions.** `cqawesome` climbed `pds` in **5 separately
   resolved steps**, the largest being **28%** of the total; `Hongjian Zhou` in ≥2 steps (largest
   43%); `Haystack` in 3 (largest 56%). No leader shows a single step carrying its direction. §4.2.
7. **`Aginglab.com` is the one clean "it was there from the start" case — and it is stuck at
   0.0989.** Their `pds` has moved **+0.0031 raw in eight days**, below the conservative leak bound
   (0.9×), while their cosine *floor* rose +0.0679. **93% of their apparent direction gain was
   amplitude.** Direction share of their whole climb: **−2%**. §2.
8. **`Haystack` confirms the two-channel split from outside.** `reach` **0.2050 at submission #3**,
   `0.2442` at #23 (+0.0392) while `pds` climbed +0.1052 and raw `mse` fell 6.8×. Ordering was
   structural on day one and never had to be bought. §4.3.
9. **The date.** The only public event in the window is the **Cell publication of the Marson/Zhu
   genome-scale CD4⁺ T-cell Perturb-seq atlas on 2026-08-28** (Crossref record created 08-28,
   deposited 08-29, CC BY-NC). Both leaders' doubling windows open **one and two days later**.
   The data itself has been public since the 2025-12-23 bioRxiv preprint, so this is a *visibility*
   event, not an availability event — but it is the only one, and 206 §6.3/§7-H1 arrived at the
   same source by four independent routes. **X-Atlas/Pisces is still 0 bytes** (re-verified today,
   `lastModified 2026-03-25T22:58:18Z`, 3 files), so it is still nobody's edge. §3.3.

> ### Verdict on Phase 3, stated against the mandate's three branches
>
> **NOT "flat and high from the start" — the opposite.** The leaders' direction arrived
> **inside the window the board can read**, in **many steps**, and it **doubled in five days**.
> That is the mandate's third branch inverted: the direction channel is not a wall, it is the
> channel that has moved fastest on this board in the last week, and two mutually independent
> hobbyist-shaped accounts both walked it.
>
> **This argues AGAINST spending the remaining 48 days only on D/E/F robustness.** It is the
> strongest evidence yet for 206's **H1** (the CD4 dose, shipped at effective fraction 0.0528
> because 062 selected it on `pds`), and it is now supported by a *trajectory* rather than by the
> `k_eff` arithmetic. 210's `cd4f034` and 207's dose sweep are the right arms.

---

## 1. Phase 0 — the instrument, and its five known-answer arms

### 1.1 What separates the channels

The whole decomposition rests on one asymmetry between two members:

| member | invariance | proof |
|---|---|---|
| **`pds`** | **exactly** invariant to the emitted order and call set | `rc_hybrid` and `rc_callset` reproduced a15's `pds_cosine` to **sixteen significant figures**, twice, while `reach` moved +0.0487 and +0.0200 raw (185, 202 BOARD-RESULTs) |
| **`pds`** | invariant to a global amplitude up to a measured leak | a10 → a11 cut α 1.0796 → 0.838 and moved raw `pds` by **+0.001503** while 197 measured the cosine **flat** (0.07554 → 0.07550) |
| **`mse`** | **not** invariant to amplitude | the same arm moved raw `mse` by **−0.442426** |

So `pds` is the amplitude-free, ordering-free readout of a real change in the prediction, and
`mse` is the amplitude-contaminated one. Their **ratio** is the discriminant:

```
dcos_pds = d(raw pds)  / 0.6479          amplitude-free, ordering-free
dcos_mse = -d(raw mse) / 2.074495        197's floor difference, amplitude-CONTAMINATED
rho      = dcos_pds / dcos_mse           1 = pure direction, 0 = pure amplitude, <0 = geometry
```

`2.074495` is 197's identity constant, not a fit. `0.6479` is the decomposition's **one fitted
constant**, taken from a14 → a15 (raw `pds` +0.003285 against 197's measured Δcos +0.00507).
§5 sweeps it; no conclusion turns on it.

The **amplitude significance threshold is adaptive**: `sigma_pds_eff = max(0.0004,
0.003397·|Δ raw mse|)`, where `0.003397 = 0.001503/0.442426` is the leakage measured on a10 → a11.
That single change is what lets the instrument tell a15's direction step from a11's amplitude step.

The **budget partitions the observed Δoverall exactly**. A counterfactual six-tuple is built by
moving `Δpds` through the direction channel's own member ratios, mapped through the board's
*clamped* scaling function (`mse` is clamped to scaled 0 above raw 0.98792, `nmae` floored at −6),
and the residual on `fid`/`reach`/`jac` is ORDERING while the residual on `mse`/`nmae` is split
between AMPLITUDE and GEOMETRY by the sign of the `nmae` residual against the `mse` residual
(amplitude: |Δnmae/Δmse| = 0.0036 measured; geometry: −0.056 measured on a19).

### 1.2 The confusion matrix — 5 of 5

| truth | DIRECTION | DIR+AMP | ORDERING | AMPLITUDE | GEOMETRY | NULL |
|---|---:|---:|---:|---:|---:|---:|
| **DIRECTION** (a14→a15) | **1** | 0 | 0 | 0 | 0 | 0 |
| **ORDERING** (rc_hybrid, rc_callset) | 0 | 0 | **2** | 0 | 0 | 0 |
| **AMPLITUDE** (a10→a11) | 0 | 0 | 0 | **1** | 0 | 0 |
| **GEOMETRY** (a15→a19) | 0 | 0 | 0 | 0 | **1** | 0 |

| arm | Δ`pds` | Δ`mse` | Δ`nmae` | Δ`reach` | ρ | called | known answer |
|---|---:|---:|---:|---:|---:|---|---|
| a10 → a11 (α 1.0796→0.838) | +0.00150 | −0.44243 | −0.00160 | −0.01007 | **0.01** | AMPLITUDE | 197 §1: cosine flat, −0.00004 |
| a14 → a15 (CD4 residual) | +0.00328 | −0.00897 | −0.00151 | +0.00056 | **1.17** | DIRECTION | 197 §1.2: +6.8% cosine |
| a15 → a19 (whitener ridge) | −0.00053 | +0.09116 | −0.00511 | +0.00263 | 0.02 | GEOMETRY | LOG 08-31: `mse` −0.0959 scaled, `nmae` **improved** |
| a15 → rc_hybrid | −0.0000000 | −0.0000000 | −0.00003 | **+0.04875** | — | ORDERING | `pds`/`mse` = a15's to 16 s.f. |
| rc_hybrid → rc_callset | **0.0000000** | **0.0000000** | −0.0000000 | +0.01996 | — | ORDERING | bit-identical `pds`/`mse` |

### 1.3 The budget on the same five arms (overall units; the four buckets sum **exactly** to the observed Δ)

| arm | direction | ordering | amplitude | geometry | Σ | observed |
|---|---:|---:|---:|---:|---:|---:|
| a10 → a11 | +0.00096 | −0.00295 | **+0.01392** | +0.00000 | +0.01192 | +0.01192 |
| a14 → a15 | **+0.00395** | +0.00000 | −0.00027 | −0.00000 | +0.00367 | +0.00367 |
| a15 → a19 | −0.00064 | +0.00133 | −0.00000 | **−0.01421** | −0.01352 | −0.01352 |
| a15 → rc_hybrid | −0.00000 | **+0.00934** | +0.00001 | +0.00000 | +0.00935 | +0.00935 |
| rc_hybrid → rc_callset | +0.00000 | **+0.00835** | +0.00000 | +0.00000 | +0.00835 | +0.00835 |

**Magnitude check** against 197's independently inverted cosines:

| arm | instrument | 197 | error | verdict |
|---|---:|---:|---:|---|
| a10 → a11 | +0.00232 | −0.00004 | **1.02σ** of its own leak floor | PASS — and the residual *is* the leak the model predicts |
| a14 → a15 | +0.00507 | +0.00507 | 0.00σ | **calibration, not a test** (0.6479 is fitted here) |
| a15 → rc_hybrid | −0.0000000 | +0.000001 | 0.00σ | PASS |
| rc_hybrid → rc_callset | +0.0000000 | 0 | 0.00σ | PASS |

### 1.4 🔴 The out-of-sample validation — on strangers' rows, with nothing fitted to them

For every arm-to-arm pair where **both** endpoints are informative (raw `mse` < 1.023611) and
`pds` sits in the band where the rate is calibrated, the two readouts are independent. Theory says
DIRECTION arms must land at ρ = 1 and AMPLITUDE arms at ρ = 0. Measured on 34 such pairs across
the whole board:

| subset | n | ρ median | range | theory |
|---|---:|---:|---|---:|
| **DIRECTION-labelled** | 13 | **+0.944** | [+0.17, +4.56] | **1.000** |
| **AMPLITUDE-labelled** | 7 | **+0.010** | [−0.031, +0.115] | **0.000** |

Whole-set regression through the origin: `dcos_pds = 0.240 · dcos_mse` — i.e. **only a quarter of
the field's `mse` motion is corroborated by `pds`.** That is 206 §4's α-shrink wave, measured from
a second direction and with a number attached.

**V2 — the amplitude channel's `pds` leak, measured on the field.** Over 51 AMPLITUDE-labelled
strangers' arms with |Δ`mse`| > 0.02, `|Δpds|/|Δmse|` has median 0.0064 and **p90 0.0236**. That
distribution is an *upper envelope* (an arm that mixes a little direction lands in this bucket
too), so p90 is used below as the **conservative** leak bound — 6.9× our own calibrated 0.0034.

---

## 2. Phase 1 — the per-team channel budget

Budgets are over each team's **clean** window (pairs whose endpoints both have raw `nmae` ≤ 1.05
and `fid` ≥ 0.45 — a broken or shakedown submission is excluded, not decomposed). The four
buckets sum to the observed climb in overall units.

| team | obs | first | last | **climb** | direction | ordering | amplitude | geometry | %dir | %ord | %amp | %geo | Δ raw `pds` | cos floor now |
|---|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **cqawesome** | 10 | 08-21 | 09-03 | **+0.1590** | **+0.1477** | +0.0069 | +0.0044 | +0.0000 | **93** | 4 | 3 | 0 | **+0.1947** | **0.13211** |
| **jw_13288** | 6 | 08-27 | 09-03 | +0.0604 | **+0.0594** | +0.0071 | −0.0061 | 0.0000 | **98** | 12 | −10 | 0 | +0.0495 | 0.13117 |
| **Hongjian Zhou** | 8 | 08-22 | 09-03 | +0.1083 | +0.0814 | −0.0113 | +0.0321 | +0.0060 | 75 | −10 | 30 | 6 | +0.1257 | 0.12986 |
| **MT Zeng** | 4 | 08-30 | 09-03 | +0.0676 | +0.0386 | +0.0045 | +0.0177 | +0.0068 | 57 | 7 | 26 | 10 | +0.0411 | 0.12749 |
| **Andrew Liu** | 4 | 08-30 | 09-03 | +0.0528 | +0.0213 | −0.0045 | **+0.0308** | +0.0051 | 40 | −8 | **58** | 10 | +0.0303 | 0.12710 |
| **Garrrbage** | 8 | 08-21 | 09-02 | +0.1107 | +0.0740 | +0.0207 | +0.0362 | −0.0202 | 67 | 19 | 33 | −18 | +0.1158 | 0.11465 |
| **Aginglab.com** | 4 | 08-25 | 09-03 | +0.0589 | **−0.0012** | +0.0225 | **+0.0375** | +0.0000 | **−2** | 38 | **64** | 0 | **−0.0046** | 0.09887 |
| **Slurm Shadi** | 6 | 08-27 | 09-03 | +0.0812 | +0.0532 | +0.0073 | +0.0141 | +0.0065 | 66 | 9 | 17 | 8 | +0.0836 | 0.08118 |
| **Haystack** | 9 | 08-22 | 09-03 | +0.0784 | +0.0695 | −0.0095 | +0.0182 | +0.0003 | 89 | −12 | 23 | 0 | +0.1052 | 0.08075 |
| **us** | 11 | 08-21 | 09-03 | +0.1087 | +0.0695 | +0.0145 | +0.0181 | +0.0065 | 64 | 13 | 17 | 6 | +0.1081 | 0.07221 |

Notes that matter for reading it:

- **`Aginglab.com` is the one flat-direction team on the board.** Their `pds` has moved −0.0046
  net since 08-25 and +0.0031 since their first informative row. Everything they have gained since
  is `mse` (amplitude, 64%) plus one `reach` step (`GeroAI_v9 → v11`, raw `reach` 0.1672 →
  **0.2674**). It is not classifiable as ordering — `nmae` moved −0.0294 and raw `mse` −0.86 in the
  same arm — and their `reach` has drifted back to 0.2570 since, against a field best of 0.2908
  (`luka`).
- **`Andrew Liu` and `MT Zeng` are amplitude-dominant** because they only *join* the visible board
  on 08-30, by which time the cluster's pipeline already carried its direction. Their budgets
  describe when they were observed, not when the method was built. 206 §3 is reproduced exactly.
- **`jw_13288`'s clean budget starts at `willow_34`** because `pm-r14-pds3o` (raw `mse` 101.77,
  `nmae` 2.04) is a broken row.
- **`Garrrbage`'s −18% geometry** is the `v10 → v11` pair, which moved `nmae` +0.0461 worse against
  `mse` −0.0395 better at pinned `pds` — the a19 signature.

### 2.1 `cqawesome`, arm by arm

| from → to | Δ# | Δ`pds` | Δ`mse` | Δ`nmae` | Δ`reach` | ρ | label | Δovr | dir | ord | amp |
|---|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|
| `test_water_v2` → `v3` | 1 | +0.0484 | +1.0309 | −0.0318 | +0.0732 | — | broken → working (excluded) | +0.1780 | | | |
| `v3` → `v6` | 3 | −0.0047 | −0.1459 | −0.0151 | +0.0022 | — | DIRECTION(pds) `*+v` | +0.0128 | −0.0030 | +0.0110 | +0.0048 |
| `v6` → `test_fire_v1` | 3 | +0.0099 | +2.1999 | −0.0033 | −0.0026 | — | UNRESOLVED `*+v` | +0.0206 | +0.0064 | +0.0146 | −0.0003 |
| `test_fire_v1` → `mp4` | 4 | **+0.0559** | +1.6199 | +0.0388 | +0.0152 | — | DIRECTION(pds) `*+v` | +0.0158 | +0.0358 | −0.0022 | −0.0178 |
| `mp4` → `E3a_ctrlprog` | 2 | **+0.0563** | −0.1335 | −0.0159 | +0.0278 | — | DIRECTION(pds) `*+v` | +0.0303 | +0.0360 | −0.0030 | −0.0028 |
| `E3a_ctrlprog` → `m0830b` | 4 | +0.0360 | **−4.7207** | −0.0340 | −0.0100 | 0.02 | AMPLITUDE `+v` | +0.0174 | +0.0230 | −0.0275 | +0.0219 |
| **`m0830b` → `m0902a`** | 5 | **+0.0391** | −0.1280 | −0.0091 | +0.0498 | **0.98** | **DIRECTION** | **+0.0580** | +0.0469 | +0.0131 | −0.0019 |
| `m0902a` → `m0902b` | 1 | +0.0002 | −0.0018 | −0.0003 | −0.0014 | 0.38 | NULL | +0.0011 | | | |
| `m0902b` → `m0903c` | 2 | +0.0019 | −0.0109 | +0.0011 | −0.0029 | 0.55 | DIRECTION+AMPLITUDE | +0.0030 | +0.0023 | +0.0004 | +0.0003 |
| **TOTAL (clean)** | | **+0.1947** | | | | | | **+0.1590** | **+0.1477** | +0.0069 | +0.0044 |

`*` = `pds` outside the [0.79, 0.88] band where the rate is calibrated · `+` = more than one
submission elapsed · `v` = at least one endpoint above raw `mse` 1.023611, so its floor is vacuous.

---

## 3. Phase 2 — the direction arms, with dates and bounds

### 3.1 MEASURED — both endpoints informative, so Δcosine is a difference of two 197 floors

Fourteen such arms exist on the entire board, over 707 usable arm-to-arm pairs.

| # | team | from → to | date | Δ# | **Δcos (floor)** | Δcos (`pds`) | ρ | cos from → to | Δovr |
|---:|---|---|---|---:|---:|---:|---:|---|---:|
| 1 | **cqawesome** | `m0830b` → `m0902a` | **2026-09-02** | 5 | **+0.0617** | +0.0603 | **0.98** | 0.06427 → **0.12599** | +0.0580 |
| 2 | MT Zeng | `Petrel1` → `Kestrel2` | 2026-09-03 | 2 | +0.0582 | +0.0492 | 0.85 | 0.06932 → 0.12749 | +0.0468 |
| 3 | athenen | `exp275-sparse` → `exp278-k562-g…` | 2026-09-03 | 1 | +0.0311 | +0.0382 | 1.23 | 0.00946 → 0.04061 | +0.0327 |
| 4 | Ibrahim Mansour | `CellSim` → `CellSim` | 2026-08-25 | 4 | +0.0189 | +0.0221 | 1.17 | 0.05502 → 0.07391 | +0.0155 |
| 5 | **jw_13288** | `willow_34` → `hazel_21` | **2026-09-01** | 5 | **+0.0167** | +0.0550 | 3.29 | 0.06333 → 0.08002 | +0.0342 |
| 6 | 王子木 | `G512 ExactTot…` → `G1002 XAtlas…` | 2026-08-27 | 2 | +0.0139 | +0.0289 | 2.08 | 0.00791 → 0.02184 | +0.0084 |
| 7 | Haystack | `GSA-015` → `GSA-018` | 2026-09-01 | 3 | +0.0108 | +0.0231 | 2.14 | 0.08246 → 0.09327 | +0.0150 |
| 8 | Yuzhu Shi | `Bamboo` → `Bamboo` | 2026-09-02 | 2 | +0.0099 | +0.0094 | 0.94 | 0.08541 → 0.09536 | +0.0041 |
| 9 | Patheia | `Cell Papa` → `Cell Papa` | 2026-09-02 | 2 | +0.0099 | +0.0094 | 0.95 | 0.08542 → 0.09536 | +0.0059 |
| 10 | black_handsome_egg | `啊小绳头V12` → `V19` | 2026-09-03 | 2 | +0.0099 | +0.0199 | 2.01 | 0.03906 → 0.04895 | +0.0215 |
| 11 | cqawesome | `m0902b` → `m0903c` | 2026-09-03 | 2 | +0.0053 | +0.0029 | 0.55 | 0.12686 → **0.13211** | +0.0030 |
| 12 | **us** | `a14-era` → `a23` | 2026-09-02 | 6 | +0.0050 | +0.0062 | 1.24 | 0.06789 → 0.07291 | +0.0029 |
| 13 | Aginglab.com | `GeroAI_v12` → `GeroAI_v19` | 2026-09-03 | 4 | +0.0041 | +0.0044 | 1.07 | 0.09478 → 0.09887 | +0.0017 |
| 14 | Dongmei Cui | `dc-12` → `dc-18` | 2026-09-02 | 6 | +0.0028 | +0.0399 | 14.5 | 0.04674 → 0.04949 | +0.0337 |

**Row 2 is an artefact and must not be read as a direction arm.** `Petrel1` is 206 §3.2's failed
side-arm. Against MT Zeng's own `Kestrel` line the honest pair is `Kestrel4 → Kestrel2`: Δcos floor
+0.03484 but **ρ = 0.15**, Δ`pds` +0.0034 — amplitude, exactly as 206 measured.

### 3.2 The two windows that carry the whole finding, stated endpoint to endpoint

| team | window | dates | Δ# | cos floor | Δ`pds` | conservative leak bound | ×bound | ρ | verdict |
|---|---|---|---:|---|---:|---:|---:|---:|---|
| **cqawesome** | `m0830b` → `m0903c` | 08-30 → 09-03 | 8 | **0.06427 → 0.13211 (×2.06)** | **+0.0412** | 0.00332 | **12.4** | **0.94** | **DIRECTION, amplitude excluded** |
| **jw_13288** | `willow_34` → `juniper_46` | 08-29 → 09-03 | 9 | **0.06333 → 0.13117 (×2.07)** | **+0.0495** | 0.00332 | **14.9** | **1.13** | **DIRECTION, amplitude excluded** |
| Hongjian Zhou | `pm-r14b` → `pm-r14f` | 09-01 → 09-03 | 4 | 0.09764 → 0.12986 | +0.0019 | 0.00158 | 1.2 | 0.09 | not resolvable vs amplitude |
| Andrew Liu | `AL-4746` → `AL-1466` | 09-01 → 09-03 | 4 | 0.09725 → 0.12710 | +0.0035 | 0.00146 | 2.4 | 0.18 | not resolvable vs amplitude |
| MT Zeng | `Kestrel4` → `Kestrel2` | 09-01 → 09-03 | 4 | 0.09265 → 0.12749 | +0.0034 | 0.00170 | 2.0 | 0.15 | not resolvable vs amplitude |
| **Aginglab.com** | `GeroAI_v11` → `v19` | 08-26 → 09-03 | 5 | 0.03096 → 0.09887 | **+0.0031** | 0.00332 | **0.9** | **0.07** | **not resolvable — 93% of the floor gain is amplitude** |
| Haystack | `GSA-015` → `GSA-023` | 08-30 → 09-03 | 7 | 0.08246 → 0.08075 | +0.0039 | 0.00008 | 46 | −3.48 | direction FLAT within ±0.007 |
| us | `a14`-era → `rc_callset` | 08-30 → 09-03 | 9 | 0.06789 → 0.07221 | +0.0033 | 0.00021 | 15.5 | 1.17 | DIRECTION, +0.0043 |

**This is the correction to 206 §8 item 1.** The cluster's direction *is* readable, through the
account that shrank first. `jw_13288` reached raw `mse` 0.8922 on **08-29**, two days before
`Hongjian Zhou`, `Andrew Liu` and `MT Zeng`; its subsequent rows record the cluster's direction
doubling. What is *not* readable is anything before 08-29, when their raw `mse` ran 2.24–101.77.

### 3.3 The date, and what our own history holds there

| date | the board | us |
|---|---|---|
| **2026-08-28** | — | a10, a11 (α 1.0796 → 0.838) |
| 2026-08-29 | **`jw_13288` `willow_34` — the cluster's first informative row, cos floor 0.06333** | a12 (tiered spread), a13 |
| **2026-08-30** | **`cqawesome` `m0830b` — their first informative row, cos floor 0.06427** | a14, **a15 — the CD4 residual, +6.8% cosine, effective fraction 0.0528** |
| 2026-08-31 | | a19 (whitener ridge, −0.01352), a15 restore |
| 2026-09-01 | `jw` `hazel_21` 0.08002; PM cluster shrinks (`pm-r14b` 0.09764, `AL-4746` 0.09725, `Kestrel4` 0.09265) | a21, a22 |
| 2026-09-02 | **`cqawesome` `m0902a` 0.12599**; `jw` `cedar_90` | a23 (mass fix) |
| 2026-09-03 | **`cqawesome` `m0903c` 0.13211; `jw` `juniper_46` 0.13117** | abc-knee, **rc_hybrid** (+0.00933), **rc_callset** (+0.00844) |

**The public event, verified today, read-only:**

- **`Genome-scale perturb-seq in primary human CD4⁺ T cells…`, Cell, DOI
  `10.1016/j.cell.2026.08.002`.** Crossref record **created 2026-08-28**, deposited 2026-08-29,
  issued 2026-08, licence **CC BY-NC**. This is 206 §6.3's Marson/Zhu atlas — the lineage-matched
  donor for context A (Jurkat, exp 041), 284–291/300 panel coverage, which we hold on disk and
  ship at effective fraction **0.0528**.
- **X-Atlas/Pisces remains an empty repository**: `usedStorage: 0`, 3 files
  (`.gitattributes`, `LICENSE.md`, `README.md`), `lastModified 2026-03-25T22:58:18Z`. 206's check
  reproduces exactly five days later. Xaira's public description of Pisces does list *resting and
  CD3/CD28-activated Jurkat T lymphoma cells* — an H1-shaped lever if it ever lands — but it has
  not landed.
- No other release, preprint or repository event was found in the 08-28 → 09-03 window; 206's
  arXiv/bioRxiv/HF/GEO/GitHub sweep of 09-03 is not superseded.

**How much weight this carries.** The coincidence is of *dates*, and the data was already public
(bioRxiv 2025-12-23, CZI Virtual Cell Models portal). What is new is that (a) both leaders'
doubling windows open within 24–48 h of the formal publication, (b) the only board-published price
for that source is GISL's **+0.01325 overall = +0.0363 cosine at 197's rate**, about half of the
+0.068 both leaders gained, and (c) our own +6.8% at 5% weight is the *only* resolvable cosine gain
in eight of our board arms. Four independent pointers, one source. **This is a hypothesis with a
date, not a proof.**

### 3.4 BOUNDED and UNREADABLE arms — stated, per the mandate

- **Bounded (7 arms):** start vacuous, end informative. The endpoint cosine is a *floor*; the move
  itself is **not measurable**. Largest: `peng yu` `Sun20260903-c → MODEL` (09-03, raw `mse` 5.560 →
  0.8971, cos ≥ 0.06097); `Aubree` `V5 → V9` (09-02, ≥ 0.05569); `Ibrahim Mansour` (08-22,
  ≥ 0.05057); `MetaLearner` (08-27, ≥ 0.04080); `huanglinn` (08-30, ≥ 0.03739);
  `PumpkinChicken` (08-25, ≥ 0.03139); `Saber` (08-31, ≥ 0.00904).
- **Unreadable (178 arms):** end with raw `mse` ≥ 1.023611. **The board carries no cosine
  information for these at all** — the floor is vacuous and the `pds` → cosine rate is calibrated
  only inside `pds` ∈ [0.79, 0.88], which only 6 of the 178 satisfy. They are reported as raw
  `pds` moves and are **not** converted to cosine anywhere in this document. This class contains
  every arm of `cqawesome`'s first 18 submissions and every arm of `Hongjian Zhou`'s first 10.

### 3.5 A by-product: the ORDERING channel is still almost entirely unused

Under the strict definition (`pds`, `mse` and `nmae` all pinned at the seed floor while a DE member
moves — the exact `rc_hybrid`/`rc_callset` signature), **only 9 arms exist in 707 usable board
pairs, and two of them are ours.**

| team | arm | date | Δovr | Δ`reach` | Δ`jac` | Δ`fid` | Δ`pds` |
|---|---|---|---:|---:|---:|---:|---:|
| **Vu Mai** | `reprank-tau200` → `repequal-restore` | 09-02 | **+0.0251** | −0.0057 | +0.0003 | **+0.0463** | **0.0000** |
| **us** | a15 → `rc_hybrid` | 09-03 | +0.0094 | **+0.0526** | +0.0001 | −0.0001 | −0.0008 |
| **us** | `rc_hybrid` → `rc_callset` | 09-03 | +0.0084 | +0.0200 | **+0.0102** | −0.0001 | 0.0000 |
| IsoCell | `m4` → `m7` | 08-23 | +0.0030 | +0.0187 | −0.0001 | −0.0009 | −0.0002 |
| Happy | `HP-027` → `HP-030` | 09-03 | +0.0000 | +0.0157 | −0.0043 | −0.0017 | 0.0000 |
| Vu Mai | `repequal-t2000` → `reprank-tau200` | 09-02 | **−0.0251** | +0.0057 | −0.0003 | **−0.0463** | 0.0000 |

*(6 of the 9 shown; the other three — `Yuzhu Shi`, `IsoCell` `m2→m4`, `Slurm Shadi` — are all
negative or under +0.001.)*

**`Vu Mai` ran a profile-preserving emission A/B on the board and threw the result away.** Three
consecutive submissions at `pds` = 0.756373 and raw `mse` = 0.9181, **bit-identical**: a
rank-ordered emission cost them **−0.0464 raw `fid`** while gaining +0.0057 `reach`, and they
reverted. That is a measured price for the channel's `fid` leg from outside our own instruments,
and it is the risk our construction avoided (`rc_hybrid` moved `fid` by **+0.00019**).

---

## 4. Phase 3 — the three tests that decide it

### 4.1 T1 — did anyone start with a high direction? **No. Nobody. Ever.**

Over the 424 teams with at least one usable observation, only 16 were already informative on the
earliest row we hold for them with a cosine floor above 0.03. Ranked:

| team | earliest observed row | # | cos floor then | `pds` then | cos floor now |
|---|---|---:|---:|---:|---:|
| Leo | 2026-09-01 `(CPA+pcqdq4)` | 12 | **0.08616** | 0.8186 | 0.08798 |
| ITHACA ZERO | 2026-09-01 `M42-hy` | **4** | **0.08614** | 0.7910 | 0.08614 |
| huhuhu | 2026-09-01 `(cellflow+8c8vwp)` | 12 | 0.08611 | 0.8188 | 0.08899 |
| GZ YAI | 2026-09-01 `(STACK+x0hneb)` | **4** | 0.08576 | 0.8179 | 0.08601 |
| Xinping Yu | 2026-09-01 `v25e` | 4 | 0.07862 | 0.7776 | 0.07200 |
| Zhenxing Fang | 2026-09-03 `v32a1g` | **2** | 0.07303 | 0.8061 | 0.07303 |
| jw_13288 | 2026-08-29 `willow_34` | 6 | 0.06333 | 0.7838 | **0.13117** |

- **Highest cosine anyone has ever carried on their earliest observed row: 0.08616.**
- **Restricted to submissions #1–#6: 0.08614** (`ITHACA ZERO`, 2026-09-01).
- Both figures are **65% of cqawesome's 0.13211**, and four of the top five are the
  `lingshu`/`CPA`/`STACK`/`cellflow` wrapper cluster 206 §5.4 already separated from the leaders.
- All four leaders are observed from an early submission number — `Hongjian Zhou` from **#1**
  (`effect-transfer-v1-k562gw`, 08-22, `pds` 0.7208, raw `mse` 16.23), `Garrrbage` from #1,
  `cqawesome` from #2, `Haystack` from #3, `Aginglab.com` from #4 — so for them "we missed their
  early history" is not available as an explanation.

### 4.2 T2 — one step, two, or many?

| team | clean pairs | resolvable `pds` steps | Σ positive Δ`pds` | top-1 share | top-2 | top-3 | reading |
|---|---:|---:|---:|---:|---:|---:|---|
| **cqawesome** | 8 | **5** | +0.1994 | **28%** | 56% | 76% | **MANY steps (incremental)** |
| **Hongjian Zhou** | 6 | 2 | +0.1259 | 43% | 77% | 98% | **MANY steps (incremental)** |
| **Haystack** | 8 | **3** | +0.1189 | 56% | 74% | 86% | **MANY steps (incremental)** |
| MT Zeng | 3 | 3 | +0.0695 | 54% | 100% | — | two steps |
| Slurm Shadi | 5 | 3 | +0.0978 | 69% | 100% | 100% | two steps |
| jw_13288 | 3 | 2 | +0.0495 | 72% | 99% | 100% | one step *(see caveat)* |
| Garrrbage | 4 | 2 | +0.1161 | 76% | 95% | 100% | one step *(see caveat)* |
| Andrew Liu | 3 | 1 | +0.0303 | 89% | 100% | 100% | one step *(see caveat)* |
| **Aginglab.com** | 3 | **1** | **+0.0031** | 90% | 100% | — | **FLAT since first observation** |
| us | 10 | 3 | +0.1115 | 79% | 95% | 99% | one step (v3 → the whitening operator) |

**The caveat is load-bearing and I will not hide it: "one step" is unfalsifiable at this snapshot
cadence.** The board is sampled 17 times over 14 days, so a pair separated by 5 submissions could
be five small steps. **Only "many steps" is a positive finding**, because it requires several
*separately observed* increments; and only "flat" is a positive negative. So:

- **Demonstrably incremental:** `cqawesome` (5 resolved steps), `Hongjian Zhou` (2 + 2 unresolved
  gaps), `Haystack` (3).
- **Demonstrably flat:** `Aginglab.com`.
- **Nobody is demonstrably a single-step acquisition**, and no team's largest observed step carries
  more than 90% of its climb.

### 4.3 T3 — `Haystack`, the natural experiment

| date | # | model | `pds` | raw `mse` | cos floor | `reach` | `fid` | `jac` | overall |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 08-22 | 3 | `GSA-003` | 0.6861 | 5.8099 | vacuous | **0.2050** | 0.5065 | 0.0261 | 0.0923 |
| 08-25 | 9 | `CDX-J75A-G8` | 0.6934 | 5.9992 | vacuous | 0.2142 | 0.5065 | 0.0258 | 0.0932 |
| 08-26 | 10 | `GSA-009` | 0.7142 | 3.6025 | vacuous | 0.2366 | 0.5073 | 0.0254 | 0.1137 |
| 08-28 | 13 | `CDX-R28B` | 0.7810 | 2.4212 | vacuous | 0.2450 | 0.5074 | 0.0250 | 0.1472 |
| 08-30 | 16 | `GSA-015` | 0.7875 | 0.8525 | **0.08246** | 0.2408 | 0.5053 | 0.0233 | 0.1671 |
| 09-01 | 19 | `GSA-018` | 0.8025 | 0.8301 | 0.09327 | 0.2488 | 0.5066 | 0.0230 | 0.1822 |
| 09-02 | 21 | `GSA-020` | 0.8014 | 0.8291 | 0.09379 | 0.2494 | 0.5067 | 0.0230 | 0.1824 |
| 09-03 | 23 | `GSA-023` | 0.7914 | 0.8561 | 0.08075 | 0.2442 | 0.5059 | 0.0229 | 0.1714 |

**Raw `reach` 0.2050 → 0.2442 (+0.0392) while raw `pds` went +0.1052 and raw `mse` fell 6.8×.**
Their third submission already read `reach` 0.2050 — just under the shrunk cohort's 0.2106,
and above everything we produced until `rc_hybrid` twelve days later. Their direction has been
**flat within ±0.007 cosine since 08-30**.

**This separates the two channels cleanly from the outside.** A team can have ordering as a
structural property from day one *and* climb direction steadily; the two do not trade. On the live
board's own conditional (`reach = +0.3465·pds − 0.0667`, residual sd 0.0387, n = 341 usable rows):

| team | `reach` residual z |
|---|---:|
| Aginglab.com | +1.00 |
| **Haystack** | **+0.95** |
| Slurm Shadi | +0.84 |
| **us (after `rc_callset`)** | **+0.70** |
| cqawesome | +0.56 |
| Hongjian Zhou | +0.55 |
| Garrrbage | −0.45 |

§8.1 recorded us at **−1.92σ** on the field's conditional before `rc_hybrid`. We are now **+0.70σ**,
above the leader. The ordering deficit is closed and the remaining headroom in that channel is
`Aginglab`'s +1.00σ and `Vu Mai`'s `fid` finding (§3.5), not another `reach` arm.

### 4.4 The verdict, and why it is the opposite of what the mandate hoped for

The mandate offered three branches. The data picks none of them cleanly and the composite is
sharper than any of them:

| branch | verdict |
|---|---|
| *"one or two large steps → acquisition or method switch, with a date"* | **NO** for cqawesome (5 resolved steps, largest 28%) and Hongjian Zhou. Unfalsifiable-but-unsupported elsewhere. |
| *"many small steps → tuning of something they already had"* | **YES, and the channel budget names the knob: 93% and 98% of the two leaders' clean climbs are DIRECTION, with ordering ≤ 12% and amplitude ≤ 3%.** They are not tuning α and they are not tuning the emitter. They are buying bulk cosine. |
| *"flat and high since the earliest resolvable arm → they started with something we do not have"* | **NO for the leaders. YES for `Aginglab.com` — and Aginglab is stuck at 0.0989, 1.25× ours, not 1.7×.** |

**The decisive fact is a negative one that no amount of private data explains: five days ago the
board's two best directions were BELOW ours.** Floor against floor, on the same instrument:
`cqawesome` **0.06427** (08-30) and `jw_13288` **0.06333** (08-29) against our own a15 floor of
**0.07221** — 11% and 12% *under* us. Four days later they are 0.13211 and 0.13117. A private atlas is a build-once asset. A team holding one
does not spend thirteen days climbing from cosine ≈ 0, and does not sit at *our* cosine on
2026-08-30 and at twice it on 2026-09-03.

**206's verdict (METHOD/PUBLIC-DATA, not private data) is upheld and its evidence is now direct.**
206 argued it from `k_eff` 10–14 exceeding the ~11 CRISPR-perturbed human lines that exist, and
from three mutually independent pipelines being up there. This pass adds the trajectory: **the
direction is not a stock anyone brought, it is a flow that two independent parties produced
inside a five-day window we can see, at a rate of **+0.068 cosine** — against a net **+0.0037**
across our own nine inverted board arms (197 §1: 0.07559 at a8 → 0.079249 at `rc_hybrid`), of which
the CD4 residual alone contributed **+0.00507**. They moved 18× what we have moved.**

---

## 5. Sensitivity — the one fitted constant

`E_DIR["pds"] = 0.6479` is the only number fitted anywhere in the instrument. `0.686` is the
board's own median over the 13 DIRECTION-labelled strangers' arms (ρ median 0.944), i.e. an
independent estimate that lands 5.9% above ours.

| team | ×0.70 | **base 0.6479** | board median 0.686 | ×1.30 |
|---|---:|---:|---:|---:|
| cqawesome | +0.1801 (113%) | **+0.1477 (93%)** | +0.1435 (90%) | +0.1302 (82%) |
| jw_13288 | +0.0771 (128%) | **+0.0594 (98%)** | +0.0571 (95%) | +0.0499 (83%) |
| Hongjian Zhou | +0.0965 (89%) | **+0.0814 (75%)** | +0.0795 (73%) | +0.0733 (68%) |
| Haystack | +0.0826 (105%) | **+0.0695 (89%)** | +0.0678 (86%) | +0.0624 (80%) |
| Garrrbage | +0.0875 (79%) | **+0.0740 (67%)** | +0.0723 (65%) | +0.0668 (60%) |
| Slurm Shadi | +0.0628 (77%) | **+0.0532 (66%)** | +0.0520 (64%) | +0.0480 (59%) |
| MT Zeng | +0.0533 (79%) | **+0.0386 (57%)** | +0.0367 (54%) | +0.0307 (45%) |
| Andrew Liu | +0.0257 (49%) | **+0.0213 (40%)** | +0.0208 (39%) | +0.0190 (36%) |
| **Aginglab.com** | **−0.0009 (−2%)** | **−0.0012 (−2%)** | **−0.0012 (−2%)** | **−0.0013 (−2%)** |
| us | +0.0822 (76%) | **+0.0695 (64%)** | +0.0679 (62%) | +0.0627 (58%) |

**No ordering of teams changes and no verdict changes over a ±30% sweep.** `Aginglab.com` is
−2% under every setting. The two headline windows in §3.2 do not use the constant at all — they
are differences of two 197 floors, with `pds` used only as a *sign* corroboration.

---

## 6. What this corrects and extends in 206

**Corrects.**

1. **206 §8 item 1 — "The board records that they had it by 09-01 and nothing about where it came
   from."** Half wrong. It is right for `Hongjian Zhou`, `Andrew Liu` and `MT Zeng`. It is wrong
   for the cluster, because `jw_13288` shrank on **08-29** and its rows carry the cluster's
   direction from cos floor **0.06333 → 0.13117**, ρ = 1.13, amplitude excluded at 14.9× the
   conservative leak bound. The cluster's ablation account is the intelligence leak.
2. **206 §1's "teams with raw `mse` < 1.0 on 2026-09-03T0030Z = 44."** The snapshot gives **41**,
   which is 197 §2.5's figure. The informative count (54 of 612) reproduces exactly.
3. **206 §3.3's "the component that moved is the amplitude, and only the amplitude"** is confirmed
   for the two joins it examined, and generalised: over the 34 both-endpoints-informative in-band
   pairs, only **24% of the field's `mse` motion is corroborated by `pds`** (§1.4's regression
   slope 0.240).

**Extends.**

4. **The channel budget did not exist before this pass** and is the deliverable: §2.
5. **`Aginglab.com`'s +0.0679 of cosine floor since 08-26 is 93% amplitude** — a case where 197's
   floor alone would have reported a large direction gain and `pds` refutes it. This is the
   sharpest available demonstration of why the two-member decomposition is necessary.
6. **The ordering channel is measured to be nearly unused by the field** (§3.5, 9 strict arms in
   707 pairs) and `Vu Mai`'s reverted A/B prices its `fid` leg at −0.0464 raw.

---

## 7. What I could NOT resolve, and why

1. **Whether `cqawesome`'s 08-30 → 09-02 doubling was one arm or five.** Five submissions
   (#19 → #24) elapsed across a 2.5-day gap in our snapshot coverage (08-30T15:55 → 09-02T06:30).
   The *whole-trajectory* shape says incremental (5 resolved steps overall, largest 28%), but the
   inside of that specific window is not observed. **A snapshot every 6 h from now on would close
   this**; the board API is free and unauthenticated.
2. **Anything about any leader before their first informative row.** `cqawesome`'s first 18
   submissions, `Hongjian Zhou`'s first 10, `Garrrbage`'s first 14 and `Haystack`'s first 15 all
   sit above raw `mse` 1.023611, where 197's floor is vacuous **and** the `pds` → cosine rate is
   uncalibrated (their `pds` is 0.49–0.79, outside the [0.79, 0.88] band). For that whole region
   I report raw `pds` moves and refuse to convert them to cosine. This is 197 falsifier 1 binding
   exactly as designed.
3. **Whether the 2026-08-28 Cell publication caused anything.** It is a coincidence of dates
   against a dataset that was already public. I did not and cannot establish causation from
   outside, and I did not run any query that could reveal what we are doing.
4. **The pds ↔ cosine map outside [0.79, 0.88].** It is certainly nonlinear (`pds` saturates) and
   there are **zero** both-endpoints-informative pairs below `pds` 0.79 anywhere on the board to
   calibrate it with. Every conversion in this document is inside the band or explicitly refused.
5. **Which knob inside the direction channel the leaders turned.** The budget says *direction, not
   ordering / amplitude / geometry*, at 93–98%. It cannot say *donor*, *estimator* or *weight*.
   The date coincidence and 206's four independent pointers name a candidate; nothing on the board
   confirms it.
6. **`Hongjian Zhou`'s and `Garrrbage`'s single-vs-many step question.** Both have unobserved gaps
   of 4–5 submissions around their largest steps. See the T2 caveat — "one step" is unfalsifiable
   at this cadence and I have not claimed it for anyone.
7. **`cqawesome`'s identity.** Untouched. This pass was about their trajectory, not about them,
   and 206 §5.3's negative stands.
8. **Per-context decomposition.** Still unavailable — the board reports the three-context average
   only, so 206 §8 item 9 is unchanged. If the leaders' edge is concentrated in context A (which
   is what H1 predicts) the board cannot show it.
9. **Whether the `Aginglab.com` `GeroAI_v9 → v11` `reach` jump (0.1672 → 0.2674) is an ordering
   arm.** It is not classifiable: `nmae` moved −0.0294 and raw `mse` −0.86 in the same arm, so it
   fails the strict pinned-profile test and the three channels are confounded in it.

---

## 8. Hand-offs

- **→ coordinator.** The mandate's "clean negative" did not happen. **The leaders' direction moved,
  recently, fast, and by 2× in five days, and the board measured it on two mutually independent
  pipelines.** The remaining 48 days should NOT go only to D/E/F robustness. 206's **H1** — the
  CD4 dose, shipped at effective fraction 0.0528 because 062 selected it on `pds`, the member 197
  proved is least sensitive to direction — is now backed by a trajectory as well as by arithmetic.
  210's `cd4f034` (built, +13.83% forward PLOLO cosine at f = 0.34 per 212 §4.3) is the arm.
- **→ whoever owns board intel next.** Pull the leaderboard **every 6 hours**, not every 12–24.
  Item 1 in §7 is unresolvable only because of snapshot cadence, and it is the single most
  informative unobserved interval on the board. The API is free.
- **→ whoever owns the ordering channel.** It is closed for `reach` (we are +0.70σ on the field's
  own conditional, above `cqawesome` and `Hongjian Zhou`) but `Vu Mai` measured **−0.0464 raw
  `fid`** for a rank-ordered emission at a bit-identical profile, i.e. `fid` is movable by ±0.046
  raw = ±0.026 overall inside the channel where `pds`/`mse`/`nmae` are pinned. 197 says we *lead*
  on `fid`; that is a member worth defending, not trading.
- **→ `AGENTS.md` §8.2.** Two updates. (a) *"the leading direction is ≥ 0.13117"* is current, but
  the sentence should record that **the leader's direction was 0.0643 on 2026-08-30** — the gap
  is four days old, not structural. (b) The `pds` invariance used throughout this pass
  (**exactly** ordering-invariant, amplitude-invariant to 0.0034·|Δ`mse`|) is a reusable
  instrument and belongs next to 197's floor.
- **→ §14 discipline.** Reconfirmed at cost to everyone but us. This experiment reconstructed the
  cosine trajectory of the board's two leading pipelines out of nothing but their published six
  metrics, and found the intelligence leak in a cluster whose four accounts publish **zero**
  descriptions — it was their *ablation account shrinking two days early*. Submission timing and
  order are themselves a public channel. Ours publishes nothing, and should keep not doing so.

---

## 9. Files

| file | what it does |
|---|---|
| `panel215.py` | builds the deduplicated 2,208-observation panel from 17 snapshots; 5 known-answer gates (scaling map vs `docs/BOARD-LEDGER.md`, our own two ordering points, 186's `pm-r14d` row, 197's cohort counts) → `work/panel215.json` |
| `chan215.py` | the instrument: channel response vectors, the adaptive `pds` threshold, the exact budget partition |
| `phase0215.py` | the five known-answer arms, the confusion matrix, the magnitude checks → `work/phase0215.{json,log}` |
| `phase1215.py` | every leader trajectory decomposed arm by arm; per-team channel budgets; direction-arrival profile → `work/phase1215.{json,log}` |
| `phase2215.py` | board-wide sweep (707 usable pairs, 280 teams); the out-of-sample ρ validation; measured / bounded / unreadable direction arms; the informative frontier → `work/phase2215.{json,log}` |
| `phase3215.py` | T1 (did anyone start high), T2 (concentration), T3 (Haystack + the `reach ~ pds` conditional) → `work/phase3215.{json,log}` |
| `valid215.py` | V2 (the amplitude leak measured on 51 strangers' arms) and V3 (amplitude exclusion for each leader) → `work/valid215.{json,log}` |
| `sens215.py` | the ±30% sweep of the one fitted constant, and `cqawesome` row by row → `work/sens215.{json,log}` |
| `work/lb-live.json` | the public leaderboard, fetched read-only 2026-09-04 (644 entries) |

**Inputs read, never written:** `data/leaderboard/*.json`, `docs/BOARD-LEDGER.md`, `docs/LOG.md`,
`experiments/186-intel-0903/board/trajectories.md`,
`experiments/185-rank-controlled-emission/BOARD-RESULT-rc_hybrid.md`,
`experiments/202-callset-rider/BOARD-RESULT.md`, `experiments/{197,206,212}/RESULTS.md`.
**Two network reads, both public and read-only:** the leaderboard API; Crossref
`10.1016/j.cell.2026.08.002` and the Hugging Face dataset API for `Xaira-Therapeutics/X-Atlas-Pisces`.
No production code was touched, no submission was built, `cell-eval2` was never run, `vcc` was
never run, and nothing was written outside `experiments/215-leader-decomposition/`.
