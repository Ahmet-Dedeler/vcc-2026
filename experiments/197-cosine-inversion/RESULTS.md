# 197 — Inverting the board for the leader's DIRECTION. The answer is **≥ 1.6×**, not 1.28×, and after `rc_hybrid` it is **110% of the remaining gap**.

**Date:** 2026-09-03 · algebra on 096's identity plus board arithmetic · **no `vcc`, no
`cell-eval2`, no build, no submission, nothing written outside
`experiments/197-cosine-inversion/`.** Peak RSS 1.47 GB (`gate197.py`, the emitter gates); the
rest is < 1.0 GB. One network read: the public leaderboard API, saved to `work/lb-live.json`.

> **`rc_hybrid`'s row landed mid-run and adjudicated this experiment's own pre-registration.**
> Predicted raw `mse` 0.8738050, board 0.8738045 — **error 4.9e−7** against a ±4e−3 band; predicted
> cosine 0.079248, inverted 0.079249 — **error 2.5e−7**. §5 falsifier 4. Every number below that
> refers to "our displayed entry" has been refreshed onto it (`rchybrid197.py`).

---

## Verdict in eight lines

1. **The inversion IS identifiable, but only from below, and only for 54 of 612 teams** (63 of 630
   on the live board). For a
   stranger the identity has three unknowns (`A`, `claim`, `X`) and one equation. Two of them
   close: the amplitude optimum sits exactly on the ρ cliff (096 (b)), and `claim ≤ K` because
   `PRED_TRACE_CAP_K = 1.0`. Both close in the *same* direction, so what comes out is a **floor**:
   `cos ≥ (1.023611 − raw) / 2.074495`, vacuous for any `raw ≥ 1.02361`.
2. **AGENTS.md §8's "at least 1.28×" is CONFIRMED as a construction and STALE as a number.**
   Re-derived here from scratch: Aginglab's raw 0.8270 → floor 0.09478 → 1.281× 096's 0.0740.
   The inequality points the way 096's late correction said. But Aginglab has not led since
   2026-08-28.
3. **The board's leading direction is ≥ 0.13117 against our measured 0.079249, i.e. ≥ 1.66×
   ours**; the overall leader is at ≥ 0.12821 = **1.62×**. Five teams sit within 4% of the top.
   Correcting for the bound's measured tightness on our own rows (×1.09) puts the true figure near
   **1.8×**. **The sentence in §8 should read "at least 1.6× ours".**
4. **Our own cosine is a measurement, not a bound: 0.079249** on the displayed entry
   (`rc_hybrid` = a15's profile), with the whole nine-arm trajectory spanning 0.0742–0.0792.
   096's 0.07556 is reproduced exactly on a8.
5. **The inversion passes its own consistency test, and now on the board.** Arms that changed
   only the emitter or only the amplitude invert to the same cosine to **0.03–0.13%**; the one
   apparent exception (abc-knee, −4.28% against a23) is the `C_p` proxy, and 187's board-calibrated
   `f = (1, 1, 0.9825)` collapses it to **−0.15%**. `rc_hybrid` closes it live: a byte-different,
   31%-sparser file with digest-identical panel totals inverts to **the same cosine to 2.5e−7**
   while moving `reach` by the largest single-member delta this project has recorded.
6. **One profile change in nine board arms measurably raised the direction: the CD4 residual,
   +5 to +7%** (a14 → a15). The tiered spread and the mass fix moved it by less than the
   instrument's own ±2% row-to-row systematic.
7. **The leader's cosine corresponds to ~10.9 effective independent donor lines (9.6–12.5) against
   our 3.09,** by 159's `1/cos²(k) = 66.5 + 1920/k`. **It is not reachable by donor count** — 166
   proved only four panel-covering lines exist, and the best equal consensus of all four is
   board-equivalent **0.0873**, i.e. **1.46× short** of the leader's floor. It is also **not
   reachable by donor count for them**: whatever they have is not a fifth cell line.
8. **What it is worth: +0.0189 to +0.0199 overall on `mse` alone** (exact through the identity),
   and **+0.03696 across `mse` + `pds` + `nmae` — 110% of the +0.03370 that is left of the gap
   after `rc_hybrid`.** We are now *ahead* of the leader on `fid` and `jac`. **Everything that
   remains between us and first place is direction.**

---

## 0. Known-answer gates — all eight PASS (`gate197.py`, `work/gate197.log`)

Rule 10's corollary. If any of these fails, the instrument is not the one 096/164/169/171/187
used and everything below is void.

| gate | expected | measured | source | verdict |
|---|---:|---:|---|---|
| `claim(spread 0.25, φ0)` | 0.11860 | **0.11860** (−0.00%) | 096 `claim_grid2.log` | **PASS** |
| `claim(spread 1.00, φ0)` | 0.19216 | **0.19438** (+1.16%) | 096 `claim_grid.log` | **PASS** (tol 3%) |
| `claim(spread 0.25, φ4)` = a8 | 0.16478 | **0.16652** (+1.06%) | 096 `board_solve` | **PASS** |
| `claim(a15 as shipped)` | 0.24581 | **0.24581** (+0.00%) | 164 §5 | **PASS** |
| a15 raw `mse` / `score_mse` | 0.88228 / 0.11112 | **0.88228 / 0.11112** | 164 §5, 169, 171 | **PASS** |
| `C_p`-free LS solve on a8+a10+a11 | `X` 0.019240, `ΣD` 0.229141 | **0.019240 / 0.229141** | 187 [KA-i] | **PASS** |
| — its three residuals | — | **+6.0e−5 / −6.1e−5 / +0.2e−5** raw | | |
| **out-of-sample**: a10 raw predicted before it was built | 1.3567 | **1.3567** (board 1.3524, err **+0.0043**) | 096 falsifier 1 | **PASS** |
| **out-of-sample**: a11 raw predicted before it was built | 0.9021 | **0.9021** (board 0.9100, err **−0.0079**) | 096 falsifier 2 | **PASS** |

The 3% tolerance on the `claim` gates is 096's own documented spread between its two grids on one
policy (0.116800 vs 0.11860). The two out-of-sample rows are the ones that make the decomposition
credible rather than merely self-consistent: both were written down before the arm existed and
both landed inside ±0.008 raw.

`PRED_TRACE_CAP_K` was read out of `cell_eval2.metrics.delta` rather than assumed: it is **1.0**,
so `claim ≤ Σ_p w_p C_p = K = 0.246529` exactly, which is the hinge of the whole field inversion.

---

## 1. Our own nine rows, inverted (`invert197.py`, `phase3197.py`)

The identity, per context (187's correction to 096, which averaged `min(claim, B)`):

```
MSE_c − 1 = ( A_c + Kc_c − min(claim_c, B_c) − 2 X ) / ΣD          raw = mean_c MSE_c
```

For our own arms `A_c`, `B_c`, `Kc_c` and `claim_c` are all measured on the shipped artifacts, so
one board equation gives `X` exactly, and

```
cos  =  X / ( sqrt(ΣD) · mean_c sqrt(A_c) )
```

is **exactly invariant to a global amplitude rescale** (X ∝ s, A ∝ s²), which is what makes the
trajectory a check rather than a fit.

| arm | α | board raw `mse` | `Ā` | `claim` | ρ_min | `X` | **cos** | what changed |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| a8 | 1.0796 | 1.362667 | 0.282948 | 0.16676 | 1.0000 | 0.019247 | **0.07559** | — (096's own 0.07556) |
| a10 | 1.0796 | 1.352415 | 0.282948 | 0.16914 | 1.0000 | 0.019233 | **0.07554** | emitter only (φ4 → s0.15/φ16) |
| a11 | 0.8380 | 0.909989 | 0.170798 | 0.16916 | 0.9775 | 0.014935 | **0.07550** | amplitude only |
| a12 | 1.0108 | 0.901349 | 0.248147 | 0.24592 | 0.9418 | 0.017693 | **0.07420** | tiered spread + α |
| a14 | 0.9820 | 0.882775 | 0.234208 | 0.24592 | 0.8889 | 0.017183 | **0.07418** | amplitude only |
| a15 | 0.9820 | 0.873805 | 0.233158 | 0.24615 | 0.8811 | 0.018317 | **0.07925** | **CD4 residual** |
| a23 | 0.9820 | 0.872353 | 0.231516 | 0.24594 | 0.8783 | 0.017930 | **0.07785** | mass fix (minnorm) |
| abc-knee | per-ctx | 0.874491 | 0.247189 | 0.24594 | 0.9997 | 0.017727 | **0.07452** | per-context α only |
| **rc_hybrid** | 0.9820 | 0.8738045 | 0.233158 | 0.24652 | 0.8811 | 0.018317 | **0.079249** | rank-controlled emitter, **panel totals identical to a15**, `claim` at 100.0% of the cap |

`a15`'s terms come from its **own delivered panel totals** — 196 streamed
`chal-totals-a15.npz` off `submissions/191-a15_local.h5ad`, which 191 §6 proved byte-identical to
the shipped a15 package — with `claim` from 191's own per-context measurement on the same file.
a23/abc-knee come from 180's `knee180.json`; a8/a10/a11/a12 from 096's `*-stats.npz`.

### 1.1 The consistency check — the mandate's central test

**Arms that changed only the emitter, or only the amplitude, must invert to the same cosine.**

| family | members | cos | relative spread | what varied |
|---|---|---|---:|---|
| v11 | a8, a10, a11 | 0.07559 / 0.07554 / 0.07550 | **0.13%** | emission (φ4→s0.15/φ16), then α 1.0796→0.838 |
| tiered | a12, a14 | 0.07420 / 0.07418 | **0.03%** | α 1.0108 → 0.982 |
| mass-fix | a23, abc-knee | 0.07785 / 0.07452 | **4.28%** | per-context α only |
| **a15 / rc_hybrid** | a15, rc_hybrid | 0.079248 / 0.079249 | **0.0003%** | the emitter only; per-gene panel totals digest-identical on 900/900 blocks |

The a15/rc_hybrid row is the mandate's decisive test and it was executed **on the board**, not on
our own terms: two byte-different files (19,049,370,200 vs 13,202,460,120 bytes, 31% sparser) whose
per-gene panel totals are digest-identical on 900/900 blocks scored `mse` and `pds` identical to
sixteen significant figures, so the inverted cosine agrees to **2.5e−7** — while `reach` moved raw
0.18117 → 0.22992, **+0.05494 scaled, the largest single-member move this project has recorded.**
A change that moves the ordering channel that hard moves the inverted direction by zero. That is
§8.1's two-channel split verified from the outside, and it is the cleanest evidence available that
the inversion reads direction and nothing else.

The v11, tiered and a15/rc_hybrid families all pass at a level far below anything the experiment
claims. **The mass-fix pair does not, and it is the most informative row in the table.** abc-knee differs from a23 in nothing but its
per-context amplitude, so a 4.3% cosine move is either a bug or an artefact of a term the
inversion holds fixed. It is the latter, and 187 had already named it:

| arm | hypothesis | `claim` A / B / C | branch A/B/C | cos |
|---|---|---:|---|---:|
| a23 | proxy `C_p` | 0.23375 / 0.23792 / 0.26613 | B / B / B | 0.07785 |
| abc-knee | proxy `C_p` | 0.23375 / 0.23795 / 0.26613 | B / B / B | 0.07452 |
| a23 | 187's `f = (1, 1, 0.9825)` | 0.23375 / 0.23792 / 0.26147 | B / B / B | **0.07785** |
| abc-knee | 187's `f = (1, 1, 0.9825)` | 0.23375 / 0.23795 / 0.26147 | B / B / **claim** | **0.07773** |

**−4.28% → −0.15%.** The reference's own jackknife on context C is 1.75% below the NTC-construct
proxy; abc-knee shipped C's α 0.9% above C's *true* cliff and paid for it, and the inversion —
which holds `claim` at the proxy — charges that overshoot to the cosine. This is the instrument's
own systematic and it is **one-signed**: a `C_p` proxy that is too high can only make our own
direction look worse than it is.

A second, independent read on the same systematic: a12's inverted cosine sits 1.7% below the v11
family's. That is exactly 187 [KA-ii]'s recorded a12 model residual (−0.0028 raw) converted into
cosine units (0.0028 × ΣD / 2 / √(Ā·ΣD) = 1.8%). **So the instrument's row-to-row systematic is
0.03–0.13% within a family and ≈ ±2% across families, rising to ±4% when an arm sits on a
mis-estimated cliff.** Every conclusion below is sized against that.

### 1.2 What actually moved our direction

| change | cos | Δ | verdict |
|---|---:|---:|---|
| tiered spread (a11 → a12) | 0.07550 → 0.07420 | −1.7% | inside the ±2% systematic — **not resolvable**, and 187 attributes it to a12's own model residual |
| **CD4 residual (a14 → a15)** | 0.07418 → **0.07925** | **+6.8%** | **the only resolvable improvement in eight board arms** |
| mass fix (a15 → a23) | 0.07925 → 0.07785 | −1.8% | inside the systematic — not resolvable |
| per-context α (a23 → abc-knee) | 0.07785 → 0.07452 | −4.3% | **artefact**, see §1.1 |

The CD4 number corroborates two independent earlier reads from different anchors: 170 gate 8
(a15's `X` 0.017766 against 0.01714 rescaled pre-CD4, **+3.7%**) and 187 (**+4.9%**). All three
are positive and of the same order; the spread between them is the anchor choice, not the sign.

**Our current aggregate transfer cosine is 0.0785 (a15 0.0792, a23 0.0779), against 096's
published 0.0740 / 0.07556.** The rise is the CD4 residual, and it is the only thing that has
raised it since August.

---

## 2. The field (`invert197.py` §2, `assump197.py`)

### 2.1 The inversion, and why it can only be a floor

At its own amplitude optimum a submission sits exactly on the ρ cliff (096 (b): below the cliff
`#348`'s budget refunds the whole differential amplitude, above it the least-squares slope takes
over, so the minimum is the kink), where

```
MSE_min = 1 + Kc/ΣD + q·ε − 2·cos·sqrt( q·(1+ε) ),      q = claim/ΣD,  ε = (A−B)/B ≥ 0
```

`ΣD`, `Kc` and the cap `K` are properties of the **withheld reference** and are therefore the same
for every team; `claim` is a property of the submission and is bounded above by `K`. `cos` is
**decreasing in q**, so the cap yields the **floor**, not a mid-range estimate — and a team that is
*off* its optimum scores worse than `MSE_min`, which also pushes its true cosine up. Both unknowns
therefore close in the safe direction:

```
cos  ≥  ( 1 + Kc/ΣD − raw ) / ( 2 sqrt(K/ΣD) )  =  ( 1.023611 − raw ) / 2.074495
```

with `Kc = 0.005410`, `ΣD = 0.229141`, `K = 0.246529`, and `ε = 0` taken as its own minimum
(a positive common mode raises the floor). The bound survives per-context optimisation exactly:
`mean_c √(K_c/ΣD) ≤ √(K/ΣD)` by Jensen.

**It is vacuous for any raw ≥ 1.02361.** That is falsifier 1, and it is not a technicality.

### 2.2 Tightness — calibrated on eight rows where we know the answer

| arm | raw | mean s/s_c | floor | **measured** | floor / measured |
|---|---:|---:|---:|---:|---:|
| a8 | 1.36267 | 1.297 | −0.16344 | 0.07559 | **VACUOUS** |
| a10 | 1.35241 | 1.288 | −0.15850 | 0.07554 | **VACUOUS** |
| a11 | 0.90999 | 1.000 | 0.05477 | 0.07550 | 0.726 |
| a12 | 0.90135 | 1.001 | 0.05894 | 0.07420 | 0.794 |
| a14 | 0.88277 | 0.972 | 0.06789 | 0.07418 | 0.915 |
| a15 | 0.87381 | 0.970 | 0.07221 | 0.07925 | 0.911 |
| a23 | 0.87235 | 0.969 | 0.07291 | 0.07785 | 0.937 |
| abc-knee | 0.87449 | 1.000 | 0.07188 | 0.07452 | 0.965 |

**Violated on 0 of 8.** For arms sitting on their cliff the floor recovers **91–97%** of the true
cosine; for a11, which sits 10% below its cliff in `claim` units, it recovers 73%.

**a8 and a10 are the decisive rows for falsifier 1.** Their measured cosine, 0.07559, is the
*highest* of any v11-family arm — and the field method reports "no information" for them. A raw
`mse` of 1.36 is **not** evidence of a poor direction; it is evidence of an un-shrunk amplitude.
Any fit that let a clamped or un-shrunk team enter at its raw value would invert the sign of the
inference. They are excluded here explicitly, not down-weighted.

### 2.3 Testing the second assumption rather than asserting it

**(a) Our own arms sweep the assumption across a 1.34× span in `s/s_c`** (a8 at 1.297 down to a15
at 0.970) and the loss is monotone in exactly the ρ-cliff shape: on the cliff, 3–9%; 10% below,
27%; 30% above, everything. That is a direct measurement of what the assumption costs, not a
plausibility argument.

**(b) The field is measured, not assumed, from its own submission history.** Walking all 15
snapshots on disk for the last step of each top team:

| team | Δ raw `mse` | Δ raw `pds` | reading |
|---|---:|---:|---|
| Hongjian Zhou | −0.0067 | +0.0006 | still descending `mse`, `pds` flat → **not at optimum** |
| cqawesome | −0.0028 | −0.0011 | still descending → **not at optimum** |
| Aginglab.com | −0.0085 | +0.0028 | still descending → **not at optimum** |
| Vivai | −0.0116 | +0.0001 | still descending → **not at optimum** |
| Andrew Liu | +0.0008 | +0.0001 | flat → consistent with being at its optimum |
| Ahmet Dedeler | +0.0021 | +0.0001 | flat → consistent (187: global α is closed at 0.982) |

Four of six leading teams are still moving `mse` downward with `pds` standing still — 096's own
signature of a team descending its amplitude curve. Those teams are **above** their optimum, so
their observed raw exceeds their own minimum and their true cosine is **higher** than the floor.
**The assumption fails in the safe direction on the teams where it fails at all.**

### 2.4 Sensitivity of the answer to everything it rests on

| perturbation | leader floor | × ours | change |
|---|---:|---:|---:|
| as computed (`ε = 0`) | 0.12749 | 1.623 | — |
| `ε` = our own 0.0052 | 0.12985 | 1.653 | +1.9% |
| `ΣD` +5% / −5% | 0.13008 / 0.12484 | 1.656 / 1.589 | ±2% |
| cap `K` +10% / −10% | 0.12155 / 0.13438 | 1.547 / 1.711 | ∓5% |
| `Kc` +20% | 0.12976 | 1.652 | +1.8% |

Every perturbation that could plausibly be wrong leaves the ratio above **1.55×**, and the one
that moves it most (`ΣD`) is the best-determined constant in the whole experiment — three of our
own board points fit it to ≤ 6e−5 raw.

### 2.5 The cohort

| cut | teams (612-team snapshot, 2026-09-03T0030Z) |
|---|---:|
| raw `mse` < 1.0 | 41 |
| raw `mse` < 1.02361 — **the informative cohort** | **54** |
| raw `mse` < 1.2 — the mandate's "shrunk cohort" | **106** |
| `pds` > 0.80 | 51 |
| `score_mse` > 0 | 45 |

*Note on the mandate's figure:* it says the shrunk cohort at raw `mse` < 1.2 is 49 teams. On this
snapshot it is 106. 49 is 186's `reach`-pair cohort; the nearest `mse`-based cuts are 41
(raw < 1.0), 45 (`score_mse` > 0) and 51 (`pds` > 0.80). The board has grown from 87 teams on
2026-08-21 to 612 here and 630 on the live read.

**The unshrunk 558 are bounded separately, and the bound is nothing.** With `A` free the identity
`2X = A + Kc − min(claim,B) − (raw−1)ΣD` admits `cos = X/√(A·ΣD)` arbitrarily close to zero for
large `A` at any observed raw. So a raw above 1.02361 constrains the direction not at all — which
is precisely what our own a8/a10 rows demonstrate from the inside. There is no useful upper bound
either (`cos ≤ 1` is never binding).

### 2.6 The board, inverted

Top of the informative cohort, 2026-09-03T0030Z:

| # | team | overall | raw `mse` | **cos floor** | × ours | raw `pds` | `nmae` | `q_w ≥` | `reach` |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | Andrew Liu | 0.20834 | 0.75914 | **0.12749** | 1.62 | 0.83286 | 0.92969 | 0.53515 | 0.2490 |
| 2 | Hongjian Zhou | 0.20776 | 0.75992 | 0.12711 | 1.62 | 0.83252 | 0.92908 | 0.53546 | 0.2454 |
| 3 | **cqawesome** (overall leader) | **0.21790** | 0.76044 | 0.12686 | 1.61 | 0.86839 | 0.92688 | 0.53656 | 0.2592 |
| 4 | Garrrbage | 0.19960 | 0.78578 | 0.11465 | 1.46 | 0.85143 | 0.97265 | 0.51368 | 0.2110 |
| 5 | Yuzhu Shi | 0.18882 | 0.82578 | 0.09536 | 1.21 | 0.81348 | 0.92152 | 0.53924 | 0.2061 |
| 7 | Aginglab.com | 0.21122 | 0.82698 | 0.09478 | 1.21 | 0.81994 | 0.88921 | 0.55540 | 0.2621 |
| — | **us (abc-knee)** | 0.17404 | 0.87449 | 0.07188 | — | 0.83679 | 0.94317 | 0.52841 | 0.1773 |

On the **live** board (630 teams, ~12 h later) the leading direction has moved to **0.13117**
(`jw_13288`) = **1.67× ours**, with five teams above 0.127.

---

## 3. Three routes to one latent quantity — and they disagree (`phase3197.py`)

| team | raw `mse` | **A**: cos floor | **B**: cos ← `pds` | **C**: `q_w ≥` | × ours (A) | × ours (B) | × ours (C−0.5) |
|---|---:|---:|---:|---:|---:|---:|---:|
| us (abc-knee) | 0.87449 | 0.07188 *(measured 0.0785)* | 0.10768 | 0.52841 | — | 1.000 | 1.000 |
| Andrew Liu | 0.75914 | 0.12749 | 0.10589 | 0.53515 | **1.62** | 0.98 | 1.24 |
| Hongjian Zhou | 0.75992 | 0.12711 | 0.10574 | 0.53546 | **1.62** | 0.98 | 1.25 |
| cqawesome | 0.76044 | 0.12686 | 0.12204 | 0.53656 | **1.61** | 1.13 | 1.29 |
| Garrrbage | 0.78578 | 0.11465 | 0.11433 | 0.51368 | 1.46 | 1.06 | 0.48 |
| Aginglab.com | 0.82698 | 0.09478 | 0.10002 | 0.55540 | 1.21 | 0.93 | 1.95 |

- **Route A** — the `mse` identity, §2. Says **≥ 1.6×**.
- **Route B** — 159 §5.4's synthetic cos → `pds` curve, inverted on the HCT116 surface. Says
  **≈ 1.0–1.1×**.
- **Route C** — 187 §5's exact identity `nmae = 1 − 2 q_w W_hi − (2q_lo − 1) K_lo + K_hi`, read at
  its `K ≤ 1` bound so that `q_w ≥ (1 + (1 − nmae))/2`. Says **≈ 1.25–1.3×** on sign-accuracy
  excess above chance.

**Do not average these. The disagreement is the finding, and it has a mechanism.**

`pds` is a **rank** statistic and our whitening buys retrieval without buying bulk cosine — AGENTS.md
§4.2 measured exactly that (macro cosine 0.0435 → 0.0295, **−32%**, while macro PDS rose +9.5%),
and 159 §5.5 measured the same thing from the other side (at cosine 0.038 the isotropic curve gives
pds 0.637; production gets 0.826). So Route B is not measuring the same object as Route A: it is
measuring direction *after* our geometry has already converted it. **`pds` is the one member where
we are not direction-limited**, and it is the one route that says the leader is close.

Route C sits between them, and its residual against Route A is itself informative. Through the
Gaussian sign law `2q − 1 = (2/π)·arcsin(cos)`:

| team | cos floor | ⇒ `q_w ≥` implied | `nmae`'s own floor | reading |
|---|---:|---:|---:|---|
| us | 0.07188 | 0.52290 | 0.52841 | consistent |
| Andrew Liu | 0.12749 | 0.54069 | 0.53515 | `nmae` **short** of what `mse` implies |
| cqawesome | 0.12686 | 0.54049 | 0.53656 | `nmae` **short** |
| Aginglab.com | 0.09478 | 0.53022 | 0.55540 | `nmae` **ahead** of what `mse` implies |

The three top-`mse` teams' `nmae` is *behind* their own direction, and Aginglab's is *ahead* of
theirs. That reproduces 187 §5 independently and from a different statistic: **the leader has no
`nmae`-specific mechanism; Aginglab does.**

### 3.0 The three routes over all 54 teams, not just the top (`routes197.py`)

| pair | Spearman ρ | p |
|---|---:|---:|
| A (`mse` floor) vs B (`pds` curve) | **+0.810** | 1.2e−13 |
| A (`mse` floor) vs C (`nmae` `q_w`) | **+0.658** | 6.4e−08 |
| B (`pds` curve) vs C (`nmae` `q_w`) | +0.720 | 8.3e−10 |
| A (`mse` floor) vs raw `reach` | +0.701 | 3.6e−09 |

**The three routes order the field the same way. Their disagreement is entirely a disagreement
about SCALE**, which is what makes it interpretable rather than a sign that one of them is broken.

The per-team ratio `B/A` — retrieval extracted per unit of bulk direction — is the whitening
effect measured on other people's submissions. Restricted to the 32 teams whose floor exceeds
0.05 (it explodes as `A → 0`):

| | `B/A` |
|---|---:|
| **us (abc-knee)** | **1.50** — rank 6 of 32 |
| cohort median | 1.19 |

**We extract ~26% more `pds` per unit of bulk cosine than the median shrunk team.** That is
AGENTS.md §4.2's whitening effect measured against the field for the first time rather than
against our own LODO, and it is why `pds` cannot be used to price the direction gap: **the member
where we are strongest is the member least sensitive to what we are short of.**

### 3.1 A correction to 159 §5.4's known-answer row

159 validated its cos → `pds` curve by writing: *"our board-implied aggregate cosine is 0.0722 …
the curve at cos = 0.072 predicts raw pds 0.743; the board reads our raw pds at ≈0.73–0.74."*

**The board has not read our raw `pds` at 0.73 since 2026-08-21.** `0.72793` is the **v3** entry
(raw `mse` 1.683). Every entry from 08-22 onward reads 0.816–0.837, and the a15-era rows that
159's cosine came from read **0.83599**. The known-answer row compared a September cosine against
an August `pds`. Re-run honestly, the curve **under**-reads our `pds` by ~0.09 at our cosine —
which is not a failure of the curve but the whitening effect above, and it is why Route B and
Route A part company. Recorded because 159 §5.4's "a synthetic construction reproduces a live board
number to 0.01" is currently cited as that curve's licence.

---

## 4. The target (`phase3197.py` §3, `target197.py`)

### 4.1 In effective donor lines

159 §5.5's law, `1/cos²(k) = 66.5 + 1920/k`, fitted on the HCT116 challenge-panel surface with
R² = 0.99998. The board ↔ surface conversion is the ratio our **own single prediction** spans
between the two, in which the surface factor cancels:

| board↔surface anchor | leader cos (surface) | **k_eff** | ceiling (board) | leader % of ceiling |
|---|---:|---:|---:|---:|
| a15 (×2.078) | 0.06136 | **9.6** | 0.25480 | 50.0% |
| a23 (×2.041) | 0.06246 | **10.1** | 0.25031 | 50.9% |
| abc-knee (×1.954) | 0.06525 | **11.4** | 0.23959 | 53.2% |
| 159's own ×1.8931 | 0.06734 | 12.5 | 0.23215 | 54.9% |

**The leader's direction is worth ~10.9 effective independent donor lines (9.6–12.5) against our
3.09.** 159 estimated 10.3 for the then-leader from a stale anchor; the number survives its own
re-derivation.

### 4.2 Is that reachable by donor count? No — for us or for them

| quantity | surface cos | board-equivalent | k_eff |
|---|---:|---:|---:|
| production (us) | 0.03814 | 0.0779–0.0792 | 3.09 |
| best equal consensus of **all four** panel-covering lines | 0.04278 | 0.0873 | 4.00 |
| **leader's floor** | 0.0613–0.0673 | **0.1275** | **9.6–12.5** |
| context-independent `k → ∞` ceiling | 0.12263 | 0.2350–0.2548 | ∞ |

166's census, over every donor table on disk: `orion_hek293t` 294/300, `orion_hct116` 290,
`kolf2` 273, `k562_gwps` 269, `h1_2025` 20, everything else **0**. Four lines, and 166 searched
Hugging Face, GEO (258 series), Europe PMC (229 papers), BioStudies/ArrayExpress, the Arc Virtual
Cell Atlas and IGVF (6,902 measurement sets) and found no fifth, with a structural reason: panel
coverage is a library-size gate with nothing in between, and common-essential libraries read
exactly 0/300 because the panel excludes common-essentials by construction.

**So the leader's floor is 1.46× the best consensus the four donor lines in existence can produce.
No re-weighting of what we hold reaches it, and no fifth line is available to them either.**

### 4.3 What else the arithmetic permits — and what it rules out

The leader sits at **~51% of the context-independent `k → ∞` ceiling**, so — and this is the load-
bearing negative — **their direction does not require a context-adapted model.** It is inside what
a context-independent predictor can in principle do. That rules out the most expensive
explanation.

What the arithmetic **permits**:

1. **A better estimator of the shared component from the same four lines.** The `1/k` law's `k_eff`
   is a statement about *equal-weight consensus*; it is not a ceiling on estimators. 159 measured
   that production already converts cosine into `pds` far better than the consensus does
   (pds 0.826 at cos 0.038 vs 0.721 at cos 0.043) — the converse trade is available too. This is
   exactly what 193 (whitener/eigenvalue map), 194 (donor re-estimation) and 195 (retrieval
   transform) own.
2. **Non-donor information that is not a cell line at all**, and therefore does not enter `k`:
   pretrained expression models (`STATE` and `Lingshu-Cell` appear verbatim in board model strings
   at ranks 11–13), sequence/graph priors, DepMap essentiality as an effect-size prior. 096 already
   measured STRING `kg_density` correlating with our per-perturbation accuracy at ρ ≈ +0.19 on both
   surfaces — and that it is at least partly a magnitude proxy.
3. **A cross-line transfer operator learned on a non-CRISPRi substrate.** 184 proved the operator
   is real, large and steeply `k`-dependent on Tahoe-100M (cos +0.4635 at k=49 against a null of
   +0.0389) and **absent** in CRISPRi donors (arm ≡ control ≡ null to three decimals). A team with
   Tahoe or a comparable multi-line drug atlas has a channel we measured and closed on our data.
4. **Amplitude they have not spent.** Four of the six top teams are still descending their own
   `mse` curve (§2.3b), so part of their remaining headroom is calibration we have already banked.

What it **rules out**:

- **More donor lines.** No fifth panel-covering CRISPRi screen exists (166).
- **Line adaptation from CRISPRi donors.** Closed on 50 lines by 184, in the strong form.
- **Global or per-context amplitude.** Closed by 187: a15/a23 sit on the joint `mse`+`nmae`
  optimum to four decimals, and the entire per-context pot is +0.0011.
- **`pds` as the instrument for any of this.** §3 above: it is a rank statistic and it is the one
  member that does not see the gap.

### 4.4 What closing it is worth

`mse` alone, computed exactly through the identity at `rc_hybrid`'s own `claim` (100.0% of the
cap) and our own optimal α, against its board `score_mse` of 0.1200566:

| cos | raw `mse` | `score_mse` | Δ overall vs `rc_hybrid` |
|---:|---:|---:|---:|
| **0.079249 (us, measured)** | 0.86830 | 0.12583 | +0.0010 |
| 0.09887 (Aginglab's floor) | 0.82740 | 0.16888 | +0.0083 |
| **0.12821 (overall leader's floor)** | **0.76627** | **0.23322** | **+0.0189** |
| **0.13117 (best direction on the board)** | **0.76009** | **0.23972** | **+0.0199** |
| 0.14298 (best floor × 1.09, the tightness correction) | 0.73549 | 0.26561 | +0.0243 |
| 0.25032 (`k → ∞` ceiling, board scale) | 0.51219 | 0.50062 | +0.0634 |

And the whole remaining gap, decomposed exactly from the board's own scaled members, on the
**live** board with `rc_hybrid` displayed:

| member | us (`rc_hybrid`) | cqawesome | Δ/6 | share of the +0.03370 gap |
|---|---:|---:|---:|---:|
| `mse` | 0.12006 | 0.24229 | **+0.02037** | **60.5%** |
| `pds` | 0.74587 | 0.81459 | +0.01145 | 34.0% |
| `nmae` | 0.09429 | 0.12509 | +0.00513 | 15.2% |
| `reach` | 0.17031 | 0.20106 | +0.00513 | 15.2% |
| `fid` | −0.02976 | −0.07137 | **−0.00693** | **−20.6%** |
| `jac` | −0.00025 | −0.00896 | −0.00145 | −4.3% |
| **total** | **0.18342** | 0.21712 | **+0.03370** | 100% |

**`mse` + `pds` + `nmae` — §8.1's pure-direction members — are +0.03696 against a gap of
+0.03370, i.e. 110% of what is left.** `reach` adds +0.00513 more (rc_hybrid took two thirds of
the `reach` gap that stood twelve hours earlier), and `fid` and `jac` run the *other* way: **we are
ahead of the leader on both, by 0.00838 combined.**

> **Bracket for closing the direction gap: [+0.0189, +0.0370] overall.**
> Lower bound exact through the identity on `mse` alone; upper bound the three pure-direction
> members. **The upper bound exceeds the entire remaining gap to first place.**

This is the single most consequential number in the experiment, and it is a change of régime:
before `rc_hybrid`, direction was 81% of a +0.0439 gap and ordering was the rest. After it,
**direction is the whole thing.** §8.1's allocation frame stands, but its two channels are no
longer both live: the ordering channel has been spent down to +0.005 and `fid`/`jac` are already
in our favour.

The board's own conditional means, over the informative cohort at 00:30Z, corroborate it:

| quartile of cos floor | n | cos | overall | raw `pds` | `nmae` | `reach` | `fid` |
|---|---:|---:|---:|---:|---:|---:|---:|
| Q1 | 14 | 0.10031 | 0.19135 | 0.82536 | 0.93144 | 0.2208 | 0.5052 |
| Q2 | 14 | 0.07315 | 0.15792 | 0.80751 | 0.94031 | 0.2029 | 0.4864 |
| Q3 | 13 | 0.04161 | 0.13872 | 0.78949 | 0.95182 | 0.1706 | 0.4986 |
| Q4 | 13 | 0.00484 | −0.21433 | 0.61714 | 1.27781 | 0.1026 | 0.1876 |
| us, **abc-knee** (as displayed at 00:30Z) | | 0.07855 | 0.17404 | 0.83679 | 0.94317 | 0.1773 | 0.5038 |
| us, **rc_hybrid** (as displayed now) | | 0.07925 | 0.18342 | 0.83599 | 0.94422 | **0.2299** | 0.5037 |

At 00:30Z our `pds` beat even Q1's mean while our `reach` sat below Q2's at a cosine between Q1's
and Q2's — §8.1's ordering deficit, seen from outside our own instruments for the first time.
**`rc_hybrid` closed it: our `reach` is now above Q1's mean while our cosine is still between Q1
and Q2.** The only member left where we sit below our own cosine's cohort is `mse`, and that is by
construction, because `mse` *is* the cosine. (Q4 is contaminated — mean `nmae` 1.278, i.e. broken
submissions — so the usable slope is Q1 vs Q2: +0.0123 overall per +0.01 cosine, cross-sectional
and confounded.)

---

## 5. Falsifiers

**1. A clamped or un-shrunk team carries no information — excluded explicitly.** 558 of 612 teams (567 of 630 on the live board)
have raw `mse` ≥ 1.02361 and enter nothing. Our own a8 and a10 are the proof that this is a real
exclusion and not a formality: their measured cosine, **0.07559**, is the highest of any
v11-family arm, and the field method reports nothing for them. A high raw `mse` means an
un-optimised amplitude, not a poor direction.

**2. The byte-identical resubmission (`5QfX14axLAG8FF2rSDJe`) inverts to the identical cosine.**
`docs/LOG.md` 2026-08-31 records `score_mse` `0.1200566151966586` against a15's
`0.1200566151966586` — identical to sixteen digits. Inverting the recovered scorer map gives raw
`0.873790` against a15's recorded `0.873805` (gap 1.5e−5, the affine map's own 2.1e−5 fit
residual, not scoring noise). Identical `score_mse` ⇒ identical per-context raw ⇒ identical
`A`, `B`, `claim`, `Kc` ⇒ identical `X` ⇒ **cos = 0.079248, both times, exactly.** The falsifier
also certifies that nothing in the inversion reads a quantity that varies between two runs of the
scorer.

**3. AGENTS.md §8's `1.28×` — adjudicated.**

> **CONSTRUCTION: CONFIRMED.** Independently re-derived here. The `#247` cap bounds `claim` from
> above and therefore bounds the implied cosine from **below**; the inequality points the way
> 096's late correction said it does, and the earlier "0.97×–1.28×, plausibly zero gap" reading
> had it backwards. The floor is now additionally validated on eight of our own board rows where
> the cosine is *measured* rather than bounded: valid on 8/8, and 0.91–0.97× the truth on the four
> rows near the cap.
>
> **NUMBER: STALE, and stale in the direction that matters.** 096 computed it against
> Aginglab.com at raw 0.8270 and our then-cosine 0.0740. Reproduced exactly: floor
> `(1.023611 − 0.8270)/2.074495 = 0.09478` (096: 0.09469), ratio `1.281×` (096: 1.28×).
> Aginglab has not led since 2026-08-28 and is now 7th by implied direction. **The board's
> leading direction is ≥ 0.13117 against our measured 0.079249.**
>
> **The sentence in §8 should read "at least 1.6× ours" — 1.62× against the overall leader,
> 1.66× against the board's best direction, and ~1.8× once the bound's measured tightness on our
> own rows (×1.09) is applied.** The strategic reading it supports is unchanged in kind and larger
> in degree: calibration is closed (187), and first place needs the direction.

**4. `rc_hybrid` — pre-registered blind, and CONFIRMED on the board two hours later.**
191 falsifier 1 proved its per-gene panel-total digests are identical to a15's on 900/900 blocks,
so its `A`, `B` and profile direction had to be a15's exactly; only `claim` moved, and only in
context B, where `min(claim, B)` was already on the `B` leg. Written down in `assump197.py` before
the row existed:

> **PREDICTION: rc_hybrid's board raw `mse` = 0.873805, and it inverts to cos = 0.07925.
> FALSIFIED** if the raw differs by more than ±0.004 (096's own out-of-sample band on a10/a11),
> or if the cosine lands outside [0.07608, 0.08242].

| quantity | pre-registered | board (`nQUJZhZb3MW2uFJHBX5e`, 11:18 UTC) | error | band | verdict |
|---|---:|---:|---:|---:|---|
| raw `mse` | 0.8738050 | **0.8738045054832962** | **−4.9e−7** | ±4e−3 | **CONFIRMED** |
| inverted cos | 0.079248 | **0.079249** | **+2.5e−7** | ±0.0032 | **CONFIRMED** |
| `score_mse` | a15's 0.1200566151966586 | **0.1200566151966586** | 0 | — | identical, 16 digits |
| `pds_cosine` | a15's 0.835987 | **0.8359866220735785** | 0 | — | identical |

Two orders of magnitude inside the band on both, from a prediction with no free parameter. And the
arm that produced it moved `reach` raw 0.18117 → **0.22992**, +0.05494 scaled — **the largest
single-member move in the project's history** — while moving the inverted direction by nothing.
That is the strongest single piece of evidence in this experiment that the inversion reads
direction and only direction, and it is the mandate's digest-identical consistency test settled on
the board rather than on our own terms.

---

## 6. What this closes, corrects, and hands on

**Closes.**
- The inversion is identifiable — as a floor, for 54 of 612 teams (63 of 630 live), with a measured
  tightness of 0.91–0.97 for teams near the cap. It is not "not identifiable"; it is one-sided.
- Our own cosine is now a measurement on nine rows with a 0.03–0.13% within-family consistency
  check and a live board confirmation to 2.5e−7, not a single solve. **0.079249 today.**
- The leader's direction is bounded from below at 1.6× ours (≈1.8× once the bound's own tightness
  is applied) and is worth **+0.0189 to +0.0370** overall. Every sensitivity leaves the ratio
  above 1.55×.
- Reaching it by donor count is closed: it needs ~11 effective lines, four exist, and 166 proved
  no fifth is available to anyone.
- **After `rc_hybrid`, direction is 110% of the remaining gap to first place.** The two-channel
  frame of §8.1 stays correct and stops being a *choice*: the ordering channel is spent down to
  +0.005 and `fid`/`jac` are already in our favour.

**Corrects.**
- **AGENTS.md §8** — "at least 1.28×" → **"at least 1.6×"**. Construction confirmed, number stale.
  This is the sentence that currently governs allocation.
- **AGENTS.md §8.1's allocation frame** — still correct in structure, but the ORDERING channel it
  lists as one of two live levers is now largely banked (`rc_hybrid`, +0.0549 scaled `reach`), and
  the residual `reach` gap is +0.0051 against a direction gap of +0.0370. Any further split of
  effort should reflect that.
- **159 §5.4** — its known-answer row compares a September cosine against an August `pds`
  (0.72793 is the v3 entry). The curve does not in fact reproduce our board `pds`; it under-reads
  it by ~0.09, and the gap is the whitening geometry, not curve error.
- **`scripts/analyse_submission.py`** — asked in the mandate. Its `cos = (r² + 1 − mse)/(2r)` gets
  the *shape* right (`mse` is quadratic in amplitude, so it inverts to a cosine) and that instinct
  is what this whole line of work rests on. It is wrong in three ways and 096 already flagged two:
  it omits `Kc` and `ρ·claim` entirely (on this panel `claim/ΣD` is 0.52–1.08, larger than the
  quantity being solved for); it uses the *plug-in* real amplitude rather than the corrected `ΣD`,
  which 096 measured attenuates by 39% under an 8× reference thinning while the deconfounded form
  does not move at all; and its `--r-assumed 0.92` default is a free parameter where the constrained
  optimum makes `r` determinate. Its LODO-slope half is unaffected and still useful. **Superseded
  by (★); do not use its cosine output.**

**Hands on.**
- **The instrument.** `cos ≥ (1.023611 − raw)/2.074495` is a one-line read on any board row, valid
  for raw < 1.02361 and tight to ~0.93 for a team near the cap. Use it on every new snapshot; it
  costs nothing and it is the only direct read on the field's direction we have.
- **The one thing that has ever raised our cosine is the CD4 residual, +5–7%.** Eight board arms,
  one resolvable improvement. Whatever 193/194/195 propose should be priced against that bar.
- **`fid` and `jac` run the *wrong* way at the top of the board** (−0.0059 and −0.0014 of the gap:
  we are *ahead* of the leader on both). Any arm that trades them for `reach` is trading a member
  we lead for one we trail, and this decomposition prices it.
- **For whoever re-runs this after `rc_hybrid` lands:** `uv run python gate197.py && uv run python
  invert197.py && uv run python phase3197.py && uv run python assump197.py && uv run python
  target197.py`. Total ~14 min, peak RSS 1.4 GB. Add rc_hybrid's row to `BOARD_RAW` in
  `invert197.py` and the §5.4 prediction adjudicates itself.

---

## 7. Files

| file | what it does |
|---|---|
| `gate197.py` | the eight Phase-0 known-answer gates, including the four emitter `claim` gates on the real challenge panel and 096's two out-of-sample predictions → `work/gate197.{json,log}` |
| `invert197.py` | our eight pre-`rc_hybrid` board rows inverted, the digest-identical consistency check, the `rc_hybrid` prediction, and the field floor over all 612 teams → `work/invert197.{json,log}` |
| `phase3197.py` | the exact per-context cosine, the `C_p`-proxy systematic, the three-route cross-check, the effective-donor translation, the §8 adjudication → `work/phase3197.{json,log}` |
| `assump197.py` | `s/s_c` per arm vs floor tightness, the sensitivity table, the field's own amplitude trajectories, falsifier 2 → `work/assump197.{json,log}` |
| `routes197.py` | the three-route rank agreement over all 54 informative teams and the per-team `B/A` retrieval-per-direction ratio → `work/routes197.{json,log}` |
| `target197.py` | `mse` at a given cosine, the binned board price of direction, the member-by-member gap decomposition, the live-board read → `work/target197.{json,log}` |
| `rchybrid197.py` | adjudicates this experiment's own `rc_hybrid` pre-registration against the landed board row, and refreshes the target onto it → `work/rchybrid197.{json,log}` |
| `work/lb-live.json` | the public leaderboard, fetched read-only 2026-09-03 (630 teams, `rc_hybrid` included) |

**Inputs read, never written:** `experiments/096-mse-floor/work/{a8,a10,a11,a12}-stats.npz`,
`ctrl-stats.npz`, `cp-proxy-ntc.npz`, `board_solve.npz`, `v11-stats.npz`;
`experiments/096-mse-floor/{msekit,decompose,mech_alpha_chal,mech_a12_audit}.py`;
`experiments/180-def-readiness/work/knee180.json`;
`experiments/191-rc-hybrid-build/work/claim191.json`;
`experiments/196-emitted-depth/work/chal-totals-a15.npz`; `data/leaderboard/*.json`;
`src/vcc2026/donor_transfer.py` (the production emitter, called not reimplemented);
`cell_eval2` primitives. No production code was touched and nothing was built or submitted.
