# 224 — Can ANY trained model beat the donor mean on held-line transfer? **Yes — the smallest one, and only that one.**

**Date 2026-09-04.** Deterministic (seed `20260904`; `20260902` for the 159 reproduction). **No
`vcc`, no `cell-eval2`, no submission, no `submissions/` write, no file modified outside
`experiments/224-learned-transfer/`, no cloud instance launched, $0 spent.** Peak RSS **9.6 GiB**
on one early run (diagnosed and fixed to 4.2 GiB, §9); **8.2 GiB** across everything after that.

---

## 0. Verdict, in nine lines

1. **YES, and it is the cheapest rung on the ladder.** A **per-gene diagonal** — one free scalar per
   gene, fitted by inner leave-one-training-LINE-out on `(donor consensus → held-out training line)`
   pairs — beats the equal-weight donor mean by **+5.16% relative transfer cosine** on the primary
   7-line surface (95% CI `[+4.63, +5.72]`, **7 of 7 held lines positive**) and by
   **+8.9% at production's `k = 3`** on the **300-gene challenge panel**, at **Δ`pds` +0.0118**.
   The pre-registered gate was **+2.9%**.
2. **Every higher-capacity rung is worse.** Feature ridge replacing the consensus: **−38.6%**
   (132 reproduced). Feature ridge stacked on it: **+7.8% on held perturbations, −9.1% on held
   lines** (195's wall reproduced, quantitatively). A full 25M-parameter learned gene-space
   operator: **+11.5% cosine at −0.032 `pds`** — it buys cosine by inflating the decoy. An
   MLP/factorised-bilinear over (target embedding × recipient basal): **+4.4%, no better than the
   plain diagonal it contains**, with every recipient-conditioned parameter driven to zero.
3. **The slope in `k` is negative or flat on every rung that clears anything** — the diagonal peaks
   at `k = 3` (+7.80%) and decays to `k = 6` (+5.16%). **So the stop rule fired and no GPU was
   rented.** Capacity is measurably not the binding constraint: a 5,070-parameter diagonal beats a
   25M-parameter operator on the endpoint that matters.
4. **All four Phase-2 falsifiers pass.** Gene-shuffle null **−7.2%**; moment-matching control (same
   estimator without the target correspondence) **−3.5%**; target-label derangement **−25.7%**;
   wrong-recipient null on the only recipient-conditioned rung reads **better than the arm**
   (184's exact signature). Leak audit: 0 panel coordinates in any gene axis, 0 duplicate line
   groups, `h1_2025` in no surface.
5. **The mechanism is identified and it is already in this repo.** The fitted gain is 189's per-gene
   cross-line reliability prior, recovered by a completely different estimator:
   **r(gain, 189's `theta`) = +0.81 / +0.90** on the two surfaces, r(gain, my own cross-line
   reliability statistic) = **+0.94**, and the gain vectors fitted for **different held lines
   correlate at r = 0.95**. **189 measured this statistic in the ORDERING channel, where it bought
   nothing. Nobody had put it in the DIRECTION channel.**
6. **It is orthogonal to whitening.** At production's `gamma = 0.65` the gain adds **+9.63%** cosine
   and **+0.0124 `pds`** *on top of* the whitened consensus — more than it adds to the raw one.
7. **The panel surfaces P4/P5/P6 read ~0 and that was POWER, not substance.** 7,312 per-gene gains
   cannot be fitted from 73–120 panel targets. Fitting on **1,758 NON-panel targets** and evaluating
   on **232 panel targets of a held line** recovers **+8.3% ± 0.3 across five seeds**.
8. **004's negative does not replicate here and 132's does.** 004 fitted the same estimator for a
   **single-donor** predictor and applied it to one; fitted for the **consensus** and pooled over
   inner folds it is worth +5–9%. 004's own variant, re-run, reads **+3.92%** — real but weaker.
9. **This is an effect-level measurement, not a board price.** It answers exactly the question asked
   (held-LINE transfer cosine against the donor mean). It does **not** run through
   `profile_from_effect` + the scorer read-back, so it is not a rule-13 arm. §10 says what is.

---

## 1. Phase 0 — the substrate (`substrate224.py`, `surfaces224.py`, `panel224.py`)

**17 donor-shaped tables on disk** (`delta.npy` + `perts.txt` + `genes.txt`), with panel coverage:

| line | group | targets | genes | panel/300 | role |
|---|---|---:|---:|---:|---|
| `orion_hek293t` | HEK293T | 17,977 | 38,606 | 294 | production donor |
| `orion_hct116` | HCT116 | 17,767 | 38,606 | 290 | production donor |
| `cd4_stim48hr` / `cd4_rest` | CD4T | 284 / 283 | 10,113 | 284 / 283 | production residual source |
| `kolf2` | KOLF2 | 273 | 18,106 | 273 | 159 §4: anti-aligned, own blob |
| `k562_gwps` | K562 | 9,787 | 8,248 | 269 | production donor |
| `feng_ipsc` | FENG_IPSC | 182 | 18,533 | 182 | 223: fails its control at every dose |
| `dld1_lowmoi` | DLD1 | 4,715 | 18,533 | 67 | 223's live source |
| `telohaec` | TELOHAEC | 2,167 | 17,464 | 65 | 105 + 223 NO-GO |
| `gse295214_h1` | H1 | 2,354 | 18,533 | 64 | Pan 2025 hESC |
| `h1_2025` | H1 | 255 | 18,080 | 20 | **leaks h1-test 85/85** |
| `myllia_a549` | A549 | 175 | 18,533 | 14 | |
| `frangieh_melanoma` | MELANOMA | 235 | 23,712 | 5 | |
| `k562_essential` | K562 | 2,021 | 8,563 | 0 | second K562 library |
| `rpe1` / `nadig_jurkat` / `nadig_hepg2` | RPE1/JURKAT/HEPG2 | 2,204 / 2,270 / 2,160 | 8,749 / 8,881 / 9,623 | 0 | training-only |

**Pairwise shared-target counts** are in `work/substrate224.json`. The load-bearing ones: the six
Replogle/Orion/Nadig tables share **1,800–2,250** targets pairwise; `gse295214_h1` shares
1,297–2,225 with them; `dld1_lowmoi` 457–4,387; **`kolf2` shares 0 with `rpe1`, `nadig_jurkat`,
`nadig_hepg2` and `k562_essential`** (which is why it can be an evaluation line but not a fit line);
CD4 shares 0 with the same four.

**Surfaces.** One TABLE per LINE (`k562_essential` dropped as a second K562 library, kept only as the
reliability anchor `rho`). Every line mean-removed over its own targets first; own-target and all
300 panel coordinates dropped from the gene axis.

| surface | LINES | targets | genes | note |
|---|---:|---:|---:|---|
| **S5** | 5 | 1,592 | 5,002 | 159 §5.2's set |
| **S7** | **7** | 968 | 5,070 | **primary — most distinct lines** |
| S8 | 8 | 219 | 5,681 | + DLD-1 |
| P4 / P5 / P6 | 4 / 5 / 6 | 240 / 232 / 147 **panel** | 7,312 / 6,406 / 6,406 | 159 §5.3's set |
| **crosspanel** | 6 fit / 5 eval | **1,758 non-panel fit → 232 panel eval** | 5,304 | §5 |

### 1.1 Known-answer rows — all pass (`known224.py`, `work/known224.json`)

Reproduces 159 §5.2 from `data/donors/*/delta.npy` alone, on the six-line intersection
(**1,592 × 5,206** — 159 reported 1,592 × 5,206).

| quantity | 159 | 224 |
|---|---:|---:|
| consensus-of-`k` cos, `orion_hct116`, k = 1..5 | 0.1065 / 0.1382 / 0.1579 / 0.1718 / 0.1822 | **0.10647 / 0.13818 / 0.15794 / 0.17181 / 0.18223** |
| consensus-of-`k` cos, `k562_gwps`, k = 1..4 | 0.1132 / 0.1492 / 0.1728 / 0.1900 | **0.11319 / 0.14923 / 0.17277 / 0.18999** |
| consensus-of-`k` cos, `nadig_jurkat`, k = 1..5 | 0.1232 / 0.1602 / 0.1841 / 0.2013 / 0.2146 | **all five to 4 dp** |
| `1/cos^2 = A + B/k`, `cos_inf` HCT116 | 0.253 | **0.2515** (R² 0.999927) |
| `cos_inf` K562 | 0.303 | **0.2995** (R² 0.999926) |
| same-LINE control, K562-gwps ← K562-essential | pds 0.9468 / cos 0.2627 | **0.9468 / 0.26270** |
| cross-line comparator, K562-gwps ← HCT116 | pds 0.7839 / cos 0.0897 | **0.7839 / 0.08968** |
| `kappa` = `cos_inf`/√rho | 0.591 | **0.5843** |
| `tau` = `kappa²` | 0.349 | **0.3414** |

A fourth known-answer row fell out of §6: re-fitting the **whitening** family on this surface
reproduces AGENTS.md §4.2's signature exactly — macro cosine **−32.8%** (§4.2 measured −32%) with
`pds` **+0.0459**.

---

## 2. The gate and the stop rule (pre-registered, `PREREGISTRATION.md`, committed `2618c79` before any fit)

```
GATE      relative cosine gain on held LINES >= +2.9%   (AGENTS.md §8.4 corrected shape-only bar,
                                                         conservative end of +2.5-2.9%)
          AND a paired target bootstrap excluding 0
STOP      the decision variable is d(gain)/dk, not the level.
          slope <= 0 and gate not cleared -> STOP, report the negative, do NOT rent a GPU.
SEPARATE  held-perturbation and held-line are always reported apart (195 cleared on one and
          failed on the other; a pooled number hides the wall).
```

---

## 3. Phase 1 — the ladder on the primary surface S7 (`ladder224.py`, `work/analyse-S7.txt`)

7 LINES, 968 targets, 5,070 genes; `T_fit` 484 / `T_eval` 484; hyperparameters selected inside
`T_fit` on a holdout that is **both a different training line and a disjoint target block**.

**Everything at `k = 6` (all six other lines as donors), paired target bootstrap, 4,000 draws:**

| rung | params | rel. cosine % | 95% CI | Δ`pds` | held-PERT % | slope in `k` | verdict |
|---|---:|---:|---|---:|---:|---:|---|
| **R0** donor mean | 0 | — | — | — | — | — | baseline `cos` **0.20910** |
| **R1** feature ridge, **replaces** consensus | 50 → G | **−38.55** | `[−42.18, −34.78]` | −0.1313 | −16.00 | −1.70 | **132 reproduced** |
| **R2** feature ridge **stacked** | 50 → G | **−9.12** | `[−11.41, −6.78]` | −0.0418 | **+7.79** | −1.98 | **195's wall** |
| **R2b** feature ridge on the residual | 50 → G | 0.00 | — | 0 | 0 | 0 | `b = 0` selected |
| **R3** learned gene-space operator | ≤ 25.7 M | +11.46 | `[+9.13, +13.94]` | **−0.0323** | +38.56 | −0.83 | cosine up, **retrieval down** |
| **R3b** operator blended with R0 | ≤ 25.7 M | +15.97 | `[+14.74, +17.21]` | **−0.0088** | +30.90 | −1.61 | same, milder |
| **R4** **per-gene diagonal** | **5,070** | **+5.16** | `[+4.63, +5.72]` | **+0.0105** | +7.74 | −0.36 | **CLEARS** |
| R4 · gene-shuffle null | — | −5.36 | `[−5.46, −5.26]` | −0.0007 | — | — | null fails ✓ |
| R4 · moment-match control | — | −3.64 | `[−3.86, −3.43]` | −0.0047 | — | — | control fails ✓ |
| R4 · fitted from a SINGLE donor (004's) | 5,070 | +3.92 | `[+3.23, +4.69]` | +0.0149 | +7.30 | −0.65 | weaker |
| R4 · ORACLE, fitted on the held line | 5,070 | +12.22 | `[+11.17, +13.29]` | +0.0066 | — | −1.11 | **leaky ceiling** |
| **R4a** recipient-basal gain (**184's shape**) | 12 | +0.41 | `[+0.34, +0.48]` | −0.0002 | +0.44 | +0.03 | **184 reproduced** |
| **R4b** basal-similarity donor weights (010's) | 1 | +2.81 | `[+2.28, +3.36]` | −0.0051 | +2.97 | +0.62 | below gate, P(>gate) 0.376 |
| **R5** MLP / factorised bilinear, rank 2/8/32 | ~10⁴ | +4.36 / +4.43 / +3.79 | — | +0.0073 | — | — | **no better than R4** |

### 3.1 Cosine vs number of training lines — the decision variable

Relative cosine gain (%), S7, macro over held lines and donor subsets:

| rung | k=2 | k=3 | k=4 | k=5 | k=6 | slope |
|---|---:|---:|---:|---:|---:|---:|
| R0 donor mean, absolute `cos` | 0.14437 | 0.16535 | 0.18518 | 0.19858 | 0.20910 | — |
| **R4 per-gene diagonal** | **+6.05** | **+7.80** | +6.87 | +5.97 | **+5.16** | **−0.36** |
| R4 diagonal, Δ`pds` | +0.0187 | +0.0167 | +0.0129 | +0.0116 | +0.0105 | — |
| R3 operator | — | +13.87 | +14.14 | +13.11 | +11.46 | −0.83 |
| R1 feature replace | −31.65 | −34.48 | −35.88 | −37.70 | −38.55 | −1.70 |
| R2 feature stack (LINE) | −1.16 | −4.10 | −6.03 | −7.98 | −9.12 | −1.98 |
| R2 feature stack (**PERT**) | **+79.51** | +37.85 | +22.38 | +13.41 | **+7.79** | — |
| R4a basal gain | +0.18 | +0.51 | +0.45 | +0.39 | +0.41 | +0.03 |

**The diagonal's advantage peaks at `k = 3` and decays.** The *absolute* gain is flatter
(+0.00874 → +0.01290 → +0.01271 → +0.01185 → +0.01080); it is the baseline that rises. On S8
(8 lines, 219 targets) the relative slope is **+0.785**/line and the gain is +4.82% at `k = 7`, so
the sign of the slope is surface-dependent while its magnitude is small everywhere — **there is no
surface on which more lines make the learned model dramatically more useful, which is the negative
184 found for CRISPRi and the opposite of what it found on Tahoe's 50 lines.**

> **A methodological correction I have to record.** On **S5** alone (5 lines, `k` capped at 4) the
> diagonal's slope reads **+0.77**/line and looks like a clean positive. It is an artefact of a
> truncated `k` range: the turnover is at `k = 3–4` and S5 cannot see past it. **A slope in `k`
> measured over fewer than ~5 values of `k` is not evidence of anything.** This is why S7 is the
> primary surface and why the S5 number is not the headline.

### 3.2 The one bug that would have inverted this experiment

The first implementation selected the operator's `(rank, lambda)` by holding out one inner training
**line** while keeping the same targets. That is contaminated: with `N < G` a near-interpolating
ridge memorises the target, and the held fold's row for target `t` is the mean of `k−2` of the same
lines that produced the fit rows for the same `t`. It selected **full rank at `lambda = 0.01`** and
read **−0.010 on the held line against a +0.188 donor mean**, i.e. a total collapse that looked like
a real negative result. Selecting on **a different line AND a disjoint target block** fixes it and
the same configuration reads **+0.175**. Recorded in `core224.py`'s docstring.

---

## 4. Rung 5 in detail — capacity, and the recipient (`mlp224.py`)

`P_hat[t,:] = P[t,:] * ( g0 + Σ_r a_r(e_t) · c_r(x_g) ) + Σ_r u_r(e_t) v_r(x_g)`, with `g0` a free
per-gene gain (so the model **contains rung 4 exactly**), `e_t` the PCA-50 target embedding, and
`x_g = [b_A, b_D, b_A − b_D]` the per-gene recipient context. Adam on the macro-cosine loss,
early-stopped on the same honest line+target holdout.

| rank | k=2 | k=4 | k=6 |
|---:|---:|---:|---:|
| 2 | +5.10% | +4.08% | +4.36% |
| 8 | +3.87% | +4.85% | +4.43% |
| 32 | +4.42% | +4.55% | +3.79% |

**Rank does nothing** (±0.6% across a 16× capacity range) and the whole family lands **at or below
the plain diagonal it contains**. The conditioned heads converge to their zero initialisation:
in **30 of 84** configurations the arm and its wrong-recipient null are **bit-identical**.

---

## 5. The panel — power, not substance (`crosspanel224.py`)

Run directly, the panel surfaces read nothing:

| surface | panel targets | best rel. cosine | 95% CI | P(>gate) |
|---|---:|---:|---|---:|
| P4 (4 lines) | 240 | **−0.04%** (k=3) | `[−1.08, +1.08]` | 0.000 |
| P5 (5 lines) | 232 | **+1.77%** (k=4) | `[−0.14, +3.84]` | 0.130 |
| P6 (6 lines) | 147 | **+0.45%** (k=5) | `[−0.92, +1.58]` | 0.000 |

Two readings: **(a)** the panel deliberately excludes common-essential genes (099 §1.1) so the
reliability ranking does not hold there, or **(b)** 7,312 per-gene gains cannot be estimated from
73–120 rows. S8 settles half of it immediately — 219 **non-panel** targets, the same thinness, and
it clears at +4.82%. The direct test: **fit the gain on 1,758 NON-panel targets, evaluate on 232
PANEL targets of a held line** (fit ∩ eval = **0**, asserted).

| k | R0 `cos` | + gain | rel % (5 seeds) | Δ`pds` | gene-shuffle | moment-match | **label-derange** |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0.01736 | 0.01919 | **+10.64** (10.41–11.03) | +0.0100 | −7.2 | −3.8 | −25.8 |
| 2 | 0.02405 | 0.02635 | **+9.56** (9.29–9.88) | +0.0112 | −7.2 | −3.7 | −25.8 |
| **3** | 0.02900 | 0.03157 | **+8.86** (8.57–9.14) | **+0.0118** | −7.2 | −3.6 | −25.8 |
| 4 | 0.03301 | 0.03575 | **+8.30** (7.99–8.57) | +0.0118 | −7.2 | −3.5 | −25.7 |

Bootstrap at `k = 4`: **+8.57% `[+5.19, +12.36]`, P(> gate) = 0.999**; ≥ 0.998 on all five seeds.
Per held line at `k = 4`: `k562_gwps` **+8.16**, `orion_hct116` **+7.46**, `orion_hek293t` **+7.89**,
`kolf2` +28.23, `cd4_rest` +2.79 — **5 of 5 positive**, and **+7.94% excluding `kolf2`** (which 159
§4 shows is anti-aligned and carries its own blob, and whose R0 is 0.013 so its relative gain is
cheap). Note the **absolute** gain *rises* with `k` here (+0.00185 → +0.00274) while the relative
gain falls.

**The strictest possible fit pool.** `rpe1`, `nadig_jurkat`, `nadig_hepg2` carry **0 of the 300**
panel targets, so a gain fitted on them alone cannot have seen a panel response in any line at any
point. It reads **+3.26% `[+1.13, +5.58]`, P(>gate) 0.628** at `k = 4` — weaker (3 fit lines instead
of 5, and the two shallowest tables) but still positive with all three nulls failing. **Treat +3.3%
as the conservative floor and +8.3% as the estimate at the fit pool actually available.**

---

## 6. Does whitening already collect it? No (`whiten224.py`)

Production ships a **whitened** consensus. The shipped operator's gene axis is
`union(donor coverage) − panel` and is **not recorded in its sidecar json** (AGENTS.md §8.3 is
exactly about this), so it cannot be restricted to this gene axis; instead the **same family** was
re-fitted on this surface from the same kind of basis (donor effect rows, held line excluded), which
is what §4.2 and `build_union_whitener.py` prescribe. Evaluated on held-line **panel** targets:

| `gamma` | `cos` (k=3) | `cos` + gain | rel % | `pds` | `pds` + gain | Δ`pds` |
|---:|---:|---:|---:|---:|---:|---:|
| 0.00 (raw) | 0.02900 | 0.03130 | +7.93 | 0.6660 | 0.6768 | +0.0108 |
| 0.30 | 0.02404 | 0.02645 | **+10.02** | 0.7002 | 0.7156 | **+0.0154** |
| **0.65 (production)** | 0.01949 | 0.02137 | **+9.63** | 0.7119 | 0.7243 | **+0.0124** |
| 1.00 | 0.01614 | 0.01748 | +8.35 | 0.7037 | 0.7122 | +0.0085 |

**The gain is larger on top of whitening than without it**, and it is largest at `gamma = 0.3–0.65`
— exactly where production sits. Whitening's own signature reproduces §4.2 (`cos` −32.8%, `pds`
+0.0459), which is the known-answer row for this table.

---

## 7. Mechanism — it is 189's `theta`, in the other channel (`mechanism224.py`)

| question | measurement |
|---|---|
| Is the gain a property of the gene or of the fold? | Gain vectors fitted for **different held lines** correlate at **r = 0.9525** `[0.918, 0.969]` (S7; 0.9178 on S5) |
| What predicts it? | per-gene **cross-line reliability r = +0.938** (Spearman +0.965); `sd_line/sd_consensus` **−0.941**; basal expression +0.689; single-line sd **+0.063** |
| Is it 189's `theta`? | **r = +0.807** (Spearman +0.794) on 4,865 shared genes, S7; **+0.897 / +0.864** on S5 |
| How big is it? | gain p10 **0.129**, median **0.227**, p90 **0.342** — a **2.7× spread** across genes |

**189 built this exact statistic**, validated it at `z = +41.6` against its own permutation null on a
leave-one-LINE-out over five lines, and then used it as an **ordering key**, where it moved `reach`
by `+0.00408 ± 0.00424` and was closed. **AGENTS.md §8.1 rule 1 says the two channels have different
instruments; 189 measured `theta` in the ordering channel and nobody measured it in the direction
channel.** That is the whole of this experiment's novelty, and the +0.81 correlation is the evidence
that the two objects are the same.

**Why 004's negative does not replicate.** 004 fitted `w_j = Σ δ_target·δ_donor / Σ δ_donor²` — the
identical estimator — on a **single donor→target pair** (K562→RPE1) with no shrinkage, applied it to
a **single-donor** prediction, and evaluated on H1. Fitted instead for the **consensus** and pooled
over inner folds it reads **+5.16%**; 004's own single-donor variant, re-run inside this harness,
reads **+3.92%** — real but 24% weaker, **negative on `orion_hct116` (−1.24%) and ~0 on
`k562_gwps` (+0.53%)** where the consensus-fitted one is positive on 7 of 7. 004's conclusion ("do not spend more time on transforms
of the donor matrix without a mechanism that explains why these two failed") is satisfied: the
mechanism is the predictor mismatch plus the absence of shrinkage.

---

## 8. Phase 2 — the four mandatory falsifiers

**1. Wrong-recipient null (the one that killed 184).** R4a (basal-conditioned gain) and R5 (the MLP)
are the only recipient-conditioned rungs. R5 on S7, 84 configurations:

| | value |
|---|---:|
| arm (right recipient's basal) | **0.17535** |
| **wrong-recipient null** | **0.17634** |
| donor mean | 0.16783 |
| configurations where arm and null are **bit-identical** | **30 / 84** |

**The wrong recipient is better than the right one**, and in more than a third of configurations the
model's recipient-conditioned parameters are exactly at their zero initialisation. This is 184's
`+0.0665 / +0.0667 / +0.0663` in a different currency. R4a's own gain is **+0.41%** against a +2.9%
gate. **Recipient conditioning is worth nothing on CRISPRi — 184's finding, reproduced on 7 CRISPRi
lines after it was established on 50 Tahoe lines.**
The diagonal that *does* work is **recipient-independent by construction**, so it cannot be a
wrong-recipient artefact; its analogous test is the decoy/`pds` pair in §3, which it passes
(`decoy` moves +0.0003 while the same-target cosine moves +0.011, and `pds` **rises**).

**2. Permutation over target labels** (`--permute-labels`, whole re-runs):

| rung | honest | labels deranged |
|---|---:|---:|
| R4 per-gene diagonal, S7 | **+5.16%** | **−4.48%** |
| R4 per-gene diagonal, P5 | +1.77% | −27.50% |
| R4 diagonal, crosspanel | **+8.30%** | **−25.7%** |
| R3 operator | +11.46% | −85.38% |
| R1 feature replace | −38.55% | −96.21% |
| R4b basal donor weights | +2.81% | +2.81% (correct — it reads no labels) |

**3. Held-perturbation vs held-line, never pooled.** The wall is real and it is where 195 found it:

| rung | held PERTURBATION | held LINE |
|---|---:|---:|
| R2 feature stack | **+7.79%** | **−9.12%** |
| R2 feature stack at k=2 | **+79.51%** | −1.16% |
| R3 operator | +38.56% | +11.46% (at −0.032 `pds`) |
| **R4 per-gene diagonal** | **+7.74%** | **+5.16%** |

The diagonal is the only rung whose held-line number is close to its held-perturbation number. 195's
`+0.00904 across unseen perturbations / −0.00618 across an unseen line` is reproduced in kind by R2.

**4. Leak audit** (`leakaudit224.py`, `work/leakaudit224.json`) — computed, not asserted:

- **`h1_2025` covers 85/85 (100%) of h1-test's targets and appears in NO surface used here.** Next
  highest: `gse295214_h1` 30/85, `dld1_lowmoi` 28/85, `nadig_jurkat`/`nadig_hepg2` 21/85.
- **0 panel target coordinates in the gene axis of any of the six surfaces** (`pds_cosine` drops all
  300 — matched); **0 own-target coordinates** in any broad surface's gene axis.
- **0 duplicate LINE groups** in any surface: `k562_essential` is dropped as a second K562 library
  and `k562_gwps`/`k562_essential` never straddle a split; `h1_2025` and `gse295214_h1` are the same
  LINE and never co-occur.
- **Production inputs present**: S5/S7/S8 carry `k562_gwps` + both Orion lines; P5/P6 additionally
  carry `cd4_rest`, which is inside 222's shipped union. A fold holding `cd4_rest` out is therefore
  **not independent of the shipped artifact** and is reported separately (it is the weakest held
  line, +2.79%).
- **crosspanel fit ∩ eval targets = 0** (asserted at runtime); the strictest variant fits only on
  lines with **0/300** panel coverage.
- Never touched: `cell-eval2`, any scorer, the whitening operator on disk, `vcc`, `submissions/`, the
  CD4 residual sizing rule, or any file outside this directory.

---

## 9. Cost, and proof that nothing is left running

**No cloud instance was launched. Cost: $0.00.**

The mandate authorised a GPU **only if rungs 1–3 showed a positive slope in training lines**. The
measured slopes on the primary surface are **−1.70 (R1), −1.98 (R2), −0.83 (R3), −0.36 (R4)**; the
only positive slopes belong to rungs that never clear the gate (R4a +0.03, R4b +0.62) and to a null
(gene-shuffle +0.91). Independently, **capacity was measured directly and is not the constraint**: a
5,070-parameter diagonal beats a 25.7M-parameter linear operator on `pds` by **0.043** and matches or
beats a conditioned MLP at three ranks. Renting a GPU would have bought a more expensive confirmation
of a measured negative. **The stop rule fired and I stopped.**

`scripts/aws_reap.sh` (221's outer guard) swept **every enabled region** before and after:

| | running instances | flagged | names |
|---|---:|---:|---|
| **before** (`work/reap224-before.json`) | 4 | 0 | `vcc-cloud-score-ka2-…`, `vcc-cloud-score-ladder{6,7,8}-…` |
| **after** (`work/reap224-after.json`) | **3** | 0 | `vcc-cloud-score-ladder{6,7,8}-…` |

All are **experiment 221's** cloud scorers (its launcher's naming), all young, all at 1.0–1.1 busy
cores i.e. actively scoring, and the count **went down by one** during this experiment because one of
221's workers hit its own `shutdown -h` guard. **224 created zero instances, zero buckets, zero S3
objects.** No `224-*` tagged resource exists in any region.

**Local resources.** Peak RSS **9.6 GiB** on the first `crosspanel224.py` run — over comfort, at the
edge of the 10 GiB cap. Diagnosed (`donor_mean` built a list of `k` full arrays and then a stacked
copy before averaging; `load_effects` materialised a 17,767 × 38,606 row block twice) and fixed by
accumulating in place and chunking the memmap read: the same run now peaks at **4.2 GiB**. Highest
peak anywhere after the fix: **8.2 GiB**. No `cell-eval2` process was started at any point.

---

## 10. What this closes, what it opens, and the honest caveat

**Closed by measurement, on top of the existing closures:**

- **The leaders' climb is not explained by a high-capacity learned transfer model.** Every rung above
  the smallest one loses, and the two biggest ones lose in the two distinct ways this project has
  seen before: replacing the same-target donor signal (132) and clearing on perturbations while
  failing at the line boundary (195). 215's "iterative retraining of a learned transfer model"
  reading is **not supported by any capacity above a per-gene diagonal**.
- **Recipient conditioning on CRISPRi is dead a second time**, now from 7 CRISPRi lines with a
  wrong-recipient null that beats the arm, and with an MLP that zeroes its own conditioning heads.
  184's verdict stands and is now over-determined.
- **Basal-similarity donor weighting stays closed** (010): +2.81%, P(> gate) 0.38.

**Opened, and this is the deliverable that matters:**

> **A per-gene cross-line reliability weight on the emitted effect — 189's `theta`, applied to the
> DIRECTION channel instead of the ordering channel — is worth `+8.9%` relative transfer cosine and
> `+0.0118` `pds` at production's `k = 3`, on the challenge panel, on top of production's own
> whitening `gamma`, with a gene-shuffle null at −7.2%, a label-derangement null at −25.7% and a
> moment-matching control at −3.5%.**

Against AGENTS.md §8.4's corrected bars that is **3.1× the shape-only gate (+2.9%)** and **0.65–0.77×
the amplitude-floor gate (+11.6–13.6%)**. Using §8.2's exchange rate of **+0.00365 overall per +0.01
of cosine** on our board cosine of **0.09027**, a +8.9% relative gain is **+0.0080 absolute → ≈
+0.0029 overall**, which would be about **18% of the +0.0159 gap to `cqawesome`** — *if* the
effect-level gain survives the transfer to the board.

**The caveat, stated as loudly as the result.** This is an **effect-level** measurement on donor
tables. It is **not** scored after `profile_from_effect` + the scorer read-back, so per experiment
156 it measures an object the board never sees; it is **not** a PLOLO number, so it is not directly
priceable through rule 13; and the PLOLO→board transfer factor (**1.17–1.39**, 219) and the direction
arm's **2.51×** propagation multiplier have **not** been applied to it here because neither is
measurable on this instrument. **Nothing in this experiment justifies spending a slot.**

**The next experiment, precisely specified.** Wire the gain in as a per-gene multiplier on the
production effect matrix (a **profile-class** arm, so rule 13's `reach` toll is **0** per 169, but
`mse` moves and must go through 096's `claim` route), fit it leave-one-donor-line-out on non-panel
targets only, and price it on **all six members through PLOLO** in the challenge frame. Two things to
check there that this surface cannot see: (i) the gain interacts with `VCC_ALPHA`, because rescaling
per-gene magnitudes changes the optimal global amplitude and 211 found `claim` already at 100.0% of
its cap; (ii) it changes per-gene effect sizes and therefore the DE call set, so `fid`/`jac` must be
priced, not carried at zero — §8.4's whole lesson.

---

## Files

| file | what |
|---|---|
| `PREREGISTRATION.md` | gate, stop rule and ladder, committed (`2618c79`) before the first fit |
| `substrate224.py` / `surfaces224.py` / `panel224.py` | the 17-table census, the pairwise matrix, the six cached surfaces |
| `known224.py` | the 159 §5.2 known-answer rows |
| `core224.py` | shared primitives; the model-selection rule and why it is load-bearing |
| `features224.py` | DepMap co-essentiality (132's own artifact) + STRING; ESM-2/GenePT are **not on disk** and were not downloaded |
| `ladder224.py` | rungs 0–4 and every null, leave-one-LINE-out |
| `mlp224.py` | rung 5, torch on MPS via `uv run --with torch` (**`pyproject.toml` and `uv.lock` unmodified — verified**) |
| `crosspanel224.py` | the power-vs-substance disambiguation, 5 seeds + the 0/300 fit pool |
| `whiten224.py` | the interaction with production's `gamma` |
| `mechanism224.py` | stability, predictors, the 189 triangulation |
| `leakaudit224.py` | falsifier 4 |
| `analyse224.py` | rung × k, slopes, paired target bootstrap, per-held-line |
| `work/*.json`, `work/analyse-*.txt`, `work/mechanism-*.txt` | every number above |
| `work/reap224-{before,after}.json` | the cloud proof |
