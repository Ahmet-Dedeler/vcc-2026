# 187 — The joint amplitude optimum: `mse` and `nmae` solved together. Amplitude is CLOSED on A/B/C.

**Date:** 2026-09-03. Algebra + board arithmetic only. **No `vcc`, no `cell-eval2`, no build, nothing
written outside this directory.** Peak RSS 0.8 GB (`mse187.py`); the rest is smaller.

## Verdict in six lines

1. **The joint optimum of `mse` + `nmae` in a GLOBAL α is exactly where a15 / a23 sit: 0.982.** Any
   global increase crosses context A's ρ cliff (at α 0.9823) and `mse` loses 10× what `nmae` gains
   per 1%; any decrease loses on both members. Global amplitude is closed — for good.
2. **abc-knee is RETRODICTED**, with zero free parameters at the point of prediction, by the only
   per-context C_p hypothesis that also reproduces a12→a14: the reference's own jackknife on
   **context C is 1.5–2.0% below the NTC-construct proxy** (A and B are at the proxy). Under it the
   knee overshot C's true cliff (α 1.0386, not 1.0479) by 0.9% and paid +0.0103 raw on C, which is
   the whole `mse` miss (model −0.00208 scaled vs board −0.00226). 096's two other pre-registered
   corrections (the h1-measured +1.95% level bias; the CV-0.052 spread) are **falsified** by the pair.
3. **The per-context optimum** — A 0.9823 / B 1.0129 / C 1.0386 — is worth **+0.0011 overall vs a23
   and +0.0004 vs abc-knee** (`mse` +0.0046, `nmae` +0.0020 scaled). Composed gate: **+0.00075 against
   a +0.0089 floor → VETO**, by 12×. Including `nmae` it is still 12% of the floor. **No slot.**
4. **`nmae` is a direction problem, not an amplitude problem.** Within the `mse`-feasible amplitude
   window the collectable `nmae` gain is 0.0012 raw (abc-knee collected 0.00116 of it); the gap to
   the leader is 0.0163 raw and to Aginglab 0.054. The "k ≈ 0.75" reading is an artefact of the
   uniform-k identity: on h1-test the |t|-weighted per-gene k has median **0.12** and p99 **1.11**,
   and the board's own amplitude response is 75% of the uniform-k slope. The credit is
   `⟨sign(t), p⟩` on the confident genes (|p|-weighted sign accuracy 0.62 on h1 vs count accuracy
   0.55); moving it needs a better `q_w`, i.e. 184's territory.
5. **A uniform k increase is safe on `nmae` up to ≈ +20%** (h1: −0.0005 raw at +3%, −0.0015 at +10%,
   turnover at +29%; board: −0.00116 at +3.3%). The binding constraint is the `mse` cliff at
   +0% / +3.1% / +5.8% (A/B/C), far before `nmae` turns.
6. **D/E/F rule:** the knee computed blind from the NTC proxy is expected to gain ≈ 0 (the C-type
   proxy error costs what the knee buys); the expected-value optimum under a 1% cliff uncertainty is
   **knee × (1 − 1.25%)**, worth +0.0009 overall over a global 0.982 with the asymmetric downside
   removed. Recorded for the transition plan; it is not an A/B/C arm.

---

## 1. The two-member model, from first principles

### `mse` — the per-context ρ-cliff identity (096 (a), 170, 180)

For context `c` at amplitude factor `s = α/0.982` on a23's profile:

```
MSE_c(s) − 1 = [ s²A_c + Kc_c − min(claim_c, s²B_c) − 2 s X ] / ΣD
raw = mean_c MSE_c ;  score_mse = mean_c clip((0.987862879 − MSE_c)/0.950162566, 0, 1)
claim_c = Σ_p w_p · E[ min(Ĉ_p, C_p) ]
```

`A_c, B_c, Kc_c, Ĉ_p` are ours, measured on the shipped artifacts (096 `*-stats.npz`, 180
`knee180.json`). `ΣD = 0.229141` and `X_prod = 0.019240` come from the C_p-free LS solve on
a8+a10+a11 (reproduced here, [KA-i]). `X` for a23's profile is solved from a23's own board raw
(one equation: 0.017930). `C_p` — the reference's per-perturbation jackknife — is withheld; 096
substitutes the delete-one jackknife of one 400-cell NTC construct, and the a12 audit named two
corrections **before abc-knee existed**: a level bias (×1.0195 on h1-test) and an across-perturbation
spread (CV 0.052 on h1-test; `min(·, C_p)` is concave so spread always lowers `claim`). This
experiment adds a third: a per-context level factor `f_c`.

Assumptions: (i) `X` and `ΣD` are common across contexts (096 (e); cannot be separated from a pooled
raw); (ii) `Ĉ_p` is flat in `s` (096 measured 0.01% over `s ∈ [0.15, 0.85]`); (iii) `X` carries
linearly with `s` (180 measured κ/s deviation ≤ 7e-3, used κ where recorded); (iv) 24 measured rows
per context stand in for 300 on a23 (180's sample; a12 has 300/300).

### `nmae` — the |t|-weighted identity, heterogeneous k

The spec is `NMAE_p = Σ_S |p − t| / Σ_S |t|`. With `k_g = |p_g|/|t_g|`, |t|-weights `w_g`, and the
kink at `k = 1`:

```
nmae = 1 − 2 q_w W_hi − (2 q_lo − 1) K_lo + K_hi              (exact, gene-wise)
d nmae / d ln α = −(2 q_lo − 1) K_lo + K_hi                      (first order)
K_lo = Σ_{k<1} w k,  K_hi = Σ_{k>1} w k,  W_hi = Σ_{k>1} w ;  q_lo, q_w = k-weighted sign accuracies
```

Uniform `k` and gene-independent `q` collapse this to CURRENT_STATE's `q|k−1| + (1−q)(k+1)`. On the
board I use the measured local law on the current pipeline, `raw_nmae(α) = 0.944330 − 0.0416 ·
mean_c ln(α_c/0.982)` (a12→a14: 0.0471; a23→abc-knee: 0.0361; scaled = −Δraw / 0.5963 from our own
board map, [KA-a]), with h1-test's shape (turnover at s ≈ 1.29) as the saturation prior. Every
optimum below lies inside ±6%, where the linear law is what abc-knee measured.

`pds` is scale-invariant on this pipeline (180: PLOLO −0.000076 raw); `fid`/`jac`/`reach` are carried
at 0 with 180's twin band (reach ±0.006–0.011 scaled, fid ±0.0003, jac ±0.0002).

---

## 2. Retrodiction table (rule 10 / audit-181 rule 1)

Known-answer rows on the implementation, then the board points. "err" = model − board.

| # | row | model | board / record | err | verdict |
|---|---|---:|---:|---:|---|
| KA-i | C_p-free LS solve a8+a10+a11 → `X`, `ΣD` | 0.019240 / 0.229141 | 0.019240 / 0.229141 | ≤ 6e-5 raw | **PASS** |
| KA-ii | a12 raw under proxy / level / spread / both | 0.898549 / 0.890600 / 0.904041 / 0.894773 | 0.898551 / 0.890602 / 0.904073 / 0.894727 (A12-AUDIT §4) | ≤ 5e-5 | **PASS** ×4 |
| KA-iii | a23 → abc-knee under proxy | +0.005025 scaled | 180: +0.005025 | 0 | **PASS** (reproduces the failed prediction) |
| KA-a | our (raw, scaled) `nmae` rows, one affine map | b 1.000439, span 0.5963 | 10 rows | 5e-4 raw max | PASS (per-context spans differ) |
| KA-b | h1-test `nmae` at s = 1 from 170's DE table | 0.959881 | 180 anchor 0.959881 | 0.0 | **PASS** (target removed; without it 0.957096) |
| KA-c | h1-test `nmae` at s = 1.0315 / 1.0671 by lfc scaling | 0.959356 / 0.958833 | 180 re-emitted 0.959180 / 0.958532 | +1.8e-4 / +3.0e-4 | PASS to 3e-4; scaling under-states the gain 10% |
| **B1** | **a12 → a14** Δraw (α 1.0108 → 0.982, one variable) | proxy **−0.018544**; level −0.010595; spread −0.024037; both −0.014768 | **−0.018574** | **+0.00003** / +0.0080 / −0.0055 / +0.0038 | **proxy PASS; level, spread, both FAIL** |
| **B2** | **a23 → abc-knee** Δscaled `mse` | proxy +0.005025; level +0.005025; spread −0.001857; both +0.001309 | **−0.002259** | +0.0073 / +0.0073 / **+0.0004** / +0.0036 | **spread PASS; proxy, level, both FAIL** |
| **B1+B2** | **per-context factor `f = (1, 1, 0.9825)`** on the proxy | B1 −0.018544; B2 **−0.002076** | −0.018574; −0.002259 | **+0.00003; +0.00018** | **PASS on both** — the only hypothesis that does |
| B3 | a23 → abc-knee Δscaled `nmae` (linear board law) | +0.002238 | +0.002081 | +0.00016 | PASS |
| B4 | a23 → abc-knee overall from `mse`+`nmae` | +0.000027 | −0.000030 (board mse+nmae; +0.000684 with reach/fid/jac) | +0.00006 | PASS |
| a15 | board raw 0.873805 | anchor: `X` solved (170: 0.017766) | — | — | X-anchor, not a retrodiction; 164's "0.88228" is the identity with the a8-family `X` before the CD4 uplift (Δ−0.0085 = CD4's +4.9% on `X`, a profile change this model does not price) |
| a19 | `mse` −0.09591 scaled (overall −0.01352) | out of domain | — | — | **not retrodictable by an amplitude model**: a19 moved the operator (gamma 0.80, ridge 0.007); 170 measured the gamma part at 21% (+0.0190 raw); the ridge changes `A`, `B`, `X` and none was measured on the challenge panel. Sign only: any ridge that raises `A/X` above the cliff is charged at full `2sA/ΣD` |

**Why B1 and B2 discriminate.** a12→a14 straddles A's cliff (a14 lands on it to 0.02%) and B's
(a12 was 0.3% above it) on the a12 profile; abc-knee lands exactly on B's and C's proxy cliffs on
the a23 profile. A spread or a lower level on **B** makes a12 further above B's cliff and predicts a
larger a12→a14 gain than the board delivered (spread: −0.0055). A lower level on **C** leaves a12→a14
untouched (a12 is below C's cliff either way) and charges the knee. Single-context scans
(`cliff187.py`): `f_C` alone zeros B2 at 0.9825 with B1 at +0.00003; `f_B` alone cannot zero both
(B1 −0.0046 at the `f_B` that fixes B2); `f_A` is bounded above at 1.000 by B1 and is otherwise
unconstrained by the deltas (0.99–1.00 also fits the a12/a14 levels). The joint region that fits both
deltas to ±0.001 is `f_B ∈ [0.995, 1.0075]`, `f_C ∈ [0.980, 0.985]`, `f_A ≤ 1`.

**What the h1-measured corrections did not transfer.** The +1.95% level bias (real perturbation groups
vs 400-cell control draws on h1-test) is contradicted on A/B (B1 +0.0080 if applied) and has the
wrong sign on C. The CV-0.052 spread is contradicted on B. **The challenge reference's `C_p` sits at
the NTC proxy on A and B and 1.5–2% below it on C** — C is the context whose NTC constructs are the
most heterogeneous (sd 0.92 vs 0.38 / 0.47), so the NTC pool overstating the perturbed groups'
within-group dispersion there is physically plausible and unverifiable from what we hold.

---

## 3. Where the cliffs are, and the joint optimum

a23's profile, retrodicted model `f = (1, 1, 0.9825)`, `X23 = 0.017930`:

| ctx | `A_c` | `B_c` | `claim_c` | ρ at a23 | **cliff (α)** | proxy said | abc-knee shipped | a23 raw |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| A | 0.234830 | 0.233593 | 0.233754 | 0.9993 | **0.9823** | 0.9823 | 0.98234 (on it) | 0.87108 |
| B | 0.224778 | 0.223624 | 0.237923 | 0.9399 | **1.0129** | 1.0129 | 1.01291 (on it) | 0.87135 |
| C | 0.234941 | 0.233728 | 0.261471 | 0.8939 | **1.0386** | 1.0479 | 1.04786 (**+0.9% over**) | 0.87463 |

abc-knee per-context Δraw under the model: A −0.00005, B −0.00458, **C +0.01034** → mean +0.0019 raw
(board +0.0021).

### Global α (sum of `mse` and `nmae`, Δ vs a23, overall units)

| α | 0.95 | 0.97 | 0.98 | **0.982** | 0.985 | 0.99 | 1.00 | 1.01 | 1.02 | 1.05 | 1.10 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Δ`mse` scaled | −0.0050 | −0.0019 | −0.0003 | 0 | −0.0015 | −0.0044 | −0.0102 | −0.0161 | −0.0272 | −0.0699 | −0.1216 |
| Δ`nmae` scaled | −0.0023 | −0.0009 | −0.0001 | 0 | +0.0002 | +0.0006 | +0.0013 | +0.0020 | +0.0027 | +0.0047 | +0.0079 |
| Δ overall | −0.0012 | −0.0005 | −0.0001 | **0** | −0.0002 | −0.0006 | −0.0015 | −0.0024 | −0.0041 | −0.0109 | −0.0189 |

**The global optimum is 0.982, to 4 decimals, under both the retrodicted and the proxy model.** Below
the cliffs `mse` and `nmae` both improve with α (slopes +0.155 and +0.070 scaled per unit ln α);
above A's cliff `mse` loses 2.0 scaled per unit ln α — 30× what `nmae` gains. `nmae` never justifies
crossing a cliff, so adding it to the objective moves nothing: the joint optimum is the `mse` cliff.

### Per-context α

| arm | A | B | C | Δraw `mse` | Δ`mse` sc | Δ`nmae` sc | **Δ overall vs a23** | vs abc-knee |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| **on the retrodicted cliffs** | 0.9823 | 1.0129 | 1.0386 | −0.00435 | +0.00458 | +0.00203 | **+0.00110** | **+0.00042** |
| cliffs × (1 − 0.5%) | 0.9774 | 1.0078 | 1.0335 | −0.00360 | +0.00379 | +0.00168 | +0.00091 | +0.00023 |
| cliffs × (1 − 1%) | 0.9725 | 1.0028 | 1.0283 | −0.00285 | +0.00300 | +0.00133 | +0.00072 | +0.00004 |
| cliffs × (1 − 2%) | 0.9627 | 0.9927 | 1.0179 | −0.00135 | +0.00142 | +0.00062 | +0.00034 | −0.00034 |
| abc-knee as shipped | 0.9823 | 1.0129 | 1.0479 | +0.00197 | −0.00208 | +0.00224 | +0.00003 | −0.00066 (board: 0) |
| a23 (global 0.982) | 0.982 | 0.982 | 0.982 | 0 | 0 | 0 | 0 | −0.00068 |

**Per-context α does help once `nmae` is in the objective — by +0.0004 over what is already on the
board.** The whole per-context pot on A/B/C is +0.0011 overall vs a23, of which abc-knee banked
+0.0007 (through `nmae` and `reach`, not `mse`).

### Under cliff-position uncertainty (expected value, `mse` + `nmae`)

Cliff of each context ~ Normal in ln α around the retrodicted value, sd from the retrodiction region
(A 0.25%, B 0.04%, C 0.13%): best α = A 0.9794 / B 1.0124 / C 1.0371 (margins −0.3% / −0.05% /
−0.15%), expected +0.00111 vs a23. **Blind** (D/E/F-like, sd 1% per context because `f_c` is unknowable
before the round): best margin **−1.24% to −1.29%** below the computed knee on every context, expected
+0.00090 vs a global 0.982; **at the raw knee the expectation is ≈ 0 (A −0.00002, B −0.00009,
C +0.00018)** — the proxy's asymmetric risk eats what the knee buys, which is what abc-knee measured.

---

## 4. The composed slot gate (rule 13, audit-181 amendment: profile-class, no emission toll)

Arm: a23 rebuilt at A 0.9823 / B 1.0129 / C 1.0386.

```
PLOLO Δpds/6 × 1.011      −0.000076/6 × 1.011   =  −0.000013     (180's measurement)
096 Δmse/6 (amplitude)    +0.00458/6            =  +0.000763
predicted                                       =  +0.000750
SHIP only if > +0.0089  (amplitude class)       →  VETO   (8% of the floor)
with nmae counted                               =  +0.001089  (12% of the floor)
```

Six-member preregistration for the record — **do not build, do not submit**:

| member | predicted Δ vs a23 | vs abc-knee | mechanism / band |
|---|---:|---:|---|
| `pds` | 0.000 ± 0.001 | 0 | scale-invariant under `VCC_CD4_TARGET_F` (180 §7) |
| `mse` | **+0.00458** | +0.00666 | C back below its cliff (+0.0103 raw on C removed), B/A unchanged |
| `nmae` | +0.00203 | −0.00021 | 0.0416 raw per ln α; C's amplitude 0.9% lower than abc-knee |
| `fid` / `jac` | 0 ± 0.0003 / 0 ± 0.0002 | | twin noise (180) |
| `reach` | 0 ± 0.01 | | unresolved; abc-knee's +0.0035 is inside seed |
| **overall** | **+0.0011** | **+0.0004** | |

The builder cannot express this arm today: `VCC_ALPHA=knee` solves the knee from the NTC proxy and
has no per-context `C_p` factor and no explicit per-context α input. Realising it would need a
builder change (a `C_p` scale per context, `A:1,B:1,C:0.9825`) — a change outside this directory,
not made, and not worth making for +0.0004.

---

## 5. `nmae`: k versus q

### The uniform-k identity reading (as the mandate framed it)

`nmae = 1 − k(2q − 1)` for `k ≤ 1`; at `k = 1` it is `2(1 − q)`.

| team | raw `nmae` | credit `1 − nmae` | k at q = 0.537 | q at k = 1 | q at k = 0.75 |
|---|---:|---:|---:|---:|---:|
| ours (abc-knee) | 0.9432 | 0.0568 | **0.768** | 0.5284 | 0.5379 |
| ours (a23) | 0.9443 | 0.0557 | 0.752 | 0.5278 | 0.5371 |
| cqawesome (leader) | 0.9269 | 0.0731 | 0.988 | 0.5366 | 0.5487 |
| Aginglab.com | 0.8892 | 0.1108 | 1.497 (impossible) | **0.5554** | 0.5739 |

Read literally: the leader is "us at k ≈ 1", Aginglab needs q ≥ 0.555 even at perfect amplitude. The
literal reading is wrong for us, and the evidence is our own amplitude response:

### What the board and h1-test actually say

| quantity | h1-test (a23 geometry, BH 0.05, target removed) | board (current pipeline) |
|---|---:|---:|
| credit `1 − nmae` | 0.0401 | 0.0557 |
| `d nmae / d ln α` | −0.0211 (first-order identity: −0.0211; 098 measured −0.0179 on a12) | **−0.0416** (a12→a14 −0.0471, a23→knee −0.0361) |
| uniform-k=0.75 identity slope | | −0.0555 (board is 75% of it) |
| `K = Σ|p|/Σ|t|` | **0.199** | not observable |
| `K_lo`, `K_hi`, `W_hi` | 0.180, 0.0195, 0.0153 | `2 q_w W_hi = 0.0141` → W_hi ≈ 1.2% |
| sign accuracy: count / |p|-weighted / below-kink | 0.552 / **0.621** / 0.613 | 0.537 count (172) |
| |t|-weighted per-gene k quantiles | p10 0.012, **p50 0.123**, p90 0.486, p95 0.646, **p99 1.106** | |
| `nmae(s)` argmin / gain | s = 1.29 / −0.0027 raw | linear to ±6%; shape prior beyond |

There is no "k = 0.75". The median gene on the reference pool is predicted at **12%** of its true
magnitude, one gene in a hundred is already past the kink, and the credit is a directional inner
product `⟨sign(t), p⟩_w` (098 (d1)) concentrated in the confident genes, where we are right 62% of
the time against 55% by count. A global α cannot bring the median gene to k = 1 without pushing the
top decile 8× past it; h1's curve turns over at +29% having gained 0.0027 raw (6% of the leader gap).

### The decomposition of our deficit

| | vs leader (0.0163 raw) | vs Aginglab (0.0540 raw) |
|---|---:|---:|
| amplitude, inside the `mse`-feasible window (B +3.1%, C +5.8%) — **measured by abc-knee** | **0.0012 (7%)** | 0.0012 (2%) |
| amplitude, `mse` sacrificed (α → 1.29, h1 shape × board slope) | ≈ 0.005 (30%) — and costs −0.019 overall on `mse` | ≈ 0.005 (9%) |
| **direction (`q_w` on the reference pool)** | **≥ 0.011 (≥ 70%)** | **≥ 0.049 (≥ 90%)** |

Field check (206 published teams, `pds > 0.70`, `nmae ∈ [0.8, 1.05]`, n = 190): `nmae = 1.255 − 0.391
pds`, σ 0.027. At our `pds` the line reads **0.9276**; we sit **+0.0155** above it; of the 25 teams
within ±0.02 `pds` of ours, **22 have better `nmae`** (median 0.9291 — the mandate's +0.0141). The
leader sits **+0.0116 above the line at their own `pds`**: their `nmae` edge over us is 75% their
`pds` edge (0.0123 of 0.0163) and not a calibration edge. Aginglab is **−0.045 below** the line —
the one team with a genuine `nmae`-specific mechanism (and the only positive `fid` in the top 12).
Amplitude explains at most 0.0012 of our +0.0155 residual at fixed `pds`; **the rest is `q_w` at fixed
`pds`** — exactly the trade the whitener makes (raw-space alignment down 32%, discrimination up; §4.2
of AGENTS.md).

**Answer to the mandate's question:** the field beats us on **q**. `nmae` is a direction member
(098's headline, now with the amplitude channel priced on the board: 0.0012 raw, already banked by
abc-knee).

---

## 6. Is a uniform k increase safe? (item 5)

Charged side: genes with `k_g > 1` pay `+k_g` per unit ln α regardless of sign; wrong-sign genes above
the kink pay `|p| + |t|`.

| uniform increase | h1 Δ`nmae` raw | W_hi after | K_hi after | board |
|---|---:|---:|---:|---|
| +3% | −0.00050 | 0.0170 | 0.0219 | abc-knee +3.3% mean: **−0.00116** (better) |
| +5% | −0.00081 | 0.0177 | 0.0230 | |
| +10% | −0.00148 | 0.0200 | 0.0265 | |
| +20% | −0.00241 | 0.0250 | 0.0341 | |
| +29% | −0.00269 (turnover) | | | |

Safe on `nmae` throughout the range any `mse`-feasible arm can reach: the |t|-mass above the kink
grows 1.5% → 2.0% at +10%, and the below-kink credit dominates until +29%. The 02:15 constraint holds
in practice; the binding constraint on amplitude is the `mse` cliff at +0% / +3.1% / +5.8%, an order of
magnitude earlier.

---

## 7. What this closes, and what it hands on

- **Global amplitude on A/B/C: closed.** 0.982 is the joint optimum of the two amplitude members to 4
  decimals; a23 / abc-knee already sit on it.
- **Per-context amplitude on A/B/C: closed for slots.** The entire remaining pot is +0.0011 vs a23 and
  +0.0004 vs the displayed entry, 12× under the gate floor. Do not build the corrected knee.
- **The `C_p` proxy: now board-calibrated on A/B/C.** `f = (1.00, 1.00, 0.9825)`; the h1-measured
  level bias and spread do not transfer. This replaces the 00:00 entry's "the proxy does not
  transfer" with a measured correction — but the correction is context-specific and cannot be known
  for D/E/F before the round.
- **D/E/F transition plan (for whoever owns it):** ship `knee × (1 − 1.25%)`, not the raw knee. Under
  a 1% cliff uncertainty the raw knee's expected gain is ≈ 0 and its downside is 30× its upside per
  1%; the haircut keeps +0.0009 of expected overall over a global 0.982. Contexts whose NTC constructs
  are unusually heterogeneous (C-like, sd/mean > 5%) deserve the larger margin.
- **`nmae`: direction (184).** The amplitude channel is priced at 0.0012 raw and banked. The 0.015 raw
  we lose to teams at our `pds` is `q_w` at fixed `pds`; the leader has no `nmae`-specific mechanism;
  Aginglab does.
- **Standing instrument note:** the uniform-k identity is right as an identity and wrong as a
  measurement of k. Use the board's amplitude response (−0.0416 raw per ln α) and h1's k-quantiles.

## 8. Files

| file | what it does |
|---|---|
| `mse187.py` | the per-context identity; KA-i/ii/iii; a14 and abc-knee under the four 096 hypotheses; (bias, CV) sweep → `work/mse187.{json,log}` |
| `cliff187.py` | per-context `f_c` / `cv_c` scans and the `(f_B, f_C)` region that fits both board deltas → `work/cliff187.{json,log}` |
| `nmae187.py` | board scale map, the four amplitude pairs, h1-test k-distribution and `nmae(s)` curve, identity readings → `work/nmae187.{json,log}` |
| `joint187.py` | the joint global / per-context optimum, expected value under cliff uncertainty, the gate → `work/joint187.{json,log}` |

Inputs read, never written: `experiments/096-mse-floor/work/{a8,a10,a11,a12}-stats.npz`, `cp-proxy-ntc.npz`,
`ctrl-stats.npz`; `experiments/180-def-readiness/work/knee180.json`; `experiments/170-gamma-mse/work/de_g0.65.npz`;
`data/leaderboard/2026-09-02T2359Z.json` and earlier snapshots.
