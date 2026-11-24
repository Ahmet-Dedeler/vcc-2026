# 246 — The direction channel, audited. Five of the six named hypotheses are CLOSED or worthless; the one OPEN construction — projecting the target-shared subspace out of an un-whitened residual — is mechanistically REAL but pays nowhere on the mandated object. **Applied instead to the shipped CD4 residual it buys `pds` +0.0123 raw (3/3 cancer folds, 8σ) for −1.3 % macro cosine, at 11–30× the dose axis's exchange rate; priced +0.0023 … +0.0082 overall. BUILD-WORTHY under the shape-only gate, conditional on the board's `pds` transfer.**

**Date 2026-09-04.** Deterministic (seeds `20726` for every derangement, 218's; `20260904` (+offsets) for
every random subspace). **No `vcc`, no `cell-eval2`, no build, no submission, no slot proposal, no
production change, nothing written outside `experiments/246-direction-hypotheses/`** (the only files
outside it that this run *read* are 218's, 222's, 216's, 062's, 162's, 188's and 203's artifacts and
the CD4/H1 banks; `docs/*` and `scripts/*` files with newer mtimes were touched by 240/241/243/245 or
the lead, not by anything here — every write in these scripts goes through `lib246.WORK`).
Peak RSS **4.14–4.96 GiB** per fold process (`/usr/bin/time -l`, one fold per process, 12 runs) and
**8.53 GiB** for `six246.py` (the operator load), **`0 swaps` on all 13 runs**; memory gate at every
launch: `memory_pressure` **85–87 % free**, `vm_stat` free+inactive **30 GiB**, while 241's
`mse241.py` (3.9 GB) was live. 9 scripts, **560 scored arm×fold rows** on 162's PLOLO (4 folds).

---

## 0. The answer in fourteen lines

1. **All 19 known-answer rows pass, every one at `0.0e+00`** (§2): 218's `f0`, `self0.1`,
   `cd40.34+none0` and `cd40.34+h1PAIRCARRY0.28` per fold — so the mandated row, **macro
   +7.6803 % / `pds` +0.004887 on the three cancer folds**, reproduces to the digit; 216/222's exact
   `nmae`, precision and delivered-norm rows through `six246.py`; and the k = 0 deflation is
   bit-identical to the shipped stack.
2. **Phase 1 (§1).** Hypotheses 3 (gating beyond 10 CPM) and 4 (sign-consistent consensus) are
   **CLOSED** by 088/213/164 and 232/234 respectively, with the exact construction inside the closing
   experiment. Hypothesis 5 (per-source dose, joint `(f, g)`) is OPEN on paper and **worthless**:
   the only source that wants a higher dose is `h1_2025`, which goes **20 → 3 targets in October**
   (229 §octo), and the CD4 axis is flat to ±0.0003 (216 §6.0; re-measured here, §3.6). Hypothesis 2
   (per-target dose) is unmeasured, adds no information, and has no mechanism for its sign — OPEN-LOW,
   design given, not run. **Hypothesis 1 is the only genuinely OPEN one**: no experiment had ever
   deflated the target-shared subspace of an *un-whitened residual* — 156/175/204/209 all deflated
   the donors or the render, objects the whitener already handles.
3. **Hypothesis 1 as posed — 099/218's self-residual with its shared subspace projected out, on top of
   the shipped stack — is NOT build-worthy.** The projection converts its `pds` cost (−0.0039 at
   f = 0.10) into a gain (+0.0017 at k = 2) **but the cosine gain vanishes** (macro +2.17 % → +0.23 %,
   aggregate −0.30 %): the stack already carries the un-whitened relaxation. Best price **+0.0014**
   (k = 2, f = 0.10) against the +0.0019 gate; **NO-GO** (§3.1).
4. **The mechanism is real, measured on the bare base** (§3.2): 218's self-residual at f = 0.10 gives
   **+5.35 % cosine / `pds` −0.0097**; with its top-5 shared directions removed it gives **+2.18 % /
   `pds` +0.0019** — cosine gain 41 % retained, `pds` cost gone — while an **energy-matched random
   subspace** (63 dims, the same 25 % of energy) leaves both untouched (+5.23 % / −0.0099). **The
   `pds` cost of an un-whitened residual sits in its top-k target-shared directions; its cosine gain
   only partly does.**
5. **Where that pays on the shipped stack is the CD4 residual itself** (§3.3): the shipped CD4 residual
   enters *un-whitened* at f = 0.34, so nothing ever removes its target-shared component — its top-5
   directions carry **20 % of its energy** across 251 rows. Removing them at the shipped fraction:
   **`pds` +0.012271 raw (3/3 folds, 8.2 LOO-σ), macro cosine −1.286 %, aggregate −1.194 %,
   q_w +0.00058**. Per fold: +0.0094 / +0.0095 / +0.0180; `kolf2` apart +0.0038.
6. **Every falsifier passes** (§3.4): two **energy-matched random subspaces** (51–52 dims, the same
   20 % of energy) read `pds` **−0.0001 / −0.0050** at cosine **−3.7 / −4.1 %**; three
   **shuffled-projector** controls read ≤ +0.0013; the **split-half projector** (basis fitted on the
   other half of CD4's rows) keeps **+0.0104 … +0.0118**, so the subspace is a property of CD4, not of
   the rows it is applied to; **non-carrying rows move ≤ 0.026 pp**; leave-top-3-out retains **56 %**
   of Δ`pds` (`h1union` retained 14.6 %, 222 §9); leave-one-fold-out k selection prices
   **+0.0015 / +0.0028 / +0.0069** at `pds` transfer 0.7× / 1.0× / 2.0×.
7. **The exchange rate is what makes it an arm** (§3.5): on the same folds and base, the CD4 **dose**
   axis (216's lever) buys `pds` **+0.0022 for −2.37 %** aggregate cosine (f 0.34 → 0.28); the
   deflation axis buys **+0.0123 for −1.19 %** — **11× per unit cosine at k = 5, 30× at k = 3**.
   212's gamma axis cancels to 3.8 % of the `mse` leg; this axis does not.
8. **§8.1d — signs or magnitudes?** The removed component is mostly *unmatched* magnitude on the
   reference's DE pool: 188's exact identity gives **`nmae` −0.00046 raw (improves), e = ΔM/ΔU =
   0.25 < ½**, and **q_w rises +0.00058**. The rate card charges `nmae` −0.00011 from the cosine
   loss; the exact identity says +0.00013. Both inside ±0.0004; **reported as a magnitude bound**.
9. **Six members, never at zero** (§4): cosine legs (`mse` −0.00042, `reach` −0.00028, `nmae`
   −0.00011, `fid` −0.00010, `jac` −0.00001) = **−0.00093**; `pds` **+0.0032 / +0.0046 / +0.0091**
   at transfer 0.7× (`h1union` realised) / 1.0× / 2.0× (§8.4's `cd4f034` factor). **Total +0.0023 /
   +0.0037 / +0.0082.** The re-solved-dose variant (`svd5_resolved`, CD4 fraction re-solved to 0.396)
   reads **+0.0026 / +0.0040 / +0.0086** at cosine −0.46 % aggregate.
10. **Against the gates:** above the **+0.0019 shape-only** gate at every transfer factor observed
    (1.2× … 4.3×; §8.1b says a direction gain of +0.002 is readable on the deterministic board);
    **below the +0.0089 amplitude gate at all of them** (0.26× … 0.92×). Delivered amplitude moves
    +0.11 % (six246 `dlv/a15` 1.0011 vs 1.0033) — the arm is in `cd4f034`'s profile class (α, γ,
    λ, emission untouched; one preserve scalar), i.e. the shape-only row.
11. **The self-generated hypothesis — CD4's partial-condition rows** (18–23 rows per fold that
    production's complete-mean3 rule zeroes) — is real but small: **+0.59 % macro, `pds` +0.0018
    (3/3, 3/3)**, shuffled ≈ 0, row-matched control (the consensus on the same rows and support)
    +0.57 % / −0.0001 — so the *cosine* is not CD4-specific, the `pds` is; price **+0.0009 …
    +0.0018**, concentrated (leave-top-3-out −6.6 %). **NOT build-worthy alone.** Composed with the
    deflation (`E_cd4_svd5+partial`): `pds` +0.0140, price **+0.0032 / +0.0047 / +0.0100**.
12. **What information it adds — the mandate's filter, answered honestly:** none. It **removes** the
    T-cell-programme hub that the CD4 residual re-injects into all 251 rows — the component the
    whitener strips from the donors (209: 12.5 % → 0.14 %) and nothing strips from CD4 or H1
    because they bypass the operator (207 §5). This is exactly STRATEGY-2026-09-04 §4.5's second
    admissible form: *a thing tested in the wrong channel* — shared-programme removal was measured
    on the donors and the render (156, 175, 204, 209), never on the one object the operator never sees.
13. **`kolf2` apart, every sign holds**: `pds` +0.0017 … +0.0049 across k, cosine −2.5 … −5.8 %.
14. **Verdict: BUILD-WORTHY (conditional), one variable, shape-only class** — `B2_cd4_svd5`
    (or its re-solved twin). The condition is stated in §5: the board's `pds` transfer for a *gain*
    has one calibration point (0.70×, `h1union`) and it is the whole margin. **NOT build-worthy:**
    hypotheses 1-as-posed, 2, 3, 4, 5, the partial rows alone.

---

## 1. Phase 1 — the audit

Filter (docs/STRATEGY-2026-09-04.md §4.5): *a new direction hypothesis must say which information it
adds that the donors + CD4 + H1 do not contain, or explain why one of the tested sources was tested in
the wrong channel.* Closed list: `docs/CURRENT_STATE.md` ("Closed, do not re-open") + the RESULTS
headlines (239's `headlines.txt`), each closing experiment read at the headline.

| # | hypothesis | mechanism | information it adds | signs or magnitudes (§8.1d) | expected `pds` | status · closing experiment (cited) | PLOLO design + known-answer row | cost |
|---|---|---|---|---|---|---|---|---|
| **1** | **Specific-subspace relaxation**: 099/218's self-residual (+4.61 % cos 4-fold, −0.0078 `pds`; **+5.35 % / −0.0097 on the 3 cancer folds**, `resid218.json`) with the component shared across the panel rows projected out — per-gene centring across carried rows, then the top-k right-singular directions | the whitener suppresses the shared programme; an un-whitened residual re-adds it and pays hubness; removing only the shared part should keep the specific cosine and not pay | **none new** — changes how much held information the operator discards, on the subspace where 193/212 say discarding buys nothing | unknown a priori; **measured: the removed part is unmatched magnitude (e = 0.25), q_w up** | **+** if the cost is in the shared part | **OPEN — nobody measured THIS construction.** Neighbours: 209 deflated the *Orion-pair excess* from the **donors** pre-whitening (rank 2 = the long-gene axis, "the whitener has already taken it, ~0.1 % left", every arm negative); 175/178 per-donor rank-1 deflation of the **donors** (−2.3 … −2.7 % cos); 204 §3.2 blob deflation of the **consensus** (<1 % change, "already deflated"); 156 §6 centroid removal of the **render** (+0.00000); 193 map-shape families incl. `head^1.2` = γ 0.78 on the top half (net −0.0021 on `nmae`); 212 the γ axis (knife-edge, cancels to 3.8 %). None touched the un-whitened residual. | **Measured, §3.** KA: 218 `self0.1` per fold; k = 0 ≡ stack | 1 h |
| **2** | **Per-target dose** for CD4/H1 from the target's cross-donor agreement, replacing global f / g | rows where donors agree need less/more residual | none — agreement is held (it is the ordering key, 097/013) | a per-row mix with row norm preserved is direction-only; un-preserved it is 158's per-target *scale* (`nmae` −0.0023 … −0.0245, CLOSED) | no mechanism fixes the sign: more residual where donors agree (CD4 aligned with a strong shared component) or where they disagree (consensus noisy)? Both are fits | **OPEN-LOW, not run.** Neighbours: 158 (per-target magnitude: NO-GO), 204 (equal donor weighting is the cosine **maximum**; `a_i` not identified), 216 §4.5 (dose optimum per *context* spread 3.2 %, curve flat to 6 decimals), and the **per-group** dose edit spec'd 2026-09-04 (`docs/NEXT-BUILD-SPEC.md` §6, **+0.00099 PLOLO**) — which is per source-*group* (H1 pair vs stack4), belongs to 229's stack (**−0.00219**), and is not per target | design: `f_i = f·(1 + β·z(agreement_i))`, row norms restored to the global-f row norms, β ∈ {±0.25, ±0.5}, LOFO on β; KA β = 0 ≡ stack. Fails the filter; deferred | 1 h |
| **3** | **Recipient-basal gating beyond 10 CPM** | an effect on a gene the recipient does not express cannot be real | none — the renderer gate already IS recipient-basal gating (`gate = cpm_recipient > 10 & b_ctrl > 0`, 162 `plolo.py`; production `--min-cpm 10`) | membership of low-CPM genes is a magnitude/call-set change | — | **CLOSED.** 088b's `min_cpm` sweep on h1-test (094 §"min_cpm gate at its peak"): **8 → +0.0028, 10 → +0.00455, 15 → −0.0056, 30 → −0.0129** vs 5; 006: cosine optimum at 20 CPM on two harnesses, "above ~50 discards real signal"; **213**: the soft form (`log(p+c)`) has cosine max **1.0159×** at c ≈ 36–40 CPM but the **`pds` maximiser is the incumbent** and the composite is **+0.00103 < +0.0019**; 164 §7: `min_cpm 5` NO-GO; 213 F2: deranging the recipient basal costs **−0.070 raw `pds`** (the basal is already load-bearing). STRATEGY addendum 8 lists `min_cpm` closed | not run | — |
| **4** | **Sign-consistent consensus**: carry magnitude only where ≥ 2 donors agree in sign | contested coordinates are coin flips (232: 0.523) | none | zeroing a magnitude is a magnitude arm; and a zero counts as a WRONG sign on the DE pool | ambiguous | **CLOSED, inside 232/234.** At production arity 3, "≥ 2 of 3 agree" is **vacuous** (a majority always exists); the non-vacuous form — zero where not unanimous — is 232's **`maj_strict`**, which scores |t|-weighted sign accuracy **0.324 / 0.395 / 0.394 vs 0.538 / 0.572 / 0.571** (−0.18 … −0.21); through 234's one-sided rate (−0.004827 per −0.01 lost) that is ≈ **−0.09 overall**. 232's majority at arity 3 is worse than the mean (q_w −0.00095, **−0.00375 overall**); 234: **every pool-blind flip arm negative on all six members**, 99.5 % of the axis is off-pool | not run | — |
| **5a** | **Per-source dose for the union** (Pan-2026 H1 and `h1_2025` separately) | the two sources have different carrying-row norms and different `pds` maximisers | none (a re-tune) | — | + for `h1_2025` (monotone to 0.45), flat for `gse295214` (max at g = 0.17) | **OPEN, worthless.** From 218's `compose218.json` (cancer means on `cd4@0.34`): `gse295214_h1` alone @0.28 **+7.98 % / `pds` +0.00167**, @0.45 +12.03 % / +0.00007; `h1_2025` alone @0.28 **+4.61 % / +0.00300**, @0.45 +7.07 % / **+0.00394 still rising**; union @0.28 +10.31 % / +0.00489, @0.45 +15.53 % / +0.00432. A split (gse@0.28 + h1_2025@0.45) is at best ≈ +0.001 `pds` at equal cosine ≈ **+0.0008 overall** — and **`h1_2025` covers 20 → 3 targets on a disjoint panel** (229 §7/§octo: "nearly worthless in October"). 229 measured per-*group* doses, not this | adjudicated from 218's rows; not re-run | — |
| **5b** | **Joint `(f, g)` grid** | CD4's argmax may shift with the union on | none | — | 216: `pds` turns over at f ≈ 0.20 | **CLOSED by measurement.** 216 §6.0 (exact six-member): f = 0.28 **−0.00024**, f = 0.40 **−0.00029** vs 0.34 — flat top; 218's 2×7 grid: (0.10, 0.28) has `pds` +0.0028 above (0.34, 0.28) at **−12 % cosine**; **re-measured here on the stack (§3.6): f 0.34 → 0.28 reads −2.37 % aggregate / `pds` +0.0022, price −0.0006 … +0.0005** | `DOSE_cd4@0.28/0.20`, §3.6 | done |
| **6a** | *(self-generated)* **CD4's partial-condition rows** — targets carried by 1–2 of the 3 CD4 conditions, which `cd4_for` / production's complete-mean3 zeroes (**24/300 on the challenge panel**, 222 §11; 18–23 per fold) | information held and discarded; the `h1union` shape (an estimate on rows that had none — the only shape that raised `pds`) | **yes**: CD4 rest/stim effects on 6–8 % of the panel that no residual reaches | as a residual on empty rows, direction | + | **OPEN — measured, §3.7.** Real but small and concentrated | KA: complete rows untouched | done |
| 6b | recipient control-cell covariance / co-expression / line-specific subspace | the only context-specific data | — | — | — | **CLOSED** 155 (null, split-half r = 0.028), 160 ("unusable: the capture IS the line's mean response, the hub `pds` punishes"), 174 (null; STRATEGY-2026-09-04 §4.5 table) | — | — |
| 6c | a fifth source | — | — | — | — | **CLOSED** 225 (space swept, saturation 1.84×), 218 (alignment, not independence, predicts contribution; kolf2 below the self-residual floor), 238 (October asset dead) | — | — |
| 6d | the operator fitted on non-essential rows | — | — | — | — | **already shipped** since a15 (`basis_prune: essential`, 091 +0.002409; STRATEGY-2026-09-04 §4b retraction) | — | — |
| 6e | partially whitening the residual | a residual pushed through the operator is flat (207 `wres` +0.5 … +0.9 %) | — | — | — | the deflation IS the partial whitening that pays: full whitening kills the residual's cosine (207), top-k deflation keeps 80 % of its energy and removes only the hub — **§3** | — | — |

**Ranking of what is OPEN:** 1 (measured, real mechanism, pays on the CD4 residual) ≫ 6a (measured,
small) > 5a (worthless in October) > 2 (no mechanism for its sign, fails the filter).

---

## 2. Instrument and known-answer rows

162's `plolo.py` imported verbatim; 218's `resid218.py` / `compose218.py` imported for
`score_effects`, `self_residual`, `source_residual`, `union_rows`, `solve_fraction_carrying`,
`merged` and `shuffled`; fold caches sha-identical to 218's (`fold-k562_gwps.npz`
`c13c1cec…`). The shipped stack is `cd4@0.34` (production `solve_fraction`) + the H1 union at
`g = 0.28` (carrying-row rule) merged into **one** add with **one** global preserve (218 §7.2 =
222's builder, `repro222` max|Δ| 0). Endpoints as 218: `cos_agg` (energy-weighted aggregate), macro
cosine (`a23_cos`, per-target mean), `pds` (exact rank statistic, a23 render), q_w (`a15_sgn_w`), all
averaged over challenge frames A/B/C; `kolf2` apart throughout.

| # | row | against | result |
|---|---|---|---|
| KA-1 | `f0` (whitened base) ×4 folds | 218 `resid218` `f0` | **max|Δ| 0.0e+00** on `cos_agg`, macro, `pds` |
| KA-2 | `self0.1_bare` ×4 | 218 `self0.1` | **0.0e+00** |
| KA-3 | `cd4only` (f = 0.34) ×4 | 218 `cd40.34+none0` | **0.0e+00** |
| KA-4 | **`STACK`** ×4 | 218 `cd40.34+h1PAIRCARRY0.28` | **0.0e+00** → the mandated row **macro +7.6803 %, `pds` +0.004887** (3 cancer folds) reproduces to the digit |
| KA-5 | deflation mode `none`, dose 0 ×4 | `STACK` | **bit-identical arrays** (`np.array_equal`) |
| KA-6 | `sweep246b`/`c` re-derive `STACK`, `B_cd4_svd5`, `C_partial@g`, `self0.1_bare` in fresh processes | `sweep246` | **0.0e+00**, 4 folds each |
| KA-7 | `six246.py`: a15 f = 0.10 and `cd4f034` — `nmae` wald89, top-89 precision, delivered ratio | 216 `six216.json` | **0.0e+00** ×6 |
| KA-8 | `six246.py`: the stack — `nmae`, precision, delivered ratio, `pds` 4f / 3f | 222 `six222.json` (`h1union g=0.28`: 0.98569 / 0.60312 / 1.0033 / 0.75065 / 0.78655) | **0.0e+00** ×5 |

Fold sizes: 269 / 290 / 294 / 273 targets; CD4 complete rows 251 / 267 / 271 / 253; partial-condition
rows 18 / 23 / 23 / 20; rows CD4 does not touch at all 5 / 8 / 8 / 3 (these are the "non-carrying" rows
for family B).

---

## 3. Phase 2 — measurement

All numbers: **three cancer folds vs `STACK`**, `kolf2` apart. `E rm` = fraction of the residual's
energy the projection removes. Price = 218 §7.4's rate card on the macro cosine + PLOLO `pds` (×2.0,
`price218.price`, the headline convention of 218/222); §4 gives all three transfer factors.

### 3.1 Family A — hypothesis 1 as posed: the deflated self-residual ON TOP of the shipped stack

| arm (dose 0.10 unless noted) | rel agg | rel macro | Δ`pds` | Δq_w | +cos / +pds folds | E rm | price |
|---|---:|---:|---:|---:|---|---:|---:|
| `A_none` (218's self-residual, no projection) | +1.60 % | +2.17 % | **−0.00390** | +0.00006 | 2/3 · 0/3 | 0 | −0.0013 |
| `A_centre` (per-gene centring only) | +1.55 % | +2.08 % | −0.00383 | +0.00001 | 2/3 · 0/3 | 0.004 | −0.0014 |
| `A_svd1` | −0.74 % | +0.22 % | +0.00109 | −0.00014 | 2/3 · 2/3 | 0.157 | +0.0010 |
| **`A_svd2`** | −0.30 % | **+0.23 %** | **+0.00167** | −0.00013 | 2/3 · 2/3 | 0.187 | **+0.0014** |
| `A_svd5` | −0.48 % | +0.02 % | +0.00155 | −0.00011 | 1/3 · 2/3 | 0.245 | +0.0012 |
| `A_svd10` | −0.52 % | −0.14 % | +0.00190 | −0.00005 | 1/3 · 2/3 | 0.315 | +0.0013 |
| `A_svd20` | −0.58 % | −0.33 % | +0.00196 | +0.00001 | 1/3 · 2/3 | 0.415 | +0.0012 |
| `A_op2` / `A_op10` (operator's own top-k eigenvectors, no centring) | −0.42 / −1.06 % | +0.44 / −0.21 % | +0.00008 / +0.00049 | | | 0.148 / 0.179 | +0.0004 / +0.0002 |
| `A_rand1…20` (random k-dim subspace, 1–8 % energy) | +1.3 … +1.6 % | +1.7 … +2.1 % | **−0.0036 … −0.0050** | | 2/3 · 0/3 | | −0.0011 … −0.0024 |
| `A_*_SHUF` (target-deranged) | −0.65 … −1.21 % | −0.57 … −0.99 % | −0.0003 … −0.0111 | | 0/3 | | negative |
| dose 0.20: `A_svd2@0.2` | −0.85 % | +0.22 % | +0.00257 | −0.00031 | | | +0.0021 (macro) / +0.0013 (agg) |

**Reading.** The projection does what the hypothesis says to `pds` (−0.0039 → +0.0017, and the random
controls show it is the *direction*, not the energy: they stay at −0.004), but the self-residual has
**no cosine left to give on top of a stack that already carries CD4 @ 0.34 and the H1 union** — its
aggregate cosine is negative at every k ≥ 1 and its macro cosine is ≤ +0.23 %. `kolf2`: negative on
cosine at every k (−1.2 … −3.6 %). **NO-GO as posed**; the best price, +0.0014 … +0.0021, sits at or
under the +0.0019 gate with the aggregate endpoint negative.

### 3.2 Family D — the mechanism, on the BARE whitened base (no CD4, no H1), vs `f0`

| arm (dose 0.10) | rel agg | rel macro | Δ`pds` | Δq_w | E rm |
|---|---:|---:|---:|---:|---:|
| `D_none` = 218's `self0.1` | **+5.35 %** | +5.03 % | **−0.00969** | +0.00033 | 0 |
| `D_centre` | +5.26 % | +4.88 % | −0.00975 | +0.00032 | 0.004 |
| `D_svd1` | +1.97 % | +2.43 % | **+0.00113** | +0.00015 | 0.157 |
| `D_svd2` | +2.54 % | +2.45 % | **+0.00174** | +0.00008 | 0.187 |
| **`D_svd5`** | **+2.18 %** | +2.09 % | **+0.00185** | +0.00003 | 0.245 |
| `D_svd10` / `D_svd20` | +1.99 / +1.72 % | +1.76 / +1.37 % | +0.00176 / +0.00157 | | 0.315 / 0.415 |
| **`D_randE5` / `b`** (random, **energy-matched** to svd5: 63 dims, 25 %) | **+5.23 / +4.56 %** | +4.29 / +4.43 % | **−0.00987 / −0.01077** | | 0.250 / 0.249 |
| `D_svd5_SHUF` (deranged rows) | −0.82 % | −0.86 % | −0.00097 | | |
| `kolf2` apart: `D_none` +2.38 % / −0.0021; `D_svd5` +1.86 % / **+0.0005**; `D_randE5` +1.27 % / −0.0029 | | | | | |

**This is the mandate's hypothesis, confirmed and bounded.** The `pds` cost of the un-whitened
self-residual lives entirely in its top-1 … top-5 target-shared directions (removing 16–25 % of its
energy there flips −0.0097 to +0.0011 … +0.0019); removing the *same* energy in random directions
changes nothing. About **41 % of the cosine gain** survives the projection (the rest was the shared
programme, which the truth also contains — 153 §c). Per-gene centring alone does nothing (0.4 % of the
energy; the residual is already centred on the bank-wide mean, 204 §3.2). **The operator's own top-k
eigenvectors are the wrong basis** (`A_op*`): the residual's shared subspace is its own, not the
donors'.

### 3.3 Family B/B2 — the same projection applied to the SHIPPED CD4 residual (fixed shipped fraction)

| arm | rel agg | rel macro | Δ`pds` | +pds folds | Δq_w | E rm | price (×2 pds) |
|---|---:|---:|---:|---|---:|---:|---:|
| `B_cd4_centre` | −0.05 % | −0.10 % | +0.00008 | 2/3 | −0.00012 | 0.001 | 0.0000 |
| `B2_cd4_svd1` | −0.20 % | −1.04 % | +0.00653 | 3/3 | −0.00009 | 0.079 | +0.0041 |
| `B2_cd4_svd2` | −0.47 % | −1.09 % | +0.00759 | 3/3 | +0.00077 | 0.118 | +0.0049 |
| `B2_cd4_svd3` | −0.75 % | −1.22 % | +0.00905 | 3/3 | +0.00086 | 0.151 | +0.0059 |
| `B2_cd4_svd4` | −1.19 % | −1.38 % | +0.00902 | 3/3 | +0.00081 | 0.182 | +0.0057 |
| **`B2_cd4_svd5`** | **−1.19 %** | **−1.29 %** | **+0.01227** | **3/3** | **+0.00058** | **0.202** | **+0.0082** |
| `B2_cd4_svd7` | −1.93 % | −1.83 % | +0.01132 | 3/3 | +0.00051 | 0.240 | +0.0071 |
| `B2_cd4_svd10` | −2.72 % | −2.45 % | +0.01243 | 3/3 | +0.00006 | 0.286 | +0.0075 |
| `B2_cd4_svd15` / `svd20` | −4.93 / −7.01 % | −3.70 / −5.17 % | +0.01295 / +0.01297 | 3/3 | −0.00036 / −0.00052 | 0.344 / 0.388 | +0.0070 / +0.0059 |
| `B2_cd4_svd3_resolved` (CD4 fraction re-solved on the deflated object) | −0.17 % | −0.87 % | +0.00899 | 3/3 | +0.00102 | 0.151 | +0.0061 |
| **`B2_cd4_svd5_resolved`** (fraction 0.374 → **0.396**) | **−0.46 %** | **−0.81 %** | **+0.01226** | **3/3** | +0.00064 | 0.202 | **+0.0086** |
| `B2_cd4_svd10_resolved` | −1.54 % | −1.67 % | +0.01281 | 3/3 | +0.00002 | 0.286 | +0.0083 |
| `B_h1_svd1/2/5` (the H1 union deflated instead) | −1.0 / −1.4 / −3.4 % | −0.7 / −0.9 / −2.1 % | +0.0009 / +0.0010 / +0.0006 | | | | +0.0001 / +0.0001 / −0.0011 |
| `B_both_svd5` | −4.60 % | −3.40 % | +0.01335 | 3/3 | +0.00052 | | +0.0075 |

Per fold, `B2_cd4_svd5`: k562 −2.45 % / −2.33 % / **+0.00937**; hct116 −2.70 % / −2.42 % /
**+0.00947**; hek293t **+1.57 % / +0.89 % / +0.01797** (the arm *gains* cosine on the fold whose
donors are K562 + HCT116). `kolf2` apart: −4.06 % / **+0.00378** (svd3 +0.0017, svd10 +0.0042).

The `pds` gain rises to k ≈ 5 and then saturates (+0.0123 … +0.0130 out to k = 20) while the cosine
cost keeps growing; **k = 3–5 is the knee**. The H1 union has no such hub (77 rows, 7 shared): deflating
it costs cosine for nothing, consistent with `h1union` having *raised* board `pds`.

### 3.4 Falsifiers, all mandatory, on `B2_cd4_svd5` (and the k it neighbours)

| falsifier | result |
|---|---|
| **energy-matched random subspace** (a random k′-dim subspace of the centred row space with k′ chosen so it removes the SAME energy as svd-k; two seeds each) — 209's control, matched | k = 3 (38–39 dims, 15–16 %): `pds` **−0.0009 / −0.0044**, cosine −2.2 / −3.5 %; **k = 5 (51–52 dims, 21 %): `pds` −0.0001 / −0.0050, cosine −3.7 / −4.1 %**; k = 10 (73 dims, 29 %): −0.0032 / −0.0030, −4.6 / −4.7 %. **Removing the same energy in random directions buys no `pds` and costs 3× the cosine.** |
| **shuffled projector** (top-k fitted on a row-deranged copy, applied to the real rows) | k = 3 / 5 / 10: `pds` −0.0007 / +0.0004 / +0.0013, cosine −0.7 / −0.6 / −0.3 % — a k-dim subspace of the same rank fitted to the wrong rows removes 1–3 % of the energy and nothing else |
| **small random subspace** (`B_cd4_rand5/b`, 2 % energy) | `pds` +0.0005 / −0.0005 |
| **target-shuffled residual "at every dose"** — for family B the dose is k and the *residual* is the shipped CD4 (shuffling it destroys the arm itself, 207 F1); the row-shuffle is applied to the *projector* above. For family A the shuffled residual is negative at every (k, f): §3.1 | pass |
| **row-matched control** — not applicable (no row set changes: every carried row is deflated, no row gains a residual). Applied where it is defined: family C, §3.7 | — |
| **non-carrying rows must not move** (5 / 8 / 8 rows CD4 does not touch) | worst \|rel Δ mean cosine\| **0.018 pp** at k = 5 (0.026 pp at k = 10) — the preserve scalar only |
| **leave-top-3-out retention** (per-target Δ from `score_effects`, each endpoint on its own ordering) | `pds`: full +0.01227, leave-top-1-out +0.00999 (**81 %**), leave-top-3-out +0.00685 (**56 %**); cosine: the loss is *broad* (LT3O −0.00081 vs full −0.00049). `h1union` retained 14.6 % of its Δ`pds` under the same statistic (222 §9) and realised 70 % of it on the board |
| **split-half projector transfer** (basis fitted on a random half of CD4's carried rows, applied to the other half, both ways; two seeds) | k = 2: +0.0061 / +0.0068; k = 3: +0.0068 / +0.0069; **k = 5: +0.0108 / +0.0104**; k = 10: +0.0118 / +0.0114 — at −0.1 … −0.8 % aggregate cosine (less than in-sample, because a half-fitted basis removes 9–12 % of the energy instead of 20–29 %). **The subspace is a property of the CD4 source, not of the rows it is applied to.** `kolf2`: +0.0007 … +0.0028, all positive |
| **leave-one-fold-out k selection** (k chosen on two folds by the rate-card price, read on the third) | fixed: k = 5, 5, 10 → **+0.0015 / +0.0028 / +0.0069** at `pds` transfer 0.7 / 1.0 / 2.0×; resolved: +0.0017 / +0.0029 / +0.0070; positive on 3/3 folds at every factor |
| **`kolf2` reported apart** | `pds` +0.0017 (k = 3) … +0.0049 (k = 10 resolved), cosine −2.5 … −5.8 %: same sign, larger cosine cost on the pluripotent recipient |
| **non-vacuity** | every residual's carrying-row median norm > 0 asserted in-loop; solved fractions 0.374 / 0.371 / 0.376 (CD4), 0.224 / 0.223 / 0.225 (union) per fold |

### 3.5 The exchange rate — `pds` bought per unit of aggregate cosine given up (same folds, same base)

| axis | arm | rel agg | Δ`pds` | Δ`pds` per −1 % agg cosine |
|---|---|---:|---:|---:|
| **dose** (216's lever) | `DOSE_cd4@0.28` | −2.37 % | +0.00221 | **0.00093** |
| dose | `DOSE_cd4@0.20` | −6.19 % | +0.00401 | 0.00065 |
| **deflation** | `B2_cd4_svd3` | −0.75 % | +0.00905 | **0.0122** |
| deflation | `B2_cd4_svd5` | −1.19 % | +0.01227 | **0.0103** |
| deflation | `B2_cd4_svd10` | −2.72 % | +0.01243 | 0.0046 |
| deflation, re-solved | `B2_cd4_svd5_resolved` | −0.46 % | +0.01226 | **0.0266** |
| random, energy-matched | `B2_cd4_randE5/b` | −3.7 / −4.1 % | −0.0001 / −0.0050 | ≤ 0 |

**11× (k = 5) to 30× (k = 3) the dose axis; ~29× for the re-solved k = 5.** 212's gamma axis has
d(overall `pds`)/d(cos) = −0.3517 against d(overall `mse`)/d(cos) = +0.3655 — they cancel. Here the
`pds` leg is +0.0046 overall (×1) against a `mse` leg of −0.0004.

### 3.6 Hypothesis 5b re-measured on the stack (the CD4 dose)

`DOSE_cd4@0.28`: −2.37 % / −1.63 % macro / `pds` +0.0022, price **−0.0006 … +0.0005**;
`DOSE_cd4@0.20`: −6.19 % / −4.57 % / +0.0040, price **−0.0003 … −0.0015**. 216 §6.0's flat top
reproduces in direction; the dose axis is closed here too.

### 3.7 Family C — CD4's partial-condition rows (hypothesis 6a), carrying-row sizing, on the stack

| arm | rel agg | rel macro | Δ`pds` | folds +cos/+pds | price (×2) |
|---|---:|---:|---:|---|---:|
| `C_partial@0.17` | +0.65 % | +0.43 % | +0.00112 | 3/3 · 3/3 | +0.0011 |
| **`C_partial@0.34`** | **+1.16 %** | **+0.59 %** | **+0.00178** | 3/3 · 3/3 | **+0.0018** |
| `C_partial_eq@0.34` (rows norm-equalised to the complete rows) | +0.94 % | +0.60 % | +0.00132 | 3/3 · 3/3 | +0.0014 |
| `C_partial_SHUF@0.34` (deranged among the partial rows) | −0.23 % | −0.16 % | −0.00021 | | −0.0003 |
| **`C2_partial_ROWSELF@0.34`** — the un-whitened consensus on the SAME rows and gene support (the row-matched control) | **+1.36 %** | **+0.57 %** | **−0.00012** | 2/3 · 2/3 | +0.0003 |
| `kolf2`: `C_partial@0.34` +1.07 % / +0.0003; row-matched −0.05 % / −0.0025 | | | | | |

The **cosine** is fully reproduced by the row-matched self control (so it is "more un-whitened residual
on rows that had none", not CD4 information); the **`pds`** is CD4-specific (+0.0018 vs −0.0001) but
concentrated — leave-top-3-out retains **−6.6 %** (top-1 alone carries 50 %). Exact `nmae`
+0.00005, `e` = 0.75. **NOT build-worthy alone** (price +0.0009 … +0.0018 across transfer factors,
under or at the gate, concentrated). Composed on the deflation (`E_cd4_svd5+partial@0.34`): aggregate
−0.02 %, macro −0.70 %, **`pds` +0.01405 (3/3), LT3O retains 61.5 %**, price **+0.0032 / +0.0047 /
+0.0100** — additive to within 0.0002 of the two parts.

---

## 4. The six-member price — every member, never at zero (§8.4)

Rate card: `cd4f034`'s realised board row per +0.01 absolute cosine (218 §7.4: `mse` +0.003634,
`reach` +0.002446, `nmae` +0.000980, `fid` +0.000847, `jac` +0.000100), applied to the **macro** PLOLO
cosine delta at k_cos = 1 (the PLOLO → board factor is **not identified**: 1.000 on `h1union`, 1.384
on `cd4f034`; a factor > 1 makes the cosine legs *more* negative). `pds` from PLOLO (raw Δ / 0.452 /
6 × 1.011) at **three transfer factors**: 0.7× (`h1union` realised 70 % of its predicted `pds` gain,
`docs/NEXT-BUILD-SPEC.md` §4), 1.0×, 2.0× (§8.4's realised factor on `cd4f034`'s `pds` *cost*, k_pds
= 1.987 — what `price218.price` charges). The exact `nmae` (six246, 188's identity, wald89 pool, 4
folds) is a **magnitude bound, not a sign** (222 §8: it mis-signed its own family's board point).

| arm | macro | Δ`pds` (σ) | `mse` | `reach` | `nmae` rate-card / **exact** | `fid` | `jac` | `pds` 0.7× / 1.0× / 2.0× | **TOTAL 0.7× / 1.0× / 2.0×** | Δq_w (232 rate) | `e` | dlv/a15 | Δ prec89 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **`B2_cd4_svd5`** | −1.286 % | +0.01227 (8.2) | −0.00042 | −0.00028 | −0.00011 / **+0.00013** | −0.00010 | −0.00001 | +0.00320 / +0.00457 / +0.00915 | **+0.00227 / +0.00365 / +0.00822** | +0.00058 (+0.00023) | 0.250 | 1.0011 | −0.00086 |
| **`B2_cd4_svd5_resolved`** | −0.809 % | +0.01226 (8.2) | −0.00027 | −0.00018 | −0.00007 / +0.00031 | −0.00006 | −0.00001 | +0.00320 / +0.00457 / +0.00914 | **+0.00262 / +0.00399 / +0.00856** | +0.00064 (+0.00025) | 0.156 | 1.0002 | −0.00047 |
| `B2_cd4_svd3` | −1.220 % | +0.00905 (6.1) | −0.00040 | −0.00027 | −0.00011 / +0.00033 | −0.00009 | −0.00001 | +0.00236 / +0.00337 / +0.00675 | +0.00148 / +0.00249 / +0.00587 | +0.00086 | 0.136 | 1.0024 | −0.00097 |
| `B2_cd4_svd3_resolved` | −0.867 % | +0.00899 (6.0) | −0.00028 | −0.00019 | −0.00008 / +0.00029 | −0.00007 | −0.00001 | +0.00235 / +0.00335 / +0.00671 | +0.00172 / +0.00273 / +0.00608 | +0.00102 | 0.167 | 0.9994 | −0.00065 |
| `B2_cd4_svd10` | −2.447 % | +0.01243 (8.3) | −0.00080 | −0.00054 | −0.00022 / +0.00006 | −0.00019 | −0.00002 | +0.00324 / +0.00463 / +0.00927 | +0.00148 / +0.00287 / +0.00750 | +0.00006 | 0.405 | 0.9967 | −0.00133 |
| `B3_cd4_svd5_splithalf` (out-of-row basis) | −0.818 % | +0.01083 (7.3) | −0.00027 | −0.00018 | −0.00007 | −0.00006 | −0.00001 | +0.00283 / +0.00404 / +0.00808 | +0.00224 / +0.00345 / +0.00749 | +0.00051 | | | |
| `E_cd4_svd5+partial@0.34` | −0.699 % | +0.01405 (9.4) | −0.00023 | −0.00015 | −0.00006 / +0.00017 | −0.00005 | −0.00001 | +0.00367 / +0.00524 / +0.01047 | +0.00316 / +0.00473 / +0.00997 | +0.00040 | 0.113 | 1.0015 | −0.00103 |
| `C2_partial@0.34` alone | +0.587 % | +0.00178 (1.2) | +0.00019 | +0.00013 | +0.00005 / +0.00005 | +0.00004 | +0.00001 | +0.00046 / +0.00066 / +0.00132 | +0.00089 / +0.00109 / +0.00175 | −0.00017 | 0.750 | 1.0061 | −0.00020 |
| `A_svd2@0.1` (hypothesis 1 as posed) | +0.232 % | +0.00167 (1.1) | +0.00008 | +0.00005 | +0.00002 | +0.00002 | 0 | +0.00044 / +0.00062 / +0.00124 | +0.00060 / +0.00079 / +0.00141 | −0.00013 | | | |
| `DOSE_cd4@0.28` (the dose axis) | −1.626 % | +0.00221 (1.5) | −0.00053 | −0.00036 | −0.00014 | −0.00012 | −0.00001 | +0.00058 / +0.00082 / +0.00164 | −0.00060 / −0.00035 / +0.00047 | −0.00013 | | | |
| `B2_cd4_randE5` / `b` (energy-matched controls) | −2.59 / −3.35 % | −0.0001 / −0.0050 | −0.0009 / −0.0011 | | | | | ≤ 0 | **−0.0019 / −0.0043** (×1) | | | | |
| leave-one-fold-out, fixed k | | | | | | | | | **+0.00154 / +0.00276 / +0.00685** | | | | |
| leave-one-fold-out, re-solved | | | | | | | | | +0.00168 / +0.00290 / +0.00698 | | | | |

Notes on the members. **`mse`** is the only member that reads the cosine loss directly (197's identity;
Δraw ≈ +0.0023 at k = 5, i.e. −0.0004 overall — the rate card's number). **`nmae`**: the rate card
charges the cosine loss; the exact identity says the removed component was **unmatched** magnitude
(`e` = 0.25 < ½, U/D falls 4× more than M/D) so `nmae` *improves* +0.00013; with 222 §8's caveat the
honest statement is **|Δ`nmae`| ≤ 0.0004**. **`reach`/`fid`** are priced through the rate card
(cosine-driven) and are the members that would recover if the q_w rise (+0.00058, worth +0.00023 on
232's rate) is what the board reads instead; the ordering half readable here (top-89 precision) moves
**−0.0009**, i.e. 1/6 of what `h1union` gained. **`jac`** is pinned near zero. **Delivered
amplitude** 1.0011 (vs the stack's 1.0033 and `cd4f034`'s 1.0009) — the arm is amplitude-neutral to
0.2 %, in the profile class rule 13 keys as shape-only.

---

## 5. BUILD-WORTHY / NOT — one line each

- **`B2_cd4_svd5` (or `svd5_resolved`): BUILD-WORTHY, shape-only class, conditional** — priced
  **+0.0023 / +0.0037 / +0.0082** at `pds` transfer 0.7 / 1.0 / 2.0× against the **+0.0019** gate
  (1.2 / 1.9 / 4.3×); leave-one-fold-out +0.0015 / +0.0028 / +0.0069; **below the +0.0089
  amplitude-changing gate at every factor** (0.26 … 0.92×). `pds` ≥ 0 on 3/3 cancer folds and on
  `kolf2`; every falsifier passes. **The condition:** the whole margin at the pessimistic end is the
  board's `pds` transfer for a *gain*, which has exactly one calibration point (`h1union`, 0.70×, on
  a gain that was 99.8 % top-3-concentrated where this one is 44 %); at 0.7× the margin over the gate
  (+0.0004) is the size of PLOLO's own `pds` LOO σ in overall units. **It is also the first
  candidate that lowers cosine** (−1.2 % aggregate ≈ −0.0011 board cosine, ~5 % of `cd4f034`'s
  gain) to buy retrieval — the opposite trade to every shipped direction arm — so a board point is
  the only instrument that can price it, and the raw `pds` gap to the leader (0.836 vs 0.870) is
  exactly the member it moves: +0.0123 PLOLO → **+0.009 … +0.025 raw board `pds`** at 0.7 … 2.0×.
- **Hypothesis 1 as posed (deflated self-residual on the stack): NOT** — +0.0006 … +0.0021, aggregate
  cosine negative.
- **Hypothesis 2 (per-target dose): NOT** — unmeasured, fails the filter, no mechanism for its sign.
- **Hypothesis 3 (gating beyond 10 CPM): NOT** — closed (088b, 213, 164, 006).
- **Hypothesis 4 (sign-consistent consensus): NOT** — closed (232 `maj_strict`, 234); vacuous at arity 3.
- **Hypothesis 5 (per-source dose / joint grid): NOT** — ≤ +0.0008, and `h1_2025` evaporates in October.
- **Hypothesis 6a (partial-condition CD4 rows): NOT alone** — +0.0009 … +0.0018, concentrated; as a
  rider on the deflation it adds +0.0009 … +0.0018 with 61.5 % LT3O retention (`E_cd4_svd5+partial`:
  +0.0032 / +0.0047 / +0.0100).

**If it is built (not proposed here — no slot, no build):** one variable, `VCC_CD4_DEFLATE_K=5`
(per-gene centring over CD4's complete rows, top-5 left-singular projector `(I − U₅U₅ᵀ)` on the
centred residual, applied before `solve_cd4_fraction`/the merge — 12 lines in `lib246.deflate`);
known-answer rows: k unset → **byte-identical to `h1union`** (`f13ce40b…`), and the fold-level
`pds` +0.01227 / macro −1.286 % must reproduce through the builder as 222's `repro222` did for 218.
**October:** the projector is fitted on CD4's rows for the *panel* (like the operator, 8.3) — refit
on the D/E/F panel's CD4 rows (2 s), and the split-half result (§3.4) says the subspace transfers
across row sets; record it in `DEF-RUNBOOK` if the arm is ever shipped. The H1 union must **not** be
deflated (§3.3, `B_h1_*`).

---

## 6. Honest limits

- **One instrument, one class of truth.** All of this is 162's PLOLO with donor-library pseudobulk
  truths; 222 realised 70 % of its predicted `pds` gain and `cd4f034` 199 % of its predicted cost.
  The direction of `pds` is what these folds can say; its board magnitude is a range.
- **The cosine cost is asymmetric across folds** (+1.6 % on `orion_hek293t`, −2.5 / −2.7 % on the
  other two). The k562/hct116 folds have the Orion pair or K562+HEK as donors; a reader who weights the
  worst fold gets price +0.0006 … +0.0053 (k = 5, 0.7 … 2.0×).
- **The removed subspace is fitted on the panel's own CD4 rows** (unsupervised; no truth enters) and
  transfers across halves of them (§3.4), but whether the *same* 5 directions are the right ones for a
  disjoint October panel is not measured here — only that a refit on a different row set behaves the same.
- **The rate card rests on one board point** (218 §7.4's caveat) and prices a cosine *loss* with a
  card fitted on a gain; a `k_cos` of 1.384 would make the cosine legs −0.0013 instead of −0.0009.
- **No `cell-eval2`, no build, no h1-test** (oracle-leaked for anything on top of `h1union`, 222 §11).
- **Family A's "row-matched control" is not defined** for a residual that carries every CD4 row; the
  energy-matched random subspace (209's control, matched) plays that role and is reported at every k.
- The family-C row-matched control in `sweep246.py` was vacuous (the self-residual is masked to CD4's
  complete-row support); it was rebuilt on the partial rows' own support in `sweep246b.py`
  (`C2_partial_ROWSELF`), and that is the control quoted.

---

## 7. Files

| script | what it produces |
|---|---|
| `lib246.py` | staging (218/222's stack, verbatim imports), `deflate()` (centre / svd-k / rand-k / op-k), `merged_many()`, `cd4_partial()` |
| `sweep246.py` | per fold: KA-1…5, families A / B / C (78 arms) → `work/sweep246-<fold>.json`, `work/pertarget246-<fold>.npz`, logs with `/usr/bin/time -l` |
| `sweep246b.py` | B2 k-sweep, re-solved dose, energy-matched and shuffled-projector controls, D (bare-base mechanism), C2 (correct row-matched control), E (composed) (40 arms/fold) |
| `sweep246c.py` | B3 split-half projector transfer, the CD4 dose axis (11 arms/fold) |
| `analyse246.py` | cancer means, `kolf2` apart, concentration, non-carrying rows, LOFO k selection → `work/analyse246.json`, `work/analyse246.log` |
| `mech246.py` | D vs `f0`, exchange rates → `work/mech246.json` |
| `six246.py` | 216/222's exact six-member route on the surviving arms + KA-7/8 → `work/six246.json`, `work/six246.log` |
| `price246.py` | §4's table → `work/price246.json`, `work/price246.log` |

```bash
cd experiments/246-direction-hypotheses
for f in k562_gwps orion_hct116 orion_hek293t kolf2; do
  /usr/bin/time -l uv run python sweep246.py  --fold $f   # ~60 s, 4.9 GiB
  /usr/bin/time -l uv run python sweep246b.py --fold $f   # ~50 s
  /usr/bin/time -l uv run python sweep246c.py --fold $f   # ~15 s
done
uv run python analyse246.py && uv run python mech246.py
/usr/bin/time -l uv run python six246.py                  # 26 s, 8.5 GiB
uv run python price246.py
```
