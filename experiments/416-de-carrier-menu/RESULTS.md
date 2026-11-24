# 416 — DE-carrier menu under ballast (H1, seed-replicated): unwhitened carriers (γ 0 + CD4) are the one change that clears noise, +0.0034 predicted overall

Status: **COMPLETE** (2026-09-22 20:11–23:17Z, worker, Opus 5.5). 55 serial H1 arm runs, 54 scored (`prior_off 0` was re-run after a NaN-in-JSON write bug). No builder edits, no
`vcc`, no `cell-eval2`; everything written inside `experiments/416-de-carrier-menu/`.

**Answer.** Of 12 carrier designs, only **γ 0 carriers** (the production composition with the whitening step skipped:
unwhitened donor consensus + the frozen CD4 residual, rows norm-matched to the shipped rows) beats the shipped carriers by
more than 2 sd of seed noise. With the shipped 410b bulk (the production construction) it reads **nmae −0.0067, fid
+0.0020, jac +0.00027 raw, reach +0.0014 ± noise, pds/mse unchanged within noise**. At board slopes that is **+0.0034
overall**, of which +0.0031 comes from the three low-noise members. **Recommended board arm: 410b (or the current bulk
incumbent) with γ 0 carriers, one variable.** The mechanism is additive and measured: CD4 is worth +0.006 of DE and
un-whitening +0.003, at either level of the other. This **corrects 410 §7** ("whitened carriers stay"). That test
dropped CD4 together with the whitening, and the CD4 loss hid the whitening gain.

## Harness
`run416.py` is 410's `run410b.py` (ballast, ratio mode, pair-jackknife boost) copied, with one variable, the carrier
design, and an emission-seed replicate (`b.make_rngs(seed)`, seeds 0, 1, 2, and 3, 4 for base and g000). The jitter and
jk-boost streams are fixed as in 410. The bulk effect is the shipped cap15 matrix for every design (`--bulk base`, per the
mandate), and 410b's key-scaled bulk for the confirmation (`--bulk ps`). The carrier effect sets the profile, the emission
spread, and the rank axis and key of the carrier block. The sign prior is CLEAN4 q80 unless the design changes it.
Carrier matrices come from `prep416.py` / `prep416b.py`: 323's `stage.npz` cap15 pre-whitening design → the H1 operator at
γ → 323's frozen CD4 residual (dose re-solved against that base, as production does) → each row norm-matched to the
shipped cap15 row. Shas are in `work/carriers416.json` and `work/carriers416b.json`.

### Known answers (every one passes; the driver aborts on failure)
| row | reproduces | result |
|---|---|---|
| `ka416.py base395` | 410 `metrics-base395.json` (400's harness, no ballast) | emitted-count, panel-total and statistics digests equal; members ≤ 1.1e-16 |
| `run416.py base 0 --bulk ps` | 410 `bal_lo03_hi3_ratio_jk` (the shipped 410b construction) | all stats arrays (sums, jk, p, q, lfc) bit-equal; members Δ 0 |
| `run416.py x125 0` | 410 `de_x1p25_ratio_jk` (a carrier change with K 1–3 and a live jk boost) | arrays bit-equal |
| `run416.py pre 0` | 410 `de_unwhite_ratio_jk` | arrays bit-equal |
| `prep416.py` | γ 0.65 + CD4 composed here == 323's cap15 (`565d57ab…`) | sha equal |

410's `bal_ident_ratio` record does not reproduce with 410's committed `run410b.py`. It was made with 999,000-count
ballast cells and no 995k trim; the committed code trims 4,000 counts per target. Details in
`work/h1/known-answer-base_s0-vs-bal_ident_ratio-FAILED-trim.json`. That is why `base_s0` differs from it by pds 1.4e-4.

## Noise floor (bulk = cap15, seeds 0–2)
Pooled within-design single-seed sd (26 dof): pds 0.00060, mse 0.00005, **nmae 0.00008, fid 0.00008, jac 0.00001**,
**reach 0.0063**. **DE-sum at board slopes: 0.00119 overall**, which is almost all reach. Base alone (n = 5): reach sd 0.0061,
DE-sum 0.0011. Reach realisation noise here (0.006–0.008) is about half the ±0.013 in the mandate, because that figure was
a 400 → 399-carrier change, which reshuffles more than a seed does. The three non-reach DE members have realisation noise
near 1e-4 raw. Any design's nmae/fid/jac movement above ~3e-4 is a real property of that design on H1.

## The menu (bulk = cap15, seeds 0–2; Δ vs base, ± = the design's own seed sd; overall = Σ slope·Δ/6)
| design | what | Δpds | Δmse | Δnmae | Δfid | Δreach | Δjac | **DE-sum ovr** | bulk ovr | total ovr | paired DE sd | > 2 sd (0.0024)? |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| base | shipped carriers (γ 0.65 + CD4, q80) | 0.88950±0.00028 | 0.97502±0.00002 | 0.95212±0.00000 | 0.50902±0.00007 | 0.27495±0.00835 | 0.07474±0.00001 | 0 | 0 | 0 | — | — |
| **g000** | (b) γ 0 (unwhitened + CD4) | −0.00112±0.00000 | −0.00073±0.00008 | **−0.00683±0.00004** | **+0.00253±0.00009** | −0.00212±0.00245 | +0.00025±0.00001 | **+0.00301** | −0.00029 | +0.00273 | 0.00146 | **YES** |
| g030 | (b) γ 0.3 | −0.00070±0.00085 | −0.00011±0.00003 | −0.00515±0.00009 | +0.00149±0.00003 | −0.00460±0.00952 | +0.00019±0.00001 | +0.00147 | −0.00024 | +0.00123 | 0.00035 | no |
| q70 | (d) protect q 0.7 | 0 | 0 | −0.00059±0.00005 | −0.00079±0.00010 | +0.00775±0.01060 | −0.00002 | +0.00117 | 0 | +0.00117 | 0.00045 | no (reach noise only; non-reach −0.0003) |
| signrep | (e2) 387's sign repair, carriers only | −0.00121±0.00090 | −0.00062±0.00004 | −0.00001±0.00015 | +0.00035±0.00004 | +0.00163±0.00282 | −0.00000 | +0.00050 | −0.00034 | +0.00016 | 0.00167 | no |
| q90 | (d) protect q 0.9 | 0 | 0 | +0.00011±0.00003 | +0.00011±0.00003 | −0.00246±0.00951 | +0.00002 | −0.00042 | 0 | −0.00042 | 0.00023 | no |
| pre | decomposition: γ 0, **no CD4** (= 410 de_unwhite) | −0.00131 | −0.00095 | +0.00281 | +0.00102 | −0.01476 | −0.00028 | −0.00310 | −0.00032 | −0.00342 | 0.00155 | no |
| g090 | (b) γ 0.9 | −0.00135 | −0.00085 | +0.00782 | −0.00072 | −0.01175 | +0.00005 | −0.00475 | −0.00035 | −0.00510 | 0.00101 | no |
| w065nocd4 | decomposition: γ 0.65, **no CD4** | −0.00145 | −0.00053 | +0.00927 | −0.00219 | −0.00935 | −0.00036 | −0.00571 | −0.00044 | −0.00615 | 0.00213 | no |
| conf4 | (e1) prior may flip only CLEAN4-unanimous genes | 0 | 0 | −0.00137 | −0.00336 | −0.02296 | +0.00001 | −0.00582 | 0 | −0.00582 | 0.00077 | no |
| gm030 | (b+) γ −0.3 (anti-whitening) | −0.00126 | −0.00169 | +0.00033 | +0.00298 | **−0.04702** | −0.00014 | −0.00733 | −0.00017 | −0.00750 | 0.00198 | no |
| prior_off | (d) CLEAN4 sign prior off (K-class, unsteered) | 0 | 0 | −0.00124 | **−0.01147** | −0.02341 | +0.00002 | −0.01049 | 0 | −0.01049 | 0.00018 | no |
| gain | (c) 228's per-gene gain on the carriers | −0.00117 | −0.00045 | **+0.04091** | −0.00076 | −0.03287 | +0.00013 | −0.01782 | −0.00035 | −0.01817 | 0.00230 | no |

g000 with 5 seeds: DE-sum **+0.00314** (paired sd 0.00109), nmae −0.00687, fid +0.00253, jac +0.00025, reach −0.0015 ±
0.0025. Table source: `work/analysis416-seeds012.{json,md}` (`analyse416.py --seeds 0,1,2`); the 5-seed version is in
`work/analysis416.md`.

## Confirmation under the production bulk (410b's key-scaled bulk, `--bulk ps`, seeds 0–2)
| | pds | mse | nmae | fid | reach | jac | DE-sum ovr | bulk ovr | **total ovr** |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| base (= 410b construction; s0 is the bit-exact KA) | 0.89486±0.00035 | 0.92529±0.00008 | 0.95232±0.00004 | 0.50920±0.00003 | 0.27382±0.00402 | 0.07480±0.00000 | 0 | 0 | 0 |
| g000 carriers, Δ | −0.00033±0.00162 | −0.00040±0.00003 | **−0.00671** | **+0.00204** | +0.00136±0.00689 | +0.00027 | **+0.00337** | −0.00005 | **+0.00332** |

Under the production bulk the small pds leak seen with the cap15 bulk disappears within noise, and mse improves slightly.
The ballast falsifier holds. K distribution: base has K 1–16 (21 targets at K 1); g000 has K 2–16. Clip mass is
negligible (551 vs 634 counts in total).

## Controls and mechanism (`analyse416k.py` → `work/analysis416k.json`; seed-paired, seeds 0–2)
- **K is not the effect.** Base carriers forced to g000's per-target ballast count K: nmae/fid/jac |Δ| ≤ 1e-5, reach
  +0.001 (noise). It does account for most of the cap15-bulk pds leak (−0.00084 raw). g000 minus this K-matched base:
  nmae −0.00682, fid +0.00254, jac +0.00026 (DE non-reach **+0.00342**), pds −0.0003, mse −0.0004.
- **2×2 on the carriers: CD4 and un-whitening are additive.**

  | effect (DE-sum overall; non-reach part) | with the other factor at A | at B |
  |---|---:|---:|
  | CD4 residual (at γ 0.65 / at γ 0) | +0.00571 (+0.00395) | +0.00612 (+0.00374) |
  | un-whitening γ 0.65 → 0 (no CD4 / with CD4) | +0.00260 (+0.00362) | +0.00301 (+0.00341) |

  CD4 moves nmae −0.0093 to −0.0096. Un-whitening moves nmae −0.0065 to −0.0068 and fid +0.0025 to +0.0032.
- **Dose response in γ** (DE-sum): −0.3: −0.0073 · **0: +0.0030** · 0.3: +0.0015 · 0.65: 0 · 0.9: −0.0048. nmae is monotone
  from 0.9 down to 0 and turns up below 0. fid keeps rising below 0, but reach collapses (−0.047 at −0.3). The optimum sits
  at or near γ 0, which is also the point where the DE side uses no whitening operator at all.
- **It is sign information.** Sign accuracy of the delivered LFC on the reference's DE genes (seed mean): base 0.5614,
  **g000 0.5719**, g030 0.5686, g090 0.5568, pre 0.5630, w065nocd4 0.5501, signrep 0.5649, gain 0.5611.
- **Where it lands.** The gain sits in targets with real DE: for n_real ≥ 50 (47 targets), nmae −0.007 to −0.012, fid
  +0.002 to +0.0065, reach +0.004 to +0.014. For n_real < 50 (38 targets), nmae +0.008 and reach −0.018 (the noisy
  stratum). nmae improves on 40 of 58 targets; fid on 54 of 85.
- **Target-sampling uncertainty** (paired bootstrap over targets on seed-averaged values; relevant to transfer, not to
  realisation): DE-sum 95% CI −0.0004 to +0.0073, p(≤ 0) 0.04 under the production bulk; −0.0016 to +0.0074, p 0.10
  under the cap15 bulk. Reach's per-target variance drives it.

## Negatives (recorded with mechanism so nobody re-runs them)
- **228's per-gene gain on the carriers: −0.018.** nmae +0.041, reach −0.033. Some of the gain's multipliers are negative
  (m min −6.7), which flips signs, and the others re-rank the head. It is a bulk-cosine object; it must never touch the DE
  side.
- **Sign prior off: −0.0105** (fid −0.0115). This reconfirms 392's prior at 3 seeds.
- **Prior gated to CLEAN4-unanimous genes (conf4): −0.0058.** Flips drop from 885 to 338 per target, and fid −0.0034.
  The non-unanimous flips also carry correct base-rate signs; confidence gating throws them away.
- **q70 / q90:** inside noise. q70's +0.0012 is reach noise only (non-reach −0.0003). q80 stays.
- **Sign repair (387) on the carriers only: +0.0005**, inside noise. On the board through the bulk its DE half read +0.0014,
  also unresolved.
- **γ 0.9: −0.0048; γ −0.3: −0.0073.**

## Recommendation (one carrier change for a board arm)
**Carriers at γ 0.** Under `VCC_PERT_SCALE_BALLAST`, compose the carrier effect with the whitening step skipped and
everything else identical: stage, PerturbMatch, normrel cap, the CD4 residual (with deflation k 5 and partial rows),
and the CD4 dose re-solved against the unwhitened base, as `apply_frozen_cd4_residual` does. Norm-match each carrier row to
the γ 0.65 effect row it replaces (before the per-target key). The carriers' profile, spread and rank key come from that
matrix. The bulk, the CLEAN4 q80 prior, the classes and the call set are unchanged. In H1 terms this is exactly
`prep416.py`'s `g000`.
- **H1 price (production bulk): +0.0033 overall**, all on the DE channel: nmae +0.0018, fid +0.0011, jac +0.0001, reach
  +0.0003 (noise). This clears the board's ordering-channel sd (0.0045) only at ~0.7 sd, but the move is on nmae and fid,
  which reproduce between realisations to ~1e-4 raw locally. Discount it: H1 has over-predicted DE gains for real contexts
  before (319pA, 399). The factorial structure makes this one more credible than a single-arm reading: two independent,
  additive factors, a monotone dose response, and a sign-accuracy mechanism. A realistic board expectation is
  +0.001…+0.003.
- **October bonus:** the DE side then uses no panel-fitted whitening operator (§8.3 trap), only the CD4 banks, which
  cover 276/300 on A/B/C.
- **Natural follow-up (not run; it would be a third (e) idea):** the CD4 dose on the carriers only. CD4 is worth +0.006
  of DE at the production dose, and nobody has measured its dose response on the DE side.

## Compliance
Serial H1 arms (`drive416.sh`), each held until `memory_pressure` ≥ 40% free and free+inactive ≥ 20 GiB. The gate held
once, for 2 minutes, while another agent's whitener build ran. Peak RSS was 4.0–7.6 GB per arm (the ps runs load the H1 operator for
the key). One deviation: `prep416b.py`, a 5-second numpy job at 4.4 GB, ran while an H1 arm was live, with free+inactive
at 28 GiB. `work/h1/*.npz` and `work/carriers416*.npz` are gitignored (≤ 64 MB each); their shas are in
`work/carriers416*.json`.

## Files
`prep416.py`, `prep416b.py` (carrier matrices), `k_emit416.py` (400's emitter plus an optional prior `gate`; diff = one
block), `run416.py` (harness), `ka416.py` (known answers), `drive416.sh` + `work/queue3.txt` (serial driver, memory
gate), `analyse416.py` (menu table, noise, bootstrap), `analyse416k.py` (K-matched, 2×2, strata, sign accuracy).
Per-run outputs: `work/h1/{metrics,diag,per-target}-<design>_s<seed>[_bps][_K<tag>].*`, logs in `work/logs/`.
