# 346 — The recipient's own labels do buy real per-gene direction information; at the registered full-replacement weight the PDS cost eats it, and 341's retrieval endpoint IS the PDS member at slope 1.01

**Status: COMPLETE — NEGATIVE against the registered kill criterion, with a resolved
positive mechanism and one post-hoc operating point that is not resolved.**

September 10, 2026, 00:12–00:45 UTC. Known-line, held-target, **cross-study** development
evidence on the reused 85-target H1 panel. Not unseen-line generalization, not an A/B/C
forecast, no board conversion factor, no submission, no production change.

## One-paragraph answer for the PI
A calibration map fitted on H1's OWN public non-evaluation perturbations (Pan GSE295214,
1,965 targets) and applied to the 85 held-out H1 evaluation targets raises the design-stage
|lfc|-weighted sign accuracy from **0.6002 to 0.6350 (+0.0348, 95% CI [+0.0241, +0.0465],
48 of 58 targets improved)** — a bigger gain than 341 measured on the source frame, and it
**survives integration into the production cap15 design**. Two controls kill the obvious
alternative explanations: deranged training labels give **−0.0084**, and a target-constant
"mean H1 response" arm gives **−0.0119**, so this is per-target direction information, not a
line-level program offset. Priced on actual emitted counts, that information resolves as a
**gain on nmae, fid, mse and jac** — but the same rotation costs **pds**, and at the
registered β=1 full replacement the PDS loss (−0.0655 raw) swamps everything: H1-local
overall **−0.0149, 95% CI [−0.0332, +0.0015]**. Kill criterion 1 fires. The new mechanistic
result is that **341's retrieval endpoint is the PDS member**: across the five registered
arms, design-stage cosine rank retrieval and delivered `pds` move together at Pearson
**0.9991** with slope **1.009**. 341's retrieval guard was therefore not a proxy — it was
`pds` measured early, and rejecting on it was correct. Because sign accuracy saturates by
β≈0.6 while retrieval keeps falling, a smaller blend weight buys most of the information for
almost none of the PDS: the post-hoc **β=0.25** arm improves **five of six members with
resolved intervals** for an overall **+0.0065, 95% CI [−0.0002, +0.0124]** — positive, right
at the edge of resolution, and not corrected for having been chosen from the curve.

## Known-answer reproductions (three, all passed)
1. **Metadata.** Independently recomputed, the four donor banks named in the registration
   (K562-GWPS, Orion HCT116, Orion HEK293T, Nadig HepG2) share **1,178** targets with the Pan
   H1 study after removing the 300 A/B/C panel targets and the 85 H1 evaluation targets; their
   common canonical gene axis after excluding the training/evaluation target union is
   **5,567**; **30/85** evaluation targets have a direct Pan row. All three match
   `docs/sprint-2026-09-09/NEXT-METHOD-SCREEN.md` exactly. `work/coverage.json`.
2. **343's design-stage sign metric.** My re-implementation of `signs343.py` returns for cap15
   `design_w_all` = 0.600185097698462, `design_w_called` = 0.6492991728970661,
   `design_w_uncalled` = 0.5580069530066262, `design_u_all` = 0.569406545928883 on 58 eligible
   targets with mean gate 1119.5 genes — **absolute difference 0.0 on all four** against
   `experiments/343-nmae-decomposition/work/signs-cap15.json`. `work/signs.json`. The β-sweep
   independently reproduces the same value at β=0 (`work/beta-curve.json`).
3. **323's cap15 six-member emission.** Re-emitted from this directory through 323's
   `emit_and_score` constants: **emitted-count digest, panel-total digest, MSE legs, and every
   array in `stats-cap15.npz` (`targets, sums, jk, p, q, lfc`) are bit-identical to 323's**,
   and **every per-target value in `per-target-cap15.csv` is bit-identical** (all 85×6).
   Five of the six panel means are exactly equal; `nmae` differs by **1.110e-16**.

   *That one ulp, diagnosed rather than waved through:* cell-eval2's `dispatch_de_metrics`
   returns the per-target NMAE rows in a **nondeterministic order** (the two CSVs hold the
   same 58 values in different orders; the other four per-target members — pds, fid, reach,
   jac — came back in identical order both times, and their panel means are exactly equal).
   `np.nanmean` of those identical values reproduces 323's recorded scalar in
   323's row order and this run's scalar in this run's order, exactly, differing by 1.110e-16.
   The registered gate said "exactly"; I therefore **raised the gate after seeing the result**,
   to: exact on the counts, the totals, the MSE legs, every statistic array and every
   per-target value, plus panel means within one ulp with at least five exactly equal.
   `work/known-answer.json` records the change and its reason. The relaxation cannot
   manufacture anything: every arm contrast below is ~1e-3, seven orders of magnitude larger,
   and every member in `work/results.json` is recomputed from the per-target values in one
   canonical target order, so the ordering artifact cannot enter a contrast at all.

## What was measured

**Deviation from the registered donor set, stated up front.** The registration names
K562-GWPS, both Orion banks and Nadig HepG2 as donors. Recomputed, Nadig HepG2 has rows for
only **21** of the 85 evaluation targets, so only **20/85** evaluation targets would have a
four-donor consensus at test time and 65 of 85 design rows would be left untouched — the arm
would be a near no-op. The four-donor construction is **not usable at test time**. Primary
construction therefore uses the **three production donors** (K562-GWPS, Orion HCT116, Orion
HEK293T) — exactly 323's `DONORS`, i.e. the very donor consensus the production design is
built from, which makes the comparison sharper, not weaker.

| quantity | four registered donors | three production donors (used) |
|---|---:|---:|
| training targets (Pan ∩ donors, minus 300 panel and 85 eval) | 1,178 | **1,965** |
| covered gene axis (donor ∩ Pan ∩ canonical, minus target union) | 5,567 | 5,513 |
| …also intersected with the Pan measured mask (row-constant, 16,268 genes) | — | **5,275** |
| evaluation targets with a full donor consensus | 20 / 85 | **80 / 85** |
| DE-universe genes (of 10,779) inside the covered axis | — | **4,836** |
| median share of each cap15 row's L2 norm on the covered axis | — | **0.621** |
| gate \|lfc\| mass on the covered axis (mean over 58 targets) | — | **0.364** |

The 5 evaluation targets with no full donor consensus (**AASS, CNN3, COL6A1, HMGCS1, MED23**)
keep their unmodified cap15 row in every arm. **30 of the 85 evaluation targets have direct
Pan rows; all are EXCLUDED from every fit** (28 of them are among the 80 modified targets).

**Estimator.** Unit-RMS raw directions on the 5,275 covered genes; lines
(k562_gwps, orion_hct116, orion_hek293t, **gse295214_h1**) with H1 outer. 334's
`source_blocks` gives three inner donor recipients, each with the other two donors' mean
direction as X; the ridge is `core341.fit` verbatim (λ = number of rows = 5,895, on
RMS-scaled features; 5,275 features), with X = the three stacked leave-one-out donor
consensuses on the 1,965 training targets and Y = the Pan H1 response tiled three times.
The query is the three-donor mean direction of the 80 covered evaluation targets.

**Integration into the production design (exact).** With `C` = the 5,275 covered columns on
`effects.npz['genes']`, `d = cap15[i, C]`, `p` = the arm's prediction:

    m       = (1-β)·d/‖d‖₂ + β·p/‖p‖₂
    new[C]  = m · (‖d‖₂ / ‖m‖₂)          new[j] = cap15[i, j] for j ∉ C

Verified for every arm: off-`C` coordinates **bit-identical** to cap15, covered-axis row norm
error ≤ **1.9e-14**, β=0 returns cap15 **bit-identically**. Amplitude — a closed axis
(317, 227) — is therefore not what changes; only direction inside the covered subspace.

**Integrity checks, all passed.** Held-target poison replay: refitting with the recipient's
evaluation-target rows replaced by 1e6-scale noise leaves every fitted `gain`, `scale`, `beta`
and every prediction **byte-identical**. Subset-linearity (341's check): the mean of the three
leave-one-out subset predictions equals the prediction on the three-donor mean to
**2.09e-14**. Maximum normal-equation residual **3.91e-11**. Bootstrap MSE and budget
identities hold to < 1e-10 for every arm; **0 non-finite resamples** across all seven
contrasts. Pan's measured mask verified row-constant before use.

## Design-stage result (pre-emission), |lfc|-weighted sign accuracy on the reference gate
58 eligible targets, computed exactly as 343. Paired 10,000-target bootstrap (seed 20260909).

| arm | design weighted sign | vs cap15 | 95% paired CI | better/worse targets |
|---|---:|---:|---|---:|
| cap15 (323 design) | .600185 | — | — | — |
| source_mean β=1 | .600901 | +.000716 | [−.00617, +.00659] | 27 / 21 |
| **diag_cross β=1** (integration control) | .600901 | +.000716 | [−.00617, +.00659] | 27 / 21 |
| same_line_diag β=1 | .600901 | +.000716 | [−.00617, +.00659] | 27 / 21 |
| crossline_full β=1 | .590071 | −.010114 | [−.02208, +.00251] | 20 / 33 |
| **shuffled-label β=1** (control) | .591741 | −.008444 | [−.02675, +.00549] | 24 / 28 |
| **mean_response β=1** (target-constant control) | .588273 | −.011912 | [−.03451, +.00926] | 22 / 33 |
| same_ridge β=0.25 | .613810 | +.013625 | — | — |
| same_ridge β=0.5 | .628311 | +.028126 | [+.01802, +.03992] | 44 / 7 |
| **same_ridge β=1** | **.634995** | **+.034810** | **[+.02412, +.04648]** | **48 / 4** |

The three diagonal arms are numerically identical because a positive per-gene diagonal cannot
change a sign — as in 341. Restricted to the covered axis, where the modification acts:
cap15 .623580 → same_ridge β=1 **.706999 (+.083419)**.
`same_ridge β=1 − mean_response β=1` = **+.046722 [+.02674, +.06845]**, 40/15 targets: the
gain is per-target direction information, **not** a line-level program offset.

## Six-member pricing on actual emitted counts (85 targets, ~145 s and ≤ 4.9 GiB per arm)
Raw members. Eligible targets: pds 85, mse 85, nmae 58, fid 85, reach 82, jac 85.

| arm | pds | mse | nmae | fid | reach | jac | H1-local overall |
|---|---:|---:|---:|---:|---:|---:|---:|
| **cap15** (known answer) | .889636 | .974946 | .950836 | .497555 | .254347 | .074760 | **.280953** |
| diag_cross β=1 | .864006 | .967672 | .953437 | .499429 | .240441 | .074538 | .268688 |
| same_ridge β=0.25 † | .885014 | .958700 | .942417 | .500936 | .269792 | .074920 | **.287472** |
| same_ridge β=0.5 | .869748 | .947174 | .936270 | .504275 | .278587 | .075255 | **.288025** |
| same_ridge β=1 | .824090 | .949992 | .936875 | .506195 | .255345 | .075091 | .266011 |
| shuffled β=1 | .743697 | 1.018906 | .964640 | .497411 | .216286 | .074275 | .214534 |

† post-hoc; see below. Lower is better for `mse` and `nmae`.

Paired 10,000-target bootstrap (seed 20260909, MSE panel ratio and jackknife claim cap
recomputed in every resample), against cap15:

| contrast | overall | 95% CI | pds | mse | nmae | fid | reach | jac |
|---|---:|---|---|---|---|---|---|---|
| **same_ridge β=1** | **−.014942** | **[−.033206, +.001497]** | −.0655 **HARM** | −.0250 **GOOD** | −.0140 **GOOD** | +.0086 **GOOD** | +.0010 | +.0003 |
| same_ridge β=0.5 | +.007072 | [−.004291, +.016096] | −.0199 **HARM** | −.0278 **GOOD** | −.0146 **GOOD** | +.0067 **GOOD** | +.0242 | +.0005 **GOOD** |
| same_ridge β=0.25 † | +.006519 | [−.000150, +.012406] | −.0046 | −.0162 **GOOD** | −.0084 **GOOD** | +.0034 **GOOD** | +.0154 **GOOD** | +.0002 **GOOD** |
| diag_cross β=1 | −.012265 | [−.021825, −.002321] | −.0256 **HARM** | −.0073 | +.0026 | +.0019 | −.0139 | −.0002 |
| shuffled β=1 | −.066419 | [−.089747, −.047176] | −.1459 **HARM** | +.0440 **HARM** | +.0138 **HARM** | −.0001 | −.0381 **HARM** | −.0005 **HARM** |

"GOOD"/"HARM" = the entire 95% interval is beneficial / harmful for that member.
Two further contrasts: **same_ridge β=1 − shuffled β=1 = +.051476 [+.03165, +.07391]** (the
labels resolve strongly against random labels); **same_ridge β=1 − diag_cross β=1 =
−.002677 [−.01999, +.01317]**, made of resolved gains on mse, nmae, fid and jac against a
resolved pds loss of −.0399 — the H1 label information adds nothing NET over the plain donor
consensus at full replacement, while adding real per-member information.

## The mechanism: 341's retrieval endpoint IS the `pds` member
Cosine rank retrieval of each arm's **designed effect** against the reference's real
bulk_lognorm delta, with the panel target genes excluded exactly as PDS's
`exclusion_scope="panel"` does and PDS's midrank / (n−1) rank rule:

| arm | design retrieval | Δ vs cap15 | delivered `pds` | Δ vs cap15 |
|---|---:|---:|---:|---:|
| cap15 | .878011 | — | .889636 | — |
| same_ridge β=0.5 | .864426 | −.013585 | .869748 | −.019888 |
| diag_cross β=1 | .856723 | −.021289 | .864006 | −.025630 |
| same_ridge β=1 | .818768 | −.059244 | .824090 | −.065546 |
| shuffled β=1 | .735994 | −.142017 | .743697 | −.145938 |

Across these five arms, **Pearson 0.9991, slope 1.009, 5/5 signs agree**; per-target
correlation within each arm 0.983–0.989. So the retrieval noninferiority guard that rejected
341 was not measuring a proxy — it was measuring `pds`, the largest scaled member in this
composition (0.856 of cap15's 1.686 total scaled points). **341's rejection was correct, and
the same guard can now be priced instead of argued about.**

## The blend-weight trade-off (design stage; `work/beta-curve.json`)
| β | design weighted sign | Δ | design retrieval | Δ |
|---:|---:|---:|---:|---:|
| 0.000 | .600185 | — | .878011 | — |
| 0.125 | .607021 | +.006836 | .878431 | **+.000420** |
| 0.250 | .613810 | +.013625 | .877171 | −.000840 |
| 0.375 | .618497 | +.018312 | .870868 | −.007143 |
| 0.500 | .628311 | +.028126 | .864426 | −.013585 |
| 0.625 | .633305 | +.033120 | .851681 | −.026331 |
| 0.750 | .633550 | +.033365 | .839916 | −.038095 |
| 0.875 | .634064 | +.033879 | .828571 | −.049440 |
| 1.000 | .634995 | +.034810 | .818768 | −.059244 |

The sign gain saturates by β≈0.625; retrieval falls roughly linearly past β≈0.25. That knee
is why β=0.25 was emitted as one extra, **post-hoc** arm: it captures 39% of the sign gain for
1.4% of the retrieval loss, and it is the only arm here whose overall lower bound comes
within 0.0002 of zero. **It was chosen after the two registered weights were priced; its interval carries no
selection correction and it is not part of any kill criterion.**

## How much of the recipient's own bank the gain needs (`work/sensitivity.json`)
Design-stage sign gain after β=1 integration, refitting on random subsets (3 seeds each):

| training targets | 125 | 250 | 500 | 1,000 | 1,965 |
|---|---:|---:|---:|---:|---:|
| mean sign gain vs cap15 | +.0090 | +.0186 | +.0153 | +.0267 | **+.0348** |

Still climbing at 1,965 and not saturated. The Nadig Jurkat bank has 2,270 perturbation rows, of which
NEXT-METHOD-SCREEN counts **1,986** non-panel targets shared with the three source banks
under the same exclusions — the same order of magnitude as H1's 1,965 here, i.e. enough to
reach roughly this gain and no more. (The Jurkat-side count is quoted from that document,
not recomputed here.)

## Verdict against the PLAN.md kill criteria
| # | criterion | fired? |
|---|---|---|
| 1 | same_ridge β=1 overall paired interval does not exclude zero above +0.002 | **FIRED** — the interval is [−.033206, +.001497] and the point is negative (−.014942) |
| 2 | the shuffled-label control reproduces ≥ half of β=1's gain | **flagged, but vacuous** — there is no positive β=1 gain to share. The control is **−.066419**, resolved harmful on 5 of 6 members; it reproduces **none** of the sign gain (−.0084 at the design stage). Read this as "no gain", never as "the control explained the gain". |
| 3 | design-stage weighted sign accuracy gain < +0.015 | **not fired** — +.034810, CI [+.02412, +.04648] |

**Registered verdict: the regime as registered (full replacement of the design direction by
the same-line calibrated direction) is NOT worth building for context A.** Criterion 1 fired.

**What the PI should actually take from this, stated separately from the registered verdict:**
the criterion tested the wrong operating point, and the underlying question came back
positive. The information regime does supply real, controlled, per-target direction
information for a held-out line (+0.035 design sign, resolved; both controls negative), it
does convert into the members 343 identified (`nmae`, `fid`, `reach`, and even `mse`), and its
entire cost is a `pds` loss whose exchange rate is now measured (1.01 raw pds per raw
retrieval). At the blend weight that respects that exchange rate the net is **+0.0065 with a
95% lower bound of −0.00015** — positive, five of six members resolved beneficial, and one
target-panel's worth of noise away from resolution. That is not a clearance, and one post-hoc
β on an 85-target panel is not evidence to build a context-A pipeline on. It is a live
candidate that would need a pre-registered β, a second recipient line, and an independent
panel before anyone spends A slots on it.

## Limits
- **H1 is one line, and one panel.** 85 targets, 58 NMAE-eligible, 82 reach-eligible; one
  fixed emission/control seed; one reference gallery; one set of H1-local anchors.
  Development evidence, not an A/B/C forecast, and no board conversion factor is applied.
- **The H1-local overall clips MSE at 0**, and cap15's scaled MSE sits at the clip
  (raw .974946 vs baseline anchor .970691). Every arm here improved raw MSE, so this
  clipping did not hide harm — but on the board, where we lead MSE, the same rotation would
  be scored on a different scale. Read the raw MSE column, not the local overall, for MSE.
- **30 of the 85 evaluation targets have direct Pan H1 rows and are excluded from every fit**;
  the poison replay proves the exclusion byte-for-byte. 5 targets (AASS, CNN3, COL6A1,
  HMGCS1, MED23) are unmodified in every arm, so all contrasts are diluted by 5/85.
- **Same line, different study.** Pan GSE295214 (hESC H1) → the VCC-2025 H1 challenge data.
  Context A would additionally assume A *is* Jurkat (041/315 evidence, not certainty), so A
  carries one transport assumption more than this experiment does.
- The covered axis reaches 4,836 of the 10,779 DE-universe genes and 36.4% of the gate's
  |lfc| mass; genes outside it keep the design value, so this bounds the achievable effect.
- Ridge λ, feature set and architecture are 341's, not retuned. β=0.25 is post-hoc.
- Reused targets, correlated genes and study/line confounding remain, as in 341.

## Exact commands to reproduce
```
cd ~/Code/virtual-cell-challenge-2026
E=experiments/346-same-line-calibration-h1
export NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4 POLARS_MAX_THREADS=4 OPENBLAS_NUM_THREADS=4
uv run python $E/coverage346.py            # metadata known answer: 1178 / 5567 / 30
uv run python $E/fit346.py                 # fits + held-target poison replay
uv run python $E/effects346.py             # norm-preserving integration, invariants
uv run python $E/signs346.py               # 343 sign known answer (difference 0.0)
uv run python $E/extra346.py               # target-constant control + paired sign intervals
uv run python $E/run346.py all             # cap15 known answer, then the four registered arms
uv run python $E/run346.py recheck         # re-evaluates the cap15 gate from the receipt
uv run python $E/retrieval346.py           # design retrieval vs delivered pds
uv run python $E/beta346.py                # blend-weight trade-off curve
uv run python $E/sensitivity346.py         # training-set-size requirement
uv run python $E/knee346.py                # POST-HOC beta=0.25 effects
uv run python $E/run346.py same_ridge_b025 # POST-HOC emission
uv run python $E/report346.py              # paired 10,000-target bootstrap -> work/results.json
```
`run346.py all` was executed once before the known-answer gate was widened; the four
registered arms were then run as
`run346.py diag_cross_b1 same_ridge_b1 same_ridge_b05 shuffled_b1`.
`work/results-registered-arms.json` preserves the report as it stood before the post-hoc arm.
`work/fourdonor-*.txt` hold the four-donor metadata known answer (1,178 / 5,567);
`work/fitted-*.txt` hold the axes actually fitted (1,965 / 5,275 / 80).
Nothing outside `experiments/346-same-line-calibration-h1/` was written; no git command, no
builder, no submission, no reserve read, no download.

## Peak memory and wall time
| step | wall | peak RSS |
|---|---:|---:|
| `fit346.py` (3 ridge fits + full poison replay) | 31.6 s | 9.14 GiB |
| `sensitivity346.py` (13 ridge fits) | 25.6 s + sign pass | 9.22 GiB |
| `run346.py` cap15 (known answer) | 171.9 s emit+score | 4.82 GiB |
| `run346.py` 4 registered arms | 633.8 s | 4.85 GiB |
| `run346.py` same_ridge_b025 | 146.0 s | 4.58 GiB |
| `report346.py` | 1.3 s | 0.38 GiB |
| `coverage/effects/signs/extra/retrieval/beta/knee` | ~6 min total | < 1.5 GiB each |

Per-arm emission+scoring: 141–172 s. `memory_pressure` was checked before every arm inside
`run346.py` (gate: ≥ 50% free); observed 70–75% free throughout, with the 22 GiB production
build (pid 9220/10210) running the whole time. Wall clock 00:13:55 (first coverage output) to 00:45:05 (final results.json), 31 minutes.
The three fitted ridge operators (640 MB, `work/model-*.npz`) were deleted after use; they
regenerate deterministically in 32 s via `fit346.py`. `work/` is 350 MB and gitignored.
