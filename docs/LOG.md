# LOG

## 2026-09-01 — STANDING RULE: never spend a slot on rank

Ahmet, after the a15 restore: the goal is **first place**, not a better end-of-day row.

On 2026-08-31 a19 (whitener ridge, gamma 0.80 / lam 0.007) scored `0.160569` against a15's
`0.174087` and our displayed rank fell 12 -> 26. The day's second slot was then spent
re-uploading `a15-cd4f010.vcc` **byte-identical** to put the better number back on the board.
It worked — rank 12 restored — and it bought **zero information**, because we already knew
that package's score exactly.

**Rule, applying to every agent and every session:**

- Never spend a submission slot to restore or improve leaderboard rank.
- Never re-submit a package that has already been scored.
- If a bad arm lands, leave it displayed. The row is informational; only the last entry after
  the Oct 22 test release is prize-eligible.
- Every slot must buy information that moves us toward first. If an arm cannot plausibly do
  that, it is not worth a slot — and an already-scored package provably cannot.

Recorded in `AGENTS.md` §6 (banner) and §11.1 rule 7, and `docs/CURRENT_STATE.md`.

# Lab notebook

Running record of what was tried, what happened, and what it means. Newest at the bottom.
Every entry: **what I expected → what I observed → what changed as a result.**

Scores in this file are **local** (2025 H1 harness) unless marked LEADERBOARD.

---

## 2026-08-21 · Setup

Workspace, `uv` env, official controls bundle (630 MB), vendored `vcc2026` metric spec.
`vcc` authenticated as `<redacted-email>`, team "Ahmet Dedeler".

Verified against the real data rather than the docs:
- `gene_names.csv` order starts `TSPAN6, TNMD, DPM1` — the website mockup shows `SAMD11, NOC2L,
  …`, which is wrong. Both reference CSVs carry a header row.
- Context controls: 18,400 × 18,533, CSR float32, all-integer, median lib 20,109, ~5,973 nnz/cell.
- Log-pseudobulk Pearson between contexts: A↔B 0.742, A↔C 0.685, B↔C 0.793. Genuinely different
  tissues; this is a transfer problem, not a denoising problem.

`scripts/build_submission.py` written and verified against `vcc prep` — streams CSR into HDF5 so
a full 360k-cell submission never has to be materialised (~16 GB at 5,600 nnz/cell).

---

## 2026-08-21 · Leaderboard column mapping (correction)

**Expected:** one number per metric.
**Observed:** the leaderboard prints **(raw, scaled) pairs**, header order `pds mse jac nmae fid
reach`. Confirmed by arithmetic — scaled values average to Overall for all 10 rows, max residual
0.0002.
**Consequence:** my first reading of the field was wrong in a way that mattered. Corrected picture:

| | pds | mse | jac | nmae | fid | reach |
|---|---|---|---|---|---|---|
| top-10 scaled | 0.348–0.627 | 0.000 (all ten) | −0.025–0.008 | −0.034–0.120 | −0.171–−0.021 | 0.029–0.104 |

The board is a **`pds`-only contest**. Five of six metrics are unclaimed.

---

## 2026-08-21 · Experiment 001 — magnitude calibration is free points

**Hypothesis:** `pds_cosine` uses cosine distance so it is scale-invariant; `expr_mse` is
quadratic in effect magnitude. If so, magnitude α is a free parameter.

**Observed:**
- `PDS_p` bit-identical for α ∈ {0.01, 0.1, 0.25, 1, 10, 100}. Invariance is exact.
- `min_α mse(α) = 1 − cos²θ` at `α* = ⟨d,δ⟩/⟨d,d⟩`. Numerics match analysis to 4 dp.
- The `mse(α=1)` curve reproduces the leaderboard: cos 0.5 → 3.00 (Vivai 3.48, Sidechain 3.61,
  joshua 3.96); cos 0.3 → 10.1 (Gia Khang 6.54, p19 6.97); cos 0.7 → 1.04 (wombat 1.56).

**Airtight form of the claim:** `min_α mse(α) ≤ mse(0) ≈ 1.003` (control-copy), so any raw mse > 1
is provably not α-optimised. All ten top teams are ≥ 1.562.

**Consequence:** α (and its counter-knob, NB emission dispersion) must be swept before any
modelling work. → experiment 002.

---

## 2026-08-21 · Data acquisition

- 2025 VCC H1 (34 GB): direct public HTTPS from the Arc GCS bucket. Test + Validation landed,
  Training in flight.
- Replogle 2022 via Figshare: **blocked** — `x-amzn-waf-action: challenge`, 403 with a browser UA.
  Routed around it via **scPerturb** (Peidli et al., Nat Methods 2024, Zenodo 10044268), which
  mirrors all three Replogle sets *harmonised* and much smaller (K562_gwps 8.8 GB vs 61.3 GB).
  Strictly better: standardised `.obs` schema across studies makes them one corpus.
  Pulling all 50 files / 39.8 GB.
- Nadig 2025 HepG2 + Jurkat from GEO: works directly.

---

## 2026-08-21 · Local scoring harness

`scripts/make_harness.py` reproduces 2026 reference conditions on 2025 H1 ground truth:
subsample to exactly 400 cells/perturbation, then **binomially thin** to a median of 20,000 UMI
(2025 ships at ~54,000). Uniform thinning rather than per-cell renormalisation, to preserve the
library-size spread.

`h1-test`: 52,400 cells × 18,080 genes, **85 perturbations** (15 of 100 dropped for <400 cells)
+ 18,400 control. Post-thinning median 19,997 UMI. Depth rate 0.3684.

Next: `cell-eval2 baseline` → `prep-real-bundle` gives the (b=0, r=1) competition scale locally.

---

## 2026-08-21 · Harness anchors vs the competition's

`cell-eval2 baseline` on `h1-fast` (40 perts × 400 cells + 4,000 controls) took **6:50**,
single-threaded. Comparing our measured baseline `b` against the official 2026 values:

| metric | h1-fast `b` | official 2026 `b` (A/B/C) | agree? |
|---|---|---|---|
| `pds_cosine` | 0.5019 | 0.500 exactly | ✅ |
| `expr_mse_unbiased_capped_norm` | 0.9555 | 0.986–0.992 | ≈ |
| `de_wilcoxon_sig_jaccard` | 0.0281 | 0.021–0.037 | ✅ |
| `de_wilcoxon_lfc_nmae` | 0.9375 | 1.0009–1.0017 | ≈ |
| `de_wilcoxon_direction_reach_raw` | 0.1682 | 0.047–0.097 | ✗ high |
| `de_wilcoxon_direction_fidelity_yield_raw` | 0.1326 | 0.505–0.522 | ✗ **low** |

`pds_cosine` landing on 0.500 is a strong signal the harness is wired correctly — the spec says
the mean-response baseline carries no discriminating information once panel target genes are
removed, so every row ties and the value is exactly 0.5. We reproduce it.

**`fid` is the outlier and it is explainable.** `F_p = k / max(n_pred, n_real)`. Our harness has
4,000 control cells against the competition's 18,400, so the Wilcoxon has much less power on
*both* sides — but the baseline's own calls (`n_pred`) collapse faster than the reference's
(`n_real`), and the ratio is what `fid` reads. Expect `h1-test` (18,400 controls) to land closer.

**Interpretation rule for this repo:** the harness answers *"does change X help?"*, not *"what
will I score on the leaderboard?"* `b` and `r` are both measured on the harness, so the scale is
internally consistent, but absolute values are not transferable to the 2026 contexts. Never quote
a harness number as a predicted leaderboard score.

Gotcha logged: `prep-real-bundle --baseline` wants the baseline **prediction .h5ad** (from
`baseline --save-pred`), not the baseline output directory. The error is an opaque
`IsADirectoryError` from h5py.

## 2026-08-21 · Iteration-cost decision

Timing on this machine (M5 Max, 18 cores, all cell-eval2 work single-threaded via scanpy's
`rank_genes_groups`):

| harness | cells | perts | baseline run |
|---|---|---|---|
| `h1-fast` | 20,000 | 40 | 6:50 |
| `h1-test` | 52,400 | 85 | ~35–40 min |

`num_threads: -1` is already set in the vcc2026 preset but the process stays pinned to one core,
so the parallelism has to be **across configs, not inside a run**. `sweep.py` now builds all
prediction files first, then scores them through a `ProcessPoolExecutor`. At 6 workers a 36-config
grid is ~45 min instead of ~4 hours.

Two harnesses kept deliberately: `h1-fast` for iteration, `h1-test` for any go/no-go call. Panel
size is not neutral — `pds` is a rank *within* the panel, so a 40-gene panel is a strictly easier
ranking problem than an 85- or 300-gene one.

---

## 2026-08-21 · Harness validated against the official scale

`prep-real-bundle` on `h1-fast` gives both ends. Comparing the **replicate anchor `r`** — the
half-split self-agreement, the hardest number to fake — against the official 2026 bundles:

| metric | our `b` | our `r` | official `r` | in range? |
|---|---|---|---|---|
| `pds_cosine` | 0.5019 | 0.9501 | 0.927–0.984 | ✅ |
| `expr_mse_unbiased_capped_norm` | 0.9555 | 0.0076 | 0.028–0.045 | ≈ |
| `de_wilcoxon_direction_fidelity_yield_raw` | 0.1326 | 0.8636 | 0.795–0.832 | ✅ |
| `de_wilcoxon_direction_reach_raw` | 0.1682 | 0.9729 | 0.958–0.978 | ✅ |
| `de_wilcoxon_sig_jaccard` | 0.0281 | 0.4698 | 0.375–0.423 | ✅ |
| `de_wilcoxon_lfc_nmae` | 0.9375 | 0.3653 | 0.369–0.431 | ✅ |

**All six replicate anchors land in or adjacent to the official ranges.** The harness measures the
same quantities the competition does.

Seed noise (sd / span), our value vs the spec's: pds 1.2% (0.5–1.3), mse 0.4% (0.8–2.2),
fid 1.7% (1.5–4.8), reach 1.0% (0.2–0.9), **jac 6.7% (0.8–2.2)**, **nmae 4.6% (0.7–2.7)**.
The 40-perturbation panel makes `jac` and `nmae` noticeably noisier than the competition's 300.
**Resolution floor on `h1-fast`: ignore `jac` differences under ~0.07 of span.**

**Correction to my earlier hypothesis.** I predicted `fid`'s low baseline was a control-count
artefact and would rise on `h1-test` (18,400 controls). It went the *other* way: 0.133 → 0.0355.
More control cells give the reference more DE-detection power, so `n_real` grows faster than the
mean-response baseline's own calls `n_pred`, and `fid = k/max(n_pred, n_real)` falls. The 2026
lines must have more *correlated* perturbation responses, making the mean response a genuinely
useful directional predictor there. Net effect: our `fid` span is 0.731 vs the competition's
0.28–0.33, so **the harness understates `fid` gains by ~2.4x**. Conservative, which is the safe
direction, but do not read `fid` deltas as leaderboard deltas.

---

## 2026-08-21 · The VCC panel excludes common-essential genes

Checked VCC-panel overlap against four essential-gene screens:

| dataset | perturbations | VCC-panel hits | expected by chance |
|---|---|---|---|
| Replogle K562 **essential** | 2,058 | **0** | ~31 |
| Replogle **RPE1** (essential) | 2,394 | **0** | ~36 |
| Nadig HepG2 (DepMap essentials) | 2,394 | **0** | ~36 |
| Nadig Jurkat (DepMap essentials) | 2,394 | **0** | ~36 |
| Replogle **K562 genome-wide** | 9,867 | **272 (91%)** | — |

Zero across four independent essential panels is not chance. The 2026 panel was **deliberately
built to avoid common-essential genes** — consistent with the stated selection criterion (>80%
on-target knockdown with a measurable, non-lethal response). Knock down an essential gene and the
cells die before you can profile them.

**Consequences, and they are large:**

1. **Three downloads are dead as direct donors.** Nadig HepG2+Jurkat (13 GB) and Replogle
   essential+RPE1 contribute *zero* perturbation-matched signal. They keep secondary value for
   basal-state modelling and general response structure, but "look up what gene X does" is out.
2. **K562-gwps is the only high-coverage donor: 272/300.** That is why "K562-gwps donor transfer"
   appears in the #5 team's model string. Anyone who worked this out is using it.
3. **Gene-axis gap:** gwps carries 8,248 genes, of which only **7,091 of the 9,923 genes >5 CPM
   in context A (71%)** are present. Even a perfect donor transfer leaves 29% of the DE-scored
   genes with no signal.

**Where differentiation lives:** the 28 uncovered perturbations and the 29% uncovered genes.
Everyone with gwps gets the same 91% × 71%. The edge is in what you do with the rest — gene-gene
correlation structure, pathway priors, or a model that generalises beyond looked-up responses.

Also worth noting: RPE1 and K562-essential still give **two extra cell lines** for the ~2,050
essential genes they share with gwps. That is the only place we can directly observe *the same
perturbation in three different lines* — which is exactly the cross-context transfer function we
need to learn. Wrong genes, right structure. Keep them.

---

## 2026-08-21 · Experiments 002 and 003 — the two numbers that set the plan

**002 (α sweep, real scorer).** `pds_cosine` returned `1.11127002` at α ∈ {0.1, 0.25, 0.5, 1.0} —
identical to eight significant figures. Scale-invariance confirmed on the real scorer.
`argmax(mse)` landed at α = 0.25, exactly the `α* = ⟨d,δ⟩/⟨d,d⟩` the construction predicts.

**But my headline claim was too big.** I said α-tuning was worth +0.04 overall. It is worth
**+0.010** on this harness (avg 0.4895 at α=0.5 vs 0.4799 at α=1.0), because `fid`, `reach` and
`jac` all *rise* with α — bigger effects clear the significance bar — so shrinking α for `mse`
costs the DE metrics. `mse` is one sixth of the score and the rest pull back. Real lever, about a
quarter of the size I first claimed. Full table in `experiments/002-calibration-sweep/RESULTS.md`.

Also: `pds` raw saturated at **1.0000** at cos θ = 0.5 on the 40-perturbation panel. Discrimination
is trivial at that panel size, so round 1 measured an unrealistically easy regime.

**003 (cross-line transfer).** Replogle K562 vs RPE1, 2,187 perturbations in both lines:

- same gene across lines: cos = **0.142**
- vs the line's mean response: 0.097
- different gene (null): 0.044

Decomposing `δ_g = m + s_g`: **`m` is the scored baseline and earns exactly zero.** After removing
each line's mean response, the gene-specific part transfers at **cos 0.128 vs a null of 0.003** —
40× the null, positive for 86% of perturbations. Real, and small.

**So the naive-transfer operating point is cos ≈ 0.13, not 0.5.** Every calibration conclusion
from 002 has to be re-derived there.

**Negative result:** donor effect magnitude does not predict transferability
(Spearman 0.063, non-monotone across quintiles). That kills the obvious per-perturbation α
weighting. Best transfer is concentrated in core RNA-processing and mitochondrial machinery
(SMG5, CMTR1, GEMIN5, UPF2, TFAM, POLRMT, METTL3), which points at **pathway identity** as the
confidence signal to try instead.

**Engineering note:** scPerturb stores K562-gwps as a *dense* gzipped array — 1,989,578 × 8,248
float32 is 65.6 GB uncompressed and cannot be loaded on this machine. `build_donor_table.py`
streams 16k-row blocks and accumulates group sums with a sparse indicator matmul, holding peak
memory at ~0.5 GB regardless of dataset size.

---

## 2026-08-21 · Projecting harness results onto the official scale (two corrections)

**Correction 1 — I had the harness/competition `fid` relationship backwards.** I wrote that the
harness *understates* `fid` gains. It **overstates** them. Harness `b_fid` = 0.133, span 0.731;
official `b_fid` = 0.514, span 0.300. A raw fid of 0.58 scores 0.611 on the harness and 0.222 in
the competition — the harness is ~2.75× generous, because its baseline sits much lower so there
is more room above it. The fix is not a correction factor: **read RAW metric values out of
`agg_results.csv` and re-scale them with the official anchors.** That is now the standard
translation and it is what the table below uses.

**Correction 2 — "mse is free money" only holds at directional accuracy we do not have.**
`min_α mse = 1 − cos²θ`, and scaled mse exceeds 0 only when `1 − cos² < b`:

- official `b` = 0.989 → needs **cos > 0.105**
- harness `b` = 0.9555 → needs cos > 0.211

At the measured gene-specific transfer of **cos = 0.128**, `min mse = 0.9836`, which scales to
**0.006** — worth +0.001 overall. Essentially nothing. The "up to +0.167" figure required cos ≈ 0.7.
`mse` becomes worth having only above cos ≈ 0.3.

**The calibration lever is still large, but it lives in `fid`, `reach` and `nmae`, not `mse`.**

### h1-fast raw values re-scaled with official anchors (pds excluded — saturates at 40 perts)

| cos | α | disp | mse | fid | reach | jac | nmae | **avg5** |
|---|---|---|---|---|---|---|---|---|
| 0.15 | 0.25 | 0.10 | 0.000 | 0.225 | 0.314 | 0.125 | −0.000 | **0.133** |
| 0.15 | 0.10 | 0.50 | 0.000 | 0.054 | 0.126 | 0.116 | 0.114 | 0.082 |
| 0.15 | 0.25 | 0.50 | 0.000 | 0.115 | 0.180 | 0.125 | −0.036 | 0.077 |
| 0.15 | 0.50 | 0.10 | 0.000 | 0.168 | 0.329 | 0.130 | **−1.109** | −0.096 |
| 0.15 | 0.50 | 0.50 | 0.000 | 0.135 | 0.256 | 0.127 | **−1.185** | −0.133 |

**Leader Vivai on those same five metrics: 0.0354.**

Calibration alone moves avg5 from −0.133 to +0.133 **at identical directional accuracy** — a
0.27 swing, roughly twice the leader's entire overall score. `nmae` is the landmine: it craters
to −1.1 at α = 0.5 and is the single largest downside risk of getting α wrong.

## 2026-08-21 · Mean vs gene-specific response, across lines

| quantity | cos |
|---|---|
| mean K562 response → mean RPE1 response | **0.491** |
| gene-specific (centered) transfer | 0.128 |

The shared mean response transfers ~4× better than the gene-specific part — but it is only 10% of
a typical K562 effect's magnitude (‖mean‖ 0.72 vs mean ‖δ‖ 7.27). In RPE1 it is 57% (8.82 vs
15.42), because RPE1's panel is all-essential genes and therefore dominated by a common growth
arrest program, while gwps is genome-wide and most perturbations do little.

Consequence for the model: `m` is the easy, transferable part and it is what the scored baseline
already captures. `s` is 90% of the magnitude in a genome-wide panel, transfers at 0.128, and is
where every point comes from.

---

## 2026-08-21 · Experiment 004 — real donor transfer, and two dead ends

**Operating point measured.** K562-gwps → H1: centered cos **+0.0727** on the full gene axis,
**+0.0949** on the donor-covered genes. Positive for **100%** of the 82 covered perturbations.
Compare K562 → RPE1 at +0.128 (exp 003). H1 hESC is a harder target than a somatic line, so
0.073 is likely pessimistic for the 2026 contexts and 0.128 more representative.

**Dead end 1 — low-rank denoising.** Truncating the donor SVD loses monotonically: rank 5 gives
0.039, rank 800 gives 0.052, full rank 0.073. The spectrum is flat — 800 components explain only
51% of variance. The transferable signal is high-dimensional. **Any bottleneck architecture
destroys it**, which is a mechanism for why simple baselines beat deep models here.

**Dead end 2 — per-gene transfer weights.** Fitted on H1 train+val → 0.071; fitted on K562→RPE1
and applied to H1 → 0.067; shuffled null 0.058; unweighted **0.095**. Real information (beats its
null) but worse than uniform: cosine is about direction, and noisy per-gene rescaling distorts
direction more than it fixes magnitude.

**Conclusion: uniform full-rank donor transfer is the estimator to beat.** Remaining headroom is
calibration, not transforms of `s`. Stop tuning the donor matrix without a mechanism explaining
these two failures.

---

## 2026-08-21 · Experiment 005 — calibrating the real model, and an error in my own reasoning

First 5 configs of the real K562-gwps donor transfer on `h1-fast`, rescaled to official anchors:

| α | β | mse | fid | reach | jac | nmae | avg5 |
|---|---|---|---|---|---|---|---|
| 0.05 | 0 | 0 | **−0.862** | 0.028 | −0.044 | 0.012 | −0.173 |
| 0.10 | 0 | 0 | **−0.785** | 0.073 | −0.041 | 0.024 | −0.146 |
| 0.20 | 0 | 0 | **−0.768** | 0.076 | −0.031 | 0.007 | −0.143 |
| 0.05 | 1 | 0 | −0.763 | 0.016 | −0.045 | 0.007 | −0.157 |
| 0.10 | 1 | 0 | −0.741 | 0.034 | −0.046 | 0.008 | −0.149 |

Leader Vivai on the same five: **+0.035**. We are far worse, and `fid` is the whole story.

**The error.** I computed `α* = ⟨d,δ⟩/⟨d,d⟩ = 0.083` analytically and treated it as *the*
optimum. It is the optimum **for `expr_mse` only** — the α minimising squared error. But at
cos = 0.07, `min mse = 1 − cos² = 0.995 > b = 0.989`, so **`mse` is clamped to 0 for every α**.
I set the one free parameter to serve the single metric that cannot respond to it.

**What α actually has to do.** The DE metrics need our predicted effects to be *the right size*,
so that roughly as many genes clear the Wilcoxon/BH bar as the reference finds:

| target | α |
|---|---|
| minimise squared error (serves `mse`, which is unreachable) | 0.083 |
| **match true effect magnitude** (‖true δ‖ 5.68 / ‖s‖ 7.98) | **0.712** |

At α = 0.05–0.2 our effects are 4–15× too small, almost nothing reaches significance, `n_pred`
collapses, and `fid = k/max(n_pred, n_real)` takes its coverage penalty — hence −0.86.

**Rule going forward: when a metric is clamped, it must not drive a shared parameter.**
`mse` contributes 0 at any α we can reach, so α is set entirely by the DE metrics.

β = 1 (adding the donor's mean response) is roughly neutral-to-harmful: +0.016 at α = 0.05,
−0.003 at α = 0.10. Consistent with exp 003 — `m` is the baseline direction and earns nothing.
Keeping **β = 0**.

→ cal2 running: α ∈ {0.7, 1.0, 1.5, 2.5}, β = 0, dispersion 0.1.

---

## 2026-08-21 · Operational bugs worth not repeating

**1. Model names are published.** v1 went up as
`donor-transfer k562-gwps alpha=0.7 gene-norm cpm20, even-spread emission` — a complete recipe.
The donor and emitter were already public via other teams' model strings, but **`alpha=0.7` was
not**: no team on a 65-team board has raw `mse` below 1, so nobody else is shrinking. We
published our one unclaimed lever. Later names: random jokes, not `v3`. Still zero recipe. See AGENTS.md §14.

**2. `pgrep -f` self-matches the waiting shell.** This waiter never exits:

```bash
until ! pgrep -f "make_submission.py" >/dev/null; do sleep 60; done
```

The shell running the loop has `make_submission.py` in its own command line, so `pgrep` always
finds a match. Four waiters sat spinning after the v3 build had already finished at 19:20 —
caught by GPT-5.6 from the outside, not by me. Use the bracket trick, which cannot match the
literal text of the command that contains it:

```bash
until ! pgrep -f "[m]ake_submission.py" >/dev/null; do sleep 60; done
```

Better still, wait on the artefact rather than the process — a completed file with the right
size is the actual condition, and it does not lie:

```bash
until [ -s submissions/v3.h5ad ] && ! pgrep -f "[m]ake_submission" >/dev/null; do sleep 60; done
```

---

## 2026-08-21 · `vcc prep` peak memory, measured

`/usr/bin/time -l` on a full 360,000-cell submission (2.67e9 stored entries, 21.4 GB h5ad):

```
peak memory footprint   91,180,499,216  (91.2 GB)
instructions retired   593,764,822,508
```

**91 GB peak on a 64 GB machine.** It completed — macOS compressed-memory and swap absorbed it —
but the ratio matches the 3.65× figure a participant reported in the Challenge Discord's
`#submission` forum, and it means the density cap admits submissions that cannot be packaged on
an ordinary workstation. Our 56.3%-of-cap file is already there.

Practical consequences for us:
- Do not run `vcc prep` concurrently with a scoring sweep or a donor build. It will swap, and on a
  worse day it will be OOM-killed with a bare SIGKILL and no message.
- If we ever push density toward the cap, packaging becomes the binding constraint before scoring
  does. Budget headroom accordingly.

---

## 2026-08-21 · Experiment 017/019/020 — spectral whitening, the largest gain so far

**Expected:** nothing in particular. This started as an inversion of an old negative result.

**The reasoning.** Experiment 004 truncated the donor SVD and lost *monotonically* — rank 5
gave transfer cosine 0.039, rank 800 gave 0.052, full rank 0.073. That was read as "the signal
is high-dimensional, stop compressing it", and it closed the whole family. The stronger reading
is that truncation **deletes** the low-variance directions, so the natural next experiment is
the operation that **amplifies** them:

```
P -> P (C + lam I)^(-gamma/2),     C = B^T B / n,  B = donor effect vectors
```

The mechanism is hubness. `pds` is a retrieval metric — is this prediction nearer its own true
effect than to the other 299? — and perturbation effects are dominated by a few shared
programmes (growth arrest, ribosome, stress). Two unrelated knockdowns that both slow the cell
down look alike under a cosine, so whichever true effects sit near the dominant programme are
near everything. Retrieval and word-embedding post-processing both fix this by flattening the
space. Nobody had asked whether it applies here.

**Observed:**

| | macro PDS |
|---|---|
| incumbent (gene_norm, 4:1, no whitening) | 0.7378 |
| + whitening, gamma = 0.65 | 0.7873 |
| **+ whitening and equal donor weights** | **0.7978** |

**+0.0600**, on the harness that retro-predicted the v3 regression, consistent across all three
folds, and replicated on a second 4-line LODO with three predictors and a mostly different gene
population (+0.0433 there).

**Three things that make this credible rather than a lucky knob:**

1. **The basis behaves differently in the two ways it should.** Fitted on the 300 panel
   predictions (259 × 7,760, rank-deficient), the sweep rises monotonically to the edge of
   every grid — the signature of an artefact, since at the limit it just orthogonalises 259
   vectors and improves a ranking by construction. Fitted on all 12,000 donor perturbations,
   the covariance is full rank and there is a clean interior optimum at gamma ≈ 0.5–0.65.
2. **It survived a harness correction that was made specifically to try to kill it.** `lodo.py`
   was normalising the held-out line's effects with the same `gene_norm` it applied to the
   predictors, so the *ground truth* was a normalised object, which the scorer's never is.
   Whitening is exactly the sort of method a normalised truth would flatter. Truth switched to
   the raw held-out delta; effect unchanged.
3. **It is insensitive to the ridge** over four orders of magnitude, and the gamma optimum is a
   broad plateau rather than a spike.

**What changed as a result — and this is the bigger finding.** Experiment 010 established that
equal donor weighting is *worse than a single donor* because averaging smooths every prediction
toward the shared response. With whitening on, equal weighting is **the best** configuration and
4:1 costs 0.009. The optimum moved from "dominant primary" to "equal".

That also resolves the 010-vs-016 contradiction logged earlier. The operative variable is
neither predictor count nor gene population on its own: averaging cancels independent noise
(good) *and* concentrates predictions on the shared response (bad for a ranking metric), and
whichever dominates decides the winner. Whitening removes the second, so averaging wins cleanly.

**Consequence for how to read the rest of this notebook.** Several rejected ideas were rejected
with *the same mechanism* attached — EB shrinkage ("+10.6% cosine, −0.017 on the board, wrong
objective"), per-gene inverse-variance pooling ("over-weights whichever donor expresses a gene
highly"). Both descriptions are the shared-programme problem. → experiment 022 re-tests the
graveyard under whitening.

**Negative results from the same runs.** Rank-1 removal ("all-but-the-top") captures almost
none of what whitening does (+0.0008), so the hubness is spread across many directions rather
than concentrated in one. Per-gene modulation by the **target's** basal profile is worth +0.003
at `p^0.25` and actively harmful by `p^1.0`. Together with experiment 010's basal-similarity
failure, the control profile has now failed to be useful for anything except defining the gene
universe — which is worth stating plainly, because it is the only context-specific signal the
Challenge gives us.

**A separate operational finding.** Whitening amplifies low-variance directions, so `|predicted
effect|` becomes a *worse* per-gene confidence signal after whitening (proxy `fid` 0.5364 →
0.5242, proxy `reach` 0.0309 → 0.0190 when ordering by `|P|`). Ordering by the un-whitened
cross-donor geometric signal from experiment 013 holds up (`fid` 0.5379, `reach` 0.0288). So
the mean profile and the per-gene confidence should come from **different objects**: whitened
effects for `pds`/`mse`/`nmae`, un-whitened agreement for the DE ordering. That is exactly the
two-channel separation experiment 013 proposed, now with a reason to exist.

**v5:** `gene_norm`, no EB, three donors at **1:1:1**, whitening on the donor basis at
**gamma = 0.65**, `alpha = 1.0796` chosen so the emitted median effect norm is 4.801 —
identical to v4b — so **direction is the only change**.

---

## 2026-08-21 · Experiments 023–026 — sizing the whitening change, and three more negatives

**023 — count-level, on the metrics LODO cannot see.** Whitening reshapes the relative
magnitudes of genes inside each prediction, which is what the four DE metrics read. Real
`cell-eval2` on `h1-fast`, at α = 1.08 / CPM gate 5:

| | fid | reach | nmae | jac |
|---|---|---|---|---|
| γ = 0 | 0.5178 | 0.1538 | 0.9397 | 0.0728 |
| γ = 0.65 | 0.5121 | **0.1532** | 0.9531 | 0.0747 |

About **−0.006 overall** on the official scale, against a `pds` gain worth +0.017 to +0.023.
Net clearly positive.

**`reach` is flat, and experiment 018's proxy said it would collapse 39%.** The proxy ordered
the reference's DE genes by `|predicted effect|`; the real metric orders by adjusted p-value
from a Wilcoxon on emitted counts, which under a minimum-variance emitter is `|effect|`
measured **in units of the control's own cell-to-cell spread**. That per-gene rescaling
absorbs most of the damage. Proxy `reach` is not trustworthy even directionally — deleted from
the decision set.

**The CPM gate was costing us, six configurations out of six.** Gate 5 beats gate 20 by
+0.005 to +0.013 avg, carried by `reach` (0.1416 → 0.1532) and `nmae` (0.9635 → 0.9531),
worth roughly **+0.006 overall**. The scorer tests every gene above **5** CPM; gating at 20
made us predict *no change* on ~1,900 genes it does test, which then enter `reach`'s pool with
an arbitrary sign and `jac`'s union for free. v1 and v4b both shipped gate 20; the flag's own
docstring already said 5.

**α = 1.08 is the measured optimum** (avg 0.2306 / **0.2362** / 0.2095 at α = 0.7 / 1.08 /
1.6), decided by `nmae` (0.9623 / 0.9531 / **1.0580**). α = 1.6 buys `reach` and destroys
`nmae`. The α chosen for v5 on magnitude-parity grounds lands on it.

**024 — the operator was under-sampled.** Capping the covariance basis at 6,000 perturbations
per donor was arbitrary. Using the full catalogue moves macro PDS 0.7995 → **0.8078** and
pushes optimal γ from 0.65 to **0.85**. That coupling is itself evidence the operator is real:
a better-conditioned covariance can be inverted harder before noise in the small eigenvalues
dominates.

Two refinements that did **not** pay: per-donor whitening (0.7376, −0.062 — after each donor
is mapped by a different operator the predictions no longer share a geometry, and averaging
mixes incompatible spaces); percentile flooring of eigenvalues (+0.0016, inside noise).

The spectrum explains why the ridge never mattered: top eigenvalue 3.90, median 0.0036,
condition number 23,000, and a ridge of 0.01 — *above the median* — changes the answer by
0.0003. **The gain comes from suppressing the dominant directions, not from differentiating
among the small ones.**

**025 — a bug in what v5 shipped.** The operator covers genes *every* donor carries (7,312 of
18,533); the model predicts on genes *any* donor carries. Whitening grows the on-axis block
~4.6×, and rescaling the whole matrix afterwards silently down-weights everything off-axis by
the same factor — 1,852 of the 7,989 genes the scorer tests at CPM > 20. Rescaling only the
on-axis block back to its own prior norm is worth **+0.0061 macro PDS**, and the fold with no
off-axis genes is bit-identical, which is the control the test needed.

General form: *a transform applied to part of a feature space must have its magnitude restored
on that part, not globally.*

**026 — a first-principles correction that loses.** `log1p` is linear below ~20 CPM and
logarithmic above, so a donor's `delta` carries **that line's** compression: the local factor
`x/(1+x)` is 0.20 at 5 CPM and 0.83 at 100 CPM, and the crossover sits in the middle of the
tested gene universe. Converting to a fold change at the donor's level and back at the
target's is a change of variables, not a heuristic — and it loses (0.7937 vs 0.7978).

Dividing by `expm1(b_ctrl)` to recover a fold change amplifies noise exactly where the
compression was largest. **`log1p`'s compression is variance-stabilising, and the bias it
introduces is cheaper than the variance the correction releases.** Worth remembering the next
time a correction to the measurement model looks obviously right.

Also the **third** independent failure for the target's basal control profile — after donor
weighting (010) and per-gene modulation (017), now as the parameter of a principled change of
variables. It defines the gene universe and nothing else.

**v6 built:** equal donors, whitening on the **full 45,531-perturbation basis** at **γ = 0.85**
with `--whiten-renorm preserve`, **CPM gate 5**, α = 1.0796. Cumulative LODO: 0.7378 → 0.8078,
**+0.0700 macro PDS**.

---

## 2026-08-21 · Whitening lowers the cosine while raising `pds` — the cleanest statement of what the metric wants

**Expected:** that whitening improves alignment, since it improves the ranking.

**Observed**, measured in the raw log-pseudobulk space the scorer reads, on the same LODO
folds:

| fold | cos, gamma = 0 | cos, gamma = 0.85 |
|---|---|---|
| k562_gwps | 0.0376 | 0.0312 |
| orion_hct116 | 0.0465 | 0.0302 |
| orion_hek293t | 0.0464 | 0.0270 |
| **macro** | **0.0435** | **0.0295 (−32%)** |

while macro PDS goes 0.7378 → 0.8078, **+9.5%**.

**Whitening makes our predictions LESS similar to the truth and MORE discriminating.** It is
not improving the estimate at all — it is reallocating similarity so that a prediction
resembles its own true effect more than it resembles the other 299.

**What changes as a result.**

1. **Mean cosine is definitively the wrong objective.** This was already suspected (exp 010:
   "+20% cosine produced −0.019 pds") but this is the clean version — a −32% cosine move
   accompanying the largest score gain in the project. Any harness tuned on cosine would have
   rejected whitening outright. Every cosine-based number in entries before 2026-08-21 should
   be read as measuring the wrong thing.
2. **`mse` moves further out of reach, not closer.** `min mse = 1 − cos²` goes 0.9981 → 0.9991
   against a 0.989 baseline. The rule stands: a clamped metric must not drive a shared
   parameter.
3. **It predicts the `nmae` cost, and the count-level run confirms it.** `nmae` rewards
   alignment, so lower cosine should cost it: measured 0.9397 → 0.9531 (exp 023), −0.022 on
   the official scale. Paid willingly for +0.10 to +0.14 scaled `pds`.

There is a general lesson here about the whole field's framing. Everyone — us included until
today — has been optimising "make the predicted effect closer to the true effect". The metric
that decides 5/6 of every score on the board does not ask that question. It asks whether the
prediction is *nearer its own target than to 299 decoys*, and those two objectives can move in
opposite directions by 30% at a time.

---

## 2026-08-21 · Experiments 027–031 — a control channel that works, a lever I mis-signed, and the harness that should have existed on day one

**027 — the per-gene emission variance channel, finally testable.**

Experiment 013 proposed it; experiment 023 supplied the missing half by showing that after
whitening, `|predicted effect|` is a worse per-gene confidence signal while the un-whitened
cross-donor agreement holds up. So the two channels want *different objects* and do not
conflict: whitened effects set the mean, un-whitened agreement sets the variance.

**Version 1 was a switch, not a dial, and it was a wash.** Sending untrusted genes all the way
to gamma-Poisson moved `reach` +31% and gave it straight back in `fid` and `jac` (≈ −0.002
overall). Full Poisson noise does not *reorder* a gene, it removes it from the called set.

**Version 2 works.** `mode="quasi"` splits each gene's panel total into an evenly-allocated
part and a randomly scattered part, so `spread` is literally the fraction of Poisson noise and
**the gene's total is preserved** — the pseudobulk, and therefore `pds`/`mse`/`nmae`, cannot
move by construction rather than by approximation.

At k = 1,000 trusted genes, spread 0.15, three seeds: `reach` 0.1517 → **0.1764 ± 0.0072**
(+0.0247, **3.4σ**), `fid` 0.5108 → 0.5128 (up in all three replicates), `nmae` flat, `jac`
−0.0034. Official scale: **+0.004 overall**. The even emitter's own seed noise is 0.0005, so
all of the variability is on the new arm, which is what a sampling change should look like.

**029/030 — a lever whose sign I got backwards, and what it actually is.**

Hypothesis: `nmae` rewards restraint (predicting zero on a gene contributes exactly `|true|`;
predicting confidently wrong contributes more), so shrinking low-confidence genes should buy
`nmae` and cost `pds`. **It did the opposite** — macro PDS 0.8017 → 0.8405 and `nmae` got
worse.

Experiment 030 identified why. The confidence signal is the geometric mean of the donors'
magnitudes and the prediction is roughly their arithmetic mean, so `P·c^λ` is close to
`sign(P)·|P|^(1+λ)`. Testing them separately, the **self** power transform — which uses no
cross-donor information whatsoever — matches or beats the consensus version at every exponent,
and adding a sign-agreement gate is worst of the three. **It is magnitude sharpening, not
consensus.** `gene_norm` and whitening between them leave the predicted magnitude distribution
flatter than the truth's, and one exponent fixes it.

`pds` gain at q = 2: **+0.038 macro PDS**, on top of whitening's +0.070.

**The cost is real and the instrument for measuring it is not.** At our actual emitted
magnitude (α = 0.92 of true, inverted from v3's board `mse`), proxy `nmae` goes 0.980 → 1.047
at q = 2, giving back roughly what `pds` gains. `pds` is exactly α-invariant and `nmae` is not,
so in principle you sharpen for one and re-tune α for the other — but that argument needs the
`nmae`-optimal α, and **the proxy puts it at ~0.75 while the count-level run measures 1.08**.
The two disagree on the one number the decision turns on. Not shipped.

**031 — the harness that should have existed on day one.**

Every count-level number in this project came from `h1-fast`: 40 perturbations of **H1 hESC**
against the **2025** panel. H1 is embryonic where the contexts are somatic — a +20% transfer
cosine on H1 produced a *negative* `pds` move on A/B/C — and the 2025 panel is
common-essential enriched, exactly the class the 2026 panel excludes. So every DE-metric
conclusion so far was measured on the wrong biology *and* the wrong genes.

`scripts/build_orion_harness.py` builds `data/harness/orion-hct116.h5ad`: **132 real VCC-panel
targets in HCT116** (colorectal, somatic) with 18,400 control cells. HCT116 is one of our
donors, so it is **dropped from the donor set** — predictions come from K562-gwps +
Orion-HEK293T. Leave-one-donor-out at count level, on the right biology and the right genes.

Two build notes worth keeping. Orion runs at ~154 cells per perturbation, so the reference is
built at **150** rather than the competition's 400 — the first attempt demanded 400 and
qualified 9 targets out of 263. And scoring reads **raw** values from `agg_results.csv` and
rescales with the official anchors, per the translation rule already in this log: no local
replicate bundle, which is what turned `h1-test` into a six-hour dead end (still running, still
nothing written, now killed by supersession rather than by finishing).

---

## 2026-08-21 · Experiment 030b/c/d — killing two explanations, and a principle that survives both

I wrote down a mechanism for the sharpening gain and then measured it. It was wrong, twice.

**Explanation 1, distribution matching: dead.** The claim was that `gene_norm` and whitening
leave the predicted magnitude distribution flatter than the truth's. Measured: sd of
log|effect| is 1.1377 for the whitened prediction and 1.1955 for the truth, so the *matching*
exponent is **1.05**. The PDS-optimal exponent is **2.0–2.5**. At q = 2 the prediction is far
spikier than the truth and PDS improves anyway.

**Explanation 2, selection: also dead.** If the transform works by suppressing the unreliable
tail of small per-gene entries, hard top-k should do at least as well — and it is a live
strategy, `Jessica A` runs `K562 top128 transfer v1` and holds the only raw `mse` below 1 on
the board. Measured: top-128 → 0.8155, top-1,000 → 0.8021, top-2,000 → 0.7912 (*below* no
transform at all), against **0.8409 for the smooth q = 2.5**. Combining them is worse than the
smooth transform alone.

**What it actually is:** multiplicative shrinkage whose strength decreases with the observed
magnitude, `|e| -> |e| · |e|^(q-1)`. If per-gene noise is roughly additive and comparable
across genes, small entries have the worst signal-to-noise and should be shrunk harder. The
transform keeps every gene's weak-but-real contribution instead of deleting it.

**The principle, now with two independent instances: smooth reweighting beats hard truncation,
in every basis tried.** Across perturbations, SVD truncation loses monotonically (exp 004)
while smooth spectral flattening gains +0.070. Within a perturbation, top-k selection gains
little or loses while smooth power sharpening gains +0.038.

That also dissolves what looked like a contradiction in this notebook. Experiment 004 said
"the low-variance tail carries signal, do not truncate"; this says "the small-magnitude tail is
noise, suppress it". Both hold — a low-variance *SVD direction* is a coherent pattern shared
across perturbations, a small *per-gene entry* is one number in one prediction.

**030d — the two transforms compose.**

| gamma \ q | 1.0 | 1.5 | 2.0 | 2.5 |
|---|---|---|---|---|
| 0.00 | 0.7210 | 0.7578 | 0.7826 | 0.7932 |
| 0.65 | 0.8017 | 0.8293 | 0.8399 | **0.8415** |
| 0.85 | 0.8011 | 0.8287 | 0.8385 | 0.8409 |
| 1.05 | 0.7954 | 0.8234 | 0.8341 | 0.8363 |

Whitening alone +0.080, sharpening alone +0.072, **both +0.120** — about 79% additive rather
than redundant, which says they are reweighting genuinely different things (directions across
perturbations vs magnitudes within one). The optimal gamma drifts down 0.85 → 0.65 as q rises,
which is the small overlap showing.

**Macro PDS 0.7378 → 0.8415, +0.104 from the incumbent.** Whether any of the sharpening half
survives its `nmae` cost is a count-level question the proxy cannot answer — it puts the
`nmae`-optimal alpha at 0.75 where the real metric measures 1.08. The q = 1.0 / 1.5 / 2.0 arms
of experiment 031 are running on the somatic harness for exactly this.

---

## 2026-08-22 · LEADERBOARD — v6 scores 0.1369, rank 1 at scoring time, rank 3 on arrival

`Attention is not all you needed` — v6, equal donors + whitening γ=0.85 on the full basis +
preserve renorm + CPM gate 5, α=1.0796.

| | v3 (previous) | **v6** | change |
|---|---|---|---|
| **Overall** | 0.0825 | **0.13689** | **+0.0544 (+66%)** |
| pds raw / scaled | 0.7279 / 0.504 | **0.8163 / 0.701** | **+0.0884 raw** |
| nmae | −0.020 | +0.065 | +0.085 |
| reach | 0.055 | 0.090 | +0.035 |
| fid | −0.047 | −0.037 | +0.010 |
| jac | 0.004 | 0.002 | −0.002 |
| mse | 0 | 0 | — |

**We have the best raw `pds` on the board: 0.8163**, against Ibrahim Mansour 0.8087 and Vivai
0.7765. Every metric except `jac` improved, from a change that only touched the *direction*.

Rank 1 at the moment of scoring; two teams submitted during our four-minute upload and the
board settled at **rank 3 of 96** (Ibrahim Mansour 0.1427, Vivai 0.1380, us 0.1369). The field
moved 0.009 in one day at the top — Ibrahim went 0.0716 → 0.1427.

**The LODO harness translated well.** LODO 0.7378 → 0.8078 predicted a raw `pds` move; the
board delivered 0.7279 → 0.8163. The somatic count-level harness (exp 031) had already
predicted the direction and rough size on real metrics, and did so *before* the window opened.

---

## 2026-08-22 · `mse` is NOT structurally unreachable — and the reason is our own emitter

This repo has said since 2026-08-21 that `expr_mse` is "closed, permanently, for everyone",
because `min_α mse = 1 − cos²` needs cos > 0.105 and we measure ~0.08. **Two teams now score
positive on it**: Ibrahim Mansour raw 0.9187 (**+0.0728 scaled**) and Immovable Myhal Rock raw
0.9585 (+0.0308). Ibrahim's `mse` alone is 0.012 of overall — twice our margin behind them.

The argument was wrong because it ignored a term. Reading the spec again:

```
N_p = (1/G_p)( ||b̂_p − b_p||² − ρ·min(Ĉ_p, C_p) − C_p )
```

`Ĉ_p` is the **delete-one jackknife estimate of our own submission's sampling noise**, and the
metric *subtracts* `min(Ĉ_p, C_p)` from the numerator — a credit for noise we emit, capped by
the reference's own noise and by `ρ`, which is bounded by our across-perturbation spread.

**Our even emitter is minimum-variance, so `Ĉ_p ≈ 0` and we forfeit the entire credit.** Every
other team emitting realistic cell-to-cell noise collects `min(Ĉ_p, C_p) ≈ C_p`. That is not a
modelling gap, it is a term we opted out of when we chose the emitter for `fid` coverage.

`1 − cos²` is the floor for the *uncorrected* squared error. The scored quantity has a second
term we control directly, and the cosine bound never applied to it.

**The trade is real and already measured.** Experiment 005 found NB emission craters `fid`
(−0.77 to −0.86 scaled) because `n_pred` collapses. But `mode="quasi"` (exp 027) sets emitted
variance anywhere between the integer floor and full Poisson *while preserving each gene's
panel total*, so it can buy the jackknife credit without moving the mean profile. A global
`spread` sweep against `mse` and the DE metrics is the obvious next experiment, and it is the
first thing that could take `mse` off zero.

---

## 2026-08-22 · CORRECTION: jackknife credit is de-biasing, not a bonus

The overnight mse×spread sweep (exp 033 correction, exp 034) **refuted** the practical half of
the previous entry. Emitting more noise makes raw mse *worse* (2.26 → 2.53). `||b̂−b||²`
already contains Ĉ; subtracting ρ·min(Ĉ,C) cancels it when ρ=1 and under-cancels when ρ<1.
Minimum-variance is optimal-or-tied for mse. Ibrahim's +0.073 is a better cosine (raw mse
0.919), not an emitter trick. Giovanni's unperturbed mse 2.02 is ρ=0.

Quasi spread **survives as a DE lever** (reach/fid), not an mse lever. Exact-total repair
makes gene-panel sums bit-identical to even (verified on v6 vs rebuilt v8, 8 perts, max|Δ|=0).

**Do not submit v8 until `out_som_spread` confirms pds is unmoved on somatic data.** Rank uses
the latest entry; one slot left today.

---

## 2026-08-22 · Board operating points (113 teams)

Us 0.1369 rank 3, best pds (0.701 scaled). Gap to Ibrahim 0.006, to Vivai 0.001. Best-in-class
sum ≈ 0.18. Field is on one α frontier: larry mse~1 / dead reach, us mse 2.1 / pds lead,
Vivai mse 3.9 / best nmae, Hongjian mse 16 / best reach and wrecked nmae. Public names that
leak method: CSBL emission variant, Sidechain unshrunk+minvar, Hongjian gwps scale 1.6,
Life Code hybrid scale 1.5 + STRING fill. Don't copy; they confirm emission and magnitude
are the live surface, which we already knew.

---

## 2026-08-27 · The instrument was broken, and it had been broken for six days

Session summary. Rank at start: **11 of 363**, live entry `v15` `0.1494363`, best held
`v11` `0.1501827`. Leader Aginglab.com `0.18991`.

### 1. The scoring function, solved on 780 board entries

Union of every leaderboard snapshot: 780 distinct entries, 363 teams. `score_avg` is the
plain mean of the six scaled members (max |diff| 3.3e-16), and `scaled = (raw − b)/(r − b)`
fits four of them at `R² > 0.9999`:

| member | `b` | `r` | span | raw Δ for +0.01 overall |
| --- | ---: | ---: | ---: | ---: |
| `pds` | 0.50020 | 0.95228 | 0.452 | +0.0271 |
| `fid` | 0.51239 | 0.80964 | **0.297** | **+0.0178** |
| `reach` | 0.07902 | 0.96414 | 0.885 | +0.0531 |
| `jac` | 0.03038 | 0.39727 | 0.367 | +0.0220 |

`score_mse` and `score_nmae` are per-perturbation normalized then averaged — not linear in
the reported aggregate. Do not do arithmetic on them.

### 2. Amplitude is retired

STRATEGY Addendum 5 argued alpha should go **up**, and a v13 (α=2.0) package was built for
it. The board says otherwise: across **288 within-team consecutive submission pairs**,
`Δscore_avg` is uncorrelated with `Δlog mse` (r = −0.07, **p = 0.24**). Between teams,
amplitude explains **9.7%** of reach variance and our residual is −0.0056 (z = −0.15) — we
are already on the line. **We527** scores reach `0.2529` at raw mse `1.3577`, within 2.6% of
ours, against our `0.1815`; **Aginglab** halved amplitude and gained `+0.100` reach. v13
stays unsubmitted.

Every other scalar knob is also at or adjacent to its optimum: `min_cpm` is a sharp peak at
10 (8 → +0.0028, 15 → −0.0056), global `spread` plateaus at 0.20–0.30, call-set truncation
is negative on the board (v11 > v15 > v14), and `mse` costs 0.387 `pds` to buy 0.034.

### 3. The harness has been lying, and here is exactly how

`data/harness/orion-hct116.h5ad` has **control cells 38% deeper than its perturbation
cells** — median UMI 16,798 against 12,178, ratio **0.725**, in **132/132 perturbations**.
`build_orion_harness.py` filled controls greedily first-come-first-served, so all 18,400
came from the first ~20 batches (median NTC UMI 24,850) while perturbation cells
accumulated across all 109 (16,589 and 18,007 later). The shipped build log says it
outright: `20/109 | 0 perturbations full | 18,400/18,400 controls`.

The Wilcoxon is a rank test and CPM does not fix the per-cell zero fraction, so uniformly
shallower perturbed cells read as uniformly **down**. The reference's DE set came out
**90.8% negative**, against `0.538` on h1-test, `0.497` on Replogle K562 and `0.553` on
RPE1. Binomial thinning proves causation: invert the depth mismatch and the skew inverts to
`0.275`; equalize it and it lands at `0.541`.

**This is the missing explanation for v14.** Experiment 093 measured `+0.0202` avg6 and the
board returned `−0.00291` with `fid` at the wrong sign; its BOARD-RESULT blamed
cells-per-perturbation. The real cause is that the harness's `n_real`, its precision curve
and its direction distribution are all wrong together.

Repairs measured:

| harness | depth ratio | zero-fraction gap | frac. reference-DE negative | `n_real` |
| --- | ---: | ---: | ---: | ---: |
| greedy (shipped) | 0.725 | +0.0340 | 0.897 | 322 |
| `--sampling reservoir` | 0.987 | +0.0012 | 0.873 | **3** |
| **`--thin-controls`** | **1.000** | **−0.0015** | **0.630** | 74 |
| h1-test | 0.985 | +0.0019 | 0.500 | 222 |

Reservoir sampling fixes depth and destroys power. Thinning is the repair, and it is a
drop-in replacement — same targets, gene axis and shape, only the control counts differ.
A residual `0.630` survives that is **not** depth (unmatched batch effects), so h1-test
remains the instrument for absolute DE numbers. `scripts/check_harness_depth.py` now gates
any harness in one line.

**The donor tables are clean** — `build_orion_donor.py` has none of the greedy constructs,
measured depth ratio 1.0031, and the panel perturbations project onto the pure depth axis at
cos −0.013 with the IQR straddling zero. The model is not contaminated.

### 4. A fix I designed, tested, and threw away

Measuring `reach` on the HCT116 tables, our confidence ordering looked *anti*-correlated
with sign correctness (`Spearman +0.50`; genes we call significant right 34% of the time,
genes we do not right 83%). I built an arm around it. Then I checked whether it transferred:
h1-test's reference is **53.8% negative, not 90.8%**, so the "inversion" was our up-first
ordering multiplied by one screen's broken direction distribution. The arm was stopped
mid-build and never submitted; the retraction is in `experiments/094-reach-ordering/`.

What survives is narrower and real: the order the scorer walks is a deterministic,
pre-computable function of `|F̄_g − ½|/√(tc_g)`, and `|predicted effect|` — the key
`--keep-top` ranks on — correlates with it at **+0.071**. So that flag never had purchase
on the thing it was meant to steer.

### 5. Engineering

- **The build loop went from 3h20m to ~18 minutes**, and an unflagged build of the new path
  is **byte-identical to `submissions/v11.vcc`**. The old runner's three hours were all
  re-downloading inputs around a build whose own timer reads 16:55.
- `WHITENER_REGISTRY` + `--whitener-id` replaces a frozen-hash pin that made any operator
  A/B impossible; it still fails closed on an unregistered id or a one-byte drift.
- `--keep-spread` decouples the two ends of the emission spread channel.
- 20/20 tests pass.

### 6. Submissions

Slot 1: **`a2-prune`** — v11 with experiment 091's pruned whitening basis
(`whitener_id: union_pruned`), the one measured-positive unshipped result we held. Every
other parameter identical to v11, coverage 300/300. It had been blocked by the frozen-hash
guard. The operator is not a marginal change: λ₂ falls **1.138 → 0.497 (−56%)**, so the
common-essential perturbations were contributing a strong secondary axis that whitening was
spending its budget suppressing.

### 7. Late addition — the counterfactuals were all wrong, and over-dispersion is the lever

Six hours in, with two ordering arms built and one about to be submitted, I built the check
that should have existed from the start: `experiments/095-clean-harness/realize.py` drives
the **unmodified production `emit_counts`** with a policy's spread vector, runs the exact DE
path, and scores what the scorer would see — with a **same-policy second seed as a noise
control** (`5e-4`).

Every counterfactual evaporated. The K=2000 block promotion measured `+0.00583` reordering a
table by fiat and **−0.00021** when actually emitted. The cross-donor agreement key, the best
ordering found all day, realises at **−0.00193**. Nine policies, spanning both mechanisms and
the full useful range of K: everything that leaves `n_pred` near production is at or below
zero, and the three positives all shrink the call set 15–40%, a trade the board has already
priced twice.

The reason is mechanical: predicted p-values span **hundreds of orders of magnitude** and
`spread` moves them by a few. **`spread` changes which genes are called long before it
changes their order.** That retires the central hope of Addendum 7, written five hours
earlier the same day.

**But it pointed at the thing that does work.** If the problem is that Poisson noise is too
weak, use noise that is not capped at Poisson. Gamma-Poisson scatter with the same
exact-total repair, realised through the same path:

| policy | `n_pred` | reach | Δreach | Δfid | ΔTOTAL |
| --- | ---: | ---: | ---: | ---: | ---: |
| production | 3,934 | 0.0838 | — | — | — |
| **spread 0.25, φ=2** | **3,218** | **0.1378** | **+0.0102** | +0.0100 | **+0.0204** |
| spread 1.0, φ=0.30 | 692 | 0.1961 | +0.0212 | +0.0167 | +0.0448 |

`nmae` is invariant to four decimals across every row — the exact-total repair holding.

It is also the fix the *data* asked for. Audited against the real challenge controls, our
per-gene CV across 400 emitted cells is **0.38–0.40** where the real control cells are
**0.84–0.92**: we have been emitting cells **~2.3× too uniform** since v6, which is what made
the rank test see an artificially confident, up-biased signal in the first place.

The big-φ rows are not bankable — their `fid` gains are the harness's coverage artifact and
would invert on the board as v14's did. `spread 0.25, φ=2` is the operating point because it
buys `+64%` reach for an **18%** call-set cut, landing near `n_pred` 4,880 on the challenge
where coverage is still ~0.99. Board-anchored projection: **≈ +0.016 overall**.

**Arm `a7-nb-phi2`** — `a2-prune` plus `--tail-spread 0.25 --overdispersion 2.0`, one
variable from our best board artefact — is the second submission.

### 8. The mechanism, nailed down

Three independent measurements say the same thing about *why* over-dispersion works, and
together they are the reason to expect it to transfer.

**It aligns the p-order with our own model.** `Spearman(p_adj, |predicted effect|)` rises
monotonically with φ: `0.060` at production, `0.170` at φ=2, `0.223` at φ=4, `0.275` at φ=8.
Production's `0.060` matches the `+0.015` measured on the real challenge artefact and the
`+0.07` on the HCT116 tables — **our confidence ordering has been very nearly independent of
our own model's magnitudes**, which is the whole defect in one number.

**It makes the head genuinely purer, at every depth.** Sign purity of the first *k*
adjudicable genes of the reference's DE pool:

| depth | production | φ=4 |
| ---: | ---: | ---: |
| 1 | 0.885 | **1.000** |
| 10 | 0.696 | 0.827 |
| 20 | 0.619 | **0.775** |
| 100 | 0.539 | 0.635 |

`+0.10` to `+0.16` everywhere, not at one depth. That rules out a lucky run at the front,
which is the failure mode a metric with a 0.9 pass mark is most exposed to.

**And it is not simply "more noise is better".** φ=8 has the *best* alignment (`0.275`) and
*worse* reach (`0.109` against φ=4's `0.168`). Past φ≈4 the scatter drowns the signal faster
than the improved ordering recovers it — the same non-monotonicity experiment 093 §7 saw and
could not interpret on a confounded harness.

**What it is not.** Over-dispersion captures **14%** of the ordering ceiling: reach goes
`0.084 → 0.161` against an oracle of `0.625`, leaving `+0.0875` overall between here and a
perfect confidence signal. The oracle is invariant to the emitter — it depends only on which
genes we get the *sign* right on — so the rest is a modelling problem, not an emission one.
**The emitter can now express a good ordering; it cannot invent one.**

And realism is not the objective: the dispersion-matching φ at the median gene is `2.53`
against a measured optimum of `4.0`, and a per-gene φ fitted to each gene's own control
dispersion **loses** to a single global φ.

### 9. Both slots landed. `0.1543607`, rank 7, and one uncomfortable calibration

| | overall | Δ | rank |
| --- | ---: | ---: | ---: |
| v11 (previous best) | 0.1501827 | — | — |
| `a2-prune` — exp091 pruned whitening basis | 0.1503084 | +0.0001257 | 11 |
| **`a8-nb-phi4` — over-dispersion φ=4** | **0.1543607** | **+0.0040523** | **7 / 363** |

`a8`: `pds` **+0.000000 exactly**, `nmae` `−0.000038`, `fid` `+0.006393`,
`reach` **`+0.021398`**, `jac` `−0.003439`.

**The `pds` zero is the result to keep.** The exact-total repair survives gamma-Poisson
scatter at φ=4 — a change that quadruples per-cell variance — on an independent build
against a reference we have never seen. The emission channel is *provably* free on
`pds`/`mse`/`nmae`, board-verified twice now.

**And `reach +0.0214` is the largest single-member gain in this project's history**, against
a previous best of `+0.0137`.

#### The prediction was right in sign and wrong by 5× in size

Pre-registered `+0.0213`, delivered `+0.0041`. The harness over-stated `reach` by 2.9× and
`fid` by 6.1×, and **`jac` moved the wrong way** — `−0.0034` against a predicted `+0.0022`,
which is most of the shortfall on its own. Nothing anticipated it: `jac` *rose* with
over-dispersion on both harnesses. The coverage analysis only ever modelled `fid`.

#### Neither harness bracketed the truth, and that is the lesson

Running the same measurement on h1-test **after** submitting:

| harness | production reach | φ=4 reach | Δ |
| --- | ---: | ---: | ---: |
| orion-hct116-thinned | 0.0838 | 0.1605 | **+0.0767** |
| **the board** | **0.1701** | **0.1890** | **+0.0189** |
| h1-test | 0.1878 | 0.1830 | **−0.0009** *(noise −0.0004)* |

The gain scales inversely with how good the baseline ordering already is — which is a
coherent mechanism, since over-dispersion repairs a *badly ordered* emission. But it means
**the harness whose absolute reach is closest to the board predicted zero**, and the one
that predicted a large gain was 4× over. Either read alone was wrong.

**Run every emission arm on both harnesses. Their disagreement is the error bar.** That is
the third methodological rule this day produced, after "force it through the real emitter"
and "check the mechanism, not just the number".

---

## 2026-08-28 · The member we had written off was half our gap, and it cost one env var

**Where the day started.** Rank 16 of 385 at `0.1543607`. Overnight the field had passed us —
every team above us submitted on the 28th, we last submitted on the 27th.

### 1. The board publishes the scaled members, and we were never reading them

**Expected** to re-derive scaled scores from the anchors in
`experiments/089-production-v11/validate_authorization.py`. **Observed** that those anchors
reproduce `score_avg` under no clipping rule (max error 1.01 clipped, 358 unclipped) — and that
each board entry carries `score_pds`, `score_mse`, `score_nmae`, `score_fid`, `score_reach`,
`score_jac` directly, whose plain mean reproduces `score_avg` to **3.3e-16** across all 385
entries. **Changed:** all board comparisons now read the published fields. The board also
carries a `description` field that **131 of 385 teams populate with real method detail**,
including two teams' explicit metric-calibration runs. We had never read it.

### 2. The diagnosis: we lead the board on `pds` and score exactly 0 on `mse`

**Observed** our `score_pds` of `0.7366` is the **best of all 385 teams** — better than the
leader's `0.7089` — while `score_mse` is exactly `0`, as it is for **376 of 385 teams**. The
leader's `+0.1693` on that one member is half our entire gap to first.

`pds_cosine` is a cosine, and the rank-based DE members are magnitude-insensitive. We had
optimised hard against five members that cannot see magnitude, and the one that can said our
magnitudes were wrong.

### 3. It was not a modelling wall — it was forfeited noise credit

`expr_mse_unbiased_capped_norm` subtracts a jackknife estimate of sampling noise,
`rho*min(C_hat_p, C_p) + C_p`, from the squared error. **Observed:** our emitted cells are
**0.67x as dispersed as the reference's**, so we claimed `0.1668` of an allowed `0.2465`, and
our amplitude sat above the `rho` cliff. Written as `(A - 2X + Kc - rho*claim)/sumD`, four of
the six terms are computable from our own submission plus control cells we already hold.

**Two board anchors, from other teams' published descriptions, pinned the scale:** "Kun Lab"
submitted an explicit zero-effect floor and scored raw `1.0424`; the scored baseline is
`~0.9879`. So *predict nothing* is measurably worse than *predict the mean response*.

### 4. `a11` — one environment variable, `+0.0119`, rank 16 → 7

`VCC_ALPHA` 1.0796 → 0.838. **Predicted** raw `mse` `0.9021` and overall `+0.0089`;
**observed** raw `0.9100` (0.9% error) and **`+0.011931`**, `0.15418 → 0.16612`. `score_mse`
0 → `+0.0820`. The decomposition has now predicted two board results out of sample, a10's and
a11's, both to within 1%.

### 5. `a10` — the emission channel is exhausted

Shipped first, before the mse work landed: `spread` 0.25→0.15, `phi` 4→16, a **24% change in
the call set**, for **`-0.000177`** overall. Both harnesses had predicted a difference ~50x
larger. Nine experiments (088-095) have now gone into emission; it produced one real gain
(a8's `+0.0041`) and two consecutive results inside noise. **Closed.**

### 6. What I got wrong, in order

- **The on-target hypothesis.** Measured cross-donor transfer at `rho = 0.647` for the target
  gene's own knockdown against `0.10-0.19` for the downstream response, then estimated via a
  split-half de-noising that mean-response-plus-exact-on-target would score `mse ~ 0.910`.
  **Two board entries refute it by 800x:** "Lumos" (on-target only) reads `1.0425` against
  "Kun Lab"'s no-signal `1.0424`. The error: I de-noised both halves of a ratio the scorer
  de-noises on one side only.
- **The mid-course amendment on `pds`.** Three hours before submitting a11 I saw three points
  of an amplitude sweep, read a monotone rise in `pds`, and amended the prediction to expect
  `pds` to **fall ~0.012**, writing that the falsifier was "now expected to fire." It **rose by
  `+0.0033`**. The completed sweep shows `pds` saturating above `alpha = 0.6` — a11's move sat
  entirely inside the plateau. **The pre-registration beat the mid-course correction**, which
  is rule 2 demonstrated against myself.
- **The sweep would have been actively harmful to steer by.** Its `reach` column falls with
  amplitude; the board's rises. Its `avg6` peaks at `alpha = 0.60`, which would have cut
  amplitude far harder than 0.838 for `mse` already saturated at the clamp. Harness `n_real`
  is 92. Rule 1, again.

### 7. Method changes

- **A metric's *level* can be unreadable while a *factor* transfers.** `fid` factorises into a
  coverage leg carrying `n_real` and a precision leg that does not; `jac` does not factorise.
- **External signals transfer; self-referential ones do not.** Our directional signal is
  concentrated — the top quartile of perturbations carries **80-85%** of it and **12-29%**
  point the wrong way. Of the candidate confidence scores, both computed from our own
  prediction fail or flip sign across surfaces; `donor_agreement` is the only one positive on
  both, and the only one surviving partial correlation against effect size.
- **Agreement between two analyses is not a check if they share an omission.** Two independent
  derivations converged on shrinking amplitude ~60%; both had dropped the same two correction
  terms and the answer is ~6%. Every wrong answer in this session died to *executing*
  something, never to arguing.
- Retracted `RESEARCH-METHOD` §1's founding bracketing example: `realize.py` re-realises a
  *policy* against each harness's own controls and never scores the shipped `.vcc`. A
  prediction artifact is baseline-bound — scoring a8's `.vcc` against orion reads `mse = 343`,
  of which the control-profile gap explains 99.8%.

### 8. Where it leaves us

**Rank 7 at `0.16611503`.** `a12` (dispersion to the correction cap + `VCC_ALPHA` 1.012) is
built, probed and pre-registered at a predicted `0.1712`. Beyond it the honest ceiling for
calibration is raw `~0.870` against the leader's `0.8270`, whose aggregate directional cosine
is **at least 1.28x ours** — a late correction that reversed a more flattering earlier reading.
**Calibration buys rank 3-5. First place needs the direction.**

### 9. Addendum — my shrink-toward-the-mean idea was refuted 12/12, and the reason matters

**Expected:** our pipeline shrinks a low-confidence gene's effect toward ZERO, i.e. toward the
context's control profile. But the metric's scored 0-anchor is the context's *mean perturbation
response*, and the board pins both ends — "Kun Lab"'s explicit zero-signal submission reads raw
`1.0424` while the baseline is `~0.9879`. So predicting nothing is measurably worse than
predicting the mean response, and shrinking toward the grand mean of our own predicted effects
(James-Stein, `mbar + a*(e - mbar)`) should dominate shrinking toward zero. I sent this as a
correction to the specification while the implementation was in flight.

**Observed:** `zero` beats `mean` on **both harnesses, at every gamma, with both keys — 12 of
12**. On h1-test the mean target is worse than not shrinking at all (`+0.00286` vs `-0.00892`).

**Why, and this is the part worth keeping.** The argument had an unstated premise: that *our*
predicted panel mean approximates *the reference's* mean effect. It does not.
`cos(our grand mean, the reference's mean effect)` measures **+0.016 on h1-test and -0.066 on
orion** — no capture, and not even a consistent sign. Two independent scripts reproduced both
numbers. Second mechanism: the metric's `rho` budget is a **centred** sum of squares, so the
panel mean is exactly what it centres out; shrinking toward it inflates `(A - B)/A` from
`~0.012` to `~0.055`, and `A - B` is precisely the term the master formula drops.

**Changed:** the default is `zero`, against my own specification, with the reasoning and the
refutation recorded at the constant and in the test that asserts it.

This is the **third** hypothesis of mine the day killed — after the on-target estimate (board
refuted it 800x) and the mid-course `pds` amendment (wrong in sign). All three were
self-consistent, all three had a mechanism, and all three died to a measurement rather than to
an argument. The general form: **each one assumed a quantity we hold approximates a quantity
the scorer holds** — our de-noised residual for the scorer's ratio, a different model's
amplitude slope for ours, our panel mean for the reference's. That is worth checking directly
before building on it, because it is cheap to check and it has now been wrong three times.

**Also corrected:** `MSE-MECHANISM.md`'s intervention 5 (`1.17-1.38x`, `+0.0044`) was computed
through a `sqrt(W/Wg)` shortcut that ignores `A - B` — on the one estimator family where
`A - B` is exactly what moves. It is overstated, and the exact-identity numbers above supersede
it. A practical warning came with it: at `gamma = 1` the mean surviving amplitude on the real
300-target panel is **0.108**, so any shrinkage arm needs its amplitude re-solved (measured
`1.62x` / `1.80x`), and `VCC_ALPHA` cannot do it because PerturbMatch and CD4 are added at
fixed magnitude afterwards.

### 10. The shrinkage value, corrected three times, and the third was our own bug

The donor-agreement shrinkage has been valued four times in one evening. The sequence is worth
recording because only the last correction came from code rather than from analysis:

| valuation | raw `mse` | overall vs a12 | why it was wrong |
| --- | ---: | ---: | --- |
| published `+0.0044` | — | +0.0044 | scored through a `sqrt(W/Wg)` shortcut assuming `A ~ B` |
| corrected #1 | 0.8957 | **−0.0030** | that shortcut, done exactly: the family was net negative |
| corrected #2 | 0.8203 | +0.0102 | right target (`zero`), but a live bug inflated it |
| **corrected #3** | **0.8569** | **+0.0055** | blank-row artefact removed, board solve re-anchored |

**The bug.** `agreement_shrink_vector` ranks genes by donor agreement, descending, ties broken
by gene index. For a row whose key is **all zero** — no two donors agreeing about anything —
every gene ties, so the tie-break alone sets the order and the function returns
`exp(-gamma * g / 2000)` indexed by **position in the gene axis**: gene 0 keeps full amplitude,
gene 18,000 is annihilated, for no reason but where it sits in the array. **111 of 900 rows on
the real challenge panel.** Those rows were shrunk to nearly nothing *and scored as though it
helped*, which is what produced the `+0.0102`.

The fix is structural rather than a special case: a row with fewer than two distinct key values
has no ordering to express and must return `alpha = 1` — which also covers a constant non-zero
row, a case nobody had counted. The invariant to assert is that **a degenerate key is
bit-identical to the shrinkage being off**, which is a stronger test than checking the all-zero
row specifically.

**What this says about the method.** Three analyses agreed on `zero` over `mean` 12/12, and all
three were reading the same defective vector. The concordance was real about the *ordering* and
uninformative about the *level*, exactly as the `claim`-not-re-emitted caveat had already
warned for a different reason. **Two independent "agreement is not a check" failures, on the
same family, within an hour.** The generalisable form is narrower than "agents can share an
omission": *any* number of analyses reading the same intermediate artifact are one measurement,
however independent their algebra.

**Also corrected:** the board solve is re-anchored on a8+a11 rather than the older v11+v12 pair
— our deconfounded cosine moves `0.0741 → 0.0701`, and a12's predicted raw `mse` moves
`0.8785 → 0.8882`. And `mse-mechanism` withdrew its own `n_real ~ 4757` coverage framing on the
evidence from a11 that coverage is not binding in this range.

**Operating point for a13**, if it ships: `gamma 0.5` (the mse surface is flat across 0.5-1.5 —
total range 0.0022 raw — and 0.5 keeps 62% of the effect norm against 54%), amplitude carried
by `--effect-scale 1.4025` **not** `VCC_ALPHA`, because `VCC_ALPHA` is applied pre-whitening
with PerturbMatch and CD4 added at fixed magnitude afterwards, so using it would change the
donor:PM/CD4 ratio rather than the amplitude. Predicted `+0.0078` overall against a11, with
`n_pred` *rising* to 5,152. Falsifier that abandons the family: `score_pds` down more than
0.03.

### 11. The whitener: the largest lead of the day, closed as a negative

`experiments/097-direction` measured that dropping the whitening operator (`gamma 0.65 -> 0`)
is worth **1.495x / 2.136x** on directional accuracy — monotone on both surfaces and the
largest single number the study found. It was correctly *not* recommended, because
`WHITEN_GAMMA` was fitted for `pds`, the member we lead all 396 teams on, and nothing then
existed that could price the trade.

Gated across all six members with re-seeded twins, on top of the arm that would actually ship:

| `gamma` | h1-test net | orion net |
| ---: | ---: | ---: |
| 0 | −0.0033 | −0.0080 |
| 0.30 | +0.0027 | −0.0032 |
| 0.45 | +0.0023 | −0.0001 |
| **0.65 (production)** | — | — |
| 0.85 | −0.0021 | −0.0002 |

**No `gamma` is net positive on both surfaces.** Every setting below production gains on
h1-test and loses on orion. The exchange rate is brutal and one-sided: **`+0.01` of `kappa`
costs `0.0181` raw `pds` on h1-test and `0.0517` on orion**, so `pds` alone more than consumes
the directional gain before any DE member counts. Production's 0.65 is at or beside the
two-surface optimum.

**This was a deliberate re-opening of a closed lever**, licensed by §9's asymmetry note — a
lever closed under one objective is not closed under a new one. It cost one measurement and
returned a null. That is the rule working, not failing.

**Third time in one evening that the `mse` surface pointed at an arm the six-member gate then
refused** — the shrinkage gamma, the per-perturbation gamma, and now the whitener. At the
aggressive end of each family the partial objective is not merely incomplete but
*anti-correlated* with the total, so optimising it ungated actively selects the loser.

**Consequence for the campaign.** The direction cannot be bought from the operator. What
survives of the whole direction study is the donor-agreement family alone — a13, `+0.0215 /
+0.0040` — against a per-perturbation oracle ceiling of `1.97x` that is mostly unbuildable.
Everything else measured was either negative, unreadable, or already at its optimum. **First
place needs a source of directional signal this project does not currently have.**

One caveat recorded rather than buried: the study's `kappa` reconstruction matches the shipped
prediction at per-row cosine 0.99922 / 0.99799 but reproduces `kappa` at 1.01x / **1.19x** the
shipped value, because the aggregate cosine is ~2% of the effect's norm and a 0.2% angular
error is 10% of the signal. Every number above is a *ratio within one reconstruction* and the
`pds` column is measured directly, so the verdict holds — but those `kappa` levels must not be
reused elsewhere.

---

## 2026-08-29 · `nmae` is a direction member, and 83% of the gap is one thing

**Board.** a12 `0.16798741` (+0.0019, rank 8); a13 `0.1476855` (**−0.0203**, the largest single-arm
loss recorded). Leader `0.211224`.

### 1. a13 — the mse lever is real, the family is closed

Raw `mse` **0.9013 → 0.8493**, a `−0.052` move against a `−0.013` model — **4x the prediction**,
and within `0.022` of the leader on that member. And it lost anyway: `nmae` reached raw
**1.0185**, past its own baseline anchor (~1.0013), i.e. **worse than predicting nothing**.
`pds` `−0.0288`, `reach` `−0.0197`.

**Both harnesses under-called `nmae` by 8x** after rule 1's "take the worse surface" was already
applied. Even granting the full `mse` gain and zeroing the `nmae` cost, the rest still nets
`+0.0003` — inside noise. **The donor-agreement shrinkage family is closed**: it buys `mse` and
pays it back across `pds`/`nmae`/`reach` at about one-to-one.

### 2. The finding: `nmae` is not a magnitude member

Read exactly, `de_wilcoxon_lfc_nmae`'s denominator and gene gate are **reference-only** — the
sole submission-side quantity is `lfc_hat`. There is no normalisation lever; route closed by
reading the spec.

And our `|lfc_hat|` is **0.138x / 0.168x** the reference's on its own gate, with `|h| <= |r|` on
**99.3%** of gated genes. On that side `|h - r| = |r| - sign(r)*h` identically, so

> **`NMAE = 1 - <sign(lfc_ref), lfc_hat>`**  (residual +0.0008 / +0.0036)

**In our regime `nmae` is a direction member wearing a magnitude member's name.** The
reference's magnitude enters only through the normaliser.

### 3. What that does to the gap

| member | leader | us (a12) | gap |
| --- | ---: | ---: | ---: |
| `nmae` | 0.1838 | 0.0940 | +0.0898 |
| `reach` | 0.2066 | 0.1168 | +0.0898 |
| `fid` | 0.0030 | −0.0315 | +0.0344 |
| `mse` | 0.1693 | 0.0912 | +0.0782 |
| `pds` | 0.7089 | **0.7379** | −0.0290 *(we lead all 439)* |

`nmae + reach + fid` is **82.6% of the entire gap to first place, and all three are the same
quantity — per-gene directional accuracy.** `mse` is the remaining 30% and is now a solved,
measured, movable member. **There is one problem, not four.**

### 4. A measurement error that invalidated every prior `nmae` reading

**Every `nmae` number this project has taken was on a gene gate 2-20x shallower than the
board's**, and the amplitude response *flips sign* across that range. The challenge gates ~2,000
genes; our harnesses call 89 (orion) and 1,120 (h1-test) at BH 0.05. The board has measured
`d(raw nmae)/d(ln amplitude)` twice at **+0.0062 and +0.0097**; both harnesses read
**−0.018 / −0.051** at BH 0.05 — the wrong sign. h1-test reproduces the board only at BH
0.35-0.45; **orion never does, at any depth.**

So at challenge depth `nmae` wants amplitude **down**, the same direction `mse` does — and `mse`
charges ~10x what `nmae` pays per 1% of amplitude. **Amplitude is closed both ways.**

Also recorded: the `nmae` bracketing that motivated the study holds on a12
(0.9189 / **0.9448** / 0.9630) and **fails on a13** (board 1.0189, outside both harnesses by
0.065). One arm of two. It must not be leaned on.

### 5. Closed by measurement this session

Per-decile magnitude calibration (loses to a transferred flat scalar in 3 of 4 cells — it is the
1.02x cell-line shape, not the 1.17x model shape); per-perturbation reallocation (exactly
`pds`-free, oracle ceiling only +0.004..+0.009, buildable +0.0004); a constant directional bias
(**structurally inexpressible** — gated genes are 99.4% of profile mass and renormalisation
annihilates it, measured ±0.0000); expression tilt; de-offset.

The whitener re-opened on `nmae` and closed again for a better reason: 097 had recorded "`nmae`
reads as zero, disagreeing in sign", which is **retracted** — `nmae` improves at every depth on
both surfaces. But `pds` costs **3.4-8.7x** the `nmae` gain, so `gamma` stays at 0.65.

### 6. Where this leaves the campaign

Calibration is spent. The remaining `+0.0432` is **83% one problem**: per-gene direction. The
donor-agreement key is the only confidence signal that has ever transferred across our surfaces,
and it captures about a fifth of a `1.97x` per-perturbation oracle ceiling. **Winning needs a
better predictor of which genes move and which way — not another metric lever.**

---

## 2026-08-30 · Experiment 103 — residual-only donor-axis correction is a near-positive NO-GO

**Expected:** experiment 102's residual-only diagnostic suggested that subtracting each real
donor's projection on its own non-panel donor-mean response axis could improve the exact staged
a15 direction. The production-shaped KOLF confirmation was frozen with exact a15 reconstruction,
20 source-label nulls, source-derived clusters, and a conservative promotion gate before KOLF
was opened.

**Observed:** exact staged-a15 baseline PDS `0.620408`; residual-only candidate `0.632124`;
delta **`+0.011716`** over 273 targets and 17,806 score genes. The target-cluster bootstrap
95% interval was **`[-0.003138, +0.027365]`**. The largest of 20 frozen nulls was
**`+0.011757`**, just `0.0000404` above the observed gain; add-one one-sided p = **`0.0952`**.
All provenance, reconstruction, score-axis, PM-support, and source-coverage checks passed.

**Consequence:** **NO-GO.** The effect's sign and size are scientifically interesting, but the
preregistered uncertainty and null gates both fail. Do not build counts, submit it, or tune the
projection using this KOLF result. A future confirmation requires a different held biological
context and a newly frozen protocol.

Operationally, experiment 085's required 1 GiB synthetic preflight passed (24 checks; peak
0.255 GiB), but its dependency checksum found drift in three shared production files. The
expected hashes match commit `2d5cd1d` exactly. That fail-closed result is correct: use an
isolated checkout of the pinned commit for the one-time launch rather than blessing current
HEAD by changing hashes.

---

## 2026-08-30 · a14 board result and experiment 085 source-gate result

### a14 — the per-context amplitude correction worked

**Expected:** raw `mse 0.9013 -> 0.8813`, `score_mse +0.0210`, total
`+0.0036`, with small positive pds/nmae and flat DE members.

**Observed:** raw `mse 0.882775`, `score_mse +0.019452`, and overall
**`0.1704242443`**, `+0.00243683` over a12. `score_pds` rose `+0.000729`;
`score_nmae` fell `-0.002002`; fid/reach/jac were slightly negative.

**Consequence:** the per-context cliff model is validated on its main member
and the amplitude correction stays. One registered falsifier fired: the
simple monotone-amplitude model for nmae is wrong. a15 is not another
amplitude arm, but its own directional falsifiers must govern the second slot.

### Experiment 085 — PerturbMatch displacement cannot route PerturbMatch

**Expected:** if the arm-C/arm-D angular displacement contained transferable
information, one monotone stump would improve both HCT116 and HEK293T by at
least `+0.0015` macro and beat the full-selection permutation null.

**Observed:** exact experiment-064 reconstruction passed on both donors, but
the selector chose `always_D`. The best nonidentity rule (`high:q=0.20`) was
`+0.000432` on HCT116 and `-0.001255` on HEK293T, macro `-0.000411`.
The selected identity rule had permutation p `1.0`; null p95 was `0.001810`.

**Consequence:** **NO-GO before KOLF.** The runner did not create the global
access marker and did not authenticate or open a KOLF file. Close this
single-feature routing family. The run used the exact pinned source from
commit `2d5cd1d`; current HEAD was deliberately rejected because three shared
dependency hashes had changed.

---

## 2026-08-30 · experiments 104-107 and STATE MPS preflight

- **Experiment 104: NO-GO.** The direct-credit adapter found positive
  continuous-credit and PDS signal, but H1 top-20% sign precision fell
  `0.0414` and the frozen permutation gate failed on H1 and HCT116.
- **Experiment 105: extraction repaired.** Target eligibility now counts only
  guides that pass the five-cell guide floor. The sealed TeloHAEC extraction
  retains 2,167 targets, covers 65/300 panel perturbations, zeros all 300 panel
  coordinates, and passes five tests. The model stage is paused until the
  strict outer-context and target-coordinate masking protocol is implemented.
- **Experiment 106: active priority.** The official STATE repository is pinned
  at `9bbfe78a434a55205e4de834e1ea99f85f7a3add`; the notebook and dependency
  lock are hashed; 37 upstream tests pass. An opt-in MPS placement patch leaves
  CUDA/CPU behavior unchanged. Both `state_sm` (162.6M parameters) and `state`
  (213.6M parameters) passed full-dimensional MPS forward/loss/backward/update
  preflights without CPU fallback. The official support archive is being
  acquired and will be CRC32C/SHA-256 verified before extraction.
- **Why STATE moved to the front:** the 447-team board snapshot has `GZ YAI`
  rank 2 at `0.1860290893` under the public name `STATE baseline`; a14 is rank
  10 at `0.1704242443`. STATE's member profile is complementary, though the
  componentwise endpoint envelope is only `0.193235`, so first place requires
  blend synergy or another directional source.
- **Experiment 107: preregistered.** GSE249595 is a high-MOI Jurkat CRISPRi
  screen whose released guide library overlaps 299/300 VCC targets. Because
  each cell carries several guides and only 375 targeted-primer genes are
  high-sensitivity, the experiment forbids target-carrier pseudobulks and
  preregisters joint additive guide deconvolution plus held-guide and
  outer-context gates. Metadata files are locally hashed; the raw archive has
  not yet been downloaded.

No challenge submission was made. a15 remains built and unsubmitted.

---

## 2026-08-30 · Experiment 082 valid run — ESM2 residual rejected before held truth

Run `20260830T042050Z` completed successfully on an AWS `r7i.8xlarge`; all
output artifacts were recovered and hash-verified, and the instance was
stopped. Staging peaked at 4.44 GiB RSS and evaluation at 0.52 GiB.

All 60 nonidentity ESM2-kernel residual candidates reduced PDS in every
outer-donor × inner-teacher cell. The best was rank 32 / ridge 1 / beta 0.02
at macro delta `-0.005120` with worst-cell delta `-0.005718`. The frozen
selector therefore chose identity and failed the non-panel pre-gate. Held
HCT116, HEK293T, and KOLF panel truth was never opened. No counts or challenge
submission were produced.

Conclusion: direct linear transfer from ESM2 target similarity into the
arm-D residual is a clean NO-GO under the leakage-safe protocol.

## 2026-08-30 · Experiment 106 held-screen blockers repaired

An independent audit found two execution blockers and three provenance gaps
before the fixed-step STATE run completed. The control builder passed a tuple
to PyArrow's `columns=` argument; checkpoint tree hashing included bytecode
files that the snapshot intentionally omitted; the control count, held
operator, and common target axis were not all revalidated end-to-end.

The runner now uses a PyArrow list, hashes the same source-file set it copies,
requires exactly 18,400 controls, binds each held operator plus metadata
through inference and scoring, and pins the exact sorted 132-target HCT116 and
HEK293T axes and their 120-target intersection. Scoring independently rejects
any allowlist count or hash drift. **51 lightweight tests**, Ruff, Python
compilation, shell syntax, and an exact real-axis hash calculation pass.

---

## 2026-08-30 · STATE closed; Jurkat protocol repaired before expression access

### Experiments 106 and 113 — fixed `state_sm` checkpoint is a NO-GO

The pinned step-20,000 checkpoint completed and froze successfully. On the
120-target HCT116/HEK293T held screen, pure STATE PDS was
`0.50252 / 0.50028` against incumbent `0.80917 / 0.79454`; the smallest
10% direction blend lost `-0.02574` macro PDS.

Experiment 113 tested whether this was merely out-of-context transfer failure.
It was not. On H1, the model's own training-domain context, four deterministic
control banks produced PDS `0.4940–0.5133`; pooled PDS was `0.50448` against
incumbent `0.86345`. The checkpoint raises mean cosine without retaining
perturbation identity. No challenge build or submission was made.

### Experiment 110 — arbitration stopped before two expensive retrains

The perfect-information HEK293T disagreement-routing ceiling was `+0.00686`,
below the preregistered `+0.015` minimum. The first scaffold also applies only
on a truth-derived significant-gene mask unavailable at deployment. Both
conditions independently force NO-GO.

### Experiment 112 — control-label correction is an algebraic no-op

The legacy Orion predicate incorrectly includes `NEGR1` as a control, but this
adds one common vector to every donor effect row. Donor decomposition and
whitener fitting subtract the row mean, annihilating the translation exactly.
No expression or held truth was opened.

### Experiment 107 — acquisition resumed and estimator made identifiable

The manifest-bound 3.75 GB GSE249595 archive download began and safely resumed
after one interruption using the frozen remote validator and partial-file
sidecar. Before any transcriptome coordinate value was read, the unidentifiable
guide-held-out ridge proposal was replaced by a target-level, untreated-only,
channel-held-out ridge protocol in `preregistered_protocol.json`.

The source paper reports more than 586,000 retained cells, a median of 13
guides per cell, deMULTIplex hash assignment, and binomial guide calling at
BH-adjusted `p <= 0.05`. The frozen protocol follows those methods, reserves
activated cells for diagnostics only, validates on six disjoint Nadig Jurkat
galleries, and permits only masked post-whitening residual integration with
bit-exact fallback. Ordinary fourth-donor insertion and unvalidated
context-A-only routing are forbidden.

At 07:46 UTC the live board had 459 entries. a14 remained the latest and best
submission at `0.1704242443`, rank 11 overall and rank 1 on scaled PDS. a15
remains built and unsubmitted.

---

## 2026-08-30 · Pan residuals and A549 sources

### Experiment 117 — Pan H1 direct residual

The repaired one-shot evaluator exposed only the frozen primary arm on the
seven authorizing surfaces and the separate all-target diagnostic on H1.
Every source, mapping, axis, truth file, and output was hash-bound.

The candidate was positive on all seven PDS surfaces and exceeded every frozen
target-label null:

- HEK293T `+0.00005783`;
- six-gallery Jurkat macro `+0.00233185`;
- equal-context macro `+0.00119484`;
- bootstrap lower 95% `+0.00083821`;
- null maximum `+0.00043153`.

It still failed the preregistered `+0.003` macro gate, and HEK293T directional
credit was exactly zero. **NO-GO.** The H1 same-context diagnostic was much
larger (`+0.01947` PDS), reinforcing that the main limitation is
cross-context transfer rather than source extraction.

### Experiment 118 — PDS-null cross-panel signal

The selected CD4 residual changes only off-diagonal coordinates among the 300
VCC target genes. Every non-panel coordinate and target-self coordinate is
zero, and every row has exact zero panel sum, so the construction is
mathematically invisible to PDS when applied in mass-conserving linear-mean
space.

On the one-shot 64-target Pan confirmation, broad signed credit transferred:
`+0.02376578`, bootstrap 95% `[+0.01145745, +0.03719858]`, above both frozen
credit nulls. The top-ranked coordinates did not transfer: weighted sign
precision was `0.35479`, below the target-label null p95 `0.40375`. Direct
unit-scale application also violated nonnegativity, while every other
PDS-invariance and mass-conservation check passed.

**NO-GO for the frozen candidate.** The family remains scientifically
interesting only through a new source-only multi-context estimator with an
independently validated confidence ordering and positivity rule.

### Experiments 119–121 — A549

Two independent A549 sources were acquired and compacted:

- Altair: 20 measured VCC targets and 7,843 exact canonical response genes;
- Myllia: 14 QC-passing VCC targets and 16,477 transferable non-panel genes.

Their 20-target, self-coordinate-free, panel-free cross-source comparison was
a strong negative: macro PDS `0.45263`, null p95 `0.575`, bootstrap interval
`[0.33947, 0.56974]`, and negative directional credit in both directions.
Ordinary direct A549 donor insertion is closed.

At 16:30 UTC the board had 474 entries. a14 remained rank 11 at
`0.1704242443`; the leader was `0.2112235154`. No submission was made.

---

## 2026-08-30 · Experiment 107 terminal result and source-model queue

At 16:55 UTC the live board had 476 entries. `a14` remained rank 11 at
`0.1704242443`; Aginglab.com remained first at `0.2112235154`. Our scaled PDS
remained above the leader's but was now rank 3 at `0.738614`; `pop` had taken
the member lead at `0.755882`. `a15` remained built and unsubmitted.

### Experiment 107 — high-MOI Jurkat deconvolution: NO-GO

The immutable v2 evaluation completed under
`experiments/107-jurkat-high-moi/work/strict-evaluation-results-v2/`.
The six-gallery nested PDS gain was `+0.0042456336`, with target-cluster
bootstrap lower 95% `+0.0012727146` and permutation maximum
`+0.0011594203`. This establishes a small non-random source signal.

The frozen promotion decision is nevertheless **NO-GO**. One of six held
galleries was negative (`-0.0003567447`), H1 directional credit changed by
`-0.0001648071`, and HEK293T directional credit changed by
`-0.0001034474`; `promotion.passes=false`. The result manifest SHA-256 is
`8f5b8e748b2613239d3d7f0e9e295441e10f8c6e962865b88286312fbbc2d595`.
No fraction, mask, or target subset may be tuned against this outcome.

### Experiment 122 — panel-response atlas: source-only gate passed

The frozen source-only phase selected `equal_mean_agreement`. Signed credit
was positive on K562 (`+0.0220050`), Orion (`+0.0529608`), and CD4
(`+0.0149301`); both frozen null p-values were `0.004975`. All 300 targets
were supported across 89,541 nonzero panel pairs, all exact PDS-null
invariants passed, and 11/11 tests passed. This earns one independent
confirmation but does not itself authorize production use.

DLD1 is not currently a safe confirmation surface. An independent audit found
bypassable two-transfer provenance, insufficient receipt binding, non-atomic
publication, mutable-manifest self-authentication, and a non-exclusive
matrix-access lock in experiment 116. Its frozen reliable-row count is also
zero, so no direct-donor held run is authorized. The active scout is comparing
HCASMC and primary-macrophage sources instead.

Experiment 080 remained active on its isolated AWS worker; Experiment 111's
frozen source-only fit and Experiment 123's exact SLIM reproduction remained
active. No challenge submission, count build, commit, or push was made.

---

## 2026-08-30 · Experiments 116, 123, and 124

### Experiment 116 — DLD1 source valid; broad held route NO-GO

Independent verification compared all 33,372,300 finite canonical author
values bit-for-bit with zero mismatches. The compact artifact and split
retrieval are real: Low1-to-Low2 PDS `0.94346450`, Low2-to-Low1
`0.95070104`, macro `0.94708277`, null p95 `0.56037992`.

The frozen held-context gate still fails. Median split cosine is only
`0.06963`, the reliable-target fraction is zero, and the 2,536 measured
response coordinates cap broad-axis coverage at 13.91%. No held context is
authorized under experiment 116. Any future panel-only DLD1 use requires a
new sparse-response preregistration.

### Experiment 123 — exact SLIM/STRING reproduction: NO-GO

The exact source-only reproduction found genuine but insufficient STRING
signal. Macro PDS gain was `+0.00119048` against the frozen `+0.003` gate:
K562 `+0.00182073`, Jurkat `+0.00287115`, RPE1 `-0.00112045`.
The PDS bootstrap 95% interval was `[-0.00121310, +0.00364374]`.
Signed-credit gain was positive (`+0.00292760`), but top-20 precision gain was
only `+0.00165254` and its bootstrap crossed zero. H1 and Orion held outcomes
were not opened.

### Experiment 124 — panel-simplex integration: retracted NO-GO

An independent audit found that Experiment 122's direct effect-space
PDS-null claim does not survive production exponentiation and whole-profile
normalization. All 300 source residual rows have nonzero sums, and naive
addition changed at least one non-panel abundance in every one of the 900
A/B/C context/target profiles.

A separately frozen within-panel simplex transport appeared to repair the
mechanism at profile level, but an independent audit found two decisive
defects. First, quasi emission floors every gene total independently:
redistributing panel mass changed the sum of panel floors in 282/300 A rows,
289/300 B rows, and 278/300 C rows at dose 1.0, so the scorer's all-gene
normalization still moved non-panel log profiles by roughly `2.1e-6` to
`3.0e-6`. Second, the transform used `q*exp(r)` even though the source
residual is in `log1p(50000*p)` units. The preserved recovery files also fail
to reproduce the historical freeze hashes recorded in their receipts.

Experiment 124 is therefore **NO-GO**, not a mechanical pass. Experiment 128
is the clean replacement: use the exact source-unit transform and preserve the
mutable panel's summed emitted integer budget, all non-panel totals,
target-self, and the grand total exactly. Production remains unauthorized.

### Experiment 128 — integer-budget transport: mechanical PASS

The clean replacement passed all required invariants across contexts A/B/C,
all 300 targets, four diagnostic doses, three fixed median-library values, and
the exact replayed production library pick for every context/target: 14,400
complete evaluations.

The transport applies Experiment 122 in its frozen
`log1p(50000*relative_abundance)` units, clips the rare negative preliminary
means, and redistributes only the already-existing mutable-panel integer
budget by deterministic largest remainder. Every non-panel emitted integer
total, target-self total, full panel total, grand total, and
scorer-normalized non-panel linear/log/effect value remained bit-identical.
At dose 1.0, median post-integerization sign retention was 0.959 / 0.938 /
0.952 for A/B/C.

Freeze SHA-256:
`67df1b5c5e4f1960fc4ec7a1cc06d7f0754b485443b4e1b7478daefee3ac1426`.
Result SHA-256:
`660967723f49fde57815294b60f7f3e9820932d5ac269b78ab5888dc4b750019`.
This is a mechanical authorization only; production use remains a NO-GO until
the exact transformed candidate passes independent biological confirmation.

The source scout found 27 exact VCC targets in IGVF HCASMC and 12 in a primary
macrophage Perturb-seq source. Neither reaches the proposed 40-target floor;
HCASMC is the strongest orthogonal secondary confirmation. Experiment 125,
a nested hierarchical random-effects panel posterior, was started as the next
source-only model improvement.

No challenge truth, `cell-eval2`, count build, submission, commit, or push was
used.

### Experiment 111 — DepMap context program: NO-GO

The frozen source-only fit selected `rr4_cr4_l1_b1_cap5`. Its tune delta was
`+0.01034858`, but the untouched nested-test delta was only `+0.00490196`,
below the preregistered `+0.010` promotion gate. Five of six batches were
nonnegative and the worst batch stayed above `-0.003`, but the target
bootstrap lower 95% was `-0.00095316`, the geometry delta missed `+0.005`,
and the observed result did not exceed the context-label null p95
`+0.00680828` (`p=0.1683`). The gain was concentrated in `apollo31`
(`+0.01960784`).

The result and artifact-manifest SHA-256 values are
`3b5855b4a69042af4958711f4ed4475cea9009550fb3a2d2259a4aa7ff44aaaa`
and `a02bd4967ac5eda9ec7cd171d3c332482fdfe973b7375ae0701c303da2056648`.
All 37 tests and lint checks passed. Held-context access is not authorized.

---

## 2026-08-30 · Experiments 129–132 and the current frontier

### Experiment 129 — Pan response-neighborhood transfer: NO-GO

The source-only learner used Pan H1 only to retrieve biologically similar
training targets. It never added the query target's own Pan response, and
RPE1/Jurkat teacher outcomes were kept strictly opposite-context.

The identity fallback won the frozen 100-target selection. Nonidentity
candidates improved signed directional credit in some cells but reduced
complete-gallery raw PDS by `-0.0110` to `-0.0209` macro. The untouched
55-target source-test split was not opened. No held HCT116, HEK293T, KOLF,
Pan-panel, DLD1, HCASMC, challenge outcome, counts, or submission was used.

### Experiment 130 — TeloHAEC confirmation closed without value access

Experiment 130 required an immutable Experiment 125 `GO`. Experiment 125 is a
terminal NO-GO, so the candidate-binding amendment remains absent and the
TeloHAEC panel-value gate can never open under this experiment number. The
metadata-only freeze remains preserved; no TeloHAEC expression value was
read.

### Experiment 131 — differentiated-cell metadata audit

HUES8 ESC, HUES8 definitive endoderm, and WTC11 hepatocyte each retain
91/300 VCC targets after QC and expose broad response axes of 9,519, 9,584,
and 11,138 genes. WTC11 cardiomyocyte retains 91 targets in the canonical QC
summary, 90 at the 100-cell floor, with a 7,260-gene response axis.

The broad canonical H5MU and calibrated gene-level outputs are not anonymously
downloadable. Synapse requires an authenticated certified account with
READ/DOWNLOAD permission. A public cardiomyocyte log2FC table exists, but it
belongs to the sparse portal family already tested by Experiment 081 and is
not sufficient evidence for another acquisition or experiment.

### Experiment 132 — DepMap response-decomposition Ridge: NO-GO

The exact upstream PCA-50/Ridge-10 target-feature model was evaluated as a
replacement for the same-target source mean on 613 frozen selection targets.
Its macro raw-PDS delta was **`-0.204710`** across K562-GWPS,
K562-essential, RPE1, and Nadig Jurkat. Macro signed-credit delta was
`-0.002369`; three of four contexts were negative. Real DepMap features beat
permuted features, so they carry information, but the learned replacement
throws away much stronger observed same-target cross-context signal.

Identity was selected and the untouched 410-target source test stayed sealed.
Future context-feature work must be residual-only and preserve the incumbent.

### Current frontier

Experiment 080 remains active on its isolated AWS worker. At 18:06 UTC its
service was healthy after 5.6 hours, the worker had 244 GiB available memory
and 195 GB free disk, staging had passed all 24 tests, and held truth remained
closed. No S3 terminal artifact existed yet.

The proposed margin-cone transform is mathematically deployable but not yet a
frozen experiment: it still lacks an independently authorized directional
proposal and a fresh source-target split. A basal-conditioned residual model
is also not yet runnable because the fresh WTC11 outcome lacks a valid basal
control vector and the remaining fresh cross-context target budget is too
small.

No challenge submission, `cell-eval2`, count build, commit, or push was made.

---

## 2026-08-30 · Experiments 114, 126, 127, and the H1 transfer queue

### Experiment 114 — DepMap/STRING neighbor replacement: terminal NO-GO

The frozen 41-row weak-support replacement family was evaluated only on its
authorized 35-target HEK293T raw-PDS diagnostic. The exact incumbent scored
`0.70420168`; STRING-only scored `0.56134454`, DepMap-only `0.54873950`, and
the primary weight-1 blend `0.58403361`.

The primary delta was **`-0.12016807`**. Its 100-draw target-label null p95 was
`-0.11932773`, so it did not even exceed the negative null threshold. The
result closes leave-target-out co-essentiality replacement of a14's
weak-support rows. H1 was not opened and should not be opened for this family.

### Experiment 126 — HCASMC acquired but still opaque

The exact 1,531,381,676-byte author SCEPTRE gzip was downloaded and published
read-only. Its compressed MD5 exactly matches the published
`4d3067b3e934842c88b9ea52d74780f9`; local SHA-256 is
`b87452c2db1f508f5013bc15dff8088828ab931b35492d77ed6aea569fb5d9b9`.
Metadata establishes 27 exact VCC targets. The gzip was never decompressed,
tested, opened, or parsed, and its value-access boundary remains closed.

### Experiment 127 — panel-only ceiling

Across H1 and the two repaired Orion development surfaces, mutable panel
coordinates contain only **2.283%** of reference DE genes and **2.284%** of
absolute reference LFC mass. Perfect panel correction improves incumbent raw
NMAE by at most **0.024665** on the equal-surface oracle average. On the one
surface with an exact prediction DE table, the separate metric-specific raw
ceilings are `+0.030948` fidelity, `+0.025908` reach, and only `+0.000579`
Jaccard. These optima are incompatible and cannot be summed.

The panel response family may still supply a small orthogonal gain, but about
97.7% of the relevant directional signal lies outside the panel. It cannot
close the leader gap alone.

### Experiment 129 — broad H1 source-only learner, with a corrected boundary

Experiment 115's verified GSE295214 compact donor supplies 2,354 H1 CRISPRi
effects on 16,290 covered canonical genes and 64 VCC targets, including 57
rows absent from the old H1 donor. It overlaps 1,348 RPE1, 1,389 Nadig-Jurkat,
2,216 K562-GWPS, 2,207 HCT116, and 2,225 HEK293T target labels, so a
high-powered non-panel source-only learner is possible.

This is not a fresh direct-transfer source. The same compact donor powered
Experiment 117's direct residual and supplied Experiment 118's 64-target Pan
confirmation. Experiment 129 was therefore redirected before any result:
direct residual integration, target masks, fractions, and confidence rankings
from those experiments are closed and may not be retuned. The experiment may
continue only with a materially distinct model selected and tested on
source-only non-panel splits.

The independent-source census found no public dataset satisfying all of:
at least 40 VCC targets, broad response coverage, and no consumed outcome.
IGVF TeloHAEC is the best actionable panel-only confirmation at 65 targets
and 17,464 response symbols. Its broad direct-donor route failed Experiment
105, but its panel coordinates remain unopened. HCASMC remains a second,
more orthogonal but smaller 27-target confirmation.

No `cell-eval2`, count build, challenge submission, commit, or push was made.

### Experiment 125 — hierarchical panel posterior: NO-GO

The immutable source-only result verified all 16 files. The selected
hierarchical posterior increased top-20% truth-weighted sign precision by
`+0.01217457`, with bootstrap 95% `[+0.00640310,+0.01791314]`, but it failed
both complete-procedure precision null thresholds (`+0.01351644` response
label and `+0.01552727` target label).

The decisive failure was broad signed credit: macro delta
**`-0.01504990`**, with all three held-family deltas negative
(`-0.01214611` K562, `-0.02608242` Orion, `-0.00692117` CD4) and bootstrap
95% `[-0.02007071,-0.01014371]`. All 1,800 transport/PDS invariant checks and
100% support retention passed, so this is a biological/statistical rejection,
not a transport failure.

Artifact manifest SHA-256:
`b3efbd1db8a39ab933c4e4fe0bd998917216866df6157441a33a0e1560cdc63a`.
Result SHA-256:
`635701eefa74271e9209489b7593c045a22b573576c5a6e4319c7e7f99621010`.
No independent confirmation is authorized for this candidate.

---

## 2026-08-30 · Experiment 135 TeloHAEC atlas confirmation

Experiment 135 is a terminal **NO-GO**. The independently reviewed evaluator
opened the previously unaggregated TeloHAEC panel responses once for the exact
Experiment 122 `equal_mean_agreement` atlas, using the staged-a15 comparator,
pooled controls, one all-300 `default_rng(0)` library replay, and the exact
Experiment 128 integer-budget transport.

Signed credit improved by `+0.004066`, but missed the frozen `+0.005` gate and
had bootstrap lower 95% `-0.000456`. Top-20% magnitude-ranked sign precision
worsened by `-0.005245`, with bootstrap lower 95% `-0.020571`. Median
transport sign retention was `0.890`, below the required `0.900`.
`promotion_passes=false`.

All integer and scorer-space invariants passed, all 10 output artifacts
verified byte-for-byte, and the output is read-only. No truth matrix,
per-target outcome, null-score sample, or bootstrap sample was persisted.

- implementation freeze SHA-256:
  `0f6a218d590853258f1456163089c60581b773fa2ab0a4064eead4ad1426a743`
- result SHA-256:
  `b6fc074904eaf11b6f954a9ba21d9425e1b67c096338d4c8ca9b8412dc0e99b5`
- artifact-manifest SHA-256:
  `0aaf2a96bfd6345c1fb5cae61fb607e8767b6ce012221cd561d28c9c05d33d81`

Do not rerun, rescale, subset, tune, integrate, or submit this atlas. Do not
build an a15-plus-panel candidate from it.

Experiment 080 remains the only active heavy evaluation. At 2026-08-30
18:52 UTC, its AWS systemd service was still active, PID 8384 was using about
798% CPU, the service cgroup used 10.9 GiB, 244 GiB memory remained available,
and 195 GB disk remained free. No terminal S3 artifact existed yet.

No challenge submission, count build, `cell-eval2`, commit, or push occurred.

### Experiment 136 — DepMap residual margin screen

Experiment 132's rejected Ridge prediction was tested as a residual around the
same-target source mean on the already-consumed 613-target selection split.
Every nonzero global dose lost PDS in every scored context. The smallest dose,
`0.0025`, changed macro PDS by `-0.0002590` while changing signed credit by
`+0.0001216`; dose `0.20` gained `+0.010110` signed credit but lost
`-0.021078` PDS.

A source-rank-preserving per-target rule lost `-0.020302` macro PDS. Exact
nearest-decoy-margin preservation reduced to identity. Even an invalid
held-outcome oracle reached only `+0.001926` macro PDS and remained negative
on K562. The untouched 410-target split remains sealed and is not authorized
for this proposal.

Result SHA-256:
`bbb26d91e4e75d72c88b0480a3876c0f594eefb0882dd006c77311e616040973`.
Runner SHA-256:
`b3885ae1dfeaa070f7d7ade59a1e24f969893d923d78d70f0b18a788cb189bf7`.

### Experiment 080 — listwise adapter terminal result

The frozen AWS evaluation completed successfully and all artifacts were
recovered before the worker was stopped.

The source-only pre-gate was genuinely positive. It selected rank 8,
`beta=0.02`, temperature `0.03`, regularization `1e-4`, with non-panel macro
PDS gain `+0.0233758` versus a complete-selection null p95 `+0.0212175`.
Every source-only pre-gate passed.

The held-panel result was a terminal NO-GO:

- HCT116 raw-PDS delta: `-0.000688396`;
- HEK293T raw-PDS delta: `+0.000643501`;
- primary macro delta: `-0.000022448`;
- KOLF confirmation delta: `+0.000686813`;
- cluster-bootstrap 95%:
  `[-0.001863159,+0.001803298]`;
- permutation-null p95: `+0.000119721`;
- `promotion_gate.passed=false`.

The adapter preserved target coverage and stayed below the 10-degree angular
limit, but failed the nonnegative HCT116 gate, `+0.003` macro gate, positive
bootstrap lower bound, permutation threshold, and cluster-dependence gate.

Recovery verified 331 files total, all 328 manifest-governed artifacts, all
303 model files, nine query files, three cluster files, the selection freeze,
the permutation plan, and every pinned implementation hash. The EC2 instance
`i-08300b395f88f6f55` is stopped with state code 80.

- result SHA-256:
  `9b12cdcb440081bd3910c72d00ba3e41f1a137b2024119c1c1a258c6ddb8ccb2`
- selection-freeze SHA-256:
  `8a7335cf21ff6c576c5ce47eaf3a85d14d6e212a62e3680c8cddcf57204f8524`
- artifact-list SHA-256:
  `971602b00c44dc2bf7747ff34c7ff94131f50c3a947417aa09fc0bc40d228312`

Do not integrate, retune, or run count-level confirmation for Experiment 080.

### Experiment 137 — high-MOI Jurkat tangent cone: terminal NO-GO

The exact a14 incumbent and frozen Experiment 107 high-MOI Jurkat source were
combined only through preregistered shortest-arc moves capped at
`0°, 1.5°, 3°, 6°`. The source-only operator and 394-target candidate were
built without opening Nadig truth. The first development process encountered
a non-contiguous byte-view assertion after opening only the authorized 315
development rows, before emitting a score or writing a result artifact. A
numerically blind `ascontiguousarray` repair was documented, bound, and passed
two independent read-only audits before the same locked run resumed.

Development selected `6°`, but nested macro raw PDS changed by
**`-0.000921659`**. Fold deltas were
`+0.000256, -0.000256, -0.005632, +0.001024, 0.000000`; the 20,000-draw
cluster-bootstrap 95% was `[-0.00430108,+0.00593958]`, and the observed result
did not beat the 99-permutation maximum `+0.000256016`. Signed credit was
slightly positive overall (`+0.000429369`) but negative in one fold.

All exact-identity, exact-fallback, and angle invariants passed.
`promotion_passes=false`. The untouched 79-target test remains sealed; no
candidate freeze, count build, scorer run, submission, commit, or push was
made.

- development summary SHA-256:
  `8e7a898367c2a9930a5c72ab361cb4ba86823853c7e14d6ed454c52ee23515d2`
- development arrays SHA-256:
  `bbfb9dda83b190e14933ba4448d8a8e2f0aa72f75d0e53c601d537f840f011b3`
- runtime amendment SHA-256:
  `2551a15548cabee6a5156494d851e140ee22603531a082f58eff95382d5615df`

### a15 — CD4 residual transfer confirmed on the board

Entry `KpoZt5Rie0EjgMaZQECE` published under the opaque model name
`The Fern Filed An Appeal`.

- overall: **`0.1740865087`**, rank 9 of 482 at the 2026-08-30 20:54 UTC
  read;
- delta versus a14: **`+0.0036622644`**;
- scaled PDS `+0.007256`, MSE `+0.009441`, NMAE `+0.002266`,
  fidelity `+0.002190`, and reach `+0.000826`;
- scaled Jaccard was effectively flat at `-0.0000066`.

All four preregistered falsifiers passed. The submitted package
`submissions/a15-cd4f010.vcc` has SHA-256
`55f7bd41bfc4d751bcc125d42e25b40c7aed474589257c6db4b575afe26aac59`.
This is the current live best entry and must be preserved.

### Experiment 138 — DLD1 replicated sparse-sign gate: terminal NO-GO

The frozen development procedure selected `k=200, beta=0.05`. Nested macro
raw-PDS improved by only **`+0.0000995142`** and signed credit by only
**`+0.0000602495`**. Both changes exceeded their frozen null thresholds, but
were roughly 30 and 50 times below the `+0.003` practical promotion gate.
RPE1 PDS was negative in two folds and PDS positive-gain concentration
slightly exceeded its frozen maximum.

`promotion=false`. The untouched 410-target test remains sealed. Do not
rerun, retune, authorize test access, or integrate this family.

### Experiment 139 — high-MOI Jurkat donor-simplex routing: terminal NO-GO

The frozen router used 138 high-MOI Jurkat landmark coordinates to choose
convex weights among five allowed same-target source contexts. It selected
exact identity in every outer development fold and on the complete
613-target development set, so nested PDS and signed-credit deltas were
exactly zero.

`promotion=false`. The untouched 410-target test remains sealed. Do not
rerun, retune, authorize test access, or integrate this family.

- development result SHA-256:
  `c72ca68b3fd52792971e872959faa02f97c23c9f846235fae0f8f1a44c595049`;
- development statistics SHA-256:
  `e9dee51616a2c1dd10a5de8e5a0288daea6a3de7480318f1e36b6bf1fb034325`;
- development artifact-manifest SHA-256:
  `73025af13f32f3991b433c125515e7912e1fe9ae7bf37478ba07ef00cc6a079c`.

---

## 2026-08-30 22:25 UTC — Experiment 140 batch-matched controls: terminal NO-GO

The final v5 protocol corrected production precision, family-level eta
eligibility, nested/fixed inferential gate parity, source-hash enforcement,
and preparation provenance before any development outcome access. All 16
tests and Ruff passed. The source freeze rehashed 373 inputs; preparation
peaked at 9.77 GiB and identity at 4.55 GiB.

All three production donors supported all 613 development targets. Frozen
global controls reconstructed exactly. Exact candidate controls excluded 154
HCT116 and 233 HEK293T NEGR1 cells that the legacy `^neg` rule had treated as
controls. Eta-zero identity passed at raw, centered, pre-whitening, and
post-whitening stages with the expected 21,610-row K562-held and 29,394-row
non-K562 `union_pruned` bases.

The result was decisively negative:

- fold eta selections: `0, 0, 0.25, 0, 0`;
- final eta: `0`;
- nested macro PDS delta: `-0.0000879634`;
- nested macro signed-credit delta: `-0.0000017071`;
- nested PDS bootstrap 95%:
  `[-0.000281220, +0.000073152]`;
- nested signed-credit bootstrap 95%:
  `[-0.000016813, +0.000013334]`;
- fixed candidate: exact identity;
- frozen Orion-only diagnostic: PDS `-0.0000852978`, signed credit
  `-0.0000046390`;
- frozen K562-only diagnostic: PDS `-0.0000168819`, signed credit
  `+0.0000029312`;
- checks passed: `5/30`;
- decision: `NO-GO`.

The 410-target confirmation split remains sealed. No production integration,
count build, `cell-eval2`, submission, commit, or push was made.

- result SHA-256:
  `4160c8d319a0894b1784d0c98436c0d05e8c48ce67044a9810048de1f6f7a3c1`;
- development arrays SHA-256:
  `955036d541d9b8a01120c9859c1d61a8631e7c5be6ec3516b0932599f05472d6`;
- artifact-manifest SHA-256:
  `e40654cfe609e3b24aeb495f8efbd1ddcc0f2ea18e41b3a1fad33b932427f6f6`.

---

## 2026-08-30 22:52 UTC — Experiment 141 Pan coordinate calibration: terminal NO-GO

The fixed `diag_rho050` primary learned a nonnegative gene-coordinate
reliability map for the same-target Pan residual under five-fold target
cross-fitting and strict destination-family teacher isolation.

The source-only surrogate result was large in absolute terms:

- primary macro PDS delta: `+0.0119714643`;
- primary macro signed-credit delta: `+0.0076031793`;
- all three biological families were positive on both endpoints.

But the preregistered distinctness comparison failed. The scalar direct-blend
control reached `+0.0118417405` PDS and `+0.0076232622` signed credit, leaving
the diagonal primary only `+0.0001297238` PDS better and slightly worse on
signed credit. The required incremental PDS was `+0.001`.

This is evidence for Pan's already-known direct residual, not for a new
coordinate-calibration mechanism. Experiment 117 has already consumed and
closed scalar Pan interpolation, so the result cannot be promoted or retuned.
The conditional null stage did not run and the 410-target split remains
sealed.

- result SHA-256:
  `d6386c4a1281d0cfb31ea5c7c079b5a81fbf5068c328cf25a3eae07fd1f7b874`;
- arrays SHA-256:
  `35090a88eb2ab1be2c788f9f3421821587e7064a49f0be4273d35700fe08d22b`;
- development-freeze SHA-256:
  `adfc59ec0918f77b37c490e297660b3d81b115035efb9edfd80231903c2ea2ae`.

No `cell-eval2`, count build, challenge submission, commit, or push occurred.

---

## 2026-08-31 18:45 UTC — Experiment 142 batch-disjoint reliability weighting: terminal NO-GO

The final development evaluator repaired six control-plane defects before any
613-target outcome access: frozen-axis hashes are enforced, Experiment 061 is
no longer dynamically executed, failed lock acquisition cannot remove another
process's lock, prepared artifacts are bound to their creating authorization,
truth bytes are verified only after the irreversible access lock is recorded,
and concentration uses the preregistered exact-zero rule.

All 45 combined tests pass. The post-run audit validates both authorizations,
the prepared artifact, the development result, and the recorded access lock
without reopening truth arrays.

The frozen result is exact identity:

- fold eta selections: `0, 0, 0, 0, 0`;
- final eta: `0`;
- nested macro PDS delta: `0.0`;
- nested macro signed-credit delta: `0.0`;
- fixed PDS and signed-credit deltas: `0.0`;
- promotion checks: `2/59`;
- decision: `NO-GO`.

Postmortem fixed arms showed that reliability improved signed credit while
slightly harming retrieval:

- `eta=0.25`: PDS `-0.000020`, signed credit `+0.001183`;
- `eta=0.50`: PDS `-0.000450`, signed credit `+0.001855`.

The non-promotable support controls selected `eta=0.5` and were positive but
below the `+0.003` gates:

| control | fixed PDS | fixed signed credit |
| --- | ---: | ---: |
| target cell count | `+0.00196141` | `+0.00213555` |
| target total UMI | `+0.00199606` | `+0.00208881` |
| occupied-batch count | `+0.00191654` | `+0.00209870` |
| pooled-effect norm | `0.0` | `0.0` |

The positive controls are highly redundant and their largest PDS family gain
is Jurkat. They do not authorize proxy selection or production integration.

The 410-target confirmation split remains sealed. No reserved confirmation
expression, count build, `cell-eval2`, challenge build, submission, commit, or
push occurred.

- evaluator SHA-256:
  `dee877786a5cf83819ce1f954ebc82f63179a36fa31e031fb47c45f2bfcc5c1a`;
- prefreeze authorization SHA-256:
  `204d3bdccf1d11b8c1b9a21942f3888156aad6631fc860cb527e9a3b286f333a`;
- development authorization SHA-256:
  `898b035ff508aeafc8684595af2cd3f0c5d2dbf779e6ff24cdf6aefe44734e3e`;
- development result SHA-256:
  `46715620f4fac830be615a7bf82bfd30da16610e7816a2a94dd29fdb6e5206c4`;
- result manifest SHA-256:
  `a2d730451497b33aab3eeb9181a159fe72f2732cce75e62d6974d3f9ddee44e6`.

---

## 2026-08-31 — Experiment 143 normrel production replication: terminal NO-GO

Experiment 143 tested one preregistered, leakage-safe analogue of Experiment
097's `normrel^-1.5` target-specific donor weighting. The source-only prepare
used only Experiment 142's frozen 613-target by 4,992-gene artifacts; the
numerator and denominator used the same allowed rows and axis, PerturbMatch
remained embedded inside the K562 composite slot, and CD4 was absent because
zero development targets had complete three-condition support.

After a fresh no-truth audit, a stage-specific development authorization was
created with every control-plane and truth-ledger file hash-bound. The guarded
development run wrote the irreversible outcome-access lock before verifying
or loading truth bytes. Peak observed process-tree RSS was `0.17 GiB`.

The result was a terminal NO-GO:

- macro raw-PDS delta: `+0.0027912833`;
- paired-bootstrap PDS 95%:
  `[+0.0011221364, +0.0045882510]`;
- macro signed-credit delta: `-0.0009390308`;
- paired-bootstrap signed-credit 95%:
  `[-0.0016508036, -0.0002102718]`;
- family PDS: K562 `-0.0005291132`, RPE1 `+0.0020258239`,
  Jurkat `+0.0068771391`;
- family signed credit: K562 `-0.0006122562`, RPE1 `-0.0016307352`,
  Jurkat `-0.0005741010`;
- PDS positive-gain family concentration: `0.7724551`;
- checks passed: `3/12`;
- decision: `TERMINAL_NO_GO`.

The small PDS gain was statistically positive but missed the frozen `+0.003`
practical gate, harmed K562, and was too Jurkat-concentrated. Signed direction
was significantly worse and negative in every family. Do not retune, subset,
blend, confirm, integrate, build, score, or submit this candidate.

The 410-target values were not loaded or interpreted; their containing donor
files were streamed bytewise only for preregistered SHA-256 verification
after the access lock. Reserved confirmation expression remains sealed. No
count generation, `cell-eval2`, challenge build, submission, commit, or push
occurred.

- development authorization SHA-256:
  `94fe544a3cb5b00b532a16c5c46a7aabf91a2bfa8467e0979fced2970597d3c3`;
- result SHA-256:
  `7f1bc6fbce304cd7fbd556e385478918d39c3e4e3f95ecbd9e288381ae012e89`;
- result manifest SHA-256:
  `93f4b09333a84b771b5c8c8ecd0981f292b18f5b40f21b2cc7a08431311f2576`;
- outcome-access lock SHA-256:
  `a5ef3e87891607ced8bae138586fcd6957bb36b9d8d90a67c6603f64da3b901e`.

---

## 2026-08-31 — Experiment 145 spectral trace-bin attenuation: terminal NO-GO

**Expected:** the gamma-0.65 whitener might over-amplify its low-variance
spectral tail. Preserving the incumbent row norm while attenuating the bottom
half of cumulative source variance could improve PDS and signed direction
together.

**Observed:** after a failed-audit v1 was preserved and corrected into a
hash-bound v2, one guarded 613-target evaluation found a monotone trade:
positive beta improved signed credit but harmed PDS. At beta 0.25 the fixed
macro changes were `-0.00005109` PDS and `+0.00986897` signed credit; by beta
1.00 they were `-0.00457143` and `+0.03406376`. Nested selection used beta
0.25 in one fold and identity in four, yielding `-0.0000684160` PDS and
`+0.0020382723` signed credit. The final selector chose identity.

**Consequence:** the low-variance tail is not merely amplified nuisance; it
contains target-discriminative retrieval signal. The exact attenuation family
is closed. Its frozen null stage and 410-target confirmation remain sealed.
No counts, scorer, build, submission, commit, or push ran.

The corrected run used selected-row memory-mapped truth access only, peaked at
`0.7933 GiB` process-tree RSS, passed 20 tests and independent audit, and
wrote a terminal `5/17` gate result. Full record:
`experiments/145-spectral-trace-bins/RESULTS.md`.

---

## 2026-08-31 — Experiment 146 context-gated sign repair: source-only NO-GO

The source-only evaluator tested one fixed, energy-capped HGB sign-flip policy
under six leave-one-experimental-batch-out folds crossed with five whole-target
folds. Availability-stratified whole-target permutations preserved every
authenticated missing context-target pattern.

The primary result was positive but scientifically unusable:

- macro signed-credit delta: `+0.0000119081`;
- experimental-batch bootstrap 95%:
  `[+0.0000006788, +0.0000235618]`;
- whole-target bootstrap 95%:
  `[+0.0000000924, +0.0000246646]`;
- macro raw-PDS delta: `-0.0000770895`;
- minimum context raw-PDS delta: `-0.0009259259`;
- positive contexts: `4/16`;
- contexts at or above 65% weighted flip precision: `4/16`;
- permutation-null p95: `+0.0003156613`;
- positive-gain concentration: `63.76%`;
- FDR `0.10` signed-credit delta: `-0.0000083846`.

It failed eight mandatory gates, including context breadth, flip precision,
matched-control superiority, permutation-null separation, FDR robustness, and
gain concentration. An independent read-only audit verified the artifact
hashes, gate calculations, leakage isolation, and terminal `NO_GO` decision.

No challenge control expression, development or held 613/410 outcomes,
`cell-eval2`, counts, build, or submission was opened or created. Full record:
`experiments/146-context-gated-sign-repair/RESULTS.md`.

---

## 2026-08-31 — Experiment 151 relaxation-consensus selector: terminal NO-GO

**Expected:** Experiment 149's dual-moment gain was large on HEK293T
(`+0.012919`) and Experiment 150's uniform transfer lost on HCT116
(`-0.005045`). If the loss were concentrated where HEK and HCT disagree about
the relaxation direction, weakening lambda on the low-consensus half of targets
should recover it. The selector was preregistered and hash-sealed **before**
Experiment 150's result existed.

**Observed:** the gated candidate lost to the incumbent by `-0.004344`
(20,000-draw paired target bootstrap 95% `[-0.009484, +0.000669]`,
`P(delta > 0) = 0.04565`), and lost to Experiment 150's own uniform candidate by
`-0.000393`. Per-target PDS was exactly unchanged everywhere. NMAE `+0.006315`
and Jaccard `+0.014645` improved; fidelity `-0.003736` and reach `-0.043628`
did not. Four of ten gates failed across both comparisons.

**Consequence:** target-level cross-context consensus does not localise where
the dual-moment emitter hurts. Uniform lambda-1 lost `-0.0451` scaled reach;
halving lambda across the entire low-consensus half recovered only `+0.0015` of
it. The reach damage is diffuse across targets, so the whole-transform HEK293T
gain has now failed to transfer to HCT116 in both uniform and target-gated form.

Only the new candidate arm was scored; Experiment 150's frozen incumbent and
uniform frames were reused, so no comparator was recomputed. Peak scorer RSS
`16.781 GiB`, minimum free memory `66%`, elapsed `4508.66 s`, final scorer
process count `0`, no parallel scorers. No challenge build or submission was
created. Full record:
`experiments/151-relaxation-consensus-selector/RESULTS.md`.

This result does not close the gene-wise family — Experiment 152 gates the same
relaxation *within* each target rather than between targets.

---

## 2026-08-31 — Experiment 152 gene-wise consensus ramp: terminal NO-GO, family closed

**Expected:** Experiment 151 showed that target-level cross-context consensus
does not localise where the dual-moment emitter hurts. Gene-wise gating *within*
each target was the remaining independent axis, preregistered and hash-sealed
before Experiment 150's result existed.

**Observed:** the gene-wise candidate lost to the incumbent by `-0.003235`
(20,000-draw paired bootstrap 95% `[-0.008406, +0.001811]`,
`P(delta > 0) = 0.10510`), failing four of six Comparison A gates. Against the
best frozen comparator (Experiment 151) it gained `+0.000716` with
`P(delta > 0) = 0.954` — real, but under the preregistered `+0.001` margin.
Per-target PDS was exactly unchanged. Reach `-0.040483` again carried the loss
while NMAE `+0.007151` and Jaccard `+0.014893` improved.

**Consequence — the family is closed, and the reason is quantitative.** Three
arms now differ only in relaxation dose, with mean lambda over all 132 targets:

| exp | placement | mean lambda | delta |
| --- | --- | ---: | ---: |
| 150 | uniform | `1.0000` | `-0.005045` |
| 151 | target-wise | `0.7727` | `-0.004344` |
| 152 | gene-wise | `0.6170` | `-0.003235` |

Least squares gives `delta = -0.004606 * mean_lambda - 0.000539` with
`R^2 = 0.9449` and residuals an order of magnitude below the bootstrap
intervals. **At zero dose the fit predicts `-0.000539` — doing nothing.** The
transform's harm on HCT116 is linear in how much of it is applied, and no
forward dose reaches the `+0.006` gate. Two independent placement schemes failed
to beat that trend, so cross-context relaxation consensus is not a useful
localiser at either granularity.

Peak scorer RSS `16.615 GiB`, minimum free memory `61%`, elapsed `4593.37 s`,
final scorer process count `0`, no parallel scorers. No challenge build or
submission was created. Full record:
`experiments/152-gene-consensus-ramp/RESULTS.md`.

---

## 2026-08-31 — Post-hoc: the reach overwrite law closes the dual-moment direction

Source-only diagnostic over the already-unblinded 149 and 150 frames, run after
152 closed its family. Full record: `analysis/dual-moment-reach-closure.md`,
reproducible via `analysis/reach_signflip_diagnostic.py`.

The emitter replaces rather than refines the incumbent's directional ranking.
Per-target `reach_delta` regressed on incumbent reach breaks even at `+0.0173`
(HEK293T) and `+0.0163` (HCT116) — within 6% across two unrelated cell lines —
and the candidate's reach is nearly independent of the incumbent's (slope
`+0.037` on HEK). The apparent context-dependence is a difference in input
distribution, not response: HEK's incumbent reach averages `0.0134` with 68/132
targets at zero, HCT's averages `0.1008`.

The fix is unavailable. Incumbent reach does not transfer
(`spearman = -0.033`, conditional zero-rate `0.290` vs base rate `0.292`), and
no source feature predicts it (`|rho| <= 0.16`). Cross-context gating moves the
HCT reach delta only `-0.0377 -> -0.0277`. A greedy per-target oracle on held
outcomes yields `+0.007207` — the family's ceiling sat just above the `+0.006`
gate it had already missed.

No sealed data was opened, no gate changed, no build or submission created.

---

## 2026-08-31 — Experiment 154 confidence-key screen: terminal NO-GO, ordering channel closed

Deterministic screen, no emitter and no seeds. Asked one falsifiable question:
does any available signal isolate a subset of our predictions with >=90% sign
accuracy, per perturbation? Answer: no.

Per-perturbation reach over 127 perturbations, median pool 70: `|pred_lfc|`
`0.10716`, best buildable key (geo-magnitude x donor concordance) `0.12473`,
oracle `0.62797`. The best key captures **3.4%** of the oracle headroom, with a
median pure head of 4 genes against the oracle's 44.

Three durable facts. **Donor concordance alone is anti-useful** (`0.05657`, worse
than `|pred_lfc|`, 42/127 perturbations at `k* = 0`) and only helps multiplied by
magnitude — the magnitude *is* the confidence signal. **"More donors" is not an
available move**: of seven donor matrices, panel coverage is `k562_gwps` 132,
`orion_hek293t` 132, `h1_2025` 12, and `k562_essential`/`nadig_jurkat`/`rpe1`
**zero**; concordance rests on two donors. **Pooled purity is a trap** — pooled,
`|pred_lfc|` reaches 1.0000 top-1% purity, but pooling cherry-picks across
targets while `reach` is scored per perturbation.

Reproduces 095 §6 on an independent implementation, so the instrument agrees with
prior work and only the conclusion is new. Consequence: `reach` is
**model-limited, not order-limited**. Do not build further ordering keys, per-gene
`phi` confidence maps, or spread policies aimed at `reach`. The remaining headroom
is per-gene directional accuracy. Full record:
`experiments/154-confidence-key-screen/RESULTS.md`.

No sealed data opened, no `cell-eval2` run, no build, no submission.

---

## 2026-08-31 — a19: the whitener ridge, one parameter nobody had ever swept

**The lever.** `src/vcc2026/whiten.py:79` computes `w = (evals + lam) ** (-gamma/2)`.
`lam` was threaded through `apply_v9_whitening`, recorded in every meta as
`whiten_lambda`, and left at its default **0.0 in every build this project has
shipped**. With `lam = 0` the union basis's smallest eigenvalues (~1e-11) are
amplified without bound, which is exactly why AGENTS.md 4.2 records that
whitening makes `|predicted effect|` a WORSE per-gene confidence signal. A ridge
floors those directions, decoupling the two jobs `gamma` does at once:
suppressing dominant shared programmes (which buys `pds`) and amplifying the
near-null tail (which costs per-gene sign). That single mechanism retro-explains
both 153 §8's monotone reach-vs-gamma and 088's `gamma = 0` collapse.

**Arm.** a15 exactly, plus `gamma 0.65 -> 0.80` and `lam 0.0 -> 0.007`. Screened
in exact scorer space on two held surfaces and selected by **worst-surface**
value, not by average — positive on `pds`, `reach` and `nmae` on both, with no
member regressing on either.

**Anti-artefact check.** Head purity rose on both surfaces
(`0.6646 -> 0.6845`, `0.5625 -> 0.5697`) while `reach` and `pds` also rose. The
artefact signature is purity DOWN while reach goes UP, so this is its opposite.
Two caveats on the record: median `k*` is unchanged (`3.0`, `2.0`), so the gain
is entirely in the mean and tail-weighted; and `k*=0` worsens `22 -> 23` and
`29 -> 32`, deterministic rather than seed noise. If the board underdelivers on
`reach`, those two lines are the first place to look.

**Two builder defects found and fixed.**

1. `assert_frozen_contract()` compared `EXP062.WHITEN_GAMMA` against the main
   `WHITEN_GAMMA`; they were equal only by coincidence in every shipped build, so
   sweeping the main gamma tripped the guard. Decoupled with a pinned
   `CD4_WHITEN_GAMMA = 0.65`, which is also the correct semantics — the CD4
   residual must not move as an unmeasured second variable.
2. **`--whitener-id` silently defaults to `union`, but a15 shipped
   `union_pruned`.** Three whitening operators exist on disk with near-identical
   names. The first three arms were built on the wrong basis and would have been
   two-variable changes sold as one. Caught only by diffing the finished meta
   against a15's and finding two config differences where one was intended. Now
   AGENTS.md 11.1 rule 10.

**Verification before the slot was spent.** Config diff vs a15: exactly
`2 of 46` fields, both intended. Operator sha `e2f65ecf...` identical to a15's.
Operator-scale falsifier: `after_all_median_norm` `4.22930 -> 4.23539`, a
**0.144%** drift against a 5% threshold, so alpha keeps its meaning. The ridge
also *reduced* operator amplification `25.74 -> 22.03` with the renorm
compensating — independent confirmation of the flooring mechanism.

**The projected gain is NOT supported, and the record must say so.** The
(gamma, lam) screen that selected this cell un-whitened experiment 104's frozen
baseline and re-whitened it at each cell. That baseline was built with a
**fourth** operator — the leakage-safe held-HCT116 whitener
(`experiments/062-cd4-post-whiten/cloud-output/whiteners/held_orion_hct116.npz`,
sha `7cecbf7f...`), recorded in `104/work/hct116/manifest.json` but absent from
its `input-audit.json`. The screen un-whitened with a different operator, so its
`U_q` is not the raw effect but a composition of two operators, and the sweep
does not measure "production rebuilt at (gamma, lam)". The identity check at
`(0.65, 0.0)` reproduces the baseline exactly by orthonormality, which is
precisely why it could never detect this.

What survives is narrower and worth stating exactly:

- **The mechanism is confirmed on the shipped operator**, independently of that
  screen, by this build's own meta: `after_operator_median_norm 25.7379 ->
  22.0309` on `union_pruned`. The ridge does reduce amplification on the real
  path. That is the intervention working; it is **not** evidence the score
  improves. `gamma = 0` also changes the operator in a well-understood way and
  088 measured it as a collapse.
- The family replicated in sign and monotone shape across two genuinely
  different operators (max eigenvalue `8.036` and `5.895`), with the optimal
  `lam` tracking the spectrum in the predicted direction.
- `0.80 / 0.007` is the smallest operator move of any candidate considered.

**Re-measured correctly, and the projection is restored.** The screen was rerun
with each surface un-whitened and re-whitened using the operator its own baseline
was actually built with (062's held-HCT116 / held-HEK293T). a19 is **positive on
all six members on both surfaces**, worst-surface `+0.00195`. The best
all-positive cell is `0.90 / 0.015` at `+0.00215` — a gap of `+0.0002`, far
inside what one board point resolves, so shipping `0.80 / 0.007` costs nothing
measurable.

**Preregistered board prediction: `~ +0.002` overall**, carried by `pds` and
`reach` both slightly positive with `nmae` slightly better (HCT `+0.00049` pds /
`+0.00109` reach; HEK `+0.00152` / `+0.00018`). **If the board lands materially
above `+0.004` or below `-0.001`, the local instrument is wrong about this lever**
— and that finding would be worth more than the arm.

**The optimum did not move across three genuinely different operators.** 095's
held-HCT (max eigenvalue `8.036`), `union_pruned` (`5.895`), and each surface's
own held operator all put the all-positive region at **gamma `0.80-0.90`, lam
`0.007-0.015`, worth `~+0.002` worst-surface**. So the (gamma, lam) family is
*not* operator-specific in the way the `n_real` regime problem is, and it should
**not** be filed next to that rule. Tuning it on a harness operator is
legitimate; only the absolute magnitude needs the caveat below.

**Standing methodological consequence.** A screen run on a leakage-safe operator
measures a different object than the one that ships, and they disagree *by
construction* — the held-safe operator excludes the very surface being scored.
There is no single number available for both. Report on the held-safe operator
and state the gap; never treat a harness-operator (gamma, lam) result as a
production claim.

That sweep was run and comes back clean: **eight distinct operators exist across
`experiments/`, all resolvable on disk, no orphans** — `5608e807` (`union`),
`e2f65ecf` (`union_pruned`, ships), `5cf3899d` (095-local), `7cecbf7f` /
`97e486e0` / `13d988fc` (062 held-safe HCT116 / HEK293T / K562), and `028645b9` /
`5684a5b1` (106's **pruned** held-safe pair). The only unresolved `operator`
strings in the corpus are prose in preregistration documents ("essential-pruned
held-HCT union operator") rather than paths — record a path and a sha, never a
description.

**Use 106's pruned held-safe pair for any future screen.** It is pruned like
production's `union_pruned` *and* excludes the scored line, making it the closest
available analogue to the shipped basis for a held surface — strictly better than
062's. It could not be used here only because un-whitening must use whatever
operator did the whitening, and 104's baseline used 062's.

Submitted as **"The Ridge Objects To Infinity"**, entry `QGpNVMwzzdSlhUoQtRdp`.

## BOARD RESULT: terminal NO-GO, `-0.01352`, and the failure is instructive

| member | a15 | a19 | delta |
| --- | ---: | ---: | ---: |
| `pds` | `0.745871` | `0.744688` | `-0.00118` |
| `mse` | `0.120057` | `0.024149` | **`-0.09591`** |
| `nmae` | `0.094239` | `0.102690` | `+0.00845` |
| `fid` | `-0.030397` | `-0.026398` | `+0.00400` |
| `reach` | `0.115368` | `0.118338` | `+0.00297` |
| `jac` | `-0.000619` | `-0.000053` | `+0.00057` |
| **overall** | **`0.174087`** | **`0.160569`** | **`-0.01352`** |

Rank 12 → 26. Predicted `+0.002` with a preregistered falsifier at `-0.001`; delivered
`-0.01352`, so the falsifier fires decisively.

**Every DE member moved as predicted.** `reach` beat prediction (`+0.00297` against
`+0.00109`/`+0.00018`), `nmae` and `fid` and `jac` all improved, `pds` was near-flat and
missed its sign by `0.0017`. The local two-surface screen was **right about everything it
measured**.

**It lost on `mse`, which the screen never measured.** No mse column existed in the
(gamma, lam) grid on any operator at any cell. A ridge on the whitening operator changes
the eigenvalue floor and therefore the effect vector's *shape*, and `mse` reads the profile
directly — so the one member with no local estimate was the one carrying the whole result.

**Standing consequence, now AGENTS.md 11.1 rule 12.** An arm must be priced on **all six
members**, or must state explicitly which are unmeasured *and why they cannot move*. That
exemption is legitimate for `spread`/`phi`/`keep_top` arms — the quasi emitter's exact-total
repair guarantees `pds`/`mse`/`nmae` invariance by construction, and v14 verified `pds` moved
by exactly `0.000000`. It is **not** legitimate for an operator change, which has no such
guarantee. a19 was certified against a three-member instrument as though the six-member
guarantee applied.

**a15 restored on slot 2**, byte-identical (`55f7bd41...` verified against the recorded
remote hash), returning the board to `0.1740865087`.

---

## 2026-08-31 — a15 restored, and the board scorer is certified deterministic

Entry `5QfX14axLAG8FF2rSDJe`, a **byte-identical** resubmission of `a15-cd4f010.vcc`
(sha `55f7bd41...`, verified against its recorded remote hash) after a19's `-0.01352`.
Rank restored to **12**, overall back to `0.1740865087`.

**Free known-answer row: the scorer is deterministic.** Resubmitting identical bytes
reproduced every member to float representation:

| member | original a15 | restore |
| --- | --- | --- |
| overall | `0.1740865087184745` | `0.17408650871847453` |
| `pds` | `0.7458708290005253` | `0.7458708290005253` |
| `mse` | `0.1200566151966586` | `0.1200566151966586` |
| `nmae` | `0.09423877808887045` | `0.09423877808887053` |
| `fid` / `reach` / `jac` | identical | identical |

Only `nmae` differs, in the 16th significant digit — float accumulation order, not scoring
noise.

**Consequence: every board delta is signal.** There is no run-to-run scorer variance to
subtract, so a board difference of any size is real and attributable to the artifact. This
retires a standing worry — the project has repeatedly had to ask whether a local effect was
inside noise, and the same question does not arise for board points. a19's `-0.01352` was
entirely real, and any future arm's delta can be read at face value.

Worth spending one slot on once, and it cost nothing extra: the restore was already required.

---

## 2026-08-31 — a21: the resampling emitter, built and probed (queued for the UTC reset)

**The lever.** `donor_transfer.py:299` sets `L = int(min(median(lib_sizes), max_counts))` —
one fixed library size for all 400 emitted cells, while the challenge's own controls span
`3,275-52,420` (CV 0.46). Emitted CPM therefore lives on a `49.7 CPM` lattice the controls do
not, and against the scorer's Mann-Whitney that is **not null under the null**: measured in
closed form on contexts A/B/C, ~19% of gated genes clear BH at **zero predicted effect** and
~88% of those are called UP. The builder already samples 400 real library sizes
(`make_submission_pm_cd4.py:1977`) and `quasi` discards them.

**The arm.** New `emit_counts(mode="resample_mix")`: draw real control cells, apply the
predicted per-gene fold change multiplicatively, renormalise each cell to **its own** library,
then repair per-gene totals **after** mixing. `VCC_RESAMPLE_X=0.14` keeps the top 14% of genes
by `|predicted effect|` on the parametric law and resamples the rest.

**Repair ordering is load-bearing.** The two halves carry different per-cell library structure,
so splicing genes between them breaks the pseudobulk pinning: rel L2 is `0.0155` repairing
before the mix and `0.0015` after. A 1.55% profile drift is the exact channel that cost a19
`-0.0959` on `mse`.

**Guards, all passed.** `mode="quasi"` verified **bit-identical** to the pre-change module at
phi 0 and phi 4 (161,937 / 152,647 nonzeros matching exactly) — the frozen-implementation hash
was re-registered deliberately on that evidence, following the 2026-08-27 `overdispersion`
precedent. Per-gene totals exact at every mix fraction. Emitted-pseudobulk rel L2 `0.001387`
against a `0.002` ceiling. Config diff vs a15: **1 of 46** fields. Operator sha matches. 36
builder tests green. Default `VCC_RESAMPLE_X` unset leaves the mode unreachable.

**Independent anchor reproduction.** At `x=1.0` (all parametric) the emitted per-gene CV is
`0.668` of the real controls', reproducing this repo's documented `0.67x` from 2026-08-28 and a
peer session's jackknife `0.656`; at `x=0.0` it is `0.963`. Three implementations agreeing.

**Probe on the built artifact** (`probe_de.py`, real challenge controls, 12 perts/context):

| ctx | a15 `n_pred` | a21 `n_pred` | a15 %UP | a21 %UP |
| --- | ---: | ---: | ---: | ---: |
| A | 4640 | 1259 | 60.8 | 82.1 |
| B | 4459 | 999 | 63.3 | 89.3 |
| C | 4162 | 1164 | 63.3 | 84.2 |
| mean | **4420** | **1141** | ~62 | ~85 |

`n_pred` 1,141 sits between the two candidate `n_real` values (~240 / ~2,500), so `fid`
discriminates them by roughly `0.28` raw.

**The decile table is direct evidence the mechanism works on a real package**, not in a model:

```
a21 ctx A, frac UP by confidence decile (1 = most confident)
  0.935  0.306  0.150  0.137  0.152  0.178  0.214  0.253  0.307  0.392
```

The discontinuity sits exactly at the `param_mask` boundary — x=0.14 of ~10k gated genes is
~1,400 against `n_pred` 1,141 — so the surviving call set is essentially the parametric head
and **the resampled tail went genuinely quiet**. The artifact is also 16.9 GB against a15's
20.1 GB, and the packaged `.vcc` 3.21 GB against 3.69 GB.

**The up-rate rising 62% -> 85% is not the risk it first appears.** BRIEF lines 232-234 count
`k` over the genes we call that are **adjudicable**, not over the reference's DE set. With
`n_pred >> n_real` most of our calls land on reference-*non*-significant genes, which
`CURRENT_STATE:741` records as drifting **54.8% up**. Under no skill, expected precision is
`u*q + (1-u)*(1-q)`, so `u: 0.62 -> 0.85` moves precision `+0.015 to +0.027` — the *opposite*
direction — and a15's raw `fid` of `0.5029` sits on that no-skill line. The unfavourable case
(calls concentrating on DE-set genes where `q ~ 0.40`) bounds the cost at `-0.006`. So the
confound is `-0.006 to +0.027` against a `0.28` coverage signal: a **14x margin**. This is
explicitly not the a19 shape, where the unmeasured member was seven times the predicted effect.

Both of 2026-08-31's slots are spent (a19, then the a15 restore). a21 is prepped and queued
for the first slot after the 00:00 UTC reset.

---

## 2026-09-01 — The board scorer has NO run-to-run variance

Stated separately from the a15 restore that produced it, because it is a durable fact about
the instrument and must not be read as a justification for that restore. The restore was a
mistake — see the slot rule in `AGENTS.md` §6 — and this was an accidental byproduct of it.

Resubmitting byte-identical bytes (`a15-cd4f010.vcc`, sha `55f7bd41...`) reproduced every
scored member exactly:

| member | first scoring | second scoring |
| --- | --- | --- |
| overall | `0.1740865087184745` | `0.17408650871847453` |
| `pds` | `0.7458708290005253` | `0.7458708290005253` |
| `mse` | `0.1200566151966586` | `0.1200566151966586` |
| `nmae` | `0.09423877808887045` | `0.09423877808887053` |
| `fid` / `reach` / `jac` | identical | identical |

Only `nmae` differs and only in the 16th significant digit — float accumulation order.

**Consequence: subtract no noise allowance from any board delta.** Every local instrument this
project uses needs one — the over-dispersion arms carry ~`0.012` of `reach` seed noise, and a
(gamma, lam) screen changed *sign* under re-derivation on 2026-08-31. The board carries none.
A board difference of any magnitude is signal and is attributable to the artifact.

a19's `-0.01352` was therefore entirely real, and any future arm's delta can be read at face
value without asking whether it clears a noise floor.

This is a one-time measurement. It is now known, so it can never again be a reason to
re-submit a scored package — and `scripts/check_submission_novelty.py` refuses to.

---

## 2026-09-01 — a21 resampling emitter: `-0.02146`, and it settles `n_real`

Entry `5Zvx3QDyb6nJMKXS0h5C`, "Every Cell Gets Its Own Library". Overall `0.1526293`
against a15's `0.1740865`. Rank 38.

| member | a15 | a21 | delta |
| --- | ---: | ---: | ---: |
| `pds` | `0.7458708290005253` | `0.7458708290005253` | **`0.00000000`** |
| `mse` | `0.1200566151966586` | `0.12005661519665851` | **`0.00000000`** |
| `nmae` | `0.094239` | `0.030025` | `-0.064214` |
| `fid` | `-0.030397` | `-0.059861` | `-0.029464` |
| `reach` | `0.115368` | `0.096886` | `-0.018482` |
| `jac` | `-0.000619` | `-0.017201` | `-0.016583` |

### 1. The emission design is validated exactly

`pds` and `mse` are unchanged **to sixteen significant digits** while the entire cell channel
was rebuilt from real control cells. The repair-after-mixing order did precisely what it was
designed to do: emitted-pseudobulk rel L2 `0.0014` locally, and bit-identical profile members
on the board. Any future emission arm can rely on this construction to leave the profile
channel untouched.

### 2. H_B is CONFIRMED — the challenge's `n_real` is small

This is what the slot bought. With `n_pred` probed at `1141`:

```
H_A (n_real ~2500): coverage 0.456 -> fid raw ~0.224 -> scaled ~ -0.966
H_B (n_real ~240):  coverage 1.000 -> fid = precision -> scaled ~ -0.078
OBSERVED                                                          -0.0599
```

The observation sits on H_B's prediction and **16x away from H_A's**. Coverage did not
engage, so it was already saturated at `n_pred` 1141 — and therefore at every larger value.

**Consequence: the call-set channel is closed.** There is no coverage to recover, `jac` is at
its ceiling for our `n_pred`, and trimming or inflating the call set cannot buy `fid`. This
retires a question that produced three wrong board predictions (v14's `fid`, v15's `fid`,
a8's `jac`) and that two sessions spent 2026-08-31 arguing over. The `~2,670` figure in
`093-spread-truncation/BOARD-RESULT.md` was circular — it inverted board `fid` using a
*harness* precision — and it is now refuted directly.

### 3. Why it lost: the cell channel is real, movable, and negatively signed for us

`nmae` fell `-0.0642`, the largest single-member move this project has ever engineered.
`nmae`/`fid`/`reach`/`jac` read the **unweighted mean of per-cell CPM**, not the pseudobulk.
`donor_transfer.py:299`'s fixed `L` welds that cell channel to the profile channel; resampling
unwelds it, and every emitted cell then carries its own real library size.

So the channel exists and the lever is large. Its **sign** is the problem. The payoff scales
as `(2*accuracy - 1)`, and our pooled per-gene sign accuracy is **0.522** — a coin flip. a21
applied the lever to the resampled **tail** (86% of genes, sign accuracy ~0.51), i.e. to
precisely the slice carrying no signal, so it amplified noise.

**The untested inverse:** top-decile sign accuracy is **0.7359**, giving `(2*0.736 - 1) = 0.47`
against the pooled `0.044` — a 10x better multiplier. Applying the cell-channel change
**head-only** is the arm this result points at, and it is the same conclusion the independent
Aginglab analysis reached from the metric algebra.

Per the slot rule, a15 was **not** restored. a21 stays displayed.

---

## 2026-09-01 — a22 head-only resample: the cell-channel sign law is CONFIRMED

Entry `EFf5zr3H1PEbJW0imuQo`, "The Head Gets Its Own Library". Overall `0.1715689`,
rank 19. Deltas against a15:

| member | a21 (tail resampled, acc ~0.51) | a22 (head resampled, acc 0.7359) |
| --- | ---: | ---: |
| `pds` | `0.00000000` | `0.00000000` |
| `mse` | `0.00000000` | `0.00000000` |
| **`nmae`** | **`-0.064214`** | **`+0.001520`** |
| `fid` | `-0.029464` | `-0.004447` |
| `reach` | `-0.018482` | `-0.010656` |
| `jac` | `-0.016583` | `-0.001523` |
| overall | `-0.021458` | `-0.002518` |

### The confirmed law

`nmae` reads the **unweighted mean of per-cell CPM**, and its payoff under a cell-channel
change scales as `(2*acc - 1)` in our per-gene sign accuracy. Two arms, same machinery,
opposite halves of the gene set:

```
a21  resampled the tail, acc ~0.51, multiplier 0.044  ->  nmae -0.064214
a22  resampled the head, acc 0.7359, multiplier 0.47  ->  nmae +0.001520
```

**The sign flipped as predicted.** This is the first mechanism this project has predicted in
advance and then confirmed on the board with a matched pair. It also confirms the two-channel
reading: `donor_transfer.py:299`'s fixed `L` welds the cell channel to the profile channel,
and unwelding it moves `nmae` while `pds` and `mse` stay bit-identical — which they did, to
16 significant digits, for the **third** consecutive arm.

### Why it still loses: reach and the ordering

`nmae` gained `+0.0015` but `reach` lost `-0.0107` and `fid` `-0.0044`, so the arm is net
`-0.0025`.

The mechanism is clear and it is a genuine tension rather than a bug. `reach` walks the
reference's DE pool **in our confidence order**, and that order comes from the emitted
p-values. Resampling a gene changes its emitted distribution, which changes its p-value,
which perturbs its rank. a22 resampled precisely the **head** — the top of the very ordering
`reach` reads — so it bought `nmae` calibration by adding noise to the ranking that `reach`
depends on.

So the cell-channel lever cannot be aimed at the head without disturbing the head's order.
`nmae` wants the change where sign accuracy is high; `reach` wants that exact region left
undisturbed. Within this family those requirements are in direct opposition, and the
magnitudes say `reach` wins the argument: `-0.0107` against `+0.0015`, a 7x margin.

### Consequence

Do **not** sweep `x` in this family expecting a positive cell. At `x -> 0` the arm converges
to a15 with zero change; at `x = 0.10` it is `-0.0025`; at `x = 0.14` inverted it is
`-0.0215`. The derivative at the origin is negative and the `reach` term dominates the `nmae`
term everywhere it has been measured.

A successor would need to move the emitted **magnitudes** without moving the emitted
**p-value order** — the two are coupled through the same emitted counts, so it is not obvious
that any emission-side change can do it. That is the honest open question this pair leaves.

Both 2026-09-01 slots spent (a21, a22). Per the slot rule a15 was not restored; a22 stays
displayed at rank 19.

---

## 2026-09-02 — a23: the renderer was undoing the geometry, and every screen since exp 017 scored the wrong object

**The finding.** `profile_from_effect` renders the log-space effect as `expm1(b_ctrl + gated)`
and divides by its own sum. The scorer reads back `log1p(TS*p) - b_ctrl`. That division is the
identity **iff** `sum_g clip(expm1(b_ctrl + h_g), 0) == TS`. The control profile satisfies it
exactly by construction (`sum(expm1(b_ref)) = 50000.0` verified, control round-trip
`||b_recov - b|| = 0.0`), and so does every reference perturbation — **the truth lives on that
manifold.** Our effect does not: the CPM gate, the alpha scaling, the whitening and the CD4
addition each break it.

So the renderer applies its own correction, and it lands as a shift concentrated on the
highest-expressed genes that is nearly the same direction for every target:

| | HCT116 | HEK293T | challenge A/B/C |
| --- | ---: | ---: | ---: |
| mass ratio | 1.00181 | 1.00168 | 1.00134 / 1.00135 / 1.00127 |
| render residual, % of row norm | 3.96% | 4.53% | 2.83% / 3.35% / 2.94% |
| row-row cosine of the residual | **0.303** | — | — |
| row-row cosine of the intended effect | **0.0040** | 0.0035 | — |

**The contamination is 75x more mutually-aligned than the signal it is added to.** `pds` is a
retrieval metric, hubness is what it punishes, and whitening's entire job (§4.2) is to make
those rows mutually dissimilar. The renderer was putting the shared blob back.

**Consequence for the instrument, which is worth more than the arm.** `identity/intended` beats
`identity/delivered` by **+0.00324 / +0.00382** raw pds. Every effect-level screen in this
project since experiment 017 — every gamma sweep, every donor comparison, every LODO fold —
scored the effect *before* the renderer touched it. The board scores the effect *after*. They
are different geometry.

**The fix.** `expm1` is convex with `S(0) = TS`, so `S` has no second root along `{gamma*h}` —
scaling cannot satisfy the constraint and the correction must move off the ray. But it is one
scalar, so the direction is free: `dS/dh_g = exp(b_g + h_g) =: w_g`, and the minimum-L2
correction is `eps = -(S - TS) * w / ||w||^2`, confined to the CPM gate. Newton with a
minimum-norm step, `VCC_MASS_FIX=minnorm`, default off and bit-identical when unset, 40 tests
green. Cost 1.76% / 2.29% of the row norm; the read-back becomes the identity to `1e-12`.

**Measured, both held surfaces, all six members in our favour:**

| | HCT116 | HEK293T |
| --- | ---: | ---: |
| Δ pds | +0.00035 | +0.00278 |
| aggregate cosine | 0.02285 -> 0.02596 (**+13.6%**) | 0.00900 -> 0.01404 (**+56.0%**) |
| sign accuracy on the reference DE pool | 0.5311 -> **0.5417** | — |
| head purity | 0.6751 -> 0.6844 | — |
| nmae (lower better) | 0.9501 -> 0.9479 | — |

Of the **318** (pert, gene) pairs in the reference DE pool whose predicted sign moved,
**68.9% land on the correct side** against 0.500 chance. The fix removes a systematic bias.

**Three internal controls, and they are the reason to believe it.**

1. `massfix-const` satisfies the **identical** scalar constraint with a uniform log shift and
   **loses** 0.0056/0.0031 pds and 4.7%/18.8% cosine. Same constraint, different direction,
   opposite sign of result — geometry, not bookkeeping. Kept in the shipped code as that control.
2. `minnorm` beats the *infeasible* `intended` upper bound on HEK293T (+0.00278 vs +0.00237).
   The feasible object is closer to the truth than the one we designed, which is what you expect
   if the truth's own manifold is the constraint set.
3. Seven correction directions (`weighted.py`) all land at +14% / +65% cosine, spread smaller
   than the between-surface disagreement. **There is no placement to tune** — which is the
   signature of a mechanism rather than a fit. Do not spend a slot on another member.

**It also explains a competitor's published board failure.** `centroid` alone, pre-render, is
worth exactly **+0.00000** on HCT116. Duy's board description reports cross-target centring at
"+0.041 scaled perturbation discrimination post-hoc on the pseudobulk delta" but "a seed-matched
local render returns only -0.002". Same mechanism: **the renderer undoes a pre-render geometry
transform.** Geometry applied *after* the constraint is imposed does stack
(`massfix + centroid d1.0`, worst-surface +0.00237) and is held as the clean follow-up.

**The anti-share family is spent for us.** Five teams report gains from removing the shared
component of the 300 predicted rows (AgingLab_BI `L <- L - mean_t(L)`, pds 0.7573->0.7961
through the real scorer; Heisenberg rank-1, 0.6518->0.6829; Zhou Li; Cameron S; Lingkuan Meng).
Our panel matrix already has 0.35% shared energy and mean pairwise cosine +0.0009, against the
0.0747 Heisenberg's operator was fixing — `Donor.decompose()` already mean-removes and
production runs `beta = 0`. **But AGENTS.md §4.2 point 1's stated REASON for banning
panel-fitted operators is wrong on the metric's own algebra**: `k_p` ranks the 300 *measured*
effects against one fixed prediction, so decorrelating predictions cannot improve a ranking "by
construction". The real hazard is a rank-deficient 259x7,760 covariance. The ban is over-broad
as written; the lever is nonetheless empty for us.

### A verification that failed, recorded because the failure is the lesson

The first artifact check recovered `delta_hat` from the emitted counts and tested
`sum(expm1(b_ctrl + delta_hat)) == TS`. a23 read **exactly 1.00000000** — and it is worthless:
`expm1(log1p(TS*p)) = TS*p`, which sums to `TS` for *any* normalised profile. It was caught only
because the same check was run on **v9**, a pre-fix build, which also read exactly 1.00000000.
Rule 10's corollary, operating as written.

The real check rebuilt the intended effect matrix independently (`post_whiten_sha256` confirmed
identical to a23's `b2426e15…`) and compared the residual the artifact carries against both
candidate corrections: it points at `eps_minnorm` on **12/12** probed targets, mean cosine
**+0.700** against **+0.338** for the renderer's implicit shift.

### Also measured tonight

**Per-target amplitude is a confidence problem, not a magnitude problem.** `mse` is a ratio of
sums and `pds` is exactly invariant to a per-row scale, so per-target amplitude is a free lever.
Oracle per-target scaling takes aggregate cosine 0.0260 -> 0.0413 (HCT116, **+59%**) and
0.0140 -> 0.0274 (HEK293T, **+95%**). But **oracle MAGNITUDE only** — knowing the true
`||delta_p||` exactly and nothing about the per-target cosine — is worth **+0.000011 /
+0.000005**, i.e. nothing. And **18.9% / 27.3% of targets have an optimal scale <= 0**: our
direction is anti-correlated with the truth for a quarter of all targets. This refutes the
COMPASS-derived hope that magnitude, being the more transferable half, is the exploitable half.
`corr(s_p*, ||d_p||)` is only +0.13. Experiment 158 is screening build-time keys.

**Cross-line transfer, measured directly in the scorer's space** (120 targets shared by the two
Orion lines, all target coordinates dropped): SAME target in a different line, cosine
**+0.1881**; DIFFERENT targets across lines, **+0.1790**. Target-specific excess **+0.0091**.
Within-line different-target cosine is +0.2476 (HCT116) and +0.4175 (HEK293T), so the truth
itself is dominated by a shared programme. Energy split over the same targets: **60.1% shared,
39.9% line-specific**, and a perfect *context-independent* predictor caps at cosine 0.77 to each
line. This reproduces, from cosine, the +0.021 sign-accuracy excess measured independently on
2026-09-01 — two statistics, one conclusion. Both are attenuated by 150-cell sampling noise and
the attenuation is not uniform; experiment 159 is disattenuating it.

### a23 as submitted

`submissions/a23-massfix.vcc`, sha256 `525574729da0e20f84e3c4cb56bbe55c5a0f1e6216e780dcbc6ee875dff20445`.
Config diff vs a15: **one semantic field** (`mass_fix: "minnorm"`), plus two null schema keys
(`resample_x`, `resample_invert`) that postdate a15's meta. Operator `union_pruned`,
`columns_sha256` and `metadata_sha256` identical to a15's. Whitening audit values match a15 to
**1e-16** — a15 was built on AWS Linux x86_64 and a23 on macOS arm64, so the shas differ by
architecture float noise; a21 and a22 were also macOS builds and reproduced Linux-built a15's
`pds` and `mse` to sixteen significant digits on the board, which settles that this cannot move a
member.

**Pre-registered prediction (`A23-PREDICTION.md`, amended on the challenge-frame measurement
before the slot was spent): +0.004 overall, range +0.003 to +0.011, landing ~0.178 against
a15's 0.1740865087.** Falsifiers: raw `mse` must fall below 0.8680; scaled `pds` must stay in
[-0.001, +0.005]; overall below a15 is a NO-GO and would mean an "all six positive on two
surfaces" profile-channel screen is untrustworthy — which would be worth more than the arm.

### 2026-09-02 06:20 UTC+3 — gamma re-opened as an `mse` lever (exp 170 launched)

Assembled from three existing facts: a19's board point (`mse` −0.0959 vs `pds` −0.0012 for
gamma +0.15, i.e. 81× the sensitivity), 097's measurement that lowering gamma raises `kappa`,
and the observation that 097 closed gamma without an `mse` column. The veto on lowering came
from the Orion pair, which 168d discredited for sizing. Naive back-extrapolation reaches raw
`mse` 0.78 at gamma 0.50 (leader is 0.7767) — treated strictly as a hypothesis. Exp 170 sweeps
gamma ∈ [0.35, 0.80] at `lam = 0` on h1-test through 096's `claim` route with the a15 audit
values as a known-answer gate, and prices all six members. Hold stands.

### 2026-09-02 06:25 UTC+3 — 161 interim: residual R² = 1.7%, not predictable

Five-line surface, three known-answer rows pass. The held-line residual is explained at under
2% by recipient control-profile features and consensus×context interactions; delivered gain is
`pds` +0.005 with cosine falling (hubness pattern, not prediction). Agent died pre-RESULTS.md;
watchdog resuming for the write-up.

### 2026-09-02 06:35 UTC+3 — CORRECTION to the 06:20 entry, and the board re-read

**The gamma re-opening rested on a false premise and I killed exp 170.** 097 §9h *does* price
`mse` across gamma on both surfaces. Direct h1-test slope `d(raw mse)/d(gamma) = +0.0947` against
a19's implied +0.6075 — **a19 overstates 6.4×, which quantifies the `lam` confound: the ridge did
~6/7 of a19's −0.0959 `mse` loss, not the exponent.** Gamma's h1-test optimum is +0.0027 at 0.30,
below rule 13's +0.0089 operator gate. Recorded as sub-threshold live, not closed.

**Board re-read (581 teams, leader `cqawesome` 0.2169, us 25th at 0.1734).** Per-member rank:
`pds` **5**, `mse` 23, `nmae` 68, `jac` 169, `fid` 176, **`reach` 192**. Reach is 37.3% of the gap
and worth +0.0162 to match the leader. Launched exp 171 on the hypothesis that our significant
calls take pole position in the reach walk carrying zero directional signal, dragging the
cumulative `P(k)` to chance at the front where it cannot be undone.

### 2026-09-02 06:45 UTC+3 — reach: the deficit is a DEFECT, quantified against the field

Field regression over 581 teams: `raw_reach = −0.0751 + 0.3577·raw_pds`, corr +0.831, resid sd
0.0409. We sit **−1.22σ below the line**, bottom 6.9%. All four teams with `raw pds` ≥ ours reach
≥ 0.2492 against our 0.1742 — and two of them (Slurm Shadi 0.883, cqawesome 0.762) hold our raw
`mse` or better, so **amplitude/shrinkage is not the explanation.** Board intelligence written to
`experiments/171-reach-prefix/INBOX-board-intelligence.md`, including Evidoa's published
single-knob amplitude probe (MAG 0.77→1.66, reach 0.2272) which 171 must treat as a control rather
than the hypothesis. Confirmed our own public row is clean: model name is a joke, description is
`None`.

### 2026-09-02 06:50 UTC+3 — the unifying pattern: our DE CALL SET is the broken subsystem

Members that read our call set (`fid` 176, `jac` 169, `reach` 192) versus members that do not
(`pds` 5, `mse` 23, `nmae` 68), out of 581. No overlap; the boundary is exactly whether the metric
passes through our emitted cells' Wilcoxon calls. Field-best headroom on the three call-set members
is **+0.0339 overall** — 0.1734 → 0.2073 — with no improvement to the direction. Retro-explains
a22 (`reach` −0.01066), a23 (`reach` = 182% of the loss) and 093's board inversion as one broken
subsystem rather than three accidents. 171's mandate broadened from "reach" to "diagnose and repair
the DE call set", with reach as the largest and best-instrumented symptom.

### 2026-09-02 06:55 UTC+3 — 161 CLOSED NO-GO; h1-test validated on the DE members; the defect measured

161: cross-line residual not predictable. Context-only model reads `R²` −0.0024 and is
indistinguishable from its permutation null; adding recipient features *lowers* cosine (T 0.0606 >
X 0.0592 > ALL 0.0550); the surviving gain is a zero-learning `sign(C)·|C|²` re-shaping worth 3.1x
the learned model. Fifth route into the 159/160 wall. Its q=2 result handed to 169 with the
target-independent-energy discriminator, since 168d closed that family on h1-test.

**169's known-answer row incidentally validated h1-test on the DE members:** `KA_a15_raw` reads
reach **0.18283** and fid **0.50701** against our live board **0.17422** and ~0.505 — within ~0.009
raw. And it measured the defect directly: **our call set's sign accuracy is 0.53435**, precision
0.50701. The genes we call DE are barely distinguishable from coin flips while the profile they
come from ranks 5th in the world. Both numbers written to 171's INBOX 2.

### 2026-09-02 07:00 UTC+3 — the call-set mechanism solved, and a decoupling nobody has exploited

Two-population solve of 169's `n_pred` 6034 / `sign` 0.53435 gives **~900 informative genes at
0.730 accuracy plus ~5,100 at chance** — and 168b (0.698–0.731 for the top 5%) and a22's board probe
(top decile 0.7359) agree independently. Our signs are fine; 85% of the call set is noise.

`fid`'s yield cap `min(1, n_pred/n_real)` is what pins `keep_top` at 8000 despite the builder's own
docstring recommending 850, and it re-reads 093's board inversion as a yield-cap effect.

**But `reach` has no yield term, and the order inside our significant block is set by our own
p-values, which the per-gene spread sets.** So one parameter serves two objectives that want
opposite things. The arm is **spread contrast at fixed call count** — keep `n_pred` ~6000 for `fid`,
drive the ~900 informative genes to far smaller p-values so they lead the reach walk. Written to
`experiments/171-reach-prefix/INBOX-3-the-decoupling.md` with its falsifier: if our p-value rank
tracks expression rather than `|effect|`, the arm cannot work and that becomes the headline.

### 2026-09-02 07:20 UTC+3 — falsifier run on disk artifacts: p-value order is EXPRESSION-set

Inside our call set `Spearman(−log10 padj, |effect|)` = **−0.391**, vs expression **+0.485** (a15,
h1-test, no leak: 095 artifact used k562_gwps + orion_hek293t only). Re-ordering the reference's DE
pool by `|H|` lifts reach 0.18283 → **0.23533 (+29%)**; "sig-first then |H| inside" and "|H| over
called genes only" both read 0.21914 exactly — the call set is right, the order inside it is wrong.
+0.0525 raw is our entire deficit below the field pds→reach line. Mechanism: Wilcoxon `z ∝ μh/σ`
runs ~5x across expression at constant spread, and `keep_top` tiers by `|effect|` so cannot touch it
— which is why 093/095/a15 realisations all failed. Fix is an expression-compensating spread
schedule; must be priced on `mse` first (a15 at 99.7% of the dispersion cap). Written to 171's
INBOX 4 with a three-number reproduction.

### 2026-09-02 07:35 UTC+3 — feasibility of the expression-compensating schedule, measured

From `_quasi_block`: `σ²_ours = s·m + φ·s²·m²`, crossover `m* ≈ 1/(φs²) ≈ 200 CPM` at a15's
setting, matchable floor ≈ 12–16 CPM. Ordering the reference's DE pool by `|H|` only above 12 CPM
(production order below, sig-first block preserved) recovers **73.4%** of the inside-block ceiling:
reach 0.18283 → 0.20949, +0.0267 raw ≈ +0.0049 overall as a proxy lower bound. Only 9.3% of the
pool is below the floor. Corrected a direction error in INBOX 4: `s` must INCREASE with `m` (add
noise to high-expression genes to pull the plateau down); `φ > 0` is required for a plateau to exist.

### 2026-09-02 07:45 UTC+3 — 165 CLOSED NO-GO; h1-test blind to reach-on-profile

No head-protected mass-fix placement moves h1-test reach; protection costs the profile gain
monotonically. Instrument finding: a23 (`minnorm`) reads +0.0004 reach on h1-test vs −0.0079 on the
board — h1-test is validated for reach on emission-law changes only. Caveat appended to 171's INBOX 4.

### 2026-09-02 07:55 UTC+3 — oracle ordering measured at 0.646; exp 172 launched

On h1-test, our own signs perfectly ordered give reach **0.64553** against production 0.18283, the
`|H|` key 0.23533 and the leader's 0.26063. The ordering gap (~0.46 raw ≈ +0.085 overall) splits
into delivery (171, emission schedule) and key quality (172, learned sign-confidence key, LODO over
the three production donors, no gene/target identity, scored only as `reach(key=c)` + top-X%
accuracy on h1-test). No composite of `|H|` beats `|H|`. Calibration: `|H|` reaches 0.235 at pool
Spearman 0.069 because its information is top-concentrated — 172 must optimise precision at the top.

### 2026-09-02 08:05 UTC+3 — learnable oracle 0.621; design-sign accuracy 0.537, top decile 0.624

Emission preserves the design sign 94.5% of the time on the reference pool, so the H-oracle
(0.62119) is nearly the full oracle (0.64553): the headroom is learnable in principle. But the
excess over chance is thin (0.537 overall, 0.624 in the top `|H|` decile, bottom half at chance).
Written to 172's INBOX with the decile table; label should be design-sign correctness.

### 2026-09-02 13:55 UTC+3 — 169/171/172 read (all finished ~07:50); all NO-GO

171: "sig calls first" is a BH no-op; spread schedules closed (|z| U-shaped in s, rank-sum reads
median shift; my expression-compensating fix doubled the expression dependence, control beat it);
keep_top 8300 is mse-free (+0.01313 knee headroom) but nets +0.0003. 172: learned keys 0.221–0.243
vs |H| 0.235, inside noise; not learnable from AgingLab's features; whitening costs 0.011 pooled sign
accuracy. 169: sharpening +0.0195 pds, reach toll falsified, but nmae −0.0129 on h1-test (a13
pattern). Hold intact; no slot spent. Remaining problem: per-gene sign accuracy of the predictor.

### 2026-09-02 14:35 UTC+3 — 173 sign repair (coordinator-run): NO-GO at net ≈ +0.002

Whitening-flipped signs (15.7% of the pool) are anti-informative (45.6% right vs cons 54.4%).
Repairing them: fid +0.00175 and nmae +0.00116 resolved; cos +7.8%; pds −0.0042 raw; reach inside
seed noise and negative. Net ≈ +0.002 against the +0.0089 profile gate. Same trade curve as gamma;
composable rider only. Random-flip null collapses everything, so the targeting is real.

### 2026-09-02 14:50 UTC+3 — 174 co-expression: NULL (signal = top-2 PCs = shared programme)

`−sign(r(t,g))` over h1-test controls reads 0.516 on the pool; removing 2 PCs takes it to 0.489.
It is the proliferation/stress axis, not regulation. Recipient control cells closed at the pairwise
level too.

### 2026-09-02 15:05 UTC+3 — 175: p53 donor matching refuted; Orion pair ~0 sign info over K562

K562 alone 0.567 vs equal mean 0.570 vs Orion-only 0.5425 on the H1 pool. p53 targets: p53-WT donor
worst. Per-target oracle donor routing ceiling 0.603. Same pds-vs-sign trade as 097/172/173.

### 2026-09-02 15:35 UTC+3 — synthesis of the day + exp 176 (CMap/LINCS) launched

Every lever today trades pds against the rest; production is at the pds-optimal point; the leader
has a better direction outright. Processing cannot supply it (004/156 close denoising both ways;
161/174 close the recipient's cells). Only new information can. LINCS/CMap — never examined here —
is the one large cross-line perturbation source left. Exp 176 prices it on h1-test with 175's sign
instrument and 173's pds pipeline. Hold stands; no slot spent today.

### 2026-09-02 18:55 UTC+3 — 176 CMap shRNA: NO-GO (coordinator-run; agent auth expired)

CMap shRNA consensus signatures read 0.534 (landmark) / 0.522 (inferred) sign accuracy on the H1
pool vs K562's 0.562/0.568; added to the trio ≤ +0.001 overall. Mechanism: median own-target
knockdown only 0.5–0.9 SD, 25–34% of targets not knocked down. CRISPR-KO modality pending.
**Claude Code needs `claude login` — seven agent launches failed on expired OAuth.**

### 2026-09-02 19:20 UTC+3 — 176 CMap CRISPR-KO measured: landmark-only signal, coverage-capped; LINCS closed

KO 21-line mean adds +0.019 sign accuracy on 867 landmark genes (real transfer) but −0.011 on 7,101
inferred genes; net +0.0006 on covered pairs; ≈ +0.0007 pooled on the challenge after coverage.
NO-GO both modalities. Public cross-line perturbation data is now closed as a route (with 166).

### 2026-09-02 19:30 UTC+3 — Claude Code auth back; 177 (resample emitter) handed to an agent

Auth failed 14:34–18:38 (nine launches/tests), works again at 19:21. Coordinator's first 177 pass hit
a self-inflicted depth confound (real 20k cells rescaled to the harness's 100k convention made
`rs_all` call 10,776/10,779 genes); agent relaunched with the fix (emit resample arms at real depth,
add a quasi-at-real-depth control, zero-effect null first) and the delivery design (symmetric
log-normal jitter schedule so z tracks |H| on an emitter with no median artefact). Motivation: Slurm
Shadi holds our pds and mse with 1.48x our reach, so reach is not bounded by the profile at our
quality — something in emission delivers it, and 171 closed delivery only through our own emitter.

### 2026-09-02 19:35 UTC+3 — 178 (composed riders) launched alongside 177

Gamma + sign repair + K562-heavy weights as one production arm, all six members, against the
+0.0089 profile gate. Composition additivity is the scientific question; a clean NO-GO closes the
rider family. Watchdog tracks 177 and 178.

### 2026-09-02 19:45 UTC+3 — Orion pair shares a protocol axis; per-donor rank-1 deflation +0.0073 sign

HCT116/HEK293T top PCs cos 0.872, K562 0.23–0.31; HEK293T top-1 = 19.7% variance. Deflating k=1 per
donor: consensus sign accuracy 0.5697 → 0.5770, k≥2 hurts. Handed to 178 as a fourth switch.

### 2026-09-02 19:55 UTC+3 — the Orion axis identified: long-gene / nuclear-fraction technical axis

Top-30 loadings are all 0.3–1.7 Mb intron-rich genes (AUTS2, GPC6, NRXN3, PRKG1, WWOX…), orthogonal
to H1's truth (cos ≤ 0.13). A sequencing artefact of the X-Atlas libraries, 19.7% of HEK293T's
variance, inherited by the consensus. Addendum to 178's INBOX with the challenge-controls check.

### 2026-09-02 20:20 UTC+3 — the Orion harness blob IS the long-gene axis; whitening already removes it from H

HEK293T truth: 36.3% of effect energy on the axis (blob cos 0.832); HCT116 13.8%; H1 <2%. Shipped
whitened H carries 0.11% (raw consensus 9.8%). Mechanism of the harness warning and 168d's rule
recorded in AGENTS.md. 178 corrected: deflation acts only via the whitener refit and the sign source.

### 2026-09-02 20:35 UTC+3 — none of A/B/C is a donor line (checked)

Log-pseudobulk Pearson of each challenge context's control profile against donor `b_ctrl`: K562
0.33–0.42, HCT116 0.62–0.67, HEK293T 0.61–0.65, Jurkat 0.36–0.59, RPE1 0.33–0.40, H1 0.69–0.74
(contexts with each other 0.65–0.76). A same-line match would read ~0.9+; none does. Zero-shot is
genuine. Hints only: A×Jurkat 0.587 vs 0.36/0.39 for B/C (A plausibly hematopoietic → production's
CD4 component matters most there); H1's uniform lead is the shared Arc/10x Flex platform.

### 2026-09-02 20:50 UTC+3 — 170 finished (kill failed): gamma TERMINAL NO-GO, 0.65 is a local optimum

Raw challenge mse V-shaped with vertex at 0.65; every other gamma net negative on all six. Mechanism:
ρ cliff — lower gamma raises A/B past claim, refund saturates. 097's h1-test mse estimate assumed
the arm sits on the cliff; withdrawn. Design rule: useful levers must raise X without raising A
(norm-preserving). 178 redirected: gamma dropped from the composite.

### 2026-09-02 22:15 UTC+3 — 177 and 178 read

178: rider family CLOSED (all negative; gains compose sub-additively, mse costs super-additively;
sign repair's pds cost is 6x larger on PLOLO than on h1-test). 177: **the ordering defect is solved
by the emission law** — a15 calls 38–45% of genes significant at zero effect on the challenge
contexts, resampling calls 0.0–0.2%; 81–110% of the |H| reach counterfactual realised with pds
pinned by digest and mse exactly 0. Arm's sign depends on challenge n_real: +0.006 to +0.021 at ~240
(four independent lines including our own board jac), catastrophic at 093's circular 2,670. Rule 13
vetoes on its flat emission-law penalty; agent + coordinator argue the gate is miscalibrated for
this resolution regime. Decision presented to Ahmet under the hold.

### 2026-09-02 22:25 UTC+3 — 179 launched: build hy_kt8000, verify, dry-run, NO SUBMIT

Package for the resample-tail emitter (a23 + `--emitter resample_mix --resample-keep-top 8000`).
Gates: 46+ tests, rule-10 bit-identity of a23's audit, exactly-two-line meta diff, 900/900
pseudobulk digest identity vs a23's package, zero-effect null on the built path, n_pred probe
(prediction ≈ 2,410 vs 4,424), `vcc prep --dry-run`. Decision to spend a slot remains Ahmet's.

### 2026-09-02 22:40 UTC+3 — CORRECTION on 177's expected value

The +0.006 for hy_kt8000 was a single-n_real reconstruction. debench on h1-test's real n_real
distribution: jac −0.0142 raw (−0.0064 overall) from the top tercile only; fid −0.001; reach +0.002.
Board jac (0.030 vs h1-test 0.068) says the challenge tail is ~half as heavy → hy_kt8000 ≈ −0.0006
overall, a wash. The slot's value is the regime information, not score. Range −0.005 to +0.006.

### 2026-09-02 22:50 UTC+3 — per-perturbation K for the resample family is closed

Spearman(||H_p||, n_real_p) = −0.036, (n_pred_p, n_real_p) = −0.085: predicted strength and our call
count carry no information about the true DE count, so the jac loss cannot be targeted. Only a global
K remains; hy_kt8000 ≈ neutral. 179 building for the information slot; decision is Ahmet's.

### 2026-09-02 23:00 UTC+3 — 177 vs 164 armA cross-check

Resampling the tail = 164's tail-φ0 arm on every DE member, minus 164's −0.0185 mse cost (real cells
keep the dispersion credit). Confirms the mechanism and the wash.

### 2026-09-02 23:20 UTC+3 — within-batch Orion rebuild ruled out before launch

Batch control offsets on the long-gene axis are large (sd 4.57) but each perturbation spans ~115
batches (~1.3 cells each), so they average to ~0.43 vs the observed per-perturbation 3.20. The axis
is a real per-perturbation cell-state shift in Orion cells, not batch mismatch. Deflation/whitening
are the right tools; no source-level fix.

### 2026-09-02 23:35 UTC+3 — D/E/F transition plan written into CURRENT_STATE

Builder currently refuses non-A/B/C contexts (line 2617 + hardcoded bundle hashes); alpha is a global
constant while the mse knee is per-context and computable from the new controls (C_p from NTC
jackknife); donor coverage must be re-audited for the new 300 targets; the board is dark Oct 22–Nov 5.
Seven-item checklist recorded; items 1–2 need code + rehearsal well before Oct 22.

### 2026-09-02 23:45 UTC+3 — per-context knee alpha validated arithmetically

A/B/C sit at ρ 0.992/0.940/0.882 under the global 0.982; knee alphas 0.987/1.013/1.047; ≈ +0.0008
overall now, pds-invariant; explains 170's V (global optimum = contexts straddling the knee). D/E/F
robustness item; implementation after 179 releases the builder.

### 2026-09-02 22:45 UTC+3 — hy-kt8000 package built; coordinator's rule-10 meta diff passes

`submissions/hy-kt8000.h5ad` (19.7 GB) built 22:37 with `submission: false`. Against a23's meta: all
54 `effects.*` audit fields identical (incl. `effects.final_sha256`, operator sha `e2f65ecf…`),
config differs in exactly `emitter=resample_mix` and `resample_keep_top=8000`; alpha 0.982, mass_fix
minnorm, keep_top 8000, tail 1.0, keep_spread 0.25, φ 4, min_cpm 10, union_pruned all unchanged.
nnz/cell 6607 → 6837 (51.8% of the density cap). Agent still running: 900/900 pseudobulk digest
identity, built-path zero-effect null, n_pred probe, `vcc prep --dry-run`. NOT submitted.

### 2026-09-02 23:50 UTC+3 — 179 done: hy-kt8000 built, 7 gates (6 as written + mechanism), NOT submitted

Package sha256 e2902589…; pds/mse pinned by 900/900 pseudobulk identity on the artifacts; board n_pred
3,528 vs a23's 4,474 (−21%, not the −46% forecast); tail false-call rate 33–49% → 0.06–1.9%; net
−0.003 to +0.001 overall; rule 13 vetoes. Decision remains Ahmet's; both of today's slots unused.

### 2026-09-02 23:30 UTC+3 — 180 launched: D/E/F readiness (items 1–2 of the transition plan)

Generalise the builder past the hardcoded A/B/C adapter with a re-pinnable frozen-input table and a
synthetic D/E/F rehearsal (A/B/C bit-identity preserved); implement `VCC_ALPHA=knee` (per-context
`s* = √(claim/B)` from the context's own controls), validate against 0.987/1.013/1.047 and ρ = 1.00
on A/B/C, price all six, build the A/B/C package, dry-run only. Hold stands.

### 2026-09-03 01:30 UTC+3 — 180 done: builder D/E/F-ready, knee alpha validated, vcc prep gotcha found

Contexts from manifest; `scripts/frozen_inputs.json` with re-pin; synthetic D/E/F rehearsal PASS with
counts byte-identical to a23's; A/B/C bit-identical (179 fields). `VCC_ALPHA=knee`: 0.98586 /
1.01275 / 1.04568 from the bundle alone, budget ratio ≈1.000, mse +0.000838 overall, pds −0.000076
raw (PLOLO), rule 13 vetoes as an A/B/C arm → D/E/F insurance (asymmetric cost: −0.030 at 30% low,
−0.090 at 10% high). `vcc prep` needs `--contexts D,E,F`. Runbook in 180 §10. Nothing submitted.

### 2026-09-03 01:20 UTC+3 — 181 launched: adversarial audit of every 2026-09-02 closure

Read-only. A fresh agent tries to break each closure (surface, instrument, seeds, arithmetic,
contradictions between documents), runs cheap re-computations against on-disk artifacts, and ranks
each claim SOLID / PROBABLY RIGHT / WEAK / CONTRADICTED with what depends on it. Specific targets: the
pds-vs-accuracy synthesis, the board partition, 177's 81–110% and the 0.46 scaling, the ρ-knee
algebra, the long-gene identification (gene-length correlation), 175 vs 178, 176's coverage
arithmetic, 172's LODO leakage, the n_real inference, and the 6× h1-test-vs-PLOLO pds gap.

### 2026-09-03 03:00 UTC+3 — audit 181 read: hy_kt8000 is DEAD (a21 is its board point, reach −0.0185)

Every NO-GO of the day holds; the one recorded win is refuted. a21 (`resample_x 0.14` = hy_kt≈2600)
was scored 2026-09-01: reach −0.0185, nmae −0.064; h1-test predicted +0.050 reach at that dose.
h1-test is not validated for emission-law reach (157's pair was a spread-policy change); rule 13's
toll = mean(a21, a22) exactly. hy_kt8000 corrected to −0.009…−0.004: DO NOT SUBMIT. Knee alpha GO
for D/E/F (board-validated ±1%). Long-gene axis confirmed by gene length (ρ −0.71). corr fixed to
+0.701. n_real restated as a distribution. New rule in AGENTS.md: retrodict the family's board point.

### 2026-09-03 03:02 UTC+3 — 182 / 183 launched

182: pretrained Arc State Transition checkpoints as a new sign source (106/113 only tested a
tiny local `state_sm`). 183: norm-preserving unique-component boost (raise X not A; audit 181
reopen). Both Claude Code, no submit.

### 2026-09-03 03:00 UTC+3 — abc-knee SUBMITTED and published (Sep-2 slot)

Ahmet lifted the hold ("ur the lead scientist… do use it"). Submitted
`submissions/abc-knee.vcc` at 23:28 UTC as `We regret to inform the genes`, description empty.
Entry `N3ui4hCkKx9bxExG2pxE`. Score `0.1740437` (rank 24/614) vs a23 `0.173360` (+0.00068) vs
a15 `0.174087` (−0.000043).

**The slot's question was answered: per-context knee does not transfer.** Δ`mse` scaled −0.00226
against a prereg of +0.0050 — falsifier fired. `pds` invariant. D/E/F default = global alpha,
not `VCC_ALPHA=knee`. `hy_kt8000` not submitted. No Sep-3 slot spent.
Write-up: `experiments/180-def-readiness/BOARD-RESULT.md`.

### 2026-09-03 02:15 UTC+3 — nmae landscape formula recorded

`nmae(k,q) = q|k−1| + (1−q)(k+1)`; optimum k = 1 worth 2(1−q); slopes −0.074 below / +1 above at
our q = 0.537; our board 0.944 = k ≈ 0.75. Explains why sharpening/a13/a21 died on nmae (overshoot).
Design constraint for the reopened pds lever: never inflate DE-pool magnitudes past k = 1.

### 🔴 CORRECTION to the entry above (187, 2026-09-03)

**"our board 0.944 = k ≈ 0.75" is an artefact and must not be used.** The uniform-`k` identity does
not describe the scored quantity: on h1-test the **|t|-weighted per-gene `k` has median 0.12 and p99
1.11**, and the board's own amplitude response is 75% of the uniform-`k` slope. 187's verdict:
**"`nmae` is a direction problem, not an amplitude problem"** — within the `mse`-feasible amplitude
window the collectable `nmae` gain is 0.0012 raw against a 0.0163 raw gap to the leader. The design
constraint the entry draws from it (never inflate DE-pool magnitudes past `k = 1`) happens to remain
correct, for a different reason: 188's exact identity has slope +1 above `k = 1`.

### 2026-09-04 06:55 UTC+3 — three slots scored in two days, 0.1740865 → 0.2049216, rank 5 of 656

`rc_hybrid` **+0.0093325** (emitted p-order at bit-identical pseudobulk) → `rc_callset` **+0.0084448**
(the call set, `n_pred` 4,480 → 2,599) → **`cd4f034` +0.0130557**, the **first arm ever to move the
direction channel** (`VCC_CD4_TARGET_F` 0.10 → 0.34; implied cosine 0.07228 → **0.09018**). Only five
teams are above 0.20 and one above 0.22. **World rank per member: `jac` 1st, `mse` 10th, `reach`
10th, `pds` 7th, `nmae` 59th, `fid` 129th — `cd4f034` improved every one.**

**The gap to 1st is +0.01939 and 82% of it is cosine** (197's inversion: ours 0.09018, `cqawesome`
0.13356, 1.481×). We are behind on all three direction members and **ahead on all three ordering
members**. `AGENTS.md` **§4.1c** records this and retires §4.1/§4.1b, whose claims ("nobody is
shrinking", "`mse` stays 0 for us") described an August board — the top 12 now sit at raw `mse`
0.730–0.856 and we score +0.1593 on it.

**MT Zeng holds the board's best direction (cosine 0.14138, above `cqawesome`'s) and ranks 8th**,
because their ordering members throw it away. A point of cosine is worth more to us than to anyone
on the board. **Do not trade ordering for direction.**

### 2026-09-04 06:52 UTC+3 — `h1union` SUBMITTED (second Sep-4 slot)

Entry `iUziqE1gewfcPKzYUDaB`, `The Nucleus Declined To Comment`, description empty, md5-verified at
2,824,785,920 bytes. 222 returned **BUILDABLE** — 3 known-answer rows, 5 falsifiers, default path
byte-identical to `cd4f034`, CD4's solved fraction unmoved, so the arm is single-variable.

**Band reissued before the score existed** (`DECISION-TREE-AMENDED.md`): the original +0.0055–0.0105
came from the 4-fold Δ`pds`, inflated by `kolf2` — an iPSC fold gaining 2.3× the next because both it
and the added H1 source are pluripotent, while A/B/C are all cancer lines. Honest 3-cancer rate gives
**+0.0010 to +0.0045, central +0.0030**, after discounting for 222's concentration finding
(Δ`pds` retains **14.6%** under leave-top-3-out; the cosine gain retains **72.7%**).

### 2026-09-04 — six routes closed, two hypotheses left

**216** killed the "better data" reading at byte level: 12 byte ranges match, `LastModified`
2026-05-28, **zero objects moved on that bucket during the leaders' doubling window**; it also put
the CD4 dose turnover at `f* = 1.018` and confirmed `f = 0.34` is the argmax over a 31-point grid
(the cosine argmax costs **−0.01195** overall and reaches only 0.799× the leaders). **218** killed
independence. **220** bounded the entire ordering channel: being the world's best `reach` team is
worth **+0.0077**, less than one direction arm. **223** resolved the sizing bug as a **measurement**
bug, not a shipping bug — **0 of 21 experiments void** — and found **DLD-1 live at +0.0025 marginal
over the union**. **187** closed `nmae` as an amplitude lever. And mining the `description` field of
all **2,306** scored entries closed the public-text route: 764 are non-empty and **not one above
implied cosine 0.09 has one** — the whole top of the direction channel keeps it empty, as we do.

**Two readings of 215's many-small-steps shape survive, both running:** **224** (iterative retraining
of a learned transfer model; gate +2.9% relative cosine on held **lines**, slope in `k` as the
decision variable) and **225** (public sources added one at a time; exhaustive census, saturation as
the decision variable).

Infrastructure: `scripts/supervisor_loop.sh` daemonises agent supervision (it was Monitor-driven and
therefore not actually guarding anything); 221's `scripts/aws_reap.sh` sweeps 17 regions on a
two-reading rule after finding that a percentage-only threshold would reap a busy single-threaded
scorer; 24 GB reclaimed from sha-verified duplicate packages.

### 2026-09-04 14:45 UTC+3 — LEAD CHANGE. The coordinator session was stopped; a fresh Fable 5.1 session took over

Ahmet stopped the previous coordinator (it had silently switched to Opus 5 at 06:36 local; its own
audit, 239, found the summary layer drifting while the archive stayed reliable). Handover state:
board `h1union` 0.2118310 rank 3/656, both 09-04 slots spent, 241 (ablation re-price) running, 242
(guard audit) finished, cloud bundle instance `i-08028907eeba5a980` building the h1-test real bundle
under a 15 h cap. First actions: committed the finished-but-untracked overnight archive
(`c44d845`); armed the preregistered 237 slot — `submissions/229-stack6-g035.vcc`, sha
`7a3c2d3e…`, fires at 2026-09-05 00:00 UTC via `scripts/arm_submit.sh` (kill
`experiments/237-gate-validation/submit-armed.log.pid` to veto); started `scripts/lead_monitor.sh`;
launched 243 (reconcile the PLOLO→board cosine factor, 1.384 vs 1.000), 245 (a cited READ-FIRST
built from the archive only), 246 (audit + measure the best OPEN direction hypothesis on top of the
shipped stack); wrote 244 (full Oct-22 rehearsal on the fixed driver, essential-free pessimistic
panel) to launch when 241 frees memory.

### 2026-09-04 15:10 UTC+3 — 245 landed: `docs/READ-FIRST.md`, and the summary layer is now subordinate to it

A 12 KB map built only from primary sources with a citation on every number; every subagent's system
prompt now says to read it before `AGENTS.md`, and a banner on `AGENTS.md`/`CURRENT_STATE.md` says
READ-FIRST wins where §8.x disagrees. Its contradiction list (15 rows, `docs/READ-FIRST-contradictions-
2026-09-04.md`) includes: §8.4's 1.80× and 2.0× pricing constants vs the two measured board points; §8.2's
"exactly two levers" both closed; the propagated `cd4f034` typo 0.2049216 (fixed to 0.2049196 in four
files); the h1-test→board `reach` factor 1.290 vs 0.769 on the real scorer; and that `vcc status` now
returns `not_found` for every entry but the current one, so the `BOARD-RESULT.md` files are the only
record of earlier scores. Supervisor relaunches now default to Fable (were Opus).

### 2026-09-04 15:40 UTC+3 — 243: the PLOLO→board cosine factor is a class constant, and the floor mis-reads amplitude

0 of 48 conventions admit one scalar within ±15% of the three board points; amplitude is excluded by
streaming each package's own `A_c`/`B_c` (measured cosines 0.0792 → 0.0979 → 0.1073). Two classes:
re-dose of occupied rows k 1.36 (a14→a15 1.315, `cd4f034` 1.405); first source into empty rows
0.94–1.00 (`h1union`). ARM6 re-priced on direction legs at +0.0022…+0.0038 with its −1.24%
amplitude charged; 229's −0.00219 differs only by the `nmae` identity leg, which is exactly what
tonight's slot tests. Reading table pre-committed in `experiments/237-gate-validation/ADDENDUM-243.md`.

### 2026-09-04 15:55 UTC+3 — 241: the October component ranking survives the renderer change; PerturbMatch and the CD4 row re-decided

200's ten rows reproduce exactly on its own renderer, and on the shipped renderer no component moves more
than one rank or more than its error bar, so the ordering is safe to quote again. PerturbMatch is
−0.002178 at the shipped dose, outside the gate but sub-gate in both channels (direction 0.74× its gate,
ordering 0.17 sd): ship it, stop citing −0.001615. The CD4 row was 2.1–2.5× understated (structural zero
on h1-test). `DEF-RUNBOOK.md` §3.1 updated.

### 2026-09-04 16:30 UTC+3 — 246: the first arm that BUYS `pds` — deflate the CD4 residual's target-shared top-5 directions

Five of six named direction hypotheses are closed or worthless; the open one (project the target-shared
subspace out of an un-whitened residual) pays nothing on the self-residual but, applied to the shipped
CD4 residual (which enters un-whitened, so nothing ever strips its T-cell-programme hub), buys raw `pds`
+0.01227 on 3/3 cancer folds (8.2σ) for −1.2% cosine. Energy-matched random subspaces buy no `pds`;
shuffled projectors ≤ +0.0013; split-half projector keeps +0.0104–0.0118; leave-top-3-out retains 56%.
Six-member price +0.0023 / +0.0037 / +0.0082 at `pds` transfer 0.7 / 1.0 / 2.0×; shape-only class.
BUILD-WORTHY conditional on the board's `pds` transfer. Candidate for the second 2026-09-05 slot.

### 2026-09-04 16:00 UTC+3 — 247: the Oct-22 path is guarded, each guard observed firing; CD4's dose is panel-invariant off A/B/C

`def_postbuild.py` now asserts the operator by registry pin, file hash and symbol-axis digest (it passed
20/20 on a swapped operator before; 31/31 on the shipped meta, FAIL on the swap now); row counts compared;
the census is fail-closed; the driver catches all six launch forms and refuses a stale operator on a
new panel; `VCC_CD4_SIZING=carrying` + `VCC_CD4_TARGET_G=0.3443655278649641` (0.34 × A/B/C's
carrying/all-rows ratio 1.01284) is wired for non-A/B/C bundles, delivering ×1.000 of A/B/C's per-row
dose where the legacy rule delivers ×1.18/×1.16 on the two pessimistic panels. Flag off is byte-inert:
the shipped `.h5ad` `8df4dadd…` rebuilt through the changed driver and builder. Novelty ledger tracked,
backfilled, and `arm_submit.sh` now refuses a scored package. Launched 249 (build the CD4 deflation
arm) and 244 (the end-to-end rehearsal) on the released builder.

### 2026-09-04 17:15 UTC+3 — 244: the Oct-22 sequence ran end to end on the fixed driver; executable only with the coverage branch taken explicitly

All eleven stages on 238's disjoint, essential-free bundle (K562 146/300, CD4 union 186/300) with nothing
exported: operator re-fit, CD4 re-extracted (182/300, over-dose factor 1.1702 bit-exact with 247's),
H1 union sized on 28 carrying rows (all-rows rule 0.0), `def_postbuild` 32/32, `vcc prep --dry-run
--contexts D,E,F` clean, 29.7 min launch-to-artifact. But the runbook's literal command stops at step 3
(047's CD4 gate, 186 < 290) and its §5 response cites a file only the relaxed extraction writes; the
prep stage also refused once under the 30 GiB gate. Eight runbook lines to change (244 §8). The skip-
branch hole 240 named is closed by 247 and was observed closed here.

### 2026-09-04 17:45 UTC+3 — 249: `cd4svd5` BUILT and armed for the second 2026-09-05 slot

`submissions/249-cd4svd5-f034.vcc` sha `6e8b41c6…` (2,824,683,520 B). Flag-off byte-inert (staging
`abeb0749…`, package `8df4dadd…` rebuilt through the changed builder); the staged arm reproduces 246's
resolved row to 0.0e+00 on 16/16 fold×arm rows; energy-matched random control through the same path
buys no `pds`; rule-10 diff exactly `cd4_deflate_k` plus the 11 pre-declared CD4-frame scalars;
non-carried rows equal the incumbent times one scalar to 1e−16; `n_pred` +0.00%; amplitude −0.41%
(predicted board raw `mse` 0.8191); prep clean. Builder commit `be029c8`. Armed via
`scripts/arm_second_slot.sh` (pid in `experiments/249-cd4-deflate/submit-armed.log.pid`): fires as
"You can do it too.", description empty, after the 237 entry is terminal and after 00:30Z.
Launched 248 (second union residual in the builder; DLD-1 rider) on the released builder.

### 2026-09-04 18:25 UTC+3 — 248: the second union residual is in the builder (inert); the DLD-1 rider does NOT clear the gate

`VCC_UNION2_SOURCES`/`VCC_UNION2_G` (commit `2632b92`) express 225/229's CUM construction; flags unset
rebuild the shipped package byte-for-byte (`8df4dadd…`, `def_postbuild` 32/32). Through the builder the
rider `dld1_lowmoi@0.35` reproduces 225's CUM1 row to 0.0e+00 (120/120 rows), but prices +0.00216 under
243's class rule (1.14× the shape-only gate, below the 1.5× required) with 82% of its `pds` gain in three
targets (leave-top-3-out keeps 25%); the exact `nmae` identity puts it at +0.00107. No package built.
The DLD-1 rider is closed as a slot arm; the capability stays for October source stacks.

### 2026-09-04 18:50 UTC+3 — 250: the deflation transfers where CD4 covers the panel; October must PIN the A/B/C basis, never re-fit it

On four disjoint 300-target panels (233's `real`/`pess`, 238's `efp`/`efr`) the day-of re-fit at k = 5
keeps the `pds` mechanism (+0.0148 / +0.0045 / +0.0035 / +0.0099) but costs 4–14× A/B/C's cosine, so it
prices ≤ 0 on both low-coverage draws (CD4 174 / 186 of 300). The PINNED A/B/C basis (`U_5` sha
`f06e0262…`, the object 249 built) beats the re-fit on 4/4: `pds` +0.0153 / +0.0037 / +0.0073 / +0.0107
at −0.1 … −2.0% cosine, price +0.0061 / +0.0009 / +0.0014 / +0.0030, never negative; any foreign basis
does the same — a basis fitted on other rows removes the shared hub and nothing else. k = 5 sits below
the `pds` plateau off A/B/C (knee 10–15) but the price is flat past 5–10. Census trigger: CD4 carried
rows < 250/300 → the re-fit variant reads ≤ 0. Launched 251 to put the pinned-basis path in the builder.

### 2026-09-04 19:55 UTC+3 — 251: the pinned deflation basis is in the builder (inert); October can ship A/B/C's V_5 on any panel

`VCC_CD4_DEFLATE_BASIS=V5_abc` (commit `c321f76`; artifact `data/donors/cd4_deflation_basis/V5_abc.npz`
sha `ba7ade3f…`, U_5 `f06e0262…` = the shipped `cd4svd5` audit). Flags unset: staging `abeb0749…` and the
package `8df4dadd…` reproduce; the re-fit flag still gives 249's `a28d5be2…`; on A/B/C the pinned path
equals the re-fit to 2.4e−15 (numeric, not byte: two associations of one projector); 250's pinned rows
on `real`/`efr` reproduce at 0.0e+00 (24/24); the driver rehearsal on the ef-bundle passes 40/40 postbuild
and prep dry-run with the captured-energy diagnostic 0.1334. `def_postbuild` accepts the new keys.

### 2026-09-04 20:15 UTC+3 — 253: READ-FIRST v2 installed; the remaining summary-layer contradictions annotated

`docs/READ-FIRST.md` rebuilt from primary sources after today's twelve results (11,963 bytes, 176 citations
checked). Of 245's fifteen contradiction rows, eight were resolved by today's edits; the rest are now either
corrected in the primary (`cd4f034`'s "pending" ledger row → +0.0130557) or annotated (221's 0.769× `reach`
factor next to the 1.290× claim; a superseded-in-part banner over `AGENTS.md` §8.1–§8.2).

### 2026-09-04 20:50 UTC+3 — 252: the shared PLOLO harness is guarded; every published number reproduces at 0.0e+00

`plolo.py` (`3decf80`) and `064/run.py` (`e2f8d71`): operators are loaded by registered name with sha and
axis checks (the three near-identical files are told apart; a one-ULP tamper raises); the fold cache is
keyed on its inputs; an all-zero truth block reads NaN instead of 0.0 sign accuracy; the on-axis preserve
raises at a zero median; and the silent `f_per_unit = 0.0` below 50% coverage now raises by default with
`sizing="carrying"` and `sizing="legacy_all_rows"` reachable by name. 14 guards shown firing against the
pre-252 byte copy; 97 known-answer rows (216, 218, 222, 233, 246, 249) at 0.0e+00; the k562 fold rebuilt
from scratch bit-identical. One existing script (`trapprice233.py`) needs `check_axis=False` to run.

### 2026-09-04 20:52 UTC+3 — joke names swapped on the two armed 09-05 slots; packages and times unchanged

Killed pids 26557 / 80453 (neither had submitted). Re-armed:
- 00:00Z `229-stack6-g035.vcc` as `I'll just keep moving forward` (pid 73835)
- 00:30Z `249-cd4svd5-f034.vcc` as `You can do it too.` (pid 73837), still waits for the first entry to be terminal
Description empty on both (§14).

### 2026-09-04 21:05 UTC+3 — first-slot name finalised as "I'll just keep moving forward." (with the period), pid 73250; second slot unchanged ("You can do it too.", pid 70577)

### 2026-09-04 21:20 UTC+3 — 254: `cd4svd5+partial` BUILT (`eea5ed18…`), inert flag committed; the partial rows alone are under the gate

`VCC_CD4_PARTIAL_G` (commit `32c5219`): flags unset rebuild `8df4dadd…`; 246's `C2_partial@0.34` and
`E_cd4_svd5+partial@0.34` reproduce through the builder at 0.0e+00 (104 rows); the package's increment over
`cd4svd5` is the partial row itself, +0.0009 … +0.0018 — under the +0.0019 gate at ≤ 1.0× `pds` transfer, so
it is NOT a slot candidate on its own (its `PREREGISTRATION-DRAFT.md` waits for the realised factor).
Branch tree for `cd4svd5`'s result pre-committed in `experiments/249-cd4-deflate/DECISION-TREE.md`.
Nothing heavy is justified overnight beyond the two armed slots; the next launches depend on their readings.

### 2026-09-04 23:05 UTC+3 — 255: `cd4svd10` BUILT (`f860632a…`) for branch A of the cd4svd5 tree; not a cleared arm on price

k = 10 on top of `cd4svd5`: 246's `B2_cd4_svd10_resolved` rows reproduce to the digit; every gate of 249's set
passes against `cd4svd5`; amplitude −0.46%; prep clean. PLOLO prices the increment over `cd4svd5` at
−0.0005 … −0.0002 (`pds` +0.00054 on 1/3 folds) — the package exists only so a branch-A board point can buy
the k-response October needs (250: plateau at k 10–15 off A/B/C). Draft preregistration in its directory.

### 2026-09-05 00:15 UTC+3 — vcc CLI is 0.2.0 (found upgraded at 22:35 local); verified before tonight's slots

Arc's release mail: cancel command, prep memory halved, no change to submission validity. Verified under
0.2.0: `--json whoami` still carries `"can_submit": true`, `--json status` `"is_terminal": true`; `submit`
keeps `-m`, `-d`, `--wait` (no default timeout; new `--wait-timeout`, `--skip-limit-check` = a daily-allowance
pre-check now exists); `prep --contexts` present; `vcc cancel` exists. Known-answer rows: `prep --dry-run` on
`249-cd4svd5-f034.h5ad` clean at 15.9 GB peak (was 22–27), and on 244's D/E/F rehearsal artifact with
`--contexts D,E,F` clean at 16.0 GB. `vcc skill install` re-run. Because of the new allowance pre-check the
first slot now fires at 00:02Z instead of 00:00Z, and `arm_submit.sh` retries once after 90 s.

### 2026-09-04 23:45 UTC+3 — 221's h1-test real-bundle job hit its 15 h cap and FAILED (S3 `STATUS`: `FAILED 1 (completed=0)`); not relaunched

The `prep-real-bundle` on h1-test did not finish inside the 900-minute guard; the instance self-terminated
at 20:27Z as designed (~$8 of credits). Not relaunched: since `h1union` shipped, every arm on top of the
incumbent carries `h1_2025`, which covers 85/85 h1-test targets, so a scored h1-test bundle could only
calibrate pre-`h1union` artifacts. Revisit only if a non-leaked local six-member instrument is needed.

### 2026-09-05 02:35 UTC+3 — New coordinator session took over (the 1M-context session was closed); tonight's two slots audited, evidence capture and result reading fixed BEFORE they fire

Both waiters survived the session close (pids 41699 first slot 00:02Z, 77532 second; re-parented to launchd).
Packages, hashes, names and times unchanged. Verified against the installed `vcc` 0.2.0 source and live calls:
`submit --wait` prints only the SCALED six (`vcc/_scores.py`), `--json status` carries the raw members, and every
non-current entry id returns `not_found` (re-confirmed on `QSe3JlNUYZ2BzvlcsL14`) — so the 237 entry's raw `mse`,
the number 243's class reading needs, would have vanished when the 249 entry published. Fixes, none touching the two
RUNNING scripts (bash reads scripts incrementally; `arm_submit_after.sh` is exec'd later, so it was safe to rewrite):
`scripts/capture_board_entry.py` (validated, atomic capture to `data/board-entries/<id>.json`; h1union captured at
23:23Z); `scripts/arm_submit_after.sh` captures the prior entry FIRST and fails closed without it, and captures its
own entry after scoring; `scripts/read_board_entry.py` reads captures or ids, applies only the named prereg
(`--prereg 237|249`), reads the `nmae` leg DIRECTLY off `score_nmae` against −0.004783 instead of inferring its sign
from the overall bin (a mocked overall +0.0020 / nmae −0.0048 entry was previously declared a mis-sign), and applies
the 249 tree with precedence E>D>A>B>C plus a HOLD branch. Pre-result addenda: `E237/ADDENDUM-2-READING.md`,
`E249/DECISION-TREE.md` §Addendum (k = 10 not automatic on A; 255 priced its increment negative). `caffeinate` for
4 h; supervisor loop has no candidates. Decision on the candidates: keep both — k5 is the stronger arm (3/3 folds,
controls, split-half); ARM6 is an information purchase on the `nmae` gate, its readings now member-resolved.

### 2026-09-05 03:35 UTC+3 — 237 board point: `ARM6` 0.2081907 (rank 4), Δ −0.0036403; the `nmae` identity is VALIDATED at scale; both records captured before the API dropped them

Fired 00:02:21Z, published 00:28:02Z (entry `ydAemVdWmsTx9xhDWcou`). Direct `nmae` read: −0.00716 overall vs the
identity's −0.004783 — same sign, 150 % of the magnitude, so the P237 §4 "mis-spent" falsifier (within 20 %) does
not fire: the identity is a conservative lower bound on the `nmae` loss for a six-source union, and 228's
per-gene gain and ARM6 stay dead on that leg. Direction legs excluding `nmae` +0.00352, inside 243's band; raw
`mse` 0.81227 vs 243's point 0.81244, realised k_cos 0.962 (new-source bin) — the class model's fourth point.
`pds` +0.00248 overall (raw +0.00666, 1.9× `h1union`'s own step) and `mse` +0.00088 are the union's real gain;
`nmae` is the whole loss. Open hypothesis, not a candidate: an `nmae`-neutral rendering of the same directions,
priced on the validated identity. Evidence: `data/board-entries/ydAemVdWmsTx9xhDWcou.json` (both capture paths
succeeded, no race), `E237/BOARD-RESULT.md`, `reading-ydAemVdWmsTx9xhDWcou.txt`, `SUBMITTED_LEDGER.json`,
`BOARD-LEDGER.md`. The second slot (`cd4svd5`) fired 00:30:21Z after the capture was on disk.

### 2026-09-05 04:20 UTC+3 — 249 board point: `cd4svd5` 0.2160944 (rank 3 of 688), Δ +0.0042634 — the new incumbent; mechanism transfers at 1.2×, cosine cost 3.8× PLOLO's read; October ships the pinned basis

Fired 00:30:21Z after the 237 record was on disk (the new capture gate), published 00:54:57Z, entry
`Xv0ayvuoHpsKYWv2iu42`, record `data/board-entries/Xv0ayvuoHpsKYWv2iu42.json`. Raw `pds` +0.01492 (k_pds 1.22; leg
+0.00545 overall, the largest `pds` step we have bought), raw `mse` +0.00454 worse (falsifier 2 fires: board cosine
−1.735 % vs PLOLO agg −0.4617 % → re-shape class k_cos 3.76), the four DE members flat as priced. Branch D by the
letter of the tree's precedence; A on the evidence (falsifier 1 not fired by 3.7×, both A conditions hold, and
RB §3.1's independent 250 trigger — "if board raw `pds` rises, ship the pinned basis" — fired). BR249 §3 records both
and decides: October ships `VCC_CD4_DEFLATE_BASIS=V5_abc`; k = 10 is dead on the realised factor; `+partial` needs its
own transfer reading; every re-shaping arm re-priced on 3.76×/1.22×. Gap to first +0.00845 (was +0.01272), to
second (Aginglab 0.2193913) +0.00330. READ-FIRST v2.2, CURRENT_STATE header, BOARD-LEDGER, RUNBOOK §3.1 updated.
Pending: the driver wiring for the pinned basis (R251 §7) with 240's row 1 re-run against `6e8b41c6…`.

### 2026-09-05 04:55 UTC+3 — 256: the CD4 deflation wired into the October driver by panel; 240's row 1 re-run byte-identical; `vcc` 0.2.0 changes the `.vcc` digest of a byte-identical h5ad, novelty guard closed

`scripts/def_pipeline.py` now defaults to `VCC_CD4_DEFLATE_K=5` on A/B/C (`ABC_ENV`) and `VCC_CD4_DEFLATE_BASIS=V5_abc`
elsewhere (`DEF_ENV`), asserts both off the artifact, and diffs against the 249 meta; `scripts/def_postbuild.py` names
`cd4svd5` as the incumbent and requires the deflation by panel. KA row 1: the driver's defaults, nothing exported,
rebuilt `aba27758…` on 12,648,807,800 B — MATCH (1383 s, 18.30 GiB), postbuild 34/34; the h1union meta fails on
exactly the new row (guard demonstrated). The re-prepped `.vcc` is `7f8336fb…`, not the scored `6e8b41c6…`: every
value inside equal, the 0.2.0 CLI writes `indptr` int64 and adds `meta.json`. So `check_submission_novelty.py` now
also refuses on the prediction digest (sibling meta's `output.sha256`), ledger backfilled 11/12, four demonstrations
in R256 §5. Remaining docs edit: RB:202-220's deflation row as SHIPS.

### 2026-09-05 05:10 UTC+3 — 257: a softer CD4 deflation DOES NOT CLEAR; hypothesis (b) of BR249 closed on A/B/C

Measurement-only subagent (25 min, 91 turns). Retrodictions first: the recipe on the realised factors (k_pds 1.22,
k_cos 3.76) reads `cd4svd5` vs `h1union` at +0.00472 / `pds` +0.00552 / `mse` −0.00083 (board +0.00426 / +0.00545 /
−0.00080); `k5 s=1` reproduces 246's row at 0.0e+00; `k10` prices −0.00212 (255's increment). Then the 4×3 s×k grid
and five single directions: every cell loses `pds` against the incumbent; the whole −0.46 % cosine arrives in the last
quarter of the removal, so the best cell (k=5, s=0.75) is +0.00055 conservative / +0.00023 point, 0.29× the gate, and
negative at k_cos 2.15 or 1.36. Structure: direction 1 (8 % of the centred energy) carries 75 % of the `pds`;
direction 4 costs +0.38 % cosine on 3/3 folds and buys none — the projector without it is 0.39× the gate (post hoc).
Random energy-matched controls buy no `pds` at any s. READ-FIRST v2.4. 258 (union dose on cd4svd5) still running.

### 2026-09-05 05:25 UTC+3 — 258: the six-source union DOES NOT CLEAR on `cd4svd5`; the family is closed on this base. Nothing is priced for 09-06 — the slots are held

Measurement-only subagent (36 min, 107 turns). R1 reproduces 229's `six229.json` at 0.0e+00; R3 retrodicts BR249
12/12; R2 (ARM6 by member) lands only with k_pds at the ARM6-realised value (2.48 A/B/C mean, 3.79 frame A), the
identity ×1.5 and the amplitude leg — and the DE-three rate is 0.0425 per unit Δcos, not 243's 0.24. Screen: 31 cells
(g 0.10–0.45 × subsets), 0/31 clear at conservative factors; best H1 pair + `dld1_lowmoi` g=0.35 at +0.00119 (0.63×
gate, `pds` 0.8 σ, leave-top-3-out −66 %); the full union −0.0027 … −0.0051 at every dose (the identity's leg is
2.5–3× its gain; `igvf_cardiomyocyte` carries the `nmae` cost). The deflation already bought a third of the union's
`pds`. An `nmae`-neutral rendering prices +0.0014 conservative / +0.0034 point on this base — under the gate.
Together with 257: no candidate clears for 09-06; per the standing rule the slots stay unspent rather than bought on
a guess. READ-FIRST v2.5. The defaults-only D/E/F rehearsal (256) started behind 258.

### 2026-09-05 05:45 UTC+3 — 256 KA 6: the D/E/F rehearsal with NOTHING exported passes end to end (postbuild 41/41; the meta carries the pinned basis, captured fraction 0.1334)

251's recipe minus its one export: the driver's own defaults built the ef-bundle with `VCC_CD4_DEFLATE_BASIS=V5_abc`
and no re-fit key (`[arms]` OK on all eight rows), postbuild 41/41, prep dry-run clean, builder restored to
`3615d20f…`. 22.4 min build at 19.42 GiB, prep 10.96 GiB. The October path carries the deflation by default.
Overnight ledger: 256 wired + rehearsed; 257 and 258 closed their families; nothing priced for 09-06.

### 2026-09-05 06:15 UTC+3 — 259: the union without igvf does not clear; the IGVF attribution is confirmed directly; the second-residual construction is worth ×1.7–1.8 the single union

Measurement-only subagent (25 min). R1 reproduces 258 at 0.0e+00 (4 × 338 rows); R2 retrodicts ARM6 by member
(−0.00367 vs −0.00364) and cd4svd5 (+0.00470 vs +0.00426); R3 reproduces 248's second residual at 0.0e+00 through the
union2 path. 21 cells + 9 controls, three coherent scenarios: 1/21 clears the pessimistic gate (B3 g₂=0.45, +0.00211,
outside the grid) and 0/21 clear all three conditions — the `pds` that remains after removing IGVF is DLD-1's three
K562 targets plus one Jiang target on hct116, negative on hek293t, leave-top-3-out −45 … −292 %. IGVF back: Δ`nmae`
+0.0031/+0.0084/+0.0142 raw for +0.0004 `pds` (attribution confirmed). DLD-1 alone as a second residual on cd4svd5:
`pds` +0.00217 (3/3, 1.5 σ), pessimistic +0.00094, retention 15 %. READ-FIRST v2.8. 260 (feasibility) reports GO for a
bounded instrument calibration (raw K562-GWPS file already on disk, native depth 11.3 k UMI, ≤ 250 cells/target,
K562-free fold model); not launched.

### 2026-09-05 06:40 UTC+3 — 261: deflating the H1-union residual does NOT transfer the CD4 mechanism — the union has no hub. Fourth clean negative of the night; nothing priced for 09-06

Measurement-only subagent (21 min, 67 turns). R1 (k = 0 == cd4svd5, 44/44 keys per fold), R2 (BR249 retrodicted
+0.00472), R3 (the same code on CD4 reproduces 246 at 0.0e+00) all pass. The union's centred spectrum is flat
(5.2 % leading vs CD4's 8.0 %); each projector costs −10.9 … −13.1 % aggregate cosine per unit energy (CD4: −2.3 %)
and buys +0.0011 … +0.0015 `pds` (0.7–1.0 σ) — pessimistic −0.0017 … −0.0081, and even the optimistic scenario is
≤ +0.00004. Split-half and shuffled projectors buy as much as the in-sample one: the gain is generic, not a
property of the source. Re-dosing at g = 0.35 nets −0.0008 … −0.0001 after the identity's `nmae` charge. Closed at
every k and dose tested. READ-FIRST v2.9. 260 (feasibility) still finishing.

### 2026-09-05 06:45 UTC+3 — 260: a K562-GWPS raw-count DE diagnostic is NO-GO on coverage, measured; the night's ledger

Feasibility subagent (50 min, 94 turns; one step peaked at 8.25 GiB against a ~6 GB ask, reported). The file is on
disk, streams in 1.7 min, the panel/control depth ratio is 1.0035 (PASS without repair), the K562-free model is
PLOLO's fold, the whole diagnostic would be ≈ 1 agent-day and ≈ $2 — and it still cannot bound what the top-89 proxy
misses: the scorer's own Wilcoxon recipe on the actual K562 cells gives the panel targets a median of 5 significant
genes; 31 % clear `nmae`'s 10-gene gate; 400 cells and depth do not change it. Ledger of the night (2026-09-05
00:02–03:45Z): two board points read and captured (ARM6 −0.00364, `cd4svd5` +0.00426 → 0.2160944, rank 3), the
`nmae` identity validated at scale, the October driver wired and rehearsed with the deflation by default, the novelty
guard keyed on the prediction, five measurement experiments (257, 258, 259, 261 negatives; 260 NO-GO), every one
retrodicting the board first. Nothing priced for 09-06; slots held. READ-FIRST v3.0.

### 2026-09-05 07:15 UTC+3 — 262: where `cd4svd5` loses `pds` — a fold-consistent set of fully-supported targets whose donor consensus points away from the held line; no hub left; the marginal band is the lever

Diagnostic subagent (30 min, 74 turns; KA rows at 0.0e+00 incl. the per-target arrays). Per-target `pds` is
bimodal; the bottom quartile carries 71–79 % of the deficit and the bottom deciles overlap 9–12 of ~27 across the
three cancer folds (Spearman +0.70…+0.76). One covariate explains it: base·truth cosine (ρ +0.875, AUC 0.05);
`pds` is a threshold in it (+0.095 per +0.01 inside (0, 0.05], ≈ 1 above 0.03, below chance below 0). Partition of the
deficit: base points away 49–58 % (17.5 % of targets, cross-line agreement ≈ 0), the residual breaks a positive base
12–18 %, hub 3–9 %, small truth 3–9 %, render 9–12 %. No hub remains after the deflation (hubness ρ ≤ 0.12).
Observable selectors cannot find the below-chance group. Pilots: agreement-keyed CD4 dose +0.00045 (0.3 σ) at
−3.46 % cosine (dead); per-target donor choice has oracle headroom but CD4 as arbiter picks wrong. Three
hypotheses, none priced: H1 dead by pilot; H2 (donor selection with an independent arbiter) → 263; H3 is a reading
rule, now in READ-FIRST §3. READ-FIRST v3.1.

### 2026-09-05 08:00 UTC+3 — 263: H2 closed for every observable arbiter; the oracle headroom is +0.046 `pds`; one production-realisable signal on the only three-donor fold — the majority rule

Subagent (43 min, 78 turns; KA-0/R1/R3/KA-4 at 0.0e+00). Per-target ORACLE selection between the two whitened donor
rows: +0.0459 `pds` (3 cancer folds; +0.049 / +0.047 / +0.042), 3.7× the deflation, 65 % on 262's anti-aligned
group; the un-whitened oracle captures only 44 %. 78 observable selection cells: every 3-fold mean negative, agreement
with the oracle 44–61 % (coin flips), random picks cost −0.003 … −0.024, reverses within thousandths of their arbiters;
0/130 cells clear. On `kolf2` (three production donors) the MAJORITY rule (keep the two agreeing raw rows) reads
+0.0121 `pds`, +2.4 % cos_agg, shape-only, retention 60 %, reverse −0.0097, random −0.0049, not reproducible by any
global donor weight; priced +0.0051 pessimistic — one fold, so not clearable under the 3/3 rule. READ-FIRST v3.2.
Next: 264 puts `kolf2` in as the third donor on the cancer folds to see whether the mechanism holds on cancer held lines.

### 2026-09-05 08:35 UTC+3 — 264: the majority rule does not hold on the cancer folds; H2 closed entirely; the kolf2 signal retracted by the correct control. No submission on 09-06

Subagent (30 min, 71 turns; KA-0, R1, R2, R3 at 0.0e+00 — the prepared kolf2 donor rows are bit-identical to the
truth rows 263's fold consumed). Letter construction (kolf2 at weight 1 through the two-donor operator, which
whitens kolf2 to 46–90× its raw norm): the rule −0.0195 `pds` vs the incumbent, 0/3. Amplitude-matched: −0.0107,
0/3; vs the three-donor mean −0.0051; every gain is dropping kolf2 (junk removal), every drop of a production donor
loses. Adding kolf2 as a donor helps on no fold. On kolf2 itself the +0.0121 stands but dropping one RANDOM donor
reads +0.0068 (56 % of the arm; 263 had used keep-one-random as the control). The oracle drop-one ceiling is +0.043 …
+0.048; the majority rule captures −11 … −25 % of it. H2 closed for the majority rule too. Seven measurements since
the board points (257, 258, 259, 261, 262, 263, 264), all negative or diagnostic, every one retrodicting or
reproducing at 0.0e+00 first. Decision: the 09-06 slots are not spent. READ-FIRST v3.3.

### 2026-09-05 16:05 UTC+3 — Second external review acted on: gap corrected (+0.01182, leader `m0905b` 0.22791909 at 01:05Z), sigma semantics fixed, 264's verdict gate fixed, H2 closure scoped, the guard made content-bound; 265/266 launched

The "σ" quoted overnight divided per-fold Δ`pds` by 162's 0.001494, the leave-one-arm-out RMS of PLOLO→board `pds`
prediction error — not a sampling σ; READ-FIRST §3 now says so and requires paired target-bootstrap intervals and the
regression / fitted-calibration / independent-validation distinction. 264's `analyse264.py` compared its pessimistic
price with the shape-only gate for an amplitude-changing cell; fixed to the cell's own gate, re-run, verdict unchanged
(−0.52 vs +0.0089). 264 rejects the majority rule IN ITS CONSTRUCTION (`kolf2` through an operator fitted without it), not
the three-production-donor rule — READ-FIRST v3.4 words it so; 265 (running) gives the `kolf2` +0.0121 its paired
bootstrap and a frequency-preserving permutation test over drop assignments (GPT's re-analysis: 95 % ≈ [−0.003,
+0.028]). 266 (running) re-runs 260's Wilcoxon with gemgroup- and depth-matched controls, the confound 260's pooled
control pool left open. The novelty guard now computes the prediction's CONTENT digest inside the container
(`zstd -d` stream + h5py/anndata over the fields the scorer reads) and refuses a scored prediction behind any sidecar;
backfill of the ledger's on-disk packages running (`E256/work/backfill-content.log`). READ-FIRST v3.5. No submission
on 09-06 stands; the reset is a quota reset, not a deadline.

### 2026-09-05 16:10 UTC+3 — 266: batch-matched controls do not change 260's NO-GO; the novelty backfill's descriptor leak fixed

266 (23 min): controls drawn from each target's own gemgroups (composition TV 0.39 → 0.005) and UMI-matched within
gemgroup (pool CV 0.440 → 0.378) leave the 50 candidates' median significant-gene count at 5.7 / 6.3 vs 6 pooled and
17–18 / 50 clearing the 10-gene gate vs 17; the KA row reproduces 260's counts to the digit. 260's NO-GO stands.
The content-digest backfill stalled after five packages: `mkstemp`'s descriptor was never closed, so each unlinked
13–19 GB scratch file stayed allocated until the process exited (119 GiB free again on exit); fixed (bdf35ef),
backfill re-run for the seven remaining packages. READ-FIRST v3.6.

### 2026-09-05 16:20 UTC+3 — 265: the majority result under proper statistics — p = 0.011 on kolf2 (one fold, uncorrected for the screen), at chance on 3/3 cancer folds; H2's majority variant closes on its own numbers

Subagent (22 min; 25/25 regression rows; GPT's re-analysis reproduced to MC error). Per-target `pds` is not additive
through the pipeline (the global preserve scalar in `Operator.apply`), but totals cancel to 0.04 null-SD on kolf2, so
the permutation ran on a labelled additive proxy. kolf2: bootstrap +0.01209 [−0.0035, +0.0278]; frequency-preserving
permutation p = 0.011 (z +2.3 vs null SD 0.0068); 264's three random seeds were an upper-tail draw (P 0.0047). Cancer
folds (matched): p 0.68 / 0.12 / 0.50; letter: the arm is significantly worse than its own control. 263's two-donor
cells: noise-sized. Recommendation adopted: no re-fitted three-donor validation. READ-FIRST v3.7.

### 2026-09-05 16:05 UTC+3 — content digests backfilled for all 12 scored packages; the guard is content-bound end to end

`scripts/backfill_ledger_content.py` (third run, after the descriptor and import fixes): 12/12 ledger records now
carry `content_sha256` (18–70 s per package, scratch clean, 119 GiB free). A re-prep of any scored prediction under any
CLI version, behind any sidecar, is refused by content.

### 2026-09-05 16:40 UTC+3 — The long game: research programme written, four feasibility workstreams launched in parallel

After the funder's brief (beat the challenge AND produce scientifically beneficial results; days of parallel
experiments are fine) and the second review's programme proposal, `docs/RESEARCH-PROGRAMME-2026-09.md` sets the
organising question (which responses are conserved, which depend on the recipient's state, when is a transfer
prediction unreliable), two evaluation tracks kept apart, binding rules (experimental unit, reserved confirmation
data with a ledger, refit inside the split, strong baselines and wrong-recipient controls, multiplicity, three scaling
axes, serial heavy jobs), and four workstreams with feasibility mandates: 267 replicate-validated oracle advantage,
268 Broad-16 context-residual map, 269 leakage-safe utility router, 270 data reconnaissance. RESEARCH-METHOD §1
scoped to the challenge score. READ-FIRST v3.8.

### 2026-09-06 04:30 UTC+3 — Programme round 1 read (the coordinator was not re-invoked when its watcher fired; the four subagents finished 14:00–15:30Z on 09-05): 267 DISAPPEARS, 268/269 closed at their constructions, 270 a weak positive to confirm; 271 launched

267 (77 min): the donor-selection oracle's +0.046 is the held reference's own sampling noise — chosen on half A and
scored on half B it reads −0.018 [−0.027, −0.009] vs the equal-weight mean (retained −0.32), coin-flip donor agreement
between halves (53–59 %); the reference's split-half residual cosine at 70–95 cells is +0.047 / +0.028 / +0.028 by
batch (32–34 % of Orion targets ≤ 0); only the most-reproducible tercile carries selection information and averaging
still wins there. The programme's first scientific result: response ESTIMATION must improve before any selector.
268 (51 min): Broad-16's line-specific component is reproducible (11.5 %, 4× null, 79 % across guides) but rank-27 and
not predictable from 2–3 basal programmes under LOLO (cosine +0.018 vs +0.001 wrong recipient, zero gain); reserved
lines JVE127/MG63/SLR23/UMRC3 (sha in the ledger); panel overlap 0. 269 (64 min): a utility router at its shuffled
controls on held context + held family; 066's same-target lookup is a PLOLO leakage path. 270 (38 min): GPSAdb 2.0
covers 96/300, 50 replicated; the frozen sign-stability feature AUC 0.89 on 4 anti-aligned targets, ρ +0.36 [0, 0.63]
— confirm on held targets; CMap at chance. Board 09-06: leader 0.22961661, gap +0.01352; the 09-06 slots unspent by
decision. Next: 271 — donor-row reliability at the donors' own depth and noise-aware estimators (reliability-weighted,
empirical-Bayes, Wiener) priced on PLOLO with a13's construction as the baseline. READ-FIRST v3.9.

### 2026-09-06 05:00 UTC+3 — 271: noise-aware donor estimation does not clear (0/10 cells); the donors are as noisy as the references; the operator has learnt its own training noise (→ 272)

Subagent (31 min, 98 turns; retrodicts BR249 at 0.0e+00 vs 257). Donor-row split-half reliability 0.09 / 0.05 / 0.05
at full depth; reliability picks the better donor 52–61 % of the time (break-even 66–70 %). Estimators: reliability-
weighted mean −0.0037 [−0.0103, +0.0028]; EB row shrinkage +0.0014 = its shuffled control, retention 31 %; Wiener
filters −0.013 pessimistic (nmae +0.0185); a13 re-implemented reproduces the board's loss. Finding kept: the fold
operator fitted on full-depth donor rows anti-correlates those rows' halves (−0.13 … −0.24, 97–100 % ≤ 0) — it has
fitted the sampling-noise realisation of its training rows; the held line's halves through it read +0.03 … +0.07.
The 400-cell ceiling: Orion-like references reach ρ 0.07–0.20 at 400 cells (√ρ 0.27–0.44), H1-like Flex 0.48 (√ρ 0.69).
READ-FIRST v4.0. Next: 272, a noise-corrected / cross-fitted whitening operator measured on PLOLO.

### 2026-09-06 06:15 UTC+3 — 272: un-learning the operator's training noise does not move the transfer (0/11); the whitening covariance is a noise covariance with a signal top decile, and pds lives in the noise-dominated directions

Subagent (71 min, 136 turns; KA-0 reproduces 271's anti-correlations; KA-2 the standard refit reproduces the cached
held operator). In the fitted eigenbasis the independent noise estimate exceeds the eigenvalue on 77–79 % of the
directions of every operator (production 78.7 %); those directions carry 87–92 % of the whitened prediction energy at
12–282× amplification. Wiener damping −0.029 … −0.061 pds on 3/3 folds at every κ; out-of-sample-variance whitening
c = 1: −0.0022 [−0.0040, −0.0004]; the cross-fitted operator (half-row fits) confirms the in-sample noise learning
(+0.022 / +0.007 out-of-sample vs −0.142 / −0.243 in-sample) and moves pds by +0.0001 … +0.0022 with intervals across
zero. Reading: the target-specific signal (5–10 % of a donor row's energy) is spread over the same low-variance
directions as the sampling noise, so the retrieval metric needs their amplification; denoising by direction removes
the signal with the noise. Construction closed; READ-FIRST v4.1. Next: 273 — a reliable-reference instrument (h1-test,
Flex, 400 cells, ρ 0.48, with the H1 union OFF so it is not leaked) to re-read the noise-limited verdicts.

### 2026-09-06 06:50 UTC+3 — 273: the reliable-reference instrument VALIDATES against the real scorer (5.6e-17); the Flex-400 regime is prediction-limited; the three PLOLO verdicts blur; the guide construct is a real unit of variance

Subagent (36 min, 87 turns). The per-target statistic (bulk log1p at TS 50,000, real control both sides, target
columns removed, cosine, midrank ties) reproduces cell-eval2's per-target `pds_cosine` on 221's cloud-scored artifact
to 5.6e-17; the union-OFF composition reproduces 170/221's effect to 1.2e-17; CD4 re-extracted for the 85 targets.
h1-test split-half residual reliability median +0.243 at 200/200 (0/85 ≤ 0; Spearman–Brown 0.39–0.50 at 400) vs
Orion +0.028; guide-construct halves +0.014 (2/5 ≤ 0) — the construct is a unit of variance and the challenge truth is
one construct per target. The knee at 0.02 cosine is the rank statistic's, not noise's. On the incumbent (union off):
12.9 % below chance, 48 % marginal, 39 % saturated; median transfer cosine 0.04. Re-reads: operator damping, EB
shrinkage and the CD4 deflation all BLUR as `pds` effects (intervals across zero on 85 targets) while every cosine
effect is larger than on PLOLO. READ-FIRST v4.2. Next: 274 (the union's directions with an nmae-neutral renderer —
the one competition candidate with a measured positive direction) and 275 (the findings write-up).

### 2026-09-06 07:20 UTC+3 — 275: docs/FINDINGS-2026-09.md written (447 numbers verified); one material correction — h1-test's "0.481 measured at 400" is fitted

The findings note stands on its own: seven findings with unit, n, interval and limit, the established/not table with
the preregistered confirmation each would need on the reserved data, and a methods appendix. 275 found three
headline/table discrepancies and reported them without fixing: (1) 159's h1-test 400-cell reliability 0.481 is an
extrapolation (measured 0.330 at 200/200; linear 0.478, power 0.458, Spearman–Brown 0.39–0.50) — READ-FIRST corrected;
(2) R273 §4.1 mislabels the floor family's production energy (94.9 / 83.1 / 67.7 % belong to c 0.25 / 0.5 / 1);
(3) R272's cross-fit upper bound rounds to +0.0056, not +0.0057. READ-FIRST v4.3. 274 still running.

### 2026-09-06 07:40 UTC+3 — 274: the union's direction cannot be rendered without its magnitude cost (0/39 buildable, 0/11 oracle); the last measured positive direction on cd4svd5 is closed

Subagent (47 min, 108 turns; 216/216 KA rows per fold, ARM6's identity leg to 3.5e-08, both board points retrodicted).
Four renderer families (pool-restricted magnitude match, whole-row rank transplant, second-residual, DLD-1 subset) at
three doses: every buildable cell prices below the pessimistic gate and below its un-rendered parent; on the scored
DE set direction and magnitude are one object. READ-FIRST v4.4. The competition track has no priced candidate on
this base; the programme's next stage is data integration for the organising question (see the next entry).

### 2026-09-06 08:50 UTC+3 — 276 (on Opus, after Fable's [bio] false positive): the conservation map reproduces and predicts transfer under LOLO; confirmed on the reserved HEK293T half-B truths (claim C1, ρ +0.361 [+0.257, +0.458], PASS)

276 died on Fable in one minute with `Details: [bio]` (the deterministic safeguard; AGENTS.md's recorded fix) and was
relaunched on Opus 5 (35 min, 103 turns; six KA rows at the digit; statistic frozen and hashed before any cosine).
Conservation over 19 lines: split-half ρ +0.273 (cells) / +0.435 (disjoint line sets); LOLO prediction of per-target
transfer ρ +0.387 (PLOLO), +0.239 (Flex-400), AUC 0.694 for the anti-aligned group; over two-donor agreement +0.066
[+0.004, +0.133] on one of three readings. The coordinator preregistered claim C1 in the confirmation ledger and made
the single access to the reserved half-B truths: n = 293, ρ +0.361 [+0.257, +0.458] — PASS. F8 added to
docs/FINDINGS-2026-09.md. READ-FIRST v4.5. Next: 277 (rendering on the identifiable low-conservation rows, priced).

### 2026-09-06 09:20 UTC+3 — 277 (Opus, second [bio] false positive on Fable): rendering the low-conservation rows differently does not clear; 276's competition-facing framing closed

Subagent (33 min, 113 turns; 28 KA rows, 25 at 0.0e+00, BR249 retrodicted, F4's Flex-400 cosine gains reproduced to
0.1 %). The un-deflated CD4 residual on the bottom-30 % conservation class reads −0.0043 [−0.0092, +0.0006] `pds` on
PLOLO; no-CD4 on the class and the soft floor on the class do not clear either; the mirror and random-set controls
behave as controls. The competition track has no candidate on this base; the scientific track's next axis is depth.
Mandates in this family now launch on Opus directly (two consecutive [bio] false positives, AGENTS.md's fix).
READ-FIRST v4.6.

### 2026-09-06 09:50 UTC+3 — 278 (Opus): the depth axis on disk — 214/300 targets have ≥ 400-cell somatic measurements, 152 only in primary CD4 T cells; no download moves it

Subagent (26 min, 121 turns; six KA rows). The panel's deep measurements are the CD4 T-cell source (GWCD4i) we
already use as a residual; only 62 targets reach 400 cells in an immortalised line; 86 targets are below 400 in every
somatic source. Orion is fully extracted (7,943,468 cells), Replogle's raw is the on-disk cells, Feng cannot reach 400,
Pisces (Jurkat) is unreleased. The replicate unit worth acquiring is GWCD4i's per-donor cell-level data (4 donors × 3
states) — the biological-donor reproducibility of a deep perturbation response, which no fold on disk can measure.
READ-FIRST v4.7. 279 (ZNF claim, preregistered) still running.

### 2026-09-06 10:15 UTC+3 — 279 (Opus): the preregistered zinc-finger claim is CONFIRMED on 17,106 held genes — and its honest scope is a family axis from core machinery to cell identity, state dependence rather than measurement

Subagent (40 min, 84 turns; prereg hashed before any held number; six KA rows, one mis-specified KA recorded and
replaced by the per-target check it could make, Spearman +0.97–0.98). Raw Δmedian −0.0105 [−0.0127, −0.0084], AUC
0.564; matched −0.0064 [−0.0117, −0.0008]; 0/100 random families; reproduced within one dataset and on the balanced
core (AUC 0.604). Not ZNF-specific (homeoboxes 0.603, bHLH 0.574; ZNF ranks 13–21/93 families). The deficit grows with
measured reliability (AUC 0.517 → 0.614) and survives reliability correction (74 %): state dependence. F9 in the
findings note; the 17,106-gene map is a by-product. READ-FIRST v4.8. 280 (GWCD4i donor replicates) running.

### 2026-09-06 11:05 UTC+3 — 280 (Opus): across four human donors a Flex perturbation response reproduces at +0.054, not H1's +0.243; the depth axis is closed on this disk; F1/F4's Flex-400 regime is an H1 property

Subagent (66 min, 123 turns; prereg hashed 06:57:35Z; 7 KA rows incl. the study's own DESeq2 donor pairs as a
cross-estimator, Spearman +0.70–0.85; 1.4 GB of panel rows extracted from GEO, nothing under data/). Across-donor
median +0.054 [+0.043, +0.062] at ~250 cells/side (29 % ≤ 0); donor is a real but marginal unit (+0.0033, 1/30 of cell
sampling, = the guide construct); three donors beat one by +0.0067 at matched budget; depth does not close the gap
(+0.054 at 400/side, H1's +0.243 needs ~4,500 cells/target); response magnitude, not cells, predicts reproducibility.
F10 added; READ-FIRST v4.9 records that 273's "prediction-limited" reading is H1-specific and the challenge's somatic
references are more plausibly in the +0.05 regime. Next: 281 (what predicts a gene's conservation, on the
17,106-gene map) and 282 (activation-state dependence within donors, on 280's extracted bank).

### 2026-09-06 11:35 UTC+3 — 282 (Opus): the preregistered across-state index was a construct artefact; deconfounded, a response keeps ~half its across-donor agreement when the state changes; the prediction from cross-line conservation is NOT confirmed; the reproducible part of a small pseudobulk is donor × guide

Subagent (24 min, 59 turns; prereg hashed before any across-state number; 5/5 KA). The frozen index read −2.68 —
same donor + same guide across states agree 5.6× better than across donors; the crossed design assigns that to the
construct. Deconfounded SD median +0.49 [+0.33, +0.61] (74 targets). ρ(SD, cross-line conservation): preregistered
+0.262 (wrong sign), deconfounded −0.144 [−0.403, +0.126] — not confirmed; conservation predicts the across-donor
ceiling (ρ +0.301), not across-state agreement. F11 added (the preregistered statistic retracted, the deconfounded
one reported). READ-FIRST v5.0. 281 (what predicts conservation) still running.

### 2026-09-06 12:00 UTC+3 — 281 (Opus): conservation is family-predictable (rho +0.81 over 51 held families) and only partly gene-predictable (rho +0.32, no better than reliability); F12; programme synthesis after rounds 1–3 written

Subagent (41 min, 107 turns; prereg hashed 08:15:18Z). READ-FIRST v5.1. The synthesis (RESEARCH-PROGRAMME §Synthesis)
lists what is established (F1–F12, C1), what is closed, and the next investments: a different cell type or stronger
responses for a resolving reference, a second construct per target in any new data, the write-up, and a held
competition track. The one remaining preregistered confirmation to run is F5 on the reserved Broad-16 lines (283).

### 2026-09-06 12:15 UTC+3 — 283 (Opus): claim C2 on the reserved Broad-16 lines — C2a FAIL by the frozen rule (1/4 lines at 268's size), C2b PASS, C2c PASS; both reserves are now spent

Subagent (16 min, 52 turns; claim written and hashed 08:46:12Z, ledger entry 08:47Z, first reserved byte 08:48Z, one
access to 08:58Z; 2.42 GB of reserved-line members fetched by range request into work/, never under data/). F5's
substance holds on held lines (the line-specific residual is detectable on all four against per-line nulls; the shared
response transfers; basal state does not predict the residual); its size at 268's threshold did not reproduce (1/4).
Recorded as a FAIL on C2a in the ledger and as an amendment to F5 in the findings note. READ-FIRST v5.2. Rounds 1–3
complete; nothing running; both reserves spent; the competition track held.

### 2026-09-06 16:10 UTC+3 — Heartbeat: all 24 programme experiments (260–283) finished and recorded; nothing running; the board moved — leader 0.23619 (gap +0.0201), second 0.22652 (gap +0.0104)

Both competitors above us improved during the two days we held (leader +0.0116 since 09-04, second +0.0071 since
09-05). Our incumbent is unchanged at 0.2160944 by decision: no priced candidate exists on this base (READ-FIRST §5–§6,
RESEARCH-PROGRAMME §Synthesis). READ-FIRST v5.3.

### 2026-09-06 16:40 UTC+3 — Third external review acted on: four verified corrections to the record (280's depth extrapolation, 282's state ratio, C1's provenance, 281's family split); C2 wording; the programme continues with 284–286

Verified in the sources: 280 thinned counts, not cells, and paired the donor reliability with the guide cohort's
cell count — the "4,500 cells" ceiling is retracted and the depth axis reopened; 282's "half the signal" ratio
conditions on 49 positive-denominator targets and the absolute difference +0.0080 [−0.0026, +0.0186] includes zero —
downgraded; C1's half-B cells were inside 276's full-truth evaluation and 267's scoring — relabelled split-half
robustness; 281 blocks on the first HGNC group — rerun as 284 with all memberships excluded. C2c reads "no detected
improvement for two frozen models", C2a's P = 0.11 is a miss, not more. READ-FIRST v5.4, FINDINGS amended (F5, F8, F10,
F11, F12), ledger C1 relabelled. Launching 284 (family split rerun + target-cluster bootstraps), 285 (the apples-to-
apples reliability benchmark on cell-level data: cells vs UMIs, guide / donor / study units, matched strata, agreement
AND challenge metrics), 286 (a cross-product signal-covariance estimator with a full-rank residual path, held
guides/assays → held contexts → production rendering). Slots stay unarmed; the board leader is 0.23619 (gap +0.0201).

### 2026-09-06 16:42 UTC+3 — Round 4 launched on Opus: 284 (family split rerun + remainder-vs-noise test), 285 (matched reliability benchmark with a bounded CD4 cell-level extraction), 286 (cross-product signal-covariance estimator, priced)

Mandates in `experiments/284-conservation-family-split/`, `experiments/285-reliability-benchmark/`,
`experiments/286-crossproduct-estimator/`; programme synthesis amended (the "external-only investments" claim
withdrawn; item 5: modelling capacity is not exhausted); Round 4 table appended. Slots unarmed; nothing priced.

### 2026-09-06 17:30 UTC+3 — 287 (novelty audit) lands in 22 min; FINDINGS v1.2; 288 magnitude check; programme rules 9–10

287 (Opus): NEW 4 / EXTENDED 11 / KNOWN 7 / CONTRADICTED 2 over F1–F12, C1, C2, with 40 references and four read-only
checks on the CD4 source's own supplementary tables (KA: its random controls at chance). Acted on: F10's surviving
"cells are not the binding constraint" sentence deleted; F5's headline narrowed to "not predictable from basal
EXPRESSION by a linear model" and the Broad-16 paper's genotype/HPV predictions recorded; F1a, F1c, F6, F10a, F11a
re-framed as re-measurements with their sources; F2 separated from Wang 2026's positive training-selection result;
F3a's mechanism cited as textbook, F3b kept as the new statement; F9's independent-estimator corroboration added;
F12a's +0.811 held until 284. 288 (coordinator, `experiments/288-f8-magnitude-check/`): F8 within magnitude tertiles —
at chance in the weakest tertile, partial ρ given magnitude +0.26 / +0.25 / +0.23 vs raw +0.41 / +0.39 / +0.37; no
gradient by reliability tertile. Programme rule 9 (depth-condition pass rules on multi-line sets) and rule 10 (audit
novelty before claiming it). 289 to follow: 268's LOLO with the paper's state features.

### 2026-09-06 17:40 UTC+3 — 289 launched (Opus): 268's LOLO rerun with the Broad-16 paper's own state features (p53 status, HPV-E6, lineage, proliferation / p53 programmes, arm CN if reachable) on the 12 working lines; the paper's two mechanistic predictions tested on our residual object; reserved lines untouched

### 2026-09-06 17:45 UTC+3 — 284 lands (Opus, 29 min): the family-split leak was real and did not inflate 281; the gene-level remainder is REPRODUCIBLE across disjoint source sets — "mostly measurement" withdrawn; FINDINGS v1.3

KA: 281's primary cell at 0.0e+00, 279's held map at 0.0e+00, 200-shuffle null (0/200). All-membership split ρ +0.334
[+0.216, +0.355] (281 +0.320); family axis under strict folds +0.805 [+0.679, +0.883] over 78 families; ZNF deficit
57 % predicted (84 % was a property of the leakier folds). Remainder test: two halves of the 28 sources sharing no
source and no line, remainder ρ +0.328 [+0.241, +0.373] vs ceiling +0.382, null +0.044 ± 0.018 (0/200), 20/20
partitions, +0.277 after reliability and magnitude removed, +0.322 on off-axis genes. F12 amended (three items);
READ-FIRST §6; programme Round 4 table extended with 287 and 289. Running: 285, 286, 289.

### 2026-09-06 17:50 UTC+3 — 290 launched (Opus): response-derived features from one half of the sources predicting the other half's conservation remainder (284's reproducible per-gene component); frozen predictions P1–P4 (proliferation / stress engagement and genericness → less conserved)

### 2026-09-06 17:50 UTC+3 — 289 lands (Opus, 1 h 40 min): the Broad-16 paper's own state features do not predict the residual transcriptome-wide (0/25 cells), but its P1 mechanism reproduces on our data and the LOLO finds it blind; 287's C-1 RECONCILED; FINDINGS v1.4

KA: 268's (c1) at 0.0e+00 (788 JSON entries, 984 × 45 CSV cells); the p53 dependency proxy validated on 1,150 DepMap
models (TP53 vs MDM2 effect r −0.74). S1 predicted-residual cosine +0.0110 / +0.0091, not separated from the wrong-
recipient control; full-response Δcos −0.0588 / −0.0236, worst in the deepest lines (ρ −0.72 with UMIs per
perturbation). P1: SKGII +0.10, C4I +0.07 top 2 of 11 on the p53 direction, permutation P = 0.018, E2F direction
opposite; UBE3A best-predicted of 82 targets under S1 with HPV as carrier, 12/12 folds. P2: wrong direction on 12
lines. Sign gating on the p53 axis +0.053 [+0.015, +0.089] over shuffled. Lesson recorded: low-cardinality state
features pass the shuffled-label control trivially — score against the wrong recipient. Running: 285, 286, 290.

### 2026-09-06 18:15 UTC+3 — 290 lands (Opus, 47 min): the reproducible remainder is carried by response GEOMETRY measured in other lines (ρ +0.317, at the ceiling; annotation −0.057); the preregistered mechanism inverted — broad, spectrally central responses are MORE conserved; F13; FINDINGS v1.5

KA: 279's map at 0.0e+00, 284's maps / remainders / ceiling / remainder at 0.0e+00, 200 shuffles (0/200, z +15). P1 PASS
(+0.3174 [+0.2897, +0.3914]), P2 PASS (+0.3948 over annotation), P3 FAIL sign inverted (+0.1184), P4 FAIL sign inverted
(+0.1147); replication partitions and 40/40 cells positive; without Orion +0.308; random matched gene sets recover 93 %
of the model; DNA binders the systematic exception. 284 §5.2's r_bar sentence corrected (four banks, all in half 1).
Running: 285, 286. Next: 291 (F13 on Broad-16's shared fraction and GWCD4i's published CD4-vs-K562 conservation).

### 2026-09-06 18:10 UTC+3 — 291 launched (Opus): F13 on two surfaces outside the 28-bank map — Broad-16's 12 working lines (disjoint lines, Cas9) and GWCD4i's published CD4-vs-K562 conservation (1,773 genes × 3 states); claims C3 / C4 frozen and logged in the ledger before any comparison

### 2026-09-06 18:35 UTC+3 — 291 lands (Opus, 23 min): C3 FAIL (low power, 82 targets), C4 FAIL (does not beat annotation) — F13's SIGN confirmed on two surfaces outside the map, its "beats annotation" headline withdrawn off-map; the map itself transfers (+0.27 on the GWCD4i authors' CD4-vs-K562 statistic); FINDINGS v1.6

KA: 290's cell +0.317385 (1.2e-07), 279's map, 287's random-control row, 268's decomposition at ≤ 8.2e-08. Claims C3/C4
logged in the ledger 15:13:32Z before any comparison. Broad-16: model +0.210 [−0.013, +0.426], gen_pc20 +0.295
[+0.070, +0.507], annotation +0.300, C_held +0.210. GWCD4i Stim48hr: model +0.203 [+0.130, +0.252] (z +8.4), partial
given the authors' magnitude +0.061 [+0.030, +0.156], gen_pc20 +0.186, breadth +0.128, annotation +0.270 (not beaten),
C_held +0.274; same pattern in Rest and Stim8hr. DNA-binder exception not reproduced (TFs below the model). Running:
285, 286.

### 2026-09-06 18:50 UTC+3 — 286 lands (Opus, ~8 h): the cross-product signal-covariance estimator DOES NOT CLEAR (0/12; held CD4 donors null; PLOLO best −0.0080 pds, −0.0207 pessimistic; Flex reference same sign); the mechanism survives an unbiased basis; F3 amended; FINDINGS v1.7

13 KA rows (263 / 257 / 271 / BR249 ≤ 2e-11; the moment identity 5e-17; Ŝ_shuf leading eigenvalue 7–19× smaller). Ŝ's
top-100 directions sit at principal cosines 0.52–0.70 to the fitted operator's — a different object — and 78–80 % of
directions have λ⁺(Ŝ) below the noise level while carrying ~81 % of the incumbent's whitened prediction energy; Ŝ's
top-100 carry 1 %. Every arm helps the anti-aligned group and pays more in the marginal band, shuffled control
included. Agent reported exceeding the memory cap (23 GiB, fixed to 4.7 GiB in-session; machine ≥ 73 % free) and a
pycache written into E265 (harmless). Nothing build-worthy. Running: 285.

### 2026-09-06 18:40 UTC+3 — 292 launched (Opus): declare reserve R3 — the unused single-cell perturbation datasets on disk, inventoried from metadata only, with a frozen access rule in the ledger (review item C.3)

### 2026-09-06 19:05 UTC+3 — 285 lands (Opus, 2 h 09 min): the matched reliability benchmark — cells and UMIs one axis (Spearman–Brown holds on real ladders), the guide construct the dominant irreducible term (ρ 0.09–0.41 CRISPRi; scorer pds to chance on a construct swap), batch free, a different study ρ 0.79, a different line 0.35–0.45, donor ≈ guide in CD4 and flat in the budget; CD4 cell extraction not admissible; F14; FINDINGS v1.8

KA: 267's medians at 0.0e+00, 273's +0.242920 recomputed at 0.0e+00, 159's 0.330039, 280's +0.0539 (1.5e-3), union
identities at 0. Two centring artefacts found by impossibility checks and corrected uniformly (thinned sides centred on
full-depth means; identity-varying sides). Two jobs at 11.8–12.0 GiB (serial, machine ≥ 73 % free). The competition
hypothesis F14 raises (the challenge's truths are construct responses; a same-library source carries what PLOLO cannot
see) is recorded unpriced; 293 checks the library-sharing fact. Running: 292.

### 2026-09-06 19:20 UTC+3 — Library fact checked: the 2026 panel shares 25/300 targets with the 2025 H1 library; no source on disk carries the 2026 constructs — the construct ceiling binds every zero-shot entry (FINDINGS v1.9); 293 launched (Opus): the ceiling in the challenge's own regime on the 2025 library's two-construct targets

### 2026-09-06 19:05 UTC+3 — 292 lands (Opus, 17 min): reserve R3 declared from metadata only — eight untouched scPerturb files, 174 held-map genes pooled, once-per-claim access rule in the ledger; 62 datasets USED / 42 UNUSED inventoried; agent breached its 6 GB cap (11.5 GiB obs scan, recorded). 295 launched (Opus): claim C5 on R3 — F8 off-map and F13's sign, one pooled access, file as the bootstrap block, modality as a stratum

### 2026-09-06 19:10 UTC+3 — 294 launched (Opus): the scientific write-up — manuscript draft, reproducible figures from saved tables, a claims ledger; sourced from FINDINGS v1.9, the confirmation ledger and 287's provenance verdicts; 293 / 295 integrated when they land (pending markers otherwise)

### 2026-09-06 19:45 UTC+3 — 295 lands (Opus, 20 min): claim C5 on reserve R3 — C5a PASS (the map predicts transfer into eight untouched datasets, ρ +0.502 [+0.153, +0.566] file-block; 0/200 shuffles; survives reliability, magnitude, within-file ranks), C5b PASS (F13's sign, `gen_pc20` +0.371 [+0.064, +0.415]); the map does not beat annotation; all eight files SPENT. 293 lands (Opus, 30 min): the challenge-library construct ceiling — ρ_guide +0.075 [−0.010, +0.288] on 17 two-construct targets, bimodal by knockdown (+0.687 where both constructs work, −0.003 where the weaker fails; two-construct targets the library's weakest, P = 0.0015); a failing construct still gives a reproducible transcriptome; zero-shot ceiling cos 0.21–0.57 vs the incumbent's 0.04 — the ceiling caps everyone and does not bind today (F14's sentence corrected); h1-test and the 2026 controls at ≈ 20 k UMI/cell vs 54 k in the 2025 files. FINDINGS v2.0. 296 launched (Opus): the knockdown stratification on 285's public libraries (is the construct term guide efficacy?).

### 2026-09-06 20:15 UTC+3 — 297 lands (Opus, 18 min): audit extension — NEW 3 (284's remainder, study-beats-construct, geometry-vs-identity ablation), EXTENDED 7, KNOWN 5 (F13's sign and F14's construct term are published: Wang et al. 2026; the CD4 source's own column); 293's "one truth in five" withdrawn (two constructs vs 0.8–4.7 % published off-target rate); FINDINGS v2.1; the draft's three provenance markers filled and its citation duties added (Wang 2026, Jolliffe 1982, Hartman 2026; tensions T6–T8)

### 2026-09-06 20:25 UTC+3 — The manuscript draft published as a private artifact (https://claude.ai/code/artifact/b7545115-c967-4dd0-b670-7329a788f142): pandoc-rendered PAPER.md + figure gallery (PNGs embedded) + the claims ledger with verdict chips; build script `experiments/294-writeup/build_html.py`. Pending: 296 → F14 wording → republish, then notify Ahmet.

### 2026-09-06 20:35 UTC+3 — 296 lands (Opus, 33 min): guide efficacy explains the construct term PARTLY (Jurkat cleanly; K562 at the library-relative cut with magnitude confounded; HepG2 not separated; Broad-16 no split; primary CD4 a powered null); cross-gene null correction — the CRISPRi construct terms are ≈ 0.03–0.08 null-corrected, only Cas9's 0.79 survives; failing constructs reproducible at scale, no shared vector signature; failure rates 19–27 % by the libraries' own calls. F14 amended (FINDINGS v2.2); the draft's F14 paragraphs, abstract, plain summary, §4.3 and T7 updated; CLAIMS C-50 amended, C-51–C-53 added. Round 4 complete.

### 2026-09-06 20:40 UTC+3 — Draft republished (v2, Round 4 complete) and sent to Ahmet; 298 launched (Opus): the reviewer-facing checks on saved outputs — C5's single-file jackknife, F13's hallmark-only fit vs random matched sets, T8's study-vs-construct in raw cosine, T7's seed scan of the 2025 library's failing spacers

### 2026-09-06 20:50 UTC+3 — 298 lands (Opus, 15 min): C5a holds without Adamson (+0.268), C5b marginal without it (file-block bound +0.006; gene bootstrap through zero); hallmark-only fit +0.158 vs random matched +0.099 (identity predicts alone; redundant given geometry — F13 narrowed); study vs construct in raw cosine +0.168 vs +0.011 on the same 112 perturbations (T8 settled); the seed scan is not testable on disk (no spacer sequences, no genome; 13/34 arms fail). FINDINGS v2.3; draft and claims updated (C-54); page republished.

### 2026-09-06 20:55 UTC+3 — Fourth external review (20:30 Ankara): both 09-06 slots unused; nothing ready; proposes ONE bounded competition pilot tonight — the whitener fitted on RELIABILITY-SELECTED training rows (row count, donor proportions, axis, coverage fixed; matched-random and reversed controls; per-fold) — and zero slots unless it clearly clears. Agreed: 299 launched (Opus; hard stop 21:30Z; magnitude-selected rows added as the third control; the Flex re-read and 272's spectrum row included). Slots stay unarmed; the decision follows the reading. (The review also noted a provider safety error in this session at 20:37 Ankara — the session recovered; no work was lost, the 298 fold-in landed at 20:50.)

### 2026-09-06 23:10 UTC+3 — 299 lands (Opus, 2 h): the reliability-selected whitener DOES NOT CLEAR (−0.0675 pds, 0/3, below every control; a row-count effect below full rank, null at full rank); post-hoc rank-preserving weights give a monotone dose curve with the OPPOSITE sign (tilting away from reliable rows +0.00094 3/3 shape-only, mirror control passed, same sign on Flex) — priced 0.74 × the gate with retention 47.8 %: not a candidate. DECISION: both 09-06 slots left unspent by rule; nothing armed. F3 amended (FINDINGS v2.4); READ-FIRST §5 / §6; the draft's §3(f) and §3(g) updated and republished.

### 2026-09-06 23:15 UTC+3 — 300 launched (Opus, ≤ 3.5 h): the spread mechanism behind 299's reversed tilt, as three preregistered arms at the incumbent's row count (the tilt dose inside the shape-only band; its interaction with 062's prune; a direct low-leverage-row up-weighting) with mirror and random-weight controls, the price, and the Flex reference with an interval. A clearing arm would be a 09-07 candidate; nothing is armed.

### 2026-09-07 01:20 UTC+3 — Board 09-06 21:53Z: leader `cqawesome` 0.23893 (+0.00274 in a day; gap +0.02284), second 0.22652 unchanged, `cd4svd5` 3/740. READ-FIRST §1 updated. 300 still running.

### 2026-09-07 02:10 UTC+3 — 300 lands (Opus, 2 h 50 min): the reversed tilt as a mechanism — monotone, 3/3 at every realisable dose, all three controls pass, NO dose clears (the shape-only band edge at d ≈ 1.1–1.2; price and amplitude rise together); spread is 72 % collinear but a quarter of the effect; the prune interaction unresolved; 062's prune not established on either instrument; the Flex interval contains zero. The operator-composition thread is CLOSED. F3 amended (FINDINGS v2.5); READ-FIRST §5 / §6; programme: the competition track's state after Round 4 written down; the draft's §3(f) / §3(g) updated and republished. Nothing armed for 09-07; no new agent launched overnight — the next investments need new data or a new mechanism.

### 2026-09-07 02:25 UTC+3 — 301 launched (Opus, ≤ 2.5 h): metadata acquisitions for the write-up — the 2025 sgRNA library spacer table, promoter sequences, the JoungZhang ORF id table (a future reserve candidate) — and, if the first two land, T7's seed scan of the 13 failing constructs.

### 2026-09-07 02:55 UTC+3 — 301 lands (Opus, 50 min): the 2025 library's spacers found (Replogle 2022 Table S1; all 34 arms), CAGE-TSS promoters fetched (KA 68/68; a GRCh38 attempt rejected by its KA), JoungZhang TF-ORF ids resolved (3,136/3,369; R4 candidate noted in the ledger, X never read), and T7's seed scan run: no off-target enrichment among failing constructs (4/13 vs 4/21, P = 0.68) — the construct term is knockdown failure. FINDINGS v2.6; the draft's T7 updated; republished.

### 2026-09-07 06:35 UTC+3 — 302 launched (Opus, ≤ 4 h): claim C6 on reserve R4 (the JoungZhang TF-ORF atlas), the single access — does overexpression reverse knockdown in the same hESC line (249 TFs also in the map's H1 CRISPRi screen; 148 with ≥ 100 ORF cells), and do conservation (C_held) and knockdown geometry (gen_pc20) predict the mirroring? Preregistered in the ledger before any X read.

### 2026-09-07 06:55 UTC+3 — 302 lands (Opus, 30 min): reserve R4 declared and spent on claim C6 — VOID by the prereg's KA-E clause (the atlas's own count columns are inconsistent; X holds half-integers); descriptively, overexpression and knockdown responses of 148 TFs in hESC are ORTHOGONAL at the atlas's depth (median −0.0007 vs null +0.0008; ORF reliability +0.025 caps any cosine at 0.11); by-products: the SB depth law for ORF overexpression (c ≈ 9,600 cells), a +0.021 isoform-construct term, gen_pc20 predicting the overlap's size (+0.38) not sign. F15; FINDINGS v2.7; the draft updated (§3(d), §4.4) and republished. All four reserves spent.

### 2026-09-07 08:25 UTC+3 — Coordinator stall and handover note. After 302's landing was recorded and Draft v7 republished (03:58Z), this session drafted the Round 5 mandates 303 (claim C7 on the Nakatake 2020 rows, declared reserve R5 at 04:20Z, unread) and 304 (F15 novelty audit), then stalled without launching anything; the ledger's R5 heading had been written twice (the truncated copy removed 05:25Z). At 05:13Z Ahmet asked the Codex session to take over; it wrote HANDOFF-MANDATE.md for 303 (pre-access only) and 304, a MANDATE.md for 305 (T7 expression test), and launched 304 through scripts/claude_subagent.sh at 05:18Z — pid 1508 exited at once with a 0-byte log (the watcher reported DONE with no result). State at 05:25Z: nothing running; R5 unread; slots for 09-07 unarmed; board unchanged (leader 0.23893, cd4svd5 third of 741). Who coordinates Round 5 is Ahmet's call; asked.

### 2026-09-07 08:33 UTC+3 — Codex takeover is active; 303/304/305 verified running

Ahmet's instruction in the Codex task is explicit: "continue wherever it left
off ... take over." Codex is the coordinator. Recovery commit `4289412` preserves
the interrupted drafts, repairs the execution map and gives each worker a
disjoint mandate. The ledger recovery was included in the predecessor's
`c18d15c` handover commit and is preserved.

The detached startup attempt produced no events; it was not completion.
Persistent-shell relaunches now have live processes AND actual trace activity:
303 43 tool calls, 304 49, 305 28 at the 08:33 check, no terminal event yet.
303's metadata and loader checks reproduce the source hash, 100-TF cohort and
15,047 saved cross-source cosines; R5 remains locked and unread. 305 has streamed
the three already-used 2025 files serially (reported peak 3.79 GiB); its results
are pending known-answer checks and coordinator audit.

`scripts/round5_status.py` distinguishes launch failure, interruption, running,
pre-access report and terminal result awaiting audit; 8 regression tests pass.
300's scientific scope/units review is in
`experiments/300-spread-tilt/COORDINATOR-REVIEW-2026-09-07.md` (NO-GO unchanged).
305's primary empirical ranks lack the resolution for 11-test Holm correction;
its nominal alternative is descriptive only, as recorded in
`experiments/305-t7-expression-test/COORDINATOR-DECISION-2026-09-07.md`
before the coordinator reviewed candidate expression outcomes.

Live board/CLI check at approximately 08:29: incumbent
`Xv0ayvuoHpsKYWv2iu42` still published, 0.2160944022, rank 3 of 744; leader
0.2389319449; limit not reached, zero pending uploads. No prediction, submission
or public artifact was changed. The old "nothing running" statement above is
superseded by this checkpoint.

### 2026-09-07 09:05 Ankara — Codex completes Round 5

Recovered 303–305 workers completed. 304/305 were independently reconstructed
and reviewed; integration commit `e9220e0` preserves the original reports and
records the inferential corrections. C6 remains VOID. 305 has two nominal
primary repression observations, zero full-11 Holm discoveries; the post-hoc
Poisson-binomial headline is not accepted as confirmatory evidence.

303's complete corrected protocol/code was committed in `45788c6` before
the one R5 read at **06:01:00Z**. Forty guard and 23 analysis tests pass,
seven loader regression rows reproduce, and the memory gate passed (86% free,
29.62 GiB free+inactive, 1510.56 MiB swap free). Validity passes; R5 is spent.
**C7a PASS; C7b/c/d FAIL.** Median −.02633 [−.05203, −.00949],
41 negative versus 13 positive .05 tails, four-cell Holm p=.001999.
`scripts/audit_c7.py` reproduces the full cosine matrix at difference zero
and all four p values from the saved product without reopening the source.
Full interpretation: `experiments/303-second-overexpression-library/COORDINATOR-RESULT.md`.

All five takeover workers finished; original wrapper PIDs no longer exist.
No new experiment or automation launched. Live CLI still reports cd4svd5
published at rank 3, 0.2160944021991198. Nothing armed, built, submitted,
cancelled or changed in production. Findings v2.9 and the local manuscript
are updated; no external republication. Resume from the completed checkpoint
in `docs/TAKEOVER-2026-09-07.md`, not the historical launch notes.

### 2026-09-07 16:52 UTC+3 — Draft v8 republished (Claude session). The private manuscript artifact now carries Codex's Round 5 revision of PAPER.md / CLAIMS.md (C7a measured on R5; F15's stronger readings withdrawn per 304; T7 unresolved per 305); the three "not republished" markers replaced with the republication stamp (13:48Z). Watcher re-armed after its Monitor died at 12:21Z (exit 144). Nothing running; slots unarmed; board unchanged (leader 0.23893; cd4svd5 third of 762).

### 2026-09-07 18:30 UTC+3 — Round 6 launched (Claude coordinator, three Opus workers, after verifying Codex's Round 5 takeover)

Resumed from `docs/TAKEOVER-2026-09-07.md` at HEAD `89d38f7`. Verified before launching: no worker or research
process was live (every `claude-run.pid` dead); two orphaned heartbeat processes from dead sessions killed; the
working tree carries only leaderboard snapshots and old run logs; board re-read from the 12:53Z snapshot — leader
`cqawesome` **0.238932** (flat since 09-06 21:53Z, ~15 h), `cd4svd5` **0.2160944** third of **762**.

**Slots: both 09-07 slots are unspent and unarmed. No candidate exists, so nothing is armed** (the standing rule:
a slot buys information toward 1st, never a better row).

Launched, disjoint write scopes, all `CLAUDE_SUBAGENT_MODEL=opus --effort max`:
- **306** (≤ 3 h) reconnaissance — what data and what published method exists today that did not on 2026-08-30
  (X-Atlas/Pisces, Arc releases, GEO/scPerturb/Zenodo since 08-01, preprints), plus what is publicly knowable about
  the 0.22–0.24 methods, plus one falsifiable check against data already on disk.
- **307** (≤ 3.5 h) the T7 repression leads done as 305's coordinator review specified: ONE aggregate
  candidate-enrichment statistic with the **construct** as the unit, against a joint null preserving seed/PAM,
  promoter architecture, GC, mappability and expression, with the candidate-specific exclusion 305's null omitted.
  Preregistered as C8 before outcomes; "unresolvable at this resolution" is an allowed verdict.
- **308** (≤ 4 h) the competition-facing arm: weight each donor row by **the on-target knockdown its own construct
  achieved**. The mechanism is ours and measured — F14 (construct dominant) + 296 (that term is efficacy) + 293
  (**a failing construct still gives a reproducible transcriptome**), which is exactly why every
  reliability/norm/agreement weighting tried (016, 097a, 262/263, 299) cannot see it. Decisive control: reliability
  weights. 097 route (d) closed |kd| as a prediction target, not as a weight.

Four candidate ideas were checked against the record and NOT launched because each is already closed at its scope:
conservation as a lever (276 §4 limit 4, 263), the utility router (269), shrink-toward-panel-mean (097, 12/12), and
Arc's pretrained STATE (182, NO-GO). Recorded in the programme note's Round 6 block so Round 7 does not re-propose
them. Heartbeat re-armed as a persistent Monitor in this session. Nothing built, submitted, armed or downloaded.

### 2026-09-07 19:20 UTC+3 — Two verifications and one strategic reading, while Round 6 runs

**Production integrity after the multi-runtime handoff: CLEAN.** Codex's five Round 5 commits touch `scripts/` only
to ADD `round5_status.py`, `test_round5_status.py`, `audit_round5.py`, `freeze_c7.py`, `audit_c7.py`. `git diff
1f50291..HEAD -- scripts/ src/` lists those five files and nothing else; no production builder, pipeline, postbuild,
census, preflight or identity file changed, and the working tree is clean under `scripts/` and `src/`. The October
sequence (244's rehearsal, 235/240/247's guards, 256's pinned `V5_abc` wiring) stands as rehearsed.

**The competitive reading, recomputed on the 09-07 12:53Z board and written into READ-FIRST §1.1.** The gap is still
transfer cosine and it is now sized: ours **0.09728**, `cqawesome` **0.14456** (197's inversion under §4.1c's
validity guard), so cosine explains **76 %** of the +0.02284 gap. **Clearing the +0.0089 gate needs +25 % relative
cosine; matching the leader needs +48.6 %.** We are 3rd of 762 while holding only the **9th-best cosine in the top
16** — 4.1c's "we convert cosine better than anyone" has sharpened, not changed. The leader gained on **all six**
members over 33 submissions (a real model improvement, not metric gaming) and has been flat since 09-06 21:00Z.

**What that rules out, with 225's number.** The public-source channel SATURATES at ~1.8× (asymptote +0.0040;
thirteen sources deliver +0.00420) — *"+4.6–5.5 % of board cosine against `cqawesome`'s +106 % climb. A factor of
twenty."* So "acquire another public donor line" cannot close a +25 % gap; only data that changes the CONSTRUCT term
(the challenge's own library in a public screen, a second construct per target, the held lines themselves) is worth
paying for, which is what 306 is scoped to find. **The honest bar: nothing in this repo has ever produced more than
+5.16 % relative cosine (224, killed by `nmae` in 228).** Any future proposal must say how it clears an order of
magnitude more than our best-ever arm.

Also checked and closed on the way (all with measured numbers, so Round 7 does not re-propose them): CMap as donors
(176 NO-GO, and 270's cross-line sign stability does not track `pds`); the line-specific component on more lines
(289 needs lines we do not have; CMap's responses do not transfer); the unassigned-cell depth lever (278: the 71 %
unassigned cells are Feng's iPSC files, not our donors); CD4 depth (already at its measured argmax, `f*` 1.018, 216).

### 2026-09-07 19:10 UTC+3 — 306 lands (Opus, 1 h 20 min): NO new data exists, and the frontier is not using data we lack

**Part A — zero acquisitions clear today.** Between 2026-08-30 and now, **no perturbation dataset touching the 300
targets, the six challenge lines, or the challenge's own library became downloadable.** X-Atlas/Pisces is still an
empty HuggingFace repo (`usedStorage 0`, three files, unchanged **165 days** since 2026-03-25; README still "Coming
Soon" for the HCT116/HEK293T re-release, the HepG2/iPSC validation subsets and **Jurkat Resting + Active**); X-Cell
still has no weights and its preprint is unrevised since 2026-03-20. GEO (six queries), ArrayExpress/BioStudies
(22 hits, newest a bacterial CRISPRi set), Zenodo (426 records) and scPerturb (no release since 2024-08) return
nothing qualifying. Arc has published nothing new: the *Cell* front-matter piece (2026-08-26) and the STATE journal
version (2026-08-31) carry no data, and the timeline is confirmed unchanged — **Oct 22 test release, Nov 5
deadline**. The only post-freeze item is an ASD cortical-neuron CRISPRi preprint (bioRxiv 2026-09-02, 1,408 genes,
hESC-derived neurons) with **no recoverable accession**. Two conditional triggers recorded, nothing downloaded:
(1) Pisces uploading its Jurkat split, (2) an accession resolving for the ASD screen. **Do not download on a
preprint alone.**

**Part B — the frontier is silent, and it stands on our own data.** Of 769 entries, 292 write a description, but
**zero of the two entries ≥ 0.220 do**, and the best-scoring entry that describes any method is **rank 42 at
0.17543** (an emission recipe, not biology). Every rival who does describe a method names **only sources already on
this disk** — Replogle K562-GWPS, X-Atlas/Orion HCT116/HEK293T, GSE314342 primary CD4 (the same study as our GWCD4i
donor), the 2025 VCC H1 context — and at least four teams have independently reached a near-verbatim description of
the incumbent's own construction (whitened multi-source unit-norm fusion, common-response centring, count
transport), scoring **0.136–0.150 against our 0.2161**. No public repository is near the frontier (the best
published trajectory is 0.0184 → 0.0766 → 0.0649). **Whatever the top two are doing, they are not doing it with
data anyone has published.** Treat every description as an unverified self-report.

**Part C — a rival's public claim checked against our own files.** `Testh` (0.14954) asserts the 2025 H1 context
covers "21/300 panel targets". Our settled 25/300 **reproduces exactly** (KA-1 PASS; panel 300, 2025 splits
150/50/100, union 300; the 25 genes enumerated with file digests). No defensible reading gives 21 (bank-as-built 20,
library ceiling 25, per-split 13/4/8). **Nobody has more of the 2025 H1 context than we do**, and 293's ceiling
question stays open because no public screen uses the 2026 library.

**Board, live at 15:48:55Z (769 entries):** leader `cqawesome`/`m0907h` **0.2393600**, submitted 15:29Z **today** —
they are still moving; ours 0.2160944, rank 3, gap **+0.0232656**.

**What it changes.** Round 4's recommendation was "new data OR a new modelling idea". The data half is **blocked by
fact** until Pisces uploads. Combined with 225's measured saturation of the public-source channel (+5 % of board
cosine), the conclusion is now doubly supported: **the leader's advantage is method, not data.** Round 6's weight
falls on 308.

### 2026-09-07 19:35 UTC+3 — The gap trajectory, read off the full snapshot series (the 3-h snapshotter is alive; 15:53Z captured, 770 teams)

| window | leader | us | gap | note |
|---|---|---|---|---|
| 09-02 06:30Z | 0.21687 | 0.17336 | +0.04351 | we were 25th of 581 |
| 09-05 03:13Z | 0.22792 | **0.21609** | **+0.01182** | our narrowest gap; `cd4svd5` published 09-05 00:55Z |
| 09-07 15:53Z | 0.23936 | 0.21609 | **+0.02327** | 3rd of 770 |

**Read it correctly.** Over the whole window the gap *narrowed* (+0.0435 → +0.0233) because **we** climbed
+0.00793/day while the leader climbed +0.00417/day. But **our score has been constant at 0.2160944 since
2026-09-05 00:55Z** — by rule, since no arm has cleared the gate — while the leader has added **+0.01144 in the
last 2.5 days (+0.0046/day)**. So the *recent* motion is the gap widening at ≈ +0.0011/day, and 306's "still
widening" is correct for that window.

**Do not extrapolate linearly.** A naive linear projection puts the leader at 0.42 by Oct 22, which is impossible:
the **empirical ceiling — field best on each of the six members, 0.28883 today — is an upper bound assembled from
different teams and is not attainable by any single submission.** At their recent rate the leader would meet even
that bound in roughly a month, so the honest expectation is that they plateau before the Oct 22 test release.
The decision-relevant quantity remains ours: **+25 % relative transfer cosine to clear the gate.**

### 2026-09-07 19:15 UTC+3 — 307 lands (Opus, 1 h 15 min): the strongest scientific result of the round, and it aborts itself honestly

**C8 (aggregate candidate-repression enrichment, construct as unit): UNRESOLVABLE-AT-THIS-RESOLUTION.** The design
the 305 review asked for was executed: ONE verdict cell, the **targeting construct** as the experimental unit, all
**34** targeting constructs (83 construct–gene pairs, no construct dropped for lacking a strong candidate), and a
joint null that resamples *candidate sets* matched on the measured off-target hit probability under 301's exact
rule, promoter GC, TSS count and density, a mappability proxy and basal expression — with every library-targeted
gene, the construct's own target and its own spacers' genes excluded. Attainable floor p = **1.0e-5**, three orders
below α (305's floor was 0.0696; that is the defect this design fixes).

**The estimate.** `S` = **+0.786** standardised units per construct, construct-level bootstrap 95 %
**[+0.387, +1.166]** (n = 34); joint-null p **1.0e-5** (the floor); Δ = **−0.0790 log1p-CP50k [−0.1366, −0.0310]**,
i.e. candidate genes sit ≈ **7.6 % below** their matched non-candidates. It survives two structurally different
nulls (real-spacer resampling p = 0.0015; wrong-recipient construct permutation p = 2.0e-4) **and the design's best
idea — a partner-construct specificity test: the identical candidate genes read at wrong-recipient level (+0.041
[−0.295, +0.369]) in the cells of the partner construct that knocks down the SAME intended gene with DIFFERENT
protospacers.** So the signal follows the **protospacer**, not the intended target.

**Why it is nonetheless not a PASS.** Both pass-rule conditions were met (`p_joint` 9.9999e-06; max
leave-one-construct-out p 1.99998e-05) but **2 of 7 known-answer rows failed**: KA-4's no-signal row carries a
construct-independent bias of **+0.0683** (rule ≤ 0.050; 5,000 derangements give +0.0717 ± 0.0022 — real, and 9 %
of the effect), and KA-6's matching missed on **one axis of nine**, promoter-window count SMD **+0.2905** (rule
≤ 0.250). The frozen abort clause governs and **the standard was not weakened after the fact** — the discipline
305 failed is the discipline 307 kept. Reported as a bounded estimate, not a verdict.

**What it would and would not establish, stated up front:** even a clean PASS would not establish dCas9-KRAB
off-target binding. Expression is observational; the 34 constructs differ in on-target knockdown (−3.04 … +0.03),
cell number (95 … 3,415), composition and handling; there is no biological replicate of any construct. The
partner-construct test excludes intended-target biology and generic candidate-gene properties, narrowing what is
left to **something carried by the protospacer sequence** — seed-mediated binding OR a construct artefact
correlated with spacer identity, which this data cannot separate.

**A bug worth recording.** 307's spacer calibration found an **Excel-escaping defect in Horlbeck's supplementary
gene labels (`="SYMBOL"`)** that had stopped the own-target exclusion firing and inflated the per-spacer L = 12 hit
rate from **0.107 to 0.597**. It was caught by a known-answer row before use and the corrected pool was used
throughout. **301's published T7 numbers are unaffected** — 307's KA-2 reproduced 301's selection exactly (failing
4/13, working 4/21, NTC 4/31).

**→ 309 launched** (19:20Z): claim C9, the calibration repair — bias-immune scoring for KA-4 and promoter-window
count added to the matched axes for KA-6, at the SAME α and the SAME single verdict cell, with C8's estimate
reported alongside whatever happens. **The bias repair moves the estimate DOWN by ≈ 9 %, so it makes passing
harder**; that asymmetry is stated in the preregistration before running, and a second known-answer failure aborts
again.

### 2026-09-07 19:18 UTC+3 — 308 lands (Opus, 45 min): the |kd|-weighted-donor mechanism is REFUTED at this construction

**308 (|kd|-weighted donors): DOES NOT CLEAR.** Best arm `soft:abs:p1.0:dis` `pds` **+0.00077 [−0.00114, +0.00283]**
3/3 folds, priced **−0.00145 pessimistic vs +0.0089**; **0 of 16 arms clear, 51 cells screened; 64/64 known-answer
rows pass** (all 0.0e+00 bar three at ≤ 5e-12; fold caches byte-identical to 262's).

The mandate's premise was checked and held: the feature is real, build-time observable, defined for 226/290/294 of
300 targets, **absent from the incumbent's base by construction**, **not** a norm proxy (ρ +0.03) and only weakly a
reliability proxy (ρ +0.04…+0.19) — it carries information the existing keys do not. The acting set is large
(**106 of 300** targets disagree between donors), so this is not a power failure.

**It fails on the decisive control and on its own mechanism.** |kd| never beats reliability weighting by a
separated margin and sits *behind* 097 route (a)'s `normrel^−1.5` (`gate:abs − normrel1.5` −0.00579
[−0.01042, −0.00132]). The one 3/3-positive cell is **indistinguishable from randomly re-weighting the same rows**
(entropy-matched random +0.00101 / +0.00054 / +0.00100, two of three LARGER than the arm). And the measured reason:
**on the disagreeing set the higher-|kd| donor's row is not the one closer to the held truth** — dropping the
failing donor moves more targets DOWN than UP in all six (cut × fold) cells, by 1.5–2.4×. The variance cost of
averaging two rows instead of three exceeds the bias benefit of removing a failed construct. `normrel1.5` was the
only positive row in the whole screen (+0.00259, 3/3, 0.55× the gate) and it does not clear either.

**097 route (d)'s distinction was verified numerically** in the installed `cell_eval2` 0.16.0 (KA-C): the scorer
drops each row's own-target coordinate from every member. |kd| as a *prediction target* was already closed; |kd| as
a *weight* was genuinely open, is now measured, and is closed at this construction.

### 2026-09-07 19:45 UTC+3 — 309 lands (Opus, 14 min): claim C9 PASSES the calibration repair, and the repair moved the estimate DOWN

**C9 (aggregate candidate-repression enrichment, construct as unit, calibration-repaired): PASS.**
`S` = **+0.71293 [+0.32419, +1.07351]** standardised units per targeting construct (n = 34), joint-null
**p = 1.99998e-05** (attainable floor 1.0e-5), max leave-one-construct-out p = 8.99991e-05, **9/9 known-answer
rows**. Δ = −0.078205 log1p-CP50k [−0.135970, −0.030017]. **C8's +0.786 [+0.387, +1.166] stands unretracted beside
it.** Peak RSS 1.25 GiB, 20 s of compute.

**The repairs and their direction.** (i) Every construct's score is centred on **its own wrong-recipient
distribution** — the same candidate genes, same matched pools, read in the cells of every other targeting construct
of its split (15/9/7 recipients); this is C8's partner-construct test generalised from one wrong recipient to all.
(ii) Promoter-window count added to the hard strata (10 basal deciles × 5 `pi_off` quintiles × 3 window classes =
150 cells). **The preregistration stated before running that repair (i) works against the effect, and it did:
+0.78797 → +0.71293, −9.5 %**, without narrowing the null (sd 0.17316 vs 0.17173), so p rose from the floor to
2.0e-5 — still 2,500× below α. **No threshold was loosened, then or now.**

**What survived.** 26/34 constructs positive (exact binomial p = 0.0015); the real-spacer resampling null
+0.67242, p = 0.00149 at its floor; and the specificity test unchanged — own **+0.78797** vs the partner construct
that knocks down the SAME gene with DIFFERENT protospacers **+0.04296 [−0.29799, +0.38107]**, difference
**+0.74500 [+0.22896, +1.25766]**, paired Wilcoxon p = 0.00459. **The signal follows the protospacer, not the
intended target.**

**The multiplicity, stated by the experiment itself.** C8 and C9 are **two verdict cells on one hypothesis and one
data set** — not an independent replication. Holm over two attempts needs the smaller p ≤ 0.025, which 2.0e-5
satisfies, but the real multiplicity is analytic: C9's design was chosen after seeing C8's diagnostics and the
repaired estimate was **predictable in advance from C8's own post-hocs** (P5 +0.71399). What C9 buys is that the
estimate is not an artefact of the two calibration defects — nothing more. It does not upgrade an association into
a mechanism, and it does not establish dCas9-KRAB off-target binding.

**Recorded as F16 (`docs/FINDINGS-2026-09.md` v3.0)** with its limits in the body: observational expression, no
biological replicate of any construct, one shared control pool across the three splits, and a partner test that
narrows the cause to something carried by the protospacer without separating seed-mediated binding from a
spacer-correlated construct artefact. The decisive next experiments are wet-lab: a **sequence-mutant control**
(silent seed mutations that abolish the predicted match while preserving the on-target site; ~24 constructs × ~500
cells, one 10x lane, order **$3–6 k**) or a **dCas9 occupancy assay** (**$5–10 k**). **→ 310 launched** to fold
Round 6 into the manuscript.

### 2026-09-07 19:50 UTC+3 — October config path re-verified today, both branches, against 256's rehearsal values

Not a full build (256's is 20 GiB / 22 min and nothing in the build path has changed since it ran on 09-05); the
proportionate check is the driver's `--dry-run` env composition, which is what would silently drift. Both branches
reproduce **256's KA 5 to the digit**, and `git status scripts/` is empty:

| branch | measured today | required (R256 KA 5) |
|---|---|---|
| **A/B/C** (`data/challenge`) | `VCC_CD4_DEFLATE_K` **'5'**; `VCC_CD4_DEFLATE_BASIS` **unset**; `VCC_CD4_SIZING` unset; ALPHA 0.982, CD4_TARGET_F 0.34, H1_UNION_G 0.28; 11 stages | k = 5, no basis |
| **non-A/B/C** (238's ef-bundle, contexts D/E/F — **the October case**) | `VCC_CD4_DEFLATE_BASIS` **'V5_abc'**; `VCC_CD4_SIZING` **'carrying'**; `VCC_CD4_TARGET_G` **0.3443655278649641**; no `VCC_CD4_DEFLATE_K` | basis `V5_abc`, sizing `carrying`, no k |

So the pinned-basis branch — the one that fires on the real October panel, which will not be A/B/C — is wired and
carries the same `cd4_target_g` 256 recorded. Combined with 308's independent evidence that the scoring path is
stable (64/64 known-answer rows at 0.0e+00, fold caches byte-identical to 262's), the environment shows no drift.
**A full ef-bundle rehearsal remains the pre-October item**; it is deferred, not skipped, and should run once
closer to Oct 22 or after any production edit. Scratch deleted; nothing written outside the scratchpad.

### 2026-09-07 20:05 UTC+3 — 310 lands and Draft v9 is republished to the same artifact URL

**310 (Opus, editorial only; no measurement, nothing published by the agent).** Round 6 is now in the manuscript:
**§3(g) is a new results subsection** for F16 carrying C8 (`UNRESOLVABLE`, self-aborted on 2 of its own 7 KA rows)
and C9 (**PASS**, `S` +0.71293 [+0.32419, +1.07351], p 1.99998e-05, 9/9), with the bias repair's **−9.5 %** stated
as declared-in-advance, the partner-construct test as the load-bearing design element, and "two tests of one
hypothesis on one data set — **not an independent replication**" in the body rather than a footnote. **§4.5's T7 is
rewritten**: 301's negative stands (reproduced exactly by 307's KA-2), 305 produced no accepted discovery, and F16
supplies the aggregate measurement — with one sentence on why 301's negative and F16 are not in conflict. §3(h)
carries **308's closure with its mechanism** and **306's reconnaissance negative** about the field. §4.1 records
that F16 postdates both audits and carries **no provenance verdict**, and that the phenomenon itself is published —
what is ours is the construct-as-unit test with this null and this control. `CLAIMS.md` gains C-55/C-56/C-57;
`REFERENCES.md` gains Horlbeck and the Excel-escaping defect. **140 provenance entries** in `work/verify310.log`,
each naming source file, line and a verbatim anchor; **seven candidate readings were considered and rejected.**

**Coordinator check before publishing.** The draft cited 305's attainable Holm floor as **0.0506** in one place and
**0.0696** in another; both are correct and cite different families (305's own 8 *testable* pairs vs the
coordinator's binding 11-test family, since untestable hypotheses cannot be claimed rejected). Both exceed 0.05, so
the conclusion is identical either way — no error. The "local; not republished" marker and the build masthead were
then corrected to **Draft v9** before the build, since they were about to become false.

**Republished to the same URL** (https://claude.ai/code/artifact/b7545115-c967-4dd0-b670-7329a788f142), 1,612 KB,
26 TOC entries, all 8 figures. Ahmet's link is unchanged.

### 2026-09-07 20:27 UTC+3 — 311: the FULL October rehearsal passes end to end, today

Not the config dry-run — the real 11-stage driver on 238's ef-bundle (contexts D/E/F), every `VCC_*` unset so the
driver's own defaults had to supply every arm. **rc 0 on both invocations; `def_postbuild` 41/41 checks pass, 0
FAIL.** Peak RSS 20.05 GiB (build), 14.90 GiB (prep). `config.cd4_deflate_basis` = **V5_abc**, `cd4_deflate_k`
**absent** (rule 10 satisfied), `deflation.mode` **pinned**, `captured_frac_of_centred` **0.1333778923231684** —
256 recorded 0.1334, so the pinned basis reproduces to the digit. **Builder restored to
`3615d20f8b535c…` = intended**, `git status scripts/` empty. The October path is verified on today's environment,
not merely assumed from 09-05.

### 2026-09-07 20:35 UTC+3 — 312: the per-target donor weighting is WIRED, and a slot is going to be spent on it

**Why, stated plainly: the slot policy was wrong.** Only the last entry after Oct 22 is prize-eligible, so every
slot before then is prize-irrelevant AND expires nightly. Holding them for an arm that clears a deliberately
pessimistic gate let four rot (09-06 ×2, 09-07 ×2) and bought nothing. The gate decides what to BUILD; it must not
decide whether to spend an expiring resource. Ahmet pushed on this and he is right.

**The arm.** `VCC_NORMREL_EXP` — per-target donor weights, w_j(t) ∝ (‖row_j(t)‖ / that donor's own panel median)^exp
at exp = **−1.5**, renormalised over the donors carrying t, and the row then rescaled to the uniform-weight row's
norm so the arm is **shape-only** (263's `:nm`; AGENTS 8.1d — signs pay, magnitudes do not). This is 097 route (a)
a6 (+0.0008–0.0009, "the only arm positive on both surfaces") and 308's `normrel1.5` (**+0.00259
[+0.00000, +0.00527], 3/3 folds, priced +0.00104**, the only positive row in a 51-cell screen).

**097's note that the builder already took a per-target weight dict was WRONG** — `weights` is
`Mapping[str, float]`, one scalar per donor for every target. But the wiring is still small, because
`build_v9_pre_whiten_effects` already calls `model.effect` once per target, so a per-target dict needs **no change
to `MultiDonorModel`**. 87 lines in the builder, 3 in the driver's arm registry.

**Known-answer rows, all four PASS** (`experiments/312-normrel-arm/test312.py`, synthetic 3-donor model):
KA-1 the unset path runs; **KA-2 exponent 0 == unset at max|diff| 0.000e+00** — an unset build is bit-identical, so
the default path cannot change; **KA-3b row norms preserved to 2.05e-16** — genuinely shape-only; KA-4 the donor
whose row is 5× its own median is down-weighted 0.452 → 0.059 by hand-computed weights. KA-3b first read
1.36e-07; a float64 control proved the logic **exact** (error 0.0) and the residual purely float32 representation,
so the norm-match was moved to float64 rather than the threshold being loosened.

Registered in `ARM_BOARD_DELTA` and written into the meta `config.normrel_exp` under the same "only when set"
convention, so provenance is on the artifact. Build of `the-cat-sat-on-the-modem` running on the real A/B/C bundle.
**Expected +0.001…+0.004 — it will not move us off 3rd.** The reason to spend the slot is that we have never tested
a per-target weighting on the real scorer: if the board returns ~4× the prediction, our pessimistic pricing is too
harsh and several closed arms reopen.

### 2026-09-07 21:35 UTC+3 — 312 SCORES: +0.0016245, the first board gain since 09-05 — and THREE CORRECTIONS to my own reading (GPT 6 Astra review, all three verified here)

**The entry.** `the cat sat on the modem`, entry **`Q519I2pLWRCW7PAEx5ul`**, published **0.21771888483638538**
against `cd4svd5`'s 0.2160944021991198 — **Δ +0.0016245**. Captured to
`data/board-entries/Q519I2pLWRCW7PAEx5ul.json` before anything can replace it. Preregistration
`644a43d4…` was frozen and committed before the entry existed; the result falls in its
**+0.0010 ≤ Δ < +0.0040** bin ("pricing roughly calibrated; the arm is real and small; not a route to first").
Rank **4** — we scored higher and still lost a place, because we sat static for three days.

| member | incumbent | arm | Δ |
|---|---:|---:|---:|
| `pds` | 0.77809 | 0.79039 | **+0.01230** |
| `mse` | 0.17477 | 0.18013 | +0.00536 |
| `nmae` | 0.11047 | 0.10277 | **−0.00771** |
| `fid` / `reach` / `jac` | | | +0.00088 / −0.00113 / +0.00004 |

**CORRECTION 1 — my "4.7× underpriced" claim was a UNITS ERROR and is withdrawn.** I compared PLOLO's
**raw** `pds` prediction (+0.002589) with the board's **scaled** `pds` delta (+0.012302). The like-for-like
comparison is raw-to-raw: board raw `pds` 0.850643 → 0.856154 = **+0.005511**, i.e. **2.13×** the local
prediction, not 4.75×. And **308's own overall point forecast was +0.00184 against a realised +0.001624** — a
12 % miss. **The pricing was well calibrated; the programme is NOT systematically underpriced, and no closed arm
reopens on this evidence.** The component errors (a `pds` under-prediction and an `nmae` sign error) largely
cancelled in the total.

**CORRECTION 2 — the shipped arm is NOT the arm 308 measured.** In 308's construction PerturbMatch is blended
into the K562 donor *before* weighting and norm-matching. Production (`apply_perturbmatch_k562_replacement`,
builder line ~1550) weights and norm-matches FIRST and then applies the PM correction using the **uniform**
`k562_share` — my patch changed only `build_v9_pre_whiten_effects`. **These operations do not commute.** So this
board point validates **the shipped implementation, not the local calibration**, and the local-vs-board component
discrepancies may be partly this. → **313** launched: a shared-fixture parity test. Today's winner is PRESERVED;
any parity-correct version ships as a SEPARATE arm, never as a silent fix.

**CORRECTION 3 — "shape-only to sixteen decimal places" was an INTERMEDIATE guarantee only.** The 2.05e-16 row-norm
identity holds at the pre-whitening effects matrix, which is where my test measured it. Downstream the delivered
magnitude does move: the meta's `after_all_median_norm` goes **4.229296 → 4.272970 (+1.03 %)** on the same recorded
operator hash (`before_all_median_norm` +1.16 %). Adding PerturbMatch after the norm-match changes norms, and later
normalisation does not preserve every row. **The score gain stands; the MECHANISM claim ("pure direction") does
not, and the `mse` improvement must not be attributed wholly to cosine until the emitted pseudobulk magnitude and
the `mse` correction terms are checked.**

**Also corrected:** exponent −3.0 is not "double the dose" — w ∝ r^−3 is the **square** of the r^−1.5 weight
ratios. The dose-2 build is left running (it is still a valid second point on the curve), but it is relabelled
accordingly, and the better next hypothesis is Astra's: **w ∝ max(1, r)^(−p)**, which suppresses oversized donor
rows WITHOUT up-weighting tiny ones — the current rule does both, and 143 is the precedent for `pds` rising while
signed credit falls.

### 2026-09-07 22:05 UTC+3 — 313 lands: production's per-target weighting is NOT the object 308 measured, and that EXPLAINS the `nmae` sign error

**Verdict: IS NOT.** Max |diff| **2.700** on the pre-whiten effect matrix against a float32 noise floor of
3.39e-07 (a factor of 7.9e+06); median per-target row cosine between the two outputs **0.9857**, 15/20 targets
below 0.99. Measured with **both implementations' real functions** on one shared synthetic fixture covering
PerturbMatch present/absent, partial gene coverage, one and two missing donors, a near-zero row and an oversized
row. Peak RSS 0.21 GiB; **production was not touched** (Astra's instruction, and the coordinator's).

**Three independent differences.** (i) The PerturbMatch correction is applied at the **uniform** donor share
`1/|have|` rather than the arm's own per-target share, and **outside** the norm match — the 312 patch changed only
`build_v9_pre_whiten_effects`, so nothing downstream learned about the per-target weights. (ii) The weights are
computed on a **different vector**: production's norm includes the panel columns and excludes the PerturbMatch
blend; 308's does the opposite (its prepared rows are panel-zeroed and PM-blended). (iii) Pre-existing, not caused
by 312: production applies PM to every official column PM observes, the PLOLO harness to `k562 ∩ PM ∩ non-panel`.

**The consequence that matters, and it settles CORRECTION 3.** 308's arm is **exactly shape-only** — delivered
row-norm ratio 1.000000 for every target. **The shipped arm is not**: delivered per-target norm **0.9930–1.3407**
on the full axis, **0.8062–2.1549** on the scored non-panel columns, whole-matrix norm **+6.00 %**. We
preregistered `VCC_NORMREL_EXP=-1.5` as shape-only with a predicted `nmae` leg of **−0.00095 (a benefit)** and the
board returned **−0.007705**. **An amplitude-changing object losing `nmae` is exactly what the shipped form
predicts** (AGENTS 8.1d). The one component we got wrong is now explained mechanistically, not hand-waved.

**What stands and what does not.** `Q519I2pLWRCW7PAEx5ul` **+0.0016245 stands** — it is a real board gain by the
real scorer. What does **not** stand is reading it as a validation of 308's calibration: it validates the shipped
implementation. Today's dose curve (−1.5 and −3.0) remains internally valid because both points share that one
implementation.

**The parity-correct arm is SPECIFIED, NOT APPLIED**: `normrel15-pm`, env `VCC_NORMREL_PM=blend` (default unset =
today's behaviour, so an unset build stays byte-identical), one function to change
(`_normrel_pre_whiten_effects`), the blend folded in **before** the norm match so the row is exactly shape-only,
two mechanical supporting edits, and an explicit instruction **not** to touch
`apply_perturbmatch_k562_replacement` (it is also the incumbent's PM path; changing it would move `cd4svd5` and
every other arm). Two second variables frozen rather than swept in the same slot (D5 the weight support, D7 the
column policy). Three known-answer rows before any build, including "unset reproduces `312-normrel15` byte for
byte" and "`VCC_NORMREL_PM=blend` with exponent 0 reproduces `cd4svd5` byte for byte".

**Why it is the next arm.** On the fixture the parity variant restores a delivered row norm of **1.0000000** for
every target while differing from the shipped arm at max |diff| 1.7297 (median row cosine 0.9982) — a real
re-shaping arm, to be priced as one, not a cleanup. If restoring shape-only recovers most of the −0.00771 `nmae`
leg, the same weighting would land near **+0.0029** instead of +0.0016. That is still short of the +0.0089 gate,
but it is the best-motivated arm on the table and it is **tomorrow's first slot**.

### 2026-09-07 22:25 UTC+3 — 312b: the dose curve TURNS OVER, and 314: don't chase PromoterAI — ChromBPNet is the better instrument

**312b (`a goose in the server room`, entry `1n5LFk6PSEHq15UHENIE`): 0.21645114** — Δ vs the −1.5 entry
**−0.0012677**, landing in the preregistered **Δ₂ ≤ −0.0008** bin: *the curve has turned over.* Three points on ONE
shipped implementation: exp 0 (`cd4svd5`) **0.2160944**, exp −1.5 **0.2177189**, exp −3.0 **0.2164511**.

| member | exp 0 | exp −1.5 | exp −3.0 | Δ(−1.5→−3) |
|---|---:|---:|---:|---:|
| `pds` | 0.77809 | 0.79039 | 0.79241 | **+0.00202** |
| `mse` | 0.17477 | 0.18013 | 0.17634 | −0.00379 |
| `nmae` | 0.11047 | 0.10277 | 0.09980 | −0.00296 |
| `reach` | 0.22139 | 0.22026 | 0.21791 | −0.00235 |

**`pds` rises MONOTONICALLY with |exp| while every other member falls past −1.5** — 143's pattern reproduced on
the real scorer, which Astra named as the precedent *before* this point existed. Implied cosine 0.09728 → 0.09973 →
0.09800. **The exponent axis is closed at its peak ≈ −1.5; more dose is not the lever.** Two slots bought a closed
question rather than a rank, which is the trade the standing rule asks for.

**314 — the new leader's method label, run down.** `ViRtUaL cElL builders` took first at **0.240185** with the
public model name `PromoterAI_precomputed`; until today no entry above 0.220 disclosed anything (306). **Treat the
label as a self-report: 0.240185 is the observable, the attribution is their claim.**
- PromoterAI is Illumina's (Jaganathan et al., *Science* 389(6760) 2025, doi 10.1126/science.ads7373). Input is
  20,480 bp of DNA centred on a TSS. **RELEASED BUT GATED**: code is PolyForm Strict 1.0.0 (noncommercial, and it
  forbids derivative works and redistribution), **no weights in the repo** (the PyPI wheel is 13,351 bytes, code
  only); weights and the precomputed table exist but come only via an Adobe Sign licence and a human emailing a
  link. Precomputed table is hg38 GENCODE v39, TSS ±500 bp, **SNVs only**, per variant × per transcript.
- **The decisive fact, read out of `architecture.py`, not argued:** the released score is
  `tanh(mean over axis=(1,2) of (f(alt) − f(ref)))` — the mean **collapses the track axis**, so all 498 ENCODE
  tracks average into one number. **The published score is cell-type-blind by construction.** The steelman is real
  (CRISPRi efficacy tracks CAGE-TSS proximity and accessibility, which is what it was trained on, so a per-gene
  aggregate is a plausible "promoter disruptability" prior) but it supplies **nothing for held-out-cell-line
  transfer — the actual difficulty** — nothing for the ~17,999 non-target coordinates, and it is confounded with
  free CpG/GC/DNase features and likely redundant with knockdown depth we already observe.
- **The reframe: PromoterAI is the weakest and most encumbered member of its own family.** **ChromBPNet is strictly
  better here — MIT, downloadable today, 1,512 ENCODE4 models over 408 biosamples, base-resolution promoter
  accessibility PER CELL TYPE.** Borzoi is Apache-2.0 with predicted expression; scooby is MIT and
  single-cell-resolved; AlphaGenome weights went downloadable 2026-01-28 (click-through, non-commercial).

**Why this matters against 306's clean negative.** 306 established there is no new *perturbation* dataset to buy.
This is a different category: a **regulatory model that supplies per-cell-type chromatin features for the HELD
line**, which no perturbation screen provides. F5/289 measured that basal EXPRESSION does not predict the
line-specific residual and concluded a stage 2 needs more LINES; per-cell-type accessibility is a richer feature
than the two expression PCs that failed, and 408 biosamples is a far larger line axis than the 11–12 we had.
**This is the first competition-track direction all day that is not already closed at its scope.**

### 2026-09-07 22:50 UTC+3 — 315: NO-GO on per-cell-type promoter accessibility, and the NO-GO is more useful than a GO would have been

**315 (per-cell-type promoter accessibility for the held line): NO-GO** — identity **A = JURKAT, B = HELA,
C = CAL33** (unanimous under 8/8 feature-rank variants, but the runner-up is a same-lineage relative in all three
and the consensus margin is only 1.0–1.5 of 25 — a **ranked match with moderate confidence, NOT a fact**);
coverage **1/12** Broad-16 lines by ChromBPNet model *and* by measured DNase/ATAC alike; Gate 3 **could not be
run**, and the count is the answer.

**The two reasons, both structural, neither fixable by a better estimator.**
1. **ChromBPNet is one model per ASSAYED biosample; it does not generalise to an unassayed cell type.** For a
   genuinely held-out line the feature is not weak — it is **undefined**. Held context C (CAL33) has zero
   ChromBPNet models, zero DNase-seq and zero ATAC-seq on the ENCODE portal today; A and B have them. Coverage of
   the held contexts is **2/3**, and October's D/E/F identities are unknown before the round, so 3/3 cannot be
   assumed. **This is the same defect that makes the leader's PromoterAI label unusable (314: cell-type-blind by
   construction) — arriving from the opposite direction.**
2. **The only surface where the claim could be TESTED has 1 of 12 lines.** A leave-one-line-out over one line has
   no training row and no held row. Of our own donors only K562, HCT116 and Jurkat carry ENCODE accessibility.
   289's harness was verified to be ready (KA-1 reproduces its frozen LOLO exactly) — **only the feature is
   missing.**

**And yet the biology is real, which is why this closure is worth having.** On a 4× larger panel than Broad-16 —
49 cell lines with both ENCODE DNase-seq and a CCLE profile — basal expression explains only **21.0 %
[19.5, 26.7]** of the line-specific variation in promoter accessibility, and the expression-**orthogonal**
remainder reproduces between two independent DNase-seq experiments of the same line at **+0.444 [+0.335, +0.552]**
against a wrong-recipient control of **−0.025**, **17/18 lines positive**. So promoter chromatin does carry
line-specific information basal expression does not — F5/289's negative was about *expression*, and that stands
unchallenged — but the information cannot be delivered to an unassayed held line. **The constraint is
availability, not biology, and it applies to every entrant equally, including the leader.**

**Verification note (coordinator).** 315's self-reported wall clock ("≈ 2 h 45 min") is **wrong** — the run
started 19:25:54Z and its last event is 19:46:25Z, about 20 minutes. The *substance* is verified: 69 ENCODE
narrowPeak files (125 MB) are on disk as claimed, and KA-1's reproduction of 289 (`S1.g` cos_res +0.010972
[+0.006007, +0.016234], Δcos_full −0.058756) matches R289's published +0.0110 [+0.0060, +0.0162] / −0.0588 to the
rounding. **A reporting defect, not a data defect** — but agents' self-reported timings are not evidence and are
not to be quoted.

### 2026-09-07 23:25 UTC+3 — The new leader is NOT winning on transfer cosine, and that contradicts §4.1c's "there is ONE route"

`ViRtUaL cElL builders` took first at **0.240185** (`PromoterAI_precomputed`, 12 submissions, 17:06Z). Their member
profile is **the opposite shape to `cqawesome`'s and to ours**, and the inversion is **in domain** for them
(`score_mse` +0.07252 > 0, §4.1c's validity guard satisfied):

| raw metric | us | LEADER | |
|---|---:|---:|---|
| `expr_mse_unbiased_capped_norm` | 0.81671 | **0.99153** | higher = worse; they sit near the α=1 baseline |
| `de_wilcoxon_lfc_nmae` | 0.93795 | **0.87166** | lower = better |
| `de_wilcoxon_direction_reach_raw` | 0.27428 | **0.37115** | +35 % over us |
| `pds_cosine` | 0.85615 | 0.86832 | |

**Scaled, leader − us:** `mse` **−0.10761**, `nmae` **+0.11257**, `reach` **+0.10990**, `pds` +0.02657,
`fid` +0.01386, `jac` −0.02050 → total **+0.02247**. **Their implied transfer cosine is 0.01546 against our
0.09973 — six times WORSE — and they beat us by +0.0225.**

**The reframe, and it is uncomfortable.** §4.1c concluded "**there is ONE route, not two: transfer cosine**", and
227 priced the `nmae` credit at +5.0 % relative aggregate cosine on h1-test / +12.3 % on PLOLO against a ±0.2 %
gate — "failed by 25–60×". **That pricing assumed we must KEEP the cosine.** This entry paid exactly that price:
it gave the cosine away (0.0155) and bought `nmae` + `reach` (+0.223 between them) for a `mse` cost of −0.108.
**Net +0.0225. The trade our own instrument said was catastrophic is what is currently in first place.**

**Status: HYPOTHESIS, not a finding.** It rests on ONE public entry, a self-reported model name, and a
cross-team inference — and §4.1c contains a **retracted** amendment (04:35 → 05:10 on 2026-09-04) from exactly this
kind of reasoning, where a clamp manufactured a phantom anti-correlation. What is solid: the six member columns
and the raw metrics are the board's own numbers, and the inversion guard is satisfied. What is NOT established:
that OUR pipeline can move along that tradeoff, or that their operating point is reachable from ours.

**What it changes tonight:** nothing shipped. **What it changes tomorrow:** after 316's parity arm, the next
feasibility question is no longer "more dose" or "a better weight shape" — it is **"can our own pipeline trade
`mse` for `nmae` + `reach`, and where does the total peak?"** That is a sweep on the emission/calibration axis we
stopped exploring when §4.1c declared the single route, and 227's own numbers are the map of what it costs. Queue
it as 317 with the operating point, not the arm, as the variable.

### 2026-09-08 01:20 UTC+3 — CORRECTION 4 (mine, self-caught): "the leader wins with 6× worse transfer cosine" is NOT supported, and 197's inversion has a second, unwritten validity limit

Two hours ago I logged that `ViRtUaL cElL builders` win at 0.240185 with an implied transfer cosine of **0.01546**
against our 0.09973 — six times worse — and that this contradicts §4.1c's "there is ONE route". **The direction
half of that claim does not survive its own arithmetic.**

| | us | leader |
|---|---:|---:|
| **raw `pds_cosine`** — the scorer's OWN per-target statistic | 0.85706 | **0.86832** |
| implied cosine from 197's `mse` inversion | 0.09800 | **0.01546** |

**Both cannot be true.** A team whose per-target rank statistic is *better* than ours does not have a transfer
direction six times worse. The scorer's own measurement wins; **the inversion is the weaker instrument and it is
what breaks.**

**The unwritten limit.** §4.1c's guard says "only invert a row whose `score_mse` is strictly positive", and theirs
is +0.07252, so it is **formally in domain** — which is exactly why I trusted it. But 197's relation
`cos ≥ (1.023611 − raw_mse)/2.074495` was fitted **within our own emission family**. Raw `mse` is a function of
emission scale and dispersion as well as direction, so a team emitting differently can carry a high raw `mse`
(theirs 0.99153) with a perfectly good direction. **A positive `score_mse` makes the inversion in-domain; it does
not make it applicable ACROSS emission families.** That limit is now written down; it was not before, and §4.1c
already carries one retraction (2026-09-04 04:35 → 05:10) from misapplying this same instrument.

**What survives, and what 317 is actually testing.** The member columns are the board's own numbers and stand:
leader − us is `mse` **−0.10761**, `nmae` **+0.11257**, `reach` **+0.10990**, `pds` +0.02657 → **+0.02247**. So the
*shape* of their advantage is real — they pay `mse` and buy `nmae` + `reach`. What is **withdrawn** is the causal
story that they achieved it by *giving up direction*. Their direction is fine. **The trade is emission/calibration
against `mse`, not direction against the DE members** — which makes 317's sweep (overdispersion, spread,
tail-spread, keep-top/spread, amplitude, all six members scored) exactly the right experiment, and makes it MORE
likely to find something, since it need not sacrifice the direction we are good at.

**Two facts that bound the enthusiasm.** (i) The DE-traded shape is **not a cluster**: 1 of the top 25 and 9 of
780 published entries have it, and the leader is the best of those 9 — so it is close to a one-off at the top.
(ii) Their own trajectory is 12 submissions of violent variance — −0.96554 (09-02), 0.05051 (09-06 03:14), 0.11352
(09-07 03:53, with `fid`+`jac` summing ≈ −0.80), then 0.24018 (09-07 18:53). **They are sweeping an unstable
surface, not executing a clean method**, which is consistent with a calibration/emission axis and inconsistent
with a model advantage.

### 2026-09-08 01:50 UTC+3 — The leader is NOT in the "never α-shrank" cluster; they sit just INSIDE the `mse` clamp, and that gives 317 a target region

Profiling all 9 published entries with the DE-traded shape (better `nmae` AND `reach`, worse `mse` than us), on
**raw** members only — no inversion, after Correction 4:

| team | raw `mse` | raw `nmae` | raw `reach` | raw `pds` | overall |
|---|---:|---:|---:|---:|---:|
| **ViRtUaL cElL builders** | **0.9915** | **0.8717** | **0.3711** | 0.8683 | **0.24018** |
| nailong king | 3.8828 | 0.9005 | 0.2771 | 0.8203 | 0.18540 |
| brocade | 3.8832 | 0.9130 | 0.2763 | 0.8203 | 0.18151 |
| The_Way | 7.3742 | 0.9313 | 0.2861 | 0.8130 | 0.17738 |
| luka | 3.8673 | 0.9222 | 0.2908 | 0.8162 | 0.17324 |
| CYH | 3.5197 | 0.9070 | 0.2769 | 0.7858 | 0.17017 |
| Evidoa.com | 2.0993 | 0.9229 | 0.2811 | 0.7948 | 0.16363 |
| kriagr | 2.6532 | 0.8936 | 0.2768 | 0.7578 | 0.15682 |
| 王爱媛 | 1.8126 | 0.8894 | 0.2762 | 0.7917 | 0.13735 |
| **us** | 0.8203 | 0.9396 | 0.2721 | 0.8571 | 0.21645 |

**Eight of the nine are §4.1c's known cluster** — *"teams who never α-shrank"* — raw `mse` **1.81–7.37**, far past
the **1.02361** clamp, so their `score_mse` is 0 and they land at 0.137–0.185. **The leader is not one of them.**
Their raw `mse` is **0.9915**, i.e. **just inside the clamp** (3 % under it), so they still collect
`score_mse` +0.0725 — and their raw `reach` **0.3711** is **28 % above the best of the other eight** and 36 % above
ours, with the best raw `nmae` on the board.

**What that gives 317, concretely.** They are not trading direction away and they are not a broken submission: they
are sitting at an **operating point near the `mse` clamp edge** where `reach` and `nmae` are maximal but `mse`
still scores. **Our raw `mse` is 0.8203 — well inside.** There is room between 0.82 and ~0.99 that we have never
explored, because every amplitude/emission sweep we ran (211, 216, 034, 094) was judged on cosine or on `pds`, not
on the six-member composite. **317's sweep should be read against a target region of raw `mse` ≈ 0.95–1.00, and
the composite — not any single member — is the criterion.** Past the clamp is a cliff, not a slope: the other
eight prove it costs 0.03–0.08 overall.

**Still bounded:** one entry defines the good point, their 12 submissions swing wildly (−0.966 → 0.051 → 0.114 →
0.240), and `PromoterAI_precomputed` is a self-report. This is a **target region for our own sweep**, not a recipe.

### 2026-09-08 03:25 UTC+3 — 317 lands: the operating-point trade DOES NOT EXIST as a competition lever, and §4.1c is re-established on far stronger evidence

**317: DOES NOT EXIST.** Instrument: the real `cell-eval2` 0.16.0 (`--preset vcc2026`) on eleven disposable EC2
workers, **$17.18 total, $1.56 per six-member score, all eleven terminated, nothing left running**. **KA-1 PASS
bit-for-bit** — 221's recorded `rc_callset` six-member row re-scored on a fresh worker, five members exact, max
|diff| 4.44e-16 — and 10/10 workers gated their rebuilt `base` sha256 against 221's board-anchored artifact.

**The emission axis is empty.** Best single point `ts070` **+0.000131** overall = **7 %** of the +0.0019
shape-only gate and **1.4× a re-seed**; its only moving member's interval includes zero; the five-point *oracle*
(each member's best, unreachable) reaches 23 %. **094's "growing the call set cannot pay" is confirmed on the real
scorer, six members at once, for the first time** — `kt16000` is the sweep's worst arm at −0.007338, carried by
`jac` −0.009270.

**The amplitude axis: the trade is REAL and runs 28.6–40.8 : 1 AGAINST us.** `amp120` buys +0.001340 across the
five non-`mse` members and pays **−0.038293** on `mse`. **Even counting `mse` FREE** — the most generous possible
reading, granting the premise that a team at the clamp pays nothing more — the best single point reaches
**+0.002374 (27 % of the +0.0089 amplitude gate)** and the five-point oracle **+0.005190 (58 %)**. **The whole
axis cannot clear the gate even for free.** And `mse` is not free: our live `score_mse` is 0.18013, so the entire
pot on that side is +0.03002 overall, every point of which the measured exchange rate spends at a loss.

**The board's own history agrees, over 2,754 distinct submissions / 1,961 consecutive same-team pairs.**
Direction-stable pairs (|Δ raw `pds`| < 0.002, so most model rewrites are removed): **−0.17583 overall points per
+1.0 raw `mse` given up** (median). Teams that **gave up** `mse`: mean Δ`score_avg` **−0.01475**, only **33.9 %
improved**. Teams that **bought** `mse`: **+0.01538**, **83.2 % improved**. **Selling `mse` is what loses on this
board; buying it is what wins.**

**So the leader's advantage is NOT a calibration choice we can copy.** Their raw `reach` 0.3711 against our 0.2721
does not come from an emission or amplitude setting — we swept both and the surface is flat or adverse. It comes
from a **different prediction**. That closes the operating-point hypothesis (my 21:46Z framing) cleanly, and with
Correction 4 it closes the whole "second route" line: there was no route, and the appearance of one was two
instrument errors of mine stacked.

**Decision-tree consequence:** every "trade EXISTS" branch in `docs/NEXT-SLOT-DECISION.md` is now dead. Slot 2 goes
to **Astra's `max(1, r)^(−p)` cap variant** if 316 pays, and to a re-think if 316 does not.

### 2026-09-08 03:35 UTC+3 — 316 SCORES 0.2169733: the preregistered IMPORTANT NEGATIVE fires, and the decomposition is better than a win

**Entry `2K5ppxkowMSj2ugoZqYb`** (`a duck reviews the diff`), first slot of 09-08, submitted at the reset,
captured immediately. **0.21697328748387812**, Δ **−0.0007456** against the shipped −1.5 entry — inside the
preregistered **Δ ≤ −0.0002** bin: *the parity-correct form is WORSE.*

**The decomposition, and it is the useful part.** Against `cd4svd5` (0.2160944):

| arm | overall | vs incumbent |
|---|---:|---:|
| shipped (weighting **+** the accidental per-target norm distortion) | 0.2177189 | **+0.0016245** |
| parity (weighting **alone**, exactly shape-only) | 0.2169733 | **+0.0008789** |
| ⇒ the distortion's own contribution | | **+0.0007456** |

**Both arms beat the incumbent.** The per-target weighting mechanism is **real and worth +0.00088**; the accidental
distortion adds a further **+0.00075** on top. My 22:05Z reading — "if restoring shape-only recovers the `nmae`
leg, the same weighting lands near +0.0029" — was **half right and the half that failed is the informative half.**

**313's mechanism was RIGHT about `nmae`, and outweighed.** Parity recovers `nmae` **+0.00339** of the −0.00771 the
shipped arm lost, exactly as an amplitude-changing object losing `nmae` predicts. But the same distortion was
**buying** `mse` **+0.00399** and `pds` **+0.00225** — together more than it cost.

**Two independent experiments converged tonight.** 317 measured, over 1,961 consecutive same-team board pairs,
that **buying `mse` is what wins here** (mean +0.01538, **83.2 %** of teams improve; selling it: −0.01475, 33.9 %).
The distortion was buying exactly that. It was not a bug getting lucky — it moved the direction the board pays for.

**The tension that is now the live question.** 317 swept **uniform** amplitude and found ×1.20 *pays* `mse`
−0.038293 (worse). Our distortion raised the whole-matrix norm ~6 % and *bought* `mse` +0.004. **Opposite signs, so
these are different operations:** a uniform scale versus a **per-target redistribution** of magnitude (delivered
row norms 0.806–2.155 where parity holds them at 1.000). AGENTS 8.1d's "magnitudes do not pay" was established for
the uniform axis. **A per-target magnitude rule is a different object and has never been tested.**

**Decision, per `docs/NEXT-SLOT-DECISION.md`'s own STOP branch: no third weighting variant ships.** The displayed
entry is now the weaker of our two arms; **it is NOT re-submitted to restore rank** (the standing rule). Slot 2 of
09-08 stays open pending a feasibility test of the per-target magnitude hypothesis, which is post-hoc on one point
and confounded (the shipped arm also applies PerturbMatch at the uniform share), so it needs isolation before it
earns a slot.

**Also caught, and worth recording as a near-miss.** 317's `amp085` shows +0.014556 on h1-test. **It is not an
opportunity** and 317 says so: h1-test must be scored with the H1 union OFF (oracle-leaked, 85/85 targets), so its
prediction is materially weaker than the board's (raw `mse` 1.016 vs the board's 0.81671) and a weaker direction
wants a smaller amplitude. **h1-test cannot locate the board's amplitude optimum; 187 and 211 locate it, at 0.982,
where we ship.** That is the third instrument-domain trap of the day, and this time the experiment caught it
before I did.

### 2026-09-08 04:10 UTC+3 — Fable session (269c6ee5) standing down to read-only. Its 30-min cron heartbeat (d9d03f57) is deleted so it cannot launch or edit while the Opus 5 session coordinates; its 3-hour watcher (completion events + board line) stays armed and only reports. No launches, no edits to READ-FIRST/FINDINGS/ledger from this session since 13:49Z 09-07. Live board at 00:36Z: leader cqawesome 0.24661, second ViRtUaL cElL builders 0.24018, Aginglab 0.22652, our live entry 2K5ppxkowMSj2ugoZqYb 0.21697 rank 4 of 781; 318 running (Opus session).

### 2026-09-08 04:45 UTC+3 — Two operational notes: a parallel session, and the leader's +0.00725 is one member only

**A second session was writing this repo.** A **Fable 5.1** session ran in parallel through the night and committed
some of the same artifacts this session was interpreting (`c1860fa`, 317's `RESULTS.md`, two minutes before this
session's `85e203b` LOG entry on the same experiment). It then **stood down to read-only** (`8272bc4`) and deleted
its cron heartbeat, explicitly to let the Opus session coordinate. **The real risk of two coordinators — both
spending slots — did NOT occur:** exactly one slot was spent on 09-08 (`2K5ppxkowMSj2ugoZqYb` at 00:01:08,
`submission_count` 34), and one remains. **Rule for any session resuming here:** check
`ls -t ~/.claude/projects/-Users-ahmet-Code-virtual-cell-challenge-2026/*.jsonl | head -3` and `git log -3` for a
foreign `Co-Authored-By` before launching, editing shared docs, or arming a submit; if another session wrote inside
2 h, stay read-only unless Ahmet is addressing you directly.

**`cqawesome` retook first: 0.23936 → 0.24661 (+0.00725), and it is ONE member.** `m0907h` → `m0908a`:

| member | before | after | Δ |
|---|---:|---:|---:|
| `fid` | −0.05836 | −0.01833 | **+0.04003** |
| `jac` | −0.00273 | −0.00002 | +0.00271 |
| `reach` | 0.23263 | 0.23358 | +0.00095 |
| `pds`, `mse` | — | — | **exactly 0.00000** |

Raw `pds_cosine` 0.88350 and raw `mse` 0.72339 are **unchanged to five decimals**. **Their prediction did not
move; only the DE call set / emission did.** That is the channel this repo treats as closed (094, 231, 232, 234).

**It does not reopen for us, and the number says why.** They moved `fid` from −0.05836 to −0.01833; **ours is
already −0.01686 — still better than their new value.** They closed their own deficit; they did not find a lever
above where we sit. A team with a `fid` deficit can collect +0.0072 from the ordering channel; we have no such
deficit to collect. Gap to first is now **+0.02964**.

### 2026-09-08 05:50 UTC+3 — October readiness: all three of 233's items are handled, verified tonight; and the endgame reframed

While 318's scorers run, I checked whether 233 §7's three October-specific items — the ones "that pay only in
October" — were ever executed. **All three are, and I verified two of them tonight directly:**

1. **The October-four sources: registered, correctly never selected.** `make_submission_pm_cd4.py:971` carries
   them ("*THE OCTOBER FOUR. Registered, never selected.*"). 233 priced them at +5.5–17.7 % relative cosine each
   in October against zero today — but **238 re-priced them and the asset largely evaporates**: 70.7–90.4 % of the
   cosine gain and **10 of 10 negative `pds`** come from **common-essential targets**, and Arc's panel contains
   **zero essential genes** (8.1g). On an essential-free panel the robust value is **16 %** (+0.000257 vs
   +0.001607). So "registered, never selected" is the right state, not an oversight.
2. **Carrying-row sizing** — verified tonight in the driver dry-run: on a non-A/B/C panel `VCC_CD4_SIZING` is
   `carrying` and `VCC_CD4_TARGET_G` is 0.3443655278649641, matching 256's rehearsal to the digit. Without it the
   October sources deliver identically zero.
3. **The `axis_sha256` guard and the operator refit** — 235's guard is in place (A/B/C axis `edf6ba88…`), and 311's
   **full 11-stage rehearsal passed tonight, 41/41 postbuild, builder restored to its intended hash.**

**The endgame, stated with numbers rather than mood.** The gap to first is **+0.02964**. Our best single arm ever
is **+0.0016** and a typical paying arm is **+0.0009**; at two slots a day to Nov 5 with a generous hit rate, the
whole remaining slot budget is worth single-digit thousandths. **Knob-turning cannot close this gap — that is
arithmetic, not pessimism.**

**But the board we are ranked on is not the board that pays.** Oct 22 releases a new test panel (D/E/F, new cell
lines) and 233 measured what a disjoint panel costs: **93 % of cosine retained but only 64 % of `pds` and 60 % of
ordering**, projecting 0.2118 → **0.137–0.119**, *and the loss is carried by the TARGETS, not by coverage*. That is
an explicit **lower bound**, because a random draw's transferability ceiling is 1.3–2.0× below Arc's and Arc's
selection protocol applies again in October. **Everyone drops; relative standing depends on whose construction is
panel-robust, and 233's own guidance was "buy breadth, not a new mechanism" — coverage explains 78 % of the cosine
loss.** The one thing no local instrument can touch is the **context** change (new cell lines), which is the same
constraint 315 hit from the chromatin side.

**Consequence for allocation:** October readiness needs no further work — it is done and re-verified. The marginal
day is better spent on breadth/robustness and on the scientific track than on chasing +0.001 arms, and the slots
remain worth spending only because they are free and they buy information.

### 2026-09-08 07:10 UTC+3 — 318: NO-GO, not separable from its own shuffles — and CORRECTION 5 (mine): I quoted FIXTURE numbers as delivered production values

**318: NO-GO, and NOT SEPARABLE FROM ITS SHUFFLED CONTROL.** The replayed per-target magnitude profile, isolated
at fixed direction and fixed whole-matrix norm, buys **+0.000093 overall = 1.0 % of the +0.0089 gate**, at
**z = −0.06** against three entropy-matched shuffles of its own multipliers — **the best shuffle (+0.000543) beats
the real profile.** Six times the measured redistribution reads **−0.001101**: the axis is negative and steepens
with dose. Instrument: the real `cell-eval2` 0.16.0 on six disposable workers; preregistration hashed twice, both
digests **before** the first worker launch. **The per-target magnitude axis is CLOSED.**

**CORRECTION 5 — the numbers I put in the LOG, in AGENTS.md §6.2b and in 318's own mandate were FIXTURE numbers.**
"Delivered per-target row norms **0.806–2.155**, whole matrix **+6.0 %**" came from `316/work/test316.json`, which
measures on **313's synthetic fixture — 20 targets × 60 genes**, a third of the columns "panel" against
production's 1.62 %. 318 measured the same quantity on the **real** thing: the delivered log-effect of every
(context, target) block of the two 12.6 GB predictions that produced the board points, 900 blocks, all 300 panel
targets × 3 contexts:

| | fixture (what I quoted) | **real panel (318)** |
|---|---|---|
| ratio range | 0.806 – 2.155 | **0.9321 – 1.2349** |
| whole-matrix norm | **1.0600** | **1.0085** |
| sd(log m) | — | **0.032885** |

**The distortion is 6× smaller in spread and 7× smaller in mean than I stated.** AGENTS.md §6.2b is corrected in
this commit. **This is the same failure as Corrections 1–4: a number measured in one domain quoted as if it were
another.** Fifth of the day, third that is mine alone.

**And "exactly shape-only" is wrong at the delivered level too.** 318 measured that the **parity arm redistributes
delivered magnitude MORE** (sd(log) **0.0485**) **than the distortion it was built to remove** (0.0329) — because
the norm match is enforced one stage before an operator whose per-row gain spans **0.591–1.286**. So 316's
`4.44e-16` identity is real *at the pre-whitening matrix* and does not survive to the emitted prediction.
**Whenever an invariant is claimed here, the stage it was measured at must be named, and every stage after it can
undo it.**

**The premise itself also fails:** **half of the +0.00075 the distortion earned is `pds`**, and per-target
magnitude **cannot move `pds`** (it is scale-invariant per target). So the +0.00075 was never a magnitude effect
to isolate.

**What survives.** The profile IS structured — reproducible across the three contexts at Pearson **0.9992** with
ICC ≈ 1.000, and dose-coherent between exp −1.5 and −3.0 at Spearman **+0.9931**. It is a deterministic function of
the target. **It just does not pay, and it does not beat a shuffle of itself.**

**Slot 2 of 09-08 therefore goes to Astra's `max(1, r)^(−p)` cap variant** — the one weighting hypothesis still
untested, which separates "suppress oversized donor rows" from "up-weight tiny ones". Patch staged at
`experiments/312-normrel-arm/apply_cap_patch.py`; the builder is now free.

### 2026-09-08 08:35 UTC+3 — 319 WINS (new best 0.2186597) and 320 deflates F16 to a re-measurement

**319 (`the parrot approved my expenses`, entry `xuDelEik7zGWMXdNu0R4`): 0.21865973944346428 — our BEST arm ever.**
Δ **+0.0009409** over the shipped `plain` arm, landing in the preregistered **Δ ≥ +0.0005** bin: *suppression alone
beats the combined rule.* Against `cd4svd5` it is **+0.0025653**, **1.58×** the previous best (+0.0016245).

| member vs `plain` | Δ |
|---|---:|
| `nmae` | **+0.00475** |
| `reach` | **+0.00193** |
| `mse` | +0.00088 |
| `pds` | −0.00208 |

**Astra's hypothesis is confirmed and it is specific:** `r^exp` conflates *suppressing* oversized donor rows with
*up-weighting* undersized ones; clamping `r ← max(1, r)` keeps only the suppression and gains, giving up a little
`pds` to buy `nmae` and `reach`. **The up-weighting half was amplifying small noisy rows** — which is what 312b's
turnover past −1.5 was. Implied cosine 0.09728 → 0.09973 → **0.10014**. **The preregistered consequence: the dose
axis REOPENS for the cap form only**, and the next arm pushes `cap` to a larger |exp| where `plain` turned over.
Credit for the idea is the review's, not this session's; I was sweeping the dose harder when the rule needed
splitting.

**320 — the novelty audit F16 never had: NEW 0 / EXTENDED 1 / KNOWN 4 / CONTRADICTED 0.** Reading and citation
only, 6/6 of its own known-answer rows reproducing 301/307/309's published numbers.

| item | verdict | prior report |
|---|---|---|
| seed-mediated CRISPRi off-target transcriptional effects exist | **KNOWN** | **Rohatgi 2024** (*Cell Genomics* 4(11):100693) — screens + RNA-seq **+ dCas9 ChIP-seq occupancy**, 98–1,141 reproducible peaks per guide, **in our constructs' own parent library**, two years earlier; Hartman 2026 across three genome-wide Perturb-seq sets |
| a 12-nt seed + NGG promoter match predicts repression | **KNOWN** | Rohatgi's and Hartman's rule, window and 12-nt threshold — **both evaluated it predictively and both report low precision** |
| construct-as-unit + a joint null over covariate-matched candidate SETS | **EXTENDED** | a recombination of Cazares 2024's aggregate seed-shift test and Davis 2023's matched-null method; new only in using a *measured* 209,070-spacer propensity as the covariate |
| the partner-construct specificity control | **KNOWN** | **Jackson 2003** — 23 years old, the founding observation of the field, and already used for exactly this purpose on CRISPRi Perturb-seq by Hartman and as the discovery device by Rohatgi |
| aggregate finds what per-pair cannot | **KNOWN** | Goeman & Bühlmann 2007; the standard motivation for aggregation testing |

**"Nothing biological survives as new. F16 is a re-measurement of a known phenomenon on a new sample — the
challenge's own 2025 library. That is worth having, and it is worth saying that it is what it is."** The audit also
finds a **statistical** overreach of ours: the headline p is a **gene-sampling** p, not a construct-level one,
which undercuts the "construct as the experimental unit" framing, and it quotes two `PAPER.md` sentences that
claim the partner-construct control as ours. Both must be weakened.

**The practical consequence, and it is a correction I owe Ahmet directly.** I recommended a **$3–6 k
sequence-mutant control** as the decisive next experiment. **Rohatgi 2024 already demonstrated the mechanism at the
level of dCas9 BINDING, with occupancy data, on the same parent library.** The mechanism is not open; our proposed
experiment would largely re-demonstrate it. **That recommendation is withdrawn pending a re-read of Rohatgi's
design** — it is no longer obvious that any wet-lab spend here buys something the literature lacks.

### 2026-09-08 09:05 UTC+3 — Draft v10 republished: every overclaim 320 named is weaker or gone

**321 (editing only, nothing published by the agent).** **23 passages changed** across `PAPER.md` (16),
`CLAIMS.md` (6) and `REFERENCES.md`, plus the build masthead — and the agent flagged the masthead itself as a
**fourth file outside the three it was authorised to edit**, which is the disclosure the rule exists for.
Its own summary: ***"Not one number was added that strengthens a claim."***

- **F16 is now presented as a RE-MEASUREMENT** with its provenance table in the body (NEW 0 / EXTENDED 1 /
  KNOWN 4 / CONTRADICTED 0).
- **The partner-construct control is credited to [JACKSON03] and is no longer claimed as ours** — 7 references in
  `PAPER.md`, using **320's wording rather than the agent's own**, which is how a "correction" avoids quietly
  re-softening into another overclaim.
- **The headline p is stated as a GENE-SAMPLING p in four places**, so the construct-as-unit design is described
  without the p-value being asked to carry it.
- **The $3–6 k sequence-mutant recommendation is withdrawn in the text**, because Rohatgi 2024 already
  demonstrates the mechanism at the level of dCas9 binding in this library's parent collection.
- Six references added with DOIs (Rohatgi 2024, Hartman 2026, Jackson 2003, Cazares 2024, Davis 2023,
  Goeman & Bühlmann 2007); §3(h) carries Round 7's five scored entries including the new best **0.2186597**.

**Republished to the same URL** (https://claude.ai/code/artifact/b7545115-c967-4dd0-b670-7329a788f142), 1,627 KB,
26 TOC entries, all 8 figures. Ahmet's link is unchanged. **FINDINGS is at v3.1** with the same re-grade.

**The day's shape, stated once.** The best competition result of the programme (**+0.0025653**, entry
`xuDelEik7zGWMXdNu0R4`, 0.2186597) came from an external review's hypothesis, not from this session's sweeping;
and the programme's best scientific finding was deflated to a re-measurement by an audit this session
commissioned. **Both were external checks catching what would otherwise have been carried forward** — one into
wasted slots on the wrong axis, one into a manuscript that claimed a 2003 control and a 2024 phenomenon as new.

### 2026-09-08 11:55 UTC+3 — 322: breadth is NOT October insurance — and CORRECTION 6 (mine), I had the refuting number in front of me

**322 (broad-source October insurance): NOT WORTH IT.** `$0` spent, 1 h 15 min of a 5 h box, no cloud, no build.
`kolf2`'s Δ`pds` contribution is **−0.01017 today** and **+0.00000 in October as the bank is built** — 0/3 folds
positive at the amplitude-matched dose on **all three panels in both constructions**. Today's cost is ≈ **−0.0038 /
−0.0028 overall** by 233's own route. The cosine half *does* retain (+0.84 % → +1.53/+2.00 % relative, 210 %) —
**but cosine is the channel 233 measured at 93 % retention, i.e. the one that was never in danger.**
**Breadth is insurance for the half that is not at risk.**

**CORRECTION 6 — my mandate's premise was inverted, and the number that refutes it was already in two of our own
documents, one of which I had read hours earlier.** I wrote: *"`kolf2` covers 273 of 300 panel targets — far
broader than `dld1_lowmoi` (67)."* Both halves are true and **they are about different quantities**:

| | `kolf2` | `dld1_lowmoi` | `gse295214_h1` |
|---|---:|---:|---:|
| coverage of TODAY's 300 | **273** | 67 | 64 |
| perturbations the bank contains **at all** | **273** | **4,715** | 2,354 |
| coverage on 233's `real` / `pess` disjoint draw | **0 / 0** | 67 / **73** | 57 / 34 |

**`kolf2` is a targeted screen OF the panel — all 273 of its rows ARE panel rows — so on a disjoint panel its
coverage is exactly zero.** `dld1_lowmoi`, which I named as the *narrow* contrast, is **17.3× broader as a source**
and one of only three banks whose coverage *rises* on a disjoint draw. **`docs/DEF-RUNBOOK.md:342` records
`kolf2` 273 → 0 → −273 in a table, and `docs/READ-FIRST.md` §7 states "`kolf2` 273→0" in a line I read at 20:17Z
last night while checking October readiness.** This is not an instrument-domain subtlety like Corrections 1–5.
**I had the number and did not use it.**

**What the experiment produced anyway, and it is worth more than the arm would have been.** Over all 21 banks on
disk, the predictor of October coverage is **not** today's panel coverage:

| predictor of coverage on a disjoint draw | ρ (`real`) | ρ (`pess`) |
|---|---:|---:|
| coverage on **today's** panel | +0.476 | +0.482 |
| **log₁₀(perturbations in the source)** | **+0.940** | **+0.933** |

**Source SIZE, not present panel overlap, is what survives the panel change** — and that is the axis on which any
future October acquisition or stack decision should be judged. It also explains 233's own attribution for free:
CD4 collapses because its coverage does; the H1 union does not move because it is sized on carrying rows.

**Consequence:** no breadth arm is proposed, `kolf2` stays unselected, and READ-FIRST §7's coverage line should be
read as the *primary* October instrument rather than a footnote.

### 2026-09-08 12:20 UTC+3 — Correction 6, sharpened: the runbook contains the rule I broke, and the follow-on was already run

`docs/DEF-RUNBOOK.md`:350 states, in bold: ***"never rank a source by today's overlap again."*** That is exactly
what 322's mandate did. The rule was written after 233 found `kolf2` and `feng_ipsc` rank #3 and #5 by today's
overlap and **#28 and #29 in October**. Line 342 tabulates the collapse; READ-FIRST §7 repeats it.

**And the follow-on I was about to launch is already done.** "Large sources with 0/300 coverage today should matter
in October" is **233 §5** — `rpe1`, `nadig_jurkat`, `nadig_hepg2`, `k562_essential`, priced at +5.5 % to +17.7 %
relative cosine, dose-monotone, all eight rows clearing the shape-only bar. **238 then re-priced them down** (their
gain is common-essential and Arc's panel has none; robust value 16 %), which is why they are *registered, never
selected* — verified last night. **Nothing new to run. Two greps replaced an experiment.**

**Guard recorded in memory:** in a repo with 300+ experiments, grep `docs/` and the relevant RESULTS.md for a data
property *before* writing a mandate premised on it, and check whether a rule about it already exists.

### 2026-09-09 — Codex continues the interrupted storage/documentation/scientific review, without subagents

**Expected:** the cleanup was fully lossless, the current-state page was a
reliable short guide, and the latest negative experiments closed broad
modeling directions. **Observed:** current artifacts are intact, but those
three broader descriptions overstate the evidence.

- Live check at **04:27Z**: 319-cap remains the best recorded score
  **0.21865973944346428**, now **5/821**. Leader **0.2620924638410216**;
  gap **0.0434327243975573**. Current published columns, not an inferred
  cosine, carry the decomposition in `docs/review-2026-09-09/state.json`.
- September 8 used two slots (316 and 319); no September 9 entry was
  recorded. Authenticated CLI says `can_submit=true`, and its read-only
  limits method says `limit_reached=false` (no numeric remaining count).
  **No submission, schedule or worker was started.** Existing snapshot
  process left alone. Ahmet's no-subagent/no-automatic-reset direction
  overrides the older cadence advice.
- Storage: hashed **21** retained containers, read **20** distinct
  payloads, reproduced all **12** existing ledger package/content bindings,
  and verified **17** historical source-content digests. The eighteenth
  source (255) had no usable historical digest; that gap remains explicit.
  All **42** pinned core inputs match the selected operator and frozen
  inputs. The large CD4 source separately matches `c355f535…`.
- **Four missing scored entries backfilled** into the novelty ledger:
  312/312b/316/319, now **16** entries. Their original captures' null SHA
  fields were not silently converted into historical server attestations.
- Recovery helper hardened: literal paths, whole-index validation,
  missing/duplicate payload rejection, optional semantic hash checking,
  no-clobber atomic publication and cleanup on failure. No user artifacts
  deleted; only audit-owned temporary decompressions were removed.
- The old sweep's **“never submitted” rule was invalid**: the ledger
  omitted historical scored families as well as recent entries. Exact
  uploaded-byte identity and exact historical rebuilds remain unverified.
  No `.vcc` backup object was found in the checked project S3 bucket.
- READ-FIRST replaced with a short live-state/reading map; the old page
  is preserved **byte-identically** under `docs/archive/`, hash
  `b9b25ee661b77a7918df545527b1406695d79cf032998656f188fea6761f5d81`.
  Current-state/history distinctions and missing 319 board result added.
- Scientific review corrects **PDS ≠ own-target cosine**, cross-family
  MSE inversion, pre-whitening vs delivered norm identity, overly broad
  317/318 closures, KOLF extraction vs source breadth, and the old F16
  novelty/wet-lab narrative. 319d3 remains unsubmitted; its performance
  bins stay frozen, but magnitude/noise mechanism attribution is limited.

**What changes:** prioritize a current-family, held-safe validation design,
then one controlled downstream-response modeling hypothesis; maintain a
separate independent-evidence science track. Details and stop conditions:
`docs/SCIENTIFIC-REVIEW-2026-09-09.md`. Storage findings:
`docs/STORAGE-AUDIT-2026-09-09.md`. Contributor-tier privacy question:
`docs/AGENT-DATA-SHARING-2026-09-09.md` (model-specific terms unverified;
no provider switch or data contribution enabled).

Verification: recovery regressions, takeover/novelty regressions, package
content audit, pinned-input hashes, byte-identical archive, git connectivity
and whitespace checks. The production builder/driver, frozen preregistration,
confirmation ledger and manuscript were not changed.

### 2026-09-09 — Experiment 323: current-family H1 validation completed, small-effect selection unresolved

**Expected:** a held-safe current-family adapter and all-six-metric evaluation
could test whether H1 resolves the recent donor-weighting choices.
**Observed:** the historical H1 adapter retained 80 excluded target coordinates
inside its dense whitening operator. Refit the operator from the same three
non-H1 sources on the proper H1 axis. Separately repaired evaluator arithmetic
to match Scanpy's normalization and ties; production emission is unchanged.
The historical six-metric regression meets every frozen tolerance.

Six actual-count arms, locally rebuilt official H1 anchors, and 10,000 paired
target resamples are complete. **4/4 historical board point signs agree;
0/4 paired 95% intervals resolve the differences.** The large target-shuffle
control passes. MSE's panel-level ratio and noise-credit cap are recomputed in
each resample. **Decision: NOT_QUALIFIED_AS_STANDALONE_SELECTOR.**
These are reused H1 development measurements, not A/B/C score estimates or
independent-context confirmation. H1 union remains OFF.

**What changes:** retain the repaired panel for regression and large-effect
rejection, not small-arm selection. Next is source-only feasibility for
independently repeatable downstream-response residual labels, with explicit
novelty checks before another model. Primary record:
`experiments/323-current-family-validation/RESULTS.md`.

Verification: 10 tests passed; all three frozen input manifests rechecked;
best 319 package SHA unchanged. Authenticated status at **08:19:43 UTC**:
same published entry and **0.21865973944346428**, rank 6 (rank 5 at 08:09:50).
No submission, production edit, paid compute, new reserve read or push.
No research worker remains; existing snapshot process left alone. Ahmet's
subsequent takeover request authorized updating the existing paused hourly
continuation in this task: solo research, no automatic submissions or rollback.

### 2026-09-09 — Experiment 324: crossed residual-label feasibility, metadata only

**Expected:** the historical 50 panel targets across multiple studies and
lines might not supply independently repeated same-context labels.
**Observed:** seven metadata known-answer checks reproduce 270's inventory.
The panel has zero literal same-two-accession × two-context rectangles;
same-method repeats cover only six targets / ten target-method-contexts
(nine with stable Cellosaurus IDs). Publication, sample/control and construct
independence are not recorded. Across all metadata, same-method rectangles
cover only MESP2 and RUVBL1, neither on the challenge panel; their provenance
is also unresolved.

**What changes:** no response-label pilot is earned by these counts. This
does not close unbalanced multi-study models or broader response prediction.
Next is one official-metadata provenance pass on the ten same-method panel
pairs. In particular KHSRP could become a publication-level crossing if its
two ENCODE experiments share a properly established study grouping.
Do not mistake accession identity for publication or construct independence.

Seven tests pass, including duplicate-edge, alias and method-mixing
falsifiers and a synthetic shared-estimate-noise check. No response matrices,
reserved outcomes or new downloads read; no prediction/scorer, production
change, submission or paid resources. Protocol frozen at `1cbda12b`, primary
result `experiments/324-residual-label-feasibility/RESULTS.md`.
