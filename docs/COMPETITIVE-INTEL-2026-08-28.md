# Competitive & literature intelligence — 2026-08-28

Scope: (1) who is winning the 2026 Virtual Cell Challenge and what their metric signatures
say they are doing, (2) the published state of the art in perturbation-response prediction,
(3) the 2025 challenge, (4) AI-for-science agent methodology.

**Confidence convention used throughout.**
`VERIFIED` = computed from `data/leaderboard/lb_*.json` or read directly in a primary source.
`INFERRED` = my reasoning from verified facts, stated as such.
`REPORTED` = a web page / repo / self-reported `model_name` string says so; I did not confirm it.
Everything a competitor writes about their own method (the `model_name` field) is self-reported
and adversarially unreliable — treated as `REPORTED` throughout.


---

> ## ⚠️ ADDENDUM 2026-09-07 — READ THIS BEFORE ANY SECTION BELOW
>
> This file is **ten days stale** and two of its central readings are now falsified by measurement. The current
> picture is `docs/READ-FIRST.md` **§1.1** (recomputed from the board's own six member columns) and
> `experiments/306-recon-oct/RESULTS.md` (reconnaissance, 2026-09-07). What changed:
>
> - **§1.3's "the leader's rise is 100 % the `mse` member" is DEAD.** `cqawesome` went 0.216873 (09-02) →
>   0.238932 (09-06 21:00Z) → **0.2393600** (09-07 15:29Z) gaining on **all six** members, not one:
>   `pds` +0.0355, `mse` +0.0405, `nmae` +0.0160, `fid` +0.0118, `reach` +0.0222, `jac` +0.0064 over 33
>   submissions. That is a real model improvement, not a metric-gaming signature.
> - **The gap is transfer cosine, and it is sized.** Ours **0.09728**, theirs **0.14456** (197's inversion under
>   §4.1c's validity guard) — cosine explains **76 %** of it. Clearing our gate needs **+25 %** relative cosine;
>   matching the leader needs **+48.6 %**. Our best-ever arm was +5.16 % (224) and `nmae` killed it (228).
> - **The frontier is silent and it is not using data we lack** (306): of 769 entries, **zero** of the two ≥ 0.220
>   describe a method; the best entry that describes one is rank 42 at 0.17543. Every rival who does describe a
>   method names **only sources already on this disk** and scores 0.136–0.150. At least four teams have
>   independently reconstructed the incumbent's own recipe.
> - **There is no new data to buy.** X-Atlas/Pisces has been an empty repo for 165 days; GEO, ArrayExpress, Zenodo
>   and scPerturb return nothing qualifying since the 2026-08-30 freeze. With 225's measured saturation of the
>   public-source channel (+5 % of board cosine), **the leader's advantage is method, not data.**
>
> Everything below is kept for its arithmetic (the scoring map, the emission identities, the account analysis),
> which still holds. Do not act on its competitive conclusions.

## TL;DR — the five findings that matter

1. **The mse baseline is not "the mean response", it is "emit the control".** Arc's spec and
   scorer source both say `raw = 1` is a control-emitting prediction. **Our 1.3627 means our
   profiles are 36% further from the truth than the control is.** No realized constant-vector
   submission on the board reads near 0.99 — they read 1.16–2.71 — because the `ρ` cap denies
   noise credit to zero-spread submissions. (Part 0)
2. **The diagnosis is calibration, and the fix is mostly to RAISE emission dispersion.**
   *(This item was rewritten after experiment 096 measured the mechanism; my original geometry
   here was wrong by an order of magnitude — see the retraction at the head of Part 0.)*
   The metric's noise credit `ρ·claim` is capped at the reference's own dispersion (`#247`), and
   **we sit at 0.67× of it, with the cap firing on 0 of 900 rows** — leaving 0.0798 of claim
   unclaimed, worth **0.349 of raw mse**. Because the predictions' common mode is only 0.84% of
   their energy, the across-perturbation budget `B` tracks the amplitude almost exactly, so
   `mse(s)` has **no interior optimum**: it falls until `ρ` saturates and rises after. The
   optimum is the cliff itself, `s_c = √(claim/B)`, measured at `VCC_ALPHA ≈ 0.83` for our
   current emitter. **Raising dispersion raises the cliff too**, so both levers point the same
   way. Measured recipe: **raw mse 0.8785, `score_mse` +0.115, +0.0192 overall** (≈ rank 3–4),
   for a **6%** amplitude cut. Adding per-gene shrinkage keyed on *donor agreement* (Miller et al.)
   takes it to **+0.0236 (transferable) to +0.0289 (optimistic)** — the optimistic end lands below
   the leader's 0.8270. Ship order: amplitude-to-cliff first (one env var, zero build risk), then
   dispersion, then shrinkage — never dispersion and shrinkage together, since both move `B` and
   the result would be unattributable.
3. **The leader is ahead on calibration AND on direction — the second part corrected late.**
   Their rise is visible on the board as almost pure mse (v11→v12: **+0.02139 of a +0.02131 jump**,
   with `reach` actually falling). But inverting their 0.8270 correctly puts their aggregate cosine
   at **≥ 0.0947, at least 1.28× ours** — the `#247` cap bounds `claim` from above and therefore
   bounds their implied `cos` from *below*. My earlier "0.97×–1.28×, plausibly zero gap" was
   arithmetic that assumed uncollectable noise credit; it is retracted in Part 0.
   **Consequence: closing our credit gap reaches raw ≈ 0.870 against their 0.8270 — roughly 4th on
   this member, not 1st.** The cheap interventions are still worth +0.0192; they are no longer a
   route to the top of the board on their own.
4. **mse is a gate, not one of six levers.** An entry that never collects mse cannot exceed
   ~0.208 even if it matched the board maximum on all five other members, and realistically caps
   near 0.176 (AgingLab_BI). The leader is already past that. We hold the board's best `pds`
   with **zero headroom left there**. (§1.5)
5. **The literature already solved this, and there is a published formula.** Miller et al.'s
   per-gene shrinkage `μ_ID = α ⊙ μ_TD + (1−α) ⊙ μ_all`, `α = 1 − p_DEG`, strictly dominates a
   global scalar. COMPASS (same six cell lines as VCC 2026) shows the **shared-response magnitude
   coefficient transfers across lines at r = 0.66 while the gene-specific residual transfers at
   r = 0.22** — magnitude is ~3× more transferable than direction detail, and we are monetising
   only the hard half. (§2.5)

**What I got wrong, recorded deliberately.** My first pass modelled this as pure geometry
(`raw = r² − 2rc + 1`), inverted the leader's score to `c ≥ 0.416`, and recommended shrinking to
0.36–0.40× amplitude. The omitted noise terms are larger than the signal being solved for, so
that inversion was off ~4.5× and the recommendation inverted the correct sign on dispersion. The
algebra was self-consistent and reproduced the board's ordering, which is exactly why it read as
confirmed. Full retraction at the head of Part 0.

---

# PART 0 — the correction that reframes everything else

The framing in `docs/BOARD-ANALYSIS-2026-08-28.md` — that the mse baseline `b ≈ 0.993` is
"the context's mean perturbation response assigned identically to every perturbation" — is
**wrong**, and the correct reading makes the member far more attackable.

`docs/reference/vcc2026-metrics-FULL.md` §2 states it outright (`VERIFIED`, primary source):

> A prediction that emits the control unchanged has expected numerator equal to the expected
> denominator, so no skill is 1.

So **raw = 1 is "emit the control"**, not "emit the mean response". Our 1.3627 means our
predicted profiles sit **36% further from the truth than the control does**. We are not slightly
behind a clever baseline; we are behind *doing nothing*.

## Why no team's constant-vector submission actually reads ~0.99

Because of the `ρ` cap in the spec (`VERIFIED`, §2): the sampling-noise correction a submission
may claim is bounded by `B`, the submission's **own across-perturbation spread**. A submission
with identical rows for all 300 perturbations has `B ≈ 0`, hence `ρ ≈ 0`, hence **no noise credit
at all**. The spec says as much: a submission with real signal has `ρ = 1` and is unaffected.

This is confirmed by four teams who submitted exactly those baselines (`VERIFIED`, board data):

| team | self-described model | mse raw | pds_cosine |
|---|---|---:|---:|
| #297 carp | `Baseline 0: NTC-copy` | 2.7054 | 0.5000 |
| #248 em2.bio | `mean_realized: one realized-mean block per context, identical across all 300 perturbations` | 1.8417 | 0.5000 |
| #276 Amy Francis | `v0 baseline (context mean + NB noise)` | 1.7135 | 0.5006 |
| #261 Rhowan | `context-mean NB baseline v0` | **1.1630** | 0.5148 |
| — | **ours (Wildflower)** | **1.3627** | 0.8318 |

`pds_cosine = 0.5000` is exact chance, so those rows genuinely are identical across perturbations.
The 0.68 spread between em2.bio (1.84) and Rhowan (1.16) for nominally the same idea is
attributable to emission/rendering choices alone.

## RETRACTED — my geometry was wrong, and the correction inverts the recommendation

**What I originally wrote here** (and propagated to the TL;DR, §1.9, §1.10 and Part 5): treating
the metric as `raw = ‖δ̂ − δ‖²/‖δ‖² = r² − 2rc + 1`, inverting the leader's 0.8270 as `1 − c²` to
get `c ≥ 0.416`, and concluding that we should **shrink our delta to ≈0.36–0.40× of current
amplitude**, with a flat basin over `s ∈ [0.30, 0.50]` worth ≈ +0.028 overall.

**That is wrong.** Experiment 096 has since measured the mechanism end-to-end against
`cell_eval2` (write-up: `experiments/096-mse-floor/MSE-MECHANISM.md`), and my simplification
dropped two terms — the constant `Kc` and the noise credit `ρ·claim` — which are **larger than
the quantity being solved for** (measured at 0.02 and 0.52–1.08 of `sumD`). Inverting a formula
whose omitted terms dominate produced a `c` that was off by ~4.5× and a shrink recommendation
that was off in *magnitude* by an order of magnitude. I am recording it rather than deleting it
because the failure mode is instructive: **the algebra was self-consistent and reproduced the
board's ordering, which is exactly why it read as confirmed.**

### The measured identity

```
MSE − 1  =  ( A − 2X + Kc − ρ·claim ) / sumD          verified against cell_eval2 to 9e-16
A = Σ_p w‖δ̂_p‖²      X = Σ_p w⟨δ̂_p, δ_p⟩      B = across-perturbation budget
```

**The key structural fact, and the reason my version failed:** `A − B` is *exactly* the squared
norm of the predictions' **common mode**, and it is only **0.84% of `A`**. So `B ≈ A` — the
across-perturbation spread is essentially all of the prediction energy — and under a global scale
`s` the noise credit tracks the energy almost exactly, until `ρ` saturates.

### The optimum is the `ρ` cliff itself; there is no optimum below it

- **Below the cliff** (`ρ < 1`): `ρ·claim = s²B`, so the numerator collapses to
  `s²(A−B) − 2sX + Kc`. Because `A−B` is 0.84% of `A`, this is **decreasing in `s`** all the way
  out to `s = X/(A−B) = 7.9`. Shrinking always *hurts*.
- **Above the cliff** (`ρ = 1`): the credit is constant, numerator is `s²A − 2sX + Kc − claim`,
  minimised at `s = X/A = 0.066`. Shrinking always *helps*.
- **Therefore the minimum is exactly at `s_c = √(claim/B)`.**

Measured on h1-test with full re-emission and re-jackknife at every grid point: **argmin at
`s = 0.42` (mse 0.9637), `ρ = 0.963` there, `ρ` hits 1.000 at `s = 0.43` and the curve turns
up.** Predicted cliff 0.4290. My "the argmin will be biased high" was right; my "the true
optimum lies below it" was wrong — nothing lies below it.

### Our cliff is much higher than I estimated, so the fix is a 6% cut, not a 60% one

My `α ≈ 1/√3 = 0.577` used the spec's "factor of about 3", which is measured **on the reference
contexts, not on our submission**. Measured `B/claim` from our own shipped `.vcc` files:

| submission | `B/claim` | cliff `s_c` | ⇒ `VCC_ALPHA` |
|---|---:|---:|---:|
| v11 (φ=0) | 2.37 | 0.65 | 0.70 |
| a8 / Wildflower (φ=4) | 1.68 | 0.77 | 0.83 |
| a10 (s0.15 φ16) | 1.66 | 0.78 | 0.84 |

### The real lever is the opposite of what I said: RAISE emission dispersion

`#247` caps the claim at the reference's own `C_p`. **We are under that cap, not over it.**
Measured: reference `C_p` = 14.54 / 14.77 / 16.37 on contexts A/B/C; our emitted `Ĉ_p` at φ=4 is
**10.30 — only 0.67×**. The cap fires on **0 of 900 rows**. We are leaving **0.0798 of claim on
the table, worth 0.349 of raw mse.**

Since `s_c = √(claim/B)`, raising `claim` *also raises the cliff*, so **less amplitude shrink is
needed.** Both effects point the same way. My §1.10(b) had this exactly backwards.

### The master formula, and what it says about the leader

At the cap:

```
MSE_min  =  1 + Kc/sumD − 2·cos·√(claim/sumD)
```

Solving `X` and `sumD` from two of our own board points differing only in `VCC_ALPHA`
(v11: 1.5801 at α=1.0796; v12: 0.9640 at α=0.40) gives `X = 0.018846`, `sumD = 0.228572`, hence
`MSE_min = 1.0237 − 2.077·cos` for our emitter. **Out-of-sample falsifier**: a8/Wildflower was
not in the fit and differs in emission and the pruned whitener — predicted 1.3671, observed
1.3627, **+0.32%**.

Our **deconfounded** cosine is therefore `X/√(A·sumD) = 0.0740` — **not 0.42**.

**RETRACTED — my sensitivity table was wrong, and it flattered us.** I originally wrote that
`sumD` is a property of the data while `Kc` and `claim` are properties of the submission, so
inverting the leader's 0.8270 with *our* `claim` assumes their emitter matches ours. That framing
is right. **The table I built from it was not**, and the error reversed the conclusion.

I took the coefficient 2.077 as "same claim as ours" and scaled it *up* to model a leader nearer
the `#247` cap. But **2.077 was already the at-the-cap value**: `(2.077/2)² = 1.0785`, matching
the measured cap `claim/sumD = 1.07856` exactly. My two rows below it implied `claim/sumD` of
**1.609 and 1.867 — both above the physical cap**, i.e. noise credit no submission can collect.
I also had the monotonicity backwards: `coef = 2√(claim/sumD)` is *increasing* in claim, so
`cos = target/coef` is *decreasing* in it. **The cap is a ceiling on claim, therefore a FLOOR on
`cos` — not a mid-range estimate.**

Corrected, with `2·cos·√(claim/sumD) = 1 + Kc/sumD − 0.8270 = 0.19667`:

| leader's emitter | `claim/sumD` | coefficient | implied `cos` | vs our 0.0740 |
|---|---:|---:|---:|---:|
| φ=0, like our v11 | 0.51947 | 1.4415 | 0.1364 | 1.84× |
| φ=4, like our a8 (current) | 0.72958 | 1.7083 | 0.1151 | 1.56× |
| **at the `#247` cap (the maximum)** | 1.07856 | 2.0771 | **0.0947** | **1.28×** |

**The leader's aggregate cosine is AT LEAST 0.0947 — at least 1.28× ours — and higher if their
emitter falls short of the cap. There is no branch in which their direction is no better than
ours.**

Robustness to the one remaining transferred constant (`Kc` is submission-specific, and I am
using ours): across a range twice as wide as the measured 0.02, the floor stays at
**1.13×–1.45×** ours. **The direction of the conclusion is robust to `Kc`; the exact multiple is
not.**

**The material consequence: calibration alone does not reach them.** Our ceiling at `cos = 0.0740`
with claim at the cap is `1.02367 − 2.0771×0.0740 = 0.8700`, against their observed **0.8270** —
and that shortfall holds under every `Kc` in range (0.8463 / 0.8700 / 0.8963 vs 0.8270). Closing
the credit gap takes this member to roughly **4th, not 1st**.

**So "their lead is calibration, not direction" — which I authored and which propagated into both
documents — is withdrawn.** Their lead is *partly* calibration, which is the cheap part we can
take; the residue is real directional accuracy that we do not have. Raising `cos` moves from
"the expensive option" to **the only route to first**.

### The corrected recipe and its price

096's measured proposal: top-8000 genes at (spread 0.15, φ 16), tail at (spread 1.0, φ 4) — which
brings `claim` to the cap at only **1.16× reference dispersion** — then `VCC_ALPHA` 1.0796 → 1.012.

| | predicted |
|---|---:|
| raw mse | **0.8785** |
| `score_mse` | **+0.115** |
| overall | **+0.0192** → ≈ 0.1736 |
| `reach` cost | ≈ −0.001 |

That is ≈ rank 3–4, not rank 1. **Sober, measured, and still by far the largest single lever
available to us** — and it costs a 6% amplitude cut rather than the 60% I originally claimed.

---

# PART 1 — who is winning, and what their signatures say

Source: `data/leaderboard/lb_20260828T1325.json`, 385 teams, 2,010 total submissions, plus the
nine earlier snapshots in `data/leaderboard/`. All arithmetic `VERIFIED`.

## 1.1 The scoring map, solved exactly

The board publishes both raw and scaled members. Fitting `score_mse = a·raw + c` on the six
well-behaved positive-mse teams reproduces them to ~1e-16:

```
score_mse = clip( (0.987863 − raw) / 0.950163 , 0, 1 )      a = −1.052451,  c = +1.039678
```

so the baseline anchor is `b = 0.987863` and the implied "perfect" anchor is `0.0377`. Exact fits:
Aginglab.com, Ibrahim Mansour, BioFarmer, IsoCell, black_handsome_egg, enough1.

**Per-context clamping is confirmed.** Exactly 9 of 385 teams have `score_mse > 0`; **zero** have
`score_mse < 0`. Three of the nine sit *above* the affine line and **none sit below**:

| team | mse raw | score_mse | residual vs affine |
|---|---:|---:|---:|
| We527 | 1.3464 | 0.0532 | **+0.4306** |
| Immovable Myhal Rock | 0.9906 | 0.0059 | +0.0087 |
| nanobot | 0.9918 | 0.0014 | +0.0055 |

One-signed residuals are the signature of `clip(·,0,1)` applied *per context* then averaged
(Jensen). Corroborating: `anchor_version = "vcc2026-valA-r4+vcc2026-valB-r4+vcc2026-valC-r4"` —
three separate per-context anchor sets. Under that reading We527 has **one context at raw ≈ 0.836
and two at ≈ 1.60** (mean 1.346, checks out). The same one-signed residual pattern appears at the
1e-3 scale in `pds`, `reach`, `jac`.

**Upgraded from inferred to PROVEN.** We527's two scored submissions have *different* raw values
and *identical* scores to 15 significant figures:

| snapshot | model | `expr_mse_unbiased_capped_norm` | `score_mse` |
|---|---|---|---|
| 08-27 16:17 | `STATE baseline` | 1.357686703495**0901** | 0.05319847662453**766** |
| 08-28 13:25 | `STATE baseline` | 1.346449203277**8235** | 0.05319847662453**7694** |

The raw moved by 0.0112 and the score did not move at all — while athenen at raw **1.003978**
scores exactly **0**. No global map can be locally constant across that interval, strictly
positive there, and zero at a *lower* raw. Per-context clamping is the only structure that fits.
It further implies the raw change happened **entirely inside the two clamped contexts** and their
one scoring context was bit-identical between submissions.

Practical consequence, unchanged from `BOARD-ANALYSIS`: **beating the baseline in ONE context
collects a third of the member with no penalty in the other two.** One context at raw 0.83 is
worth +0.0092 overall (≈ rank 5); all three at 0.83 is worth +0.0277 (≈ rank 1).

Other members, robust affine fits (`VERIFIED`): `score_pds = 2.2151·raw − 1.1078` (baseline raw
0.5001, i.e. chance); `score_fid = 3.3621·raw − 1.7230` (baseline 0.5125); `score_reach =
1.1292·raw − 0.0893` (baseline 0.0791); `score_jac = 2.7232·raw − 0.0826` (baseline 0.0303);
`score_nmae = −1.6617·raw + 1.6626` (baseline 1.0005) **with a hard floor at −6.0000** — several
teams are pinned there from unnormalised counts.

## 1.2 The top 16, z-scored against the top-120 field

Scaled members, z-scored; `+` = strong for this part of the board, `−` = weak.
**Snapshot: `lb_20260828T1325.json`.** A later pull the same day (16:46) has the board moved —
notably Mike Chen (VIVAI) up to #6 having *cleared the mse floor* at `mse_raw` 0.9183; see §1.11.
Family D is therefore no longer a stable description of that cluster.

| # | team | avg | pds | mse | nmae | fid | reach | jac | family |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---|
| 1 | Aginglab.com | 0.21122 | +1.66 | **+7.31** | **+2.18** | +1.16 | +1.79 | +0.22 | **A — complete** |
| 2 | AgingLab_BI | 0.17637 | +1.61 | −0.19 | +1.70 | +1.08 | +1.48 | +1.11 | B — broad, no mse |
| 3 | We527 | 0.17054 | +1.22 | +2.17 | +1.59 | +0.79 | +0.77 | +0.87 | A — partial |
| 4 | Dewei Hu | 0.17017 | +0.80 | −0.19 | +1.34 | **+1.95** | **+2.36** | +0.75 | C — DE-reach |
| 5 | Ibrahim Mansour | 0.16890 | +1.66 | **+6.72** | +0.44 | −0.27 | −0.36 | +0.13 | A — mse-only |
| 6 | Cell-tics | 0.16248 | +1.43 | −0.19 | +0.99 | +0.72 | +1.22 | +0.85 | B |
| 7 | Steve | 0.16104 | +1.68 | −0.19 | +1.26 | +0.44 | +0.72 | **−0.93** | D — VIVAI |
| 8 | cc | 0.16022 | +1.71 | −0.19 | +1.20 | +0.58 | +0.55 | **−1.09** | D — VIVAI |
| 9 | Kate P | 0.16014 | +1.71 | −0.19 | +1.30 | +0.32 | +0.64 | **−1.24** | D — VIVAI |
| 10 | Jurassic Park | 0.15944 | +1.05 | −0.19 | +1.51 | +0.65 | +0.88 | +1.64 | B |
| 11 | Mike Chen | 0.15940 | +1.66 | −0.19 | +0.84 | +0.67 | +0.83 | −0.61 | D — VIVAI |
| 12 | VIVAI Inc | 0.15852 | +1.66 | −0.19 | +0.95 | +0.56 | +0.68 | −0.65 | D — VIVAI |
| 13 | UM | 0.15794 | +1.47 | −0.19 | +0.68 | +0.61 | +1.37 | −0.91 | D — "Boom" |
| 14 | Weijing Zhao | 0.15781 | +1.49 | −0.19 | +0.58 | +0.67 | +1.40 | −0.87 | D — "Boom" |
| 15 | GISL | 0.15753 | +0.93 | −0.19 | +1.01 | +1.11 | +1.21 | +1.53 | B |
| **16** | **us (Wildflower)** | **0.15436** | **+1.94** | −0.19 | +0.58 | +0.21 | −0.08 | +0.62 | **E — pds-extreme** |

### The four families

**Family A — magnitude-calibrated (Aginglab.com #1, Ibrahim Mansour #5, We527 #3).**
The only teams collecting `mse`, and the only ones above 0.168. Aginglab.com is the sole
*complete* signature on the board: best mse, best nmae, 3rd-best reach, positive fid.
Ibrahim Mansour is the pure form — high mse and pds, **negative** fid and reach — i.e. a
well-calibrated mean profile with no DE-callability engineering at all. `INFERRED`: Ibrahim
is solving the amplitude problem and nothing else; Aginglab solved amplitude *and* the DE
members.

**Family B — broad, no-mse (AgingLab_BI #2, Cell-tics #6, Jurassic Park #10, GISL #15).**
Positive on all five non-mse members, `mse = 0`. AgingLab_BI at 0.17637 is essentially the
**ceiling of the "ignore mse" strategy** — see §1.5.

**Family C — DE-reach maximisers (Dewei Hu #4).** Extreme fid (+1.95) and reach (+2.36) with
only mid pds. Reaches 0.170 on the Wilcoxon members alone.

**Family D — direction maximisers with negative jaccard (VIVAI ×5, "Boom" ×2).** A distinctive
fingerprint: high pds, high nmae, **negative jac**. `INFERRED`: they push discrimination hard
and under-call significance, paying the jaccard floor for it.

**Family E — us, alone.** The board's single most pds-extreme signature (+1.94, and the outright
best raw `pds_cosine` 0.8318 of all 385), zero mse, and — unlike family D — *positive* jaccard
and neutral reach. We are the only entry in the top 16 that is simultaneously pds-extreme and
not paying a jaccard penalty. That is a real asset; it is just attached to a broken amplitude.

## 1.3 The leader's rise is 100% the mse member

From our own snapshot history (`VERIFIED`):

| snapshot | model | avg | pds | mse raw → score | nmae | reach |
|---|---|---:|---:|---|---:|---:|
| 08-25 | GeroAI_v9 | 0.15446 | 0.7255 | 1.8219 → **0.0000** | 0.1314 | 0.1000 |
| 08-27 | GeroAI_v11 | 0.18991 | 0.7082 | 0.9594 → 0.0410 | 0.1784 | 0.2128 |
| 08-28 | GeroAI_v12 | **0.21122** | 0.7089 | **0.8270 → 0.1693** | 0.1838 | 0.2066 |

The v11→v12 step decomposed over the six published `score_*` fields (`VERIFIED`, reproduces the
observed `score_avg` delta to a residual of −3.1e-17):

| member | v11 | v12 | Δ | contribution (/6) |
|---|---:|---:|---:|---:|
| **mse** | 0.0410 | 0.1693 | **+0.1283** | **+0.02139** |
| nmae | 0.1784 | 0.1838 | +0.0054 | +0.00090 |
| pds | 0.7082 | 0.7089 | +0.0007 | +0.00012 |
| jac | −0.0044 | −0.0042 | +0.0002 | +0.00003 |
| fid | 0.0034 | 0.0030 | −0.0005 | −0.00008 |
| **reach** | 0.2128 | 0.2066 | **−0.0062** | **−0.00104** |
| | | | **total** | **+0.02131** |

**+0.02139 of a +0.02131 jump is mse. All five other members net −0.00007, and `reach` *fell*.**

This is the strongest single piece of evidence in this document, because it is an **external
replication on someone else's build**: the mechanism experiment 096 measured internally — buy the
mse member, pay essentially nothing elsewhere, `reach` if anything slightly down — is visible in
the leader's published scores, from the board alone, with none of our instrumentation involved and
nothing fitted.

**What it says about the leader's directional accuracy.** Their `pds_cosine` raw moved
0.819621 → 0.819937 (+0.0003) while `mse` raw went 0.9594 → 0.8270. `INFERRED`, with one caveat
that matters: `pds_cosine` is a **discrimination** statistic — the rank of the true perturbation
among candidates under cosine distance — *not* the bulk directional agreement `X/√(A·sumD)` that
the master formula's `cos` denotes. A rank statistic can sit flat while the bulk cosine moves,
because ranks saturate. So flat `pds` is strong but **not conclusive** evidence that `cos` was
fixed.

The conclusion survives anyway, by a second route: under `MSE_min = 1 + Kc/sumD −
2·cos·√(claim/sumD)`, improving 0.9594 → 0.8270 **at fixed `cos`** requires raising
`√(claim/sumD)` — i.e. claiming more noise credit. That is the **high-claim branch** of the
sensitivity table in Part 0, and the high-claim branch is precisely where the leader's inferred
`cos` is *lowest* — but "lowest" is **0.0947, still 1.28× ours** (Part 0 retraction; my earlier
0.0775/0.0720 figures assumed claim above the physical cap). **So the correct reading is: they
out-calibrated us *and* they are ahead on direction. This decomposition shows the calibration move
they made between v11 and v12; it does not show the directional edge they already had at v11.** And on 08-25 they were at 0.15446 — *statistically indistinguishable from where we
are now* — with the same broken mse. **They were us three days ago.** Their pds went *down*
(0.7255 → 0.7089) when they fixed mse, which is exactly what the geometry predicts if the fix was
an amplitude/calibration move (cosine is scale-free, so a pure rescale should not move pds; the
small drop suggests a rescale plus a change to the emitter). Under the corrected mechanism
(Part 0) the likeliest reading is that they moved onto the `ρ` cliff and raised their claim.

Ibrahim Mansour shows the same story over 15 submissions, monotonically: mse raw
0.9187 → 0.9095 → 0.8703 → 0.8348 → 0.8397. IsoCell has held 0.9197 since 08-25 without moving.

## 1.4 Aginglab.com vs AgingLab_BI — related, and the answer is instructive

`VERIFIED`: **identical `organization` string "AgingLab", different `team_id`s**
(`n2IRzh…` vs `cderAl…`), different member lists (`Suan Cho` vs `Sungwon Jeon, Hyojung Ryu`),
different model naming (`GeroAI_v12` vs `pinned-plus-gram-whitening`). They are the only two
teams on the board sharing that org string.

`INFERRED`, and this is the useful part: **their signatures are the same shape except for mse.**
Both have pds ≈ 0.818/0.817 raw, nmae 0.889/0.905, fid 0.5134/0.5126, reach 0.2621/0.2502,
jac 0.0288/0.0321 — extremely close. The single difference is mse: 0.8270 vs 2.7357. Two
accounts from one organisation running near-identical pipelines where **only one has the
amplitude fix**. That is the cleanest natural ablation on the entire board, and it prices the
fix at **+0.0349 overall** — from rank 2 to rank 1 by a margin larger than the gap between
ranks 2 and 12.

Note also that AgingLab_BI's mse raw of **2.7357** is close to carp's NTC-copy 2.7054 — `INFERRED`,
that is what a strongly *over*-amplified or `ρ`-penalised submission reads.

## 1.5 The ceiling of the "ignore mse" strategy

`VERIFIED`. Taking the board-wide maximum of each scaled member:

| member | board max | held by |
|---|---:|---|
| score_pds | 0.7366 | **us** |
| score_mse | 0.1693 | Aginglab.com |
| score_nmae | 0.1838 | Aginglab.com |
| score_fid | 0.0313 | f15 (#157) |
| score_reach | 0.2818 | jw_13288 (#282) |
| score_jac | 0.0166 | xiaoyan chen (#25) |

- Frankenstein ceiling (board max of every member): **0.23655**
- …with `score_mse` pinned to 0: **0.20833** — and the best actually-achieved no-mse entry is
  AgingLab_BI at 0.17637.
- Leader: 0.21122. Us: 0.15436.

**So an entry that never collects mse cannot exceed ~0.208 even if it simultaneously matched the
best team on the board at all five other members**, and realistically caps near 0.176. The leader
is already past that ceiling. `INFERRED`: **mse is not one of six levers, it is the gate.**

Our own headroom, per member:

| member | ours | leader | board max | gap→leader | gap→max |
|---|---:|---:|---:|---:|---:|
| pds | 0.7366 | 0.7089 | 0.7366 | −0.0276 | 0.0000 |
| **mse** | 0.0000 | 0.1693 | 0.1693 | **+0.1693** | +0.1693 |
| nmae | 0.0943 | 0.1838 | 0.1838 | +0.0895 | +0.0895 |
| fid | −0.0287 | 0.0030 | 0.0313 | +0.0317 | +0.0600 |
| reach | 0.1243 | 0.2066 | 0.2818 | +0.0822 | +0.1574 |
| jac | −0.0004 | −0.0042 | 0.0166 | −0.0039 | +0.0169 |
| total /6 | 0.15436 | 0.21122 | 0.23655 | +0.05686 | +0.08219 |

We already hold the board maximum on pds. There is literally zero headroom there.

## 1.6 The `fid` member is nearly dead

`VERIFIED`: `de_wilcoxon_direction_fidelity_yield_raw` has a chance level of exactly 0.5000
(every constant-row submission reads 0.4991–0.5029) and the **board-wide maximum over 385 teams
is 0.5222** (#157 f15). See §1.9 for *why*: at `n_pred ≥ n_real` the coverage term saturates and
`fid` reduces to pure directional precision, so these numbers say the whole field — us at 0.5041,
the leader at 0.5134 — is calling direction barely better than a coin flip. The leader is at 0.5134. Best-5: 0.5222, 0.5211, 0.5155, 0.5140, 0.5134.
Nobody has moved this member meaningfully. Our 0.5041 costs us −0.0287 scaled, and the *entire*
available upside is +0.0600 scaled = +0.0100 overall, of which realistically ~+0.005 is
reachable. **Deprioritise.** (`jac` is similar: board max 0.0166 scaled.)

## 1.7 Coordinated multi-account submission — confirmed, but not "a hundred teams"

The Discord claim that "Vivai has a hundred teams submitting back to back" is **not supported by
the board**, but a smaller, clearly deliberate version of it **is**.

`VERIFIED` — **there is a submission budget of 16.** The maximum `submission_count` over all 385
teams is exactly 16, with 15 teams at 15–16. That is a cap, not a coincidence.

`VERIFIED` — **the VIVAI experiment counter runs continuously across five separate `team_id`s:**

| team_id | display name history | submissions | model |
|---|---|---:|---|
| `Rg2VUo…` | "VIVAI Inc" (08-23, 08-25) → **renamed "cc"** (08-28), member `Arshia1` | 12 | `m8` → `vivai-m15` → `m47` |
| `Gek0pu…` | "VIVAI Inc" (new account, 08-28), member `Arshia Alaiee` | 15 | `vivai-m43` |
| `kVDhWZ…` | "Mike Chen" | 3 | `vivai-m29` → `vivai-m39` |
| `T1QSoS…` | "Steve", member `Nanoxd`, org "P" | 3 | `vivai-m44` |
| `UnUEbG…` | "Kate P" | 3 | `probe` → `vivai-m46` |

One counter — m8, m15, m29, m39, m43, m44, m46, m47 — across five accounts, **36 submissions
against a 16 cap**, and the original account was renamed while a fresh one took the old name.
Three of the five accounts (Steve, Kate P, Mike Chen) have exactly 3 submissions each and appear
only at the top of the board.

`VERIFIED` — their signatures are a tight, isolated cluster. Scanning all 385 entries by
z-scored 6-raw-vector distance to the VIVAI centroid: the five VIVAI accounts sit at
d = 0.029, 0.094, 0.117, 0.129, 0.150; **the next nearest team on the whole board is at 0.451.**
They submitted within a 14-minute window on 08-28 (Mike Chen 00:32, Steve 00:36, Kate P 00:46).

`VERIFIED` — **the same rename-and-clone pattern appears in a second cluster.** team_id `32ClOT…`
was named "OpenPert" on 08-25 and is now "luka" (org DenseAI, model `PertMind_Pro`); a *different*
team_id `Y4Wk2y…` now carries the name "OpenPert" (model `PertMind_Pro_Gram`). Ranks 39 and 17.

`VERIFIED` — a third likely pair: **#13 UM (`BoomLi`, org "Macau") and #14 Weijing Zhao (`PAMI`)**
have raw-vector distance 0.0266, the closest non-trivial pair in the entire top 30; #31 Mengran Li
runs `BoomAI`. Shared org strings elsewhere: `INFevo` ×3, `ZHEJIANG UNIVERSITY` ×3, `AgingLab` ×2,
`Soongsil University` ×2, `Central South University` ×2.

**Scale check** (`VERIFIED`): the board grew 353 → 385 teams (+32) between 08-27 and 08-28, and
only **four** accounts new since 08-27 landed in the top 60. Total board submissions 1,702 → 2,010.
A hundred coordinated accounts would show as a hundred rows; they do not. **Conclusion: ~5
coordinated VIVAI accounts evading a 16-submission cap by ~2.25×, plus at least two other
2-account clusters. The "hundred teams" figure is not corroborated by the leaderboard.**

Caveat (`INFERRED`): the board shows only each team's *best* entry, and account renames are
visible only because we snapshot over time. Both facts mean this analysis undercounts. Continue
snapshotting `team_id` → `team_name` mappings; renames are the tell.

## 1.8 What the `model_name` field leaks about the field's methods

All `REPORTED` — self-declared strings, potentially misleading, but the distribution is informative.

**Foundation models are losing, badly.** Every self-declared foundation-model run sits at or near
the bottom:

| # | team | model | avg | mse raw |
|---:|---|---|---:|---:|
| 204 | flyingfog won | `v8_scgpt_embed_beta1` | −0.002 | 11.93 |
| 217 | NotToday | `state-4model-addressing-v2` | −0.015 | 1.92 |
| 312 | Lottery | `scGPT-MM-Replogle-v1` | −0.263 | 191.2 |
| 343 | Formosa Bear | `STATE sm 10000 steps, ESM2 pert features, trained on 2025 H1` | −0.480 | 131.9 |
| 350 | The Interventionists | `scGPT four-context baseline v2` | −0.858 | 873.3 |
| 362 | cccccc | `STATE_Replogle_k562_flow_baseline` | −0.968 | 1748.9 |
| 363 | QLZZZ | `state_sm_test` | −0.978 | 950.8 |

**The single exception is #3 We527, whose model is named literally `STATE baseline`** and which
scores 0.17054 with a positive mse. `INFERRED`: either We527 found the normalisation/rendering that
the other STATE runs are missing (note the others' mse raws of 131–1749 are unnormalised-count
signatures), or the name is not descriptive. Worth watching; do not assume STATE is the answer.

**The dominant successful family is transfer from public Perturb-seq.** Recurring across dozens of
mid-board entries: **Replogle K562 genome-wide Perturb-seq (GWPS)** as the delta source, plus
**X-Atlas/Orion (HCT116)**, **Nadig HepG2/Jurkat**, **Adamson**, **RPE1**, and the **2025 VCC H1**
data. Representative strings: `replogle-k562-delta-transfer`, `K562-GWPS continuous cross-line
transfer + expression gate + NB synthesis`, `cross-line transfer (Nadig HepG2/Jurkat + Harmonizome
K562-GWPS fallback)`, `multi-context joint train (h1+k562+adamson+rpe1), LOCO-selected, 3-seed
ensemble`.

**Amplitude/calibration language is everywhere in the mid-board**, which corroborates Part 0:
`alpha=0.35`, `post_scale 0.6`, `EFFECT_SCALE=0.25`, `magnitude-scaled transfer, alpha=3.0`,
`trial6_uniform_magcal_s1.6`, `amp1.8`, `scale=1.50 after DE-yield audit`,
`magnitude-shrunk after fid −0.47 at x1.5`, `component-rescaled`. The field has found the
amplitude knob; most are tuning it against the *DE* members and therefore setting it too high.

**"Gram" appears in two independent top-20 entries** — AgingLab_BI's `pinned-plus-gram-whitening`
(#2) and OpenPert's `PertMind_Pro_Gram` (#17), plus `K562 delta + count-Gram v2` (#121).
`INFERRED`: whitening the perturbation-by-perturbation Gram matrix is circulating as a technique,
almost certainly a **pds/discrimination** trick (decorrelating predicted deltas makes each
perturbation more distinguishable) rather than an mse one — consistent with AgingLab_BI having
excellent pds and zero mse.

---
## 1.9 The DE members under the corrected recipe

`VERIFIED` from `docs/reference/vcc2026-metrics-FULL.md` §§3–5, then `INFERRED`.

Read the three Wilcoxon set/rank members' definitions carefully and they depend on the submission
**only through two things**:

1. **the sign of each predicted `lfc`** (`fid` scores sign correctness; `reach` scores sign purity
   down the submission's own ranking), and
2. **the effect-size-to-dispersion ratio per gene**, which is what sets the Wilcoxon p-value and
   therefore which genes are "called" (`jac` is the Jaccard of the called sets; `fid` penalises
   under-yield; `reach` puts significant calls first in its ordering).

The spec is explicit that the test is rank-based and magnitude-blind — it says the statistic is
unmoved by how large an outlying cell is, only by where that cell falls in the order. `reach`'s
ordering is `significant calls first → ascending predicted adj-p → ascending p → descending |lfc|`,
every tier of which is invariant under a positive scalar.

**Therefore** (`INFERRED`): the three Wilcoxon members are invariant to a positive rescaling of
the predicted delta *provided each gene's effect-to-dispersion ratio is held fixed*. `pds_cosine`
is a cosine and is scale-free, so it is invariant too.

**This section survives 096's correction, but its importance shrank.** The measured amplitude move
is `VCC_ALPHA` 1.0796 → 1.012, a **6% cut**, not the 60% I originally derived — so the DE members
were never seriously at risk from the amplitude leg. 096 prices the whole recipe's `reach` cost at
**≈ −0.001**.

**The live question is now the opposite one, and it has been measured.** The recipe *raises*
emission dispersion (top-8000 genes at spread 0.15 / φ 16, tail at spread 1.0 / φ 4, reaching
1.16× reference dispersion). I hoped the added dispersion could be spent where nothing was
callable anyway. **096 measured that and it is false for our submission** (`mech_callability.py`,
real path on challenge context A): the tail policy **churns 8.7–14% of the call set**, and
`n_pred` goes *up* +13.5% (4,772 → 5,417), not down.

### fid's coverage term is already saturated — so `fid` *is* our directional precision

`VERIFIED` from spec §3: `F_p = k/max(n_pred, n_real) = (k/n_pred)·min(1, n_pred/n_real)`, and
*"calling more genes than the reference found buys nothing, since the coverage term saturates."*
With `n_pred` 4,772 ≥ `n_real` ≈ 4,757 **the coverage term is already 1**, so

> **our `fid` raw of 0.5041 is not a coverage problem — it is our directional precision, and it is
> 0.41% better than a coin flip on the genes we call.**

That reframes §1.6: `fid` is not "nearly dead" because the metric is broken, it is nearly dead
because *nobody, us included, can call direction better than chance at this depth.*

Same exercise on `jac`: solving `J = |∩|/(n_pred + n_real − |∩|) = 0.0302` gives
**|∩| ≈ 279 genes — only ~5.9% of the reference's significant set is in our call set.**

### The whole dispersion risk collapsed to one number — which has now been measured at chance

Because coverage is saturated, the churn's cost depends *only* on the directional precision of the
~1,060 new calls. **096 measured it on both harnesses** (`mech_precision.py`, exact DE path against
real controls, 25 perturbations per surface):

| | h1-test | orion-thinned |
|---|---:|---:|
| `n_pred` a10 → tiered | 4,198 → 6,210 | 3,598 → 4,745 |
| precision `k/n_pred` | 0.5110 → 0.5062 | 0.5565 → 0.5456 |
| **precision of the NEW calls** | **0.4946** | **0.5131** |
| precision of the DROPPED calls | 0.4938 | 0.5132 |
| jaccard | 0.0751 → 0.0814 | 0.0121 → 0.0113 |

**The two surfaces bracket chance from either side, and the dropped calls are at chance too — the
churn is chance-in, chance-out.** The pessimistic branch required the new calls to be *actively
anti-correlated* with the reference; neither surface shows it, and there was never a mechanism for
it. **Risk resolved: ≈ −0.0005 overall, not −0.028.**

### Why the board will pay less than the harnesses, and the closed form

`INFERRED`, mine, and it refines the point estimate in the favourable direction. When churn is
chance-in/chance-out the above-chance mass `E = k − 0.5·n_pred` is **conserved**, and only the
denominator grows:

```
p_after − 0.5  =  (p₀ − 0.5) · n_before / n_after
Δscore_fid     =  3.362 · (p₀ − 0.5) · (1 − n_before/n_after)
```

| surface | `p₀` | `n` before → after | closed form | observed |
|---|---:|---|---:|---:|
| h1-test | 0.5110 | 4,198 → 6,210 | 0.50744 | 0.5062 |
| orion | 0.5565 | 3,598 → 4,745 | 0.54284 | 0.5456 |

The two harnesses bracket the closed form (h1 below, orion above).

**But do not promote "the harnesses bracket the truth" to a prior — and my attempt to refute it
was itself confounded.** The sequence is worth recording in full, because it is the third time in
this campaign that a confident, self-consistent argument turned out to rest on a mismatch.

I built this table to test 096's proposed rule (*"they bracket dimensionless quantities and
mechanisms; they do not bracket levels"*):

| quantity | h1-test | orion | board | brackets? |
|---|---:|---:|---:|---|
| a10 hit rate into the reference DE set | 0.0996 | 0.0123 | 0.0585 | yes |
| `expr_mse` raw **level** | 2.0551 | 1.9405 | 1.3627 | no — both high |
| deconfounded **cosine** | 0.0519 | 0.0325 | 0.0740 | no — both low |

and argued the dimensionless half was refuted, since the cosine is dimensionless and does not
bracket. **Two corrections landed on that, one of which is fatal to my own argument.**

**(i) My cosine row was mismatched statistics.** The originally-quoted 0.043/0.018 are *plug-in*
cosines `X/√(A·R)`; 0.0740 is *deconfounded* `X/√(A·sumD)`. Like-for-like the gap is much smaller
than I made it look (plug-in: 0.0432 / 0.0176 / 0.0511; deconfounded: 0.0519 / 0.0325 / 0.0740).
Both statistics still read below on both surfaces, so the direction survives — but I was arguing
against a number that had been mismatched.

**(ii) I proposed reference sampling noise as the shared systematic. Measured, and refuted.**
`mech_noise_bias.py` thins h1-test's reference 8× (400 → 50 cells/pert): the plug-in cosine
attenuates **39%**, but the **deconfounded cosine does not move** (±5%, no trend) and `sumD` is
flat to 3% — because `sumD = R − K − Kc` has the reference's sampling variance subtracted off
by construction. The correction does exactly what it was designed to do. My mechanism explains
one row of the table and not the other, which is why it needed measuring rather than reasoning.
(`n_real` is not in either cosine at all — it is a DE call-set quantity — so it could not have
been the channel even in principle.)

**(iii) The real shared systematic is provenance, and it invalidates the whole table — including
the row that bracketed.** The harness predictions come from `experiments/031-somatic/run.py`
(donor transfer + whitening); the board numbers come from `scripts/make_submission_pm_cd4.py`.
`VERIFIED`: 363 lines with zero mentions of PerturbMatch or CD4, versus 1,809 lines with 26 and
97 respectively. **They are different models**, and the harnesses additionally target held-out
cell lines while production targets the challenge's own contexts.

So every row above compares *model A on a harness* to *model B on the board*. Model provenance
varies simultaneously with measurement surface — a textbook confound of exactly the kind
`RESEARCH-METHOD` §4 exists to catch. **The table is not 1-for-3 evidence; it is n = 0 evidence.**
096's rule is not refuted and mine is not supported. Both are untested.

### The criterion I proposed was also wrong — and I mislabelled an inference as VERIFIED

I proposed: *"bracketing evidence is valid only when the identical prediction artifact is scored
on every surface,"* and cited `RESEARCH-METHOD` §1's a8 case as the founding clean example, on the
grounds that the identical `.vcc` had been scored on both harnesses and the board.

**That premise is false.** What I actually checked was that `a8-nb-phi4.vcc` exists in
`submissions/` and was submitted. I *inferred* the harness half and labelled the whole claim
`VERIFIED`. `experiments/095-clean-harness/realize.py` consumes a **profile** and drives the
production `emit_counts` against **each harness's own controls** — it never reads a `.vcc`. So
what ran on the harnesses was the same *emission policy* re-realised into a different artifact.
**Same policy, different artifact.** `RESEARCH-METHOD` §1's founding example is therefore an
instance of the confounded case, not the clean one, and cannot support the criterion I built on it.

This is the fourth self-consistent error in this thread and the only one that broke my own stated
confidence convention. Recording it as such.

### And the criterion is ill-posed anyway: artifacts are baseline-bound

096 ran the test. `orion-thinned`'s 132 targets are a 100% subset of the challenge's 300 on a
bit-identical gene axis, so the shipped `a8-nb-phi4.vcc` *can* be scored against the orion
reference. The result:

| arm | mse | A | deconf cos |
|---|---:|---:|---:|
| a8 (shipped), ctx A, scored on orion | **342.97** | 68.32 | 0.0243 |
| a8 (shipped), scored on the board | **1.3627** | 0.2830 | 0.0741 |

mse 343 against 1.36 — it is not measuring the surface, it is measuring the **cell line**. A
prediction artifact is **baseline-bound**: `b̂_p` is the challenge context's control profile plus
an effect, so scoring it elsewhere measures the control gap. Confirmed to 0.1% — the per-gene
control gap `‖b_ctrl(chal) − b_ctrl(orion)‖²/G` times the 132 rows explains **99.8–99.9%** of the
measured `A` on all three contexts. The cell-line baseline is **548× the per-row perturbation
signal**.

> **Second clause, and it is the binding one:** cross-surface scoring of an identical artifact is
> meaningful only when the surfaces **share a baseline** — same cell line, same control cells.
> Otherwise the identical artifact measures the distance between cell lines, not between scorers.

**No two of our three surfaces share a baseline, so the bracketing question cannot be settled with
anything we hold.** Both rules are recorded closed-until-new-evidence rather than either left
standing.

### The one construction that would settle it — and the trap already documented inside it

The only available design that holds the baseline fixed while varying the surface is **one
reference panel split into disjoint halves, with the same submission scored against both**. That
is exactly the spec's replicate anchor, and Arc already implements it: `cell_eval2/anchor.py`,
`_disjoint_halves`, five seeds off one pinned base seed.

**If we build it, `anchor.py`'s own docstring documents the trap** (`VERIFIED`, primary source).
`control_source` is forced to `"pred"` and is *not* optional — under `"real"`, scoring half B
against half A computes **both** halves' log fold-changes against half A's control, which

> shares the control's sampling noise between the two quantities whose agreement is being measured
> and biases the anchor upward

Arc measured it: **0.5–2.3% optimistic on `lfc_nmae`, consistently, on all 12 comparisons**, and
`ceiling.py` measured the same mistake turning a **0.54 split-half reliability into 0.74**. The
code enforces it with an assertion (`_assert_disjoint_controls`), and `control_source` is pinned
to `"pred"` in `_inner_config`. **Any split-half bracketing test we build must use disjoint
controls, or it will produce a confidently wrong answer of a size Arc has already quantified.**

**The clause I originally left out, which makes it worse** (`VERIFIED`, same docstring): an
upward-biased anchor

> biases every submission's score DOWNWARD, uniformly, with nothing in the output that looks wrong

**Silent and uniform** — the exact class of error this whole thread kept producing. A test whose
documented failure mode is *"manufactures the agreement you are testing for, invisibly"* deserves
the flag before it is built, not after.

### Why this is also a load-bearing check on Part 1, and it passes

Every scaled figure in this document rests on the recovered map `b = 0.987863`, span `0.950163`,
implied replicate `r = 0.037700` (§1.1), which is fitted **against the shipped anchors**. Had
those anchors been built with shared controls, the span would be biased and every scaled number
here would be uniformly wrong with nothing visible to catch it. They were not: disjoint controls,
assertion enforced, `control_source` pinned.

Two independent cross-checks that the recovered map is sound:

- **My fitted `r = 0.0377` sits inside the spec's independently published replicate range of
  0.028 – 0.045** (`VERIFIED`, metrics-FULL §2 "Practical values"), and the fitted `b = 0.987863`
  sits inside its published baseline range of 0.986 – 0.992.
- The `vcc2025` research agent recovered `BASELINE 0.9879 / CEILING 0.0377` independently from
  archived challenge material, matching the fit to 1e-4.

**And the robustness argument for Part 1 specifically:** a span bias is a single multiplicative
constant applied to every submission's score on that member. So **within-member ordering is
preserved exactly** — "we hold the board's best `pds`", "9 of 385 cleared the mse floor", "the
leader's rise is entirely mse" are all invariant to it. What *would* be sensitive is the
**cross-member gap decomposition** (§1.5's "mse is 50% of our gap to the leader"), since that
compares magnitudes across members with different spans. That decomposition is the one number in
Part 1 to re-derive if the anchors are ever revised.

### The generalisable lesson

**`RESEARCH-METHOD` §1 should read "the metric's *level* is unreadable", not "the metric is
unreadable".** Before declaring a member unmeasurable, ask whether it **factorises** into an
`n_real`-carrying leg and one that does not:

- `fid = (k/n_pred) · min(1, n_pred/n_real)` — the coverage leg carries `n_real` and is
  unreadable; **the precision leg carries none and transfers.**
- `expr_mse_unbiased_capped_norm` looks like it needs the reference entirely, but written as
  `(A − 2X + Kc − ρ·claim)/sumD`, **four of its six terms come from our own file plus control
  cells we already hold** — only `X` and `sumD` needed the board. That factorisation is what made
  experiment 096 possible at all.
- `jac`'s union carries `n_real` and does **not** factorise, so its magnitude stays unreadable —
  096 quotes only its sign (+0.0062 / −0.0008, i.e. zero within noise). The reason it is
  irreducible: `n_real` sits *inside the union, added to* our own set rather than multiplying it,
  and no rearrangement separates them.

Across the four members `RESEARCH-METHOD` §1 names, the split is **one clean, one partial, one
impossible, one not applicable**: `fid` factorises cleanly; `reach` factorises *partly* (the
purity-versus-normalised-depth curve is readable — the "ordering-driven part" §1 already allows —
while the `1/n_real` normalisation is not); `jac` not at all; `nmae` is not affected by §1 in the
first place. **The rule earns its keep by telling them apart, not by rescuing all of them.**

## 1.10 The `ρ` cap — measured, and I had the sign backwards

`VERIFIED` from spec §2 and, since, measured end-to-end by experiment 096.

```
ρ = min(1,  B / Σ_q (1/G_q)·min(Ĉ_q, C_q) )
B = Σ_g Σ_{p∈R_g} (1/G_p)·( b̂_{p,g} − weighted-mean-over-p of b̂_{·,g} )²
```

**(a) The hazard was real, and half my reading of it was right.** As amplitude falls, `B → 0`,
`ρ → 0`, and the submission forfeits its `min(Ĉ, C)` credit — which does *not* shrink, because
`Ĉ` is the jackknife scatter of the *submitted cells*. The board proves the size of it: carp's
NTC-copy reads **2.7054**, not 1.0. I was right that a naive 1-D sweep's argmin is biased high and
that `ρ` must be instrumented. I was wrong that the true optimum lies below the cliff — **there is
no optimum below it** (Part 0). And my `α ≈ 1/√3 = 0.577` used the spec's "factor of about 3",
which is measured on the *reference* contexts; our own measured `B/claim` is 1.66–2.37, putting
our cliff at `VCC_ALPHA` **0.70–0.84**.

**(b) RETRACTED — I claimed reducing emission dispersion buys headroom. It is the reverse.**
My reasoning was that a smaller `Ĉ` raises `B/Σ min(Ĉ,C)` and so pushes the cliff down toward the
optimum. That is arithmetically true and strategically backwards, because it assumed we were *over*
the `#247` cap. We are **under** it: reference `C_p` = 14.54 / 14.77 / 16.37 on A/B/C against our
`Ĉ_p` = **10.30** (0.67×), and the cap fires on **0 of 900 rows**. Cutting dispersion destroys
credit we are already failing to collect. **Raising `Ĉ` toward `C_p` is the single biggest free
lever** — worth 0.349 of raw mse — *and* it raises the cliff, so less amplitude shrink is needed.

**(c) The one part of the original that stands.** `B` subtracts the per-gene mean over
perturbations before squaring, so it measures only the **perturbation-specific** spread. 096
measured the complement directly: the common mode is **0.84% of `A`**. That is what makes `B ≈ A`,
which is precisely why there is no interior optimum — and it means shrinking the shared component
is nearly free in `B` while doing almost nothing to the total energy either. The two-parameter
shrink I proposed is therefore not the lever I thought it was; the dispersion cap is.

## 1.11 The leaderboard has a public API with a `description` field we are not capturing

**Action item, highest practical value in this document.**
`https://virtualcellchallenge.org/api/leaderboard` returns all entries including a free-text
**`description`** field. **133 of 386 teams have filled it in with real method detail** —
hyperparameters, data sources, ablation results, even leaked configs. Our snapshots in
`data/leaderboard/` have `"description": null` on every row, so **our scraper is dropping the
single most informative field on the board.** Fix the scraper and re-pull.

(Our own entry's description is empty, so we are not leaking. Consider keeping it that way.)

### What the descriptions say — the three that matter

**(1) VIVAI leaked their amplitude, and it is 0.35** (`VERIFIED`, rank 6 on the 16:46 pull):

```
CALIB_S=1.00  SUB_CURVE_ALPHA=3.00  ALPHA=0.35  SUB_UNSHRINK=0.60
SUB_EMIT=calib  SUB_DROPOUT=1.0;  npred 3224
```

and that submission (`m48`, Mike Chen) reads **`mse_raw = 0.9183`, `score_mse = 0.0732`,
`pds_cosine = 0.8235`, overall 0.16518** — they cleared the floor. Their own earlier `m43`
(VIVAI Inc) reads `mse_raw = 3.8334`. **Same pipeline, same cluster, amplitude dropped, mse
raw falls from 3.83 to 0.92 while pds stays at 0.82.** That is a within-team single-variable
confirmation of Part 0, from a competitor, at the exact amplitude the geometry predicts.

Their call-set size, `npred 3224`, is directly comparable to our a10's ~4,757 — **we call ~48%
more genes than a team that is beating us.**

**(2) AgingLab_BI states outright that they trade cosine for pds** (`VERIFIED`, rank 2):

> Same pinned emission, plus joint Gram whitening of the 300 predicted rows against each other
> (alpha 1.0) at a two-way p99 match. […] **It pays cos to buy pds, the opposite of the
> constant-vector trap.**

So "Gram whitening" is a **pds** technique that *knowingly* sacrifices cosine accuracy — and
their `mse_raw` is 2.7357, the price. This closes §1.4: the AgingLab.com/AgingLab_BI difference
is that #1 does not run the whitening, or runs it without paying the cosine. **Do not chase Gram
whitening.** We already lead the board on pds; the technique buys the one thing we do not need
and costs the one thing we do.

**(3) Zhou Li (#49) published the mechanism behind the `fid` floor** (`VERIFIED`), and it is
a better explanation than anything in our own notes:

> Of the 5,601 genes it means to assert, its cells make only 2,226 significant. What separates
> the two groups is not confidence: their directional accuracy is 0.516 against 0.518, and the
> genes that fail carry the larger predicted fold change. **It is expression.** Genes that clear
> the test sit at a median 107 counts per million in control, those that fail at 38.7. At this
> depth a 38 CPM gene averages under one count per cell, so its column is almost all ties at
> zero, and a rank test cannot separate what is not there.

`INFERRED`: this explains why the whole board sits at chance on `fid` (§1.6) — the callable
universe is bounded by control expression depth, not by model quality.

**⚠ My follow-on inference — "allocate calls by control CPM, not by predicted |lfc|" — was
measured and does NOT reproduce on our submission.** 096 ran it on challenge context A: our
*called* genes sit at median control CPM **48.1** and our *uncalled* genes at **87.6** — expression
predicts callability in the *opposite* direction for us. Our callability is set by our own effect
sizes and spread policy, not by depth. Genes at control CPM ≥ 40 are only 62.7% of our called
genes, so a CPM-keyed tail policy would over-disperse more than a third of them.

Zhou Li's observation is almost certainly true *of their pipeline* — it is a well-argued
measurement — but it is a **submission-specific** property, not a property of the data, and it
does not transfer to ours. A good reminder that a competitor's published diagnosis is evidence
about their model, not about the task.

### Other useful disclosures

- **#184 "To infinity and beyond"** quantifies the noise floor: *"A 400-cell NTC draw carries
  ‖delta‖~2.79 of pure sampling noise; v1's 0.5x shrinkage put the transferred effect at SNR
  0.77, so the scorer measured mostly noise (pds_cosine 0.638, fid 0.427 below the 0.5124
  chance point)."* — an amplitude that drops below the sampling-noise SNR is measuring its own
  noise. An independent confirmation, from a competitor's board data, that the amplitude axis has
  a hard lower guardrail — consistent with 096's finding that `mse(s)` has no interior optimum
  below the `ρ` cliff.
- **#108 Qigong Li** swapped emission from bootstrap-resampled real control cells to a
  parametric gamma-Poisson with per-gene dispersion fitted from that context's NTCs, and reports
  it *"strictly dominates bootstrap on all six metrics: pds +0.035, fidelity yield +0.048,
  overlap +0.025, nmae −0.012 (better), reach +0.003, **mse flat**."* — emission family moves
  five members and **not** mse. `INFERRED`: mse is driven by the mean profile, not the emission
  family — which is what makes the amplitude fix orthogonal to the emission work.
- **#231 Toshi** explicitly avoids stacking noise: *"bootstrap w/o replacement + apply delta +
  library-size restore + stochastic rounding (no Poisson resample on top of real bootstrap
  variance)"* → `mse_raw` 1.4290.
- **#249 Biomni-TusoAI** reports *"Full 5-fold LOCO +0.2229"* locally and scores **−0.0514** on
  the board. A ~0.27 local-to-board divergence, from a competent team (Kexin Huang / Biomni,
  Alistair Turcan / TusoAI, CMU). Cautionary data point for anyone trusting a local harness.
- **#71 Sidechain (Saber Hafezqorani, UBC/BC Cancer)** publishes his entire model series openly
  at `saberhq.com/sidechain/models/` and `github.com/saberhq/sidechain`.
- **#17 OpenPert's description is a paper title**: *PertMind: Eliciting Emergent Biological
  Reasoning in LLM via Reinforcement Learning on Cellular Perturbation Data*
  (`arXiv:2608.16419`, Zhejiang University + Tencent; `github.com/shapsider/PertMind`, ~115
  stars; HF `tzcfly/PertMind`). RL over Tahoe-100M with gene/pathway/format rewards on Qwen3-4B.

### Cross-team amplitude evidence is confounded — do not over-read it

I extracted a self-reported amplitude/alpha/effect-scale from 31 descriptions and correlated it
against `mse_raw`: **r = +0.165**, i.e. essentially nothing. Median `mse_raw` is 1.75 for
α ≤ 0.5, 4.18 for 0.5 < α ≤ 1.5, 3.26 for α > 1.5. `INFERRED`: each team's α multiplies a
*different* base delta with a different norm, so α is not comparable across teams. **Only
within-team single-variable moves are evidence** — and those (VIVAI m43→m48, Aginglab v9→v12)
both point the same way.

---

## 1.12 Arc's scorer is PUBLIC and we should be reading it

`VERIFIED`. **`github.com/ArcInstitute/cell-eval2`** — MIT licensed, 13 stars, last pushed
2026-08-22. It is the actual 2026 scoring implementation: `src/cell_eval2/catalog.py`,
`scales.py`, `anchor.py`, `norm.py`, `ceiling.py`, `metrics/{de,delta,direction,discrimination}.py`,
plus `configs/vcc2026.yaml`. Our `docs/reference/vcc2026-metrics-FULL.md` is the prose spec from
`docs/vcc2026_metrics/`; **the source, the `CHANGELOG.md` and the `AUTHORS` file are additional
and are not in our tree.**

Confirmations straight from the source (`VERIFIED`):

- `catalog.py` — `expr_mse_unbiased_capped_norm` carries
  `Scoring(scored=True, direction="lower", anchor=0.0, penalty="boxcox", clamp_low=0.0,
  clamp_high=1.0)` with `agg="ratio_of_sums"` over
  `numerator="expr_mse_unbiased_capped" / denominator="expr_distance_unbiased"`. The comment
  states it is *"the ONE scored vcc2026 member that clamps to [0, 1], on BOTH scales"*, and
  *"The metric is capped so that 1.0 = predicting the control, so a value past the baseline is
  already no-skill."* Part 0 confirmed at the source.
- The same comment names the cause of reading above 1: *a prediction whose own correction
  exceeds k times the reference's — **fewer cells OR more dispersed ones** — reads above 1,
  because #247's cap refuses a correction the reference does not earn.*
  **⚠ My inference from this sentence was backwards, corrected by measurement.** I read it as
  "our emission is *over*-dispersed and that is inflating our mse". 096 measured the opposite:
  our `Ĉ_p` is **0.67× the reference `C_p`** and the cap fires on **0 of 900 rows** — we are
  *under*-dispersed and forfeiting credit. The clause bites submissions above the cap; we are
  below it. The correct reading of our 0.9640 → 1.3627 regression is amplitude
  (`VCC_ALPHA` 0.40 → 1.0796) crossing the `ρ` cliff, not dispersion.

  A related disclosure in the same source explains why *our* emitter gets a better deal than most
  of the board: the exact-total repair **pins our group aggregates**, so our emission noise does
  not land in `A` and we can claim `min(Ĉ, C)` for free. The teams whose realized constants read
  1.16–2.71 (§Part 0) do not pin theirs and pay for it. That is the `#348` residual the scorer's
  own docstring discloses.

`AUTHORS` records externally-reported issues and their rulings (`VERIFIED`). Two are directly
relevant:

- **#278** — *"`expr_mse_unbiased_capped_norm` still reads within-group cell dispersion."*
  Ruled docs-only at the time; residual recorded as **#294** and parked until after the
  competition; *"re-measured a week later on a live submission as the largest single member
  contribution in the scored set, and the measured exploit was closed in #348 (PR #353);
  **#294's general class stays open.**"*
- **#279** — *"`direction_reach_raw`'s no-skill point is a function of the reference's `N_conf`
  distribution, not the constant the scale registry assumed, and `k*` is decided at the head of
  a ranking the submission controls, which makes it **exploitable** rather than merely
  miscalibrated."* Measured: a no-skill arm collects `reach_raw` 0.042 / 0.080 / 0.092 on the
  three contexts. *"The exploit half is open."*

`INFERRED`: #279 says `reach` is decided by the **head** of a ranking we control, and that a
no-skill arm already collects 0.042–0.092 (our raw is 0.1890). Combined with Zhou Li's finding
that callability is set by control CPM, the `reach` lever is *ordering*, and the ordering should
put high-control-CPM genes first. That is consistent with exp094/095's direction and gives it a
mechanism.

**Recommendation: vendor or track `ArcInstitute/cell-eval2` and read `metrics/delta.py`
(`mse_unbiased_capped`) directly rather than working from the prose spec.** It is the scorer.

---
# PART 2 — the published state of the art

## 2.1 The headline: "you cannot beat the mean" is the field's central result

This literature is large, recent, and directly load-bearing. Confidence per row.

| paper | venue | finding |
|---|---|---|
| Ahlmann-Eltze, Huber, Anders, *Deep-learning-based gene perturbation effect prediction does not yet outperform simple linear baselines* | **Nat Methods 22(8):1657–1661 (2025)**, [doi:10.1038/s41592-025-02772-6](https://www.nature.com/articles/s41592-025-02772-6) | 5 foundation models + GEARS + CPA vs "no change"/"additive"/"mean"/linear. **None consistently beat the mean.** `VERIFIED` |
| Kernfeld, Yang, Weinstock, Battle, Cahan, *A comparison of computational methods for expression forecasting* | **Genome Biology 26 (2025)**, [doi:10.1186/s13059-025-03840-y](https://doi.org/10.1186/s13059-025-03840-y) | 11 datasets, 9 methods. *"The mean or median baseline was almost always the top performer."* Best mechanism section in the field. `VERIFIED` |
| Csendes, Sanz, Szalay, Szalai, *Benchmarking foundation cell models for post-perturbation RNA-seq prediction* | **BMC Genomics 26:393 (2025)**, [PMC12016270](https://europepmc.org/articles/PMC12016270) | Train-Mean beats scGPT/scFoundation on Pearson-delta on 4/4 datasets. `VERIFIED` |
| Wong, Hill, Moccia, *Simple controls exceed best deep learning algorithms…* | **Bioinformatics 41(6):btaf317 (2025)**, [PMC12202205](https://pmc.ncbi.nlm.nih.gov/articles/PMC12202205/) | The "CRISPR-informed mean". `VERIFIED` |
| Wei et al., **scPerturBench** | **Nat Methods 23(2):451–464 (2026)**, [doi:10.1038/s41592-025-02980-0](https://doi.org/10.1038/s41592-025-02980-0), code [bm2-lab/scPerturBench](https://github.com/bm2-lab/scPerturBench) | 27 methods / 29 datasets. Numbers below recomputed from the authors' raw CSVs. `VERIFIED (recomputed)` |
| Wu et al., **PerturBench** | [arXiv:2408.10609](https://arxiv.org/abs/2408.10609) | Rank metrics + **mode-collapse** diagnostics; *"simpler architectures are generally competitive"*. `SECONDHAND` |
| Wenteler et al., **PertEval-scFM** | ICML 2025 / PMLR v267, [bioRxiv](https://www.biorxiv.org/content/10.1101/2024.10.02.616248v2.full.pdf) | On Replogle K562 the *best* scFM beats the mean by **0.16% relative**; UCE and Geneformer are *worse* than the mean on Norman. `VERIFIED` |
| Heidari et al., *Evaluating Single-Cell Perturbation Response Models Is Far from Straightforward* | [bioRxiv 2026.02.14.705879](https://www.biorxiv.org/content/10.64898/2026.02.14.705879v1.full) | Replicate Pearson-delta **ceiling is only 0.4–0.7**, sometimes <0.3. `VERIFIED` |
| Vollenweider & Bühlmann, *Signal, Bounds, and Baselines* | bioRxiv 2026.04.20.719650 | 7 datasets; DL *"often still fails to meaningfully surpass simple linear baselines"*. `SECONDHAND` |
| Mao et al., *Benchmarking virtual cell models for in-the-wild perturbation response* | [arXiv:2604.27646](https://arxiv.org/abs/2604.27646) | Embedding-based (non-random) splits collapse performance. `VERIFIED` |
| *A Systematic Comparison…* (13 methods, 25 datasets, 24 metrics) | [bioRxiv 2024.12.23.630036](https://www.biorxiv.org/content/10.1101/2024.12.23.630036) | *"Models shared a **conservative bias**, with fine-tuned foundation models compressing variance and underestimating synergistic effects."* And: **"no method consistently generalized across cell types."** `VERIFIED` |

## 2.2 WHY the mean wins — the mechanism, not the headline

Kernfeld et al. name two non-exclusive causes (`VERIFIED`, quoted from the paper):

**(a) Bias–variance. This is the dominant cause, and it is statistical, not biological.**

> The mean of the control samples is unbiased, but given limited sample sizes, it has high
> variance. The mean of all training data will be biased by perturbation effects but will have
> lower variance. If perturbation effects are generally small and noise levels generally high,
> the mean of all samples may be a better estimator.

Their control experiment is the one to remember: **in simulated data containing only noise and
no true perturbation effect at all, the mean baseline yields observed-vs-predicted log-fold-change
correlations > 0.6.** A Pearson-delta of 0.6 is consistent with zero signal.

**(b) A real shared biological axis.** They confirm that LFC(A/ctrl_A) correlates with
LFC(B/ctrl_B) for *different* perturbations with *independent* controls, "far from zero" on four
datasets, and Enrichr identifies the shared axis as **stress / p53 / apoptosis / differentiation**
(p53 effectors on replogle4; erythroid differentiation — 3 hemoglobin subunits, 2 glycophorins,
SLC25A37, BLRVB — on replogle2/K562). Their conclusion: *"by learning the mean of the training
data, perturbation prediction models could appear to display some level of biological insight,
even without including any correct mechanisms."*

### The quantitative picture

**How much is shared** — median pairwise correlation between DE profiles of *different*
perturbations (Csendes Fig 2D, `VERIFIED`): **Adamson (UPR-focused) 0.662**, **Replogle K562
(genome-wide) 0.117**. Replogle's own background distribution of gene-pair profile correlations
is median r = 0.10 (vs 0.61 within a CORUM complex). Csendes show an **inverse relationship**
between a dataset's median intra-correlation and the headroom any model has over Train-Mean.
`INFERRED`: VCC 2026 is a curated 300-gene panel, so its shared fraction is probably between
these two poles — closer to Replogle than Adamson, but not as low.

**How many perturbations do nothing** — Replogle 2022 K562 genome-wide (`VERIFIED`): permuted
energy-distance test finds **2,987 / 9,608 (31.1%)** perturbations produce a significant
transcriptional phenotype (controls 1.9%); Anderson–Darling with >10 DEGs gives 30.5%.
**≈69% of perturbations are transcriptomically silent.** Upstream cause: CRISPRi knockdown often
just fails — only ~40–50% of targeted genes are effectively knocked down.

**How small the real effects are** — TRADE (Nadig et al., *Nat Genet* 57:1228–1237 (2025),
[PMC11244993](https://pmc.ncbi.nlm.nih.gov/articles/PMC11244993/), `VERIFIED`): a typical
perturbation affects **~45 genes** (essential genes >500) out of ~18,000 — everything else *is*
the mean. FDR-significant genes capture only 13–36% of true transcriptome-wide impact. And the
number that reframes everything: **median correlation of raw log-fold-change point estimates
between two replicate K562 experiments is 0.16; after deconvolving sampling noise, 0.90.**
Measured deltas are almost all noise; the truth underneath is highly reproducible.

**Is it a signal, metric, or generalisation problem?** All three — but for **MSE specifically the
metric result settles it** (Miller et al., `VERIFIED`):

> On Replogle22, the technical duplicate — a real re-measurement of the same perturbation,
> i.e. the physical ceiling of any predictor — scores **worse than the mean baseline on MSE**
> for **~95% of perturbations.**

Because the duplicate is built from half the cells and carries more sampling noise than a
population mean: *"Metrics that treat all genes uniformly reward this noise reduction, even though
it is biologically irrelevant."* They name the two artifacts **control bias** and **signal
dilution** (setting the top-25 perturbation-specific DEGs — ≈0.3% of genes — exactly equal to
truth barely moves unweighted MSE).

**And the sentence that describes our exact situation** (Ahlmann-Eltze, `VERIFIED`):

> Unlike the L2 distance, the Pearson delta metric does not penalize predictions that are
> consistently too small or too large in amplitude and, thus, prioritizes correct prediction of
> the direction of the expression change.

Substitute `pds_cosine` for "Pearson delta" and `expr_mse` for "L2 distance" and that is
§Part 0, written a year early. Note the field's *typical* error is **under**-shooting (models
collapse toward the mean under MSE training); ours is over-shooting. We have the rarer problem,
and it is the easier one.

## 2.3 What the few methods that DO beat the baselines are doing

Five mechanisms; only two are architectural.

1. **Use the perturbation's identity directly.** Wong et al.'s **CRISPR-informed mean**: predict
   the training mean of perturbed cells for every gene *except* the target, set the target to 0.
   Beats GEARS by 0.08 Pearson-delta (p = 9.3e-4) and scGPT by 0.11 (p = 8.1e-6); on the rank
   metric **4.7–213.9× better than GEARS**. `VERIFIED`.
   **⚠ Closed off for us:** VCC 2026 excludes each perturbation's own target gene from both legs
   of every member. Arc quantified the exploit: an adversary predicting control everywhere plus
   truth at its target *"was taking 6–11% of the scaled range before the exclusion was applied,
   because that gene is the largest-moving one in 57–66% of perturbations."* `VERIFIED` (spec).
2. **Cross-context transfer through a latent basis.** Ahlmann-Eltze's *only* consistent winner:
   a linear model `argmin_W ‖Y − (GWPᵀ + b)‖²` with the perturbation embedding **P pretrained on
   the other cell line**. *"Pretraining on the single-cell atlas data provided only a small
   benefit over random embeddings, but pretraining on perturbation data increased predictive
   performance."* `VERIFIED`.
3. **Functional priors over gene identity.** RF + Gene Ontology beats every foundation model
   (0.739/0.586/0.480/0.648 Pearson-delta vs Train-Mean's 0.711/0.557/0.373/0.628 on
   Adamson/Norman/Replogle-K562/RPE1). A plain **k-NN over a biological knowledge graph "beats
   almost all methods on out-of-distribution perturbation prediction"**
   ([arXiv:2606.08816](https://arxiv.org/abs/2606.08816)). `VERIFIED`.
4. **Weighted loss to stop mode collapse.** WMSE-trained GEARS: test WMSE improves **2.5×**,
   median R²_w(Δ) goes **−0.068 → 0.555**. `WMSE = Σ wᵢ(μᵢ − μ̂ᵢ)²`, weights from |t-score vs
   *Rest*| (not vs control), min-max normalised, squared, summed to 1
   ([arXiv:2506.22641](https://arxiv.org/abs/2506.22641)). `VERIFIED`.
5. **Model the distribution, not the mean** — STATE, PRiMeFlow. See §2.5; the re-benchmarking is
   unkind.

## 2.4 Cross-cell-line vs unseen-perturbation — and it matters which one we are

Recomputed from scPerturBench's raw tables (`VERIFIED, recomputed`):

| scenario | trainMean baseline | best method | **headroom** |
|---|---:|---:|---:|
| **Unseen context, known perturbation** (top-100 DEG) — **our task** | 0.597 | trVAE 0.631 | **+0.034** |
| Unseen perturbation, known context (top-100 DEG) | 0.402 | CPA 0.638 | **+0.236** |
| Unseen perturbation, all 5,000 genes | 0.264 | CPA 0.497 | +0.233 |

**The exploitable headroom over a trivial mean is ~7× smaller in our setting.** Only 3 of 13
methods beat trainMean cross-context; on 2 of 12 datasets trainMean beats everything. True o.o.d.
cross-context costs the average published method −17.6% vs an i.i.d. split — and biolord /
SCREEN / baseMLP / baseReg don't degrade at all, **because they were never using context
information in the first place.**

**The ceiling**, from TRADE — median transcriptome-wide-impact correlation for the same
perturbation across lines, 1,660 perturbations significant in all four (`VERIFIED`):

| HepG2↔RPE1 | K562↔Jurkat | Jurkat↔HepG2 | K562↔HepG2 | Jurkat↔RPE1 | K562↔RPE1 |
|---:|---:|---:|---:|---:|---:|
| 0.75 | 0.74 | 0.69 | 0.64 | 0.56 | **0.40** |

Structure follows lineage (K562+Jurkat hematopoietic/p53-mutant/suspension; HepG2+RPE1
epithelial/p53-WT/adherent). A two-component mixture splits perturbations **56% broadly-conserved
/ 44% cell-type-preferential**. Qi & Chapfuwa (bioRxiv 2026.05.08.723625) prove the condition:
effects transport **iff the columns of the causal influence matrix for the perturbed genes are
invariant** — distributional similarity is *not* sufficient. `VERIFIED`.

`INFERRED`, and it matters: **the board's flatness is not our failure, it is the task.** The
top-16 spread 0.154→0.211 against a Frankenstein ceiling of 0.237 is exactly what a +0.034-headroom
task looks like. This is an argument for harvesting the *mechanical* members (mse, nmae) rather
than trying to out-model the field on biology.

## 2.5 Method-by-method: what each would actually buy us

### Arc's own STATE — `VERIFIED` from source, not the paper
[`github.com/ArcInstitute/state`](https://github.com/ArcInstitute/state), **CC BY-NC-SA 4.0**
(non-commercial), 651 stars. Paper [bioRxiv 2025.06.26.661135v2](https://www.biorxiv.org/content/10.1101/2025.06.26.661135v2), *Cell* 2026 in press.

- **Outputs a set of cells, not a mean.** `StateTransitionPerturbationModel` consumes a set
  (`cell_set_len: 512`) and emits a matched predicted set. `predict_residual: True` → basal +
  delta. A `FinetuneVCICountsDecoder` maps latent → counts with a **learnable `read_depth`**.
- **Magnitude:** yes, structurally. Default loss is `energy` (`geomloss.SamplesLoss`), i.e.
  energy distance / MMD between predicted and observed cell sets — not scale-invariant. But it
  penalises magnitude *jointly with spread*, in whatever space `embed_key` selects, which is a
  much weaker constraint than squared error on the scored log-pseudobulk profile.
- **Unseen contexts:** the paper's own stated limitation is *"we have not tested its usability on
  entirely held-out datasets for which no cell contexts have been seen during training"* — which
  is precisely the 2026 task.
- **The baselines ship in the repo**: `perturb_mean.py`, `context_mean.py`, plus
  `celltypemean.yaml` / `pseudobulk.yaml` / `globalsimplesum.yaml`.
- **No MSE/MAE-type headline metric appears in STATE's own results.**

**Verdict: STATE helps direction members. It is not a magnitude-calibration tool, and it is
non-commercially licensed.** Consistent with the board: every self-declared naive STATE run sits
at the bottom (§1.8).

### The three papers that are actually our problem, solved

**★ COMPASS — Liang & Singh (Duke), [bioRxiv 10.64898/2026.08.03.742643](https://www.biorxiv.org/content/10.64898/2026.08.03.742643v1), 6 Aug 2026.** `VERIFIED`, full text.
The model *is* the missing piece:

```
z_cp  =  β_cp · u_c  +  r_cp
```

`u_c` = unit vector along cell line c's **mean perturbation direction**; `β_cp` = a **scalar
magnitude coefficient per perturbation**; `r_cp` = gene-specific residual, orthogonal to `u_c`
by construction. Evaluated on **2,270 CRISPRi perturbations measured in each of six cell lines**
(Replogle K562/RPE1, Nadig HepG2/Jurkat, X-Atlas/Orion HCT116/HEK293T) — **the same six lines as
VCC 2026.**

**The number that should change our plan: leave-one-cell-line-out, the shared-response
coefficient β transfers at Pearson r = 0.66, while the gene-specific residual transfers at
r = 0.22.** Their sentence: *"the shared-response coefficient thus transfers considerably more
reliably than the gene-specific residual."* Cross-context: COMPASS-H 0.34 de-biased Pearson-Δ /
0.23 cosine-PDS gain, vs GEARS ≤0.28/0.01, scGPT ≤0.28/0.00, training-mean 0.19–0.32/0.00.
Code promised at `github.com/rohitsinghlab/compass`.

**Read-across (`INFERRED`): in a held-out cell line, magnitude is ~3× more transferable than
directional detail — and we are currently monetising only the hard half.** This also gives the
two-parameter decomposition a principled form — split `β` (shared magnitude) from `r`
(gene-specific residual) and treat them separately. Note §1.10(c): the shared component is only
0.84% of our prediction energy, so a `β`/`r` split is a *direction* lever here, not the mse lever
I originally took it for.

**★ Miller et al. (Shift Bioscience), [bioRxiv 2025.10.20.683304](https://www.biorxiv.org/content/10.1101/2025.10.20.683304v1.full).** `VERIFIED`.
Besides the ~95% duplicate-loses-to-mean result, they give **a published per-gene shrinkage
formula** — the literature's direct answer to "match the mean while keeping perturbation-specific
direction":

```
μ_ID  =  α ⊙ μ_TD  +  (1 − α) ⊙ μ_all       with   α = 1 − p_DEGs
```

Per gene: weight your perturbation-specific prediction by (1 − its DE p-value), the context mean
by the p-value. Genes you are confident about keep your signal; the other ~99% fall back to the
mean and **inherit the mean baseline's MSE exactly**. This is DESeq2's moderated-LFC shrinkage in
different clothing, and it is **strictly better than a single global scalar**, because the
MSE-optimal scale differs per gene by confidence. `INFERRED`: this is the highest-value single
idea in this document after Part 0 itself.

**★ Shoeibi & Yousefi, [arXiv:2608.00152](https://arxiv.org/abs/2608.00152), 31 Jul 2026 — run on VCC H1 data.** `VERIFIED`.
Four scalar magnitude features of the input (`m₁` mean|x|, `m₂` RMS, `m₃` ℓ∞, `m₄` std, z-scored)
predict the log Anderson–Darling distance from NTCs:

| model | R² |
|---|---:|
| RF on `[x; m̃]` | **+0.370 ± 0.004** |
| **linear regression on the 4 scalars alone** | **+0.208** |
| strongest `x`-only classical model | +0.190 |
| deep MLP encoder | +0.082 |
| **direction-only (row-normalised `x`)** | **+0.062** |
| shuffled-magnitude control | +0.071 |

**Magnitude-only predictors transfer zero-shot to Replogle RPE1 and K562-essential; direction
predictors do not.** Row-normalising collapses the best model from +0.190 to +0.062.

### Direction and magnitude are nearly the same variable
Raju, [arXiv:2604.16642](https://arxiv.org/abs/2604.16642) (2026), `VERIFIED`. Per-perturbation
coherence `S_p = (1/n_p) Σᵢ (dᵢ·d̄)/(‖dᵢ‖‖d̄‖)` vs effect magnitude, Spearman ρ:

| Norman | Adamson | Dixit | Papalexi | **Replogle** | pooled |
|---:|---:|---:|---:|---:|---:|
| 0.950 | 0.951 | 0.843 | 0.945 | **0.977** (n=1,832) | **0.974** (n=2,285) |

*"Directional coherence is largely a restatement of effect magnitude."* `INFERRED`, and directly
usable: **our per-perturbation cosine confidence is already a magnitude estimator.** We have the
input to a per-perturbation β without training anything.

### Why scaling up helps discrimination and destroys MSE — the derivation
Liu, Zhang, Du, Zhao, Wang (**VCC 2025 3rd place, team Outlier**),
[arXiv:2511.16954](https://arxiv.org/abs/2511.16954). `VERIFIED`. For a scaled prediction `cΔ̂ᵢ`,
the L2 ranking difference is `(‖Δᵢ‖² − ‖Δᵢ'‖²) − 2c(Δ̂ᵢᵀΔᵢ − Δ̂ᵢᵀΔᵢ')`, so **as c→∞ the L2-PDS
converges to cosine and the L1-PDS to a magnitude-weighted sign-cosine** — discrimination rises
to a plateau while MSE grows without bound. On identical predictions they measure L1/L2 PDS ≈ 0.50
(chance) vs cosine PDS ≈ 0.80, and **norm-matching does not fix it**. Also: *"Total effect
magnitudes are intrinsically difficult to predict: they depend on guide RNA efficiency,
experimental design, noise levels… and even subtle preprocessing choices."*

**Their entire 3rd-place contribution was a global linear scaling factor selected by CV.**

### The rest, compressed

| method | outputs | magnitude? | unseen context? | source / licence |
|---|---|---|---|---|
| **X-Cell / X-Cell-Ultra** (Xaira, 4.9B, diffusion LM, X-Atlas/Pisces 25.6M cells) | control→perturbed refinement; headline is Pearson-Δ on logFC | **direction-first**, no MSE/MAE reported | yes, zero-shot + test-time adaptation | [bioRxiv](https://www.biorxiv.org/content/10.64898/2026.03.18.712807v1); `Xaira-Therapeutics/X-Cell`, CC BY-NC-SA |
| **PRiMeFlow** (Altos; VCC-2025 Generalist winner) | **single cells in count space** | **yes, explicitly** — inference post-processing is `expm1 → median-normalise → round → median-normalise again` | pretrained atlas, fine-tuned per context | [arXiv:2604.13986](https://arxiv.org/abs/2604.13986); [`altoslabs/primeflow`](https://github.com/altoslabs/primeflow), academic non-commercial |
| **PerturbMap** | source→recipient ridge experts + local low-rank base | **explicitly calibrates magnitude**; raw-copy vs calibrated-copy is a control; scored on full-effect MSE | that is the whole paper | [arXiv:2607.28090](https://arxiv.org/abs/2607.28090) |
| **Response decomposition** (Molina & Zhang, Altos) | splits global / pert-specific / line-specific / interaction | **the global component is low-dimensional and inferable from control expression**; specific components are not | yes, unseen lines; **simple linear/MLP beat SOTA architectures** once aligned to the right component | [bioRxiv](https://www.biorxiv.org/content/10.64898/2026.07.24.740459v1) |
| **SLIM** | **640 trainable params.** Bilinear extension of Ahlmann-Eltze + 64-d STRING, closed-form ridge → predicted **mean**, then retrieves training cells and **rescales each gene to match the predicted mean** | **yes, by construction** | within-dataset only | [bioRxiv](https://www.biorxiv.org/content/10.64898/2026.08.07.743481v1); `RasmussenLab/SLIM`, MIT |
| **TxPert** | multi-KG latent transfer (STRING, GO, Recursion PxMap) | direction-oriented | **explicitly tested on unseen lines**; approaches split-half reproducibility | [Nat Biotech 2026](https://www.nature.com/articles/s41587-026-03113-4); `valence-labs/TxPert` |
| **Lingshu-Cell** | masked discrete diffusion over quantised expression tokens (281 bins), all ~18k genes | yes, via the token grid | not demonstrated | [arXiv:2603.25240](https://arxiv.org/abs/2603.25240) |
| **scDFM** | distributional flow matching | MSE on velocity + optional multi-σ MMD | population-level, not context transfer | `AI4Science-WestlakeU/scDFM`, MIT |
| **Tabular FMs (TabICL/TabPFN)** | in-context pseudobulk regression, no biology-specific architecture | pseudobulk | cross-cell-type in-context | Palla/Royer (CZI), [bioRxiv](https://www.biorxiv.org/content/10.64898/2026.06.28.735106v1) — *"consistently out-perform specialized baselines"* |
| **GEARS / CPA / chemCPA / biolord / scGen / CellOT / PRnet / sams-vae** | latent shift or mean profile | mostly mean-level, MSE-trained → **under**-disperse | mostly unseen-perturbation, not unseen-context | GEARS MIT, CPA BSD-3, scGen GPL-3, CellOT BSD-3 |
| **scGPT / Geneformer / scFoundation / scBERT / UCE / Nicheformer / C2S** | **embeddings only** | **no** | n/a | scGPT MIT, scFoundation Apache-2.0, scBERT GPL-3, UCE MIT |

`INFERRED`, and worth stating plainly: **every foundation model outputs an embedding, so the
prediction — and therefore all the magnitude calibration — lives in a decoder head you supply.
None of them helps `expr_mse` on its own.** The board agrees (§1.8).

Counter-evidence worth reading before dismissing FMs: GenBio AI's *Foundation Models Improve
Perturbation Response Prediction* ([bioRxiv](https://www.biorxiv.org/content/10.64898/2026.02.18.706454v1),
`genbio-ai/foundation-models-perturbation`) benchmarks 600+ models, finds many FMs beat PCA and
that *"many 'advanced' perturbation prediction methods (e.g. GEARS, latent diffusion, flow
matching) underperform FMs"* — **and it ships a `cross_context/` directory with direct comparisons
to STATE in our exact setting.**

## 2.6 Which member each family buys

| family | `expr_mse` | `de_wilcoxon_lfc_nmae` | direction members (`pds`, `fid`, `reach`, `jac`) |
|---|---|---|---|
| Cosine/rank-trained, KG methods, FM embeddings + Pearson-Δ heads, X-Cell | **no** (scale-invariant objectives) | weak | **yes** — their whole value |
| STATE (energy/Sinkhorn on cell sets) | partial (not in the scored space) | partial | yes |
| Generative-in-count-space with renormalisation: PRiMeFlow, Lingshu-Cell, eVAE, scDFM | **yes** — the field's biggest MAE wins come from here, via **post-processing not architecture** | yes | yes, but they trade along the eVAE "distribution richness" axis |
| **Explicit magnitude decomposition: COMPASS, PerturbMap, response-decomposition, SLIM, Miller's interpolated duplicate** | **yes, by construction and cheaply** | **yes** | preserved |
| Mean/context-mean baselines | defines the anchor | ~0 | chance |

**One more metric-design result worth internalising** (Bhattacharya et al., eVAE,
[bioRxiv](https://www.biorxiv.org/content/10.64898/2026.06.15.732063v1), `VERIFIED`): a decoder-head
ablation exposes **one design axis — the richness of the predictive distribution at inference —
that splits standard metrics into two groups, "variance-sensitive and mean-sensitive, which move
in opposite directions along the axis."** That is exactly our `pds_cosine`-vs-`expr_mse` tension,
isolated to a single tunable knob, and it is the same knob as §1.10's emission dispersion.

---
# PART 3 — the 2025 challenge, and what actually transfers

## 3.1 Final results

Recovered from the archived leaderboard JSON (Wayback snapshot `20260207134815` of
`virtualcellchallenge.org/leaderboard`, 125 entries), matching Arc's announced winners.
`VERIFIED`. Winners announced 6 Dec 2025 at NeurIPS; 5,000+ registrants, 1,200+ teams submitting,
300+ final.

| # | team | model | score | DES | PDS | **MAE** | prize |
|---:|---|---|---:|---:|---:|---:|---|
| 1 | BM_xTVC (BioMap) | xTrimoSCPerturb | 0.3401 | 0.3488 | 0.8720 | 1.026 | $100k |
| 2 | XLearning Lab (Xi Peng, Sichuan U.) | X | 0.3279 | 0.3555 | 0.8505 | 1.069 | $50k |
| 3 | Outlier (UChicago/Dartmouth/HKU) | TransPert | 0.3259 | 0.3611 | 0.8444 | **4.219** | $25k |
| 4 | CaiLab | PertsPredictor | 0.3115 | 0.2669 | 0.8744 | 1.114 | — |
| 6 | SIMON | **"Naive"** | 0.3053 | 0.3547 | 0.8175 | 1.111 | — |
| Generalist | **cleopatra / Altos Labs** | "go with the flow" = **PRiMeFlow** | 0.2869 (avg rank) | — | — | **0.086** | $100k |

Note **a model named "Naive" finished 6th**.

Methods (`VERIFIED`/`SECONDHAND` as marked):
- **BioMap (1st)**: scFoundation-derived + xTrimo protein embeddings + **explicit DEG-frequency
  and mean-expression features**; **fused loss = PDS + DES + small-weight MAE**. No code, no
  paper — LinkedIn announcement only. `SECONDHAND`.
- **XLearning Lab (2nd)**: FCN on **pseudobulk**, ESM-2 embeddings, UMI-count indicator;
  **residual learning of the delta with MSE loss**. Nothing else findable. `SECONDHAND`.
- **Outlier (3rd)**: **pure statistics** — pseudobulk log-fold-change matrix + Wilcoxon p-values
  only, plus **a global linear scaling factor tuned by CV**. Their own slide: *"Simply using
  observed effect in other cell lines as predictions… already outperform existing approaches."*
  Theory published as [arXiv:2511.16954](https://arxiv.org/abs/2511.16954). Full TransPert
  manuscript still forthcoming. `VERIFIED`.
- **Altos (Generalist)**: **code is open** — [`altoslabs/primeflow`](https://github.com/altoslabs/primeflow),
  paper [arXiv:2604.13986](https://arxiv.org/abs/2604.13986). Flow-matching in full gene space,
  U-Net velocity field, classifier guidance, pretrained on ~7M cells. README has an explicit
  *"Reproduce PRiMeFlow Results in the Arc Virtual Cell Challenge"* section. `VERIFIED`.
- **Turbine AI (15th)** published a [full post-mortem](https://blog.turbine.ai/p/how-did-we-get-a-regression-model):
  plain regression, no pretraining, external Perturb-seq (Replogle, Nadig) as features —
  *"a simple model first described in 1970… ultimately landed 15th."* `VERIFIED`.

## 3.2 The organizer post-mortem

[Arc's wrap-up](https://arcinstitute.org/news/virtual-cell-challenge-2025-wrap-up), `VERIFIED`:

> **"Almost all models performed worse than baseline on MAE… this effectively excluded MAE from
> competitive optimization."**

> "hybrid models combining deep learning with statistical features outperformed pure neural
> networks, and strategic loss function design mattered as much as architecture."

> "Purely AI-based approaches did not consistently outperform statistical baselines on DES and MAE."

Arc also acknowledged the community finding that *"L1-based PDS is inherently scale-sensitive."*

**Our failure mode is the field's failure mode, at scale — and it repeated.** In 2025 nearly the
whole field was worse than the mean baseline on the expression-error member. In 2026,
**`score_mse` is exactly 0 for 376 of 385 teams (97.7%)**. Same failure, relocated.

## 3.3 Baselines and the 2025 exploit

Arc never published official baselines either year. Ingenix.ai
[measured them](https://ingenix.substack.com/p/the-virtual-cell-challenge-towards-80b) (`SECONDHAND`):

| | mean baseline | replicate ceiling | DepMap-similarity baseline |
|---|---:|---:|---:|
| DES | 0.106 | 0.427 | 0.22 |
| PDS | 0.514 | 0.943 | 0.62 |
| MAE | 0.027 | 0.010 | 0.04 |

BM_xTVC captured 83% of PDS headroom and 76% of DES headroom — but its MAE was **38× worse than
the mean baseline**; Outlier's was **156× worse**. **Altos was the only top entry with a sane MAE
(0.086, 3× baseline) — which is precisely why it won the well-roundedness prize.**

**And 2025 was openly gamed.** A team submitted a model literally titled *"Leaderboard Hacker"*
and hit rank 3 on the validation board with MAE 14.13; a model named *"random-baseline"* ranked
14. The [exploit writeup](https://gmdbioinformatics.substack.com/p/arc-virtual-cell-challenge-has-the):
variance inflation + log1p transforms multiplied scores, and *"applying these transformations to
random data returns better scores than using the top models."* Root cause: the MAE penalty was
**capped at zero** — Arc confirmed *"submissions were not penalized for scoring worse than the
mean."* `VERIFIED`.

## 3.4 Why the metrics changed, and what transfers

Arc's stated rationale ([2026 announcement](https://arcinstitute.org/news/virtual-cell-challenge-2026), `VERIFIED`):
*"no single metric captures model quality, and… a scoring function with a narrow surface invites
optimization against the metric rather than the biology."* `cell-eval2`'s docs give a 1:1 swap map:
`pds_l1 → pds_cosine`, `expr_mae → expr_mse_unbiased_capped_norm`,
`de_wilcoxon_overlap → de_wilcoxon_sig_jaccard`, plus three net-new DE members.

**Does NOT transfer — actively inverted:**
- **Magnitude inflation to boost PDS.** L1-PDS asymptotically approaches sign-cosine as you scale
  up (Liu et al. measured L1-PDS ≈ 0.50 vs cosine-PDS ≈ 0.80 on identical predictions).
  **`pds_cosine` is scale-invariant by construction — that entire attack surface is gone.**
  Turbine's "magnitude matters more than direction" is now false. And the *cost* of inflation
  moved from free to fatal, because `expr_mse` is quadratic in scale.
- **The log1p / variance-inflation hack.** Closed in `cell-eval2` #348.
- **DES padding.** 2025 DES rewarded predicting many common DEGs; 2026 `sig_jaccard` is symmetric
  and **punishes** padding (Arc measured a skilled arm losing 0.26 of span for ten extra calls per
  target). `INFERRED`: this is a direct argument against a10's ~4,757-gene call set versus VIVAI's
  3,224.
- **The target-gene gift.** Excluded from both legs in 2026; was worth 6–11% of the scaled range.

**DOES transfer:**
- **Hybrid > pure end-to-end.** Stated independently by Arc, BioMap and Turbine.
- **Cross-cell-line effect transfer as a strong prior.** Outlier's core result — reusing observed
  effects from other lines beat foundation models — maps directly onto 2026's zero-shot
  multi-context framing, and is what most of the 2026 mid-board is doing (§1.8).
- **Summary-level / pseudobulk sufficiency.** Outlier and Turbine both won on pseudobulk +
  Wilcoxon. 2026 is **more** Wilcoxon-centric (4 of 6 members), so this gets stronger.
- **A CV-tuned global scale factor was worth 3rd place.** Still true, and now worth more.

**⚠ One 2025 lesson that must be *un*learned:** "the expression member is a free-fire zone, you
cannot lose on it." That was true in 2025 (clamped at zero, everyone below baseline). It is *half*
true in 2026 — `expr_mse` is still clamped at 0, so being at 1.36 costs the same as 1.001 — **but
`de_wilcoxon_lfc_nmae` is floored at −6, not 0, so magnitude inflation now bleeds there,
uncapped and invisibly.** Our nmae raw is 0.9442 against a 1.0005 baseline; the leader is at
0.8892. `VERIFIED` from the board.

---

# PART 4 — AI-for-science method

Framing note: `docs/RESEARCH-METHOD.md` is already a strong, hard-won version of this. Its eight
rules independently rediscover much of what follows (pre-registered predictions with falsifiers,
"a null whose mechanism was never verified is not a result", the harness-confound and oracle-leak
incidents, "understand the metric's shape before optimising it"). **What follows is what to *add*,
not a replacement.** The three highest-value additions are marked ★.

## 4.1 ★ AIDE — the ML experiment loop, MIT, ~700 lines
[`WecoAI/aideml`](https://github.com/WecoAI/aideml) — **MIT**, © 2024 Weco AI. Files:
`aide/agent.py`, `aide/journal.py`, `aide/utils/config.yaml`.

**The search policy** (`Agent.search_policy`), copyable as logic:

```
if len(journal.draft_nodes) < num_drafts:      return None            # draft a fresh root
if random() < debug_prob:
    debuggable = [n for n in journal.buggy_nodes
                  if n.is_leaf and n.debug_depth <= max_debug_depth]
    if debuggable:                             return choice(debuggable)
if not journal.good_nodes:                     return None            # back to drafting
return journal.get_best_node()                                        # greedy
```
Defaults: `steps: 20`, `k_fold_validation: 5`, `num_drafts: 5`, `debug_prob: 0.5`,
`max_debug_depth: 3`.

**Three rules worth adopting outright:**

1. **No metric printed ⇒ buggy ⇒ worst possible score.** In code, not prose:
   `node.is_buggy = response["is_bug"] or node.exc_type is not None or response["metric"] is None`,
   then `if node.is_buggy: node.metric = WorstMetricValue()`. A run that does not emit a number
   can never win.
2. **Metric direction is fixed once, globally.** One call at campaign start picks the metric name
   and `lower_is_better`; it is then injected into every subsequent task description, and
   `_warn_on_direction_conflict` logs any node that disagrees. Kills mid-campaign metric drift —
   directly relevant given six members that move in different directions.
3. **Atomicity.** Every improvement proposes exactly one change, and it must be atomic so its
   effect is experimentally attributable. (Already `RESEARCH-METHOD` §2 in spirit; AIDE makes it
   a schema constraint.)

**Prompts are nested dicts, not strings** — `{"Introduction", "Task description", "Memory",
"Previous solution", "Instructions": {...}}` with `Instructions` merged from reusable blocks, so
implementation guidelines cannot drift between draft/debug/improve modes.

**The Journal is the experiment ledger.** `Node`: `code, plan, step, id, ctime, parent, children,
_term_out, exec_time, exc_type/info/stack, analysis, metric, is_buggy`; derived `stage_name ∈
{draft, debug, improve}`, `is_leaf`, `debug_depth`. `generate_summary()` emits per good node
`Design: <plan>` / `Results: <analysis>` / `Validation Metric: <value>` and deliberately excludes
code — that is the `experiments/exp0XX` rollup, generated rather than hand-maintained.

**From AI-Scientist-v2** ([`SakanaAI/AI-Scientist-v2`](https://github.com/SakanaAI/AI-Scientist-v2)),
two structural ideas worth lifting (**licence: "AI Scientist Source Code License v1.0", RAIL-derived
— adapt the schemas, paraphrase the prose, do not vendor**):

- **A 4-stage campaign state machine** with caps: `initial_implementation` (20 iters) →
  `baseline_tuning` (12; **"DO NOT change the model architecture from the previous stage"**) →
  `creative_research` (12) → `ablation_studies` (18; **same datasets as stage 3**). The stage-2
  and stage-4 constraints hold the confound surface fixed — exactly the discipline
  `RESEARCH-METHOD` §9's "closed under one objective is not closed under a new one" needs.
- **A fixed on-disk artifact contract.** Every run must `np.save` a dict of one exact shape to a
  fixed filename, so cross-experiment aggregation is possible at all. Our analogue: a fixed
  `results.json` schema across all `experiments/NNN-*/`.
- **Multi-seed evaluation as a first-class node type** (3 seeds), with seed nodes flagged out of
  the best-node competition.

## 4.2 ★ Curie + RD-Agent — rigor as schema validation, not prompting

### Curie — [`Just-Curieous/Curie`](https://github.com/Just-Curieous/Curie), **Apache-2.0**, [arXiv:2502.16069](https://arxiv.org/abs/2502.16069)

**The experimental-plan schema (`curie/formatter.py`) — adopt this wholesale** as the
`PREDICTION.md` front-matter:

| field | type |
|---|---|
| `hypothesis` | str |
| `constant_vars` | List[str] — held fixed |
| `independent_vars` | List[str] — intentionally changed |
| `dependent_vars` | List[str] — measured |
| `controlled_experiment_setup_description` | str |
| `control_group` | List[Dict] — **concrete values, not prose** |
| `experimental_group` | List[Dict] — one per condition |
| `priority` | int > 0 |

The load-bearing piece is a hard `@model_validator(mode="after")` called `control_group_has_vals`
that raises when the list is empty. **The agent physically cannot register a plan with no control
group — the tool call fails.** That is the difference between a rule and a gate.

**Two validators:**
- *Static* (`prompts/llm-verifier.txt`): recursive checks over every script in the chain —
  workflow addresses the question; `cat` every dependency recursively; **any hardcoded or mock
  data**; **any placeholder values**; **every declared variable including constants is actually
  used**, searched recursively; results file must exist. Records `{is_correct, log_message}` and
  is **forbidden from modifying files**.
- *Dynamic* (`nodes/exec_validator.py`): re-runs the control script in a clean env up to 3 times
  and hands the analyzer **all N runs concatenated**. Distinguishes non-zero exit, timeout (not
  an error), and **succeeded-but-no-results-file → `is_correct=False`**.

Also: **tiered write access** — a technician may append results only to its own partition and
cannot touch unrelated plan sections. Applied here: an experiment may not retroactively edit
another experiment's record *or its own hypothesis*.

### RD-Agent — [`microsoft/RD-Agent`](https://github.com/microsoft/RD-Agent), **MIT**

The strictest "did this actually support the hypothesis" gate found anywhere
(`scenarios/data_science/dev/prompts.yaml`, key `exp_feedback`). Six steps; **1 and 2 can kill a
numerically-better result**:

1. Submission-format check. Fail ⇒ do not replace best; reasoning must begin `[Submission format error]`.
2. **Alignment with the task** — is there a mismatch between validation and test performance?
   Sub-checks: exact match between the validation metric and the official metric; consistent
   prediction methodology; no fold-specific strategy applied inconsistently; **claims of unreliable
   validation must cite concrete evidence from the code, not assumptions**. **Two-tier severity:**
   *severe* (likely to invert the val→test trend — leakage, distribution shift) ⇒ not aligned;
   *mild* (scaler fit on full train but applied consistently) ⇒ aligned, but the bias must be
   written into the reasoning. **This is exactly the harness-vs-board problem in
   `RESEARCH-METHOD` §1, formalised.**
3. Compare results, ensemble-vs-ensemble only.
4. **Tie-break on code quality when scores are similar**, in priority order: less overfitting /
   no leakage → best practice → interpretability → efficiency.
6. Acceptable only if 1, 2 and 4 all pass.

Two more rules to steal:
- `ExperimentFeedback.from_exception(e)` forces `decision=False`. **A crashed run is never evidence.**
- **The novelty rule:** score a hypothesis 1/10 if it repeats or is substantially similar to a
  previously attempted one — *successful or failed* — unless the prior attempt died on a trivial
  implementation bug and this is the correct implementation of the same core idea. That single
  rule kills the "re-run the same thing with cosmetic changes and call it new" loop, and it is
  the enforcement `RESEARCH-METHOD` §9's closed-list currently lacks.

## 4.3 ★ The co-scientist supervisor — a campaign controller with numeric triggers

Google's *Towards an AI co-scientist*, [arXiv:2502.18864](https://arxiv.org/abs/2502.18864)
(paper text is Google's copyright — **paraphrase**). Agents: Supervisor, Generation, Reflection,
Ranking (Elo, new hypotheses start at 1200), Proximity, Evolution, **Meta-review**.

The **meta-review loop** is the most adaptable idea: it synthesises all reviews and tournament
debates into recurring critique patterns, and that text is **appended to the Reflection and
Generation prompts on the next iteration**. No fine-tuning — the system learns by accumulating
its own review history in context. Applied here: every ~5 experiments, produce a "recurring
criticisms of our own experiments" block and prepend it to the next design prompt.

Reflection's **deep verification review** is the other one: decompose the hypothesis into
constituent assumptions, break each into sub-assumptions, **decontextualise** them, and evaluate
each independently — so the reviewer cannot rescue a bad sub-assumption using the hypothesis's own
framing. Failures are *fatal* if the assumption is fundamental, *non-fatal* if peripheral.

**The open reimplementation has the real prompt text:**
[`conradry/open-coscientist-agents`](https://github.com/conradry/open-coscientist-agents), **MIT**,
LangGraph, prompts as `.md` files in `coscientist/prompts/`.

**`supervisor_decision.md` is the single most directly liftable artifact in this report** — a
complete campaign state machine with numeric thresholds. Six actions (`generate_new_hypotheses`,
`evolve_hypotheses`, `expand_literature_review`, `run_tournament`, `run_meta_review`, `finish`),
routed on observed state (hypothesis counts, Elo trajectory, **average pairwise cosine similarity
of hypotheses**, **number of distinct clusters**):

- **generate** when total < 8–10, or avg cosine similarity > 0.85 (too similar), or median Elo < 1300
- **evolve** when ≥4 hypotheses exceed Elo 1300 **and** cosine similarity < 0.85 (diversity exists,
  so refining won't over-optimise)
- **meta-review** at ≥4 newly ranked since last, always at ≥10, and **on a performance plateau**
- **finish** only when ≥3 exceed Elo 1400, max/median Elo has plateaued over ≥3 meta-reviews,
  **and the most recent action was a meta-review**
- explicit exploration schedule: <12 hypotheses = explore, 12–25 = balance, 25+ = exploit

**Adapted for this campaign: swap Elo for the board/harness score and hypothesis cosine-similarity
for diversity of experiment configs, and you have a principled "what should exp097 be?"
controller with numeric triggers instead of vibes.** Note the plateau trigger would have fired
here — nine experiments in the `reach` member while `mse` sat untouched.

Other prompts in that repo worth reading: `deep_verification.md` (*"Do not be unnecessarily
charitable"*, and **"Do not pass a final judgement"** — judgement is deferred to the tournament so
the verifier cannot shortcut to a verdict); `tournament.md` (first instruction: **disregard any
numerical scores in the reviews, they may not be comparable across reviews**);
`meta_review_tournament.md` (identify **bias patterns in the review process itself**).

## 4.4 Claude-native artifacts, drop-in

**Negative finding:** `anthropics/skills` contains **nothing** on research methodology.
The material is third-party.

- **[`theam/limina`](https://github.com/theam/limina) → `skills/experiment-rigor/SKILL.md`**,
  **Apache-2.0**. Best single Claude skill for this. Two things to add to our `PREDICTION.md`
  template today:
  - a **`Shortcut Risks`** field: *what could make this result look good without the underlying
    capability improving?*
  - the **negative-result taxonomy**. Before writing "X doesn't work," classify into exactly one
    of: `valid negative` (method-valid test failed the threshold) / `invalid test` (setup didn't
    expose the capability) / `implementation failure` / `insufficient signal` / `trade-off
    failure`. **Only `valid negative` licenses REJECTED; everything else is INCONCLUSIVE plus a
    redesign.** This is `RESEARCH-METHOD` §3 with teeth, and it is exactly the discipline the
    φ=0 bug and the "amplitude closed" ruling both needed.
  - a 6-question **method-validity gate** answered before the hypothesis is finalised, including
    *"could failure be explained by implementation, preprocessing, evaluation, or dataset mismatch
    rather than the method itself?"*
- **[`MarcoMorandin/thesis-with-context`](https://github.com/MarcoMorandin/thesis-with-context)** —
  a *working* Claude Code ML-research campaign. **No licence file → design pattern only, do not
  copy text.** Structure worth copying: a `/register-experiment <id> "<hypothesis>"` slash command
  that writes a registry row, creates a config containing **only the delta from base**, and cuts a
  branch; `knowledge/ablations.md` as the registry (`ID | Hypothesis | Config | Branch | Status |
  Result`, every run registered **before** it launches); an `experiment-reviewer` subagent that
  emits a ✓/✗/? checklist and **blocks the run on any ✗**; a `result-aggregator` that audits
  result JSONs before any table is trusted and is **instructed to report only, fix nothing**; and
  **`.claude/hooks/guard-data-readonly.py`, a `PreToolUse` hook that `sys.exit(2)` on any write
  under the dataset path — the data of record is read-only by hook, not by promise.**
  Their registry includes exactly the negative controls we need analogues of: *"Shuffled-frames
  control (vision actually read)"*, *"Mismatched-frames control"*, *"Vision-only upper bound"*.
- **[`yogsoth-ai/de-anthropocentric-research-engine`](https://github.com/yogsoth-ai/de-anthropocentric-research-engine)**,
  **Apache-2.0**, 900+ skills. Three are directly applicable:
  - **`artifact-detection`** — the shortcut-learning probe. Stage 1 is a **partial-input
    baseline**: run with the key input removed. Threshold: **a gap of <10 points between partial
    and full input indicates severe artifacts.** Our analogue: a **perturbation-blind** baseline
    (context mean, identical rows) and measure how much of our score it already captures. We can
    read this off the board today — em2.bio's mean_realized scores −0.04975 against our 0.15436,
    so the gap is real; but the same test per-member is worth running.
  - **`circular-validation-audit`** — a pre-build gate on any validator, with a **non-circularity
    matrix**: rows = load-bearing claims, columns = the validator's ground-truth construction
    choices, cells = `CIRCULAR` / `INDEPENDENT` / `PARTIAL`. Green-light only if every load-bearing
    claim has at least one INDEPENDENT test. Directly applicable to the two local harnesses.
  - **`ablation-design`** — budget rule: N components ⇒ N+1 runs minimum, with mandatory
    full-system and minimal-system anchors.
- **[`ARA-Labs/Agent-Native-Research-Artifact`](https://github.com/ARA-Labs/Agent-Native-Research-Artifact)
  → `skills/rigor-reviewer/SKILL.md`**, **MIT**. A post-hoc auditor scoring write-ups 1–5 on six
  dimensions. Best items: the **name-deletion test** (strike your component names from a claim; if
  no transferable relationship survives it is a league table, not an insight); **dead-end failure
  modes must be specific** ("divergence after 1000 steps due to gradient explosion", not "didn't
  work") and **a tree with zero dead-ends is suspicious**; findings carry an `evidence_span` that
  **must be an exact verbatim quote from the artifact**.

## 4.5 Structural patterns that beat any prompt

- **METR RE-Bench** (MIT): `assets/torch_rule_enforcer.py` monkeypatches ~20 torch ops to raise
  under an env var; `heldout_setup.py` constructs the held-out set **outside the agent's
  directory** — the split is a filesystem boundary. Their measured result on deterrence is the
  point: a clarification that copying is not allowed produced **no** improvement, while an explicit
  disqualification clause dropped reward hacking from 10/10 to 3/10. **Treat prompt-level
  deterrents as brittle and back them structurally.**
- **`openai/mle-bench`** (MIT): `validate_submission()` validates by **actually running the real
  grader** against private answers and catching the exception, rather than reimplementing a format
  check — so the validator can never drift from the grader. Our analogue: **validate against
  `ArcInstitute/cell-eval2` itself** (§1.12), never against a reimplementation.
- **Weco CLI**: the evaluation harness lives **outside the agent's edit surface, by construction**
  (`weco run --source <editable file> --eval-command "<fixed harness>"`).
- **`Future-House/robin`** (Apache-2.0): a `Step` has `parallel: N`, so the same analysis runs N
  independent times and a consensus step meta-analyses them, plotting **agreement rate** (all / >50%
  / <50%). *"How many of N independent analyses agreed?"* is a cheap, strong self-deception check —
  and it is the generalisation of `RESEARCH-METHOD` §1's "run every arm on BOTH harnesses".
- **Anti-pattern, flagged deliberately:** AgentLaboratory's `mlesolver.py::get_score` uses an LLM
  "professor agent" to score 0–1 how well the code followed the plan. That measures *plan
  adherence*, not held-out performance. Do not copy it.

## 4.6 The consolidated checklist, deduplicated

`[H]` = make it a **hard gate** (schema validation, hook, or code), not a prompt.

**Pre-registration**
1. `[H]` Plan declares `hypothesis`, `independent_vars`, `dependent_vars`, **and `constant_vars`**.
2. `[H]` Non-empty `control_group` with concrete values.
3. What single decision does this resolve? If it cannot change a decision, do not run it.
4. Confirm / reject / inconclusive thresholds written **before** seeing results. *(already §2)*
5. `[H]` Novelty check against the closed-list: is this substantially a previously attempted
   hypothesis? If so it is not new — unless the prior attempt died on a trivial bug.
6. **Shortcut Risks**: what could make this look good without the capability improving?
7. `[H]` Exactly one variable changed; config expressed as a **delta from base**, so the diff *is*
   the experiment.

**Setup validity, before trusting any number**
8. Does the offline metric **exactly** match the scored metric — formula, normalisation,
   aggregation, gene axis? *(the `n_real` problem, already §1)*
9. `[H]` Did any transform get **fit** on data including the evaluation rows? Classify *severe*
   (can invert the trend) vs *mild*, and record which.
10. `[H]` Is the held-out set a **filesystem boundary**, not an intention?
11. `[H]` Is there a runtime guard that raises on a forbidden data touch?
12. Are all declared variables — **including constants** — actually used, checked recursively?
13. Any hardcoded values, mock data, or placeholders in the path that produced this number?
14. `[H]` No metric printed ⇒ buggy ⇒ worst score.
15. `[H]` Circularity: does the validator's ground-truth construction *assume* the claim it tests?

**Result validity, before calling it a gain**
16. `[H]` Re-run end-to-end from clean state ≥2–3 times; report all runs side by side. *(already §7)*
17. Is the delta larger than seed-to-seed spread? **A single-run quantitative claim is flagged by
    default.**
18. Report the **delta against the named baseline**, not an absolute score.
19. Is the improvement in the *pre-registered primary* metric, or did the metric get switched
    after seeing results?
20. **Partial-input control:** how much does a mechanism-blind baseline already capture?
21. **Negative controls:** shuffle / mismatch / zero the input the model is supposed to be using,
    and confirm performance degrades.
22. **Name-deletion test:** strike your component names from the claim — does a transferable
    relationship survive?
23. Tie-break for similar scores: less leakage → fewer overfitting risks → better practice →
    interpretability → efficiency. Not "the one I built last."

**Negative results and bookkeeping**
24. `[H]` Classify before writing "X doesn't work": `valid negative` / `invalid test` /
    `implementation failure` / `insufficient signal` / `trade-off failure`. **Only the first
    licenses REJECTED.**
25. A sub-component that did not run leaves its part of the hypothesis **unverified**, not refuted.
26. `[H]` Any run that raised is `decision=False`. A crash is never evidence.
27. Were previous "fixes" actually effective, or only partially — is the issue still open?
28. Dead ends documented with a **specific** failure mode. A log with zero dead-ends is a red flag
    on the log.
29. `[H]` Raw per-run metrics in machine-readable files with code revision, data version, seeds —
    never only in prose.
30. `[H]` Tiered write access: an experiment may append to its own record and may **not**
    retroactively edit another's, or its own hypothesis.
31. Numbers in the write-up are **cited from** the results file, never retyped.
32. Every ~5 experiments, **meta-review** your own past reviews for recurring criticisms and
    prepend that block to the next design prompt.
33. Judge artifacts and code, not the agent's own narrative.

## 4.7 Licences

| source | licence | guidance |
|---|---|---|
| `WecoAI/aideml` | MIT | copy freely |
| `SakanaAI/AI-Scientist(-v2)` | AI Scientist Source Code License v1.0 (RAIL-derived) | **adapt schemas, paraphrase prose**; distribution requires shipping the full licence and propagating §3.2 restrictions |
| `Just-Curieous/Curie` | Apache-2.0 | reuse with notice; schema field names safe |
| `microsoft/RD-Agent` | MIT | copy; rewrite metric-specific text |
| `conradry/open-coscientist-agents` | MIT | copy freely |
| Google co-scientist paper | Google copyright | **paraphrase only** |
| `theam/limina` | Apache-2.0 | reuse with notice |
| `ARA-Labs/Agent-Native-Research-Artifact` | MIT | copy |
| `yogsoth-ai/de-anthropocentric-research-engine` | Apache-2.0 | reuse with notice |
| `openai/mle-bench`, `METR/RE-Bench` | MIT | copy |
| `MarcoMorandin/thesis-with-context` | **none** | design pattern only |

---
# PART 5 — what I would do, in order

Ranked by (expected points) × (confidence) ÷ (cost). Items 1–3 are all the same finding.

1. **Fix the leaderboard scraper to capture `description`.** Zero risk, one field, and 133 teams
   are publishing their methods into it. `https://virtualcellchallenge.org/api/leaderboard`.
   (§1.11)
2. **Raise emission dispersion to the `#247` cap. This is the lever, and it is free.**
   Our `Ĉ_p` is 0.67× the reference `C_p` and the cap fires on 0 of 900 rows — we forfeit 0.0798
   of claim, worth **0.349 of raw mse**. 096's measured recipe: top-8000 genes at (spread 0.15,
   φ 16), tail at (spread 1.0, φ 4), reaching the cap at only 1.16× reference dispersion.
   *(This replaces my original item 2, which recommended a 2-D sweep on the premise that we were
   over-dispersed. We are under-dispersed; see the Part 0 retraction.)*
3. **Then set amplitude to the cliff, `s_c = √(claim/B)` — a 6% cut, not 60%.** `VCC_ALPHA`
   1.0796 → 1.012. There is **no interior optimum**: `mse(s)` falls until `ρ` saturates and rises
   after, because the predictions' common mode is only 0.84% of their energy. Raising the claim in
   step 2 raises the cliff too, so the two steps compound. Predicted: raw mse **0.8785**,
   `score_mse` **+0.115**, **+0.0192 overall** → ≈ 0.1736, rank 3–4 *on today's board* (§"The rank
   projection has a second, larger uncertainty"). Verify on both harnesses that the added dispersion
   lands only on genes that were never callable (§1.9). **Note the ceiling**: this route tops out at
   mse raw ≈ 0.870 against the leader's 0.8270, because their aggregate cosine is at least 1.28×
   ours (Part 0 retraction). Calibration is the cheap points, not the win.
3b. **Then per-gene shrinkage — now measured, and it is the best buildable lever found.**
   096 tested Miller's `μ_ID = α ⊙ μ_TD + (1 − α) ⊙ μ_all` with `α` keyed two ways, re-optimising
   the scale onto the cliff each time so `A`, `B` and `X` move together (which answers my
   objection that it lowers the cliff). **A donor-agreement key (gwps × hek) beats a self-|dh|
   key on both harnesses at every setting**; transferable estimate 1.17× on `√W`, optimistic
   1.38×:

   | | raw mse | `score_mse` | overall |
   |---|---:|---:|---:|
   | intervention 2+3 only | 0.8785 | 0.1151 | +0.0192 |
   | + Miller at 1.172× | 0.8536 | 0.1414 | **+0.0236** |
   | + Miller at 1.381× | 0.8232 | 0.1733 | **+0.0289** |

   The 1.381 end lands **below the leader's 0.8270**. The mechanism note is the generalisable
   part: **it fits no slope.** It interpolates between two profiles the model already produces,
   weighted by an *external* confidence score. A per-gene regression slope is a **cell-line
   property** and transfers at 1.02×; a confidence score is a **model property** and transfers at
   1.17×. That is why this works where the per-gene slope did not. Caveats: donor coverage was
   81/85 and 132/132 on the harnesses and must be checked against the challenge's 300 targets;
   and this is mse-only accounting — shrinking toward the panel mean moves `pds`, which reads
   exactly that geometry. (§2.5)
3c. **And treat raising the aggregate cosine as the only route to first, not as the expensive
   option.** Corrected inversion (Part 0) puts the leader's `cos` at ≥ 0.0947 against our 0.0740,
   and our ceiling with claim at the `#247` cap is raw ≈ 0.870 against their 0.8270. Every
   calibration lever above is bounded by that. Combined with the board-drift finding — the #1
   threshold is the fastest-climbing rank — **the cheap interventions buy rank 3–5 and hold it for
   about a day; only a directional gain competes for first.** The literature's concrete candidates
   are in §2.3 and §2.5, and 096's measurement supplies a **selection criterion** that reorders
   them: prefer whichever is easiest to express as a **per-(p,g) external confidence** rather than
   as a per-gene coefficient, since external confidence is the only thing that has transferred
   across surfaces in this entire campaign.

   1. **KG k-NN / functional priors over gene identity** (RF+GO beats every foundation model;
      [arXiv:2606.08816](https://arxiv.org/abs/2606.08816) reports plain k-NN over a biological
      knowledge graph beating almost all methods o.o.d.). **Best fit by the criterion**:
      neighbour consensus is *natively* a per-(p,g) external confidence — the agreement among a
      perturbation's KG neighbours at gene `g` — and is structurally the same object as donor
      agreement with KG neighbours substituting for donor screens. Cheapest test of the family.
   2. **COMPASS's β/residual split.** Medium fit — β is per-perturbation and the residual per-gene,
      so it is a decomposition rather than a confidence score. But it independently *measures* the
      property that matters: shared magnitude transfers across lines at **r = 0.66**, gene-specific
      residual at **r = 0.22** — the same asymmetry as our per-gene slope (1.02×) versus per-gene
      confidence (1.17×).
   3. **Ahlmann-Eltze's perturbation-pretrained latent basis.** Strongest standalone method in its
      benchmark, but **weakest fit by this criterion** — it yields a prediction, not a confidence,
      so it would need a separate residual model to produce a weight.

4. **Revisit the call set — but the CPM-keyed version of this is retracted.** VIVAI calls 3,224
   and beats us; a10 calls ~4,757. `sig_jaccard` is symmetric and penalises padding. **But**
   `fid`'s coverage term saturates at `n_pred = n_real ≈ 4,757`, so cutting *below* that costs
   `fid` proportionally — the two members pull in opposite directions and 4,757 is close to the
   `fid`-optimal point already. And my "allocate by control CPM" inference was measured false on
   our submission (§1.11): our called genes sit at *lower* control CPM than our uncalled ones.
   The live question is precision, not size. (§1.9, §1.11, §3.4)
5. **Stop work on `fid` and `jac`.** Board-wide maximum on `fid` raw is 0.5222 against a chance
   level of 0.5000 over 385 teams; `jac`'s scaled max is 0.0166. Combined realistic upside is
   ~+0.005 overall. (§1.6)
6. **Do not chase Gram whitening.** AgingLab_BI says in their own description that it *"pays cos to
   buy pds"* — it buys the one member we already lead the board on and costs the one we score zero
   on. Their `mse_raw` is 2.7357. (§1.11)
7. **Track `ArcInstitute/cell-eval2` and validate against it directly** rather than against a
   reimplementation — MLE-bench's rule that the validator must *be* the grader. Read
   `src/cell_eval2/metrics/delta.py`. Also read `AUTHORS` #279: `reach`'s `k*` is *"decided at the
   head of a ranking the submission controls, which makes it exploitable"*, and Arc has left that
   half open. (§1.12)
8. **Adopt three cheap methodology additions** to `docs/RESEARCH-METHOD.md`: Curie's
   `control_group` schema validator as a hard gate; limina's `Shortcut Risks` field and 5-way
   negative-result taxonomy (only `valid negative` licenses "closed"); and the co-scientist
   meta-review every ~5 experiments. The plateau trigger in `supervisor_decision.md` would have
   fired on this campaign — nine experiments in `reach` while `mse` sat at a hard zero. (§4)

## The rank projection has a second, larger uncertainty: the board moves

`VERIFIED` from the ten snapshots in `data/leaderboard/`. Experiment 096 stress-tested its
intervention table against ±10% anchor-span error and found the **decision** invariant (the mse
gain is 12–19× the sum of the reach/nmae/fid costs throughout) but the **projected rank** not:
0.17347 → rank 3, 0.17156 → rank 3, 0.16999 → rank 5. That is my §1.9 within-member /
cross-member distinction landing on their headline, and it is right.

**But anchor error is the smaller of the two uncertainties.** The board climbs:

| interval | days | #1/day | #3/day | #5/day | #16/day |
|---|---:|---:|---:|---:|---:|
| 08-23 → 08-25 | 2.46 | +0.00780 | +0.00541 | +0.00841 | +0.01261 |
| 08-25 → 08-27 | 1.87 | +0.01283 | +0.00649 | +0.00062 | +0.01058 |
| 08-27 → 08-28 | 0.88 | **+0.02421** | +0.00441 | +0.01544 | +0.01404 |

The #3 threshold is decelerating (+0.00650 → +0.00441/day). **Even at the decelerated rate, one
day of board drift is 1.3× 096's entire ±10% span-induced spread of 0.00348.** Propagating a fixed
0.17156 forward against a uniformly drifting board:

| submitted after | projected rank |
|---|---|
| 0 days | 3 |
| 1 day | 6 |
| 2 days | 7 |
| 3 days | 13 |

**Caveat, stated because the number is dramatic** (`INFERRED`): this is linear extrapolation of a
noisy rate over a short window, and the #1 figure in particular is driven by GeroAI v12's single
jump inside a 0.88-day interval rather than a smooth trend. Treat the direction as solid and the
3-day figure as an upper bound on decay, not a forecast.

**Two things follow.** First, **time-to-submission dominates anchor error in the rank
projection** — an extra day of polishing costs more expected rank than a 10% anchor error does,
which is an argument for 096's shipping order (cheapest intervention first, immediately) on
grounds independent of build risk. Second, **the #1 threshold is climbing fastest of any rank**;
the leader is pulling away rather than consolidating, so intervention 1 alone (+0.0192) closes the
gap to *today's* leader but not necessarily to the leader we will face.

## What would falsify the central claim

The corrected claim is that our mse deficit is **forfeited noise credit plus amplitude past the
`ρ` cliff**, not directional error. It predicts:

- `mse(s)` **monotone decreasing** below the cliff and increasing above it, with the argmin at
  `s_c = √(claim/B)` and `ρ` crossing 1.000 there. **Falsified by** an interior optimum strictly
  below the cliff, or by `ρ` already at 1.000 across the whole grid.
- Raising `Ĉ_p` toward the reference `C_p` **lowers** raw mse roughly as `√claim`. **Falsified by**
  raw mse failing to improve once `Ĉ/C_p` passes 1.0, or by the `#247` cap firing on a material
  fraction of the 900 rows before then.
- **Resolved, no longer open.** The dispersion churn's cost was bounded by the directional
  precision of the ~1,060 new calls; measured at **0.4946 (h1-test) / 0.5131 (orion)** — chance,
  bracketed. Dropped calls likewise at chance. **Board prediction: `score_fid` falls −0.0016 to
  −0.0027, and `n_pred` must stay above `n_real` ≈ 4,757.** Intervention 1's total DE cost is
  **−0.0015** (amplitude −0.0010, churn −0.0005) against **+0.0192**. **Falsified by** `score_fid`
  falling more than ~0.006, or `n_pred` dropping below `n_real` (which would de-saturate the
  coverage term and cost proportionally). (§1.9)
- `pds_cosine` is near-flat across the amplitude move. **Falsified by** `pds` moving materially,
  which would mean the change is not a scalar.

096's out-of-sample check already passed once: a8/Wildflower was not in the two-point fit and was
predicted at 1.3671 against an observed 1.3627 (+0.32%). Each remaining item is measurable
locally, without a submission slot.

## Proposed addendum to `docs/RESEARCH-METHOD.md` — NOT APPLIED, out of scope

My brief is read-only on the repo except this file, so I have **not** edited
`docs/RESEARCH-METHOD.md`. The 096 owner declined the same edit on the same grounds. Recording it
here for whoever owns the standing rules:

> **§2 addendum — agreement between two agents is not a check.** Two agents reading the same spec
> and reaching the same conclusion have produced **one reading held twice**, not independent
> confirmation. In experiment 096 the *first* agreement reached was the wrong one: two independent
> analyses both omitted `Kc` and `ρ·claim` and agreed on a shrink of ~0.4× where the answer is
> ~0.94×. What broke it was running `cell_eval2`'s own functions, not further argument.
> **Convergence between readers should LOWER confidence relative to a single agent who has
> executed the code.**

A second candidate, from §1's own founding example being mis-cited in this document (§1.9):

> **§1 addendum — "same policy" is not "same artifact".** `realize.py` re-realises an emission
> policy against each harness's own controls; it does not score the shipped `.vcc`. Any claim of
> the form "we ran the same thing on both surfaces" must name the artifact, not the policy. And a
> prediction artifact is **baseline-bound** — scoring one against another cell line's reference
> measures the control gap, which was measured at **548×** the per-row perturbation signal.

## Open questions I could not resolve

- **Who Ibrahim Mansour is** (#5, `CellSim`, one of only two teams beating the mse baseline
  comfortably, 15 submissions, no description). No GitHub, Scholar, LinkedIn, or preprint trace.
- **What We527 actually ran.** Model name is literally `STATE baseline` and it is rank 3 with a
  positive mse, while every other self-declared STATE run on the board is at the bottom. No
  identity, no description. `INFERRED`: probably a misleading label.
- **BioMap's and XLearning Lab's 2025 methods** beyond press-release level — no repo, no preprint,
  no non-English writeup surfaced.
- ~~Whether the `ρ` "factor of about 3" holds for our submission.~~ **Resolved by measurement:
  it does not.** Our `B/claim` is 1.66–2.37 depending on emitter, putting the cliff at
  `VCC_ALPHA` 0.70–0.84 rather than 0.577.
- ~~Whether our directional accuracy is concentrated in the perturbations that discriminate well.~~
  **Resolved by measurement, and the effect is larger than I framed it.** 096 measured the
  per-perturbation cosine distribution:

  | | h1-test | orion |
  |---|---:|---:|
  | aggregate cosine | 0.0519 | 0.0325 |
  | per-perturbation cosine, mean ± sd | +0.035 ± 0.043 | +0.017 ± 0.030 |
  | p5 / p50 / p95 | −0.005 / +0.027 / +0.111 | −0.019 / +0.010 / +0.076 |
  | perturbations with cos ≤ 0 | 11.8% | 28.8% |
  | **top quartile's share of total X** | **79.7%** | **85.0%** |

  **A quarter of the perturbations carry four-fifths of all our directional signal, and between a
  ninth and a quarter point the wrong way.** That is why a rank statistic (`pds_cosine`, where we
  lead the board) sits at the top while the bulk cosine sits near the bottom — they are measuring
  different tails of the same distribution.

  **Only one confidence score sees it**, by Spearman against the true per-perturbation cosine:

  | score (all buildable, no reference) | ρ h1 / orion | `√W` multiplier |
  |---|---:|---:|
  | `‖δ̂_p‖` (self) | −0.001 / +0.253 | 0.595× / 0.983× |
  | uniqueness `min_q(1−cos)` (self) | −0.337 / −0.102 | 0.698× / 0.915× |
  | **donor agreement (gwps × hek)** | **+0.209 / +0.242** | **1.144× / 1.035×** |
  | oracle `cos_p` (leak, ceiling) | 1.000 | 1.838× / 1.676× |

  **Both self-referential scores are useless or actively harmful and disagree in sign across
  surfaces; the external one is positive on both.** This is the third instance of one cut, and it
  now generalises past the slope-vs-confidence rule: **external signals transfer, self-referential
  ones do not.** Per-gene slope (cell-line property) 1.02×; per-gene donor confidence (external)
  1.17×; per-perturbation donor confidence (external) 1.035×. It also *predicts* which methods will
  fail — anything deriving confidence from the model's own output.

- **A board-side test of the donor-agreement hypothesis, reported as inconclusive.** I classified
  the 133 public descriptions into external-agreement methods (IVW, multi-source consensus,
  reliability weighting) vs self-referential gating (top-K by own |effect|, magnitude gates):

  | group | n | median `score_avg` | median rank |
  |---|---:|---:|---:|
  | external agreement | 23 | 0.05209 | 155 |
  | self-referential gating | 7 | 0.02236 | 184 |
  | all with descriptions | 133 | −0.00489 | 207 |

  Directionally consistent, **but I would not lean on it**: no control for team quality (teams
  writing careful methodological descriptions are probably better across the board), n = 7 in one
  arm, and crude regex classification. **The strongest caveat is structural: not one
  external-agreement team is in the top 50** — the best is #51. And decisively, **nobody on the
  board is using donor agreement the way 096 is.** Those teams use it to *build the delta*
  (blending source screens); 096 uses it as a *confidence weight* on our own predictions. Same
  signal, different object. **The board cannot test this hypothesis — the harness measurements are
  the only evidence, and they should not expect board-side corroboration.**

- ~~Does a second independent external key add to the donor key, or duplicate it?~~ **Resolved:
  it adds, and the answer is that EXTERNAL agreement is the operative property — not donor screens
  specifically.** 096 ran it using a STRING embedding asset (`data/external/slim/
  gene_string_embeddings.v0.3.h5`, 21,688 symbols × 64 features) that **experiment 087 had
  preregistered for exactly this ablation and never run**; coverage 85/85 and 131/132 targets, so
  the second key cost one script.

  | per-perturbation score | ρ h1 / orion | `√W` mult | corr. with donor key |
  |---|---:|---:|---:|
  | donor agreement (screens) | +0.209 / +0.242 | 1.144× / 1.035× | — |
  | **`kg_density`** (STRING neighbourhood) | **+0.190 / +0.187** | 0.979× / 1.029× | −0.053 / +0.304 |
  | `kg_consensus` (STRING neighbour agreement) | +0.190 / **−0.058** | 0.982× / 0.972× | +0.161 / −0.140 |
  | donor + `kg_density`, averaged | **+0.286** / +0.224 | 1.104× / 1.046× | — |

  **`kg_density` contains no perturbation data at all** and still tracks our true per-perturbation
  directional accuracy at ρ ≈ +0.19 on *both* surfaces, nearly as well as the donor key and largely
  independent of it. Two unrelated external priors see the same thing: **we are more accurate on
  well-characterised genes.**

  **⚠ My own recommendation was half wrong.** I argued for KG **neighbour consensus** on the
  grounds that it is natively a per-(p,g) confidence. That is the score that **flipped sign between
  surfaces** (+0.190 / −0.058) and should not be used. The one that worked is the cheaper and less
  clever **neighbourhood density**. Right family, wrong member — and the failure mode is the one
  this document keeps recording: the more elaborate construction was the one that did not transfer.

  096's own calibration is worth preserving: **confident about ρ on both surfaces, only weakly
  confident about the `√W` gain** (1.035× → 1.046× on the rule-1 figure — real but marginal). The
  limiting factor is the **functional form**, not the signal: a mean of two min-max rescales is the
  crudest possible combiner against a 1.68× oracle ceiling.

  **⚠ RESOLVED — the confound was real, and it splits by surface.** 096 ran the partial Spearman
  against `cos_p` controlling for true `‖δ_p‖`:

  | | h1-test | orion |
  |---|---:|---:|
  | confound strength (true `‖δ_p‖` vs `cos_p`) | **+0.739** | **+0.104** |
  | `kg_density` raw → partial | 0.190 → **−0.118** | 0.187 → +0.178 |
  | **donor agreement** raw → partial | 0.209 → **+0.144** | 0.242 → **+0.229** |
  | `‖δ̂_p‖` raw → partial | −0.001 → +0.002 | 0.253 → +0.242 |

  **The donor key survives on both surfaces; `kg_density` does not.** Where magnitude predicts
  accuracy strongly (h1, 0.739), `kg_density`'s partial goes *negative* — its signal was magnitude.
  Where the confound is weak (orion, 0.104) it survives intact. The two surfaces disagree precisely
  because the confound's strength differs between them.

  **Corrected reading:** external *is* the operative property — `kg_density` carries something on
  the surface where magnitude is not doing the work — but its apparent accuracy signal is at least
  partly and possibly *entirely* a magnitude proxy, and two surfaces cannot exclude "entirely".
  **The donor key is the only score in this experiment that predicts where our predictions are
  trustworthy after effect size is removed.** Build the combiner on it alone; admit `kg_density`
  only after the confound is excluded on a third surface; and run 087's protocol with `‖δ_p‖`
  partialled out from the start — its preregistration does not currently ask for that. A magnitude
  prior remains worth having, but it belongs in the model as an **effect-size prior** (DepMap
  Chronos) rather than in the confidence weight.

  *(One calibration note on the corroboration: 096 read their h1 confound of +0.739 as reproducing
  Raju's ρ = 0.977. Those are **different statistics** — Raju measures within-perturbation cell-to-cell
  coherence against magnitude; this measures our model's agreement with truth against magnitude.
  The direction corroborates; the numbers are not a reproduction, and 0.739 vs 0.977 is a wide gap.)*

  **The original reasoning that motivated the test** (`INFERRED`, mine): dense
  STRING neighbourhoods belong to well-studied hub genes, which are disproportionately **essential**
  — and essential genes have *large* effects. TRADE measures a typical perturbation at ~45 affected
  genes but essential genes at **>500** (§2.2), and Raju measures directional coherence ≈ effect
  magnitude at **ρ = 0.977** on Replogle (§2.5). So `kg_density` may be predicting **true effect
  magnitude** rather than our model's accuracy per se. The two are not equivalent for our purposes:
  a magnitude proxy would partly duplicate any effect-size prior and would not be independent
  information about where *our* predictions are trustworthy. **Cheap test, and both quantities are
  already on the harnesses: partial out true `‖δ_p‖` and see whether `kg_density` retains its
  correlation with the per-perturbation cosine.** If it survives, it is genuinely about
  characterisation; if not, a direct magnitude prior would serve better — **DepMap Chronos
  essentiality** is the obvious candidate, is used by at least one board team (#249 Biomni-TusoAI,
  *"DepMap 24Q2 Chronos essentiality scaling"*), and is free.

- **Whether the leader is at the `#247` cap.** It sets *how far* ahead on direction they are —
  1.28× ours at the cap, 1.56× at our current φ=4 claim, 1.84× at φ=0 — but **no longer whether**,
  since the cap bounds the implied `cos` from below at 1.28× (Part 0). Their `pds_cosine` of
  0.8199 against our board-best 0.8318 is the standing puzzle: a *rank* statistic where we lead,
  against a *bulk* cosine where we are behind by at least 28%. Those measure different things and
  the tension is real, not an error — worth resolving, because it bears on whether our directional
  advantage is concentrated in the perturbations that discriminate well.
- **Arc's official baseline scores** — never published in either year. The 2025 numbers here are
  Ingenix's independent measurements; the 2026 anchors are reconstructed from the board (two
  independent derivations agreed to 1e-4, so high confidence, but they are not published values).

## Sources

Primary: `data/leaderboard/lb_*.json`; `https://virtualcellchallenge.org/api/leaderboard`;
[`ArcInstitute/cell-eval2`](https://github.com/ArcInstitute/cell-eval2) (MIT — source, `AUTHORS`,
`CHANGELOG.md`); `docs/reference/vcc2026-metrics-FULL.md`;
[Arc 2025 wrap-up](https://arcinstitute.org/news/virtual-cell-challenge-2025-wrap-up);
[Arc 2026 announcement](https://arcinstitute.org/news/virtual-cell-challenge-2026).
Literature and repo URLs are inline throughout Parts 2–4.

*Compiled 2026-08-28. Read-only research; nothing was submitted, posted, joined, or registered.*
