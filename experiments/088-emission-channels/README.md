# 088 — Emission channels: the variance channel, and the mse normalization

**Preregistered 2026-08-23, before any arm was scored.**

## Why

v10 leads the board at `0.149633`. Its gain came mostly from `nmae` (`+0.01782` of
`+0.00524`), not pds. Against the live board we hold the **best `score_pds` of all 150
teams**, so that metric has zero headroom and every remaining point is elsewhere:

| metric | ours | best | headroom → overall |
| --- | ---: | ---: | ---: |
| pds | 0.7338 | 0.7338 (us) | **0** |
| reach | 0.1063 | 0.2051 | **+0.01646** |
| mse | 0 | 0.0826 | +0.01377 |
| nmae | 0.0909 | 0.1460 | +0.00918 |
| fid | −0.0358 | −0.0063 | +0.00493 |
| jac | 0.0026 | 0.0200 | +0.00290 |

The four DE metrics are read off a Wilcoxon on emitted counts, so the **emission model**
is the lever (`STRATEGY.md:227`).

## Track A — the variance channel (`tight_k`)

`experiments/031-somatic/run.py:124` switches spread from **global** (every gene gets the
same Poisson fraction) to **confidence-weighted**: the top-k genes by cross-donor agreement
are emitted at the variance floor, the rest get `spread`. v10 ships `tight_k = 0`, i.e.
global — every gene equally noisy, ordering advantage discarded.

exp027 measured `k=1000, spread=0.15` on h1-fast at `reach +0.0247, 3.4σ` over three seeds,
and put it at **+0.004 overall**. It was never shipped: the response was non-monotone at the
low end on single seeds, and h1-fast is the wrong biology (2025 common-essential genes,
which the 2026 panel excludes).

A third reason is recorded at `run.py:110` — global spread "buys the `expr_mse` jackknife
credit". **We score 0 on mse.** That rationale is void.

### Hypothesis

Confidence-weighted spread raises `reach` and `fid` against the `k=0` baseline, with `pds`,
`mse` and `nmae` unchanged — the quasi emitter preserves each gene's total by construction,
so the pseudobulk cannot move.

### Design

`k ∈ {1000, 3000}` × `spread ∈ {0.15, 0.25}` × `seed ∈ {0,1,2}`, plus the shipped
`k=0, spread=0.25` baseline × 3 seeds. **15 configs.** Production operator throughout:
`g=0.65, q=1.0, a=1.08, cpm=5.0`, union axis, challenge panel excluded.

Three seeds is not optional. exp027 measured the quasi emitter's own sampling sd at `0.0125`
against a `+0.0247` effect; single seeds cannot resolve that, and that is precisely how its
`spread=0.05` non-monotonicity arose.

A 2×2×3 fully-powered design was chosen over a 4×4×3 scan of the same compute. The question
is whether the channel works at all, not where its optimum sits.

## Track C — the mse normalization

The transfer-cosine derivation was retracted in `92a9c3e`. `min_α mse = 1 − cos²θ` requires
`expr_mse_unbiased_capped_norm` to be `||pred−true||²/||true||²`, i.e. predicting nothing
scores exactly `1.0`. Fitting `mse(α) = c0 − 2αA + α²B` to exp077's three alpha points gives
`c0 = 1.3388`. Three points fit a three-parameter quadratic exactly, so that fit had zero
degrees of freedom and validated nothing.

### Design

Four **new** alpha points — `a ∈ {0.20, 0.55, 0.70, 0.90}` at `g=0.65` on the HCT116
comparator, giving seven with exp077's existing `{0.35, 0.40, 1.0796}`. Four degrees of
freedom against a three-parameter model is the first real test of whether the curve is
quadratic and what `c0` actually is.

This decides whether mse's `+0.01377` is reachable, or whether the `0.0039` gap measured by
exp077 means something other than what it appears to.

## Harnesses

**h1-test is primary for Track A.** exp027 nominated HCT116, but HCT116's raw reach is
`0.01502` against the board's `0.16890` and h1-fast's `0.1519` — roughly 10x too little
signal to resolve a `+0.025` effect. h1-test matches the competition on cell counts and
predicted v8's board gain to within `0.00002`.

**Track C runs on HCT116**, because it extends exp077's alpha curve and must be measured on
the same surface as the three points it is being fitted with.

Per `docs/V10-RESULT.md`, proxies are directionally reliable and their magnitudes do not
transfer in either direction — exp065's HCT116 proxy under-predicted v10's board gain by
2.6x, while h1-test over-predicted v9's DE metrics by 4-5x. **These arms rank candidates and
fix signs. No board number is projected from them.**

## Promotion gates, fixed in advance

An arm is promotable to a submission only if:

1. Sign is consistent across **all three seeds** on h1-test.
2. Its six-metric average beats the `k=0, spread=0.25` baseline on h1-test.
3. Its neighbour in `k` or `spread` is also positive. An isolated cell is noise.
4. **pds, mse and nmae are unchanged to within seed noise.** The construction guarantees
   this; movement means an implementation bug, not a finding, and voids the arm.

Failing any gate, the arm is recorded as negative and v10 stands.

## Compute

Capacity is the binding constraint, not budget. The account's EC2 quota is **16 on-demand +
32 spot Standard vCPU**; every other family is 0. Increase requests to 256 were filed
(`L-1216C47A`, `L-34B43A08`) and went to `CASE_OPENED`, i.e. AWS support review rather than
auto-approval. The grids above are sized to 48 vCPU, not to the budget.

- Track A: one spot `r7i.8xlarge`, `--jobs 6`
- Track C: one on-demand `r7i.4xlarge`, `--jobs 3`

## Artifacts

Predictions (~2.3 GB each) and the fitted whitener (~1.3 GB) stay on the worker; only logs,
`results.json`, `summary.json` and checksums are published to
`s3://<your-s3-bucket>/output/exp088-emission-channels/<RUN_ID>/`.
