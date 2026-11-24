# 102 — donor-estimated COMPASS axis transfer

Status: **completed; preregistered primary arm is a NO-GO**.

Date: 2026-08-29.

## Question

Does separating each source response into an axial loading and an orthogonal
residual improve strict leave-one-donor-out transfer when the target axis is
estimated from source donors only?

This is an independent implementation of the COMPASS equations. No code is
copied from the unlicensed checkout under `/tmp`.

The target context contains no measured perturbations, so this is explicitly a
**donor-estimated target-axis surrogate**, not official COMPASS. The upstream
COMPASS implementation estimates `u_c` from perturbations measured in the target
cell line; VCC 2026 provides only target controls. A positive or negative result
here therefore applies only to this source-only approximation.

## Frozen outer surfaces

| surface | held truth | source donors |
| --- | --- | --- |
| `h1` | H1 test | K562-GWPS + Orion-HEK293T |
| `hct116` | clean Orion-HCT116 | K562-GWPS + Orion-HEK293T |
| `hek293t` | clean Orion-HEK293T | K562-GWPS + Orion-HCT116 |

The held donor is absent from every fitted source mean, scale, axis, beta,
residual, basal weight, and frozen prediction.

## Coordinate and row isolation

For a held surface with gene axis `G` and complete gallery `T`, define

```text
forbidden = VCC_2026_panel_genes union T
eligible  = G minus forbidden
```

The forbidden coordinates are removed before computing any source mean, unit
axis, beta, residual, norm, basal-control similarity, prediction, or PDS.
Zeroing them after a projection is not permitted because a forbidden
coordinate could otherwise alter a scalar beta and leak back into every
eligible coordinate.

Every source geometry on a surface is fitted on the same perturbation rows:
the intersection of source-donor perturbations after removing both the 300
VCC panel targets and that surface's held-gallery targets. The held targets'
source responses remain available only as the per-target inputs being
decomposed after the global source geometry is frozen. This mirrors the
challenge path, where none of the 300 prediction rows is allowed to help fit
the axis it is evaluated under.

For source `c`, on this common fit set `S` and eligible gene coordinates:

```text
m_c       = mean(q in S) delta_cq
s_cq      = delta_cq - m_c
a_c       = median(q in S) ||s_cq||_2
Delta_cp  = delta_cp / a_c
mu_c      = m_c / a_c
u_c       = mu_c / ||mu_c||_2
beta_cp   = <Delta_cp, u_c>
r_cp      = Delta_cp - beta_cp u_c
beta_bar_c = mean(q in S) <Delta_cq, u_c>
```

The common-row median centered norm `a_c` is the only donor-magnitude
equalization. It is shared by every arm.

## Frozen arms

For target `p`, averages use only source donors carrying `p`. Source coverage
is frozen once and must be identical for all arms.

- **A — incumbent centered:** `mean_c(Delta_cp - mu_c)`.
- **B — raw donor average:** `mean_c(Delta_cp)`.
- **C — residual only:** `mean_c(r_cp)`.
- **D — primary centered COMPASS:**
  `mean_c(r_cp) + mean_c(beta_cp - beta_bar_c) * u_target`.
- **E — official-form COMPASS:**
  `mean_c(r_cp) + mean_c(beta_cp) * u_target`.
- **F — basal diagnostic:** arm D with `u_target_basal` in place of
  `u_target`.

The primary target axis is donor-only:

```text
u_target = normalize(mean_c(u_c))
```

For F only, each source axis is weighted by the nonnegative cosine between the
held context's supplied control pseudobulk and that source's control
pseudobulk, using only eligible coordinates present in that source:

```text
w_c = max(0, cosine(control_target, control_c))
u_target_basal = normalize(sum_c normalized(w)_c u_c)
```

If all similarities are nonpositive, F falls back to equal source weights.
The basal weights change only the target axis; source residuals and beta
loadings remain uniformly averaged.

## Beta-permutation null

Twenty deterministic null draws independently permute each source donor's
centered beta values across the held target labels that donor carries.
Residuals, source coverage, the primary target axis, and each donor's beta
marginal distribution remain fixed:

```text
null_jp = mean_c(r_cp) + mean_c(permute_j(beta_cp - beta_bar_c)) * u_target
```

All real and null predictions are frozen to disk before the held perturbation
counts are read.

## Evaluation

The only primary metric is raw, uncentered, complete-gallery PDS:

1. aggregate the held context's real control and perturbation counts;
2. form raw control-referenced log-pseudobulk effects;
3. score all gallery targets against all gallery truths by cosine rank;
4. use exactly the eligible coordinate mask defined before fitting.

No centering, whitening, target subset, DE gate, or partial gallery is allowed
in the primary result.

Target-block uncertainty resamples whole target rows. The macro estimate gives
each of H1, HCT116, and HEK293T equal weight.

## Primary promotion gate

Arm D is a **GO** only if every condition holds:

1. `PDS(D) > PDS(A)` on H1, HCT116, and HEK293T;
2. macro `PDS(D) - PDS(A) >= 0.003`;
3. the deterministic 2,000-draw target-block bootstrap 95% lower bound is
   above zero;
4. the observed macro delta is above the beta-permutation null 95th
   percentile;
5. source coverage is unchanged between A and D on every surface.

Any failed condition is a NO-GO. B, C, E, and F are diagnostics and cannot
rescue a failed D.

## Staged, resumable runtime

Each substantive command is one local process and is capped at **24 GiB**.
The runner also takes an experiment-local nonblocking lock, writes artifacts
atomically, skips complete stages, and resumes missing null draws.

```bash
uv run python scripts/run_with_memory_guard.py --max-rss-gib 24 -- \
  uv run python experiments/102-compass-axis-transfer/run.py prepare \
  --surface h1 --max-rss-gib 24

uv run python scripts/run_with_memory_guard.py --max-rss-gib 24 -- \
  uv run python experiments/102-compass-axis-transfer/run.py freeze \
  --surface h1 --max-rss-gib 24

uv run python scripts/run_with_memory_guard.py --max-rss-gib 24 -- \
  uv run python experiments/102-compass-axis-transfer/run.py score \
  --surface h1 --max-rss-gib 24
```

Repeat those commands serially for `hct116` and `hek293t`, then:

```bash
uv run python experiments/102-compass-axis-transfer/run.py summarize
```

Synthetic tests only:

```bash
uv run python experiments/102-compass-axis-transfer/test_run.py
```

## Hard scope boundary

- Local only; no AWS.
- Never more than one experiment-102 data process at a time.
- No `cell-eval2`.
- No whitening in the preregistered primary comparison.
- No challenge build, submission, or git commit.
- Every generated artifact stays under
  `experiments/102-compass-axis-transfer/`.
