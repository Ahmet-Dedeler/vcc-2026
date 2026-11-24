# 384 — PromoterAI gene priors: gated upstream, proxies carry no reach (NEGATIVE)

Status: **NEGATIVE** (kill criterion fired on the reach leg).

## Known-answer reproduction

Imported `reach348` walk verbatim (348's directory never written):

| check | here | expected | match |
|---|---|---:|---|
| O1 current panel mean | 0.25434726 | 0.25434726355678233 | exact (0.00e+00) |
| O4 \|design\| panel mean | 0.28074544 | 0.28074543623823345 | exact (0.00e+00) |
| parquet-vs-npz truth | max\|dq\|=0, max\|dlfc\|=0 | 0 | ✓ |
| design universe coverage | 10,777/10,779 | 10,777/10,779 | ✓ |
| 189 τ (all ten pairs) | 0.3417639846 | 0.3417639845890673 | ✓ (1e-12) |
| 189 θ_g axis | 10,725 genes | 10,725 | ✓ |

Proceeded only after this passed.

## What was measured

H1 panel (85 targets, 82 scored; pool = ref q<0.05 minus own target gene,
348 conventions). Per-gene delivered sign accuracy
a_g = Σ(match∧pool)/Σ(pool) over scored targets; support floor ≥10
(1,986 genes full; 324/306 per target-half). Features (universe 10,779):
189 θ_g, log basal CPM, gnomAD LOEUF/pLI/lof_z/mis_z (v2.1.1, max-exp_lof
transcript per symbol), GTEx v8 breadth (#tissues TPM≥1) and mean-log-TPM,
CpG-island promoter flag/count (GENCODE v44 TSS ±[2000,500] vs UCSC islands).
Per feature: Spearman vs a_g / θ_g / expression; split-target held-out
decile spread on half B vs 500-permutation null; OLS a_g ~ θ_g+expr+f
(ΔR², partial r). Best feature priced as reach key with 348 `walk_reach`
verbatim + 10k paired-target bootstrap. Self-caught bug: first stats pass
divided whole-matrix matches by pool support (a_g up to 58); caught by a
range check, fixed to pool-restricted numerator, all numbers below are the
re-run. Reach code never had the bug (always match∧pool).

## 1. Inventory: no per-gene Illumina-public feature is downloadable in ≤30 min / <5 GB

| source | URL (accessed 2026-09-11) | licence | per-gene bulk? |
|---|---|---|---|
| PromoterAI code | github.com/Illumina/PromoterAI (branch `master`; paper Jaganathan, Ersaro, Novakovsky et al., Science 2025, doi:10.1126/science.ads7373) | academic/non-commercial free; commercial via AI_licensing@illumina.com | **No** — precomputed all-promoter-SNV scores gated behind an Adobe-sign licence, download link "shared via email shortly after submission" (README). Not acquirable in-session. |
| PromoterAI/SpliceAI lookup | promoterailookup.broadinstitute.org (title: "SpliceAI / PromoterAI Lookup", served from GitHub Pages) | same academic terms | **No** — interactive per-variant lookup only, no bulk table. |
| PrimateAI-3D | primateai3d.basespace.illumina.com (+ github.com/illumina/primateai-3d) | academic/non-profit download after licence request (contact primateAI_licensing@illumina.com) | **No** — 70.6M missense variants ≈1.7 GB per-variant file, licence-gated; no per-gene summary published. |
| SpliceAI precomputed | basespace.illumina.com/s/otSPW8hnhaZR (via Illumina/SpliceAI README) | PolyForm Strict (code); CC BY-NC 4.0 (models/scores) | **No** — per-variant VCFs, tens of GB. |

Derivable-per-gene PromoterAI quantities (promoter sensitivity = mean/max
|score| over promoter SNVs; strong-element counts) all require the gated
table or running the TF2.15 model on hg38 — both out of budget. Used the
mandate's stand-ins instead (downloads → `data/`, 42 MB total, all <2 GB):
gnomAD v2.1.1 LoF-by-gene, GTEx v8 median TPM, UCSC CpG islands (via
hgdownload.soe mirror; main host TLS-blocked), GENCODE v44 basic.
Coverage of universe: θ 10,725; LOEUF 10,695; GTEx 10,518; CpG 10,596.

## 2. Results — association (units: accuracy fraction; spread = top−bottom decile on held-out half B)

| feature | Spearman vs a_g | vs θ_g | vs expr | spreadA | spreadB | null (500) | z | p (1-sided) | ΔR² over θ+expr | partial r |
|---|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|
| θ_g | +0.099 | +1.000 | +0.257 | +0.005 | +0.079 | +0.003±0.042 | +1.79 | 0.036 | — | — |
| log_expr | +0.057 | +0.257 | +1.000 | −0.029 | +0.040 | −0.001±0.041 | +1.01 | 0.148 | — | — |
| loeuf | +0.064 | −0.014 | −0.194 | +0.029 | −0.001 | −0.000±0.041 | −0.01 | 0.496 | +0.0037 | +0.061 |
| pLI | −0.023 | +0.084 | +0.205 | −0.054 | +0.026 | −0.000±0.042 | +0.63 | 0.266 | +0.0015 | −0.039 |
| lof_z | −0.057 | +0.012 | +0.173 | −0.006 | +0.045 | −0.000±0.043 | +1.05 | 0.168 | +0.0022 | −0.047 |
| mis_z | −0.039 | +0.008 | +0.148 | −0.053 | +0.053 | −0.002±0.039 | +1.42 | 0.072 | +0.0012 | −0.035 |
| gtex_breadth | +0.097 | +0.152 | +0.330 | +0.003 | +0.060 | −0.003±0.035 | +1.81 | 0.026 | +0.0073 | +0.086 |
| **gtex_meanlog** | +0.103 | +0.248 | +0.392 | +0.036 | **+0.140** | −0.001±0.043 | **+3.28** | **0.000** | +0.0067 | +0.082 |
| cpg_promoter | +0.047 | +0.015 | +0.071 | +0.014 | +0.058 | −0.002±0.022 | +2.68 | 0.004 | +0.0018 | +0.042 |
| cpg_count | +0.039 | +0.007 | +0.056 | +0.023 | +0.058 | −0.001±0.023 | +2.55 | 0.004 | +0.0010 | +0.031 |

Constraint (LOEUF/pLI/z) is dead: |spread| ≤ 0.05, signs flip across
halves. θ_g on the H1 pool is +0.079 (p=0.036, n.s. after correction) —
same "real but small" as 189 §3.3. Only gtex_meanlog survives Bonferroni
(×10): spreadB +0.140, 0/500 null ≥ it (two-sided ≈0.002–0.004 → ≈0.02–0.04
corrected). Reading: it correlates +0.39 with H1 expression while H1
log_expr itself is noisy single-line basal CPM — gtex_meanlog is best
understood as precision-weighted typical expression, not promoter biology
(partial r after controlling H1 expr is only +0.082).

## 3. Results — reach pricing (units: raw reach, 82 scored targets)

| scheme | mean | Δ vs O1 | Δ vs O4 |
|---|---:|---:|---:|
| O1 current | 0.25435 | — | −0.02640 |
| O4 \|design\| | 0.28075 | +0.02640 | — |
| θ_g desc | 0.19330 | −0.06105 | −0.08745 |
| expr desc | 0.19523 | −0.05912 | −0.08552 |
| FEAT gtex_meanlog (best) | 0.20031 | **−0.05404** | −0.08028 |
| FEAT then \|design\| | 0.20031 | −0.05404 | −0.08028 |
| shuffled-feature (3 seeds) | 0.183/0.211/0.211 | — | — |

Paired bootstrap (10k): FEAT ΔO1 = −0.0540 CI [−0.0964, −0.0137],
P(>0)=0.004; ΔO4 = −0.0803 CI [−0.1208, −0.0430]. The best association
feature ranks at shuffled-feature level — 189 §5's 90%-bar mechanism
exactly: decile-mean accuracy gains never approach 0.90 purity at the
prefix head, so k* does not move.

## Verdict against the PLAN.md kill criterion

Kill = (no feature beats expression's spread by +0.02) OR (Δreach < +0.01).
Leg 1 NOT fired (gtex_meanlog +0.140 vs expression +0.040, Δ=+0.100).
Leg 2 FIRED: Δreach = −0.054, CI entirely negative. → **NEGATIVE**.
No per-gene regulatory/constraint prior tested here carries
reach-usable information; the one surviving association (cross-tissue mean
expression) is a denoised-expression effect and prices below random keys.

## Limits

H1 development panel only, not A/B/C. Per-gene support is thin (mean pool
≈8 targets/gene; halves ≈4) — the ±0.04 null band is the floor, and only
gtex_meanlog clears it. PromoterAI-native per-gene features were not tested
(gated); proxies only. GENCODE v44 TSS vs H1-study genome build mismatch
would only blur the (already weak) CpG result toward zero.

## Exact commands to reproduce

```bash
uv run --no-sync python experiments/384-promoterai-gene-priors/run384.py ka    # KA row
uv run --no-sync python experiments/384-promoterai-gene-priors/run384.py build  # data/ -> work/features384.npz
uv run --no-sync python experiments/384-promoterai-gene-priors/run384.py stats  # §2 + work/stats384.json
uv run --no-sync python experiments/384-promoterai-gene-priors/run384.py reach  # §3 + work/reach384.json
```

Read-only inputs: `experiments/323-current-family-validation/work/`,
`experiments/189-per-gene-reliability/work/prior189.npz`,
`experiments/348-reach-ordering/reach348.py` (imported). Nothing outside
`experiments/384-promoterai-gene-priors/` written.

## Peak memory and wall time

`memory_pressure` free 86% at start. ka 2.9 s, build 6.3 s, stats 2.6 s,
reach 3.6 s. Peak RSS 0.77 GiB (`/usr/bin/time -l` on stats; largest object
85×10,779 float64 ≈ 7 MB). Threads capped at 4.
