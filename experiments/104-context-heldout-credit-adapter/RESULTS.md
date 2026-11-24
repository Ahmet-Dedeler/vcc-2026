# 104 — context-held-out direct-credit adapter

**Verdict: NO-GO.**

The preregistered full rank-8 STRING/context model improved linear signed
credit and complete-gallery PDS over identity on all three held contexts.
However, it reduced top-20% sign precision on H1 and failed to beat the
response-label permutation maximum on H1 and HCT116. The full model therefore
does not satisfy the all-surface direction-safety gate and must not proceed to
count generation or a challenge build.

No AWS worker, KOLF artifact, `cell-eval2`, count generation, challenge
package, submission, or commit was used.

## Authoritative run

`work-v2/` is the authoritative run. All three folds use one implementation
and README hash.

An earlier `work/` integration run is non-authoritative. It scored HCT116 and
HEK293T before H1 exposed two implementation issues: the production-pruned
operator required its own metadata loader, and scikit-learn's incremental
weighted variance emitted negative-roundoff warnings. The implementation was
repaired before any H1 score, redirected to a fresh `work-v2/` namespace, and
all three folds were rerun. The old artifacts are retained rather than
silently deleted.

## Effect-level results

All values below are candidate minus comparator. Confidence intervals are
2,000-draw whole-target percentile bootstraps.

| surface | full − identity credit | full − agreement credit | full − agreement top-20% precision | full − identity PDS |
| --- | ---: | ---: | ---: | ---: |
| H1 | **+0.004496** `[+0.002443,+0.006895]` | +0.000154 `[-0.000987,+0.001499]` | **−0.041424** `[-0.070703,-0.013736]` | +0.002521 `[-0.000560,+0.006583]` |
| HCT116 | **+0.008110** `[+0.005921,+0.010391]` | **+0.003769** `[+0.002271,+0.005404]` | +0.018582 `[-0.007603,+0.046008]` | +0.006130 `[-0.001679,+0.014517]` |
| HEK293T | **+0.003890** `[+0.002472,+0.005575]` | **+0.002748** `[+0.001639,+0.003965]` | **+0.064843** `[+0.042013,+0.087406]` | **+0.007807** `[+0.002660,+0.013532]` |
| macro point estimate | **+0.005499** | **+0.002224** | **+0.014000** | **+0.005486** |

The candidate therefore found real continuous-credit signal. It is not a
safe confidence ordering: on H1, the genes promoted most strongly were
materially less sign-accurate than those promoted by the agreement-only
ablation.

## Permutation nulls

Each null independently permuted signed response labels within every complete
target block and sample-weight quintile, then refit the full model.

| surface | observed full − agreement credit | null maximum | empirical p |
| --- | ---: | ---: | ---: |
| H1 | +0.000154 | +0.000831 | 0.3810 |
| HCT116 | +0.003769 | +0.004579 | 0.1429 |
| HEK293T | +0.002748 | +0.002082 | 0.0476 |

Only HEK293T beats all 20 nulls. The cross-context aggregate cannot rescue the
failed all-surface gate.

## Invariants

All authoritative folds passed:

- all 300 VCC panel coordinates are zero for PDS;
- directional credit excludes only the row's own target, aside from three H1
  genes absent from the challenge's 18,533-gene axis;
- the transformation never creates a nonzero coordinate;
- no sign changes;
- row-level L2 norm restoration;
- agreement-only and full models have identical target/gene support;
- the held real donor is absent from its baseline and training source sets;
- every artifact was frozen and hashed before that surface's perturbation
  truth was read.

## Resource use

Peak guarded RSS:

| stage | H1 | HCT116 | HEK293T |
| --- | ---: | ---: | ---: |
| freeze | 12.14 GiB | 6.06 GiB | 6.48 GiB |
| score | 4.16 GiB | 2.77 GiB | 3.19 GiB |

## Interpretation

This result sharpens experiment 100 rather than overturning it. Optimizing the
continuous signed-credit target recovers a reproducible NMAE/PDS improvement,
but STRING/context features still fail as a universal replacement for
donor-agreement confidence. The H1 precision loss is large and statistically
resolved.

Do not tune beta, ridge strength, STRING rank, target subsets, or null seeds
against these outcomes. A new attempt would need a genuinely new external
confidence feature, not another loss or threshold applied to this feature
family.

The preregistered macro cluster bootstrap and count-level call-coverage checks
were not run because earlier mandatory gates had already failed. They cannot
reverse the NO-GO.
