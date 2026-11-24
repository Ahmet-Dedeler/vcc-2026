# 212 — method, gates, and what each control reads when the mechanism is absent

Audit trail for `RESULTS.md`. Every number there comes through the pipeline below.

## The question, and why it is not 170's question

170 swept `WHITEN_GAMMA` and priced it on `mse` (through 096's `claim` route) and `pds`, with
`VCC_ALPHA` **pinned at 0.982** while gamma moved the predicted effect's mean square `A` by
14%. `mse`'s numerator is

```
num_c(s) = s² A_c + Kc_c − min(claim_c, s² B_c) − 2 s X_c
```

so a gamma move at fixed `s` is a gamma move **confounded with an amplitude move**, and above
the ρ cliff the metric refuses to refund the amplitude. 170's V, with its vertex on the shipped
0.65, is that confound measured correctly.

The transfer cosine `cos = X / (√ΣD · mean_c √A_c)` is exactly invariant to a global rescale.
At the constrained optimum in `s` the amplitude cancels, and `mse` becomes a monotone
decreasing function of `cos` **and of nothing else**. §8.2 says the cosine is the entire
remaining gap. So this experiment asks 170's question on the endpoint the gap is measured in,
at each gamma's own alpha, and on the emission that is actually on the board.

## Pipeline, stage by stage

| stage | what runs | file |
|---|---|---|
| pre-whitening effects | 170's cached `candidate_pre`, sha re-asserted against a15's shipped audit | `gate212.py` |
| operator | shipped `apply_v9_whitening`, `union_pruned`, gamma swept, `lam = 0`, 15 gammas | `gate212.py` |
| CD4 residual | shipped `apply_frozen_cd4_residual` / `add_post_whiten_residual`, re-solved per gamma | `gate212.py` |
| challenge `A`/`B`/`claim`/ρ/Ĉ | 170's `chal.py`, **a15 geometry** (the incumbent's), 15 gammas | `chal212.py` |
| the two cliffs, three emissions | shipped `emit_counts(mode="quasi")` and `rc_emit.rc_hybrid_block` (± `CallSet`) | `cliff212.py` |
| **cosine (primary)** | 162's PLOLO, 4 folds, challenge frames A/B/C, both renders | `plolo212.py` |
| Δcos / Δ`pds` intervals | 204's 4,000-draw paired target-cluster bootstrap | `boot212.py` |
| six members on h1-test | 170's `h1sweep.py`, imported not copied, five new gammas | `h1sweep212.py` |
| alpha re-solve + pricing | 096's identity, `min()` evaluated per context at every `(gamma, s)` | `alpha212.py` |
| falsifiers | permuted-operator and spectrum-shuffle controls on the same instrument | `falsify212.py` |

Nothing is reimplemented. `plolo212.py` and `falsify212.py` import 162's `plolo` and 204's
`score204` and call them; the only change is that `gamma` is an argument instead of a module
constant (`score_at` is `score204.score` with `PROD["gamma"]` replaced).

## Known-answer rows (rule 10's corollary). Every one of them is 170's or 204's.

| # | gate | source | result |
|---|---|---|---|
| 1 | `candidate_pre` sha256 = `0601336f…` | a15's shipped audit | **MATCH** |
| 2 | a15 post-whitening audit, γ 0.65 / λ 0 — 8 scalars | a15's meta | **PASS**, worst rel 1.1e−15 |
| 3 | a19 post-whitening audit, γ 0.80 / λ 0.007 — 6 scalars | a19's meta | **PASS**, rel 0.00e+00 |
| 4 | the 20 shared effect arrays vs 170's cache | 170 `work/eff/` | **20/20 BIT-IDENTICAL** |
| 5 | `msekit` on the frozen 095 h1-test artifact = 2.0551018299736388 | cell-eval2's own `agg_results.csv` | `h1sweep212.py` |
| 6 | per-context `Kc` = 0.005084 / 0.005228 / 0.005919 | 096 §(e) via 170 gate 7 | **PASS**, rel ≤ 5.7e−5 |
| 7 | challenge cap = 0.24653; `claim`(γ 0.65) vs 164's 0.24581 / 170's 0.24610 | 164, 170 | `chal212.py` |
| 8 | PLOLO `equal` at γ 0.65, four folds, `pds`/`cos`/`q`/`q_w` | 204's published rows | **PASS**, rel 0.00e+00 on 11 of 11 |
| 9 | PLOLO `equal`'s `max_abs_dev_from_162_pre` per fold | 204's own audit | **PASS**, exact on 4 of 4 |
| 10 | 178's PLOLO Δ`pds` at γ 0.45 (−0.00930) and γ 0.30 (−0.02386) | 178 §7 | see RESULTS |
| 11 | Δcos bootstrap half-width ≈ 4.3e−5 / 0.20% relative | 204 §8.3 | **3.90e−5 / 0.178%** |
| 12 | 208's Ĉ/C_p: `quasi` 1.125, `rc_callset` 2.940 | 208 §2.4 | `cliff212.py` |

**Gate 4 is the one 170 could not run** and it is the strongest: the whole 300 × 18,533 effect
array at nine gammas and both surfaces, reproduced bit for bit from a cache built two days
earlier on the same machine, which pins the operator, the CD4 solve and the gene axis together.

## Rule 11 — path and sha256, never a description

| operator | path | sha256[16] | max eig | n |
|---|---|---|---:|---:|
| **production** | `data/donors/whitener_union_pruned.npz` | `e2f65ecf69b7e0b0` | 5.894570 | 17,942 |
| PLOLO fold `k562_gwps` | `experiments/062-cd4-post-whiten/cloud-output/whiteners/held_k562_gwps.npz` | recorded in `work/plolo212_ka.json` | | |
| PLOLO fold `orion_hct116` | `…/held_orion_hct116.npz` | recorded | | |
| PLOLO fold `orion_hek293t` | `…/held_orion_hek293t.npz` | recorded | | |
| PLOLO fold `kolf2` | `data/donors/whitener_union_pruned.npz` (production; `kolf2` contributed no row) | `e2f65ecf69b7e0b0` | 5.894570 | 17,942 |

The production operator's basis is **29,394 × 17,942** (193's correction; §4.2's stated
12,000 × 7,760 matches no file on disk). `optable212.py` tabulates every operator on disk.
The sha is asserted **before** the file is loaded, in `gate212.py`, `plolo212.py` and
`falsify212.py`, because 170's own warning stands: the identity check at production settings
reproduces the baseline by orthonormality and therefore cannot detect a mismatched basis.

## Controls, each with its reading when the mechanism is absent

| control | with the mechanism | without it |
|---|---|---|
| γ = 0 under the permuted operator | must equal the real operator exactly (γ = 0 is the identity, so the permutation is a no-op) | anything else ⇒ the control is not a control |
| `gene_perm` `pds` at γ 0.65 | must **collapse** — this is §8.3's Oct-22 trap made deliberate | if it does not, `pds` is not reading the operator and the screen is void |
| `gene_perm` composite anywhere | must never beat the real operator | 209's outcome: the random control was the maximiser over all 36 arms |
| `spectrum_shuffle` | same gene pairing, same eigenvalue multiset, variances detached from directions | isolates *which* eigenvalue goes with *which* direction |
| PLOLO `equal` at γ 0.65 | bit-equal to 204's published row | any drift ⇒ not 204's instrument |

## What this experiment deliberately does NOT do

* It does not refit the operator's basis (193 closed the eigenvalue map, 195 closed learned
  transforms) and does not touch the donor tables (194, 204, 209).
* It does not size anything on the Orion **pair** (168d's standing rule). The two Orion lines
  enter only as two of PLOLO's four folds, where each is a held-out truth with its own
  leakage-free operator, which is a different construction from 097's Orion harness.
* It does not use `data/donors/h1_2025` (85/85 coverage of h1-test's targets, pure leakage).
* It sweeps gamma alone at `lam = 0`. 193/198 established production ships `WHITEN_LAM = 0.0`,
  so there is no ridge; a19 moved both and is the reason its slope is unusable.
* No `vcc`, no `cell-eval2`, no build, nothing written outside
  `experiments/212-gamma-cosine/`.

## Two operational notes

1. **`h1kit.builder()` cannot express the incumbent's geometry.** It sets
   `VCC_MASS_FIX=minnorm` at module scope and *raises* if the mass fix is off, while
   `rc_callset` ships without it — the production builder itself refuses `VCC_MASS_FIX`
   together with `--emitter rc_hybrid|rc_callset`
   (`scripts/make_submission_pm_cd4.py:2755`). `chal212.py` therefore loads the production
   builder directly. 170's whole challenge-panel column is in a23's geometry; this
   experiment's is in a15's, and both are reported.
2. **162's and 204's fold caches are CLONED into this experiment's `work/plolo/`** (APFS
   `cp -c`, byte-identical) before `plolo.WORK` is redirected, so no code path can write
   into another experiment's directory even if a cache were to miss.
