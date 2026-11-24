# 473 — 410d2 (228's per-gene reliability gain m) re-priced on 446b's bulk: NO-GO. Sharpening q 1.75 already took it.

2026-09-26 23:05Z, third session, at the coordinator's request. Analysis only: no build, no builder edit, no `vcc`, no
`cell-eval2`. Preregistered in `PREREGISTRATION.md` (committed `8ae5c1ca` before any 446b-base number). Scripts:
`ka473_price414.py`, `abc473.py`, `arms473.py` (one PLOLO fold per process, `run473.sh`), `analyse473.py`.

## Answer
- **The primary arm A1 prices at −0.0008 overall vs 446b** (m clipped [1/2, 2], applied before sharpening, render-gated
  constant, both arms re-levelled; cancer folds, k 1). k 1.407 gives −0.0009. Pooled target-bootstrap interval
  [−0.0029, +0.0012]. Only 1 of 3 cancer folds is positive, and kolf2 reads −0.0045.
- The relative direction change is −0.5 % [−1.5 %, +0.7 %], against 414's +3.9 % … +5.1 % for the same m on 410c's
  unsharpened bulk.
- **GO items 1, 2 and 3 fail. Item 4 passes:** the permuted-m control loses on 3/3 cancer folds, −0.016 overall.
- **Mechanism: sharpening and m reweight the same genes.**
  - On 446b's bulk, m correlates with the gene's mean |h| at Spearman 0.43 (with control CPM also 0.43).
  - Sharpening at q 1.75 alone raises the energy share of the top-10 % m genes from 14.5 % to 24.4 %.
  - Placed before sharpening, m then acts like m^1.75 on the sharpened shape, which overshoots.
- **The shelved +0.0025 … +0.0046 (414) was a price on 410c.** The q 1.25 → 1.75 steps on the board have since spent
  it. Close 410d2 on the current bulk. A2 (unclipped) is also closed by the per-target row rule, as preregistered.

## Known answers (all before any arm number)
| # | row | expected | measured |
|---|---|---|---|
| KA-4 | 414's bounded-arm price, recomputed from 414's saved inputs (`price410d` arithmetic) | matched_overall +0.0045474 (cancer) / +0.0025158 (kolf2) | identical, \|d\| 0 |
| KA-a | input shas: stack dump `bb19f0e5…`, carriers dump `aed59b0b…`, m `da937cfa…`, ntclevel-446b `71e2e6d2…` | pinned | PASS |
| KA-b | the base's realised render energy vs 446b's measured (ntclevel) | within 1 % | **−1.5 / −1.6 / −2.6 %** (design route +4.5 / +3.8 / +3.2 %): **fails 1 %**, see below |
| KA-c | emulated ballast K vs 446b's build meta | median 6, max 16 | median 6, max 16, min 1: PASS |
| KA-414g | my gain code on 228's composition, `Kinf+gf@0.25` vs `arms414-<fold>.json` | 1e-12 | \|d\| 0 on 4/4 folds |
| KA-PP | compose_arm(defaults) == compose_pp; K3 weights == 418's W3 | bit-identical | PASS 4/4 |
| KA-420 | `stack\|K3+sh1.75` designed and realised (KMAX 16) vs `comp420-<fold>.json` | 1e-12 | \|d\| 0 on 4/4 folds |
| cliff | 187's model at X solved from 446b's board raw (0.6408807601) | reproduces the raw | X 0.044838, \|d\| < 1e-12 |

About KA-b:
- The builder's ballast runs in ratio mode: block totals = the carriers' delivered composition × p_bulk / p_de.
- 414's emulation instead clips S·p_bulk at the carriers' design totals.
- At q 1.75 (ballast K median 6, vs 4 at 417a, where 420 validated the route at 0.4–0.7 %), the two differ by ~2 % in
  level. The two routes bracket the measured energy, so every arm is read as a ratio to the base within one route. The
  price's anchor is 446b's measured energy.

## Render safety on the real A/B/C bulk (`abc473.py`; realised route, arm / base)
| arm | constant | energy A / B / C (non-panel) | max per-target (target) | SYVN1 | GO item 1 |
|---|---:|---|---|---|---|
| **A1** m clip [1/2, 2], pre-sharpen, rg | 0.6252 | ×0.939 / ×0.952 / ×0.937 | ×1.39 / ×1.49 / ×1.35 (TRPC4AP, LMBRD2) | ×1.02–1.06 | energy FAIL (safe side), rows OK |
| A2 unclipped | 0.5039 | ×0.895 / ×0.910 / ×0.903 | ×1.62 / ×1.66 / ×1.46 | ×1.09–1.41 | FAIL (energy, rows) |
| A3 clip [1/3, 3] | 0.5311 | ×0.914 / ×0.928 / ×0.911 | ×1.66 / ×1.68 / ×1.48 | ×1.09–1.15 | FAIL |
| A1post clip [1/2, 2] after sharpening | — | ×0.962 / ×0.964 / ×0.959 | ×1.21 / ×1.39 / ×1.47 | ×0.98–0.99 | energy FAIL (safe side) |
| C1 A1's m permuted | 0.7815 | ×0.983 / ×0.990 / ×0.991 | ×3.25 / ×2.99 / ×3.22 (MED25) | — | rows FAIL |

- No 410d-style blow-up in any clipped arm: the largest single-gene mass share rises at most from 0.014 to 0.016.
- **Every m arm lowers the context energy.**
  - The render-gated constant restores the effect-space energy.
  - m then concentrates it on expressed genes, where the minnorm mass fix compresses it.
  - The item-1 band was written to catch a blow-up past the cliff, so this is a failure on the safe side. The
    re-level refunds it.
  - **The re-level factors that refund A1's energy** (EV-optimal under 187's cliff sd, realised energy map):
    - ×1.049 / ×1.045 / ×1.054 against the base's own ×1.002 / ×1.005 / ×1.005;
    - i.e. 446b's factors A 1.1004 → ~1.152, B 1.1517 → ~1.197, C 1.1720 → ~1.229.
    - A1post needs ×1.023 / ×1.025 / ×1.029.

## Direction on PLOLO (`arms473.py`; realised, KMAX 16; Δ vs `stack|K3+sh1.75`)
| arm | fold | r = Δcos_agg (rel) | Δpds board-emul. (λ 0.5) | Δraw mse (re-levelled, k 1) | Δoverall k 1 [95 % target CI] | k 1.407 |
|---|---|---:|---:|---:|---|---:|
| **A1** | k562_gwps | −0.0089 | −0.00295 | +0.00373 | −0.0017 [−0.0043, +0.0009] | −0.0020 |
| | orion_hct116 | −0.0031 | −0.00345 | +0.00147 | −0.0015 [−0.0054, +0.0022] | −0.0016 |
| | orion_hek293t | −0.0024 | +0.00284 | +0.00117 | +0.0008 [−0.0033, +0.0049] | +0.0008 |
| | kolf2 | −0.0204 | −0.00822 | +0.00828 | −0.0045 [−0.0094, +0.0006] | −0.0051 |
| | **cancer mean** | pooled −0.0047 [−0.0149, +0.0065] | | | **−0.0008** [−0.0029, +0.0012], 1/3 + | **−0.0009** |
| A1post | cancer mean | pooled +0.0090 [+0.0021, +0.0161] | −0.0015 / −0.0005 / +0.0046 | | +0.0009 [−0.0004, +0.0022], 2/3 + | +0.0012 |
| | kolf2 | −0.0098 | −0.00483 | +0.00401 | −0.0025 | −0.0028 |
| A3 | cancer mean | pooled −0.024 | | | −0.0030, 0/3 + | −0.0036 |
| A2 | cancer mean | pooled −0.077 | | | −0.0070, 0/3 + | −0.0092 |
| C1 (control) | cancer mean | pooled −0.141 | | | −0.0165, 0/3 + | −0.0204 |

Price: overall = (2.2138 Δpds_bd − 1.051 Δraw mse)/6.
- Energy: 446b's measured A_c/B_c × the arm/base realised ratio (all genes but the own target, ntclevel's convention).
- X from 446b's board raw. The arm's X = X_base (1 + k r) √(A ratio).
- Base and arm are each re-levelled to their own EV-optimal per-context factor.
- Bootstrap: 4,000 target draws per fold, holding energies and factors at the point estimate.

## GO rule (preregistered) → NO-GO
| item | rule | A1 | pass |
|---|---|---|---|
| 1 | energy 0.97–1.03 every context, max row ≤ 1.5, mass share ≤ 2× | energy ×0.94–0.95 (below the cliff), rows ≤ ×1.49 | **no** (safe side) |
| 2 | Δoverall > 0 on ≥ 2/3 cancer folds and pooled r lower bound > 0 | 1/3; r lower −0.0149 | **no** |
| 3 | cancer ≥ +0.0015 and kolf2 ≥ 0 | −0.0008; kolf2 −0.0045 | **no** |
| 4 | C1 negative on ≥ 2/3 | 3/3 | yes |

**Reading:**
- A1post (the post-sharpening placement, a reference arm, not preregistered as the primary) keeps a small direction
  gain: +0.9 % cos pooled, with its interval above zero.
- It still prices at only +0.0009 (2/3 folds), below item 3's +0.0015. kolf2 is −0.0025, and its pds is not robust
  (−0.0015 / −0.0005 / +0.0046).
- It is not an arm. If a later bulk change lowers q, A1post is the placement to re-price.
- 486's public-data LOCO test of a LEARNED per-(target, gene) gate is a different object (Testh's), and this result
  does not close it.

## Environment
- `run473.sh` started each fold only once vm_stat free+inactive was ≥ 33 GiB (it read 35–36 GiB, 22:52–22:53Z, between
  the coordinator's builds).
- Four folds ran strictly serially: 13–18 s each, peak RSS 5.2–5.8 GB.
- `abc473.py` peaked at 0.95 GB, `ka473_price414.py` at 0.1 GB. Nothing was written outside this directory.
