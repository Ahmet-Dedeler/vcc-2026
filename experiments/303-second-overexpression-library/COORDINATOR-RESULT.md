# C7 result — modest observed reversal, not a universal inverse

**C7a PASS; C7b/C7c/C7d FAIL.** Data validity passes. The frozen protocol and
complete code were committed as `45788c6` before access.

R5 was read once at **2026-09-07 06:01:00 UTC / 09:01 Ankara**. All 510
source profiles pass the CPM reconstruction, mask, identity and finite-value
checks. R5 is now **SPENT**. The read took 3.9 seconds, peak RSS 0.41 GiB.
No other reserve or source row subset was opened.

## Frozen test

100 matched targets, 15,211 downstream coordinates, common-cohort centering,
2,000 uniform matching permutations and a four-cell Holm family.

| question | observation | Holm p | verdict |
|---|---|---:|---|
| C7a: antisymmetric observed matching | median cosine **−.02633**; conditional target-bootstrap 95% **[−.05203, −.00949]** | **.001999** | PASS |
| C7b: conservation predicts reversal | rho **−.12883**, opposite the hypothesis | 1 | FAIL |
| C7c: geometry predicts reversal | rho **+.07606** | 1 | FAIL |
| C7d: geometry predicts absolute overlap | rho **+.08070** | 1 | FAIL |

64/100 matched cosines are negative. **41** are below −.05 versus **13**
above +.05; tail excess .28 [.14, .41]. Both C7a components exceed all
2,000 matching-null draws in their predeclared directions, so each raw
Monte Carlo p is 1/2001, not zero. Null median-cosine mean is +.00607.

The three association tests have 97 feature-complete targets. Their
non-rejections do not establish that no predictor could work.

## Descriptive robustness, not additional discoveries

The median remains negative without centering (−.02029), under bank-wide
centering (−.03705), and under the draft's pairwise own-gene exclusion
(−.02554). These were declared sensitivities.

Own-transcript signs are positive for 84/98 OE targets and negative for
98/98 KD targets. The own coordinates are excluded from the primary score.

C7 S correlates with 302's S_g at +.426 on all 58 overlapping scored targets
and +.502 on its 45 primary targets. **Both statistics share the same KD
bank**, so this is descriptive agreement, not independent biological
replication or an OE reliability estimate.

## Interpretation

This resource supports a small **target-matched negative alignment**, which
302's different assay did not detect. A median −.026 is far from a perfect
inverse (−1), and 36 targets are non-negative. Line, assay state, timescale,
platform and measurement error remain confounded. The intervals and
permutation test are conditional on these rows and target-label
exchangeability; no family-held or study-level replication was established.

The result does not rehabilitate C6, establish why the datasets differ, or
make raw overexpression a general knockdown predictor. It is a bounded
scientific finding, not a newly cleared competition method or novelty claim.

## Independent verification and artifacts

`scripts/audit_c7.py` independently reconstructs the cosine matrix, all
2,000 permutations, feature mappings, four p values and Holm results without
importing the experiment's analysis code. Matrix difference **0.0**.
`work/COORDINATOR-VERIFICATION.json` binds the checked outputs.

- Sufficient product `work/nak303.npz`:
  `40065fadf89afbca554cc3f88944d9269a05be087e0f92b2e7138738b90ad599`.
- Frozen manifest canonical digest:
  `7298dfe5258fc1b5bafe5576b7317178fa0d7f9233cfd9041625aa452de4bc48`.
- Complete result: `RESULTS.md`, `work/stats303.json`,
  `work/symmetry303.csv`, `work/analysis303.npz`.
- Access provenance: `work/RELEASE.json`, `work/ACCESS-MARKER.json`,
  `work/access303.json`. No new source read occurred during verification.

No production change, build, submission or external publication.
