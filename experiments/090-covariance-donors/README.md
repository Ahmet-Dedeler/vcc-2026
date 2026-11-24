
## Confound check, resolved at launch

The concern flagged in the commit was that adding donors widens the whitening
axis as well as enriching the covariance, so a positive result could not be
attributed. The observed basis shapes settle it:

| arm | basis rows | axis genes | Δrows | Δaxis |
| --- | ---: | ---: | ---: | ---: |
| baseline (2 donors) | 27,764 | 17,942 | — | — |
| +nadig_jurkat | 30,034 | 17,970 | +2,270 (+8.2%) | +28 (+0.16%) |
| +rpe1 | 29,968 | 17,966 | +2,204 (+7.9%) | +24 (+0.13%) |
| +k562_essential | 29,785 | 17,949 | +2,021 (+7.3%) | +7 (+0.04%) |
| +all three | 34,259 | 17,984 | +6,495 (+23.4%) | +42 (+0.23%) |

The axis moves by at most `0.23%` while the basis grows by up to `23.4%`. The
treatment is therefore ~100x more covariance enrichment than axis widening, and
any effect is attributable to the covariance estimate rather than to whitening
more genes. This is a cleaner separation than the design intended, and it holds
because the three added donors are broad common-essential screens measured on
essentially the same transcriptome — they contribute many perturbation rows and
almost no new gene coverage.
