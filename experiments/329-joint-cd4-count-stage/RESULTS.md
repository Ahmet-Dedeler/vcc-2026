# 329 — Source-density improvement does not earn count-stage advancement

September 9, 2026. **DOES_NOT_EARN_FURTHER_DEVELOPMENT**.
Frozen construction: **7fd6e8dc**. No package, production change or submission.

The SE-aware joint posterior improves source predictive density (328), but
its current-family H1 count prediction does not establish an overall gain.
Relative to cap15, NMAE and fidelity have clearly harmful paired intervals.
This closes the registered construction, not all multivariate estimation.
Do not tune its dose, prior, normalization or deflation until it passes.

## Actual-count result

All numbers below are **joint minus cap15** on the reused 85-target H1 panel.
Intervals use 10,000 paired target draws, fixed gallery/anchors/models.

| member | raw delta | paired 95% CI |
|---|---:|---|
| PDS | +.0047619 | [-.0102241, +.0184874] |
| MSE, lower is better | +.0070755 | [-.0104847, +.0277419] |
| NMAE, lower is better | **+.0138885** | **[+.0019423, +.0334234]** |
| fidelity | **-.0016449** | **[-.0029132, -.0004718]** |
| reach | -.0158834 | [-.0408734, +.0045561] |
| Jaccard | -.0000826 | [-.0005055, +.0002975] |

H1-local scaled overall: **-.0056597 [-.0169085, +.0028336]**.
Cap15 overall .2809531; joint .2752934. These are local H1-anchor scores,
not estimates of the A/B/C board score.

All five arms have point scaled MSE **0**, so the raw MSE difference must
remain visible. The bootstrap recomputes the panel ratio and global
jackknife-credit cap, rather than averaging nonexistent per-target MSEs.

Controls:

- Joint minus diagonal overall:
  **+.0012919 [-.0047792, +.0091763]**, unresolved.
- Joint minus shuffled-fit overall:
  **+.0018493 [-.0058937, +.0113686]**, unresolved.
- Joint's raw NMAE beats both fitted controls, but both controls are worse
  than cap15; superiority to them is not baseline improvement.
- Deranging joint's target rows reduces PDS by **.3893557** (.8943978 to
  .5050420). The negative control passes.
- All three overall-advancement intervals fail; the no-resolved-NMAE-harm
  gate fails. No alternative arm is promoted.

## Construction and delivered geometry

All **772,464 eligible triples** (77 complete targets × 10,032 genes) were
estimated. Source coefficients outside eligibility retained their original
deltas. Prepared CD4 support/count arrays stayed exact: 80 carried, 77 complete,
with incomplete rows and H1-target coordinates zero through deflation.

The preprocessing rule, k=5, f=.34 and main-donor/PM/whitening/emission settings
were fixed, but estimated preprocessing and deflation objects changed as
registered. No claim of delivered amplitude invariance is made:

- Emitted non-target pseudobulk row-norm ratio min/median/max:
  **.780161 / .975730 / 1.738820**.
- Relative Frobenius change: **.334300**.
- **53,520** emitted signs changed among nonzero baseline coordinates.
  They are sign changes, not established corrections.

Source predictive density is therefore insufficient to select this transfer
construction. The experiment does not identify whether model misspecification,
source-to-recipient biology, preprocessing, or their interaction caused the
loss; no causal mechanism is claimed from this single comparison.

## Verification and resources

- Ten unit tests pass, including exact Python-float/float32 extractor
  arithmetic under NumPy 2.5.2.
- Identity conversion reproduces every compact-bank delta exactly; prepared
  CD4 arrays, baseline deflation U hash and final effect hash are exact.
- Cap15 emitted count and panel-total digests match 323; **all six raw
  metric differences are 0.0**.
- Every stage has a completed execution receipt. Guarded process exited 0;
  733.83 s total, sampled tree peak **5.13 GiB**, macOS maximum RSS
  **5,505,925,120 bytes**. No surviving scorer/guard process was found.
- Independent read-only arithmetic recomputed all five point overall scores
  from raw metrics/anchors and verified every overall contrast to 1e-12;
  separately reconstructed the primary PDS and four DE paired intervals
  directly from saved target rows.
- Input manifest SHA-256:
  `2c2c360ba249a7b3b1b3425d276b283bc68480bc30f6acbbb5126d0192e68e71`.
- Terminal result SHA-256:
  `434d6e05edf3eb5c675617203dc846df67ddd816a0500628b10df142d4b6c40b`.

No reserve, paid computation, external response download, package or submission.
329's source/fit/candidate files remain available for reproducibility.
