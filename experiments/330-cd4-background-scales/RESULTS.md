# 330 — Small favorable point estimates, unresolved; no advancement

September 9, 2026. **DOES_NOT_EARN_FURTHER_DEVELOPMENT**.
Protocol/code/input freeze: **57e12aa7**. No production or package change.

Replacing the compact panel's CD4 gene-scale estimates with estimates from
1,000 fixed non-panel source targets moves every raw metric favorably at the
point estimate. However, **every primary member interval includes zero** and
both required overall intervals include zero. It does not earn advancement.
Do not retune its background size, floor, gene set or blend to pass.

## Actual-count readout

Background1000 minus cap15, reused 85-target H1 development panel:

| member | raw delta | paired target 95% CI |
|---|---:|---|
| PDS | +.00182073 | [-.00042017, +.00406162] |
| MSE, lower better | -.00009680 | [-.00053816, +.00035309] |
| NMAE, lower better | -.00011870 | [-.00113223, +.00102443] |
| fidelity | +.00015810 | [-.00014929, +.00046549] |
| reach | +.00036808 | [-.00430690, +.00602440] |
| Jaccard | +.00003304 | [-.00001913, +.00008633] |

H1-local scaled overall: **+.00083497 [-.00030217, +.00216558]**.
Cap15 .28095315; background1000 .28178811.
All arms' point scaled MSE is zero; raw MSE remains reported.
10,000 paired target draws, fixed gallery/models/anchors; panel MSE and its
global claim cap recomputed in each draw.

Controls:

- Versus background80 overall: **+.00019213 [-.00131979, +.00160961]**.
  No established overall benefit from the larger background sample.
- Versus gene-shuffled scale overall:
  **+.00445940 [-.00185254, +.01201348]**, unresolved.
  NMAE/fidelity beat that control, but this does not establish baseline gain.
- Target derangement decreases PDS by **.3894958** (.8914566 to .5019608);
  negative control passes.
- The sample-size point condition and profile non-harm gates pass. Both
  required overall-improvement gates fail. No candidate is promoted.

This is a scoped unresolved normalization experiment, not proof that panel
composition is irrelevant or that source-scale estimates are identical.

## Source and construction checks

- 1,000 fixed source targets × 9,751 response genes × 3 conditions.
  All coefficients finite. Both panel target sets excluded from background
  targets and response coordinates.
- **768,000** already-stored coefficients (1,000 × 256 × 3) reproduced exactly
  against experiment 328's archive before interpreting broader extraction.
- Raw panel coefficients unchanged. Only registered source gene-scale
  denominators changed; panel centering and median-normalization rules stayed.
- 80 carried / 77 complete CD4 targets throughout; H1 target coordinates and
  incomplete residual rows stayed zero through deflation.
- Identity CD4 arrays, U5 hash, final effect hash, emitted-count digest and
  panel-total digest match 323 exactly. Five raw metrics reproduce exactly;
  NMAE differs only **1.11e-16**, within the frozen tolerance.
- Seven pure unit tests pass.

Final emitted non-target pseudobulk, background1000/cap15:
row-norm ratio min/median/max **.998754 / 1.000463 / 1.017154**;
relative Frobenius change **.0345001**; **5,188** signs changed among nonzero
baseline coordinates. These are descriptions, not identified corrections.

## Reproducibility and resource record

All eight execution receipts are complete. Guarded process exited 0 after
**779.84 s**, sampled tree peak **5.83 GiB**, macOS maximum RSS
**6,327,730,176 bytes**. No surviving scorer/guard was found.
An independent read-only calculation rebuilt all five overall point scores
from raw metric files and anchors and verified all overall contrasts to 1e-12.

- Input manifest:
  `e88ac0f97dcd350a0fdd8a82d984eee944d0a58a981c0b06ceb013f77ab27de1`.
- Scale archive:
  `e0e8a2afd77860ef8c4852b43241f84031af52985be1abe366c11ec34165e22d`.
- Terminal results:
  `385e0826397d36faec286277f9d9e052572733913b0d89b6fa73b78252e32738`.

No new reserve, response download, paid service, package, submission or push.
No A/B/C performance claim follows from this one reused H1 context.
