# Experiment 141 results

Status: **terminal raw-gate NO-GO**

Date: **2026-08-30**

## Bottom line

The Pan same-target residual is strongly useful on the reused 613-target
production-shaped surrogate, but the new gene-coordinate calibration is not.
Nearly the entire gain comes from the scalar/direct-blend ancestry already
tested by Experiment 117.

The fixed primary `diag_rho050` improved:

- macro raw PDS by `+0.0119714643`;
- macro signed credit by `+0.0076031793`.

However, `scalar_control` already improved:

- macro raw PDS by `+0.0118417405`;
- macro signed credit by `+0.0076232622`.

Therefore the primary's incremental PDS over the scalar control was only:

```text
+0.0001297238
```

against the preregistered `+0.001` distinctness gate. The raw gate failed and
the conditional null/bootstrap stage did not run.

## Cross-fitted results

| Candidate | Macro PDS delta | Macro signed-credit delta |
|---|---:|---:|
| scalar control | +0.0118417405 | +0.0076232622 |
| diagonal rho 0.25 | +0.0119505841 | +0.0076696225 |
| **diagonal rho 0.50 primary** | **+0.0119714643** | **+0.0076031793** |
| diagonal rho 0.75 | +0.0119919003 | +0.0074195448 |

Primary family-level PDS deltas:

- K562: `+0.0225999318`;
- RPE1: `+0.0023190353`;
- Jurkat: `+0.0109954259`.

Primary family-level signed-credit deltas:

- K562: `+0.0133846797`;
- RPE1: `+0.0077485144`;
- Jurkat: `+0.0016763438`.

All absolute effect gates passed. The sole failure was the required evidence
that coordinate calibration added materially beyond a scalar Pan blend.

## Interpretation

This experiment does not overturn Experiment 117. It explains why the reused
source-only development surface looks attractive: a learned scalar move toward
Pan is broadly positive there. But the diagonal coefficients add only
`1.30e-4` PDS and slightly reduce signed credit relative to the scalar.

The new mechanism therefore has no demonstrated value. The strong scalar
result is not promotable because direct Pan interpolation, dose tuning, and
Pan-derived confidence were already consumed and closed by Experiment 117.

Do not run the null stage, open the 410-target split, reconstruct the exact-a15
comparator, tune rho, or build a submission from this family.

## Integrity

- eight synthetic/contract tests passed;
- Ruff passed;
- two independent static reviews completed before development access;
- Pan centering used 1,331 non-gallery rows and excluded all 1,023 gallery
  targets by name, including the sealed 410;
- all Experiment 140 identity members, Pan compact-donor members, and
  development truth inputs were hash-validated;
- peak process RSS: `1.261 GiB`;
- no `cell-eval2`, count generation, challenge build, submission, commit, or
  push occurred.

Artifacts:

- result SHA-256:
  `d6386c4a1281d0cfb31ea5c7c079b5a81fbf5068c328cf25a3eae07fd1f7b874`;
- per-target arrays SHA-256:
  `35090a88eb2ab1be2c788f9f3421821587e7064a49f0be4273d35700fe08d22b`;
- artifact manifest SHA-256:
  `510cb2c9907f338420534e145edb76467bce074a197032e96afb6944e582a34e`;
- development freeze SHA-256:
  `adfc59ec0918f77b37c490e297660b3d81b115035efb9edfd80231903c2ea2ae`.

