# 468 — `VCC_DE_PRIOR_GATE`: the strength-gated programme prior in the production builder

**Status (2026-09-26 18:36Z):** builder `693fde15` committed (dc05673e) after KA-1 and KA-2 passed byte for byte.
458g is built and gated for 09-30 s2 (package `eb502099…`). `compare458g.py` is running.

## The knob
- `VCC_DE_PRIOR_GATE=q<NN>` takes NN from 1 to 99 and needs `VCC_SAMELINE` without `noprog`.
- In each context whose default block prior the same-line programme replaced, the builder first computes every target's
  design strength: the number of genes whose |gated design row| > 0.2. That row is the emitter's own `key`, computed by
  the same `emission_rank_axis` call.
- Targets below the context's NN-th percentile get the context's CLEAN4 prior (`de_sign_prior_full`) through the
  per-target prior map (`sl_prior`). A weak target that already carries a same-line per-target prior (`pt`) keeps it.
- The build prints `[prior_gate]` per context. The meta's `config.de_prior_gate` records q, tau, and per context the
  threshold, the counts and the strength quantiles.
- Unset: no code path runs and no config key is written.

## Tests (`chain468.sh`, 16:56–16:57Z)
- `test_gate468.py`, all pass:
  - unset → None;
  - VCC_SAMELINE alone → None;
  - q33 / " Q50 " / q1 with pt parse;
  - refused: q0, q100, x33, 33, q3.3, noprog, and the gate without VCC_SAMELINE.
- 451's VCC_SAMELINE tests: all pass.
- Builder suite: 56 passed, 0 failed, 1 skipped. 417's four unit tests: pass.

## Known answers
| row | env | reference | result |
|---|---|---|---|
| KA-1 | every new knob unset (446b's exports) | 446b `9559923a774eb248…` | **PASS, byte-identical** (17:32Z) |
| KA-2 | + `VCC_SAMELINE=arch1d_programme_abc.npz:6951dbc7…:0`, no gate | 458d `1d09cc91a8d754c5…` | **PASS, byte-identical** (18:03Z) |

## The arm: 458g (built 18:03–18:34Z)
458g = 458d + `VCC_DE_PRIOR_GATE=q33` (`PREREGISTRATION.md`), built by `run_arm_generic.sh` (`build_458g/run458g.log`).
- **Config diff vs 458d: exactly `de_prior_gate`.** The postbuild's single FAIL is the usual config-drift line.
- **The gate on A/B/C:** threshold 2 in every context, so a target with at most one gated design gene above 0.2
  keeps CLEAN4. That is the same cut as on H1 (strength ≤ 1; H1's q33 threshold 1.72).

  | context | programme (strong) | CLEAN4 (weak) | strength quantiles (0, 25, 50, 75, 100 %) |
  |---|---:|---:|---|
  | A | 214 | 86 | 0, 1, 3, 7, 466 |
  | B | 217 | 83 | 0, 1, 3, 7, 424 |
  | C | 217 | 83 | 0, 1, 3, 7, 459 |

- Zero-gene audit: 0 zero blocks (lfc < −5: 17; > +3: 3, max +4.11).
- Cliff check OK (A 0.9603, B 1.0282, C 1.0785).
- Artifact `9929218c…`, package `eb502099d2ff69480d5d6984188f4f8840be007f8aa2e29dab15017e0601d39c`, novelty rc 0.
- `chain0930c.sh`'s readiness check passes for both 458d and 458g, so 09-30 = s1 458d, s2 458g.

## The pair's noise model: an exact three-way split (`compare458g.py`, 18:36–18:37Z)
The script restored 446b, 458d and 458g from their packages and compared all 900 (context, target) blocks, CSR data
and indices.

| context | blocks | = 458d | = 446b | = both | = neither | gate weak / strong |
|---|---:|---:|---:|---:|---:|---|
| A | 300 | 214 | 87 | 1 | **0** | 86 / 214 |
| B | 300 | 217 | 83 | 0 | **0** | 83 / 217 |
| C | 300 | 217 | 83 | 0 | **0** | 83 / 217 |

- Every 458g block equals 458d's (the strong targets) or 446b's (the weak ones). In A, one strong target's block is
  identical in all three packages: the two priors gave it the same emission.
- **The emitter's random draws do not depend on the block prior.** So the three board entries decompose the
  programme's effect with no emission re-roll noise between them:
  - 458g − 458d is the prior change on the weak third alone;
  - 458g − 446b is the programme on the strong two thirds alone.
- The ±0.009 widening in the preregistration (for shifted draws) does not apply. The bins stand as written.
- Record: `work/compare458g.json`, which lists the weak targets per context. The restored files were deleted.
