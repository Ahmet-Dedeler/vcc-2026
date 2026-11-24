# Experiment 145 — spectral trace-bin attenuation

Decision: **TERMINAL NO-GO**

Date: **2026-08-31**

## Question

Does the incumbent gamma-0.65 whitener over-amplify the low-variance half of
the source-covariance spectrum?

The frozen family preserved every incumbent row norm and attenuated the third
and fourth cumulative-trace bins with:

```text
beta = 0, 0.25, 0.50, 0.75, 1.00
multiplier = [1, 1, exp(-beta/2), exp(-beta)]
```

## Control-plane correction

The first source-only package, `work/preoutcome-v1`, was not authorized for
development. An independent audit found that its proposed validator would
have recomputed SHA-256 over complete donor `delta.npy` containers containing
the sealed 410-target rows.

The corrected v2 package:

- preserved v1 without overwriting it;
- reproduced every candidate, plan, and null array exactly;
- validated truth containers by frozen digest binding and byte size only;
- wrote the irreversible access lock before importing truth loaders;
- loaded only the authorized 613 rows and 4,992 response columns through
  Experiment 132's memory-mapped `load_effect_rows`;
- enforced an internal and external 8 GiB process-tree limit;
- used identical observed and bootstrap selector semantics;
- passed 20 tests, Ruff, formatting, compilation, and two independent
  no-truth audits.

## Development result

The one authorized 613-target run selected identity in four of five outer
folds and on the complete development set.

| beta | fixed macro PDS delta | fixed macro signed-credit delta |
| ---: | ---: | ---: |
| 0.00 | `0.00000000` | `0.00000000` |
| 0.25 | `-0.00005109` | `+0.00986897` |
| 0.50 | `-0.00081566` | `+0.01920282` |
| 0.75 | `-0.00237279` | `+0.02740058` |
| 1.00 | `-0.00457143` | `+0.03406376` |

Nested selected result:

- macro PDS delta: **`-0.0000684160`**;
- macro signed-credit delta: **`+0.0020382723`**;
- PDS bootstrap 95%:
  **`[-0.0004686213, +0.0002572156]`**;
- signed-credit bootstrap 95%:
  **`[0.0, +0.0100185203]`**;
- final beta: **`0.0`**;
- promotion checks: **5/17**;
- decision: **`TERMINAL_NO_GO`**.

The mechanism is clear: suppressing the low-variance tail monotonically
improves signed directional credit but monotonically damages retrieval. The
tail is not gratuitous PDS noise; it contains discriminative target identity.
This exact family must not be retuned, routed, nulled, confirmed, integrated,
built, scored, or submitted.

## Safety and provenance

- process-tree peak RSS: `0.7933 GiB`;
- minimum observed system free memory: `68%`;
- untouched 410-target values: not opened;
- reserved confirmation expression: not opened;
- full truth-container hashes: not recomputed;
- frozen null outcomes: not evaluated;
- counts / `cell-eval2` / challenge build / submission: not run.

Hashes:

```text
run.py                         b524286012aab6004784673c05f053a3481a54cb0ed06b221f65ec4c3c83710b
development authorization     f995dd157be828eacb12e604a6114781525bce3ff9eaacd358afa40d41a7c82f
outcome-access lock           e55eda2da9344f7146db7b3b932a5f513dae986ff1c506f20936d06beb7f40a3
development result            34abc9bb8cea648456fe2e96da9e6d79c243ad7fabd3d45b527183a1d02d78f0
per-target/bootstrap arrays   82d185627b84f814c6a4399ad095f0606a9e32cbd51f768d7d8eaadb45bc17cb
result manifest               222f7e0ef3141fa47a4671259a02216643d583d0a83d6a69d892ab9e2f2a0d25
```
