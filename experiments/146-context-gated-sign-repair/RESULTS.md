# Experiment 146 — Results

Status: **SOURCE-ONLY NO-GO**

Evaluated: **2026-08-31**

Experiment 146 tested whether a strictly cross-fitted, context-aware classifier
could identify a very small set of downstream coordinates whose a15-core sign
should be flipped. The fixed full-feature HGB arm produced a positive average
signed-credit delta, but the effect was too small, too concentrated, and not
distinguishable from the preregistered permutation null or matched controls.

**Do not promote, integrate, build, or submit this method.**

## Primary result

| Quantity | Result |
| --- | ---: |
| Macro signed-credit delta | `+0.0000119081` |
| Experimental-batch bootstrap 95% CI | `[+0.0000006788, +0.0000235618]` |
| Whole-target bootstrap 95% CI | `[+0.0000000924, +0.0000246646]` |
| Macro raw PDS delta | `-0.0000770895` |
| Minimum context raw PDS delta | `-0.0009259259` |
| Positive contexts | `4 / 16` |
| Contexts with weighted flip precision at least 65% | `4 / 16` |
| Permutation-null p95 | `+0.0003156613` |
| Positive-gain concentration | `63.76%` |

The positive average is about 26 times smaller than the permutation-null p95.
It therefore does not establish transferable sign-repair signal.

## Mandatory gate outcome

Passed:

- batch-bootstrap lower bound above zero;
- whole-target-bootstrap lower bound above zero;
- all five held target-fold deltas nonnegative;
- macro and minimum-context PDS noninferiority;
- intervention count, energy, PDS-axis energy, and signed-L1 caps;
- target-fold sizes and metric-mask separation.

Failed:

- at least 12 of 16 contexts improve;
- at least 12 contexts reach 65% weighted flip precision;
- full HGB confidently beats zero-tail;
- full HGB confidently beats donor-only HGB;
- full HGB beats all other matched controls;
- observed gain exceeds permutation-null p95;
- all frozen FDR sensitivities are positive;
- positive-gain concentration is at most 30%.

At FDR `0.10`, macro signed-credit delta was negative
(`-0.0000083846`). Full HGB also trailed donor-only HGB in macro signed
credit (`-0.0000158517` difference).

## Safety and isolation

The completed evaluation:

- opened no challenge control expression;
- opened no development or held 613/410 outcomes;
- did not run `cell-eval2`;
- generated no counts;
- created no challenge build or submission;
- respected the preregistered batch-by-target cross-fitting and
  availability-stratified permutation design;
- peaked at `0.910 GiB` process-tree RSS.

An independent read-only audit reproduced the gate decision and verified the
evaluation artifact hashes.

## Immutable evidence

| Artifact | SHA-256 |
| --- | --- |
| `work/source-evaluation-v1/artifact-manifest.json` | `0604e64f19afce039aad110db51fd907c261e420762835cc6772d3a18e99e60d` |
| `work/source-evaluation-v1/result.json` | `a1d1a59460182874a691769a05b732f64dfc4d93e97769cf21765c0b7f9dc70d` |
| `work/source-evaluation-v1/metrics.npz` | `e2382620e783e3a2c6fc2e3fa15842806b617fbd93471ec7dc554fd422423ef7` |

The evaluated source inputs had these hashes:

| File | SHA-256 |
| --- | --- |
| `run.py` | `198f9041ba7a7f11c01b8c2ae01405271c7aa7a6c888266b8fedd163c5dbb4b8` |
| `test_run.py` | `d7e6836a1ed30913742d49d7b266dcf28ea43aa39a63a5a3ba85f101abb15750` |
| `README.md` | `a56cbf8ae6b3fcb02e2771c6fb16cbe88bcd226d0160729ecbd5fe6950f45914` |
| `preregistered_protocol.json` | `44e1e91d0a13f095a090cacebf4f40d6c13aa97ff4014cdf00a6e6bf4c2b843c` |

`README.md`, the protocol, and evaluated source files should remain unchanged
so these bindings continue to describe the completed run. This result file is
the closure record for the historical scaffold language in `README.md`.
