# 331 — local runtime feasible; no biological evidence

September 9, 2026. Protocol/code/lock frozen at **29c141ad**. Terminal
status: **ENGINEERING_FEASIBLE_ONLY**. No biological arrays were read.

The downloaded 114,324,594-byte checkpoint matches the pinned revision and
SHA-256 `0db9cb538f114e79026bf08f45f41ad8dd7ad2de2aaca9a5ca8cd3bd9748ae7a`.
CPU, four threads, one estimator, no AMP/cache/automatic download; inference
ran with Hugging Face offline mode and telemetry disabled.

| synthetic training rows | mean MSE | training-mean null MSE | shuffled-label MSE | predict / repeat seconds |
|---|---:|---:|---:|---:|
| 512 | 0.0526233 | 23.25768 | 22.51695 | 0.1227 / 0.1063 |
| 2,048 | 0.0284078 | 23.50663 | 24.26322 | 0.4494 / 0.4305 |

Both sizes pass output shape/finiteness, repeatability, mean-null,
shuffled-label and runtime gates. Repeated mean predictions are exactly
identical. The very high accuracy on this synthetic Friedman-like problem is
not a biological result or a forecast of performance on perturbation data.
MPS was not tested.

Runner elapsed 14.44 seconds; outer wall 19.55 seconds. Runner macOS maximum
RSS is **1,070,350,336 bytes**; outer `/usr/bin/time -l` recorded
1,070,383,104 bytes, approximately 1 GiB. The five-second guard sampled only
0.31 GiB and missed the brief maximum; that sample is not the peak claim.

The isolated environment contains TabICL 2.2.0, Torch 2.14.0, NumPy 2.5.3,
SciPy 1.18.1, scikit-learn 1.9.0 and huggingface-hub 1.30.0. The root
environment still has NumPy 2.5.2 and no Torch; root `uv.lock` remains
`2e458aacefa51a18d1057e09f4984357ad5340674466047350e7fb12a163d78f`.

Independent post-run verification checked every registered file, checkpoint
bytes/hash, workload-to-terminal JSON identity and MSE/R2 arithmetic. Both
unit tests pass. No terminal probe rerun was needed.

- Input manifest SHA:
  `9ebcbd8c0a336daa93c12d05a1863e07680c1f56bc03a2766248a7c33a2f7ab4`
- Terminal results SHA:
  `87d67ebd67ce51cdb7135c552caba629ac785b0fb91c2940257128e12dff7e8a`

This earns only consideration of a separately registered biological pilot
with a distinct construction and matched controls. No package, production
default, reserve access, paid API/cloud or submission.
