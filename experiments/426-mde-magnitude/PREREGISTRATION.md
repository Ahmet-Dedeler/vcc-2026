# 426a preregistration: MDE-structured DE magnitude on the board (written 2026-09-23 04:30Z, before the entry exists)

**Arm.** 426a = 417e + `VCC_DE_MDE=0.1,0.1,0.3`. That is one variable vs 417e (09-24 slot 1).
- Config diff: exactly `de_mde` + `de_sign_emitter_sha256` (the 426 emitter, whose knob-off path is bit-identical: KA1).
- Build: `run426_arm_build.sh`, builder `43bd737f…`, emitter `72be7e16…`.
- Slot: 09-24 slot 2, via `arm0924_s2b.sh`, provided 417e's DE members sit inside the band vs 417a (the DE channel
  must be pinned for this read).

**Mechanism.**
- nmae is scored only on the reference's significant genes. Given significance, |lfc| ≥ the minimum detectable effect,
  which is ~0.24/sqrt(mean control count per cell).
- Adding sign × min(0.1/sqrt(count), 0.1) to every steered gene's delivered arithmetic-mean LFC gains (2a − 1) × add per
  gene (a = our sign accuracy) while add < |lfc|. Genes below 0.3 counts/cell get nothing.
- The add moves only the A-channel: the K-class split at fixed per-gene totals, then the pass-2 |z| design as before.
  Per-gene totals, and so the bulk, are untouched. Signs are unchanged except for 0.17 % of genes (on H1, 78 of our
  ~221k calls across 85 targets).

**H1 (410b construction, seeds 0–2, paired with 416's base):**

| member | Δ raw |
|---|---:|
| nmae | −0.00960 ± 0.00022 |
| fid | −0.00028 ± 0.00004 |
| jac | −0.00014 |
| reach | −0.0031 ± 0.0022 |
| pds, mse | 0 |

That is **+0.0019 overall** at board slopes. 412's best magnitude rule (global × 1.25) was −0.004 nmae.

**Predictions for the board (426a − 417e, raw):**
- **pds: identical** to ≥ 10 significant figures (the bulk is untouched). Any difference is a bug.
- **mse:** |Δ| < 2e-4 (only the jackknife pair noise on the ballast moves).
- **nmae: −0.010 … −0.022** (H1 × 1 … × 2.3; the board read 410's pure-magnitude weld at 2.3× H1), central −0.013.
- fid, jac: |Δ| ≤ 0.001.
- reach: −0.003 expected, inside board noise.
- **overall: +0.0018 … +0.0055.**

**Reading bins (nmae):**
- **Δnmae ≤ −0.008: the MDE law transfers.**
  - It stays in every future arm.
  - Next is the dose ladder on H1 first (c 0.15 / cap 0.15 and the concordant variant), then the board.
- **−0.008 < Δnmae ≤ −0.003: partial transfer.** Our A/B/C sign accuracy on the gates is lower than H1's 0.56. Keep it
  if overall > 0; no dose increase.
- **Δnmae > −0.003: no transfer.** Close it and record the mechanism: A/B/C gate signs are too close to chance for
  magnitude to pay.
- **pds not identical: stop.** Diagnose before any further DE arm.

## Addendum (2026-09-23 05:05Z, after the build, before the entry exists)
426a built: package `f968fe54…`, config diff exactly `de_mde` + emitter sha, CLIFF identical to 417e's (0.9477 /
0.9408 / 0.9447; the bulk is untouched), novelty rc 0.
- **The MDE dose binds on A/B/C.** `mde_lambda` median 0.84 (min 0, max 1) across the 900 blocks. On H1 it was 1.0 on
  every block. The prior's class-capacity check (6,000-count margin) caps the split on many A/B/C blocks.
- Other diagnostics against 417e:
  - clip_frac median 0.094 (417e 0.032);
  - prescribed calls missed median 5 (417e 2);
  - unsuppressible calls unchanged (median 0, max 87).
- **Scaled predictions, which replace the ones above:** Δnmae −0.008 … −0.018 (central −0.011), overall +0.0015 …
  +0.0045. The bins stand.
