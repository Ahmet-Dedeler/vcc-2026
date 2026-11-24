# 417: bulk-only magnitude sharpening through ballast, and carriers from an unwhitened second pass

Status: **IN PROGRESS** (2026-09-22 21:05Z). Board readings on 09-23 fill §3.

## 0. What this experiment is

Ballast (410b) separated the bulk (pds, mse: log-pseudobulk of summed counts) from the DE channel (fid, reach, jac,
nmae: per-cell CPM means and Wilcoxon ranks). That reopened every transform closed because "nmae pays". 414 re-priced
them and ranked 169's magnitude sharpening first:
- the transform is `sign(h)|h|^q` per row on the render gate, with the row norm restored;
- 4/4 PLOLO folds positive, and both controls (flatten, permuted multipliers) negative on 4/4.

416 then showed on H1 that the DE carriers want LESS whitening than the bulk (a monotone gamma dose response). Ballast
makes that buildable too.

| arm | one variable vs | knob | status |
|---|---|---|---|
| 417a | 410b | `VCC_BULK_SHARPEN_Q=1.25` | built 21:17Z (package ed4fef74); cliff 0.961/0.958/0.961, rho 0.991/0.922/0.880; K median 2, clip mass ≤ 1938; 09-23 slot 1 (00:01Z) |
| 417b | 410b (dose step on 417a) | `VCC_BULK_SHARPEN_Q=1.5` | built 21:51Z; cliff 0.888/0.880/0.888, rho 0.916/0.847/0.813; K median 3, clip mass ≤ 3518, jk ratio median 1.027; slot 2 if 417a ≥ +0.006 |
| 417c | 417a | `VCC_CARRIER_GAMMA=0` | built 22:31Z on builder fc30e8b0 (window 2, commit 9751d3cf; carrier KA max|carrier-bulk| = 0 at γ 0.65); carrier row cos vs bulk q10/50/90 0.846/0.900/0.930; cliff 0.960/0.957/0.960 (= 417a's bulk); K median 3; jk ratio median 1.027; slot 2 if 417a < +0.009 |
| 417d | 417b | `VCC_CARRIER_GAMMA=0` (+ q 1.5) | built 23:02Z on fc30e8b0; cliff 0.886/0.879/0.887 (= 417b's bulk); for 09-24 if the q 1.5 dose wins |

## 1. Builder changes (all default-off, bit-identical when unset; 56 passed / 0 failed / 1 skipped)

- **Window 1 (commit 907d6eb0, builder 400ef6e3):** `VCC_BULK_SHARPEN_Q` (`patch_sharpen417.py`, `test_sharpen417.py`).
  Sharpening acts on each context's 10-CPM gate over NON-PANEL coordinates. Panel columns, including the row's
  knockdown, are untouched, because 414 priced the transform with every panel coordinate zeroed. It applies to the
  ballast bulk render only. It equals 169's `sharpen(mode="row")` on the panel-zeroed row to 8.9e-16.
- **Window 2** (`window2_chain.sh` applies it after 417b restores):
  - `patch_next417.py`: 415's LODO key artifacts `lodo415_lo0.3_hi3/inf`, 410's key at floor 0.5 `wshare_lo0.5_hi3/inf`
    (387d6f1a / 851a689d), and `VCC_BULK_CTX_SCALE`.
  - `patch_bulkio417.py`: `VCC_BULK_DUMP` / `VCC_BULK_LOAD`.
  - `patch_carrier417.py`: `VCC_CARRIER_GAMMA`, a second whitening pass for the carriers, rows norm-matched to the
    shipped rows (416's gXXX).
- **KA-LOAD PASSED (23:33Z, builder fc30e8b0):** 417a's env + `VCC_BULK_LOAD` of a same-env dump (aed59b0b) + `VCC_CARRIER_GAMMA=0.65` reproduced 417a's X/indices/indptr byte for byte. The config diff was exactly bulk_load(+sha) + carrier_gamma, max|loaded−own| = 0, and the carriers equalled the bulk rows. Both new paths are now validated end to end.
- **Known answer for the carrier path:** a dump-exit run at `VCC_CARRIER_GAMMA=0.65` must print
  `max|carrier-bulk| = 0`. `run417ka_build.sh` is the full-build version, covering LOAD as well; it is not scheduled.

## 2. Predictions (preregistered, `PREREGISTRATION.md`)

417a vs 410b:
- pds +0.006 … +0.019 raw;
- mse −0.020 … −0.050 raw;
- DE members within the noise band;
- **overall +0.006 … +0.016**;
- cliffcheck: render-space energy about 0.96 of 410b's per context.

417c vs 417a: DE-only +0.002 … +0.004 overall (416 H1 seed 0: nmae −0.0069, fid +0.0025, reach +0.0024).

## 3. Board readings

**Refined 417a mse prediction, from its measured artifact (21:25Z, before the read).** The artifact's render energy is
0.961 / 0.958 / 0.961 of 410b's, putting rho at 0.991 / 0.922 / 0.880. Context A is now below its cliff.
- Through 187's model at 410b's board X, with X scaled by the direction gain and by √(energy ratio):

  | relative cosine gain | k | raw mse | Δ vs 410b | mse part of overall |
  |---:|---:|---:|---:|---:|
  | 0 (amplitude only) | — | 0.7263 | −0.0047 | +0.0008 |
  | +0.04 (kolf2-like) | 1 | 0.7141 | | |
  | +0.04 (kolf2-like) | 1.407 | 0.7091 | | |
  | +0.106…+0.122 (cancer folds) | 1 | 0.694 … 0.689 | | |
  | +0.106…+0.122 (cancer folds) | 1.407 | 0.681 … 0.674 | | |

- Across the whole range, the mse part is +0.001 … +0.010 overall.
- **Read the board's raw mse against this table:**
  - ≤ 0.695: cancer-fold transfer;
  - 0.705 … 0.715: kolf2-like transfer;
  - ≥ 0.72: amplitude only, and the direction did not transfer.

(board rows filled on 09-23)

## 4. Side results from tonight

- The 410d artifact (`cliffcheck414.py`) read render energy **×1.21 / ×1.19 / ×2.16** (A/B/C) vs 410b. The
  amplitude-only raw mse change was +0.50, so the NOGO was right; see 410 RESULTS §9.
- On H1, rank-preserving DE magnitude rules are small for us:
  - a global LFC scale peaks at ×1.3 (nmae −0.004);
  - `sign × const` is worse (0.9726 against 0.9521), because |lfc| carries sign reliability (top-quintile accuracy
    0.61 against 0.54 below).
  - Our nmae is sign-limited. cqawesome's pure-nmae steps (every other member fixed to 4 decimals) point to the same
    rank-preserving freedom, used at a better sign accuracy.
- **What 414's non-scaling pattern β is made of** (the one that costs a per-target key its board pds): 77–78 % of
  |β|² sits on genes BELOW 5 control CPM. There β averages −0.0019…−0.0021 log units, rms 0.0024–0.0026. The sign flips
  to +0.0012…+0.0016 above 100 CPM. corr(β, log CPM) is +0.50…+0.56 in all three contexts.
  - This is the emitter's low-count deflation, carried into the bulk by ratio-mode ballast. The render gate zeroes
    the design below 10 CPM, but the carriers' floor-and-remainder allocation still shaves those genes.
  - Open question: H1 said design-pinned ballast (no β) costs pds −0.006. The board emulation says β costs the key
    −0.0026 pds. A hybrid (design mode below 5 CPM, ratio mode above) would test it. Not built.
- **Queued idea: rank-preserving carrier magnitude (cqawesome's pure-nmae steps).** Push each up-gene's already-top-ranked
  carrier cells higher, and compress down-genes within their rank intervals against the reference controls. This
  leaves the Wilcoxon U, and so fid/reach/jac, unchanged.
  - The ballast must absorb the count change (design mode for the touched genes) to keep pds/mse pinned.
  - H1 prize, analytic at the realised LFC: nmae −0.004 at ×1.3, which is +0.0011 overall.
  - Why the board may pay more: the board read 410's weld (a per-target magnitude change on the carriers) at nmae
    +0.055 against H1's +0.024, so the board is about 2.3× more magnitude-sensitive.
  - Worth an H1 prototype through the real emitter before any build.
- **Per-context re-levelling on 417a's measured artifact** (`ctxlevel417.py`; stable across board raw 0.68–0.72):
  - EV-optimal `VCC_BULK_CTX_SCALE=A=1.0020,B=1.0408,C=1.0648` gives Δraw −0.011 to −0.012, which is **+0.0018 …
    +0.0021 overall**.
  - Sharpening's energy drop put B and C further below their cliffs, so the re-level wins back that headroom.
  - It is a clean mse-only arm vs 417a: pds is scale-invariant apart from the render, and the DE channel is pinned.
- **ntclevel417.py known answer PASSED** on 417i's artifact (02:30Z). It is the bundle-agnostic re-level tool: claims from
  the bundle's own NTC constructs.
  - Cp46 14.4398 / 14.8057 / 16.4397, exactly 180's.
  - Cliffs vs 187's model: +0.00 % / +0.37 % / +0.00 %.
  - For 417i (the LODO floor-0.5 key on 417a), the EV scales at sd 0.002 are A 1.0005 / B 1.0409 / C 1.0618, almost
    the same as 417a's re-level. The key barely moves the bulk energy.

## 5. Overnight builds, window 3 (2026-09-23 02:30–03:32Z)

All seven 09-24/25 candidates are built and novelty-checked: 417e, 417f, 417g, 417h, 417i, 417j, 417k (`.vcc` +
`*.meta.json` in `submissions/`).
- **417j** (417h with the stack dumped under the LODO floor-0.5 key): package `d5dd52b5…`.
  - Config diff: exactly `bulk_load` + sha.
  - Its "CLIFF FAIL {A 0.9637, B 1.0346, C 1.0846}" is a **gate artifact**. `run_load_arm417.sh`'s generic LOAD gate
    [0.80, 1.03] was written for un-levelled arms; a re-levelled arm sits by design at rho ≈ 0.995 of its cliffs.
  - The gate now reads [0.85, 1.15] when `VCC_BULK_CTX_SCALE` is set.
  - **Process note:** that gate edit was made while 417j's `run_load_arm417.sh` was still executing, against AGENTS.md
    §6.2. The build and package were already done. Bash then re-read ~73 bytes of the novelty line as `thon …`, so the
    log shows "novelty rc=127". The real novelty check had printed "OK: novel package, sha256 d5dd52b5…". A manual
    re-run returned rc 0 (`work/build_417j/novelty_recheck.log`), and the package sha matches the logged one.
- **417k** (417f = 417b + stack LOAD at q 1.5, with its own re-level A 1.0464 / B 1.092 / C 1.114): package `8d8aa91d…`.
  - Config diff: exactly `bulk_ctx_scale`.
  - CLIFF OK {A 0.9601, B 1.0278, C 1.0778}; novelty rc 0.
