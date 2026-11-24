# 429 — 423's exact-U down-gene pass (`dn_k250_zs`) stacks on 426's MDE: nmae −0.00179 on top of it on H1, U exact on every block; production patch ready, applied only to copies

Status: **COMPLETE** (2026-09-23 05:50–06:58Z; worker, Opus 5.5). H1 six-member harness, plus a capacity census on the
real A/B/C blocks of the 428a package. No git, no `vcc`, no `cell-eval2`, no edits to `scripts/` or `src/`.

## Answer

- **zs stacks on the MDE, sub-additively.** On H1 (3 seeds, production bulk), zs on top of the MDE gives
  **nmae −0.00179 ± 0.00007**, i.e. **+0.00049 overall**. The other five members stay bit-identical.
  - That is 0.77 of zs's stand-alone gain (−0.00232).
  - The loss is overlap, not budget: the MDE has already raised the holder genes' |l|/|r| from 0.39 to 0.56.
- **The stacked arm, MDE + repair + zs − base: nmae −0.01139 ± 0.00020, +0.00234 ± 0.00039 overall.** The MDE alone
  is +0.00185 ± 0.00039.
  - The zero repair is an exact no-op on H1: 0 zero genes on H1 carriers, bit-identical with and without zs.
- **zs keeps U exactly on the MDE's class splits.** Across 255 MDE blocks:
  - 0 genes change their rank sum (423's verify, scorer arithmetic); 0 call flips, 0 sign flips;
  - per-gene totals, libraries and ballast rows are identical;
  - the float32 scorer path gives 0 U changes on 563 moved genes.
- **On A/B/C, zs has 0.70× H1's capacity.** Measured on all 900 428a blocks through the production call: A 0.91×,
  **B 0.14×**, C 1.04×. HeLa has few mass-holder genes. Every block passes the production asserts.
- **Predicted board effect of zs on top of the MDE: nmae −0.0013 … −0.0029 raw, +0.0003 … +0.0008 overall.** No
  other member moves. This is H1 × 0.70 capacity × the 1 … 2.3 pure-magnitude window.
- **Recommendation: ride, don't slot.**
  - Keep zs out of 428a, so the repair's read stays clean.
  - If 426a's board nmae transfers (Δnmae ≤ −0.005 raw vs 417e), zs rides every later MDE arm from the first arm
    after 428a.
  - If 426a's nmae does not transfer, zs waits (§6).
- **Production patch prepared, not applied** (shas in §5). The PI applies it with `patch_zs429.py` and verifies with
  `test_zs429.py`: 229/229 PASS on a patched copy of the current production builder.

## 1. Known answers (all reproduced; `ka429.py`, the driver aborts on a MUST failure)

| row | reference | seeds | stats arrays (sums, jk, p, q, lfc) | members max\|d\| |
|---|---|---|---|---|
| base `--bulk ps` | 416 `base_s{S}_bps` | 0, 1, 2 | bit-identical | 1.1e-16 / 0 / 1.1e-16 |
| base + zs (vendored module) | 423 `dn_k250_zs_s{S}_bps` (423 §7) | 0, 1, 2 | bit-identical | 1.1e-16 |
| MDE c 0.1 / cap 0.1 / min 0.3 | 426 `base_s{S}_bps_mde0.1cap0.1min0.3` | 0, 1, 2 | bit-identical | ≤ 1.1e-16 |
| MDE + repair −5 (SHOULD = MDE) | 429 MDE | 0, 1, 2 | bit-identical | ≤ 2.2e-16 |
| MDE + repair + zs (SHOULD = MDE + zs) | 429 MDE + zs | 0, 1, 2 | bit-identical | ≤ 1.1e-16 |

- The 1e-16 differences are the scorer's polars aggregation, as in 426's KA0.
- In every zs run, every block was also passed through 423's `rp423.postpass`, and the vendored module's output had to
  be bit-identical to it (`run429.py` asserts per block). It was, on all 765 blocks (3 arms × 3 seeds × 85 targets).
- The zero repair touched 0 genes on H1's MDE carriers (`zero_repaired` sum 0).

## 2. The table (H1, `--bulk ps`, seeds 0/1/2, Δ vs the same-seed base; overall = Σ slope·Δ/6 at pds 2.2138, mse −1.051, nmae −1.6532, fid 3.3645, reach 1.1292, jac 2.7263)

Absolute (seed mean):

| arm | pds | mse | nmae | fid | reach | jac |
|---|---:|---:|---:|---:|---:|---:|
| base | 0.89486 | 0.92529 | 0.95232 | 0.50920 | 0.27382 | 0.07480 |
| base + zs | 0.89486 | 0.92529 | 0.95000 | 0.50920 | 0.27382 | 0.07480 |
| MDE | 0.89486 | 0.92529 | 0.94272 | 0.50892 | 0.27077 | 0.07465 |
| MDE + zs | 0.89486 | 0.92529 | 0.94093 | 0.50892 | 0.27077 | 0.07465 |
| MDE + repair | 0.89486 | 0.92529 | 0.94272 | 0.50892 | 0.27077 | 0.07465 |
| **MDE + repair + zs** | 0.89486 | 0.92529 | **0.94093** | 0.50892 | 0.27077 | 0.07465 |

Deltas (mean ± sd over 3 seeds):

| contrast | Δpds | Δmse | Δnmae | Δfid | Δreach | Δjac | overall |
|---|---:|---:|---:|---:|---:|---:|---:|
| base + zs − base (zs alone) | 0 | 0 | −0.00232 ± 0.00005 | 0 | 0 | 0 | **+0.00064 ± 0.00001** |
| MDE − base | 0 | 0 | −0.00960 ± 0.00022 | −0.00028 ± 0.00004 | −0.00305 ± 0.00220 | −0.00014 ± 0.00001 | **+0.00185 ± 0.00039** |
| MDE + zs − base | 0 | 0 | −0.01139 ± 0.00020 | −0.00028 ± 0.00004 | −0.00305 ± 0.00220 | −0.00014 ± 0.00001 | **+0.00234 ± 0.00039** |
| **MDE + repair + zs − base (the stacked arm)** | 0 | 0 | **−0.01139 ± 0.00020** | −0.00028 ± 0.00004 | −0.00305 ± 0.00220 | −0.00014 ± 0.00001 | **+0.00234 ± 0.00039** |
| **MDE + zs − MDE (zs on top of the MDE)** | 0 | 0 | **−0.00179 ± 0.00007** | 0 | 0 | 0 | **+0.00049 ± 0.00002** |
| MDE + repair − MDE | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

- "0" means bit-identical members, not rounding.
- The sd on the stacked arm's overall is the MDE's reach (426: reach is within noise, discrete small-n_real jumps).
- The zs increment's sd is 0.00002 overall.
- **Interaction:** (MDE + zs) − MDE − zs = +0.00053 nmae (per seed +0.00055 / +0.00055 / +0.00050). zs keeps 0.77 of
  its stand-alone gain on top of the MDE.

**Target bootstrap** (58 scored targets, seed-averaged per-target Δnmae, 20,000 paired resamples):

| contrast | mean | 95% CI | p(≥0) | better / worse | worse by > 0.005 | top-3 share |
|---|---:|---|---:|---|---:|---:|
| zs alone | −0.00232 | [−0.00433, −0.00073] | 0.0006 | 32 / 21 | 0 | 0.60 |
| MDE | −0.00960 | [−0.01350, −0.00577] | 0.0000 | 45 / 12 | 9 | 0.23 |
| MDE + repair + zs | −0.01139 | [−0.01590, −0.00711] | 0.0000 | 46 / 12 | 8 | 0.25 |
| **zs on top of the MDE** | **−0.00179** | **[−0.00338, −0.00052]** | 0.0006 | 26 / 25 | **0** | 0.61 |

On MDE carriers the zs gain is less spread than alone: 26 targets better and 25 worse, with the top 3 targets carrying
61%. No target loses more than 0.005 in any contrast that includes zs alone or on top of the MDE; the 8–9 larger per-target
losses in the MDE rows come from the MDE. The CI stays below zero.

## 3. Does zs keep U exactly on the MDE's class splits? Yes

423's `verify` ran per target block in every zs run (`diag-*_dnk250zs.json`, 255 blocks per arm):

| check | base + zs | MDE + zs | MDE + repair + zs |
|---|---|---|---|
| genes with a changed rank sum (touched / untouched) | 0 / 0 | 0 / 0 | 0 / 0 |
| call flips / lfc sign flips | 0 / 0 | 0 / 0 | 0 / 0 |
| per-gene totals, libraries, ballast rows identical | all blocks | all blocks | all blocks |
| vendored module == `rp423.postpass` | all blocks | all blocks | all blocks |
| genes moved per block (mean) | 123 | 99 | 99 |
| genes whose tie term moved (block-internal ties) | 83 | 56 | 56 |
| max \|Δ log p\| | 1.1e-3 | 1.5e-3 | 1.5e-3 |
| float32 scorer path (`f32check429.py`, 6 seed-0 blocks): genes with U moved | 0 of 718 | 0 of 563 | — |

- U is exact, so calls and signs are unchanged. The Wilcoxon p changes only through the tie term among a block's own
  cells (two carriers at the same value stop tying), by ≤ 1.5e-3 in log. This moved no call and no fid, jac or reach
  value on H1.
- **mse.** The pass raises the block jackknife by a median of 6.9% (base) and 5.7% (MDE), max 21%. The mse numerator
  is still bit-identical in every seed: the scorer's correction is min(pred trace, k·real trace), and on H1 the cap
  binds for every target. On the board a larger capped correction can only lower raw mse or leave it (423 §1).

**Why only 0.77 on top of the MDE** (`capacity429.py` over 255 blocks; `budget429.py` on 6 blocks):

| | base carriers | MDE carriers |
|---|---:|---:|
| genes moved / block | 123 | 98 |
| Σ\|Δl\| over moved genes / block | 10.9 | 7.9 |
| moved genes in the reference DE set / block | 16.4 | 13.0 |
| their delivered sign accuracy | 0.708 | 0.719 |
| their median \|l\|/\|r\| before → after | 0.39 → 0.49 | 0.56 → 0.67 |
| Δnmae per unit Σ\|Δl\| (58 scored targets) | −2.18e-4 | −2.32e-4 |
| sub-gate budget binding (share of \|dA\| lost vs an unlimited budget) | 8.7% | 11.7% |

- The MDE has already raised these genes' magnitude (|l|/|r| 0.39 → 0.56). Their class splits leave less interval
  headroom, so zs moves 20% fewer genes and 28% less Σ|Δl|.
- Per unit of Σ|Δl|, the nmae gain is the same (−2.2e-4 vs −2.3e-4). So Σ|Δl| is a usable capacity yardstick for
  pricing A/B/C.
- The carriers' sub-gate budget is not the limit: it costs 9–12%, mostly on one target (APOE).

## 4. A/B/C rehearsal: the production call on all 900 blocks of 428a (`abc429.py`)

428a (`c7bf5d47…`) is 426a + `VCC_DE_ZERO_REPAIR=-5`: exactly the carriers zs would see in production, after ballast
and repair.
- The package was restored into `work/abc/`, read block by block, and deleted afterwards.
- Each context's scorer-order control block was built from `data/challenge/context_{A,B,C}.h5ad` with
  `src_rank_mag.scorer_control_block`, as the patched builder does.
- **The scorer-order gate equals the emitter gate on A, B and C** (0 differing genes; smoke run with `--gate`), and on
  H1.
- **All 900 blocks passed `apply()`'s production asserts** (totals, libraries, U per moved gene, ballast, no flip).
  Every block had at least one moved gene. K was read from the ballast libraries and asserted to be the last K rows.

| context | blocks | genes moved / block | pure mass-holders / block | Σ\|Δl\| / block (mean, q10–q90) | per-block median \|l\| of moved genes, before → after (mean over blocks) | capacity vs H1 MDE |
|---|---:|---:|---:|---|---|---:|
| A | 300 | 40.9 | 33.0 | 7.15 (4.00–10.66) | 0.32 → 0.54 | **0.91** |
| B | 300 | 31.3 | 12.9 | 1.09 (0.60–1.60) | 0.39 → 0.43 | **0.14** |
| C | 300 | 56.4 | 42.8 | 8.20 (5.65–11.31) | 0.36 → 0.54 | **1.04** |
| H1 MDE | 255 (6 for holders and \|l\|) | 98.5 | 46.3 | 7.86 | 0.29 → 0.33 | 1 |

- **B (HeLa) has a third of A's and C's pure holders.** Its moved genes barely move (0.39 → 0.43).
- A and C move fewer genes than H1 (41 and 56 against 98) but move each one much further (0.32 → 0.54, against H1's
  0.29 → 0.33). Their summed capacity is H1's.
- It is not the library-class spread: B's (1/q25 − 1/q75)·q50 is 0.46, the same as H1's; A's is 0.70, C's 0.64.
- **Mean capacity over the three contexts: 0.70× H1's.** This is a capacity measure only; A/B/C sign accuracy on these
  genes cannot be measured without truth.

## 5. Production patch (prepared; applied only to copies)

| file | sha256 | what |
|---|---|---|
| `src_rank_mag.py` | `ed6728327a15f039d89b05bfc8f65281f138242ea914544116d54ec99ecde638` | the module; the PI copies it verbatim to `src/vcc2026/rank_mag.py` |
| `patch_zs429.py` | `b82f742f76fa7f8779f85923ab4ef3cbb5d0bf8f9b5edf556c9003166c3cd602` | builder patch plus module install; refuses if `src_rank_mag.py` ≠ the tested sha |
| `test_zs429.py` | `8845b47596c55e29078466e352061eca28584ab575b1146d2b4fc4576791c7f3` | tests (dry-copy mode via `TEST_MODULE` / `TEST_BUILDER` / `TEST_BUILDER_ORIG`) |
| `postcheck429.py` | `e8c40472faf88fb785ef46975bc17603f9322a4a0821e7ea1058340bf5a12b38` | post-build check of a zs arm against its pass-off parent (423 §8 item 4) |

**`src_rank_mag.py`**: 423's zs path, vendored with a module docstring.
- Vendored verbatim: `cell_bounds`, `ranks_of`, `fill`, `select`, `_levels`, `_zs_step`, `_take`.
- `postpass` is zs-only. The ballast legs of `_subgate` / `settle_ballast` became `_subgate_zs`, the pool-only
  bookkeeping zs actually uses; it asserts the pool empties exactly.
- New: `scorer_control_block` (323's `exact_de323.control_block`, streamed with the same chunking), `block_lfc` (the
  evaluator's lfc), `parse_spec` (`<sel>_k<k*100>_zs`, zs only) and `apply()`, the builder entry point.
  `apply()` returns flat float diagnostics prefixed `zs_`, because the builder's `rc_diagnostics_summary` calls
  `float()` on every value; rp423's `None` fields would crash it.
- Extra asserts beyond rp423: ballast rows untouched, carriers' libraries > 0, and a K = 0 block works (no ballast
  needed).

**`patch_zs429.py`** makes five anchored edits (each anchor must occur exactly once), and it refuses:
- if the 428 patch is missing;
- if the builder already carries the 429 patch;
- if the module file exists.

Nothing is written until every anchor has matched. The edits:
1. The knob `VCC_DE_MAG=<sel>_k<k*100>_zs` (only `dn_k250_zs` is priced), placed after the zero repair's knob. It
   requires `VCC_DE_SIGN_PRIOR`, and it pins the module sha: `DE_MAG_MODULE_SHA256 = ed672832…`.
2. A guard that requires `VCC_PERT_SCALE_BALLAST`.
3. 423 §8's scorer-order control block (count / (lib / 1e6)), one per context.
   - It is built right after `rc_control`, while the count matrix is alive, and released with `rc_control`.
   - The module's sha is checked here, before any block is emitted: missing or wrong → `SystemExit`.
   - It costs +1.6 GB per context. 426a's build log shows the context stage at 17.7 GiB, set earlier in the run;
     `--max-rss-gib` is 32.
4. The pass, in the block loop after `ballast_rows` and the zero repair, before `to_csr`:
   `_rank_mag.apply(block, cells_per_pert − ballast_k, rc_scorer_control, DE_MAG)`. Its diagnostics join `rc_diag`.
5. The config key `de_mag` (`variant`, `module_sha256`, `where`), written only when the knob is set.

Unset, no code runs, no control block is built and no key is written. The pass consumes no RNG, so a knob-on build
differs from its knob-off parent only in moved carrier counts.

**Run order after ballast: zero repair, then zs.**
- (a) The repair gets first claim on each carrier's sub-gate counts. It pays one count per repaired cell from the cell's
  largest sub-gate gene. zs then budgets from what is left and cannot overdraw it. The reverse order could leave a
  carrier with nothing to pay with: zs exhausts ~20–28 carriers' sub-gate budget per H1 block.
- (b) The pass's asserts (per-gene totals, libraries, every moved gene's U, ballast rows, no flip) then certify the final
  delivered carriers; nothing edits the block after them.
- (c) zs moves only genes with a finite, nonzero delivered lfc, so it sees repaired genes as what they are after the
  repair. Their single counts are effectively frozen: a 1-count cell cannot drop to 0 without tying the control zeros.
  The repair's floor is therefore kept.

**Applied to copies only** (`work/dry/`):

| builder copy | before | after `patch_zs429.py` | module copy |
|---|---|---|---|
| `builder_head.py` = HEAD `a575cdca` = the current production builder | `986e8fad8e4a…` | `9f562de3c4a0818165deabc7cd844fde6cc8eb4eaa2a3746403bf2ecf51b78f4` | `rank_mag_head.py` `ed672832…` |
| `builder_wt.py` (the working tree at 06:00Z, carrying an uncommitted D/E/F-rehearsal edit since reverted) | `c0e3a923…` | `ccb9c040343b9e0479abe4cdd5d2e08eb973067a6d614dc42b2973b3c6a2a35e` | `rank_mag_wt.py` `ed672832…` |

Re-running the patch on a patched copy refuses ("already carries the 429 patch"), and it refuses when the module file
exists.

**`test_zs429.py` results:** 229/229 PASS on the HEAD copy, full mode, peak 7.7 GB (`work/test_head_full_final.log`);
224/224 on the working-tree copy, `--no-control` (`work/test_wt_nocontrol_final.log`).
- **Module invariants on 14 real H1 blocks** (6 base and 6 MDE production-bulk dumps with K 2–7, plus 423's two cap15
  K = 1 blocks):
  - k = 1 is a bit-identical no-op;
  - dn_k250_zs keeps per-gene totals, libraries, ballast rows, non-negativity, every gene's R1 and the calls;
  - no sign flips; untouched genes' tie terms are exact; all moved genes are down and none rises;
  - a K = 0 block works.
- **Known answer vs 423:** the module equals `rp423.postpass` bit-for-bit, block and shared diagnostics, on all 14
  blocks. `block_lfc` equals the evaluator's lfc bit-for-bit, `apply()` equals `postpass(block_lfc)`, and the
  diagnostics are flat floats.
- **Grammar:** 4 accepted, 8 refused.
- **Scorer-order control block:** rebuilt from the h1-test controls, it equals 323's `reference.npz` `keep`/`cs`/`ref_lin`
  exactly (sha256 `5c26cdf3…`, `cc3155b8…`, `8940ae3b…`). It equals the emitter gate, over 18,400 controls.
- **Builder** (dry copy at `scripts/.dry429_<pid>.py`, deleted afterwards; the test verifies the deletion):
  - unset → `DE_MAG None`; `dn_k250_zs` parses alongside `VCC_DE_MDE` and `VCC_DE_ZERO_REPAIR`;
  - 6 refusals: no ballast, no sign prior, `_pc` mode, no mode, k < 1, garbage;
  - the pinned sha equals the module under test;
  - run order: zero repair < zs < `to_csr`; the control block is built before `del control_matrix`;
  - **additions only:** deleting the patch's snippets from the patched copy restores the unpatched builder
    byte-for-byte. Every added statement is one of: a module-level knob definition (`DE_MAG_ENV`, `DE_MAG = None`,
    the module path and sha constants), code under `DE_MAG is not None`, a config entry that is `{}` when unset, or one
    of two `rc_scorer_control = None` no-ops.
  - This is the static proof that the knob off is bit-identical. The dynamic proof is the build itself: 429a's config
    diff vs 428a must be exactly `de_mag`, and `postcheck429.py` must pass.

**For the PI (in order):**
1. `uv run python experiments/429-zs-production/patch_zs429.py`, in place: builder `986e8fad…` → `9f562de3…`,
   module `ed672832…`.
2. `uv run python experiments/429-zs-production/test_zs429.py`. With the production defaults it runs in place;
   `TEST_BUILDER_ORIG` is optional.
3. Commit with those shas.
4. Build 429a = 428a's exact config + `VCC_DE_MAG=dn_k250_zs`, same seed.
5. Check that the config diff is exactly `de_mag`.
6. Restore both packages and run `postcheck429.py 429a.h5ad 428a.h5ad`. It checks that every block keeps its totals,
   libraries, ballast rows and U, and that it moves only down genes. The self-test on 18 of 428a's blocks passes, and
   injected corruptions are caught.

## 6. Priced recommendation

**The zs increment on the board** (on MDE carriers; every other member unchanged, mse ≤ 0):

| reasoning | Δnmae raw | overall |
|---|---:|---:|
| 423 §6 window as stated (H1 × 1 … 2.3) | −0.0018 … −0.0041 | +0.0005 … +0.0011 |
| **capacity-corrected (H1 × 0.70 × 1 … 2.3)** | **−0.0013 … −0.0029** | **+0.0003 … +0.0008** |

- The ×1 … ×2.3 window is a board/H1 ratio for global pure-magnitude changes; its only board reading is 410's weld,
  H1 +0.024 → board +0.054, a harm. zs acts only on down mass-holders, whose prevalence differs by context, so the A/B/C
  capacity factor (0.70; B 0.14) applies on top of it.
- **The 426a correction has not landed** (09-24 slot 2; `docs/BOARD-LEDGER.md` has no 426a row as of 07:00Z 09-23).
  When it lands, replace the window with k = Δnmae(426a − 417e) / (−0.0096). Predicted zs increment ≈ −0.00125 × k raw.
  - This slightly under-reads: the MDE delivers at `mde_lambda` median 0.84 on A/B/C against 1.0 on H1, so k carries
    that shortfall and zs does not.
  - zs's genes carry higher H1 sign accuracy (0.72) than the MDE's average gene, so they should transfer at least as
    well as the MDE, unless A/B/C's holder signs are worse.

**The stacked arm (MDE + repair + zs) vs 417e:**
- MDE: +0.0018 … +0.0055 (426's preregistration).
- repair: 0 … +0.0022 (428's estimate: nmae 0 … −0.008 raw, unmeasurable on H1).
- zs: +0.0003 … +0.0008.
- **Total: +0.0021 … +0.0085 overall.** H1 prices the stack at +0.00234 ± 0.00039; pds is bit-identical and mse is
  identical or lower.

**Ride or wait?** zs should ride with the MDE, conditional on 426a's read. It should not go into 428a.

- **426a's board Δnmae ≤ −0.005 raw (k ≥ 0.5): ride.**
  - Build 429a = 428a + `VCC_DE_MAG=dn_k250_zs` (§5). It is one variable vs 428a and deterministic: same seed, and the
    pass uses no RNG. So its nmae delta vs 428a is the pass alone, with no realisation noise.
  - Carry zs on every later DE-stacked arm, starting with the first arm after 428a.
  - If the 09-25 partner is a bulk step on 428a's carriers, zs can ride it. Its nmae delta vs 428a reads zs: bulk steps
    move nmae by ~1e-4 raw (417a vs 410b), against zs's predicted −0.0013 … −0.0029. pds stays bit-identical.
  - The caveat is mse, which zs can only lower through the capped jackknife (0 on H1). If that partner's mse read has to
    be clean at the 1e-4 level, read zs in 429a instead.
- **−0.005 < Δnmae ≤ −0.003 (0.3 ≤ k < 0.5): wait for a spare slot.** Send 429a alone when a slot has nothing better.
  Expected +0.0001 … +0.0003: a calibration read of the holder class, not a gain.
- **Δnmae > −0.003 or wrong-signed (k < 0.3): wait.** H1 then mis-prices A-channel magnitude, and zs is the same class.
  Keep the patch unapplied or knob-off.
- **Do not add zs to 428a.** Both are nmae-only changes, and 428a's single-variable read of the zero-gene repair (0 …
  −0.008) would be confounded.
- **As a slot of its own, zs is small** (+0.0003 … +0.0008 overall). It is worth a slot only as a calibration read, or
  as a free rider.

## 7. Negatives and notes

- **zs is sub-additive with the MDE**: the interaction is +0.00053, so it keeps 0.77. The two passes act on the same
  genes' arithmetic means from different ends: the MDE through the class split, zs within the control intervals. More
  A-channel magnitude on these genes does not add linearly.
- **The zs gain on MDE carriers is lumpier than alone**: 26/25 targets better/worse, top-3 share 0.61, against 32/21 and
  0.60 alone. There are still 0 targets worse by more than 0.005, and the CI is [−0.0034, −0.0005].
- **B (HeLa) is a weak context for zs**: 0.14× capacity, 13 holders per block. A D/E/F line may land anywhere in the
  0.14–1.04 range; the pass needs only the context's own controls, so it transfers mechanically.

## 8. Compliance

- **H1 runs:** 9 runs, strictly serial (`drive429.sh`, queue `work/queue_a.txt`).
  - Each started behind the memory gate; every check passed at 73–77% free and 27–32 GiB free+inactive.
  - Peak RSS 7.00–7.59 GB per run (`/usr/bin/time -l` in `work/logs/`), under the 8 GB cap.
  - 428a's production build was running when the queue started and finished at ~06:05Z. It had peaked at 19.7 GiB
    before the queue started.
- **Other jobs** (not H1 runs; several ran alongside an H1 run, with system memory at 73–77% free throughout):

  | job | peak RSS |
  |---|---:|
  | `test_zs429.py` (earlier versions) | 4.4 / 9.8 / **11.6** / 9.3 GB |
  | `test_zs429.py` (final version) | 7.7 GB |
  | `abc429.py` smoke / full | 9.4 / 5.6 GB |
  | `budget429.py` | 3.7 GB |
  | `postcheck429.py --selftest` | 4.8 GB |
  | `f32check429.py` (run alone after the queue) | 9.3 GB |

  - **Deviation:** one test run hit 11.6 GB, above the ~10 GB guidance. It read the whole h1-test matrix and copied
    the 1.6 GB control block to hash it. The final test streams only the control rows and hashes in place, and peaks
    at 7.7 GB.
- **Files written:**
  - Everything is inside `experiments/429-zs-production/`, with two exceptions:
    - the mandated transient dry copy `scripts/.dry429_<pid>.py`, deleted by the test (verified absent);
    - the 12 GB restored 428a h5ad in `work/abc/`, deleted after the census.
  - No numba or bytecode cache landed in other experiments' directories (checked with `find -newer`).
  - `NUMBA_CACHE_DIR` is `work/numba-cache`.
- **Not done:** no git, no `vcc`, no `cell-eval2`, no submission package, nothing in `submissions/` (only read: the
  428a `.vcc` and its meta). No edits to `scripts/`, `src/` or any other experiment's files. The patch was applied
  only to copies in `work/dry/`.

## Files

| file | what |
|---|---|
| `run429.py`, `k_emit429.py` (= `k_emit426.py`, sha `d52199a1…`), `drive429.sh`, `ka429.py` | harness (426's plus the zs pass after ballast and repair; vendored vs rp423 asserted per block), serial driver, known answers |
| `table429.py` → `work/table429.{md,json}` | the six-member table, bootstrap, verify census |
| `capacity429.py` → `work/capacity429.json`; `budget429.py` → `work/budget429.log` | what zs moves on base vs MDE carriers; sub-gate budget |
| `f32check429.py` → `work/f32check429.log` | float32 scorer-path U check |
| `abc429.py` → `work/abc/abc429.json`, `abc429_full.log` | the A/B/C production-call census on 428a |
| `src_rank_mag.py`, `patch_zs429.py`, `test_zs429.py`, `postcheck429.py` | the production deliverables (§5) |
| `work/dry/` | builder and module copies before/after the patch |
| `work/h1/` | `{metrics,stats,per-target,diag}-<tag>.*` for the 15 arm×seed rows |
| `work/blocks/` | 12 dumped seed-0 blocks (base and MDE carriers, first 6 targets) used by the tests |
