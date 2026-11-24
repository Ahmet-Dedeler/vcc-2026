# 489 — exact-U magnitude pass (`dn_k250_zs`), before the 09-28 riders

**Do not stop the 09-28 riders.** Nothing in `src/vcc2026/rank_mag.py` or its call site breaks a library, a gene total, the ballast, or a Wilcoxon rank sum. The 2.5× target is the formula 423 already priced, and on the mass-holder geometry it does not change a single count between k = 2 and k = 10.

One sentence in 429 is false. A repaired 1-count is frozen only when some control value sits at or below that cell. When every control sits above it, the pass moves the count onto a higher-library cell, the rank sum stays exact, and the delivered lfc becomes more negative (here −1.58 and −1.00). That cannot put the gene back at lfc −35. It is an nmae-only risk on the 24 A/B/C zero genes from 428, and this VM has no control matrices to say which of them are in that regime.

Wilcoxon in this review: `rp423._r1_tie_kernel` / `rc_emit.wilcoxon_z` (scanpy-identical midrank R1, tie-corrected z, `p = 2·Φ̄(|z|)`, BH via `scipy.stats.false_discovery_control`). Cross-check: `cell_eval2.compute_de(backend="scanpy")`, the installed CPU path (`pdex` and `gpudge` are not installed). On the synthetic block the two p-values matched (median ratio 1, max |log ratio| 0).

Reproduce: `uv run python experiments/489-cloud-rankmag-review/repro_rankmag.py` (42 checks, 0 failed; numbers in `work/repro.json`). No challenge data. Linux CPython; the pass is Python and numpy, so the macOS bash 3.2 / BSD grep traps do not apply to it.

## 1. CONFIRMED — a repaired 1-count is not universally frozen

`rank_mag.py:117-158` (`cell_bounds`), `215-247` (`_zs_step`). 429 RESULTS.md:202-204 says a 1-count cannot drop to 0 without tying the control zeros, so the repair's floor is kept. That is true only when the open interval does not contain both the 1-count and 0.

- Controls include 0, and the CPM gaps are narrower than one count at the carrier library. The repaired gene is unchanged. `lfc` stays −4.588. (`repro`, "repair with control zeros". Columns 1 and 2 of that fixture also change: column 1 is the sub-gate payment, column 2 is a second gated gene the fixture left empty, which repair fills and which has no control zero. Column 0, the gene with control zeros, does not move.)
- Every control is count 4 at library 20_000 (CPM 200). A 1-count at library 8_000 is CPM 125, strictly below every control, so 0 is inside its interval. All four low-library 1-counts move onto the high-library class. `lfc` −3.138 → −4.722 (Δ −1.585 = log2(24000/8000)). Full R1 unchanged, Δ log p = 0.
- The real `zero_repair` at floor −5, control CPM 2_000, places a 1-count in every carrier. The pass then moves every low-library 1-count onto the high class (totals held). `lfc` −4.722 → −5.722 (Δ −1.00). High cells go from 1 to 2; low cells go from 1 to 0.

Score effect: `pds` and the profile half of `mse` unchanged (totals). `fid` / `reach` / `jac` unchanged in this geometry (R1 exact, and here the tie term did not move either). `nmae` moves, because that is the point of the pass. A gene repaired to about −5 can be pushed another 1–1.6 lfc more negative if its repaired mass sits on the low library class and no control value is ≤ that class's 1-count CPM (1e6/lib, about 111 at lib 9_000). 428's 24 zeros are deep knockdowns that emitted nothing (A 6, B 11, C 7), including SOD1, H3F3B and NDUFA13. Whether those controls have a zero or a cell below ~111 CPM was not measured here (SUSPECTED that some do not). H1's repair was a no-op, so 429's stacked nmae (−0.00179) does not include this interaction.

Proposed fix, not for the packages 456 already builds (it changes the module sha the builder pins at `make_submission_pm_cd4.py:746`): after the tie override in `cell_bounds`, do not let a positive count fall to 0. See `patch_freeze_last_count.diff`. Mass-holders sit above the control max, so their `xl` is already ≫ 1 and this line does not touch the geometry 423 priced. Do not rebuild the 09-28 riders to pick it up.

## 2. CONFIRMED — the pass targets 2.5× the lfc it is handed; mass-holders never get there

`rank_mag.py:320-324`. `A_t = (ref+1e-9)·2^(k·l) − 1e-9` with `k = 2.5` and `l` the block's delivered lfc after the emitter, the MDE add, ballast and the zero repair. It is not the pre-shrink design lfc.

`rc_emit_k.py:957-965`: q80 protection writes the delivered lfc back into `tgt` for the strongest block genes, and `_t2` then multiplies every negative block target, protected ones included, by the mass-balance factor. A gene whose delivered lfc was shrunk by 0.6 is handed −0.48 when the design was −0.8; the pass aims at −1.20, which is 1.5× the design. Reproduced as a scalar in `repro` ("target is 2.5x the handed lfc").

What the counts do depends on the geometry:

- Mass-holders (the genes 423 says carry the nmae gain). Controls max out near CPM 611 and include zeros, so empty carriers are tied. Two cells at lib 9_000 and two at lib 22_000 sit above that max. k = 2, k = 2.5 and k = 10 move the same 228 counts. The mass-weighted fraction of the k = 2.5 target is 0.42 on this toy and 0.021 on H1 (423 RESULTS, production bulk). The 2.5× multiplier is not load-bearing there.
- A mild gene inside one wide gap (controls at CPM 0 and 2_000, carriers strictly inside). k binds. Base lfc −0.170. k = 2 delivers −0.342 (target −0.340, overshoot 0.0023 from the half-count round at `rank_mag.py:236`). k = 2.5 delivers −0.424 (target −0.425, short by 0.0012). k = 10 delivers −0.748 against a target of −1.699 (fraction 0.50). Full R1 identical. One gene's tie term moves; max |Δ log p| = 7.0×10⁻⁶.

Score effect of the formula, on the genes that move in production: none beyond what 423 measured, because those genes are capacity-bound below even the k = 2 target. Score effect if a slack gene exists: it is taken to 2.5× its post-shrink lfc, i.e. up to 1.5× the pre-shrink design, which can overshoot a true response. 423's feasible fraction says that is not where the nmae prize is. No change for 09-28.

The half-count round (`floor(rem/eff + 0.5)`, line 236) overshot the k = 2 target by 0.0023 lfc on the slack gene. The sign guard (lines 354-356) only skips a flip through zero. A down move cannot flip a sign, and the guard does not clip an overshoot. Negligible next to the channel noise.

## 3. CONFIRMED — `lib` and `lib+1` do not share a level

`rank_mag.py:208-212`. The comment says the +1 remainder cells share their class. `np.round(cw/cw.max(), 4)` splits them: keys for libraries 9_000, 9_001, 14_000, 22_000 are 1.0, 0.9999, 0.6429, 0.4091. `class_cells` (`rc_emit_k.py:284-318`) puts that ±1 on the first cells of the block.

On a down gene whose only two libraries are 9_000 and 9_001, both above the control max, the pass moves 48 counts (8 → 2 on the low cells, 8 → 14 on the high cells). lfc −0.433040 → −0.433100 (Δ −6.0×10⁻⁵). R1 unchanged, tie term unchanged (each class stays an equal-sized tie group). The big class pairs run first, because pairs are sorted by |Δ CPM|; this pair only spends whatever target remains, which on a capacity-bound gene is all of it. Dust on nmae. Same code 423 priced. No change for 09-28. A later cleanup is to group by the emitter's `cell_class` instead of the rounded ratio.

## 4. CONFIRMED — the asserts test the sum of carrier control-ranks, and they cannot be skipped on a returned block

`rank_mag.py:347-349` checks `sum(ranks_of)` over carriers, which is `2·(# controls below) + (# controls equal)`, not the full midrank R1 and not each cell. Lines 368-375 then require libraries, per-gene totals, ballast rows and non-negativity. A failed assert raises; the build does not emit.

Paths that `continue` (nothing moved, line 345; sign guard, line 356) do so before the write. `_subgate_zs` raises if the pool is short or leftover (`repro`: "subgate shortfall raises"). There is no return that skips the four final asserts.

On a block that moves 3 genes: per-cell `ranks_of` identical, full R1 identical (the 423 kernel), ballast identical, libraries identical, totals identical, no sign flip, k = 1 bit-identical to the input, two `apply` calls identical, and the output matches `rp423.postpass(..., "dn_k250_zs")`. `cell_eval2` scanpy on a second draw: 0 call flips, 0 sign flips, max |Δ lfc| 0.595, max |Δ log p| 0 (that draw's tie term did not move).

`cell_bounds` on 400 random trials, 8_441 cells, 134 exact ties: 0 cells whose integer interval failed `searchsorted`, 0 loose edges. Tied cells, including zeros tied with a control zero, come back frozen. A one-library-class down gene (lfc −2.49, sub-gate budget available) is left bit-identical. All-zero carriers against controls that contain zero stay bit-identical (lfc −38.2).

The tie term is not asserted. It moved on the slack gene by 7×10⁻⁶ in log p. 429 measured max |Δ log p| 1.5×10⁻³ and 0 call flips on 255 H1 blocks; those blocks are not on this VM, so that number is not re-measured here. At this size a call flip would be a surprise, and none appeared in the scanpy cross-check.

Proposed strengthening, same sha caveat as finding 1: assert `ranks_of` cell-wise, not only the sum. It would not have fired on any block in this review.

## 5. CONFIRMED negatives

- Evaluator arithmetic. `cell_eval2.de_compute._to_linear` on an in-memory counts AnnData is float64 `count / (lib/1e6)`, then `log1p`. Max abs difference against `rank_mag._normalized` + `np.log1p` is 0. The float32 rowstore path is not this path. 3_070 f64 interval edges recomputed in float32: 0 control-rank crossings. Not an exhaustive search of edges the pass can sit on; 429's H1 f32 check (0 U changes) was not re-run.
- Module sha `ed6728327a15f039d89b05bfc8f65281f138242ea914544116d54ec99ecde638` equals `DE_MAG_MODULE_SHA256` (`make_submission_pm_cd4.py:746`), and the builder checks it before the first block (lines 5575-5577). The import used here is `src/vcc2026/rank_mag.py`. The pass has no RNG.
- `dn` keeps finite `l < 0` with `A > 0` (line 319). `up` is the opposite. The builder's regex (line 749) also accepts `all`, `q50`, `q80`, `q80up`, `q80dn`. 423 closed up-genes (they cost nmae). `chain456.sh` hardcodes `VCC_DE_MAG=dn_k250_zs`. A footgun for a later hand edit, not for these riders.
- Call order (`make_submission_pm_cd4.py:5756-5780`): ballast, then `zero_repair`, then `rank_mag.apply`, and only inside the `bulk_effects` branch, which is on whenever ballast is on. MDE is applied inside the emitter (`rc_emit_k.py`, the `mde_add` / `mde_coarse` keys at lines 5749-5750), so the pass sees the post-MDE lfc. Coarse MDE is what the riders set (`VCC_DE_MDE=0.1,0.1,0.3,coarse`, repair −5 or −2). No invariant above depends on coarse versus plain. The nmae price 429 published was on plain MDE; the coarse stack's board delta is what the slot is for, and it was not re-scored here.
- Jackknife, builder formula with bulk target sum 50_000, 8 seeds of the mass-holder block: relative change +0.354 on every seed, never a decrease. Totals identical, so the profile term in `mse` is identical. 423's H1 read (cap binds, raw mse bit-identical, trace up 7–17%) was not re-run. An increase in the trace can lower raw mse only when the prediction's trace is below the real trace; it cannot raise the profile error.
- Emitter gate `count·(1e6/lib)` versus scorer gate `count/(lib/1e6)`, 2_000 random trials near a handful of counts: no gene crossed CPM 5. 429 found 0 disagreements on A/B/C and H1. A near-threshold case is still possible and was not hunted gene-by-gene.

## What would stop a rider, and why none of this does

A stop would be a move that changes U, a call, a sign, a library, a gene total, or the ballast, or a magnitude change that puts a repaired gene back near lfc −35. None of those happen. The riders stay the right 09-28 slot-2 read: exact-U nmae on top of the coarse MDE and the zero repair, which is the stack 429 did not put on the board.

## Coordinator's disposition (2026-09-27 06:50Z; brought to master from `cloud/489-rankmag-review`)
Accepted; nothing stops the 09-28 zs rider (446bj2z), as the review concludes. Findings 1 (a repaired 1-count can move when
no control sits at or below it: an nmae-only risk on 428's 24 zero genes), 3 (±1-library cells split into their own level)
and 4 (assert ranks_of cell-wise) change the module sha, so they wait for zs's 09-28 board read. If zs is adopted for
October, `patch_freeze_last_count.diff` plus the cell-wise assert go in as one KA-gated module change before the October
rehearsal (unset/zs-off byte-identical; zs-on retraced on H1 with 0 call flips).
