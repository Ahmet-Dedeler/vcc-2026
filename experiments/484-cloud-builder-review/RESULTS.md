# 484 — Production builder bug hunt

**Decision: one real silent score bug, on the DE channel of every current entry including 446b.** `clean4q80` does not keep the delivered LFC of a protected down gene. The mass-balance factor `kf` multiplies it. On the clean4q80 rehearsal `experiments/432-def-rehearsal-new-knobs/work/b2/rc_audit_b2.jsonl` (900 blocks, contexts D/E/F, the September stack plus `VCC_DE_MAG=dn_k250_zs`), the factor that actually hits those genes has median **0.605**, minimum 0.385, maximum 0.971, and is below 0.99 on **every** block. Protected up genes are left alone. Signs are kept. `pds` and `mse` are pinned by per-gene totals. The member that moves is `nmae`.

Do not patch this into the next slot as a correction of 446b. 395 and 446b were scored with the shrink in them. The patch below is a new arm: price `nmae` locally first. A raw `nmae` move of about 0.017 is the size that showed up on the board for coarse MDE (−0.0168 raw, 433a). Whether this bug is that large depends on how much of the reference |lfc| mass sits in the protected down block; that share was not measured here (no `data/harness`).

Reproduction: `uv run python experiments/484-cloud-builder-review/repro_builder_bugs.py` (log: `repro.out.txt`). Label key: **CONFIRMED** means the production lines were executed or a recorded audit was re-read. **SUSPECTED** means the mechanism is shown and the score size is not.

## 1. CONFIRMED — q80 shrinks protected down genes

`src/vcc2026/rc_emit_k.py:959-965`. Same lines in `experiments/394-protected-block-signs/k_emit394.py:933-938`, so this has been in the scored q80 path since 394. The builder comment at `scripts/make_submission_pm_cd4.py:622` and 394's RESULTS ("the strongest 20% keep their delivered LFC") describe the other behavior.

What the code does, after restoring `tgt = E0` for genes with `|E0|` above the block's 80th percentile:

```text
side = blk & (tgt < 0)
_t2(kf) = minimum(tgt * kf, -shrink_floor)  on side, else tgt
```

A protected down gene has `tgt = E0 < 0`, so it stays in `side`, and `_t2` returns `kf * E0` (the floor 0.002 only binds if the product is tinier than that). A protected up gene is not in `side`, so it stays at `E0`. `kf` is solved in `[0, 1]` so the CPM block's arithmetic mass is unchanged (`rc_emit_k.py:967-981`). The delivered LFC, once the feasible dose `λ` is applied, is `E0 * (1 - λ(1 - kf))`, not `E0` (`rc_emit_k.py:1059-1060`).

MDE does not undo it. `base_t` for a block gene is `_t2(kf)` (`rc_emit_k.py:996`), and the add is capped at 0.1 (`VCC_DE_MDE=0.1,0.1,0.3,coarse`). A median 40% cut of a top-quintile `|E0|` is larger than 0.1. `dn_k250_zs` then targets `2.5 ×` the lfc it is handed (`src/vcc2026/rank_mag.py:320-323`, selector `dn` is every down gene, `rank_mag.py:203-204`), so a pass that reaches its target amplifies the missing magnitude.

**Reproduction.** The script execs `rc_emit_k.py` lines 941–982 (dedented) on a 24-gene gated axis with `protect_q=0.8` and an all-down prior. `R(0)=+1675`, `R(1)=−915`, `kf=0.569`. Protected down genes: delivered target −0.684 and −0.598 against E0 −1.201 and −1.051 (ratio = kf). Protected up genes: ratio 1.000, not in `side`.

**Recorded kf, re-read here** (`repro.out.txt` F1b; joint `λ` from the same file):

| audit | prior | blocks | kf min / median / max | effective scale vs E0 |
|---|---|---|---|---|
| 432 `b2` | `clean4q80` (RESULTS.md configuration table) | 900, D/E/F | 0.385 / 0.601 / 0.971, all &lt; 0.99 | 0.385 / 0.605 / 0.971, all &lt; 0.99 |

`λ = 1` on 89.1% of those blocks, `0 < λ < 1` on 10.9%, `λ = 0` on none. The effective scale already folds `λ` in, and it is still below 0.99 everywhere. 429a's A/B/C audit has the same kf range (median 0.608); that file's prior name was not re-checked in this pass, so 432 b2 is the pin.

The tiny end-to-end `rc_hybrid_block_k` call did **not** place counts: `λ` came out 0 on that toy and placement then raised `rc_hybrid_k did not meet the per-cell totals`. That is a bad toy, not a second bug. When `λ > 0` the emitter asserts the placed arithmetic mean against this same target (`rc_emit_k.py:1244-1253`).

**Score.** `nmae` only, in the wrong direction: the strongest down genes in the CPM block are the ones `nmae` is most likely to be reading, and they are delivered at about 60% of the LFC the comment claims, before the ×2.5 pass. `pds` / `mse` should not move (totals are held). Ordering members can move only insofar as the class histograms move. 394 already failed to separate q80 from the full-block flip on `nmae` (−0.0012, inside the CI). That is what you see if the "protected" magnitudes are still shrunk.

**Fix.** `patches/q80_exclude_protected_from_shrink.patch`. Keep protected genes out of `side`, so `_t2` returns `E0` and their move is 0. Mass balance then shrinks only genes whose target actually changed. Unset `protect_q` is untouched. Applying it changes 446b; it needs a known-answer rebuild and a local `nmae` before any slot.

## 2. CONFIRMED mechanism, latent on A/B/C — prior gate no-op when the quantile is ≤ 0

`scripts/make_submission_pm_cd4.py:5516-5517`. Strength is an integer count. The cut is `strength < np.quantile(strength, q)` with numpy's default linear quantile. Ties at the threshold are excluded as a block.

On the built 458g gate this did what 468 wrote down: threshold 2 in every context, weak 86/83/83, which is strength ≤ 1, the same integer cut as H1's threshold 1.72. That part is not a surprise relative to `experiments/468-prior-gate/RESULTS.md`. Two sharper edges:

- A rank reading of "the weakest 33%" is 99 targets. Strict `<` at a threshold of 2 drops everyone sitting on 2. The stand-in in the script (50 zeros, 36 ones, 20 twos, …) reproduces threshold 2 and 86 weak, with 20 targets tied at 2 left on the programme.
- If the interpolated quantile is 0 (a new panel whose bottom third is all zeros), `strength < 0` is empty and the gate silently keeps the programme for every target. The script shows that no-op. Not the 446b path (the gate is unset there) and not 458g (threshold 2). This is the failure mode to refuse before D/E/F, which is teammate 481's panel but a one-line guard in code we already run.

**Fix.** `patches/prior_gate_refuse_zero_threshold.patch` raises if the threshold is non-finite or ≤ 0. It does not change a threshold of 2. It does not pull the strength-2 ties into the weak set; 468 treated that cut as intentional.

## 3. CONFIRMED residual — zero-repair misses the −5 floor on the edges

`src/vcc2026/zero_repair.py:40-51`, called at `scripts/make_submission_pm_cd4.py:5767-5772`. The repair is real: a gated gene with a zero block total no longer sits at lfc ≈ −35. The comment says the delivered lfc "lands near `floor_log2`". Against the scorer's mean over all 400 cells (`count / (library/1e6)`), a synthetic with 384 carriers, library 20_000, floor −5:

| control CPM | counts placed | scorer lfc | versus −5 |
|---:|---:|---:|---:|
| 8 | 2 | −5.000 | 0 |
| 6 | 1 (rounds) | −5.585 | −0.59 |
| 80 | 19 | −5.074 | −0.07 |
| 400 | 96 | −5.059 | −0.06 |
| 8000 | 384 (capped at `n_car`) | −7.381 | −2.38 |

The cap `n = min(n, n_car)` binds when the gene's expected count per cell exceeds `2^(-floor)` (32 at floor −5). A highly expressed total knockout is left several log2 short of −5, still far from −35. One such gene is a small `nmae` term unless many of the 24 zero blocks 428 found were in this regime; their CPMs were not on this VM.

Separate latent bug, same function: if the sub-gate cannot pay, the cell is skipped (`continue`) but the gene is still appended to `names` (`zero_repair.py:47-51`). The diagnostic says repaired, the lfc stays ≈ −35, and the "other genes unchanged" assert excludes that column so it does not fire. Real blocks have a large sub-gate, so this did not show up in 428's audit. No patch: the residual is not worth a slot next to item 1.

`rc_control.ref_lin` passed into `zero_repair` is already the gated-axis CPM, which matches `ctrl_cpm_gated[j]` for `j` in `flatnonzero(gated)`. That call is aligned. Not a bug.

## 4. CONFIRMED mechanism, low score impact — two different "CPM > 5" gates

`scripts/make_submission_pm_cd4.py:4838-4841` (`emission_rank_axis`, pseudobulk CPM `1e6 * sum(counts) / sum(libraries)`) versus `src/vcc2026/rc_emit.py:745-752` (mean of per-cell CPM). The docstring at `rc_emit.py:684-686` says a gene outside the builder gate can never be called. That is false when the two gates disagree: the gene is inside `control.keep`, its key was zeroed, and `signs[signs == 0] = 1` (`rc_emit_k.py:898-899`) forces it up.

A 200×400 toy with library–expression correlation produces genes that pass one gate and not the other (`repro.out.txt` F4). The builder's own comment (`make_submission_pm_cd4.py:5548-5549`) says the two profiles differ by an sd of 0.020 in log2 on the real controls, so only genes within about 1.4% of CPM 5 can flip. Those are the bottom of the DE universe, not the top-2200 CPM block. `fid` only scores called genes. **SUSPECTED** negligible on A/B/C. A new line with a stronger library–expression correlation is 481's territory.

## 5. CONFIRMED crash, not a silent count — ballast remainder versus the 1e6 cap

`scripts/make_submission_pm_cd4.py:982-998`. After the trim to `k * 999_000`, remainders are rotated by `off = gene_index % k`. An adversarial remainder (every gene with `index % 16 == 0` has remainder 1, nothing else) puts one cell at **1_000_086** and the assert fires. Twenty random compositions at the same trim cap stayed at or below **999_091**. Succeeded production builds did not hit the adversarial case. A new panel could crash a build rather than emit an illegal count. 481. No patch.

## 6. CONFIRMED mechanism, small — the target gene can take a call-set head slot

`rc_emit_k.py:867-871` ranks `|key|` over the whole DE universe. The scorer drops the perturbation's own target from every DE metric (`docs/reference/vcc2026-metrics-BRIEF.md`, "Gene restrictions") and drops all 300 panel genes from `pds`. If the target is the largest `|H|`, it occupies one of the 400 head slots and one real gene does not enter the prescribed set. The script shows a prescribed set of 15 containing the target and 14 genes the DE metrics can see. Yield `min(1, n_pred/n_real)` is usually already capped, so the `fid` cost is about one gene. Not worth a slot by itself.

## Checked, not findings

- Builder unit tests (`experiments/410-pertscale-revisit/run_builder_tests.py`): **56 passed, 0 failed, 1 skipped** (`test_panel_exclusion_for_operator_and_cd4_residual`, skipped by that runner as stale since 481cd092). They synthesize their own controls.
- `experiments/417-bulk-sharpen/test_bulkio417.py` and `test_carrier417.py`: pass. `test_sharpen417.py` and `test_next417.py` stop on missing `data/challenge/pert_counts.csv`.
- `experiments/468-prior-gate/test_gate468.py`: does not start; `experiments/458-arc-lab-programme/work/arch1d_programme_abc.npz` is absent (it hashes that file at import).
- `experiments/451-sameline-mode/test_sameline451.py`: parse and refusal cases pass, then `data/challenge/gene_names.csv` is missing. The blend section never ran.
- `data/challenge`, `data/donors`, `data/harness` are not on this VM. No board number was recomputed.
- Emitter `wilcoxon_z` versus scanpy `rank_genes_groups` on one 80×25 toy (int counts, scanpy then float32): max |Δz| = 5.65e-3, median 6.3e-5, on a gene whose z is 6.12. The gap from `z_cap = 1.90` to the normal 0.05 critical value is about 0.06, so this does not un-suppress the complement. Not ranked.
- `to_csr` drops `block <= 0` (no explicit CSR zeros). Per-cell counts at or below 1e6 are exact in float32. Jitter seed `1000 + row // 400` is unique per (context, target) because `row` is global. `class_cells`' remainder slice decrements the same prefix it checks. Ballast jackknife noise sums to 0 per pair, so it cannot push a library over the cap by itself.
- This pass did not re-test the macOS bash-3.2 / BSD-grep waiter traps in AGENTS.md §6.2. The bugs above are in Python that runs the same on Linux.

## What to do

1. Treat item 1 as a preregistered `nmae` arm, not as a quiet rebuild of 446b. The patch is the whole diff.
2. Item 2's guard is safe to land with the October builder; it does not change 458g.
3. Items 3–6 do not earn a slot on the evidence here.
