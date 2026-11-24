# 093 — does the emission `spread` channel actually achieve the call-set truncation?

Experiment 092 established the *counterfactual*: our submissions call ~17x more genes
significant than the reference (median `n_pred` 5340 against `n_real` 304, all 132/132
perturbations over-calling), `min(1, n_pred/n_real)` is therefore always 1, `fid`
degenerates to pure directional precision `k / n_pred`, and truncating our call set to the
top ~850 by our own adjusted p-value raises `fid` 0.5313 → 0.5906, worth **+0.0299** on the
six-metric average with `reach` and `nmae` exactly invariant.

092 achieved that by fiat — it relabelled significance in a table. This experiment asks
whether the only mechanism that could produce it for real, `emit_counts(mode="quasi")`'s
per-gene `spread`, actually does.

**Answer: yes, but it saturates at roughly a third of the way there, and that is still worth
+0.0203.** The channel cannot reach 850; the best it can do is `n_pred` ≈ 2300.

## Harness

`data/harness/orion-hct116.h5ad` (132 VCC-panel targets in HCT116, 18,400 real control
cells, control `non-targeting`, `obs['target_gene']`), scored against
`experiments/031-somatic/out_som_gamma/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad`.
Production geometry is identical — 400 cells/target, 18,400 controls, 18,533 genes — so the
numbers transfer without rescaling.

The re-emission needs nothing from the model. In the quasi emitter the only thing carried
from the profile into the block is the per-gene panel total `floor(n * L * p_g)`, and the
emitter preserves that integer exactly whatever `spread` is. So the totals are read straight
off the existing prediction as its per-group column sums, and feeding
`p := (totals + 0.5) / (n * L)` back into the **unmodified production `emit_counts`**
recovers them bit-exactly. DE then runs through `experiments/092/de_kernel.py`.

Validation of the whole path, on the production prediction:

| metric | this harness | official scorer (`out_som_gamma/results.json`) | agreement |
|---|---|---|---|
| `fid`   | 0.5312827466125417  | 0.5312814007682398 | 2.5e-6 rel |
| `jac`   | 0.034028361902305626 | 0.0340283619023056 | exact |
| `reach` | 0.014651536921379217 | 0.0146515369213792 | exact |
| `nmae`  | 0.9594589944322061  | 0.9594589939803294 | 4.7e-10 abs |

Re-emitting at the production spread 0.25 with a fresh seed reproduces it to +0.0001 avg6.

## 1. `pds` and `mse` cannot move (`s6_pseudobulk.py`)

cell_eval2 0.16.0 builds the expression comparator from RAW per-group column sums
(`prep._grouped_sums` → `prep.bulk_lognorm_means`), not from per-cell normalized values.
Checked with cell_eval2's own functions across production / global 1.0 / keep-850:

```
group column sums identical (integer):  True
bulk_lognorm_means identical:           True   max|Δ| = 0.000e+00
```

The one cell-level term in that family is `moments.jackknife_correction`, which is
*subtracted* from the raw mse. More spread raises it (`C_p` 10.49 → 15.58), so it can only
lower `expr_mse`; and the official mse is clamped and already sits at 0 from far above its
own baseline (2.175 normalized against 0.989), so it does not move either way.

Every arm in the sweep asserts `np.array_equal(block.sum(0).astype(int64), totals)` per
perturbation. None ever failed.

## 2. `|predicted effect|` is a weak stand-in for our own p-value rank (`s1`, `s2`)

Ranking the call set by `|predicted effect|` — the only confidence ordering available
before a cell is emitted — recovers almost none of the counterfactual:

| truncation to m | fid by our p_adj rank | Δavg6 | fid by \|effect\| rank | Δavg6 |
|---|---|---|---|---|
| 850  | 0.590557 | **+0.0299** | 0.551118 | +0.0041 |
| 2000 | 0.573970 | +0.0225 | 0.554682 | +0.0096 |

Top-850 set overlap between the two orderings is 0.37. A quantile-displacement key
(`|F_ctrl(v_g) − ½|`, the statistic the Wilcoxon actually reduces to under a low-variance
emitter) is no better — best +0.0142, overlap 0.10, because it saturates at 0/1 for
thousands of genes at once.

**This turns out not to matter**, because the spread channel does not *choose* a set. It
raises every gene's p-value and lets the weakest evidence fall out first, which is close to
the p_adj order by construction. That is why the measured result lands near the p_adj
counterfactual curve rather than the `|effect|` one.

## 3. The achieved `n_pred` distribution (`s11_npred.py`)

`n_real`: median 304, mean 329, range [3, 1549].

| policy | min | p10 | p25 | med | p75 | p90 | max | mean | median n_pred/n_real | over-calling |
|---|---|---|---|---|---|---|---|---|---|---|
| production (0.25) | 4616 | 4997 | 5164 | **5340** | 5539 | 5727 | 6096 | 5355 | 17.5x | 132/132 |
| global 0.40 | 3460 | 4026 | 4191 | 4402 | 4659 | 4903 | 5428 | 4435 | 15.0x | 132/132 |
| global 0.55 | 2382 | 3090 | 3314 | 3588 | 3869 | 4163 | 4913 | 3606 | 12.6x | 132/132 |
| global 0.70 | 1663 | 2397 | 2642 | 2957 | 3298 | 3606 | 4547 | 2971 | 10.4x | 132/132 |
| global 0.85 | 1143 | 1894 | 2218 | 2528 | 2932 | 3223 | 4342 | 2556 | 8.8x | 132/132 |
| **global 1.00** | 961 | 1722 | 2054 | **2398** | 2775 | 3089 | 4232 | 2403 | 8.3x | 132/132 |
| keep-850, tail 1.0 | 1092 | 1711 | 2024 | **2317** | 2658 | 2981 | 4164 | 2345 | 8.1x | 132/132 |

`spread` is clipped to 1.0 inside `emit_counts` — Poisson is the noisiest a total-preserving
scatter can be — so **2300 is the floor of this channel**, a 57% cut but still 8x the
reference. Every perturbation still over-calls, so `min(1, n_pred/n_real)` stays pinned at 1
and `fid` stays pure precision; the gain comes entirely from `k/n_pred` improving as the
weakest calls leave.

Holding a top-K tight can only push `n_pred` back *up* relative to the same tail spread, so
no policy in this family goes below the global-1.0 floor. (`keep-850` lands slightly *below*
global 1.0 at 2317 vs 2398 only because protecting 850 genes at spread 0.25 reallocates
scatter: the tail genes it does not protect are unchanged, while the protected block stops
donating its own noise to the per-cell library, marginally sharpening the rest.)

## 4. Six metrics, all 28 arms, one seed (`s4_sweep.py`, `s5_metrics.py`)

`pds` and `mse` are omitted because they are bit-identical everywhere (§1); Δavg6 divides
the four DE deltas by 6.

| policy | med n_pred | fid | jac | reach | nmae | Δfid | Δjac | Δreach | Δnmae | **Δavg6** |
|---|---|---|---|---|---|---|---|---|---|---|
| production | 5340 | 0.531283 | 0.034028 | 0.014652 | 0.959459 | — | — | — | — | — |
| global 0.25 (control) | 5352 | 0.531358 | 0.033909 | 0.015254 | 0.959464 | +0.0003 | −0.0003 | +0.0007 | −0.0000 | +0.0001 |
| global 0.40 | 4402 | 0.541427 | 0.034196 | 0.016485 | 0.959456 | +0.0338 | +0.0005 | +0.0020 | +0.0000 | +0.0061 |
| global 0.55 | 3588 | 0.553141 | 0.034606 | 0.017217 | 0.959455 | +0.0729 | +0.0016 | +0.0029 | +0.0000 | +0.0129 |
| global 0.70 | 2957 | 0.559686 | 0.035246 | 0.017168 | 0.959446 | +0.0947 | +0.0033 | +0.0028 | +0.0000 | +0.0168 |
| global 0.85 | 2528 | 0.562449 | 0.035776 | 0.018212 | 0.959462 | +0.1039 | +0.0047 | +0.0040 | −0.0000 | +0.0188 |
| global 1.00 | 2398 | 0.564459 | 0.035625 | 0.018649 | 0.959457 | +0.1106 | +0.0043 | +0.0045 | +0.0000 | +0.0199 |
| effK400 s0.6 | 3316 | 0.556298 | 0.034899 | 0.017303 | 0.959451 | +0.0834 | +0.0024 | +0.0030 | +0.0000 | +0.0148 |
| effK400 s0.8 | 2594 | 0.562802 | 0.035598 | 0.017355 | 0.959455 | +0.1051 | +0.0042 | +0.0030 | +0.0000 | +0.0187 |
| effK400 s1.0 | 2342 | 0.563985 | 0.035461 | 0.017318 | 0.959458 | +0.1090 | +0.0039 | +0.0030 | +0.0000 | +0.0193 |
| effK600 s0.6 | 3301 | 0.556581 | 0.034941 | 0.017263 | 0.959467 | +0.0843 | +0.0025 | +0.0029 | −0.0000 | +0.0149 |
| effK600 s0.8 | 2581 | 0.563258 | 0.035442 | 0.016652 | 0.959457 | +0.1066 | +0.0038 | +0.0022 | +0.0000 | +0.0188 |
| **effK600 s1.0** | 2308 | 0.565869 | 0.035488 | 0.017880 | 0.959461 | +0.1153 | +0.0039 | +0.0036 | −0.0000 | **+0.0205** |
| effK850 s0.6 | 3304 | 0.556370 | 0.034743 | 0.017037 | 0.959467 | +0.0836 | +0.0019 | +0.0027 | −0.0000 | +0.0147 |
| effK850 s0.8 | 2590 | 0.562966 | 0.035177 | 0.017038 | 0.959450 | +0.1056 | +0.0031 | +0.0027 | +0.0000 | +0.0186 |
| **effK850 s1.0** | 2317 | 0.565983 | 0.035114 | 0.017459 | 0.959463 | +0.1157 | +0.0029 | +0.0031 | −0.0000 | **+0.0203** |
| effK1200 s0.6 | 3309 | 0.556593 | 0.034581 | 0.017093 | 0.959454 | +0.0844 | +0.0015 | +0.0027 | +0.0000 | +0.0148 |
| effK1200 s0.8 | 2586 | 0.563191 | 0.034929 | 0.016958 | 0.959455 | +0.1064 | +0.0024 | +0.0026 | +0.0000 | +0.0186 |
| effK1200 s1.0 | 2325 | 0.565340 | 0.034911 | 0.017639 | 0.959463 | +0.1135 | +0.0024 | +0.0033 | −0.0000 | +0.0199 |
| effK2000 s0.6 | 3337 | 0.556562 | 0.034259 | 0.016303 | 0.959461 | +0.0843 | +0.0006 | +0.0018 | −0.0000 | +0.0145 |
| effK2000 s0.8 | 2632 | 0.562449 | 0.034610 | 0.016203 | 0.959441 | +0.1039 | +0.0016 | +0.0017 | +0.0000 | +0.0179 |
| effK2000 s1.0 | 2370 | 0.564666 | 0.034469 | 0.017251 | 0.959452 | +0.1113 | +0.0012 | +0.0029 | +0.0000 | +0.0192 |
| effP15 s0.6 | 3322 | 0.557005 | 0.034791 | 0.017084 | 0.959457 | +0.0857 | +0.0021 | +0.0027 | +0.0000 | +0.0151 |
| effP15 s0.8 | 2587 | 0.563660 | 0.035160 | 0.016637 | 0.959458 | +0.1079 | +0.0031 | +0.0022 | +0.0000 | +0.0189 |
| effP15 s1.0 | 2322 | 0.565723 | 0.035180 | 0.017160 | 0.959463 | +0.1148 | +0.0031 | +0.0028 | −0.0000 | +0.0201 |
| padjK850 s1.0 *(oracle)* | 2378 | 0.563950 | 0.035463 | 0.016795 | 0.959453 | +0.1089 | +0.0039 | +0.0024 | +0.0000 | +0.0192 |
| padjK2000 s1.0 *(oracle)* | 2690 | 0.562085 | 0.035025 | 0.015493 | 0.959433 | +0.1027 | +0.0027 | +0.0009 | +0.0000 | +0.0177 |
| tightK850 s1.0 | 2270 | 0.566492 | 0.034958 | 0.016415 | 0.959458 | +0.1174 | +0.0025 | +0.0020 | +0.0000 | +0.0203 |
| tightK2000 s1.0 | 2396 | 0.564821 | 0.033726 | 0.016395 | 0.959463 | +0.1118 | −0.0008 | +0.0019 | −0.0000 | +0.0188 |

`effP15` is the scale-adaptive "top 15% of our own call list" form; `tightK` holds the
top-K at spread 0.0 rather than 0.25; `padjK` is an oracle arm that ranks by the production
run's own adjusted p-value (not buildable — it needs the emission it is choosing for).

Three things to read off this table.

* **`nmae` is invariant to five decimal places on every arm.** Largest |Δnmae| raw is
  2.7e-5, largest scaled 0.0001. This is the whole point of raising the p-value instead of
  nulling the effect: 092 measured that flattening the tail's lfc costs `nmae` 0.9595 →
  1.0154, which is −0.093 scaled and wipes out half the gain.
* **`reach` and `jac` both improve slightly** (+0.0045 and +0.0043 scaled at global 1.0)
  rather than staying flat as the 092 relabelling counterfactual predicted. Dropping the
  weakest calls removes genes that were mostly *not* in the reference set, so the Jaccard
  intersection survives better than the union shrinks.
* **The oracle arm loses.** Choosing the kept set by our own p_adj — the ordering that wins
  the 092 counterfactual by 7x — scores *below* the plain global arm (+0.0192 vs +0.0199).
  Protecting the genes that were already most significant spends the protection where the
  spread was never going to remove them. So there is no better pre-emission ranking key
  worth chasing: the per-gene channel has essentially no headroom over the global one.

## 5. Seed repeats on the leading arms (`s9_seeds.py`, `s10_seed_table.py`)

Four independent emission seeds each.

| policy | med n_pred | fid | jac | reach | **Δavg6** | vs global 1.0 |
|---|---|---|---|---|---|---|
| effK600 s1.0 | 2320 | 0.565571 ± 0.0002 | 0.035554 ± 0.0001 | 0.017543 ± 0.0002 | **+0.0203 ± 0.0001** | +0.0006 (+3.3 sd) |
| effK850 s1.0 | 2307 | 0.565868 ± 0.0003 | 0.035185 ± 0.0001 | 0.017104 ± 0.0003 | **+0.0202 ± 0.0002** | +0.0005 (+2.6 sd) |
| tightK850 s1.0 | 2285 | 0.566075 ± 0.0004 | 0.034888 ± 0.0002 | 0.016502 ± 0.0002 | +0.0201 ± 0.0003 | +0.0004 (+1.6 sd) |
| global 1.0 | 2389 | 0.564126 ± 0.0003 | 0.035592 ± 0.0001 | 0.018369 ± 0.0005 | +0.0196 ± 0.0002 | — |
| padjK850 s1.0 | 2390 | 0.563931 ± 0.0002 | 0.035602 ± 0.0001 | 0.016899 ± 0.0004 | +0.0193 ± 0.0001 | −0.0004 (−2.2 sd) |
| global 0.85 | 2538 | 0.562466 ± 0.0002 | 0.035691 ± 0.0001 | 0.017741 ± 0.0004 | +0.0186 ± 0.0002 | −0.0010 (−4.7 sd) |

Seed noise is 0.0001–0.0003 on avg6. `effK600` and `effK850` are tied with each other and
both beat the global arm by a small but real margin.

## 6. Recommendation

**`--tail-spread 1.0 --keep-top 850`, expected Δavg6 = +0.0202 ± 0.0002.**

`K = 600` measures 0.0001 higher, which is inside noise; 850 is chosen because the K curve
is a broad plateau (400 → +0.0193, 600 → +0.0205, 850 → +0.0203, 1200 → +0.0199,
2000 → +0.0192) and 850 sits in its middle with more margin on the high side, and because it
is the value experiment 092 independently landed on. Dropping `--keep-top` entirely still
gets +0.0196 if the simpler one-knob mechanism is preferred.

Expected per-metric move: `fid` 0.5313 → 0.5659, `jac` 0.0340 → 0.0352, `reach` 0.01465 →
0.01710, `nmae` 0.95946 → 0.95946, `pds`/`mse` bit-unchanged.

That is 68% of the +0.0299 the perfect truncation would give. The missing third is the
saturation: the channel stops at `n_pred` ≈ 2300 and the counterfactual optimum is at 850.

## 7. Negative result and follow-up: the clip is the binding constraint (`s7_headroom.py`)

`spread` is clipped to 1.0 because Poisson is the noisiest scatter the current emitter can
draw. Probing past it — gamma-Poisson scatter at dispersion `phi` with the *same* exact
total repair, so the pseudobulk invariance is untouched — measures what an emitter extension
would be worth. This is deliberately outside `src/`; nothing here changes production.

| arm | med n_pred | fid | jac | reach | nmae | Δavg6 |
|---|---|---|---|---|---|---|
| phi 0.00 (= spread 1.0) | 2370 | 0.563659 | 0.035574 | 0.019147 | 0.959440 | +0.0195 |
| phi 0.05 | 2078 | 0.567117 | 0.034464 | 0.022119 | 0.959528 | +0.0215 |
| phi 0.15 | 1574 | 0.568621 | 0.033817 | 0.037831 | 0.959686 | **+0.0249** |
| phi 0.40 | 1559 | 0.543191 | 0.045035 | 0.101269 | 0.960146 | +0.0275 |
| phi 1.00 | 3786 | 0.498972 | 0.048477 | 0.063601 | 0.961093 | −0.0028 |

`phi = 0.15` is the interesting point: `n_pred` 1574, best `fid` of anything measured here,
and no metric degrading. `phi = 0.40` scores higher only because `reach` explodes (+0.097
scaled) while `fid` falls *below* the recommended arm, and `phi = 1.00` is non-monotone —
extreme zero-inflation makes genes significant in the wrong direction and `fid` drops under
baseline. So the extension is worth roughly another +0.005 at `phi ≈ 0.15`, but it needs an
`emit_counts` change plus its own validation of where the reach gain comes from. Not shipped.

## Files

| file | what |
|---|---|
| `lib093.py` | re-emission from panel totals through the unmodified `emit_counts`, the gate, and the exact DE path |
| `s1_rankcheck.py` | `\|predicted effect\|` vs our p_adj ordering; truncation counterfactual under both |
| `s2_betterkey.py` | quantile-displacement key; same counterfactual |
| `s3_pilot.py` | 12-perturbation global-spread ladder, cost and totals-exactness check |
| `s4_sweep.py` | the 28-arm emission sweep |
| `s5_metrics.py` | all six metrics per arm, cell_eval2-exact |
| `s6_pseudobulk.py` | pseudobulk invariance proved with cell_eval2's own functions |
| `s7_headroom.py` / `s8_score_nb.py` | over-Poisson probe past the clip |
| `s9_seeds.py` / `s10_seed_table.py` | seed repeats on the leading arms |
| `s11_npred.py` | the `n_pred` distributions |

## Production change

`scripts/make_submission_pm_cd4.py` gains `emission_spread_vector()` and three flags:
`--tail-spread`, `--keep-top`, `--keep-top-frac`. All three default to off, and with
`--tail-spread` unset the function returns exactly `np.full(profile.size, SPREAD)` — the
literal expression the shipped call used — so an unflagged build is bit-identical to v11.
`--keep-top`/`--keep-top-frac` without `--tail-spread` is rejected rather than silently
ignored. Two tests added to `scripts/test_make_submission_pm_cd4.py`
(17/17 pass, peak RSS 0.132 GiB).
