# 164 — the call-set artefact and `fid` precision: decompose it, then price it on all six

Date: **2026-09-02**. Measurement only. No `cell-eval2`, no build, no `vcc prep`, no
submission, no existing file modified. Everything is written inside this directory.

Decision, in one line: **NO-GO on every call-set arm.** The zero-effect artefact is real,
it is 41.8% of our call set, and its directional precision is **0.4999** — dead chance, as
the thesis predicted. Removing it is still wrong, because the entire `fid` prize for
deleting it *perfectly and for free* is **+0.00089 overall**, while the same tail block
carries **all** of a15's `mse` credit and removing it costs **−0.01852 overall**. The
channel is 20× underwater.

Two inherited claims are corrected here, one in our favour and one against.

---

## 0. What changed, honestly

* The predecessor's headline hypothesis — that the 794/958/997-gene sub-gate band is a
  **normalisation artefact polluting the call set** — is **half right and strategically
  wrong**. The band's scored `lfc` *is* a pure mean-of-ratios/ratio-of-sums artefact and it
  *is* systematically DOWN (53–64%). But it is **never called**: 0.0–0.2% of the band clears
  BH at zero effect under *every* policy, because 5–10 CPM genes have no power at 400 cells
  (median `|z|` 0.3–0.8). The band cannot pollute `fid` or `jac` because it never enters
  them. §2.
* 157 §5b's "the shipped call set is ~30% zero-effect artefact" is **too high**: with the
  correct block composition it is **26.8 / 25.5 / 23.5%**. And its "100% DOWN" was the
  Wilcoxon `z` sign; the **scored** direction is `lfc-up` **42.4 / 44.1 / 43.5%** — the
  artefact is near-balanced in the direction `fid` actually reads. §2.
* The h1-test surface **has no sub-gate band at all** — the 095 artifact was not built at
  `--min-cpm 10`, proven by a statistic quantisation cannot fake. So the h1-test
  decomposition prices production's *gated tail* sub-block only, and the band is priced from
  the challenge controls alone. §3.1.
* `pred_lfc` is **not** arm-invariant. The assertion the predecessor asked for **fails as
  stated** and is replaced by a weaker true one. §3.2.
* The `mse` instrument reproduces a8's board raw `mse` to **+0.61%** and 87% of the whole
  a8→a15 board move. That is the known-answer row, and it is what kills the arms. §5.

---

## 1. Instruments, and the three known-answer gates

| gate | expected | measured | verdict |
| --- | ---: | ---: | --- |
| `a15` median `n_pred`, h1-test | 6035 / 6032 (157) | **6035 / 6032** | bit-identical |
| `armC` (a8 flat) median `n_pred` | 3581 / 3589 (157) | **3581 / 3589** | bit-identical |
| `a15` native-BH `reach` | 0.18283 (157 `power_h1test.log`) | **0.18283** | PASS |
| `armC` native-BH `reach` | 0.20291 (157 `power_h1test.log`) | **0.20291** | PASS |
| block census | 794 / 958 / 997 band (`blocks.log`) | **794 / 958 / 997** | PASS |
| zero-effect rate, `0.25/φ4`, DE universe | 19.0 / 21.6 / 16.0% (153 §15) | **19.3 / 21.4 / 16.1%** | PASS |
| `claim(0.25, φ0)` | 0.11860 (096 `claim_grid2.log`) | **0.11860** | PASS (−0.00%) |
| `claim(1.00, φ0)` | 0.19216 (096 `claim_grid.log`) | **0.19438** | PASS (+1.16%) |
| `claim(0.25, φ4)` = a8 as shipped | 0.16478 / 0.16676 (096) | **0.16652** | PASS (−0.14% vs board_solve) |
| a8 raw `mse` | **1.3627 observed on the board** | **1.3711 modelled** | PASS (+0.61%) |

The `claim` grid tolerance is 2%, set by the documented spread between 096's own two grids
(`claim_grid` 0.116800 vs `claid_grid2` 0.11860 on the same policy = 1.5%).

**Surfaces.** h1-test (85 perturbations × 10,779 gated genes), donors `k562_gwps` +
`orion_hek293t` — `h1_2025` covers 85/85 targets and never enters. The real challenge
controls `data/challenge/context_{A,B,C}.h5ad` for everything that can be computed without
truth. `claim` on the real challenge panel via 096's cached `v11-stats` + `cp-proxy-ntc` +
`board_solve`.

---

## 2. The sub-gate band: the artefact is real, and it is invisible

`blocks.py` (predecessor) established that the scorer's DE universe (**mean-per-cell** CPM >
5) strictly contains production's spread gate (**pseudobulk** CPM > 10, `--min-cpm 10.0`);
`prod-gated & not DE` is **0** in all three contexts. `profile_from_effect`
(`make_submission_pm_cd4.py:1281`) sets the predicted effect to **exactly zero** outside that
gate, so for the band genes our scored `lfc` is a pure normalisation identity:

```
lfc_0(g) = log2( pb_cpm_g / cell_cpm_g )      ratio-of-sums / mean-of-ratios, same cells
```

`zero_lfc.py`, exact, controls only:

| ctx | DE universe | prod gate | **band** | median `lfc_0` in band | band **DOWN** | `\|lfc_0\|>0.1` |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| A | 9929 | 9135 | **794** | −0.0152 | **64.2%** | 4.7% |
| B | 9626 | 8668 | **958** | −0.0146 | **62.8%** | 11.1% |
| C | 10124 | 9127 | **997** | −0.0031 | **53.1%** | 2.6% |

**The sign hypothesis is confirmed and the magnitude hypothesis is refuted.** The artefact is
systematically DOWN, but its median size is 0.003–0.015 log2 units, and — decisively — the
band is **never called**:

| policy | band call rate at zero effect (A/B/C) | median `\|z\|` |
| --- | ---: | ---: |
| `0.25/φ4` | 0.1% / 0.1% / 0.2% | 0.3–0.4 |
| `1.00/φ4` (a15's tail) | 0.1% / 0.1% / 0.2% | 0.6–0.8 |
| `1.00/φ0` | 0.0% / 0.1% / 0.1% | 0.4–0.5 |
| `0.50/φ4` | 0.0% / 0.1% / 0.2% | 0.3 |

0–2 genes out of 794/958/997. **A 5–10 CPM gene has no power at 400 emitted cells against
18,400 controls, whatever the emission law does to it.** The band is invisible to `fid`,
`jac` and `nmae`'s gate. It is not a removable dilutant, because it is not in the call set.

### The artefact that *is* in the call set is the gated tail, and its scored direction is balanced

| block, policy | call rate at zero effect (A/B/C) | `z`-up | **`lfc`-up (what `fid` scores)** |
| --- | ---: | ---: | ---: |
| gated, `0.25/φ4` | 21.0 / 23.8 / 17.8% | 87.5 / 86.3 / 81.0% | **53.6 / 56.4 / 56.4%** |
| gated, **`1.00/φ4`** (a15's tail) | **86.1 / 82.2 / 84.4%** | **0.0%** | **42.4 / 44.1 / 43.5%** |
| gated, `1.00/φ0` | 8.4 / 11.2 / 9.1% | ~100% | 67.9 / 68.3 / 67.6% |
| gated, `0.50/φ4` | 33.9 / 31.5 / 34.5% | 3.7 / 7.1 / 7.7% | 47.6 / 49.5 / 40.5% |

This settles the predecessor's reframing **by measurement**: 157's "100% DOWN" is the `z`
sign and is irrelevant to `fid`; the scored direction is 42–44% up, i.e. near-balanced. That
is exactly why the block reads precision 0.4999 on h1-test (§3.3).

**Corrected shipped mixture** (157 §5b said 30.6 / 31.4 / 28.8%):

| policy | zero-effect calls as a share of the DE universe (A/B/C) |
| --- | ---: |
| **a15 as shipped** | **26.8 / 25.5 / 23.5%** |
| armA tail φ0 | 17.9 / 20.6 / 15.1% |
| armB tail 0.5 | 20.8 / 22.0 / 17.9% |
| armC a8 flat | 19.4 / 21.5 / 16.1% |

---

## 3. h1-test: the precision decomposition

### 3.1 The h1-test surface has no sub-gate band — proven, not assumed

`decomp.py`'s first detector (modal `pred_pb/ctrl_pb` spike) **failed**: it returned 98 genes
with control pseudobulk CPM up to 105, which is integer coincidence, not a CPM gate. Its
output is discarded. A second comparison (median `|effect|` in the band 0.173 vs 0.123 in the
10–20 CPM slice) is **also unsafe**, because the emitter floors each gene's panel total and
that noise grows exactly where the band lives.

`verify_band.py` settles it with a statistic quantisation cannot fake — the **between-target**
spread of a gene's predicted pseudobulk. A zeroed gene has the same profile for every target,
so its spread is bounded by the flooring slack; a gene with a real donor effect moves.

| block | n | median between-target spread | flooring bound | ratio |
| --- | ---: | ---: | ---: | ---: |
| **sub-gate band** | 1181 | 0.7782 | 0.0641 | **12.14** |
| pb-CPM 10–20 | 1551 | 0.5484 | 0.0301 | 18.25 |
| pb-CPM 20–50 | 2812 | 0.3920 | 0.0126 | 31.11 |
| pb-CPM > 50 | 5235 | 0.2970 | 0.0037 | 80.71 |

A zeroed band would read ≈ 1. It reads **12.14**. **The 095 h1-test artifact was not built at
`--min-cpm 10`**, so the challenge's 794/958/997-gene block has **no h1-test analogue**, and
everything in §3.3 prices production's *gated tail* only. (`realize164.py` also ranks by
`|log2(pcpm/cpm)|` over a cell-CPM-5 gate rather than production's `|effect|` over a
pb-CPM-10 gate; that is 157's board-validated convention and is kept, but it is a proxy.)

### 3.2 `pred_lfc` invariance — the assertion FAILS as stated

| comparison vs a15 | max \|Δ\| | median \|Δ\| | sign flips / 916,215 |
| --- | ---: | ---: | ---: |
| **a15 vs a15_b (same arm, re-seed)** | 9.22e−02 | 2.06e−03 | **14,185 (1.55%)** |
| armA | 6.44e−02 | 1.71e−03 | 10,527 (1.15%) |
| armB | 8.11e−02 | 1.75e−03 | 10,851 (1.18%) |
| armC | 6.00e−02 | 1.72e−03 | 10,464 (1.14%) |

`pred_lfc` is **not** invariant — 1.1% of gene directions flip between arms, because per-cell
library sums differ slightly even with panel totals repaired exactly. **But the re-seeded twin
of a15 flips MORE (1.55%) than any arm does.** So the correct statement is: *spread/phi arms
move no more direction than emitter seed noise does; they change which genes are called, and
any direction change they cause is not distinguishable from re-running the same arm.* The
strong form ("never any gene's predicted DIRECTION") is false and should not be repeated.
This is why the twin was run — the check would otherwise have had no scale to be read against.

### 3.3 The decomposition

`decomp.py` / `decomp2.py`, a15, call set ∩ adjudicable, split by the saved per-target spread
vector. Pooled and per-perturbation agree to 4 dp, so rule 8's trap is not active here:

| block | calls | share | pooled `k/n` | `mean_p(k/n)` |
| --- | ---: | ---: | ---: | ---: |
| keep (spread 0.25) | 297,869 | 58.2% | 0.5122 | 0.5120 |
| **tail (spread 1.0)** | **214,265** | **41.8%** | **0.4999** | **0.4999** |
| ALL | 512,134 | 100% | 0.5070 | 0.5070 |

The tail block's per-perturbation precision is 0.4999 with an sd of the mean of **0.0011** —
**indistinguishable from chance**. Paired keep − tail = **+0.0122, t = 4.33**. The thesis's
central empirical claim is **confirmed**: ~42% of our call set carries no directional signal.

**But h1-test `fid` must not be read directly.** h1-test's `n_real` is median 89 / **mean
765** — a heavy tail that makes `max(n_pred, n_real)` bind for some perturbations once the
call set shrinks. That is AGENTS.md rule 9 exactly. The challenge-readable statistic is
**precision** = `mean_p(k/n_pred)`, since the challenge has `n_pred` ≈ 4,424 against
`n_real` ≈ 240 for every perturbation:

| arm | `n_pred` | **precision** | h1-test `fid` (contaminated) | reach | jac | nmae | pds |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| a15 | 6021 | 0.50735 | 0.50736 | 0.17910 | 0.06829 | 0.96650 | 0.79979 |
| armA φ0 tail | 3123 | **0.51305** | 0.50460 | 0.20882 | 0.05389 | 0.96665 | 0.79881 |
| armB tail 0.5 | 4663 | 0.50900 | 0.50892 | 0.18284 | 0.06731 | 0.96659 | 0.79860 |
| armC a8 flat | 3561 | 0.51230 | 0.50793 | 0.19757 | 0.06088 | 0.96668 | 0.79895 |
| *seed spread (a15)* | | 0.00070 | 0.00070 | 0.00746 | 0.00010 | 0.00018 | 0.00070 |

Precision rises exactly as predicted; h1-test `fid` moves the *opposite* way for armA
(−0.00276) purely through the yield cap. `nmae` and `pds` move within seed noise, confirming
the exact-total repair.

---

## 4. Does it reconcile with board `fid` 0.5022? **Yes — and the reconciliation is what kills it**

Using 096's conserved-headroom form `p_after − 0.5 = (p0 − 0.5)·n_before/n_after`, **not** a
transfer of harness precision levels.

**First, validate the model on h1-test** (it could have failed — a block with negative
headroom would break it):

| arm | `n_pred` ratio | predicted precision | observed | realised |
| --- | ---: | ---: | ---: | ---: |
| armA | 1.9281 | 0.51417 | 0.51305 | **92.1%** |
| armB | 1.2913 | 0.50949 | 0.50900 | **94.8%** |
| armC | 1.6910 | 0.51243 | 0.51230 | **99.0%** |

**Then the reconciliation.** h1-test headroom is +0.00735, the board's is +0.00220 — a single
scale factor of **0.299**. If the board's tail block holds the same 41.8% share at exactly
0.5, the board's *signal* block must sit at **0.50378**. h1-test's signal block headroom
(+0.01204) × the same 0.299 = **+0.00360** against the required +0.00378. **One factor maps
h1-test's overall precision AND its per-block precision onto the board.** It reconciles.

**The ceiling this implies is the whole finding.** Delete the entire tail block, perfectly and
for free:

```
board fid  0.50220 -> 0.50378      d raw +0.00158
                                   d scaled +0.00532      d OVERALL +0.000886
```

Our gap to the board **maximum** `fid` (0.5175, `The_Way`) is +0.01530 raw = **+0.008579
overall — 9.7× the entire artefact channel.** The artefact is real, it is exactly as large as
advertised, and deleting all of it does not move us.

**The board has already tested this family four times and every shrink lost.** Same profile,
emission-only changes (`v11 → v15 → v14`, per 096's own note that v15 has v11's profile
exactly):

| artifact | challenge `n_pred` (measured) | board raw `fid` |
| --- | ---: | ---: |
| v14 | 1,563 | 0.49126 |
| v15 | 1,911 | 0.49654 |
| a8 | 3,836 | 0.5021 |
| **a15** | **4,424** (`096/work/probe-a15.json`) | **0.5022** |
| v11 | 5,961 | 0.50178 |

`fid` is flat-to-rising in `n_pred` across the whole measured range and saturates near 0.502
above ~4,000. Under either reading of `n_real` (240 → pure precision; ~2,670 → coverage) the
**sign is the same: shrinking the call set has never once raised board `fid`.** The
conserved-headroom projection is the *optimistic* bound and the board contradicts even its
sign.

---

## 5. `mse`: the tail block is the entire credit, and a15 is on the cliff

`claim164.py`, real challenge panel, per-gene `(spread, φ)` built exactly as
`emission_spread_vector` builds it (including the trap that the vector is **initialised to
`tail_spread`**, so sub-gate genes sit in the tail too).

| policy | gated | `Ĉ_p` | `Ĉ/C_p` | **`claim`** | Δclaim vs a15 |
| --- | ---: | ---: | ---: | ---: | ---: |
| flat `0.25/φ0` *(known)* | 8976 | 7.33 | 0.481 | 0.11860 | −0.12721 |
| flat `1.00/φ0` *(known)* | 8976 | 12.01 | 0.788 | 0.19438 | −0.05143 |
| flat `0.25/φ4` = a8 *(known)* | 8976 | 10.29 | 0.675 | 0.16652 | −0.07929 |
| **a15 as shipped** | 8976 | 16.53 | **1.085** | **0.24581** | — |
| armA tail φ0 | 8976 | 10.56 | 0.693 | 0.17087 | **−0.07494** |
| armB tail 0.5 | 8976 | 11.69 | 0.767 | 0.18917 | −0.05664 |
| armC a8 flat | 8976 | 10.29 | 0.675 | 0.16652 | −0.07929 |
| armB′ tail 0.75 | 8976 | 13.76 | 0.904 | 0.22279 | −0.02302 |
| armA′ tail φ1 | 8976 | 12.05 | 0.791 | 0.19512 | −0.05069 |
| a15 @ `--min-cpm 5` (spread only) | 9885 | 16.81 | 1.104 | 0.24634 | +0.00053 |

**a15's `claim` is 0.24581 against the reference cap `K` = 0.24653 — it is at 99.7% of the
ceiling.** That is what a12 was built for and it is fully banked.

At a15's shipped `VCC_ALPHA = 0.982`, `s = 0.90960` and `s²B = 0.23268`, so `ρ = 0.9466`:
**a15 is on the cliff.** Any arm that drops `claim` below 0.23268 lands in the `ρ = 1` branch
where `∂(raw mse)/∂claim = −1/ΣD = −4.375`, with **no amplitude compensation available**
(amplitude is a second variable and outside this mandate):

| arm | `claim` | raw `mse` | `score_mse` | **Δ overall** |
| --- | ---: | ---: | ---: | ---: |
| **a15** | 0.24581 | 0.88228 | 0.11112 | — |
| armB′ tail 0.75 | 0.22279 | 0.92557 | 0.06556 | **−0.00759** |
| armA′ tail φ1 | 0.19512 | 1.04663 | **0.00000** | **−0.01852** |
| armB tail 0.5 | 0.18917 | 1.07266 | **0.00000** | **−0.01852** |
| armA tail φ0 | 0.17087 | 1.15272 | **0.00000** | **−0.01852** |
| armC a8 flat | 0.16652 | 1.17175 | **0.00000** | **−0.01852** |
| a15 @ `--min-cpm 5` (spread only) | 0.24634 | 0.88228 | 0.11112 | +0.00000 |

`score_mse = clip((0.987862879 − raw)/0.950162566, 0, 1)`. **Every arm that removes a material
part of the tail block drives raw `mse` past 0.98786 and clamps `score_mse` to zero** —
surrendering the whole +0.0185 that a15 banked. This is the a19 failure mode, priced before
the fact instead of after.

---

## 6. The six-member table

`d pds = 0.000000` exactly and `d nmae ≈ 0` are **not unmeasured** — they are pinned by
`emit_counts`'s exact-total repair (rule 12's stated exemption, which covers exactly this
`spread`/`φ`/`keep_top` family), and h1-test confirms both inside seed noise (§3.3).

**At a15's shipped `VCC_ALPHA 0.982` — a true one-variable arm:**

| arm | `d pds` | `d mse` | `d fid` | `d reach` | `d jac` | `d nmae` | **NET overall** |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| **armA** tail φ0 | 0.000000 | **−0.018520** | +0.001145 | +0.002394 | −0.000333 | ~0 | **−0.015314** |
| **armB** tail 0.5 | 0.000000 | **−0.018520** | +0.000359 | +0.000301 | +0.000760 | ~0 | **−0.017099** |
| **armC** a8 flat | 0.000000 | **−0.018520** | +0.000852 | +0.001488 | +0.000570 | ~0 | **−0.015610** |
| **armD** post-hoc drop | — | — | ≤ +0.000886 | — | — | — | **not expressible** |
| **`--min-cpm 5.0`** | ? | +0.000026 † | ~0 ‡ | ? | ~0 ‡ | ? | **unpriced — see §7** |

**Known-answer row a8 → a15**, board-measured at **+0.019726** overall. The instrument
predicts `mse` +0.018520, `reach` −0.001488, `fid` +0.000056 → **+0.017088, i.e. 87% of the
observed move**, with the residual +0.002638 attributable to `jac`, `nmae` and the CD4 /
whitener confounds in that comparison. Modelled a8 raw `mse` **1.3711** against the observed
**1.3627** (+0.61%). Signs and magnitudes both reproduce.

Granting each arm its **own** re-optimised `VCC_ALPHA` (a second variable, outside mandate,
and it would need its own α screen):

| arm | `d mse` | `d fid` | `d reach` | `d jac` | **NET** |
| --- | ---: | ---: | ---: | ---: | ---: |
| armA | −0.004010 | +0.001145 | +0.002394 | −0.000333 | **−0.000804** |
| armB | −0.002952 | +0.000359 | +0.000301 | +0.000760 | **−0.001532** |
| armC | −0.004270 | +0.000852 | +0.001488 | +0.000570 | **−0.001361** |

**Every arm is negative on every accounting.** The best case in the whole family is armA with
a re-tuned amplitude at **−0.0008**, which is inside the noise of the α screen it would
require.

Method notes for the table: `d fid` is the conserved-headroom projection onto the board's own
+0.0022 headroom (the optimistic bound, and §4 shows the board contradicts even its sign).
`d reach` is the h1-test delta shrunk by **0.428**, the factor calibrated on armC against the
board's own a8↔a15 `reach` pair (0.1890 vs 0.1811) — armC is the known-answer row for reach.
`d jac` is rebuilt on the board's `n_real` 240 / `n_pred` 4,424 with the overlap scaled by the
h1-test-measured ratio, because raw h1-test `jac` is unreadable under rule 9.

---

## 7. `--min-cpm 5.0` — NO-GO, and it is not in the exempt class

The mandate raised this arm because the band's effects are currently a normalisation artefact
rather than model content. Three findings close it without a rebuild:

1. **It is not a call-set arm.** `profile_from_effect` gates the **effect**, not the spread
   (`make_submission_pm_cd4.py:1281`; `donor_transfer.py:213`). Changing `--min-cpm` changes
   the emitted **profile**, so `A`, `X`, `pds`, `nmae` and `mse` all move. Rule 12's
   exemption explicitly does **not** cover anything touching the profile. My instruments
   (which hold the profile fixed and move only `claim`) **cannot price it**, and the
   `+0.000026` in the table above is the spread-vector half only — the profile half is
   unmeasured. Saying otherwise would be the a19 error again.
2. **The mechanism that motivated it does not exist.** The band is called at **0.0–0.2%**
   under every policy (§2). Genes at 5–10 CPM are statistically invisible at 400 cells, so
   giving them real effects instead of a normalisation artefact cannot recover `fid` or
   `jac` — there is nothing there to recover. The premise of the arm is measured false.
3. **The direct prior is against it and it is not stale in the relevant way.** Exp 088
   measured 10.0 as a sharp optimum over 6 seeds (+0.00455 avg6, sd 0.00058), improving
   `reach`, `pds` and `nmae` *together*. AGENTS.md's warning that a lever closed under one
   objective may reopen under a new one is correct in general, but the new objective here is
   `fid`-as-precision, and (2) shows this flag has **no** channel to `fid`.

**Verdict: NO-GO, and it must not be shipped on the strength of this experiment's
instruments even if someone later wants it — it needs a build plus `mech_chal_stats` on the
new artifact.**

---

## 8. GO / NO-GO

### NO-GO on all of armA, armB, armC and `--min-cpm 5.0`. Do not spend a slot.

**a15 stays exactly as it is:**

```
VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.10 \
  --min-cpm 10.0 --overdispersion 4.0 --tail-spread 1.0 --keep-top 8000 \
  --keep-spread 0.25 --whitener-id union_pruned \
  --whitener data/donors/whitener_union_pruned.npz \
  --whitener-metadata data/donors/whitener_union_pruned.json
```

**Zero flags change. There is no recommended arm.**

The one-line mechanism: the tail block is doing **two** jobs, and they are priced on opposite
sides of the same ledger. As a call set it contributes 41.8% of our calls at exactly chance
precision, worth at most **+0.00089** to remove. As emitted dispersion it is the sole reason
`claim` sits at 0.24581 against the cap 0.24653, worth **+0.01852**. The second is 21× the
first. **`fid` is not the member this block should be optimised for, and `mse` already has
it.**

If anyone revisits this, the *only* survivor worth a second look is **armB′ (`--tail-spread
0.75`)** at −0.00759 — it is the single arm that does not clamp `score_mse` to zero. It is
still negative and it is still NO-GO, but it is the only point in the family with a
continuous, non-clamped price.

### Preregistered falsifiers, had this been a GO

Recorded so a future agent can reuse them rather than re-derive them:

1. `claim` on the built artifact, measured by `mech_chal_stats` + `jackknife_correction`,
   must land within 3% of the value in §5. If an arm's built `claim` exceeds 0.23268 the
   `mse` price is *not* the clamped −0.01852 and the arm must be re-priced.
2. `pds` must move by < 1e−5 and `nmae` by < 0.0005; anything larger means the exact-total
   repair broke and the arm is a profile change in disguise.
3. Challenge-controls `n_pred` probe must land within 5% of the h1-test-implied ratio
   (armA 0.519, armB 0.774, armC 0.591 × 4,424).
4. Board `fid` must **rise**. Given §4's five-point board series it almost certainly will
   not, and that is the single cheapest test of the whole thesis.

---

## 9. What I did NOT measure, and why

| not measured | why | can it move? |
| --- | --- | --- |
| `pds` and `nmae` per arm on the challenge | pinned by `emit_counts`'s exact-total repair; the pseudobulk is bit-identical across `(spread, φ)`. Confirmed on h1-test inside seed noise (Δpds ≤ 0.0012, Δnmae ≤ 0.00018, seed spread 0.00070 / 0.00018) and by v14's board `Δpds = 0.000000`. | **No** — rule 12's stated exemption |
| the **band's precision** | the challenge ships no reference truth, and h1-test has no band (§3.1). Only its size, direction and call rate are measurable. | Bounded: at ≤ 0.2% call rate it cannot move `fid` by more than ~1e−4 |
| `--min-cpm 5.0` on `A`, `X`, `pds`, `nmae`, `mse`-via-profile | changes the profile; needs a rebuild + `mech_chal_stats` on the new artifact. Deliberately not attempted rather than guessed. | **Yes** — which is exactly why it is a NO-GO here rather than a maybe |
| arm D (post-hoc drop) as a shippable arm | we submit **counts**, not calls; the call set is recomputed by the scorer. Reported as a bound only (+0.000886). | Not expressible |
| per-arm `VCC_ALPHA` re-optimisation | amplitude is a second variable and outside this mandate (161/162 own adjacent work). Priced both ways in §6 so the reader can see it does not rescue any arm. | Priced, not measured |
| `jac` raw on the challenge | rule 9 — `n_real` is in the union. Rebuilt from the board's own `n_real`/`n_pred` with the h1-test overlap ratio; the model, not a measurement. | Small either way (\|Δ\| ≤ 0.0008) |
| `keep_top` variation | one variable per run (rule 5); `keep_top` is a different axis from the tail policy and was not opened. | Untested |

### Integrity

* Ten known-answer gates, listed in §1, all passing; four of them (`n_pred` ×2, `reach` ×2)
  are bit-identical reproductions of 157 from a fresh emission.
* Every emission arm carries a re-seeded twin, and the twin is what made the `pred_lfc`
  invariance check readable rather than vacuous (§3.2) — the arm effect is *smaller* than the
  seed effect, which is only visible because the seed effect was measured.
* **Two of my own controls failed and are reported as failures, not quietly dropped**: the
  modal-ratio zeroed-gene detector (§3.1) and the `pred_lfc` invariance assertion (§3.2). The
  first was replaced with `verify_band.py`, whose failure mode is stated in advance and whose
  band ratio would read ≈1 if the band *were* zeroed; it reads 12.14.
* The `mse` numbers use the general form `1 + (s²A − 2sX + Kc − min(s²B, claim))/ΣD`, **not**
  `mech_claim_grid2.py`'s optimal-`s` shortcut, which silently re-tunes the amplitude and
  would have understated every arm's cost by ~4.6×.
* No submission, no `vcc prep`, no `cell-eval2`, no build. Peak RSS 8.8 GB (`realize164.py`,
  inherited running). No file outside this directory was modified.

## 10. Files

| file | what it produces |
| --- | --- |
| `blocks.py` / `blocks.log` | block census on the real challenge controls (predecessor) |
| `realize164.py` / `realize164.log` | h1-test DE tables for 7 arms + the per-target spread vector |
| `de164.*.npz` | those tables (85 × 10,779, real + predicted `p`/`padj`/`lfc`, `spread`) |
| `zero_lfc.py` / `zero_lfc.log` | §2 — the normalisation artefact, call rates by block, scored direction |
| `decomp.py` / `decomp.log` | §3.3 — six members, precision decomposition, arm D bound |
| `decomp2.py` / `decomp2.log` | §3.1, §3.3 — pooled vs per-perturbation, call-set anatomy |
| `verify_band.py` / `verify_band.log` | §3.1 — the between-target test that settles the band |
| `claim164.py` / `claim164.log` | §5 — `claim(policy)` on the real challenge panel |
| `price.py` / `price.log` | §4, §6 — reconciliation, the six-member table, the known-answer row |
