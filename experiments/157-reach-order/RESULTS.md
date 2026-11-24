# 157 — `reach` and the confidence ordering: is the channel really closed?

Date: **2026-09-02**. Measurement only. No `cell-eval2`, no build, no `vcc prep`, no
submission, no existing file modified. Everything is written inside this directory.

Decision, in one line: **the ordering channel is closed — NO-GO on a slot for reach
ordering.** The exp-094 inversion is gone, exp 154's closure survives on a better surface
with a better instrument, and no buildable key beats the shipped emitted order by enough to
matter. One unrelated finding came out of the work and is flagged at the end.

---

## 0. What changed, honestly

Three things, and only the third is a new lever:

1. **Exp 094's inversion headline was measured on the depth-confounded harness.** It is
   reproduced here to three decimals on `experiments/092-call-set-truncation/de_tables.npz`
   and does not survive the repair.
2. **Exp 154 was run on the one surface whose `reach` magnitudes are known not to transfer,
   against a phi=0 prediction, with `|pred_lfc|` rather than the emitted order as its
   baseline.** Redone on the reach-representative harness with a power-matching instrument
   that `reach` can actually be read through, its conclusion holds.
3. **The shipped tier spread policy puts a 74-79%-rate, 100%-DOWN, median-`|z|`-7 artefact
   on ~1,900 genes, and `reach` reads them first.** Measured in closed form on the real
   challenge controls, so this one is not harness-limited.

---

## 1. Surfaces, and which numbers on them are readable

| surface | file | perturbations | gated genes | median `n_real` |
| --- | --- | ---: | ---: | ---: |
| broken HCT116 (pre-repair) | `092/de_tables.npz` | 132 | 10,258 | 304 |
| repaired HCT116 | `095/de_thinned.npz` + this dir's `de_shipped.*` | 132 | 10,261 | 68 |
| **h1-test** | this dir's `de_h1test.*` (built here; none existed) | 85 | 10,779 | 89 |
| the challenge | — | 300 x 3 | ~9,900 | ~240 (a21) |

**AGENTS.md §11.1 rule 9 applies to every raw `reach` / `fid` / `jac` number below except
§5's.** `reach` divides by `n_real`, which is 68 / 89 / ~240 across these three.

**h1-test is validated on this member, and HCT116 is not.** Board raw `reach` against the
same two emission laws realised here at `f = 1`:

| emission law | board raw `reach` | h1-test (2 seeds) | HCT116 (2 seeds) |
| --- | ---: | ---: | ---: |
| flat `spread 0.25`, `phi 4` (a8) | **0.1890** | 0.2029 / 0.1922 | 0.1525 |
| tier8000 `tail 1.0`, `phi 4` (a15) | **0.1811** | 0.1828 / 0.1754 | 0.0659 / 0.0750 |

h1-test reproduces both board points to within 0.014 raw. HCT116 misses a15's by 0.11 —
**a factor of 2.6**. This is the standing requirement of AGENTS.md §11.1 rule 3, and only
h1-test meets it for `reach`.

---

## 2. Q1 — the inversion is GONE, and 094's headline was a harness artefact

`inversion.py`, on the reference's significant adjudicable pool. Positive Spearman means
*higher* predicted `p_adj` is *more* often sign-correct, i.e. our order is inverted.

| surface / emission | acc INSIDE our call set | acc OUTSIDE | paired `t` | inverted perts | pooled `rho` |
| --- | ---: | ---: | ---: | ---: | ---: |
| **094 as recorded** | 0.3446 | 0.8259 | −54.0 | 99.2% | **+0.5009** |
| broken HCT116, reproduced here | **0.3282** | **0.8335** | **−53.9** | **100.0%** | **+0.5007** |
| repaired HCT116, phi=0 (pre-a8) | 0.5307 | 0.6292 | −5.2 | 69.4% | +0.0494 |
| **repaired HCT116, SHIPPED a15** | **0.6329** | **0.4664** | — | — | **−0.1294** |
| repaired HCT116, a15 re-seeded | 0.6358 | 0.4664 | — | — | −0.1329 |
| **h1-test, SHIPPED a15** | **0.5437** | **0.5187** | — | — | **−0.0204** |
| h1-test, a15 re-seeded | 0.5450 | 0.5185 | — | — | −0.0213 |
| h1-test, phi=0 | 0.5323 | 0.5360 | — | — | −0.0004 |

**094's three headline numbers reproduce to 0.017, 0.008 and 0.0002 on the broken
harness, and on no other surface.** They are the depth confound documented in
`docs/HARNESS-DEPTH-CONFOUND.md`, not a property of our emitter: that file's reference DE
skews 90.8% negative, our calls skewed up, and the "inversion" is those two facts colliding.

**Under the shipped configuration the association is correctly signed on both repaired
surfaces.** `rho` is **−0.1294** (HCT116) and **−0.0204** (h1-test) against 094's +0.5009,
and the accuracy inside our call set now *exceeds* the accuracy outside it on HCT116
(0.6329 vs 0.4664) rather than being half of it. Do not cite 0.3446 / 0.8259 / +0.5009
again.

The effect size is small on h1-test (`rho` −0.02) — the shipped order is roughly
*uninformative*, not anti-informative. That is the honest state: the inversion was fixed,
and nothing replaced it.

---

## 3. Q2 — exp 154's closure survives, and the instrument it needed did not exist

### 3a. `debench.Surface(n_real_target=N)` is the wrong instrument for `reach`

It replaces every perturbation's BH set with the reference's top-`N` by p-value. That is
correct for `fid` and `jac`, whose artefact is the *size* of `n_real`. For `reach` it
destroys the statistic: h1-test's native median `n_real` is 89, so `N = 240` should be a
mild change, and instead the shipped arm's emitted reach falls **0.1828 → 0.0080**, 23x, for
a 2.7x pool change.

The cause is that `reach = mean_p(k*_p / n_real_p)` and `n_real_p` is heterogeneous. Reach
by `n_real` tercile on h1-test at native power, shipped emission:

| | small `n_real` | mid | large |
| --- | ---: | ---: | ---: |
| emitted order | **0.5011** | 0.0343 | 0.0013 |
| 154's best key | 0.5455 | 0.0568 | 0.0103 |

**98% of our reach comes from the small-`n_real` third of perturbations.** Flattening every
perturbation to a common large `N` deletes that stratum by construction, so the collapse
measures the flattening, not the power.

### 3b. Proportional power matching

`power.py` scales each perturbation's pool by a common factor instead:
`pool_p = the reference's top-round(f * n_real_p) by its own p-value`. `f = 1` reproduces the
native BH set to within p-value ties — the known-answer row required by AGENTS.md §11.1
rule 10's corollary, and it passes.

**h1-test, shipped a15 emission** (donors `k562_gwps` + `orion_hek293t`; `h1_2025` is
excluded — see §7):

| | f=0.5 | f=1.0 | f=1.5 | f=2.0 | **f=3.5** | f=5.0 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| median `n_real` | 44 | 89 | 134 | 179 | **314** | 449 |
| emitted order (incumbent) | 0.2381 | 0.1828 | 0.1507 | 0.1220 | **0.0536** | 0.0283 |
| `\|pred_lfc\|` | 0.2693 | 0.2217 | 0.1868 | 0.1383 | 0.0662 | 0.0495 |
| 154's best key (geo x concord) | 0.2610 | 0.2084 | 0.1955 | 0.1504 | 0.0631 | 0.0378 |
| ORACLE (correct signs first) | 0.6551 | 0.6455 | 0.6368 | 0.6221 | **0.5790** | 0.5742 |
| **capture % of oracle headroom** | **7.5** | **8.4** | **9.2** | **5.7** | **2.4** | **3.9** |
| overall value of the best key | +0.0059 | +0.0073 | +0.0084 | +0.0054 | **+0.0024** | +0.0040 |

The challenge's `n_real ~ 240` sits between `f = 2` and `f = 3.5` on this surface.

**Answer to the question as posed: the buildable-vs-oracle ratio does change with `N`, but
it changes DOWNWARD — 9.2% at `f = 1.5` to 2.4% at `f = 3.5`.** 154 reported 3.4% at its own
power against `|pred_lfc|`; the challenge regime on the reach-representative surface reads
2.4-5.7% against the *emitted* order, which is the harder and correct baseline. **154's
closure is not an artefact of pool size, and higher power makes the channel worse, not
better.**

The oracle is almost flat in `f` (0.655 → 0.574) because `oracle reach ~ accuracy / 0.9`:
ordering all correct signs first gives a pure prefix of length `accuracy * n_real`, which
clears the 0.90 bar at `accuracy/0.9 * n_real`. That is a pure function of per-gene sign
accuracy and carries no pool-size dependence at all — which is precisely why the gap 154
measured is a **model** gap.

### 3c. On the a8 flat emission there is no gap to capture at all

Same instrument, `a8_flat_s0.25_phi4`:

| surface | f=0.5 | f=1.0 | f=1.5 | f=2.0 | f=3.5 | f=5.0 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| h1-test capture % | 7.4 | 4.5 | 0.7 | 1.5 | 1.5 | 1.0 |
| HCT116 capture % | **−10.4** | **−8.3** | **−2.3** | **−3.4** | **−1.0** | **−1.5** |

**Under a8's emission law the emitted p-order already equals or beats every key 154
screened.** The keys are not ahead of the incumbent; on HCT116 they are behind it.

---

## 4. Q3a — the key families 154 never screened are all worse

`keys_extended.py`, repaired HCT116. Three build-time families that are entirely free and
had never been tried: the surface's own **control cells** (per-gene CPM, detection rate,
CV — the challenge hands us `context_{A,B,C}.h5ad`), the **gene-level marginal
responsiveness** of the donor panel, and products of both with `|pred_lfc|`.

| key | native BH | N=240 | overall vs incumbent @240 |
| --- | ---: | ---: | ---: |
| incumbent (emitted p-order) | 0.07702 | 0.01578 | — |
| `\|pred_lfc\|` (154 baseline) | 0.10082 | 0.02107 | +0.00100 |
| **geo x concord (154 best)** | **0.11462** | **0.03808** | **+0.00420** |
| `ctrl_cpm` | 0.05897 | 0.00966 | −0.00115 |
| `ctrl_detect` | 0.05452 | 0.00890 | −0.00129 |
| `−ctrl_cv` | 0.03484 | 0.00634 | −0.00178 |
| `\|lfc\| x log ctrl_cpm` | 0.11184 | 0.02823 | +0.00234 |
| `\|lfc\| x ctrl_detect` | 0.11921 | 0.02877 | +0.00245 |
| `\|lfc\| / ctrl_cv` | 0.10350 | 0.02138 | +0.00106 |
| donor marginal `\|delta\|` | 0.07456 | 0.01051 | −0.00099 |
| `\|lfc\| x donor marginal` | 0.11018 | 0.02335 | +0.00143 |
| geo x concord x `ctrl_detect` | 0.10875 | 0.03060 | +0.00279 |
| geo x concord x log cpm | 0.11439 | 0.02968 | +0.00262 |

**Nothing beats 154's winner.** Control-cell measurability and donor-panel responsiveness
are both real gene-level signals and neither carries per-gene *sign* information. 154's key
set was adequate; the key set was never the limitation.

---

## 5. Q3b — the one thing that does move the order, and it is challenge-regime

### 5a. The shipped tier spread block

a12/a14/a15 ship `--tail-spread 1.0 --keep-top 8000 --keep-spread 0.25 --overdispersion 4.0`:
~8,000 gated genes at `spread 0.25` and the remaining ~1,900 at `spread 1.0`, `phi 4`
throughout. `tail_sweep.py`, repaired HCT116, 132 perturbations, profile and `phi` held
fixed so **only the spread vector changes**:

| arm | `n_pred` | `reach` | Δ raw | Δ overall | top-10 purity | **top-10 drawn from the tail block** |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `a8_flat_0.25` | 2,867 | 0.15249 | — | — | 0.7729 | 0.0% |
| `a15_tier8000_tail1.0` | 4,555 | 0.06591 | −0.08658 | −0.01630 | 0.6266 | **79.7%** |
| `a15_tier8000_tail1.0_b` | 4,552 | 0.07504 | −0.07745 | −0.01458 | 0.6276 | **79.6%** |
| `tier8000_tail1.0_gatedonly` | 4,545 | 0.06786 | −0.08463 | −0.01594 | 0.6316 | 79.5% |
| `tier8000_tail0.50` | 3,020 | 0.08903 | −0.06346 | −0.01195 | 0.6721 | 59.3% |
| `tier9500_tail1.0` | 3,480 | 0.10888 | −0.04361 | −0.00821 | 0.6820 | 51.9% |

Monotone in the size and strength of the tail block, and the `gatedonly` row rules out the
per-cell library-variance channel: it is the tail genes' own p-values.

### 5b. Why — measured on the real challenge controls, so rule 9 does NOT apply

`zero_effect.py` extends 153 §15's exact-integer, exact-tie closed form to the shipped
policy. Zero predicted effect everywhere; the profile is the control pseudobulk; the only
input is `data/challenge/context_{A,B,C}.h5ad`.

| `spread`/`phi` | `Fbar` median | `\|z\|` median | significant at ZERO effect (A / B / C) | called UP |
| --- | ---: | ---: | ---: | ---: |
| `0.25` / `phi 0` (v11-era) | 0.5388 | 2.7 | 44.9% / 48.5% / 36.3% | 100% |
| `0.25` / `phi 4` (a8's law) | 0.5204 | 1.7 | **19.0% / 22.8% / 16.2%** | 88 / 86 / 81% |
| `1.00` / `phi 0` | 0.5151 | 1.1 | **7.7% / 10.2% / 8.2%** | 100% |
| **`1.00` / `phi 4` (a15's tail)** | **0.3918** | **7.4** | **78.9% / 73.9% / 76.2%** | **0.0%** |
| **a15 as shipped (mixture)** | — | — | **30.6% / 31.4% / 28.8%** | — |

The `0.25 / phi 4` row reproduces 153 §15's 19.0 / 21.6 / 16.0% — the known-answer row.

**`spread` and `phi` interact non-monotonically and the interaction flips the call
direction.** At `spread 1.0` the gamma-Poisson with `phi = 4` (gamma shape 0.25) sends most
of a low-expression gene's 400 cells to zero, so the emitted median falls below the
control's and the rank-sum reports a large, highly significant **DOWN** shift at zero
effect: median `|z| = 7.4` against the keep block's `1.7`. `reach` orders our significant
calls first by ascending `p_adj`, so **those ~1,900 genes are placed at the very front of
the ranking, by construction, and none of them carries any model signal.** That is exactly
the 79.7% figure in §5a.

Two corollaries worth recording:

* **The shipped call set is ~30% zero-effect artefact, not ~19%.** 153 §15's decomposition
  ("~19 points artefact, ~25 points effect-driven") was computed for a8's flat law and does
  not describe a12/a14/a15. Correct figure for the shipped build: **28.8-31.4%** of gated
  genes are called at zero effect, against an `n_pred` of ~4,400 out of ~9,900 gated.
* **`spread 1.0 / phi 0` is the cleanest null this family has ever produced** — 7.7-10.2%,
  better than 153 §16's varying-library estimate of 7.1% and reachable with a flag rather
  than a rebuild. Setting the tail block's `phi` to 0 while keeping its `spread` at 1.0
  takes the shipped mixture from 30.6/31.4/28.8% to **16.8/20.7/14.6%**.

### 5c. But the board does not corroborate the size, and says so directly

The one board test of this change is **confounded and points the other way**. a11 → a12
introduced the tier policy *together* with `VCC_ALPHA` 0.838 → 1.0108, and scaled `reach`
went **UP** `+0.004569` (`A12-BOARD-RESULT.md`). 094's amplitude slope (+0.013 raw reach per
doubling) accounts for +0.0035 raw = +0.0040 scaled of that on its own, so the tier
residual on the board is **≈ 0**, not the −0.016 overall HCT116 reads.

Cross-surface, at `f = 1`, tier minus flat:

| surface | Δ raw `reach` | seed spread |
| --- | ---: | ---: |
| repaired HCT116 | **−0.0820** | 0.009 |
| **h1-test** (validated, §1) | **−0.0185** | 0.008 / 0.011 |
| board, a8 vs a15 (confounded by 4 other changes) | −0.0079 | — |
| board, a11 vs a12 residual after the amplitude model | ≈ 0 | — |

**HCT116 overstates it by 4.4x against the surface that reproduces the board.** The effect
is real and its mechanism is measured on the challenge's own controls, but its board value
is between 0 and −0.008 raw `reach`, i.e. **0 to −0.0015 overall**.

---

## 6. Answers, and the recommendation

**Q1. Is the inversion still present?** No. Pooled `Spearman(pred p_adj, sign correct)` is
**−0.1294** (HCT116) and **−0.0204** (h1-test) under the shipped configuration, against
094's +0.5009; sign accuracy inside the call set is **0.6329 vs 0.4664** outside on HCT116
and **0.5437 vs 0.5187** on h1-test, both correctly ordered. 094's numbers reproduce at
+0.5007 / 0.3282 / 0.8335 **only on the depth-confounded harness**. What remains is not an
inversion but an absence of signal, plus a *volume* problem: **28.8-31.4% of gated genes are
called significant at zero predicted effect** under the shipped law, worse than a8's
16.2-22.8%.

**Q2. Does 154's closure survive?** **Yes, and it strengthens.** Its surface was reach-blind
and its instrument was wrong, and fixing both does not rescue the channel: on the validated
h1-test surface the best buildable key captures **8.4%** of the oracle headroom at native
power and **2.4%** at challenge-like power (`f = 3.5`), worth **+0.0024 overall** as a
counterfactual — and 095's realisation work already showed counterfactual key gains do not
survive the emitter (`+0.00583` counterfactual → `−0.00021` realised). On a8's emission law
the emitted order already **beats** every key on HCT116. The oracle is invariant in pool
size because `oracle reach ~ sign accuracy / 0.9`, which is the cleanest possible statement
that this is a model gap.

**Q3. Best deployable ordering, and what is it worth?** There isn't one. No key — including
three families 154 never screened — beats the shipped emitted order by more than
**+0.024 raw reach / +0.0046 overall** counterfactually at native power, falling to
**+0.013 / +0.0024** at challenge power, before realisation losses that historically take
such gains to zero. The only deployable change that moves the order materially is the
**emission law**, specifically the tier tail block, and its board-scale value is
**0 to +0.0015 overall** with `mse` unpriced.

### GO / NO-GO

**NO-GO. Do not spend a submission slot on reach ordering.** The ceiling is +0.0024 overall
at the challenge's power regime, it is a counterfactual that has never realised, and our gap
to the top ten on this member is +0.065 raw reach (+0.0122 overall) — 27x larger. Reach is
model-limited. This confirms 154 rather than overturning it.

### The one item that should not be dropped

`spread 1.0 / phi 4` is a **74-79%, 100%-DOWN, median-`|z|`-7 null violation** on ~1,900
shipped genes. That is an emission-law defect, not an ordering opportunity, and it is
measured on the challenge's own controls. It plausibly touches `fid`, `jac` and `mse` far
more than `reach`, and `spread 1.0 / phi 0` fixes the null at 7.7-10.2%. **This is outside
this experiment's mandate and is handed over, not recommended:** the tail block exists to
raise the emitted dispersion `mse`'s jackknife credit reads, so any change to it must be
priced on `mse` first (AGENTS.md §11.1 rule 12). `pds` and `nmae` are structurally invariant
across every arm in §5a — exact-total repair plus fixed per-cell `L` — and were not measured
for that reason.

---

## 7. Integrity, and one error caught

* **A donor leak was caught and fixed mid-run.** `data/donors/h1_2025` covers **85/85**
  h1-test targets, so the first h1-test key screen read `reach 0.404` for
  `geo x concordance` — pure oracle. All h1-test numbers above use
  `DONORS=k562_gwps,orion_hek293t`. HCT116 keeps `h1_2025` (12/132 coverage, a different
  cell line) exactly as 154 did, and excludes `orion_hct116`.
* **Known-answer rows, both passing.** `prod_phi0_s0.25` re-emitted here reads `n_pred`
  4,019 and reach 0.0765 against `de_thinned`'s 4,020 and 0.0770. `zero_effect.py`'s
  `0.25 / phi 4` row reads 19.0 / 22.8 / 16.2% against 153 §15's 19.0 / 21.6 / 16.0%.
* Every emission arm carries a re-seeded twin. Seed spreads: 0.009 (HCT116 a15), 0.008 and
  0.011 (h1-test). All reported effects except §5c's board reconciliation exceed them.
* `de_thinned.npz`'s prediction is `g0.65_q1.0_a1.08_spread0.25` — **`phi = 0`**, built
  2026-08-27 20:32, before a8 shipped over-dispersion. Every ordering statistic previously
  read off it, including all of 154's, describes the pre-a8 emitter.
* The h1-test reference DE table did not exist and is built here
  (`de_h1test.*.npz`, median `n_real` 89). `CURRENT_STATE.md` quotes 222 for h1-test; that
  figure is from a different build and is not reproduced by this one.

## 8. Files

| file | what it produces |
| --- | --- |
| `inversion.py` | Q1: in/out call-set sign accuracy, paired `t`, per-pert and pooled Spearman |
| `realize_shipped.py` | HCT116 DE tables for 4 emission laws through production `emit_counts` |
| `realize_h1test.py` | h1-test reference DE table + the same emission laws |
| `power.py` | proportional power matching; the §3b and §3c tables |
| `power_screen.py` | the uniform-`N` version, kept to document why it is wrong for `reach` |
| `keys_extended.py` | §4, the three unscreened key families |
| `arms.py` | per-arm reach / anatomy / inversion across surfaces |
| `tail_sweep.py` | §5a, the tail-block dose response |
| `zero_effect.py` | §5b, the challenge-regime zero-effect artefact |
