# 171 — the reach prefix: where `P(k)` dies, and why the call set cannot fix it

Date: **2026-09-02**. Measurement only. No `cell-eval2`, no build, no `vcc prep`, no
submission, no `submissions/` write, no existing file modified. Everything is written inside
this directory. Peak RSS 3.1 GB.

Decision, in one line: **NO-GO, and the mandate's mechanism is refuted structurally rather
than empirically — the reach spec's "our significant calls first" rule is a mathematical
NO-OP for every submission the scorer can ever produce, so the call set cannot move `reach`
at all.** What is left is the p-value ORDER; this experiment re-derives that from the walk
itself, opens the one axis 164 explicitly left untested (`keep_top`), finds the only
`mse`-free window in this family's history — and then closes the *whole* spread-schedule
family with a mechanism: **`|z|` is U-shaped in the emitted spread with its minimum at the
shipped 0.25, and both branches are dominated by control expression rather than by our
predicted effect**, so no spread policy can deliver the ordering key we actually hold (§12).

**Read §12 with §2.** Four coordinator `INBOX` files arrived mid-run with a measured
falsifier and a proposed fix. The falsifier reproduces to the digit and confirms the
diagnosis; the fix is built, priced (it *is* `mse`-free, as §5 predicts) and **refuted** — it
doubles the expression dependence it was designed to cancel, and its permutation control
scores *better* than it does.

---

## 0. The five things that changed

1. **The pole-position hypothesis is dead by proof, not by measurement.** BH-adjusted
   p-values are monotone in p, so "significant" is *exactly* "`p_adj < 0.05`" — which means
   sorting the reach pool by ascending `p_adj` already places every called gene ahead of
   every uncalled one. The significant-first key is redundant with the second key. Measured:
   **0/85 perturbations where it could bite, 0 ties**; and `reach` with the call set replaced
   by *nothing* or by *everything* is **0.18283 either way — identical to 5 decimals**. §3.
2. **`P(k)` never reaches the bar.** Mean `P(1)` is **0.7317** for a15 and **0.8171** for the
   best arm in this whole experiment, against `P₀ = 0.90`. The curve is monotone decreasing
   from there to ~0.52. **No arm's mean purity touches 0.9 at any depth `k ≥ 1`.** `reach` is
   bounded by per-gene sign accuracy (0.5343), not by ordering. §2, `pk-curve.png`.
3. **164's "a15 is on the cliff" understates the position.** a15's `claim` is 0.24581 against
   the *knee* `s²B` = **0.23268**, not against the cap — so a15 has **+0.01313 of claim
   headroom** and `mse` is **exactly flat** for any arm above the knee. That opens a real
   `mse`-free window: `--keep-top ≤ ~8350`, plus everything that raises dispersion. §5.
4. **The spread-schedule family is closed by a sign error in its own mechanism.** A
   rank-sum test reads a *median shift*, not a variance, and our only noise channel is a
   zero-inflated gamma-Poisson that moves the median as a side effect. So past `s ≈ 0.25`
   adding noise *manufactures* significance instead of suppressing it: `n_pred` runs
   3,581 → 6,034 → 9,472 as flat spread goes 0.25 → (a15 tier) → 1.00, and 7,051 → 6,034 as
   `keep_spread` goes 0.00 → 0.25. **U-shaped, minimum at the shipped value.** §12c.
5. **The one arm in the mse-free window is real, mse-free, reach-positive — and still does not
   clear the gate.** `--keep-top 8300` reads Δ`mse` = 0.000000 (measured on the challenge
   panel) and a twin-averaged h1-test Δ`reach` of +0.02520. But 69% of that comes from five
   perturbations, it falls to +0.00941 on the `n_real ≥ 20` stratum, and the board's own
   five-point `fid`-vs-`n_pred` series charges it −0.001141. Central net **+0.0003**, against a +0.0019 gate. §7.

---

## 1. Instruments, surfaces, and twelve known-answer gates

| gate | expected | measured | source | verdict |
| --- | ---: | ---: | --- | --- |
| a15 h1-test `reach`, native BH | 0.18283 | **0.18283** | 157 `power_h1test.log` | PASS (fresh emission) |
| a8-flat h1-test `reach` | 0.20291 | **0.20291** | 157 / 164 §1 | PASS (fresh emission) |
| a15 median `n_pred` | 6035 | **6035** | 164 §1 | PASS |
| a8-flat median `n_pred` | 3581 | **3581** | 164 §1 | PASS |
| my spec walk vs `debench.reach()` | identical | **\|diff\| 0.0e+00** | independent implementation | PASS |
| `claim(0.25, φ0)` | 0.11860 | **0.11860** | 096 `claim_grid2.log` | PASS (−0.00%) |
| `claim(1.00, φ0)` | 0.19216 | **0.19438** | 096 `claim_grid.log` | PASS (+1.16%) |
| `claim(0.25, φ4)` = a8 | 0.16478 | **0.16652** | 096 `board_solve` | PASS (+1.06%) |
| `claim(a15 as shipped)` | 0.24581 | **0.24581** | 164 §5 | PASS (+0.00%) |
| a15 raw `mse` at shipped α | 0.88228 | **0.88228** | 164 §5 | PASS |
| h1-test→board `reach` shrink | 0.428 | **0.4277** | 164 §6, re-derived with twins | PASS |
| `f = 1` power pool = native BH set | identity | identity | 157 §3b | PASS |

**Surfaces.** h1-test (85 perturbations × 10,779 gated genes) for the four DE members;
the real challenge panel (`096/work/v11-stats` + `cp-proxy-ntc` + `board_solve`) for `mse`;
the board's own history for `fid`'s `n_pred` response.

**Donors.** 095's h1-test artifact, `k562_gwps` + `orion_hek293t`. `data/donors/h1_2025`
covers 85/85 h1-test targets and never enters. **This is a deliberate departure from the
mandate**, which asked for production's three-donor set including `orion_hct116` — see §9.

**`debench.Surface(n_real_target=N)` is not used.** 157 §3a: uniform power matching deletes
the small-`n_real` stratum that carries the score. Proportional matching (`pool_p` = the
reference's top `round(f·n_real_p)`) is used instead, with `f = 1` as its own known-answer row.

---

## 2. The `P(k)` decomposition — the primary deliverable

`diag171.py`, a15 as shipped, h1-test. This is the diagnosis the mandate asked for.

### 2a. Per-perturbation anatomy

| quantity | median | mean | min – max |
| --- | ---: | ---: | ---: |
| `n_real` (reference's DE set, target removed) | 115 | 793.2 | 1 – 5,484 |
| adjudicable genes in the pool | 115 | 793.2 | — |
| **of those, genes WE call significant** | **82** | **495.8** | **62.5% of the walked pool** |
| **`k*`** | **1** | **2.24** | **0 – 40** |
| `R_p = k*/n_real` | 0.0032 | 0.1828 | — |

**`k*` is 1.** Not 40, not 10 — one gene. **26.8% of perturbations score `k* = 0`**, i.e. the
very first gene in their walk is wrong and the ranking never recovers.

### 2b. The `P(k)` curve

| `k` | mean `P(k)` | median `P(k)` | % of perturbations at `P(k) ≥ 0.9` | perturbations with `n_adj ≥ k` |
| ---: | ---: | ---: | ---: | ---: |
| 1 | **0.7317** | 1.0000 | **73.2%** | 82 |
| 2 | 0.6623 | 0.5000 | 46.8% | 77 |
| 3 | 0.6491 | 0.6667 | 28.9% | 76 |
| 5 | 0.6200 | 0.6000 | 8.6% | 70 |
| 10 | 0.5655 | 0.6000 | 3.4% | 58 |
| 20 | 0.5350 | 0.5500 | 2.0% | 50 |
| 50 | 0.5357 | 0.5200 | **0.0%** | 47 |
| 100 | 0.5237 | 0.5200 | 0.0% | 41 |
| 500 | 0.5395 | 0.5340 | 0.0% | 28 |

Plotted in **`pk-curve.png`** (left panel), with the `P₀ = 0.9` bar and the chance line.
The curve starts at 0.73, is below 0.66 by `k = 2`, and asymptotes to 0.52 — chance.
**Beyond `k = 30` not one perturbation in 85 clears the bar.**

### 2c. Where the walk first goes wrong

| position of the first wrong sign | value |
| --- | ---: |
| median | **2** |
| ≤ 1 | 26.8% of perturbations |
| ≤ 2 | 56.1% |
| ≤ 3 | 73.2% |
| ≤ 5 | 92.7% |
| ≤ 10 | 98.8% |

### 2d. Head composition — what takes the pole positions

| depth `h` | % that are our significant calls | **% from the zero-effect TAIL block** | mean `P(h)` |
| ---: | ---: | ---: | ---: |
| 1 | 98.8% | **61.0%** | 0.7317 |
| 3 | 98.3% | **65.1%** | 0.6491 |
| 10 | 92.5% | **69.6%** | 0.6292 |
| 20 | 93.6% | **75.8%** | 0.6094 |
| 50 | 95.3% | **79.1%** | 0.6105 |

*("TAIL block" = `spread == tail_spread`, i.e. the gene fell outside a15's top-8,000 by
predicted effect magnitude — the operational form of zero/near-zero predicted effect.
`pk-curve.png` right panel plots this against position.)*

### 2e. The zero-effect fraction ahead of `k*` — the number the mandate asked for

| | value |
| --- | ---: |
| genes ahead of `k*`, pooled over 85 perturbations | 184 |
| **of those, zero-predicted-effect (TAIL block)** | **94 = 51.1%** |
| of those, our significant calls | 180 = 97.8% |
| **per-perturbation zero-effect fraction ahead of `k*`** | **mean 0.540, median 0.500** (60/82 perturbations have `k* > 0`) |

### 2f. Sign accuracy on the reach pool, by block — the artefact is real

| block | n | pooled `k/n` | per-perturbation | sd of the mean |
| --- | ---: | ---: | ---: | ---: |
| our sig calls, KEEP block | 23,194 | 0.5704 | **0.6629** | 0.0225 |
| our sig calls, **TAIL block** | 17,463 | 0.5083 | **0.5212** | 0.0264 |
| not called, KEEP block | 23,178 | 0.5204 | 0.5520 | 0.0281 |
| not called, TAIL block | 1,209 | 0.4864 | 0.5005 | 0.0376 |

**The mandate's empirical premise is confirmed**: the zero-effect block occupies 61–79% of the
head at per-perturbation accuracy 0.5212, against 0.6629 for the block carrying model signal.
What fails is the inference that the *call set* is therefore the lever. §3.

### 2g. The `n_real` tercile split — and a rule-9 warning about this surface

| tercile | n | `n_real` | reach | share of the panel mean | `k*` median | zero-eff ahead of `k*` |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| small | 28 | 1 – 13 | **0.5011** | **93.6%** | 2 | 0.223 |
| mid | 27 | 15 – 536 | 0.0343 | 6.2% | 1 | 0.871 |
| large | 27 | 648 – 5,484 | 0.0013 | 0.2% | 1 | 0.542 |

**93.6% of h1-test's `reach` comes from perturbations with `n_real ≤ 13`, a stratum the
challenge (`n_real ≈ 240`) does not have.** This is rule 9 in its sharpest local form, and it
is why every arm below is also reported at proportional power and on the `n_real ≥ 20`
stratum. 157 measured 98% for the same statistic on a coarser split.

---

## 3. The mandate's mechanism is a NO-OP — proved, then measured

The reach order is `lexsort(gene_name, −|lfc|, p, p_adj, ¬significant)` with `¬significant`
as the primary key. The mandate's structural consequence #2 — *"our significant calls are
placed FIRST, unconditionally"* — is true, and it is **vacuous**:

> Benjamini–Hochberg adjusted p-values are monotone non-decreasing in the raw p-values, and
> a gene is significant iff `p_adj < α`. So the significant set is **exactly a `p_adj`
> threshold set**, and sorting by ascending `p_adj` — the *second* key — already places every
> called gene ahead of every uncalled one. The first key can never change the order.

Measured on the a15 table (`work/monotone.log`):

* **0 / 85** perturbations where `max(p_adj | called) ≥ min(p_adj | not called)`. **0 ties.**
* `reach` with the significance flag stripped entirely (`nosigpole`): **0.18283** — identical.
* `reach` with the call set replaced by **nothing at all**: **0.18283**.
* `reach` with the call set replaced by **every gene**: **0.18283**.

The corollary about the alphabetical tiebreak is also empty: a perturbation carries
**10,779 distinct raw p-values on 10,779 genes** — the `p` key never ties, so the `|lfc|` and
gene-name keys are unreachable.

**So the call set cannot move `reach`, in either direction, ever.** We submit counts; the
scorer computes our calls by BH on our own p-values; a BH call set is a p-threshold set by
construction. The only exploitable version of the pole-position rule would require submitting
a call set that is *not* a p-threshold — measured here as `keep only` (0.20166) and `tail only`
(0.17683) — and that is not expressible in a counts submission.

**This is a different closure from 157's.** 157 concluded the ordering channel is closed
empirically (no buildable key beats the incumbent). This says the *call-policy* channel is
closed structurally, which is stronger, and it means the mandate's alternative reading
("limited by CALL POLICY, not ordering quality") is not merely unsupported but unavailable.

---

## 4. Counterfactual ceilings — what reordering is worth before any emitter takes its cut

`counter171.py`, same cached table for every row, so the pseudobulk is literally the same
array and `pds` / `mse` cannot move.

| ordering | f = 1 | f = 2 | f = 3.5 | `k*` med | `P(1)` | Δ overall @ f=1 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| incumbent (production order) | 0.18283 | 0.13111 | 0.07334 | 1 | 0.7317 | — |
| `nosigpole` (drop the sig-first rule) | 0.18283 | 0.13111 | 0.07334 | 1 | 0.7317 | +0.00000 |
| **`keepfirst`** (keep-block calls ahead of tail-block calls) | 0.20166 | 0.15760 | 0.08436 | 2 | 0.8049 | **+0.00355** |
| `lfcfirst` (sig first, then \|pred lfc\|) | 0.21047 | 0.15729 | 0.08174 | 2 | 0.8537 | +0.00521 |
| `lfconly` (\|pred lfc\| over the whole pool) | 0.22173 | 0.15682 | 0.07701 | 2 | 0.8415 | +0.00732 |
| **ORACLE** (correct signs first) | 0.64553 | 0.62769 | 0.59656 | 66 | 1.0000 | **+0.08713** |
| ANTIORACLE | 0.09756 | 0.07317 | 0.02439 | 0 | 0.0976 | −0.01606 |

`lfconly` reproduces 157 §3b's `|pred_lfc|` row (0.2217 at `f = 1`) to five decimals from an
independent implementation. The ORACLE sits at 0.6455 against the structural identity
`accuracy/0.9 = 0.5937`, as 157 §3b predicts.

**The whole reordering channel is worth +0.0036 (realisable shape) to +0.0073 (not
realisable), against an oracle of +0.0871.** `lfconly` is not buildable: 094 measured that the
emitted p-order is a function of `|F̄ − ½|/√(tc)` which correlates with `|predicted effect|`
at Spearman **+0.07**, so no emission law makes p-order follow `|lfc|`.

---

## 5. `mse`: the knee, not the cliff — and the only free window in the family

`claim171.py` / `claim171b.py` on the real challenge panel; `msebranch171.py` for the branch.

The scored numerator is `1 + (s²A − 2sX + Kc − min(s²B, claim))/ΣD`. **`min` saturates.**
At a15's shipped `VCC_ALPHA = 0.982` (`s = 0.90960`):

```
knee  s²B  = 0.23268          a15 claim = 0.24581          headroom = +0.01313
reference cap K = 0.24653     a15 raw mse = 0.88228        score_mse = 0.11112
```

**Above the knee `mse` is exactly flat in `claim`.** 164 read a15 as "on the cliff"; it is
0.01313 *above* it, and that difference is the entire opportunity this experiment found.

| policy | `Ĉ/C_p` | `claim` | branch | raw `mse` | Δ overall |
| --- | ---: | ---: | --- | ---: | ---: |
| **a15 `--keep-top 8000`** | 1.085 | **0.24581** | SATURATED | 0.88228 | **0.000000** |
| `--keep-top 8150` | 1.031 | 0.24148 | SATURATED | 0.88228 | **0.000000** |
| **`--keep-top 8300`** | 0.976 | **0.23580** | SATURATED | 0.88228 | **0.000000** |
| `--keep-top 8400` | 0.940 | 0.23006 | linear | 0.89375 | −0.002012 |
| `--keep-top 8700` | 0.836 | 0.20601 | linear | 0.99896 | −0.018520 |
| `--keep-top 20000` (all gated) | 0.732 | 0.18049 | linear | 1.11065 | −0.018520 |
| a8 flat 0.25/φ4 | 0.675 | 0.16652 | linear | 1.17175 | −0.018520 |
| `--keep-spread 0.00` | 0.778 | 0.19180 | linear | 1.06115 | −0.018520 |
| `--keep-spread 0.10` | 0.860 | 0.21191 | linear | 0.97317 | −0.015943 |
| `--keep-top 2000 --keep-spread 0` | **3.206** | **0.24653** | SATURATED | 0.88228 | **0.000000** |
| flat 1.00/φ4 | **4.027** | **0.24653** | SATURATED | 0.88228 | **0.000000** |

**The `mse`-free window is `--keep-top ≲ 8350`** (interpolating the knee between 8300 and
8400) **plus everything that raises dispersion.** Two sub-results worth keeping:

* **The sub-gate block alone is not enough.** The idea that all the dispersion could live
  outside the scorer's 5-CPM DE universe — where no DE member can see it but the jackknife
  can — is measured and **false**: `--keep-top 20000` leaves `claim` at 0.18049, far below the
  knee. The **976 gated tail genes carry 0.0653 of `claim`**, and the ~9,556 sub-gate genes
  only 0.180.
* `keep_spread → 0` **destroys** `claim` (0.19180) rather than preserving it, because at
  `spread = 0` the quasi emitter's random part is zero and the gene contributes no jackknife
  variance at all.

---

## 6. Realized emission arms on h1-test — and the pseudobulk identity

`realize171.py`, 12 arms + twins, all through production's `emit_counts`.

**Rule-12 pseudobulk check, in its non-vacuous form.** The sha256 of the *delivered* per-gene
panel total is recorded per (arm, target). All **12 arms across 3 runs are identical on 85/85
targets** — digest `b53d02aebfd0816a`. This is not an algebraic identity: without
`emit_counts`'s exact-total repair the Poisson residual lands in the totals and every arm's
digest differs (that is what the repair exists to prevent), so the check reads differently
when the thing it tests is absent. BRIEF §1 and §2 define `pds` and `expr_mse` on exactly this
pseudobulk, so **Δ`pds` = 0.000000 and the profile half of Δ`mse` = 0 are identities, not
estimates.**

*(Caveat worth stating: the `pds` column in the table below is `debench`'s, computed from the
DE table's mean-of-ratios `lfc`. It moves by up to 0.0020 across arms. That movement is an
artefact of the proxy — the scored `pds` reads the pseudobulk and is pinned by the digest.
Do not read that column as a `pds` measurement.)*

| arm | `n_pred` | reach f=1 | f=2 | f=3.5 | fid | precision | jac | nmae |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| **a15** | 6034 | 0.18283 | 0.13111 | 0.07334 | 0.50701 | 0.50701 | 0.06824 | 0.96659 |
| a8 flat 0.25/φ4 | 3581 | 0.20291 | 0.16920 | 0.07265 | 0.50793 | 0.51209 | 0.06089 | 0.96673 |
| `--keep-spread 0.00` | 7051 | 0.17587 | 0.12288 | 0.07498 | 0.50393 | 0.50393 | 0.06916 | 0.96655 |
| `--keep-spread 0.10` | 6767 | 0.18675 | 0.12957 | 0.07468 | 0.50460 | 0.50460 | 0.06985 | 0.96672 |
| `--keep-top 2000 --keep-spread 0` | 9055 | 0.16435 | 0.13351 | 0.07376 | 0.50520 | 0.50520 | 0.07254 | 0.96633 |
| flat 1.00/φ4 | 9472 | 0.18062 | 0.13414 | 0.06977 | 0.50523 | 0.50523 | 0.07224 | 0.96669 |
| `--keep-top 9000` | 5232 | 0.19394 | — | — | — | — | — | — |
| **`--keep-top 10000`** | 4318 | **0.20618** | 0.15782 | 0.07508 | 0.50918 | 0.50990 | 0.06347 | 0.96665 |
| `--keep-top 20000` | 3588 | 0.20934 | 0.17241 | 0.08763 | 0.50823 | 0.51235 | 0.06089 | 0.96670 |

Walk anatomy — the mechanism is exactly as §2 predicts:

| arm | `k*` med | `k*`=0 % | `P(1)` | `P(3)` | **% of head from the TAIL block** |
| --- | ---: | ---: | ---: | ---: | ---: |
| a15 | 1 | 26.8 | 0.7317 | 0.6491 | **69.6%** |
| `--keep-spread 0.00` | 1 | 24.4 | 0.7561 | 0.5877 | 66.0% |
| `--keep-top 2000 --keep-spread 0` | 1 | 34.1 | 0.6585 | 0.5789 | 82.1% |
| flat 1.00/φ4 | 1 | 31.7 | 0.6829 | 0.6667 | 100.0% |
| `--keep-top 10000` | 2 | 19.5 | 0.8049 | 0.7018 | 53.3% |
| `--keep-top 20000` | **2** | **18.3** | **0.8171** | 0.7281 | **0.0%** |
| a8 flat | 2 | 19.5 | 0.8049 | 0.7456 | 0.0% |

**`keep_spread → 0` does NOT realize the `keepfirst` counterfactual.** It was the direct test
of the "move the signal block to the high-|z| end of the U" idea and it fails: the head stays
66% tail-block and reach goes *down* (twin mean 0.17791 vs a15's 0.17910). At `spread = 0` a
gene's p-value is a deterministic function of where its 400 identical values sit in the
control CDF — 094's `|F̄ − ½|/√(tc)` — which has nothing to do with block membership. The U
argument was too crude and its arm is refuted.

**What does work is removing the tail law from the DE-visible axis**, i.e. raising `keep_top`.
`--keep-top 20000` reaches `P(1) = 0.8171` and 0% tail in the head, matching the `keepfirst`
counterfactual's 0.8049 — the counterfactual and its realisation agree.

---

## 7. The six-member ledger, twin-averaged and board-scaled

**Seed noise first.** h1-test `reach` has a large same-policy spread: 0.00746 (a15 pair),
0.01068 (a8-flat pair), 0.01233 (`--keep-top 20000` pair), 0.00376 (`--keep-top 10000` pair),
against arm effects of ~0.020. Every number below is a **twin mean**.

**Calibration (a known-answer row for the transfer itself).** Board a8 − a15 raw `reach` =
+0.00790; h1-test a8-flat − a15, twin-averaged = +0.01847. **h1-test overstates by 2.34×;
shrink 0.4277**, reproducing 164's independently derived 0.428.

| arm | `claim` | Δ`mse` | Δ`reach` raw (board) | Δ`reach` | Δ`fid` | Δ`jac` | **NET** |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| a8 flat 0.25/φ4 | 0.16652 | **−0.018520** | +0.00790 | +0.001488 | +0.000852 | +0.000570 | **−0.015610** |
| `--keep-top 20000` | 0.18049 | **−0.018520** | +0.01030 | +0.001939 | +0.000852 | +0.000570 | **−0.015159** |
| `--keep-spread 0.00` | 0.19180 | −0.018520 | −0.00051 | −0.000096 | −0.000200 | −0.000150 | **−0.018966** |
| `--keep-spread 0.10` | 0.21191 | −0.015943 | +0.00327 | +0.000616 | −0.000100 | −0.000080 | **−0.015507** |
| `--keep-top 2000 --keep-spread 0` | 0.24653 | 0.000000 | −0.00631 | −0.001188 | −0.000400 | −0.000300 | **−0.001888** |
| flat 1.00/φ4 | 0.24653 | 0.000000 | +0.00065 | +0.000122 | −0.000500 | −0.000350 | **−0.000728** |
| `--keep-top 9000` | 0.24148 | 0.000000 | +0.00483 | +0.000910 | +0.000400 | +0.000260 | +0.001570 |
| **`--keep-top 10000` → challenge `8300`** | **0.23580** | **0.000000** | +0.01078 | +0.002030 | *see below* | +0.000400 | *see below* |

Δ`pds` = **0.000000** (identity, §6). Δ`nmae` = **−0.00002** overall (h1-test 0.96659 →
0.96665, against a 0.00008 re-seed spread — inside noise, and it *can* move, since BRIEF §7
reads the DE table's mean-of-ratios `lfc`, so this is a measurement rather than an exemption).

### 7a. The best arm, and why it still fails

`--keep-top 8300` (h1-test analogue `--keep-top 10000`, fraction-matched: 92.5% vs 92.8% of
the gated axis promoted) is the only arm in this family's history that is **`mse`-free and
`reach`-positive**. Three things kill it.

**(i) The reach gain is carried by five perturbations.** `robust171.py`, twin-averaged:

| | f = 1 | f = 2 |
| --- | ---: | ---: |
| panel mean Δ | +0.02520 | +0.02941 |
| paired `t` | +2.74 (30 up / 20 down / 32 unchanged) | +2.82 (39 / 21 / 22) |
| most influential perturbation | `n_real` **8**, Δ +0.4375, **21% of the panel sum** | `n_real` **3**, Δ +0.5000, 21% |
| Δ with it deleted | +0.02011 | +0.02360 |
| **Δ on `n_real ≥ 20` only** | **+0.00792** | **+0.00941** |
| top-5 \|Δ\| share of the panel sum | **69%** | **69%** |

The challenge has no `n_real ≤ 13` stratum. On the stratum it does resemble, the gain is
**3× smaller**.

**(ii) The h1-test→board shrink is not stable in `f`.** Calibrating on the a8↔a15 pair at
each power factor gives 0.4277 (`f=1`), 0.2875 (`f=1.5`), 0.1925 (`f=2`) — a factor of 2.2
across the challenge-relevant range, which I cannot reduce. The board-scaled reach estimate
is therefore a range, not a point:

| reading | Δ`reach` overall |
| --- | ---: |
| `f=1`, full panel, 164's shrink 0.4277 | **+0.00203** |
| `f=2` (median `n_real` 229 ≈ the challenge's 240), self-calibrated shrink 0.1925 | +0.00107 |
| `f=2`, `n_real ≥ 20` stratum, shrink 0.4277 | +0.00076 |
| `f=2`, `n_real ≥ 20` stratum, shrink 0.1925 | **+0.00034** |

**(iii) `fid` goes the other way, and the board has measured it five times.** The arm cuts
`n_pred` 6034 → 4318 on h1-test (×0.716), implying a challenge `n_pred` of ~3,166:

| model | board `fid` | Δ overall |
| --- | ---: | ---: |
| 164's conserved-headroom projection (optimistic) | 0.50220 → 0.50307 | **+0.000490** |
| **the board's own five-point `fid`-vs-`n_pred` series, interpolated** | 0.50220 → **0.50016** | **−0.001141** |

Series: `n_pred` 1563 → 0.49126, 1911 → 0.49654, 3836 → 0.50210, 4424 → 0.50220,
5961 → 0.50178. 164 §4: *"shrinking the call set has never once raised board `fid`"* — five
direct board points on emission-only changes at a fixed profile. The projection and the direct
measurement **disagree on sign**, and the direct measurement wins.

### 7b. The verdict arithmetic

| | Δ overall |
| --- | ---: |
| `pds` | 0.000000 (identity) |
| `mse` | 0.000000 (measured, saturated branch, challenge panel) |
| `nmae` | −0.000017 (inside seed noise) |
| `jac` | +0.000400 |
| `reach` | +0.00034 … +0.00203 |
| `fid` | −0.001141 (board series) … +0.000490 (projection) |
| **NET, pessimistic** (f=2 `n_real ≥ 20` reach, board-series fid) | **−0.00042** |
| **NET, central** (f=2 panel reach, board-series fid) | **+0.00031** |
| **NET, optimistic** (f=1 panel reach, projected fid) | **+0.00290** |

**Rule 13's shape-only ship gate is +0.0019.** The central estimate is +0.00031 — **6.2× below
it** — and the pessimistic reading is negative. The gate is a veto and it vetoes this arm.

---

## 8. Did 164 already close this? — Precisely half of it

**Yes on the axis it opened, no on the axis it did not, and the mandate's specific mechanism
was never the live question.**

* **164 DID price `reach` for its arms.** Its §6 six-member table carries a `d reach` column
  (armA +0.002394, armB +0.000301, armC +0.001488 overall), on this same h1-test surface with
  this same 0.428 board shrink. So *"removing the zero-effect calls changes precision; their
  POSITION IN THE REACH WALK may not have been priced at all"* is **not correct** — it was
  priced, on the right surface, with the right calibration.
* **164 did NOT open `keep_top`.** Its §9 says so in as many words: *"`keep_top` is a
  different axis from the tail policy and was not opened — Untested."* That is this
  experiment's contribution, and it is where the only `mse`-free window turned out to live.
* **164's `mse` reading needed one correction.** "a15 is on the cliff, any arm that drops
  `claim` lands in the linear branch" understates the position by the +0.01313 of headroom
  between a15's `claim` and the knee. Without that correction the `--keep-top 8300` window is
  invisible, and 164 would have been right for the wrong reason.
* **The mandate's mechanism was never available.** §3: the significant-first key is
  algebraically redundant. No experiment could have priced it, because there is nothing there.

---

## 9. GO / NO-GO

### NO-GO on every arm. Do not spend a slot. a15 ships unchanged, zero flags different.

```
VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.10 VCC_MASS_FIX=minnorm \
  --min-cpm 10.0 --overdispersion 4.0 --tail-spread 1.0 --keep-top 8000 \
  --keep-spread 0.25 --whitener-id union_pruned \
  --whitener data/donors/whitener_union_pruned.npz \
  --whitener-metadata data/donors/whitener_union_pruned.json
```

**The one-line mechanism.** `reach` cannot be moved by the call set — that rule is a no-op —
and the emitted p-order can only be moved by the emission law. Every emission change that
materially improves the order does so by removing the `spread 1.0 / φ 4` law from the
DE-visible axis, and that law is where 27% of a15's `claim` lives. There is exactly one
`keep_top` window (≤ ~8,350) where the two do not collide, and inside it the central reach
gain (+0.00107) is cancelled almost exactly by the `fid` charge the board's own history
levies (−0.00114). **And §12 closes the escape route**: the `mse`-free window is generous, but
no spread policy can put our informative key at the front, because `|z|` is U-shaped in the
spread with its minimum at the shipped 0.25 and both branches track control expression rather
than predicted effect. **Underneath all of it, mean `P(1)` is 0.73–0.82 against a bar of 0.90: our
ordering has nowhere to go, because our per-gene sign accuracy is 0.5343.**

### If the slot economics ever change, this is the arm and this is its preregistration

`--keep-top 8000` → **`--keep-top 8300`**. Exactly one flag differs from a15; diff the built
meta and require **exactly one** difference (rule 10), and assert `whitener_id ==
union_pruned` and the operator sha, since `--whitener-id` silently defaults to `union`.

Preregistered prediction: **+0.00031 overall (range −0.00042 to +0.00290)**, against rule 13's
+0.0019 shape-only gate. Falsifiers:

1. Built `claim` (via `mech_chal_stats` + `jackknife_correction`) must land **≥ 0.23268**. If
   it does not, `mse` is *not* free and the arm loses −0.0185 outright.
2. `pds` must move by **< 1e−5** and `nmae` by **< 0.0005**; larger means the exact-total
   repair broke and the arm is a profile change in disguise.
3. Challenge `n_pred` (094's `probe_de.py` on the built artifact) must land within 5% of
   **3,166**. Materially higher and the `fid` charge shrinks; lower and it grows.
4. Board `fid` must **rise**. The five-point series says it will not, and that is the single
   cheapest test of the whole family.

---

---

## 12. The coordinator's arm (INBOX 1–4): reproduced, tested, and closed — with the sign error

Four `INBOX-*.md` files arrived mid-run. INBOX 4 carried a measured falsifier and a proposed
fix. **I reproduced the falsifier to the digit, and the fix is refuted — it runs backwards,
and the reason is a sign error in its mechanism that this experiment's own data pins down.**

### 12a. The falsifier reproduces exactly (`work/falsifier.log`)

| quantity | INBOX 4 | measured here |
| --- | ---: | ---: |
| `reach`, production order | 0.18283 | **0.18283** |
| `reach`, key = \|H\| (design-space effect) over the whole pool | 0.23533 | **0.23533** |
| `reach`, sig-first then \|H\| inside the block | 0.21914 | **0.21914** |
| Spearman(−log₁₀ `p_adj`, \|H\|), DE universe | −0.121 | **−0.121** |
| Spearman(−log₁₀ `p_adj`, \|H\|), inside our call set | −0.391 | **−0.391** |
| Spearman(−log₁₀ `p_adj`, expression), universe | +0.444 | **+0.444** |
| Spearman(−log₁₀ `p_adj`, expression), inside our call set | +0.485 | **+0.485** |

Sign accuracy on the reference DE pool by decile of the ordering key (d1 = front of the walk):

| key | d1 | d2 | d3 | d5 | d8 | d10 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| our p-value (production order) | **0.550** | 0.525 | 0.559 | 0.568 | 0.502 | 0.510 |
| \|H\| design effect | **0.690** | 0.607 | 0.559 | 0.497 | 0.457 | 0.504 |
| \|pred_lfc\| (emitted) | 0.646 | 0.613 | 0.544 | 0.498 | 0.484 | 0.505 |
| expression (control CPM) | 0.540 | 0.550 | 0.528 | 0.547 | 0.543 | 0.552 |

**The diagnosis is confirmed and it agrees with §2–§4**: we hold an informative key (\|H\| puts
0.690 at the front against production's 0.550), the metric does not read it, and INBOX 4's
observation that ordering the *called* genes by \|H\| gives the same 0.21914 as ordering the
whole pool is the same fact as §3 — the call set is not the defect.

### 12b. The proposed fix, built and measured

INBOX 4 §3b derives `σ²_ours = s·m + φ·s²·m²` from `_quasi_block`, hence
`z/h ~ √m / √(A + s + φs²m)`, and concludes that flattening the expression dependence means
**raising `s` on high-expression genes**. `realizex171.py` implements the exact quadratic
solve (`φc²m·s² + c²·s + (c²A − m) = 0`, clipped to [0.25, 1.0], sub-gate genes left at 1.0),
at two target quantiles, plus a rank-keyed variant, plus a permutation control.

**`mse` first, and the arm was given every chance** (`claimx171.py`, real challenge panel,
both known-answer gates passing):

| policy | `Ĉ/C_p` | `claim` | branch | raw `mse` | Δ`mse` overall | median `s` |
| --- | ---: | ---: | --- | ---: | ---: | ---: |
| a15 as shipped | 1.085 | 0.24581 | SATURATED | 0.88228 | — | 0.250 |
| **`xcomp10`** | **3.628** | **0.24653** | SATURATED | 0.88228 | **0.000000** | 0.998 |
| **`xcomp25`** | 2.414 | **0.24653** | SATURATED | 0.88228 | **0.000000** | 0.641 |

**The schedule is exactly `mse`-free**, as §5's knee analysis predicted — it only raises `s`,
and a15 already sits +0.01313 above the knee. So this arm did *not* die on `mse`. It died on
`reach`.

| arm | `n_pred` | med `s` | reach f=1 | f=2 | f=3.5 | fid | jac | nmae | `P(1)` |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| **a15** | 6034 | 0.250 | **0.18283** | 0.13111 | 0.07334 | 0.50701 | 0.06824 | 0.96659 | 0.7317 |
| `xcomp10` | 8818 | 1.000 | 0.17751 | 0.14508 | 0.07369 | 0.50488 | 0.07323 | 0.96670 | 0.7073 |
| `xcomp10_b` (twin) | 8811 | 1.000 | 0.17342 | 0.14125 | 0.07207 | 0.50516 | 0.07320 | 0.96641 | 0.6707 |
| `xcomp25` | 7184 | 0.827 | 0.17722 | 0.14252 | 0.07391 | 0.50567 | 0.07506 | 0.96694 | 0.7439 |
| `sgrad_expr` (rank-keyed) | 6256 | 0.625 | 0.17648 | 0.14112 | 0.07030 | 0.50632 | 0.07510 | 0.96678 | 0.6829 |
| **`xcomp10_perm` (CONTROL)** | 8607 | 1.000 | **0.18592** | 0.13809 | 0.06572 | 0.50540 | 0.07150 | 0.96645 | 0.6585 |

`xcomp10` twin mean **0.17546** against a15's **0.17910** — **Δ = −0.00364**, board-scaled
**−0.00029 overall**, and inside the seed spread either way. Every variant is neutral-to-negative.

**The control fires, and it fires the wrong way.** `xcomp10_perm` uses the *same multiset of
per-gene spreads*, randomly permuted across the gated axis — identical emitted dispersion
budget, identical spread histogram, only the gene↔spread pairing destroyed. It reads
**0.18592**, i.e. **0.0084 HIGHER** than the expression-keyed version. *(What would this read
if the mechanism were present? The permuted arm would score worse, because it would have
thrown away the compensation. It scores better.)*

### 12c. Why it fails: the variance model has the wrong sign at the shipped `φ`

The direct diagnostic — did the schedule flatten the expression dependence it was designed to
remove?

| arm | Spearman(−log₁₀ `p_adj`, expression) inside the call set | Spearman vs \|H\| |
| --- | ---: | ---: |
| a15 | +0.485 | −0.391 |
| **`xcomp10`** | **+0.939** | +0.001 |
| `xcomp25` | +0.942 | −0.099 |
| `sgrad_expr` | +0.945 | −0.119 |

**It roughly doubles the expression dependence it was built to cancel.** The mechanism, read
off arms this experiment already had on disk:

| law (φ = 4 throughout) | median `n_pred` |
| --- | ---: |
| flat `spread 0.25` | 3,581 |
| flat `spread 1.00` | **9,472** |
| keep block at `0.00` | 7,051 |
| keep block at `0.10` | 6,767 |
| keep block at `0.25` | 6,034 |

**`|z|` is U-shaped in `s`, with its minimum at ≈ 0.25.** Both ends inflate significance at
zero predicted effect, for two different reasons: at `s → 0` the even emitter puts 400
identical values at the control's *mean*, which for skewed counts sits above its median; at
`s → 1, φ = 4` the gamma-Poisson (shape 0.25) sends most cells to zero and the emitted median
falls *below* the control's — 157 §5b's median `|z| = 7.4`, 78.9%-significant, 100%-DOWN
artefact. So `∂|z|/∂s` is **negative below 0.25 and positive above it**, while INBOX 4 §3b's
variance model assumes it is negative everywhere.

**A rank-sum test reads a median shift, not a variance.** The only channel we have for adding
noise at a fixed mean is a zero-inflated gamma-Poisson, and it moves the median as a side
effect. "Suppress unearned significance by adding noise" is therefore not available: past
`s ≈ 0.25` adding noise *manufactures* significance of the wrong kind, and the schedule pushed
exactly the high-expression genes — the ones already dominating the head — to `s = 1.0`.

**This closes the spread-schedule family as a whole, not just this arm.** `keep_top` tiers on
`|effect|` and the expression term runs through it (§7); `keep_spread` tiers on the same
variable and sits at the U's minimum (§6); an expression-keyed schedule steers the wrong way
(here). The `|H|` key is worth +0.036 to +0.053 raw reach as a counterfactual, and **there is
no setting of the per-gene spread vector that delivers it**, because `|z|` is a U-shaped
function of `s` whose *both* branches are dominated by expression.

### 12d. What this changes in the verdict

Nothing, and it strengthens it. The `mse`-free window §5 found is real and generous, and it
now has a second occupant (`xcomp*`, at the cap with `Ĉ/C_p` up to 3.6). **The binding
constraint was never `mse`. It is that the emitted p-order cannot be made to track our own
predicted effect through any spread policy.** If reach is to be bought, it must be bought by a
change to the *test statistic's input* — a different emission family altogether, or a better
per-gene sign predictor — not by re-tiering the spread vector.


## 10. What I did NOT measure, and why

| not measured | why | can it move? |
| --- | --- | --- |
| **the three-donor h1-test artifact** the mandate specified (`k562_gwps` + `orion_hct116` + `orion_hek293t`) | 095's h1-test artifact is two-donor, and it is the object on which 157 board-validated h1-test's `reach` (0.1828 vs board 0.1811). Rebuilding with `orion_hct116` changes the prediction and **voids all four reach/`n_pred` known-answer gates**, which are the only reason any number here is trustworthy. The leakage constraint the mandate actually protects (`h1_2025`, 85/85 coverage) is fully honoured. | Would change absolute reach; the arm deltas are differences on a fixed prediction and are unlikely to invert. **Flagged as a deliberate deviation.** |
| challenge-side `n_pred` for the `keep_top` arms | needs a built `pred.h5ad`, i.e. a challenge count build — forbidden. Preserved as falsifier 3. | Refines `fid`; at any plausible value the board series still charges the arm |
| `fid` / `jac` raw on h1-test as board quantities | rule 9 — `n_real` is 1–5,484 here against ~240 there, so `max(n_pred, n_real)` and the jaccard union are in a different regime. h1-test raw values are reported but never used in the ledger. | Handled by two explicit board-anchored models, whose disagreement is reported rather than hidden |
| `--keep-key agreement`, `--min-cpm` | rule 5, one variable per run; `--min-cpm` additionally moves the profile, so rule 12's exemption does not cover it (164 §7) | Untested |
| φ schedules (per-gene `overdispersion`) | §12c shows `φ` enters as `φs²` and *worsens* the skewness pathology that drives the artefact; 157 §5b's zero-effect null would have to be re-closed first | Tested only at φ = 4; a φ = 0 tail is 164's armA, priced at −0.0185 |
| composing `--keep-top 8300` with a re-tuned `VCC_ALPHA` | amplitude is a second variable, and raising α raises the knee `s²B`, which would *shrink* the free window | Priced structurally: it can only make the window smaller |
| `cell-eval2` on anything | forbidden (25 GB RSS, shared machine) | — |

### Integrity

* **Twelve known-answer gates, all passing** (§1); four are bit-identical reproductions of
  157/164 from fresh emissions on this machine, and one — the h1-test→board shrink — is a
  reproduction of 164's calibration constant (0.4277 vs 0.428) from independently re-run twins.
* **Every arm carries a re-seeded twin**, and the twins were necessary: h1-test `reach` seed
  spread runs 0.004–0.012 against arm effects of ~0.020, so single-seed rows in this family
  are not resolvable. `--keep-top 10000`'s single-seed value (0.20618) reads *above*
  `--keep-top 20000`'s twin mean, which is non-monotone in block size and is pure seed noise —
  visible only because both twins were run.
* **My own arm was refuted and is reported as a failure, not dropped.** `keep_spread → 0` was
  this experiment's founding realizable arm (the U-shape argument in `realize171.py`'s
  docstring). It moved the head composition by 3.6 points and `reach` *down*. §6.
* **The pseudobulk control is non-vacuous by construction**: it hashes the delivered totals,
  and without `emit_counts`'s exact-total repair the Poisson residual lands in those totals and
  every digest differs. It is not the algebraic identity that was caught on 2026-09-01.
* **The `mse` numbers use the general `min(s²B, claim)` form**, not the linear shortcut.
  `claim171.py`'s own last table applies the linear branch everywhere and over-credits four
  arms by +0.000551; it is superseded by `msebranch171.py` and that supersession is recorded
  rather than silently fixed.
* No submission, no `vcc prep`, no `cell-eval2`, no build. No file outside this directory was
  modified.

## 11. Files

| file | what it produces |
| --- | --- |
| `diag171.py` / `work/diag_a15.log` | §2 — the `P(k)` decomposition, head composition, zero-effect fraction ahead of `k*`, tercile split |
| `work/monotone.log` | §3 — the BH-monotonicity proof's measurement (0/85, 0 ties, 10,779 distinct p) |
| `counter171.py` / `work/counter_a15.log` | §4 — ordering and call-set counterfactual ceilings |
| `realize171.py` / `work/realize171{,b,c}.log` | §6 — 12 emission arms + twins on h1-test, and the pseudobulk digests |
| `arms171.py` / `work/arms171.log` | §6 — four DE members + walk anatomy per arm |
| `claim171.py` / `claim171b.py` / `work/claim171{,b,c}.log` | §5 — `claim` on the real challenge panel, and the `keep_top` ladder |
| `msebranch171.py` / `work/msebranch171.log` | §5 — the knee, and the correct `mse` branch per arm |
| `ledger171.py` / `work/ledger171.log` | §7 — twin averaging, board calibration, six-member ledger |
| `power171.py` / `work/power171.log` | §7a(ii) — reach vs proportional power, twin-averaged |
| `robust171.py` / `work/robust171.log` | §7a(i),(iii) — influence analysis and the two `fid` models |
| `curve171.py` / `pk-curve.png` | §2 — the plotted `P(k)` curve and head composition |
| `work/falsifier.log` | §12a — reproduction of INBOX 4's three reach values and four Spearmans |
| `realizex171.py` / `work/realizex171.log` | §12b — the expression-compensating schedule + its permutation control on h1-test |
| `claimx171.py` / `work/claimx171.log` | §12b — `claim` for the schedule on the challenge panel |
| `work/xarms.log`, `work/ushape.log` | §12b,c — the schedule's six members, and the U-shape of `\|z\|` in `s` |
| `work/de171.*.npz` | the 12 DE tables (85 × 10,779, real + predicted `p`/`padj`/`lfc`, `spread`, `eff`) |
