# 230 — Is rule 13's `reach` toll a STEP or does it SCALE with the profile change?

**Verdict: neither, because the question was mis-posed. Rule 13's `−0.0146` is not a step and not
a slope — it is MIS-KEYED. Its two calibration points, a21 and a22, deliver per-gene panel totals
IDENTICAL to the incumbent's on 900/900 blocks, so by the table's own digest test they sit in the
row whose toll is `0`. The row the `−0.0146` is charged to contains ZERO board points. Hamilton
does not belong in it either, and should be priced at toll `0` — which leaves 208's rider at
`+0.002203` against a `+0.0019` floor, alive by 1.16× and inside the instrument's own noise.**

*No slot proposed. Nothing built. No `vcc`, no `cell-eval2`, no production change. Peak RSS
1.67 GB (`emit230.py`), everything else ≤ 0.41 GB.*

---

## 0. Known-answer rows — 27 of them, all PASS

Nothing downstream is believable without these (`phase0_230.py`, then a gate at the head of every
later script). Full log: `work/run230.log`.

| # | row | got | want |
|---|---|---:|---:|
| KA1 | scaling map, a8 raw `reach` 0.1890 → scaled | `0.124212` | `0.124349` |
| KA2 | toll `== mean(a21, a22)` `reach` | `−0.014569` | `−0.014569` |
| KA2b | toll on overall (`/6`) | `−0.002428` | `−0.002428` |
| KA3 | 214's shipped-renderer factor | `1.8407` | `1.84` |
| KA4 | 208's rule-13 composite (208's renderer) | `+0.001197` | `+0.001197` |
| KA4b | 219's composite (shipped renderer) | `+0.002204` | `+0.002203` |
| KA4c | 208's six-member sum, as measured | `+0.002899` | `+0.002899` |
| KA4d | 219's six-member, corrected | `+0.003904` | `+0.003909` |
| KA4e | composite with the `−0.00243` toll charged | `−0.000226` | `−0.000227` |
| KA5 | per-gene panel-total digest `(A, ABCD1)` | `bb6befc01caf7493` | `bb6befc01caf7493` |
| KA6 | a21/a22/rc_* board Δ`pds` = Δ`mse` = 0 exactly | `0.0 × 8` | `0.0 × 8` |
| — | 196's reconstruction: `m=1` totals `==` artifact totals | **`1.000000` on 900/900** | `1.0` |
| — | off-gate `floor(x) == T` | `1.000000` on 900/900 | `1.0` |
| — | **208's challenge tilt, `floor`** | **`−0.012511`** | `−0.01234` |
| — | **208's challenge tilt, Hamilton** | **`+0.000339`** | `+0.00039` |
| KA-A | `rc_callset` totals `==` a15's (208 gate 14) | `900/900` | `900/900` |
| KA-B | `cd4f034` totals `!=` `rc_callset`'s (210 f3a) | `0/900` | `0/900` |
| KA-C | `cd4f034` digest `(A, ABCD1)` (214 recorded) | `2cdbc2f0c36ba28b` | `2cdbc2f0c36ba28b` |
| KA-D | 222's `effects-off == effects-pre` | `True` | `True` |
| KA-E | **222's effect matrices → `cd4f034`'s own delivered totals** | **`1.00000000`, 900/900 blocks, gene for gene** | `1.0` |
| KA-F | `emit230`: every emission law kept the totals | `24/24` | `24/24` |

`rc_hybrid`'s digest result reproduces as KA5 + KA-A: `rc_callset` (which inherits `rc_hybrid`'s
emitter) is byte-identical to a15 on all 900 blocks, and block `(A, ABCD1)` hashes to the value
191/202/208/210 all recorded.

---

## 1. The measurement that decides it

### 1.1 Rule 13's two calibration points have a profile-change magnitude of EXACTLY ZERO

Three independent lines, all measured:

**(a) The renderer's arithmetic.** `emit_counts(mode="resample_mix")` — a21's and a22's law —
computes `totals = floor(n * L * profile_p)`, the *same* expression `quasi` uses, then calls
`_repair_totals(block, totals, rng)`. `L = int(min(median(lib_sizes), max_counts))` is a function
of the picked control libraries only. The emission law never enters the totals.
(`src/vcc2026/donor_transfer.py:496–520`.)

**(b) Cells emitted through production code, on real challenge blocks** (`emit230.py`, 6 blocks ×
4 laws, real context-A control cells, one shared seed):

| law | digest `==` quasi | genes differing | max\|ΔT\| | relL1 | **emitted lib CV** | frac resampled |
|---|---:|---:|---:|---:|---:|---:|
| `quasi` (a15) | 6/6 | 0.0 | 0 | 0 | 0.0187 | 0.00 |
| `phi = 0` (a8's channel) | 6/6 | 0.0 | 0 | 0 | 0.0043 | 0.00 |
| **`x = 0.14` (a21)** | **6/6** | **0.0** | **0** | **0** | **0.3821** | 0.86 |
| **the complement (a22)** | **6/6** | **0.0** | **0** | **0** | **0.1050** | 0.14 |

**(c) The board itself.** a21 and a22 moved `pds` and `mse` by `0.00000000` — sixteen significant
digits, `docs/LOG.md` 2026-09-01. That identity is only possible if the pseudobulk, hence the
per-gene panel totals, is unchanged.

### 1.2 The magnitude of every scored arm, and the arm the answer gates

`relL1 = Σ|ΔT| / ΣT` and `TV = ½Σ|p_a − p_b|` on delivered per-gene panel totals, averaged over
all 900 challenge blocks. Sources: artifacts (`cd4f034`), production emitter (a8/a21/a22), closed
form on the renderer (a23, Hamilton), closed form on 222's own effect matrices (`h1union`, with
KA-E as the check).

| arm | vs | digest | frac moved | **relL1** | **TV** | lib CV | **board Δ`reach`** |
|---|---|---:|---:|---:|---:|---:|---:|
| a8 | a2-prune | 900/900 | 0 | `0` | `0` | 0.0043 | **+0.021398** |
| a10 | a8 | 900/900 | 0 | `0` | `0` | 0.0187 | −0.000754 |
| **a21** | a15 | **900/900** | **0** | **`0`** | **`0`** | 0.3821 | **−0.018482** |
| **a22** | a15 | **900/900** | **0** | **`0`** | **`0`** | 0.1050 | **−0.010656** |
| a23 | a15 | 0/900 | 0.290 | `1.500e−03` | `7.506e−04` | 0.0187 | −0.007933 |
| rc_hybrid | a15 | 900/900 | 0 | `0` | `0` | 0.0187 | **+0.054940** |
| rc_callset | rc_hybrid | 900/900 | 0 | `0` | `0` | 0.0187 | +0.022590 |
| cd4f034 | rc_callset | 0/900 | 0.427 | `7.452e−03` | `3.726e−03` | 0.0187 | **+0.026402** |
| h1union | cd4f034 | 0/900 | 0.213 | `2.434e−03` | `1.218e−03` | 0.0187 | +0.003639 |
| **ham (208 A)** | cd4f034 | 0/900 | 0.438 | **`8.779e−04`** | `5.768e−04` | 0.0187 | **?** |

`h1union`'s Δ`reach` is `0.22293924` (leaderboard `2026-09-04T0524Z`) − `0.219300` (210's board
result) = **+0.003639**, both scaled, both read off our own scored rows. No board row was
inverted anywhere in this experiment, so the `score_mse`-clamp trap does not apply.

---

## 2. The decisive finding — rule 13's charged row is EMPTY

The table's first row has **two conjunctive conditions**: the arm must be an **emission-law**
change *and* its digest must differ. Applied to every scored arm:

| arm | emission law? | digest differs? | rule-13 row | toll charged | **realised** |
|---|---|---|---|---:|---:|
| a8, a10, **a21**, **a22**, rc_hybrid, rc_callset | yes | **no** | emission, profile-PRESERVING | `0` | +0.0214, −0.0008, **−0.0185**, **−0.0107**, +0.0549, +0.0226 |
| a23, cd4f034, h1union | no | yes | profile-class | `0` | −0.0079, +0.0264, +0.0036 |
| — | **yes** | **yes** | **emission, profile-CHANGING** | **`−0.0146`** | **n = 0** |

- **`emission, profile-CHANGING`: n = 0.** No arm this project has ever scored is in the row the
  `−0.0146` is charged to.
- **`emission, profile-PRESERVING`: n = 6, mean +0.011506, range [−0.018482, +0.054940].** This is
  where a21 and a22 actually live — the row the table assigns `0`. The constant is the mean of two
  points drawn from a row whose own mean is **positive**.
- **`profile-class`: n = 3, mean +0.007369, range [−0.007933, +0.026402].** Also positive, and
  `−0.0146` lies outside its entire observed range.

---

## 3. STEP vs LINEAR vs SATURATING — the fit, and why it decides nothing

Fit on all 9 board points, three magnitude axes, three functional forms (`fit230.py`):

| axis | law | R² | AIC | slope | toll at ham | LOO RMSE |
|---|---|---:|---:|---:|---:|---:|
| relL1 | STEP | 0.0080 | −64.91 | −0.0041 | −0.004137 | 0.027134 |
| relL1 | **LINEAR** | **0.0290** | −65.10 | +1.5900 | +0.001396 | 0.031385 |
| relL1 | SATURATE | 0.0000 | −64.83 | −0.0003 | −0.000149 | 0.027661 |
| TV | STEP / LINEAR / SATURATE | 0.0080 / 0.0290 / 0.0000 | | | −0.004137 / +0.001834 / −0.000180 | 0.027 / 0.031 / 0.028 |
| frac_moved | STEP / LINEAR / SATURATE | 0.0080 / 0.0001 / 0.0024 | | | −0.004137 / +0.000656 / −0.002634 | 0.027 / 0.029 / 0.028 |

**The evidence for each is nil.** Best R² is `0.029`. The sd of board Δ`reach` over the 9 arms is
`0.023164`; the best LOO RMSE of any law is `0.027134` — **1.17× the sd, i.e. every fitted law
predicts a held-out arm worse than predicting the mean does.** AIC separates the three forms by
less than 0.3, which is nothing.

The reason is structural and not a sample-size problem: **six of the nine arms sit at magnitude
exactly 0, and their Δ`reach` spans 0.0734 scaled (sd 0.027).** An axis on which the calibration
points are all at the origin and the outcome varies by three times the whole toll cannot identify
a slope, a step height, or a saturation scale.

### Falsifier 1 — out-of-sample retrodiction. Every law fails.

| held out | law | predicted | realised | error |
|---|---|---:|---:|---:|
| a21 | STEP | +0.017504 | −0.018482 | **+0.035986** |
| a21 | LINEAR | +0.012575 | −0.018482 | +0.031057 |
| a22 | STEP | +0.015938 | −0.010656 | +0.026594 |
| a22 | LINEAR | +0.011262 | −0.010656 | +0.021918 |
| a23 | STEP | +0.015020 | −0.007933 | +0.022953 |
| cd4f034 | STEP | −0.002147 | +0.026402 | −0.028549 |
| cd4f034 | LINEAR | −0.030804 | +0.026402 | **−0.057206** |

Errors run 1.4× to 3.9× the entire `−0.0146`, with the wrong sign every time. **A law that cannot
retrodict a held-out board point is not a law**, and none of these can.

### The one fit that does work — and it has the wrong sign

Fit *within* the digest-differing class on its own three scored members:

`Δreach = −0.012968 + 5.3585 · relL1`, **R² = 0.9639**

The slope is **positive**: inside the class rule 13 taxes, the bigger the profile change, the
*better* `reach` did. This is n = 3 with 2 parameters (1 residual df) and it is confounded —
`cd4f034` and `h1union` are direction arms that improved the underlying prediction, so their
`reach` gain is not a property of their profile-change size. **Do not use this slope.** It is
reported because it is the only magnitude fit with any explanatory power at all, and its sign is
the opposite of the one rule 13 assumes.

---

## 4. Falsifier 2 — separating the call-set effect from the profile-change effect

The mandate asked whether the `−0.0146` might be a call-set toll misattributed to the profile.
**It is separable, and the separation is clean, because one of the two terms is exactly zero.**

- **The profile term is 0.** a21's and a22's profile change is `0` — measured three ways (§1.1).
  So **100 % of the `−0.0146` is attributable to what else they changed**, with no residual to
  argue about. This is not "cannot separate"; it is separated at an exactly-zero coefficient.
- **What else they changed was the cell channel *and* the call set together.** a21 took challenge
  `n_pred` 4,420 → 1,141 with recall 0.669 → 0.271, by resampling 86 % of the gene axis off the
  parametric law onto real control cells (emitted library CV 0.019 → **0.382**).
- **But the call set alone is not the carrier.** `rc_callset` changed the call set by *more* on the
  jaccard (0.2072, `n_pred` 4,480 → 2,599) at zero cell-channel change and realised **+0.022590**.
  So a large call-set change is worth `+0.023` or `−0.018` depending on whether the new key is
  better or worse than the old one. **Size does not carry the sign; quality does.**
- **The cell channel is the only axis on which the calibration points are extreme.** A linear law
  on emitted library CV over all 9 arms gives `Δreach = +0.017109 − 0.1041·CV`, **R² = 0.2996** —
  10× the best profile-magnitude law, with the correct sign. It is still not a law (LOO RMSE
  0.033840, worse than the mean), and with only 2 of 9 arms off the baseline CV it is largely a
  re-description of a21 and a22. **Reported as the better-keyed axis, not as a fitted toll.**

**This is the bigger finding the mandate anticipated, and it lands harder than expected: the
`−0.0146` is not a profile toll at all. But it does not re-open every arm rule 13 has charged,
because — see §2 — it has never actually been charged to a scored arm. Its only live effect is
the veto it puts on Hamilton.**

---

## 5. Hamilton is not in this family at all

Every other arm in the digest-differing class changes the **intended** profile. Hamilton changes
only its **integer realisation** — and it moves it *toward* the intended one:

| rule | relL1 vs the INTENDED `x = n·L·p` | 208's size tilt (log2) |
|---|---:|---:|
| `floor` (a15 / the incumbent) | `8.778e−04` | **−0.012511** |
| **Hamilton (208 rider A)** | **`4.597e−04`** | **+0.000339** |

**Hamilton's delivered profile is 1.91× CLOSER to the intended one and its size tilt is 37×
smaller.** `max|ΔT| = 1` on every gene, every move UP (208 falsifier 4, reproduced). Its
"profile change" is a repair of the incumbent's own rounding error, not a new perturbation — and
under rule 13's first column it is not an emission-law change either, so it fails the first of the
row's two conjunctive conditions.

---

## 6. Phase 2 — 208's Hamilton rider, priced under every law

Rule-13 composite on the shipped renderer = **+0.002203** (219 §2.2, reproduced as KA4b). Floor
+0.0019. **These prices reuse the same board points the laws were fit on, except the two marked
otherwise — stated per the mandate's third falsifier.**

| toll law | toll scaled | on overall | composite | vs floor | verdict |
|---|---:|---:|---:|---:|---|
| rule 13 as written (STEP, −0.0146) | −0.014600 | −0.002433 | −0.000230 | 0.00 | VETO |
| most punitive in-class point (a23's realised) | −0.007933 | −0.001322 | +0.000881 | 0.46 | VETO |
| digest-differing **class mean** | +0.007369 | +0.001228 | +0.003431 | 1.81 | SHIP |
| within-class LINEAR at ham's relL1 *(extrapolated below the fitted range)* | −0.008264 | −0.001377 | +0.000826 | 0.43 | VETO |
| STEP fit on all 9 (R² 0.008) | −0.004137 | −0.000689 | +0.001514 | 0.80 | VETO |
| LINEAR fit on all 9 (R² 0.029) | +0.001396 | +0.000233 | +0.002436 | 1.28 | SHIP |
| **cell-channel law — the axis a21/a22 moved** | **0** | **0** | **+0.002203** | **1.16** | **SHIP** |

**Confidence interval.** 95 % on the class-mean toll: `+0.007369 ± 0.019768` scaled =
`−0.002066 … +0.004523` on overall. The composite under that interval runs
**`+0.000137 … +0.006726`** against a floor of `+0.0019`. **The interval straddles the floor.**

For scale: the *unexplained* spread of board Δ`reach` among digest-identical arms is sd `0.027011`
scaled = **`0.004502` on overall — 2.4× the entire floor**. No toll term of any shape can be
estimated to the precision this gate needs.

**Hamilton's price, stated honestly:**

- The toll no longer kills it. `−0.0146` is charged on a mis-keyed test, its calibration points
  are in the zero row, and its own row is empty. The correct charge on the evidence is **0**.
- At toll 0 the arm reads **+0.002203 against +0.0019 — 1.16×.**
- That 0.0003 margin is roughly **1/15th** of the unpriced noise (0.0045 on overall). Rule 13's own
  §7.4 — *"stop shipping small arms"* — applies with full force. **This is a live candidate, not a
  cleared one.**

---

## 7. What this does and does not re-open

- **Rule 13's `reach` toll should be re-keyed, not re-scaled.** The digest is the wrong test for
  the wrong quantity: it fires on `cd4f034` (+0.026) and not on a21 (−0.018).
- **No previously-charged arm is re-opened,** because §2 shows the charged row has never contained
  a scored arm. The toll's only realised effect on this project has been the veto on 208's rider.
- **A profile change is not, by itself, a `reach` cost.** Three profile-class board points
  (a23, `h1union`, `cd4f034`) average **+0.007369**, and 169's independent h1-test measurement of a
  profile-reshaping arm gave `+0.000605` overall with its entire 95 % interval above the toll.
- **`reach` remains not locally priceable** (162: 2/6 signs). Nothing here changes that. What
  changes is that the *toll* is no longer a defensible substitute for measuring it.

## 8. Caveats

1. **`h1union`'s magnitude is computed, not read off its artifact.** Its totals were derived from
   222's own effect matrices; the check is KA-E, which reproduces `cd4f034`'s delivered totals from
   the same path on 900/900 blocks gene-for-gene. The two arms share `n·L` because `L` depends only
   on the picked control libraries and 222 proved the default path byte-identical to `cd4f034`.
2. **a23's magnitude is closed-form**, from `mass_consistent_effect` re-expressed on
   `e = clip(expm1(b+h), 0)`; the on-gate profile carries 196's exact U(0,1) posterior, so a23's
   relL1 has a resolution floor of order one count per gene.
3. **The within-class fit is n = 3 with 2 parameters and is confounded by arm quality.** Flagged in
   §3; not used for the verdict.
4. **Nothing here re-measures a Class-A arm on the board**, builds a package, or proposes a slot.

## 9. Pinned inputs (rule 11: a path and a sha, never a description)

229 is editing `scripts/make_submission_pm_cd4.py` live. Per `INBOX-BUILDER-WARNING.md` a byte copy
was taken and the sha recorded; **nothing in this experiment imports the live builder at runtime** —
`profile_from_effect` and `mass_consistent_effect` are reimplemented in `fit230.py` /
`magnitude230.py` and validated against the artifacts rather than against the source text.

| file | sha256 | role |
|---|---|---|
| `work/make_submission_pm_cd4.PINNED.py` | `a8785cab8a6511e9…4181ec8` | the builder as read, byte copy |
| `src/vcc2026/donor_transfer.py` | `9ffca7ecac591954…865c7dfc` | `emit_counts`, imported live by `emit230.py`; **not** touched by 229 |

- `totals = np.floor(n * L * profile_p)` occurs **4 times** in `donor_transfer.py` — once per
  emitter mode. Every mode computes the same totals from the same profile. That is §1.1(a) as a
  grep rather than a reading.
- `MIN_CPM = 5.0` is the module default; production invokes `--min-cpm 10` (visible in 229's live
  build command), and 196's challenge gate is `cpm_pool > 10.0`. This experiment uses 10, and
  **KA-E confirms it empirically** — at any other gate the effect matrices would not reproduce
  `cd4f034`'s delivered totals on 900/900 blocks.

**Operating constraints (`INBOX-MEMORY.md`):** every job run strictly serially, none backgrounded;
true peak RSS `/usr/bin/time -l` measured per script (max **1.67 GB**, `emit230.py`); no sample size
reduced and no falsifier dropped. At close, swap free was 1,308 MB with a 229 production build
live, so no further job was launched.

## 10. Reproduce

```bash
cd experiments/230-toll-magnitude
uv run python phase0_230.py       # 27 known-answer rows          -> work/phase0_230.json
uv run python magnitude230.py     # 900-block profile magnitudes  -> work/magnitude230.json   7s, 0.38 GB
uv run python emit230.py          # production emitter, real cells-> work/emit230.json       40s, 1.67 GB
uv run python fit230.py           # the table + STEP/LINEAR/SAT   -> work/fit230.json         4s, 0.99 GB
uv run python phase2_230.py       # Hamilton's price              -> work/phase2_230.json     2s, 0.41 GB
```

Full transcript: `work/run230.log`.
