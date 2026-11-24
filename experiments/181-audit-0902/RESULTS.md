# 181 — ADVERSARIAL AUDIT of everything closed on 2026-09-02

Date: **2026-09-03**. Read-only outside this directory. **No build, no `vcc`, no `cell-eval2`,
no submission, nothing written to `submissions/`.** Every number below is either recomputed here
from artifacts already on disk or transcribed from a primary source with the path given. Peak RSS
6.9 GiB (`axis181.py`, memory-mapped Orion `delta.npy` pages counted resident).

**Headline: one claim is CONTRADICTED by a board point the closing experiments never ran, and it
is the day's largest positive result.** 177's "the ordering defect is solved by the emission law"
and 179's package prediction rest on h1-test's response to the resample emission law. **The board
has already scored that exact law — `a21`, 2026-09-01 — and h1-test gets 2 of its 4 DE members
wrong in SIGN, including `reach`, the member the whole claim is about.** 177 cites `a22` (the
wrong-way arm) and never retrodicts `a21` (the right-way arm), which lost more.

Six claims came out **SOLID** under independent recomputation, three of them strengthened
(the ρ-cliff/knee algebra is now board-validated to ±1%; the long-gene axis is now a −0.71
correlation instead of a gene-name eyeball; 176's CMap NO-GO is understated, not overstated).

---

## 0. What was recomputed, and the known-answer rows that gate it

Every script in this directory reproduces a published number before it prints anything new.

| # | known-answer row | source | measured here |
|---|---|---|---|
| 1 | six scaled members average to `score_avg` on every board entry | metric spec | **0 / 581 fail** |
| 2 | our board entry `0.1733597079`, rank 25 of 581 | CURRENT_STATE | **exact, rank 25** |
| 3 | per-member ranks 5 / 23 / 68 / 176 / 169 / 192 | CURRENT_STATE | **5 / 23 / 68 / 176 / 169 / 192** |
| 4 | `reach(production order)` = 0.18283 | 157 / 171 / 172 / 177 | **0.18283** |
| 5 | `reach(key=\|H\|)` = 0.23533 | INBOX-4 / 171 / 172 / 177 | **0.23533** |
| 6 | `reach(sig-first, then \|H\|)` = 0.21914 | INBOX-4 / 177 | **0.21914** |
| 7 | `reach(ORACLE)` = 0.64553 | 171 / 172 | **0.64553** |
| 8 | median `n_pred` 6034, call-set sign accuracy 0.53435 | 169 / 177 | **6034 / 0.53435** |
| 9 | 164's call-set split: 297,869 / 214,265 calls, precision 0.5122 / 0.4999 / 0.5070, paired +0.0122 | 164 §3.3 | **all exact** |
| 10 | 175's per-donor sign accuracy .5670 / .5375 / .5308 / .5425 / .5697 on n = 39,945 | 175 | **all exact** |
| 11 | 175's PC variance shares 5.6 / 8.6 / 19.7 %, cos(HEK,HCT) 0.872, vs K562 0.23–0.31 | 175 addendum | **5.64 / 8.60 / 19.65 %, 0.872, 0.234 / 0.307** |
| 12 | 180's per-context knee alphas 0.98234 / 1.01291 / 1.04786 | 180 §5 | **0.98234 / 1.01291 / 1.04786** |
| 13 | a23 raw `mse` 0.872341 through 096's identity | 180 §6 / 170 §4 | **0.872339** |
| 14 | 172 §7.3 power table: production reach 0.18283 / 0.10953 / 0.05061 / 0.01944 at f = 1 / 2 / 3.5 / 5 | 172 | **all exact** |
| 15 | 177 §4 Spearman table: a15 +0.485 expr / −0.391 \|H\| in-set | 177 | **+0.485 / −0.391** |
| 16 | per-member 1/span fitted on the 581-entry board vs AGENTS.md §4's published spans | AGENTS.md §4 | **1.129 / 3.325 / 2.721 / 1.648 — all inside** |

Scripts: `board181.py`, `mse181.py`, `reach181.py`, `reach181b.py`, `nreal181.py`, `masks181.py`,
`fid181.py`, `mech181.py`, `axis181.py`, `genelen181.py`, `lengths181.py`, `a21back181.py`,
`pool181.py` (this last one is an inconclusive control and is reported as such in §12).

---

## 1. 🔴 CONTRADICTED — "the emission law delivers the ordering" (177), and 179's board prediction

**1. What is claimed.** Replacing a15's poisoned tail block with resampled real control cells
realises **81–110 % of the `|H|` reach counterfactual** (h1-test Δ`reach` +0.035 to +0.056 raw)
with `pds`/`mse` pinned by identity, so `hy_kt8000` nets **−0.003 to +0.001** overall and is "a
wash that buys `n_real` information".

**2. What it rests on.** One surface (h1-test, 85 perturbations), matched-seed twins, eleven
known-answer gates — **all eleven of which are harness-internal reproductions of 157 / 164 / 173 /
INBOX-4 numbers.** Not one is a board point of the resample family. The `pds`/`mse` pins are real
identities (900/900 delivered pseudobulk digests, 179 §4 — verified non-vacuous). The DE members
are counterfactual-then-realised on h1-test only.

**3. How it is wrong.** **`a21` (2026-09-01) is that family's board point and it was never used.**
`submissions/a21-resample-x014.h5ad.meta.json` differs from a15 in `resample_x: 0.14` and nothing
else (alpha 0.982, keep_top 8000, keep_spread 0.25, tail 1.0, φ 4, min_cpm 10, cd4_target_f 0.1 all
identical) — i.e. *keep the top 14 % of genes by `|effect|` on the parametric law and resample the
rest.* That is `hy_kt8000`'s mechanism at a larger dose: `round(0.14 × 18,533) = 2,595` genes stay
parametric, and since `profile_from_effect` zeroes the effect outside the pb-CPM-10 gate all 2,595
are gated genes, so **a21 keeps ≈ 28 % of the ~9,100-gene gated axis parametric — it is
`hy_kt≈2600`, sitting between 177's `hy_kt2000` (18.6 % kept) and `hy_kt4000` (37 %), which is
exactly where 177 measures its LARGEST reach gain.** Its board deltas are exact (the board scorer
has zero run-to-run variance, `docs/LOG.md` 2026-09-01). Retrodicting it with 177's own instrument
(`a21back181.py`, converted with the board-fitted 1/span):

| h1-test analogue | member | Δ raw | → scaled | **a21 realised (board, scaled)** | sign |
|---|---|---:|---:|---:|:-:|
| `hy_kt2000` | **`reach`** | **+0.04424** | **+0.04996** | **−0.018482** | **WRONG** |
| `hy_kt4000` | **`reach`** | **+0.04778** | **+0.05396** | **−0.018482** | **WRONG** |
| `rs_all` | **`reach`** | **+0.03483** | **+0.03933** | **−0.018482** | **WRONG** |
| `hy_kt2000` | **`nmae`** | **+0.00054 (better)** | **+0.00088** | **−0.064214 (worse)** | **WRONG** |
| `rs_all` | **`nmae`** | **+0.00647 (better)** | **+0.01067** | **−0.064214 (worse)** | **WRONG** |
| `hy_kt2000` | `fid` | −0.04013 | −0.13342 | −0.029464 | ok, 4.5× hot |
| `hy_kt2000` | `jac` | −0.02820 | −0.07672 | −0.016583 | ok, 4.6× hot |

**h1-test gets 2 of 4 DE members wrong in sign for this emission law**, and one of them is `reach`,
by **0.07 scaled** — h1-test says +0.050…+0.054, the board says −0.018.

**The one check that could have overturned this passes.** 177 §7 names two call-site contracts
whose violation inverts the arm — (a) raw control counts with the cells' own libraries,
(b) `control_profile` = the *pooled* control profile. At commit `9e521f5`, the tree state that built
a21/a22/a23, the builder sets `resample_ctrl = control_matrix` (the raw `context_{A,B,C}.h5ad` CSR,
never rescaled), `libraries = control_matrix.sum(axis=1)`, and
`resample_profile = totals / totals.sum()`. **Both contracts were honoured. a21 is a valid board
test of the resample emission law.**
177 §7 asserts "157 §1 / 165 §6 validate h1-test for `reach` on exactly this arm class (emission
law, fixed profile)" — but 157's validation pair is **a8 vs a15, a `spread` *policy* change inside
the quasi law**, not an emission-law change. `a21`/`a22` are the only emission-law board points and
h1-test was never checked against either.

177 §6's headline table fails the same test outright: for `rs_all` at the settled `n_real ≈ 240` it
predicts **+0.0211 overall** on the four DE members; `a21` realised **−0.0215**. Off by **+0.0426
overall, sign inverted**, on the nearest board point of its own family.

And the board says the *precision* premise is wrong too: a21 cut challenge `n_pred` 4,420 → 1,141,
deleting ~86 % of the call set including the whole chance-precision tail, and board **raw `fid` fell
0.5034 → 0.4945 (−0.0089)** while h1-test precision *rose* +0.0105. That is the sixth board
measurement (after v14 / v15 / a8 / a15 / v11) saying **shrinking the call set has never once raised
board `fid`.**

**Rule 13's −0.0146 toll is not a formality — it is exactly `mean(a21, a22) reach = −0.014569`,
two exact board measurements of the very emission law 177 proposes.** 177's counter-argument
("rule 13's reach was fitted on arms inside the seed spread") addresses PLOLO's 2/6 signs, a
*different* term, and does not touch the toll. 177 cites a22 as "the wrong-way arm" and never
mentions a21, the right-way arm, which lost 73 % more `reach`.

**4. Cheap test — RUN.** `a21back181.py`, 40 s, no new data.

**5. Verdict: CONTRADICTED.** Load-bearing on: 179's package (`submissions/hy-kt8000.h5ad`,
built, not submitted), the D/E/F transition-plan item 5, and CURRENT_STATE's claim that the
+0.0525 ordering deficit "is now a solved problem".

**Corrected six-member estimate for `hy_kt8000`** (board evidence, not h1-test):

| member | 179 §8e | 181 (board-anchored) | basis |
|---|---:|---:|---|
| `pds` | 0.0000 | **0.0000** | identity, 900/900 digests — agreed |
| `mse` | 0.0000 | **0.0000** | identity + saturation above the knee — agreed |
| `nmae` | +0.0000 | **−0.0024 … 0.0000** | a21 −0.0642 scaled at 86 % dose; a22 +0.0015 at 10 % *head*. The (2·acc−1) law (LOG 2026-09-01) says a *tail* resample at acc 0.4999 is negative |
| `reach` | +0.0009 | **−0.0031 … −0.0007** | a21 −0.0185, a22 −0.0107 scaled; flat toll or dose-scaled |
| `fid` | −0.0005 (or +0.0029) | **−0.0005** | the board's five-point `n_pred` series interpolated at 3,528: −0.00098 raw = −0.00054 overall (§7) |
| `jac` | −0.0030 | **−0.0028** | h1-test −0.0142 raw × board/h1-test jac ratio 0.44 |
| **net** | **−0.003 … +0.001** | **−0.009 … −0.004** | |

The arm is most likely a **half-point loss**, not a wash. The NO-GO stands and hardens; what
changes is that the "spend a slot for regime information" argument now costs ≈ 0.005 overall, and
a21 already bought most of that information.

---

## 2. 🟠 WEAK — the 0.46 board-scaling of DE deltas (179 §8d) double-counts

**1. Claimed.** The board `n_pred` cut is −21 % against h1-test's −46 %, so every DE delta scales
by **0.46**.

**2. Rests on.** One measured pair of `probe_de.py` counts (4,474 → 3,528 board; 6,066 → 3,304
h1-test). The `reach` number it scales was *already* multiplied by the h1-test→board shrink
0.393–0.428.

**3. Wrong because** that shrink was calibrated on **a8 vs a15**, and that pair's own board/h1-test
`n_pred`-cut ratio is **0.327** (`fid181.py`): board 3,836/4,424 = −13.3 % against h1-test
3,581/6,034 = −40.7 %. The surface's smaller call-set response is therefore *already inside* the
0.428 constant. Applying 0.428 × 0.46 = 0.197 charges the same phenomenon twice. If anything
`hy_kt8000`'s ratio (0.464) is *larger* than the calibration pair's, so the transfer should be
better than 0.428, not worse.

Independently, 177's own K-series shows Δ`reach` is **not proportional to the `n_pred` cut** at
all — `hy_kt2000` (n_pred 1,579) +0.0523, `hy_kt4000` (2,307) +0.0559, `hy_kt8000` (3,304) +0.0352,
`rs_all` (1,113) +0.0413 — non-monotone. There is no law to scale by.

**4. Cheap test — RUN** (`fid181.py` §c). Also: the 0.428 shrink is a **one-point calibration with
±50 % width** — board +0.00790 over an h1-test twin-mean difference of +0.01847 whose own SE is
±0.0065, giving a shrink range **0.32–0.66**.

**5. Verdict: WEAK.** Load-bearing on 179 §8d/§8e and on CURRENT_STATE's "every DE delta scales by
0.46". It does not change §1's verdict (the sign is wrong before any scaling), but it should not be
inherited as a rule.

---

## 3. 🟠 WEAK — the 22:40 correction's "≈ −0.0006 overall" is not reproducible from its own parts

**1. Claimed.** On h1-test's real `n_real` distribution `debench` reads `hy_kt8000` at `fid`
−0.0010, `jac` −0.0064, `reach` +0.0020 overall; board-scaling the `jac` loss puts the arm at
**≈ −0.0006**, a wash.

**2/3.** The three component numbers reproduce exactly (`reach181.py`): `fid` 0.50789 → 0.50605
(−0.00183 raw = **−0.00102** overall), `jac` 0.06830 → 0.05410 (−0.01420 raw = **−0.00644**),
`nmae` 0.96658 → 0.96615, recall 0.6612 → 0.3861; and the tercile decomposition reproduces
(low +0.00014, mid +0.00048, **high −0.04372**). But two things do not hold together:

* **The h1-test→board shrink is applied to `reach` and to nothing else.** `reach` +0.0271 raw
  becomes +0.0020 overall only after ×0.393; `fid` and `jac` are carried at 1.0. Applied
  consistently the arm reads **−0.0022** (no shrink anywhere: +0.0051 − 0.0010 − 0.0064 = −0.0022;
  shrink everywhere: +0.0020 − 0.0004 − 0.0025 = −0.0009). Neither is −0.0006.
* Halving the `jac` loss for the board's lighter tail gives −0.0032, and
  +0.0020 − 0.0010 − 0.0032 = **−0.0022**, not −0.0006. I cannot recover −0.0006 under any
  combination of the stated factors.

**4. Cheap test — RUN** (`reach181.py`). **5. Verdict: WEAK** — the *conclusion* (a wash) is right
for the h1-test instrument and the *arithmetic* is not reconstructible. It is superseded by §1
anyway.

---

## 4. 🟠 WEAK — "81–110 % of the counterfactual realised" is a two-seed artefact

**1. Claimed.** `hy_kt8000` 69 %, `hy_kt4000` 110 %, `hy_kt2000` 102 % of the `|H|` ceiling.

**2. Rests on** a base twin mean of two seeds. **3.** The base pair's own spread is
`base_mf` 0.18432 vs `base_mf_b` 0.16805 = **0.01628 raw**, which is **31 % of the whole 0.0525
counterfactual** being measured. Against `base_mf` alone the same arms realise **52 % / 91 % / 84 %**;
against `base_mf_b` alone they exceed 130 %. Twin-averaging is the right move, but with n = 2 the
SE of the base is ±0.008 = ±15 % of the ceiling, so "81–110 %" should read "roughly 55–110 %,
±15 points".

**4/5. Verdict: WEAK on the headline number, SOLID on the mechanism.** The mechanism is
seed-robust and reproduces exactly (`mech181.py`): in-call-set `Spearman(−log p_adj, ·)` moves
a15 **+0.485 expression / −0.391 `|H|`** → `hy_kt8000` **+0.249 / +0.485** → `hy_kt4000`
**+0.090 / +0.514**. The p-order really does stop tracking expression and start tracking `|H|`.
What §1 shows is that **the board does not reward it.**

---

## 5. 🟠 WEAK — 179 §8c's `fid` upside (+0.0029) contradicts 164 §4 and 171 §7a(iii)

**1. Claimed.** With `n_real ≤ 3,528` the yield cap saturates, `fid` = precision, and h1-test's
precision gain (+0.00516) carries straight over: **+0.0029 overall.**

**3. Wrong because** the board has its own five-point `fid`-vs-`n_pred` series on emission-only
changes at a fixed profile (164 §4): `n_pred` 1,563 → 0.49126, 1,911 → 0.49654, 3,836 → 0.50210,
4,424 → 0.50220, 5,961 → 0.50178. Interpolated (`fid181.py`), the a23 → `hy_kt8000` move
4,474 → 3,528 reads **0.50219 → 0.50121, i.e. −0.00098 raw = −0.00054 overall** — *opposite in
sign* to 179's +0.0029 and different by 0.0034 overall. 171 §7a(iii) charged its own `keep_top`
arm −0.001141 from this same series and called it "the direct measurement, and the direct
measurement wins". 179 does not cite it. **a8 is the direct board test of removing the tail law
(n_pred 3,836, `fid` 0.50210 against a15's 4,424 / 0.50220) and it bought −0.0001**, and a21 (§1)
bought −0.0089.

**4. Cheap test — RUN** (`fid181.py` §b). **5. Verdict: WEAK→CONTRADICTED** on the +0.0029 branch.
Load-bearing on 179 §8e's "honest range −0.005 to +0.003" and on the slot-value argument.

---

## 6. 🟢 SOLID (strengthened) — 170's ρ-cliff and 180's per-context knee alpha

**1. Claimed.** `raw mse(s)` falls to the knee `s* = √(claim/B)` and rises after it; a15 sits at
ρ = 0.99 / 0.94 / 0.88 per context; per-context knee alphas 0.982 / 1.013 / 1.048 are worth
**+0.000838** overall; gamma 0.65 is a local optimum because lowering gamma raises `A`/`B` past
`claim`.

**2. Rests on** 096's identity `MSE_c − 1 = (s²A_c − 2sX + Kc_c − min(s²B_c, claim_c))/ΣD`, with
`claim_c` capped by a `C_p` proxied from 400-cell NTC jackknives and `X`, `ΣD` board-anchored and
**not separable per context** (096 §(e) says so explicitly, 180 §6 inherits it).

**3. Two ways it could be wrong, both tested.**

*(a) Is `min(s²B, claim)` the right form, and is the knee the minimum?* Yes, and with wide margins.
Below the knee `d(num)/ds ∝ 2s(A−B)_raw − 2X_raw = 0.00266 − 0.03906 < 0`; above it
`∝ 2sA_raw − 2X_raw ≫ 0` since `X_raw/A_raw = 0.070` against `s ≈ 0.91`. The knee is the minimum
unless `X_c` is wrong by **14×** in either direction.

*(b) Is the 400-cell NTC `C_p` proxy valid?* The `manifest.json` gives `ground_truth_cells` 138,400
over 300 perturbations = **461 cells per reference group, not 400**, and `C_p ∝ 1/n`, so the proxy
could be ~13 % high — which would move the knee from α 0.982 to **0.915** and invert 180's whole
recommendation. **Tested against the board and refuted** (`mse181.py`): anchor `X` on a14's own
board point (α 0.982, raw `mse` 0.882761) and predict a12 out of sample (α 1.0108):

| `C_p` × | predicted a12 raw `mse` | board 0.901243 | residual | global argmin α |
|---:|---:|---:|---:|---:|
| 0.900 | 0.928794 | | +0.02755 | 0.932 |
| **0.868 ("461 cells")** | ≈ 0.9345 | | ≈ +0.033 | 0.915 |
| 0.980 | 0.904302 | | +0.00306 | 0.973 |
| **0.990** | **0.900941** | | **−0.00030** | **0.978** |
| **1.000 (shipped)** | **0.898778** | | **−0.00247** | **0.982** |
| 1.050 | 0.881776 | | −0.01947 | 1.007 |

**The board's own pure-α pair validates the `C_p` proxy to ≈ ±1 %** and puts the global argmin at
α 0.978–0.982, where the shipped value sits. A 13 % biased proxy is off by 0.033 raw — **22× a
genuine board move** (a23 moved 0.00146). The board also brackets the minimum directly:
α 0.838 → 0.9099, **0.982 → 0.8828**, 1.0108 → 0.9012.

*(c) Does the unmeasurable per-context `X_c`/`ΣD_c` break the +0.000838?* No — it rescales it.
`knee_gain` runs **+0.00039 … +0.00175** for `X_c` × 0.5 … × 2.0 and **+0.00072 … +0.00102** for
`ΣD_c` ∈ [0.19, 0.27] (096's own sensitivity range). The **sign is robust; the magnitude is
±2×.** The asymmetry argument (10 % long costs 0.090 scaled `mse`, 30 % short costs 0.030) is
structural — above the knee the whole `2sA` is charged, below it only `2s(A−B)` — and does not
depend on `X_c` at all.

**4. Cheap tests — RUN.** `mse181.py` reproduces 180's knee alphas to 5 dp, a23's raw `mse` to
2e−6, and the +0.000838 to +0.000843.

**5. Verdict: SOLID**, and the strongest result of the day. One bookkeeping note: 180's **shipped**
(all-46-NTC) alpha for context A is **0.982338**, not the 0.98586 that CURRENT_STATE quotes as
"computed from the bundle alone" — 0.98586 is the *096-`C_p`* variant used as a known-answer row.
The `mse` gain is unaffected (A contributes +0.000053 of the +0.005025 scaled).

---

## 7. 🟢 SOLID (strengthened) — the long-gene axis identification

**1. Claimed.** The Orion pair's shared PC1 (cos 0.872, 19.7 % of HEK293T's variance) is "the
nuclear-RNA-fraction / long-gene technical axis", identified from the top-30 loading gene names.
Now in `AGENTS.md`'s harness warning.

**2. Rests on** a gene-name eyeball of 30 genes out of 10,779. The obvious alternative reading:
long intron-rich genes are also *low-expression* genes and the positive-loading genes (PSMB3, LSM7,
NME1, UBA52, RACK1, TUBA1B…) are short *housekeeping* genes, so an **expression** axis would
produce exactly the same top-30 list.

**3/4. Tested** (`axis181.py` + `lengths181.py` + `genelen181.py`): PC1 recomputed (variance shares
5.64 / 8.60 / 19.65 %, cosines 0.872 / 0.307 / 0.234 — all KA rows pass), Ensembl genomic spans
fetched for **10,601 of the 10,779** axis genes:

| donor | ρ(loading, log span) | ρ(loading, log CPM) | **partial ρ(loading, span \| CPM)** | partial ρ(loading, CPM \| span) |
|---|---:|---:|---:|---:|
| **orion_hek293t** | **−0.714** | +0.127 | **−0.720** | +0.181 |
| **orion_hct116** | **−0.654** | +0.135 | **−0.659** | +0.178 |
| k562_gwps (control) | **−0.022** | −0.166 | **−0.022** | −0.167 |

`ρ(log span, log CPM)` on this gated axis is **−0.0008**, so expression is not a confound at all;
the mean loading is monotone across all ten span deciles (+0.0051 → −0.0170); median span is
**712 kb** for the top-30 negative loadings against **40 kb** for the axis and **8 kb** for the
top-30 positive. K562 — which 175 says lacks the axis — reads **−0.022**, exactly the null this
test needed to be non-vacuous.

**5. Verdict: SOLID.** The identification is not merely supported; it is now a −0.71 correlation
over 10,601 genes with a clean donor-level control. Keep the AGENTS.md warning.

---

## 8. 🟢 SOLID — 171's structural closures, and 164's decomposition

* **"Our significant calls first" is a BH no-op.** `reach` with the call set = production / nothing
  / everything all read **0.18283** (`mech181.py`). The call set cannot move `reach` in either
  direction. *(My monotonicity counter reads 83/85 rather than 171's 0/85 only because
  `debench` zeroes the self gene out of `pred_sig`, pushing the smallest `p_adj` into the
  "not called" group; the reach identity is the decisive form and it passes.)*
* **164's call-set decomposition** reproduces to the digit: 297,869 keep / 214,265 tail / 512,134
  total, precision **0.5122 / 0.4999 / 0.5070**, tail share **41.8 %**, paired keep−tail **+0.0122**.
* **164's `mse` closure** is confirmed by §6's board bracket: every tail arm drops `claim` below the
  knee `s²B` and lands in the linear branch. SOLID.

---

## 9. 🟢 SOLID — the board decomposition, with one number wrong

Recomputed from `data/leaderboard/2026-09-02T0630Z.json` (`board181.py`).

| member | us | leader | edge/6 | share of gap | **our rank / 581** | field best | headroom |
|---|---:|---:|---:|---:|---:|---:|---:|
| mse | 0.1216 | 0.2375 | +0.01931 | 44.4 % | 23 | 0.2375 | +0.0193 |
| reach | 0.1074 | 0.2049 | +0.01625 | 37.3 % | **192** | 0.2379 | +0.0218 |
| pds | 0.7476 | 0.8165 | +0.01148 | 26.4 % | **5** | 0.8165 | +0.0115 |
| nmae | 0.0941 | 0.1207 | +0.00444 | 10.2 % | 68 | 0.1838 | +0.0149 |
| jac | −0.0005 | −0.0096 | −0.00151 | −3.5 % | 169 | 0.0196 | +0.0034 |
| fid | −0.0301 | −0.0688 | −0.00646 | −14.8 % | 176 | 0.0227 | +0.0088 |

Every figure reproduces; edges sum to the gap 0.04351 exactly; **the 5/23/68 vs 176/169/192
partition is real.** `nmae` at rank 68 sits between the two groups but on the *profile* side of the
boundary, and the boundary is a mechanism (does the metric pass through our emitted cells'
Wilcoxon calls?), not a rank cut — `nmae` reads our `lfc` over the *reference's* gate, it never
reads our call set. **The partition holds.**

**One number is wrong.** CURRENT_STATE and `docs/LOG.md` state `raw_reach` and `raw_pds` are
"coupled at **corr +0.831**". The regression they quote reproduces exactly
(`raw_reach = −0.0751 + 0.3577·raw_pds`, residual sd 0.0410, us at **−1.22 σ**, 39 of 581 further
below the line = 6.7 %), but its correlation over the same 581 entries is **+0.7013**
(Spearman +0.647; on scaled units +0.7003). No filter of the snapshot recovers 0.831. The
conclusion is unaffected — the residual, not the correlation, carries it — but the number should be
corrected.

---

## 10. 🟠 WEAK — "`pds` vs per-gene accuracy is the unifying trade, and production sits at the `pds`-optimal point"

**1. Claimed.** Every lever found on 2026-09-02 trades `pds` against every other member, "because
`pds` rewards decorrelation and every other member rewards per-gene accuracy, and production sits
at the `pds`-optimal point of that curve."

**2. Rests on** four riders — gamma (097/170), sign repair (173), K562-heavy donor weights (175),
per-donor deflation (175) — plus whitening itself (172). **All five are un-whitening moves.** The
synthesis is a description of one family generalised into a law, and 178 §10.5 then writes it into
CURRENT_STATE as "production sits at the `pds`-optimal point *and* at the `mse` optimum of the
amplitude curve, and the two pins are independent."

**3. Three problems.**
* **169 is a same-day counterexample in its own record.** Magnitude sharpening raises h1-test `pds`
  by **+0.019452 scaled** (the largest `pds` move in the project) and dies on **`nmae` −0.01289**.
  That is `pds` *up*, per-gene member *down* — the opposite trade. So "every lever trades `pds`
  against every other member" is false as stated; what is true is "every *un-whitening* lever does".
* **170 §2 point 5 directly contradicts "production sits at the `pds`-optimal point":** on the
  *delivered*, mass-fixed object `pds` has a **seed-free interior maximum at gamma 0.50–0.55**
  (+0.00126 raw / +0.00046 scaled), not at the shipped 0.65. 178's PLOLO reads gamma 0.45 at
  −0.0093 raw `pds`, the opposite sign. The two instruments disagree about the very claim.
* The two "independent pins" are not independent in the way stated: 170 §6 says the `mse` optimum
  in gamma "is not a property of the operator alone — change `VCC_ALPHA`, the dispersion policy or
  `keep_top` and this optimum moves with them."

**5. Verdict: WEAK.** True of the un-whitening family, unproven as a law, and contradicted on the
`pds`-optimality half by 170's own delivered-space measurement. Load-bearing on "only new
information can improve the predictor", which is a much stronger claim than the evidence supports —
it was established by testing five moves of one shape.

---

## 11. 🟠 WEAK — "h1-test understates `pds` costs of un-whitening arms ~6×; PLOLO is the `pds` instrument"

**1. Claimed** (178 §8, now in CURRENT_STATE's instrument table). Sign repair costs **−0.0245 raw
`pds`** on PLOLO where 173 read **−0.0042** on h1-test.

**2/3. Three of the 6× are not instrument error.**
* 178's **own** h1-test rebuild reads sign repair at **−0.00840**, 2× 173's −0.0042, on the same
  nominal surface. Its per-row cosine to the recovered artifact is **0.672** (178 §2), so a large
  part of the gap is reconstruction, not instrument.
* `pds` is a **retrieval** statistic and h1-test ranks **85** decoys where PLOLO ranks **269–294**
  and the board ranks **299**. 162 §6 measured exactly this frame factor on an identical
  quantity: 156's `intended − delivered` reads **+0.00324/+0.00382** on the 132-decoy Orion frame
  and **+0.0173/+0.0174** on the 290-decoy challenge frame — **4.5–5.4×, for the same geometric
  change.** A 6× ratio between an 85-decoy and a 290-decoy frame is inside the documented factor.
* **The rule cannot be applied generally.** 169's h1-test `pds` gain is +0.019452 scaled =
  +0.0032 overall; ×6 would be **+0.0194 overall**, larger than the entire `pds` headroom to the
  field best (**+0.0115**). Impossible.

**4. Cheap test — ATTEMPTED, INCONCLUSIVE.** `pool181.py` subsamples h1-test's 85 perturbations to
10/20/40/60 and measures Δ`pds` for a 173-shaped sign flip. The deltas (−0.0000 … −0.0013) are
dominated by subset noise and show no clean growth. Reported as a failed control rather than
dropped — `debench.pds` reads the DE table's mean-of-ratios `lfc`, not the pseudobulk, which 171
§6 already warns is a proxy, so this was the wrong instrument for the question.

**5. Verdict: WEAK.** PLOLO *is* the right frame for board `pds` (it matches the board's decoy
count and retrodicts 6/6 signs) — but the lesson should be recorded as **"`pds` deltas scale with
the decoy pool; always read `pds` in the challenge frame"**, not as "h1-test understates by 6×".
The second form, applied to a non-un-whitening arm, produces impossible numbers.

---

## 12. 🟢 SOLID — 176's CMap NO-GO (and the coverage arithmetic is *conservative*)

**1. Claimed.** CRISPR-KO adds +0.019 sign accuracy on 867 landmark genes but coverage caps the
challenge effect at ≈ +0.0007 pooled: `82/300 targets × ~13 % of pool pairs × 0.019`.

**3. The specific worry** the audit was asked to test: if landmark genes are *over-represented* in
the reference's DE pool, the pooled effect would be larger than 13 % implies.

**4. Measured** (`masks181.py`, using 176's own `work/landmark_mask.npy`): landmarks are **8.04 %**
of the 10,779-gene DE universe and **9.53 %** of the 65,044-pair reference DE pool — enrichment
**1.19×**, and *below* the 13 % 176 assumed. The corrected arithmetic is
`82/300 × 0.0953 × 0.019 = +0.000495`, i.e. **176's cap was already generous by 40 %.**
The `82/300 × 0.13 × 0.019 = +0.000675` arithmetic is itself correct.

**5. Verdict: SOLID.** NO-GO confirmed and understated. No action.

---

## 13. 🟠 WEAK — 175's "the Orion pair adds ~nothing to sign accuracy"

**1. Claimed.** K562 alone 0.5670, the Orion pair alone 0.5425, all three 0.5697 — "the Orion pair
carries two-thirds of the weight and adds +0.0027 over K562 alone."

**2. Rests on** a common mask of **39,945** pool pairs where all three donors carry a non-zero
entry. That mask reproduces exactly here, and the accuracies reproduce to 4 dp.

**3. The mask is a subset of K562's 8,248-gene axis, and it hides a third of the pool.**
`masks181.py`:

| stratum | n pairs | K562 | Orion pair | equal(3) |
|---|---:|---:|---:|---:|
| all three carry (175's mask) | 39,945 | 0.5670 | 0.5425 | 0.5697 |
| **K562 carries NOTHING, both Orion do** | **22,707** | — | **0.5294** | **0.5294** |
| whole signable pool (≥ 1 donor) | 64,645 | — | — | **0.5555** |

**K562 can sign only 64.1 % of the signable pool.** On the other 35 % the Orion pair is the entire
consensus and reads 0.5294 — above chance and the only source there. Against a "K562 or nothing"
baseline (0.641 × 0.5670 + 0.359 × 0.500 = 0.5430) the Orion pair is worth **≈ +0.012**, not
+0.0027.

**5. Verdict: WEAK as stated.** It does not change any decision (178 killed K562-heavy weights on
`mse` and on the operator anyway, and its own measurement — whitening consumes the sign gain,
+0.0043 → **−0.0012** after the operator — is the binding one). But **CURRENT_STATE's sentence
"the Orion pair adds ~nothing" should not be carried into D/E/F donor-mix decisions**: on a third
of the pool the Orion lines are the only donors we have.

---

## 14. 🟢 SOLID / 🟠 — the `n_real` inference, and a new bound nobody had

**1. Claimed.** Board `jac` 0.0301 at `n_pred` 4,424 implies recall **0.57–0.71** if `n_real` is
191–240 (matching h1-test's 0.67) and 0.078 if it is 2,670; the weight of evidence is ~240.

**2/3. The algebra is right but it is a mean-field inversion of a mean-of-ratios metric, and it
fails on the only surface where it can be checked.** `jac` is the *mean of per-perturbation
Jaccards*. Inverting the panel mean (`nreal181.py`): on h1-test, true mean recall is **0.6685**;
the inversion returns **4.395** with the median `n_real` (89) — a 6.6× error — and **0.568** with
the mean `n_real` (765), a 0.85× error. **So the inversion is only usable with the *mean*
`n_real`, and CURRENT_STATE quotes 191–240 as a typical value.** If the challenge's distribution is
heavy-tailed like h1-test's, jac 0.0301 with mean `n_real` ~765 implies recall 0.19, not 0.67.

**4. A cheap test nobody had run, and it is a new and much tighter bound.** `reach = mean_p k*/n_real_p`
is monotone-decreasing in power, and the board reads our raw `reach` at **0.17420 (a15) / 0.17416
(a23)**. Proportional power matching on h1-test (`nreal181.py`; f = 1 reproduces the native BH
surface exactly, and f = 2 / 3.5 / 5 reproduce 172 §7.3 to 5 dp):

| f | median `n_real` | mean `n_real` | reach |
|---:|---:|---:|---:|
| 0.50 | 44 | 383 | 0.22927 |
| **1.00** | **89** | **765** | **0.18283** ← KA |
| 1.25 | 111 | 956 | 0.16329 |
| 2.00 | 178 | 1,528 | 0.10953 |
| 3.50 | 312 | 2,415 | 0.05061 |
| 5.00 | 445 | 2,950 | 0.01944 |
| **30.0** | **2,670** | 5,148 | **0.00340** |

**The board's own `reach` places the challenge at f ≈ 1.05: median `n_real` ≈ 90–110.** 093's
2,670 would read reach **0.0034 — 51× below the board's 0.1742** — which refutes it through a
member that never touches `fid`, precision or the `C_p` proxy, so it is not circular.

**A consequence that is not currently written down anywhere:** the board's `reach` (dominated by
the *small*-`n_real` stratum) says the challenge *has* one, while board `jac` 0.0301 against
h1-test's 0.0683 (dominated by the *large*-`n_real` stratum) says its heavy tail is about half as
prevalent. **`n_real` is a distribution, not the scalar 240.** Treating it as a scalar is the same
error class as rule 8, and it is exactly why 177 §6's single-`n_real` table had to be reversed at
22:40 — and why a21's "coverage did not engage at n_pred 1,141" inference (LOG 2026-09-01) is not
clean: under a heavy tail the yield cap engages for the top ~20 % of perturbations even at 1,141.

**5. Verdict: the ~240 reading is PROBABLY RIGHT as a median (I now bound it to 90–110 and confirm
191 ≤ n_real), 2,670 is CONTRADICTED independently, and the scalar framing is WEAK.**

---

## 15. 🟢 SOLID — the ordering counterfactual is real but concentrated where the challenge is thin

**Claimed.** The `|H|` ordering gap is **+0.0525 raw ≈ +0.0097 overall**, "our whole deficit below
the field line".

**Measured** (`reach181b.py`), the same counterfactual split by `n_real`:

| stratum | n perts | mean Δ`reach` | share of the panel Δ |
|---|---:|---:|---:|
| `n_real ≤ 13` | 31 | +0.09793 | **63.7 %** |
| 14 – 800 | 30 | +0.04912 | 34.2 % |
| > 800 | 24 | +0.00374 | 2.1 % |
| `n_real ≥ 20` | 50 | +0.01419 | 16.5 % |
| `n_real ≥ 100` | 41 | +0.00804 | 7.7 % |

and the same for the realised arms: `hy_kt8000` +0.03523 full panel → **+0.00937 on `n_real ≥ 20`**
→ +0.00233 on `n_real ≥ 100`; 67 % of its gain sits in the 31 perturbations with `n_real ≤ 13`.
This is 171 §7a(i)'s kill shot for `--keep-top 8300` ("3× smaller on the stratum the challenge
resembles") applied to 177's arms, which 177 never ran.

**But §14 rescues it partly:** the board's own reach level says the challenge *does* have a
small-`n_real` stratum (f ≈ 1.05), and 172 §7.3 measures the `|H|` gain holding at f = 1 (+0.0525),
f = 1.5 (+0.0513) and f = 2 (+0.0612) before collapsing at f = 3.5 (+0.0148). **So the
counterfactual survives at challenge power.** What does not survive is its *delivery* — §1.

**Verdict: the counterfactual is SOLID at challenge power; the claim that it is "now a solved
problem" is CONTRADICTED (§1).**

---

## 16. 🟡 PROBABLY RIGHT — 172's "the ordering gap is not learnable", including the leakage question

**The leakage worry** the audit was asked about: in the `orion_hct116` fold the held-out line's
platform axis (cos 0.872, §7) is shared with a *training* source (`orion_hek293t`), so
"cross-source agreement" could be predicting platform agreement, not biology; and the `k562_gwps`
fold trains on two same-atlas lines.

**Three things make the leak non-load-bearing.** (i) 172 already flags the `k562_gwps` fold as
"weakest by construction" and reports it as consistently the worst. (ii) Its §6 ablation shows the
exposed feature group — **cross-source agreement** — can be dropped with **+0.0022 fold reach /
−0.0006 h1 reach**, i.e. it carries nothing; the fold signal lives in `expression` (−0.0069) and
`donor-panel responsiveness` (−0.0067), neither of which reads the shared axis. (iii) The leak's
direction is *conservative*: it would inflate apparent learnability, and 172's verdict is that
nothing beats `|H|`.

The one real cost is **selection**: the LODO-selected class (M1, 0.22448) is *below* `|H|` while the
unselected M4 gbt (0.24287) is above it, and if the folds are mildly corrupted the selection is
too. 172 reports both and the paired bootstrap cannot separate either from `|H|`
(P = 0.074 / 0.778), so the conclusion holds. **Verdict: PROBABLY RIGHT.**

172's other finding — that its own briefed permutation null (N1, within-target shuffle) **leaks**,
reading 0.209 against the honest N2's 0.177 — is exactly the class of vacuous control this project
has been burned by twice, caught and reported by the experiment itself. Good.

---

## 17. Internal contradictions between documents

| # | contradiction | resolution |
|---|---|---|
| 1 | **097 §9h "gamma 0.30 = +0.0027" vs 170 "every gamma is net negative"** | Resolved by 170 and confirmed here: 097 priced `mse` through `MSE_min = 1 + (Kc − 2√(claim·W))/ΣD`, which **assumes the arm sits ON the cliff**. §6's board bracket confirms 170. **No live contradiction.** |
| 2 | **157/154 "the ordering channel is closed" vs 177 "the emission law delivers it"** | Both partly right: 177's mechanism is real (§4's Spearman flip) and 157's *board* conclusion is right (§1's a21 point). The reconciling statement is that h1-test's `reach` does not transfer for emission-**law** changes, only for spread-**policy** changes. |
| 3 | **164 "removing the tail is worth +0.0009 and costs −0.0185 on `mse`" vs 177 "the same DE gain with the `mse` cost removed"** | 177's coordinator addendum 2 states this correctly and it is right on `mse` (real cells keep the dispersion credit; `Ĉ/C_p` 1.094 > 1.085, above the cap). But it inherits 164's h1-test DE *gains*, which §1 shows the board does not pay. |
| 4 | **173's h1-test `pds` −0.0042 vs 178's PLOLO −0.0245 (6×)** | §11: ~2× is reconstruction (178's own h1-test rebuild reads −0.0084), the rest is the 85-vs-290 decoy frame, which 162 independently measured at 4.5–5.4×. Not a 6× instrument bias. |
| 5 | **169 "the blanket reach toll rule 13 charges is FALSIFIED (+0.00061, whole CI above −0.00243)" vs 177/179/180 all being vetoed by that toll** | 169's arm is a *profile* change, not an emission-law change; the toll's two calibration points (a21, a22) are emission-law arms. §1 says the toll is *correct* for the family it was fitted on. Rule 13 should distinguish the two classes explicitly — right now one constant does both jobs. |
| 6 | **178/CURRENT_STATE quote a15 at ρ = 0.992/0.940/0.881; 180 measures 0.9993/0.9399/0.8783** | The `C_p` estimator: 096/170 use the first 16 NTC constructs, 180's shipped code uses all 46. Context A's knee alpha is 0.9859 under the first and **0.9823** under the second. CURRENT_STATE's "0.98586/1.01275/1.04568 computed from the bundle alone" mixes the two; the bundle-alone value for A is 0.982338. |
| 7 | **179 §8c `fid` +0.0029 vs 164 §4 / 171 §7a(iii) board series −0.0005** | §5. The board series wins; 171 said so explicitly about its own arm and 179 does not cite it. |
| 8 | **CURRENT_STATE "corr +0.831" vs the snapshot's +0.7013** | §9. Arithmetic error; conclusion unaffected. |

---

## 18. Ranked table

| # | claim | verdict | what depends on it | action |
|---|---|---|---|---|
| 1 | **177: the emission law delivers the `\|H\|` ordering; 179's `hy_kt8000` is a wash** | **CONTRADICTED** | the built 19.7 GB package; the D/E/F resample plan; "the ordering defect is solved" | **REOPEN — retrodict a21 before any resample arm is ever submitted.** Do not spend a slot on `hy_kt8000`. |
| 2 | 179 §8c: board `fid` +0.0029 if `n_real ≤ 3,528` | **CONTRADICTED** (board series says −0.0005) | 179's positive branch; the slot-value argument | **RE-MEASURE** — use 164's five-point series, as 171 did |
| 3 | 093's `n_real ≈ 2,670` | **CONTRADICTED** (new, non-circular: board `reach` would read 0.0034) | 177 §6's catastrophic column | **KEEP CLOSED** — now closed by a second, independent member |
| 4 | CURRENT_STATE: `raw_reach`/`raw_pds` "corr +0.831" | **CONTRADICTED** (+0.7013) | nothing — the −1.22 σ residual carries the argument | **FIX THE NUMBER** |
| 5 | 179 §8d: scale every DE delta by 0.46 | **WEAK** (double-counts the 0.428 shrink) | 179's whole six-member table | **RE-MEASURE** |
| 6 | 22:40: `hy_kt8000` ≈ −0.0006 overall | **WEAK** (not reproducible from its parts; mixed shrink conventions) | the "wash" framing | **RE-STATE** with one convention |
| 7 | 177: "81–110 % of the counterfactual realised" | **WEAK** (base twin spread is 31 % of the effect) | the headline | **RE-STATE as 55–110 % ± 15** |
| 8 | "production sits at the `pds`-optimal point; only new information can help" | **WEAK** (one family; 170's delivered `pds` peaks at gamma 0.50–0.55; 169 is a same-day counterexample) | the entire October research direction | **REOPEN** — test one lever that is *not* an un-whitening move |
| 9 | "h1-test understates `pds` costs ~6×" | **WEAK** (frame effect, not instrument bias) | the instrument table | **RE-STATE**: read `pds` in the challenge frame; do not apply a 6× multiplier |
| 10 | 175: "the Orion pair adds ~nothing to sign accuracy" | **WEAK** (measured on the 64 % of the pool K562 can reach) | D/E/F donor-mix decisions | **RE-STATE** |
| 11 | `n_real ≈ 240` as a scalar | **WEAK** (it is a distribution; board reach ⇒ median 90–110, board jac ⇒ mean several hundred) | `fid`/`jac` pricing everywhere | **RE-STATE as a distribution** |
| 12 | 172: the ordering gap is not learnable | **PROBABLY RIGHT** (the leak exists but the exposed features carry nothing) | closing the learned-key family | **KEEP** |
| 13 | 170 ρ-cliff + 180 per-context knee alpha | **SOLID** (board-validated to ±1 % on `C_p`; knee reproduced to 5 dp) | the whole D/E/F insurance policy | **KEEP** |
| 14 | the long-gene technical axis | **SOLID** (ρ = −0.71 / −0.65 over 10,601 genes; K562 control −0.02) | the AGENTS.md harness warning | **KEEP** |
| 15 | board partition: 3 call-set members 176/169/192 vs 3 profile members 5/23/68 | **SOLID** (exact) | the day's reframing | **KEEP** |
| 16 | 171: the call set cannot move `reach` (BH no-op) | **SOLID** (0.18283 with the call set = nothing, everything, production) | closing the call-policy family | **KEEP** |
| 17 | 164: tail block 41.8 % of calls at precision 0.4999 | **SOLID** (exact) | the `fid` closure | **KEEP** |
| 18 | 176: CMap NO-GO | **SOLID** (landmarks are 9.5 % of the pool, not 13 % — the cap was generous) | closing public cross-line data | **KEEP** |
| 19 | the `\|H\|` counterfactual is worth ≈ +0.0525 raw at challenge power | **SOLID** (172 §7.3's f-series + §14's f ≈ 1.05 bound) | the remaining reach headroom | **KEEP — it is still unclaimed** |

---

## 19. The 3–5 highest-value things this audit says should happen next

1. **Make a board retrodiction a mandatory gate for any emission arm.** 177 and 179 carried
   eleven and seven gates respectively and every one was harness-internal. The one board point of
   their own family (a21) was on disk the whole time and inverts two of four members. Add to
   AGENTS.md §11.1 rule 10: *a screen whose family already has a board point must retrodict it, and
   the known-answer row must be the board's, not the harness's.*
2. **Do not submit `hy-kt8000`.** Corrected expectation **−0.009 to −0.004** overall (§1), and the
   `n_real` information it would buy is largely already bought by a21 and is now bounded much more
   tightly by the board's own `reach` (§14: median `n_real` 90–110).
3. **Split rule 13's `reach` toll into two constants.** The −0.0146 is an *emission-law* toll
   calibrated on a21+a22 (both exact board points, both negative) and it is correct there; 169
   falsified it for *profile*-class arms. One constant currently vetoes both classes, which is why
   177/179/180 all read as "vetoed by a formality" when for two of them it is the direct evidence.
4. **Re-open the `pds`-optimality claim with one non-un-whitening lever.** Every input to "only new
   information can help" is a move that undoes whitening. 170's own delivered-space `pds` maximum
   at gamma 0.50–0.55 and 169's `pds`-up/`nmae`-down result are both inside the same day's record
   and neither fits the stated law. A norm-preserving transform that raises `X` without raising `A`
   (170's design criterion) is still unexplored outside sign repair.
5. **Write `n_real` down as a distribution.** Board `reach` 0.1742 ⇒ f ≈ 1.05 on h1-test's power
   scale (median ≈ 90–110); board `jac` 0.0301 vs h1-test 0.0683 ⇒ the heavy tail is ~half as
   prevalent. Every `fid`/`jac` price in the repo currently indexes a scalar, and that is what
   forced the 22:40 reversal. Carry it into the D/E/F plan's item 5.

---

## 20. What I did NOT measure, and why

| not measured | why |
|---|---|
| anything through `cell-eval2` | forbidden (25 GB RSS, shared machine) |
| a rebuilt `pds` for sign repair / gamma in the challenge frame | needs PLOLO's fold machinery and a rebuild; §11 argues the frame factor from 162's own recorded numbers instead |
| ~~a21's build against 177 §7's two call-site contracts~~ | **MEASURED after the first draft** — verified at commit `9e521f5` (the tree that built a21/a22/a23): raw control CSR, own libraries, pooled control profile. Both contracts honoured; §1 stands |
| whether a21's `nmae` collapse is the resample law or the library unwelding specifically | the two are the same change on this path; only a build at a third `resample_x` could separate them, and builds are forbidden here |
| `hy_kt8000` at the challenge's own `n_real` | rule 9; §14 brackets it instead |
| 180's synthetic D/E/F rehearsal and byte-identity claims | a build; taken on 180's own gate evidence, which is internally consistent |
| 159's `1/k` law and effective-donor count | already flagged as unverified twice (159's own 2–3× overstatement, 166's ±0.001–0.002 correction); no new evidence added here |
| a decisive decoy-pool scaling law for `pds` | attempted (`pool181.py`) and inconclusive; reported as a failed control in §11 |

---

## GO / NO-GO

**NO-GO on `submissions/hy-kt8000.h5ad`, on stronger evidence than 179 had.** a15/a23 stay as they
are; zero flags change; the submission hold is untouched and no slot was spent by this audit.

**GO on `VCC_ALPHA=knee` as the D/E/F default** — its algebra is now board-validated (§6) and its
one soft spot (the `C_p` proxy) is pinned to ±1 % by the board's own α ladder.

**The day's science survives in its closures and fails in its one positive result.** Every NO-GO
audited here holds. The single claim the day recorded as a *win* — "the ordering defect is solved by
the emission law" — is the one the board has already refuted.
