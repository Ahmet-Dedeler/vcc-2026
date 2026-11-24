# 487 — Adversarial robustness review of the evidence behind the October recipe

Status: **in progress** (started 2026-09-26 21:55Z; resumed 21:57Z). Analysis only: no builds, no `def_day.sh`, no `vcc`,
no `cell-eval2`, no git, no edits outside this directory. Memory rule until 04:30Z 09-27: every process ≤ 4 GB peak RSS
(free+inactive read 27.6 GiB at 21:58Z, below the 33 GiB gate, so nothing bigger is started).

Interim findings are appended below as they land; the final answer will replace this header.

## Decision 1 — same-line rows at w 2 + per-target prior (interim, 22:40Z)

Script `sameline487.py` / `sameline487b.py` (saved per-target outputs of 450's H1 runs and 435's base; 112 MB RSS).

**Known answer (PASS):** from the per-target CSVs and mse legs, paired by seed with 435's SLOPE, every row of 450 §9/§10
reproduces: full w2 +0.0934 ± 0.0002, w1 +0.0837, w4 +0.0933 ± 0.0005, rows-only w2 +0.0790, flip 0.2 w0.5/w1/w2
+0.0323/+0.0316/+0.0268, flip 0.3 +0.0166/+0.0086/+0.0007. Uncovered targets contribute exactly 0 on every member.

| check | number |
|---|---|
| target bootstrap, all 85 resampled (full w2) | +0.0934 [+0.063, +0.125], P(≤0) 0 |
| stratified (30 covered / 55 uncovered fixed) | [+0.070, +0.115] |
| non-essential targets only (67, 17 covered) | +0.0652 [+0.034, +0.096] |
| split of +0.0934: covered essential (13) / covered non-essential (17) | +0.052 / +0.041 |
| w2 − w1 (3 seeds, target bootstrap) | +0.0097 [+0.005, +0.015] |
| w4 − w2 | −0.0001 [−0.003, +0.002] |
| flip 0.2 bank (≈0.66 signs): w2 − w0.5 | −0.0054 [−0.016, +0.003] |
| flip 0.3 bank (≈0.61 signs): w2 − w0.5 | −0.0160 [−0.029, −0.005] |

**Robust on H1 target resampling. Three things the record does not state:**
1. **H1's covered targets are strong responders.** Median n_real of the 30 Pan-covered targets: 1,862 (essential) and
   537 (non-essential), against 13 for the uncovered non-essential targets. Only 4 of 30 have n_real < 50 (SLC5A6,
   STUB1, DPF2, TOR1A: per-target gains ×85 = +0.07, −0.15, +0.17, +0.22). Per-target gain rises with strength
   (Spearman with n_real +0.46). A genome-wide October bank covers the panel's weak third too; H1 has almost no evidence
   there, so the coverage curve's +0.20 at 100 % coverage (sampling covered targets) is an upper figure. For scale:
   n_real < 50 holds for 54 % of h1-test's 67 non-essential targets (39 % below 10) and for 13 % of the 30 covered ones.
2. **w 2 is flat in w but on an edge in bank quality.** At Pan quality w 2 ≈ w 4 and beats w 1 by +0.010 (CI excludes 0).
   At ≈0.61 sign quality w 2 loses 0.016 to w 0.5 (CI excludes 0). The `auto` rule classes banks by own-target
   knockdown (Orion HCT116 0.911, just over the 0.9 line), which is not the bank's sign accuracy against Arc's truth.
3. **The +0.0934 is one context.** If one of D/E/F matches, the overall moves by a third of it.

4. **The October lines with high coverage have banks that are ALREADY production donors, and there the rows add ~no
   direction.** 451's coverage table: K562 300/300 and HCT116 300/300 (both "strong", so `auto` gives w 2); every
   new-information bank covers ≤ 83/300 (DLD-1, sparse table, weak class) or 2–9 (HepG2, A549, MCF7). K562-GWPS and Orion
   HCT116 are two of our three donors. `lift487.py` measures the direction lift of the blend on one statistic (centred
   cosine to an independent same-line truth; unwhitened 3-donor mix as "ours"; ratio of means, target bootstrap):

   | case | truth (lab) | targets | lift at w 0.5 | lift at w 2 |
   |---|---|---:|---:|---:|
   | H1, Pan = new information (the +0.0934 case) | Arc H1, all cells | 125 (62 non-ess.) | ×1.87 | **×2.50** [2.26, 2.80]; non-ess. ×3.03 |
   | same, Pan flipped p 0.2 (emitter price w2 +0.0268) | | | ×1.49 | ×1.69 |
   | same, Pan flipped p 0.3 (emitter price w2 +0.0007, w0.5 +0.0166) | | | ×1.27 | ×1.25 |
   | K562, GWPS if it were new information | VIPerturb (Satija) | 3,000 | ×1.60 | ×1.92 |
   | **K562, GWPS already in the mix** | VIPerturb | 2,517 non-ess. | ×1.10 | **×1.09** [1.03, 1.16] |
   | HEK293, Orion HEK293T if new information | XuCao (408 cells/pert) | 171 | ×1.22 | ×1.23 (non-ess. ×1.51) |
   | **HEK293, Orion HEK293T already in the mix** | XuCao | 171 (65 non-ess.) | ×1.01 | **×0.95** [0.90, 0.99]; non-ess. ×1.05 |

   - The ratio is invariant to independent truth noise, so VIPerturb's 48-cell depth does not bias it.
   - With the PRODUCTION composition (465's efr dump, 277 VIPerturb-covered non-essential targets, `lift487b.py`) the
     truth has almost no signal there (cosines 0.004–0.006): w 2 ×0.85 [0.54, 1.15], w 0.5 ×0.99. Underpowered, but no lift.
   - So a same-line bank that is already a donor is a donor re-weighting (equal weighting is board-validated, 204), not
     new information. At w 2 its lift sits below the flip-0.3 calibration point, where w 2 prices ≈ 0 on H1 and w 0.5
     beats it by 0.016. The Orion lab's other line (HEK293T), with the deepest independent truth we hold, loses direction
     at w 2.
   - Even as new information, cancer-line banks lift less than Pan on H1 (K562 ×1.9, HEK293 ×1.2 vs ×2.5): the ESC
     line-specific response is the strong case, as 450 §11 already said for the programme.

**What 450's K562/VIPerturb replicate supports, re-read:** that a same-line bank beats other-line banks as NEW
information (×1.9 lift on 2,517 non-essential targets here; 450 §4's 2.4× was essential targets against one other line).
It does not support Pan's size for a cancer line (×1.9 vs ×2.5 on the same unwhitened statistic). It says little about
the October case, because K562's bank is already a donor: in-mix the lift is ×1.09 unwhitened (×1.29 whitened, below).
450 §6 saw the in-mix case (+15 %) on essential targets only and did not carry it into the price.

**Verdict, decision 1: robust on H1, fragile for October.** Its H1 price survives target resampling and the
essential-free subset. The October transfer is not supported where it matters most: the likely high-coverage match
(HCT116 via Orion) gets w 2 + pt on all 300 targets with a bank that is already a third of our composition, where the
measured direction lift is ×0.95–1.1 and the H1 dose calibration says w 2 ≈ 0 or worse at that quality. The 451 bank audit
(own-target knockdown) cannot tell a new-information bank from an in-mix one.

**The per-target prior (`pt`) for an in-mix bank is not supported either (`pt487.py`).** Per target, the truth's top-100
genes by |delta| (own gene out); sign accuracy of the bank row's own sign (`pt`) against CLEAN4's pinned programme sign
(391's `psign`, 84 % down):

| truth (lab) vs in-mix bank | stratum | n | pt | CLEAN4 | all-down | pt − CLEAN4 [95 %] |
|---|---|---:|---:|---:|---:|---|
| VIPerturb K562 (Satija) vs K562-GWPS (Weissman) | non-essential | 2,090 | 0.515 | 0.532 | 0.668 | **−0.017** [−0.021, −0.014] |
| | non-essential, weaker half | 1,045 | 0.510 | 0.528 | 0.651 | −0.017 [−0.022, −0.013] |
| | essential | 410 | 0.561 | 0.533 | 0.623 | +0.028 [+0.018, +0.039] |
| XuCao HEK293 vs Orion HEK293T (Xaira) | non-essential | 78 | 0.529 | 0.491 | 0.503 | **+0.038** [+0.009, +0.070] |
| | essential | 108 | 0.566 | 0.553 | 0.612 | +0.013 [−0.015, +0.042] |

The sign of `pt − CLEAN4` flips between the two replicates (the programme is a lab signature, 450h; CLEAN4 is
Weissman's), and on non-essential targets both are within ±0.04 of chance. H1's +0.0144 for `pt` came from 30 strong
targets with a 0.764-accurate bank; nothing measured supports it on a genome-wide in-mix bank's weak targets.

**`pt` also bypasses the weak-target protection (read from the builder, `make_submission_pm_cd4.py` ~l. 5520–5555).**
The prior gate runs only when the context's programme default is on (`_prog`), and it hands CLEAN4 only to weak targets
*without* a same-line per-target prior ("weak_kept_pt"). With the October default `pt:noprog` there is no gate at all.
So an HCT116 context (Orion covers ~290/300) puts every covered target on the bank row's own signs, the design-weak
~40 % included, where 467/470 measured that weak targets want CLEAN4. The nearest H1 analogue is 467's `c4off33`
(the weak third on noisy line-own signs, CLEAN4 elsewhere): −0.0032 overall-equiv per context, all of it reach.

**Suggested rule change (analysis-only evidence; the coordinator decides):** `auto` gives w 0.5 (or rows off) when the
matched bank is one of the production donors or their same-lab siblings (`k562_gwps`, `k562_essential`, `orion_hct116`),
and keeps w 2 only for strong new-information banks (Nadig, Sunshine). For an in-mix bank, drop `pt` too (CLEAN4 stays the
block prior). The rows cannot be priced for October the way 450 priced them: H1 has no production donor from its own line.

**Settling check, run (`whitemix487.py`, 3.2 GB peak, 3 s): the same lift with a production-like WHITENED composition.**
"Ours" is 481's base-only proxy of the builder's composition (438's `Donors438.pre`: three donors, normrel cap, PM50;
the targets' own columns zeroed; whitened with the A/B/C operator `e2f65ecf` at γ 0.65 with the on-axis median preserve;
× 0.982). Known answer PASS: the A/B/C panel's pre-whitening median row norm 4.2401 = 481's `whiten481.json`. The H1
calibration is recomputed on the same construction, so the rows compare like for like:

| case (whitened ours; raw bank row, as `def_sameline` stores it) | targets | lift w 0.5 | lift w 2 [95 %] | h1-test emitter price, w 2 / w 0.5 |
|---|---:|---:|---:|---|
| H1, Pan = new information (clean) | 125 (62 non-ess.) | ×2.82 | **×4.37** [3.91, 4.96]; non-ess. ×5.42 | +0.0934 / +0.0531 (rows only) |
| H1, Pan flipped p 0.2 | | ×2.06 | ×2.82 | +0.0268 / +0.0323 |
| H1, Pan flipped p 0.3 | | ×1.69 | ×2.06; non-ess. ×2.50 | **+0.0007** / +0.0166 |
| **HEK293, Orion HEK293T in the mix** (truth XuCao) | 171 (65 non-ess.) | ×1.14 | **×1.16** [1.06, 1.26]; non-ess. ×1.24 [1.03, 1.46] | — |
| **K562, GWPS in the mix** (truth VIPerturb) | 800 non-ess. | ×1.21 | **×1.29** [1.14, 1.46] | — |

- Whitening does strip some same-line signal: the in-mix lift rises from ×0.95–1.09 (unwhitened mix) to ×1.16–1.29.
- It stays far below the weakest calibrated bank. At flip 0.3 (lift ×2.06 at w 2) the h1-test emitter already prices w 2
  at ≈ 0, with mse +0.020 raw WORSE than base, and w 0.5 at +0.0166. In-mix banks sit below that point, and their lift is
  flat from w 0.5 to w 2, so w 2 buys no direction over w 0.5 and carries the w 2 downside.
- Where the h1-test price tracks (lift − 1), at w 0.5 (flip 0.3: 0.38 of clean's lift excess → 0.31 of clean's price),
  an in-mix bank's (lift − 1) at w 0.5 is 0.08–0.12 of Pan's. That puts it near +0.004–0.006 per context at 35 %
  coverage, or ≤ +0.015 at full coverage if weak targets gained like strong ones (they did not on H1, point 1). Divided
  by three contexts: **≤ +0.005 overall for an HCT116-like context at w 0.5, and ≈ 0 or negative at w 2.**

The rule change above stands on the deepest independent same-line truth we hold (XuCao, 408 cells per perturbation),
with the production-like composition. **Remaining check, if the coordinator wants one before changing `auto`:** the
h1-test emitter price of a same-line row at the in-mix lift level. Build a synthetic bank with that lift on h1-test
(Pan's rows sign-flipped at p ≈ 0.4, whose whitened lift should read ~×1.3) and run `run450c.py` at w 0.5 / 2 with
`pt` on and off: 4 runs, ~3 min each, 7.6 GB, after 04:30Z.

Memory note: `lift487.py` peaked at 5.58 GB RSS for ~8 s (mmap pages of full-width donor rows plus a float64 copy) while
free+inactive read 29.2 GiB, above the mandate's 4 GB per-process line. Later scripts read rows one at a time (`lift487b.py`
0.33 GB).

## Decision 4 — the re-level target ρ* 0.9793 and the NTC-proxy cliff (interim, 22:40Z)

`relevel487.py` (analytic; ntclevel417's slopes: raw mse −0.31 per unit scale below the cliff, +1.78 above; one raw mse
unit in one context = 1/(0.9502 × 18) overall). **Known answer PASS:** ρ*(sd 0.01) = 0.97934, a margin of −1.044 % in scale.

| true cliff vs the NTC proxy (% scale) | shipped point vs the true cliff | raw mse per context | overall, one context | all three |
|---|---:|---:|---:|---:|
| −10 % (cliff lower) | +8.96 % past | +0.159 | −0.0093 | −0.0280 |
| −5 % | +3.96 % past | +0.070 | −0.0041 | −0.0124 |
| −2 % | +0.96 % past | +0.017 | −0.0010 | −0.0030 |
| **−0.88 % (context C's fitted f 0.9825)** | 0.16 % below | +0.0005 | −0.00003 | −0.0001 |
| 0 (proxy exact) | 1.04 % below | +0.003 | −0.0002 | −0.0006 |
| +2 % (cliff higher) | 3.04 % below | +0.009 | −0.0006 | −0.0016 |
| +5 % | 6.04 % below | +0.019 | −0.0011 | −0.0033 |
| −29 % / +29 % (481 F2: 800- or 200-cell NTC constructs) | | | | **−0.090 / −0.017** |

- Random proxy error δ ~ N(bias, sd): the regret of shipping 0.9793 against the best margin is ≤ 0.0001 per context for
  sd ≤ 1 % with no bias, ≤ 0.0003 at sd 2 %, ≤ 0.0006 at sd 3 % with a C-type bias (−0.88 %).
- If the above-cliff slope on D/E/F were 2× or 4× A/B/C's, the optimal ρ* is 0.972 / 0.966 and 0.9793 costs 0.00003 /
  0.00013 per context.
- The secant loop's ±1 % ρ tolerance (±0.5 % in scale) moves any row above by ≤ 0.0009 per context.

**Verdict: robust.** The margin covers A/B/C's worst measured proxy error (C, −0.88 % in scale), the loss is ~0.001 overall
per context per 1 % past the cliff, and no plausible slope or sd error moves the answer by more than 0.0006 per context.
The one large failure is structural and already found: 481 F2 (NTC constructs not 46 × 400) moves the proxy by ±30 %.
**Single settling check:** on Oct 22, confirm cells per `ntc_id` = cells_per_pert in each D/E/F control file (one line
of pandas), or apply 481's F2 patch, whose A/B/C known answer is bit for bit.

## Decision 6 — Arc's programme and the prior gate (interim, 23:15Z)

`prog487.py`: on H1, block-prior arms share their emission draws, so any gate rule is an exact per-target mix of 435's
CLEAN4 run and 462's ungated Arc-programme run. **Known answer PASS:** 467's gated33 equals ungated on strength ≥ 2 and
CLEAN4 on strength ≤ 1, target by target, member by member, all 3 seeds (0 mismatches).

Per strength class (467's n_gt_0.2), programme − CLEAN4 per target, h1-test's 67 non-essential targets:

| class | n | fid | reach | nmae | overall-equiv per target |
|---|---:|---:|---:|---:|---:|
| 0 | 3 | +0.026 | −0.000 | −0.019 | +0.020 |
| 1 | 23 | +0.007 | **−0.150** | −0.022 | **−0.019** |
| ≥ 2 | 41 | +0.017 | +0.003 | −0.013 | +0.012 |

Expected per-context overall-equiv (H1 terms) on each panel's class shares (0 / 1 / ≥ 2), class-stratified bootstrap:

| panel (shares) | ungated | **le1** (≤ 1 on CLEAN4) | eq0 (= q33 once > 33 % have strength ≤ 1) | le1 − eq0 |
|---|---|---|---|---|
| A/B/C (0.07 / 0.21 / 0.72) | +0.0062 [−0.001, +0.013] | **+0.0087** [+0.005, +0.013] | +0.0048 | +0.0039 [−0.002, +0.010] |
| Arc-style stand-in (0.14 / 0.26 / 0.60) | +0.0050 [−0.004, +0.013] | **+0.0073** [+0.004, +0.011] | +0.0024 [−0.006, +0.010] | +0.0049 [−0.002, +0.013] |
| efr (0.23 / 0.31 / 0.46) | +0.0041 [−0.006, +0.013] | **+0.0056** [+0.003, +0.008] | −0.0003 [−0.010, +0.009] | +0.0059 [−0.003, +0.015] |

- **le1 is the only rule whose H1 value is robustly positive on every panel,** and it degrades gently with the weak share:
  0.84× A/B/C's on an Arc-style panel, 0.64× on efr, because the weak targets carry ~0 and only the strong share
  shrinks. F4's lower delivered dose is a weak-target property (469: median λ 0.953 once they are on CLEAN4), so le1
  removes it too.
- **q33 is a different rule on the likely October panel:** it becomes eq0 and leaves −0.005 (Arc-style) to −0.006 (efr)
  on the table relative to le1 in H1 terms (CIs span 0: the class-1 reach loss rests on 4 targets).
- **What each 09-30 outcome licenses** (458g − 458d; the exact split makes 458g − 446b the programme on the strong
  two-thirds, and 458d − 458g the programme on the weak third):
  - gate transfers (≥ +0.0015) → programme + **le1** (not q33). Expected D/E/F value ≈ (458g − 446b) × 0.84 on an
    Arc-style panel (× 0.64 if efr-like).
  - null (±0.0005) → the prereg ships the ungated programme. On D/E/F the strength-1 class is 1.3–1.5× A/B/C's share, so
    any per-target cost the A/B/C read hides grows with it. le1 is still the lower-variance choice (CI above 0 on every
    panel, H1 terms). Recommend le1 there too, unless 458d − 458g is clearly positive.
  - harm (≤ −0.0005) → no gate; the ungated programme only if 458d itself read positive.
  - 458d null/harm and 458g − 446b fid ≥ +0.0015 → programme + le1, scaled as above.
- **Verdict: fragile in its rule, robust once le1 replaces q33.** The one board read decides the programme itself; its
  transfer to D/E/F depends on the gate's meaning, which le1 fixes (472, in flight). **Single settling check:** 472's
  KA-le1 (le1 rebuilds 458g byte for byte), then `DEF_DAY_PRIOR_GATE=le1` in the runbook's decision line.

## Decision 3 — sharpening q 1.75 (interim, 23:25Z)

Evidence: board steps 417a (q 1.25, +0.0070), 433k (q 1.5 levelled, +0.0046), 446b (q 1.75 levelled, +0.0023), all
bulk-only (pds/mse reproduce to 16 s.f., so the A/B/C curve is measured without ordering noise); 438's PLOLO on
October-like panels with 417a's calibration (k 0.65 mse, k_pds 0.42). Per-fold increments re-derived from
`438/work/summary438.json` (levelled, calibrated):

| panel | q 1.5 → 1.75, per fold (K562 / HCT116 / HEK293T) | mean | q 1.75 → 2, per fold | mean |
|---|---|---:|---|---:|
| abc | +0.0015 / +0.0022 / +0.0027 | +0.0021 | −0.0003 / +0.0005 / +0.0008 | +0.0003 |
| real | +0.0003 / +0.0010 / +0.0010 | +0.0007 | −0.0007 / −0.0001 / −0.0007 | −0.0005 |
| efp | +0.0013 / −0.0000 / +0.0008 | +0.0007 | +0.0002 / −0.0008 / +0.0001 | −0.0002 |
| efr | +0.0009 / +0.0013 / +0.0017 | +0.0013 | +0.0003 / +0.0004 / +0.0005 | +0.0004 |

- **Known answer (the calibrated instrument retrodicts two board steps it was not fitted on):** abc 1.25 → 1.5 +0.0044
  vs board 433k +0.0046; 1.5 → 1.75 +0.0021 vs board 446b +0.0023.
- On October-like panels q 1.75 beats q 1.5 on 8 of 9 fold-panels (the ninth is −0.00002), and q 2 is a wash (−0.0005 …
  +0.0004). The optimum is 1.75–2 everywhere; the regret of 1.75 is ≤ 0.0004 on any panel measured.
- Not in the PLOLO price: the ballast's nmae cost, which grows with q and with the re-level factor (446b: +0.0013 raw
  nmae inside its +0.0023; October factors 1.21–1.29 vs A/B/C's 1.10–1.17). It is in the board number for A/B/C and small.

**Verdict: robust.** The optimum is not panel-specific at the precision that matters, the instrument retrodicts the
board twice, and the sign of the 1.5 → 1.75 step holds on every October-like panel. No check needed.

## Decision 2 — the October stack H1 ×2 / CD4 ×1 (interim, 23:25Z)

Evidence: 438's levelled PLOLO surfaces (3 cancer folds) on abc / real / efp / efr. Regret of each CD4 dose at H1 ×2
against the best dose on each panel (levelled k = 1, from `summary438.md`; calibrated values are ~0.4–0.65× these):

| CD4 at H1 ×2 | abc | real | efp | efr | max regret |
|---|---:|---:|---:|---:|---:|
| ×0.75 | 0.0039 | 0.0001 | 0.0004 | 0 | 0.0039 |
| **×1** | 0.0017 | 0 | 0 | 0.0007 | **0.0017** |
| ×1.25 | 0.0006 | 0.0005 | 0.0004 | 0.0018 | 0.0018 |
| ×1.5 (A/B/C's) | 0 | 0.0016 | 0.0011 | 0.0031 | 0.0031 |

- The CD4 optimum is panel-specific: ×1.5 on A/B/C (3/3 folds), ×1 or lower on every October-like draw (0/3 folds
  positive for ×1.5 on efp and efr). An Arc-style panel sits between A/B/C and efr in design strength (481 F3), so ×1
  (or ×1.25, within 0.0001 of it on max regret) is the right minimax dose. Stakes ≤ 0.002.
- The H1 axis is under-dosed on essential-screen-free subsets: 438's Arc-like column prices H1 ×2 → ×3 at +0.0026 (real),
  +0.0015 (efp), +0.0005 (efr), but the step fails leave-top-3-out on efp and efr. Stakes ≤ 0.002 (calibrated ≤ 0.0015).
- 471 changes the H1 union itself (h1_2025_full); the H1 ×2 price above is for the old H1 pair.
- The 09-29 read (457dr CD4 ×1.5 vs 453dr ×1, October path on A/B/C) tests A/B/C's panel, where ×1.5 is expected to
  win (+0.0017 k = 1). A win validates the instrument on the CD4 axis. It does not license ×1.5 on D/E/F.

**Verdict: robust in expectation, panel-specific in its optimum; stakes ≤ ~0.002 on either axis.** The instrument that
would pick the day's weights exists (`dayof438.py`, 27–35 min, ≤ 8.8 GB) but is not in `def_day.sh`, and its q grid
stops at 1.5 (481 item 9). **Single settling check:** extend `dayof438.py`'s q grid to 1.75 and rerun it on efp (its
rehearsal panel); if it reproduces the chain's efp surface, run it on Oct 22 between setup and the dump.
