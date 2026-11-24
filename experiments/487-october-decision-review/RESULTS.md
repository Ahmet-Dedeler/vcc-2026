# 487 — Adversarial robustness review of the evidence behind the October recipe

**Answer (2026-09-27, analysis only; every number below is re-derived from saved outputs, each with a known-answer row).**
Ranked by expected October impact × fragility:

| # | decision | verdict | what moves the final score | single settling check |
|---|---|---|---|---|
| 1 | **Same-line rows, `auto` w 2 + `pt` (451/450)** | **robust on H1, fragile to unsupported for October** | The +0.0934 is H1 (ESC; the bank is new information; 30 strong targets). The likely high-coverage matches (HCT116 via Orion, K562 via GWPS) are banks that are **already production donors**. There the measured direction lift is ×1.16–1.29 with a production-like whitened composition, against ×4.37 for Pan on H1, and it sits below the weakest bank 450f calibrated (flip 0.3, ×2.06, where w 2 already prices ≈ 0 with mse +0.020 raw worse). Even as new information, cancer-line banks lift only ×1.62 (HEK293) / ×2.08 (K562). Realistic value: ≤ +0.005 overall at w 0.5, ≈ 0 or negative at w 2. `pt` on an in-mix bank is ±0.02–0.04 against CLEAN4 (sign flips between replicates), and `pt` bypasses the prior gate's weak-target protection. | No code change: `def_day.sh` args 5–6 = `0.5` and `noprog` (or `""` with the programme source): w 0.5 for every cancer-line bank (changes only HCT116 in practice) and no `pt`. Rehearse once (452b-style, forced HCT116). A finer per-bank patch is ready (`patches/in_mix_sameline.diff`, flag-off byte for byte). Open: an h1-test emitter run below flip 0.3 (`inmix487.py`, after 04:30Z). |
| 2 | **Prior gate q33 (468) with Arc's programme (458d)** | **fragile as q33, robust as le1** | On an Arc-style panel q33 degrades to "strength 0 only" (481 F3). On H1's non-essential targets that leaves −0.005 relative to le1 (CI spans 0). le1 is the only rule whose H1 value is positive with its CI above 0 on every panel: +0.0087 A/B/C, +0.0073 Arc-style, +0.0056 efr. It transfers 0.84× (Arc-style) / 0.64× (efr) of A/B/C's strong-target effect. | 472's KA-le1 passed (22:46Z). Set `DEF_DAY_PRIOR_GATE=le1` whatever the 09-30 bin. D/E/F value ≈ (458g − 446b) × strong-share ratio. |
| 3 | **Stack H1 ×2 / CD4 ×1 (438)** | **robust in expectation, panel-specific optimum** | CD4 ×1 has max regret 0.0017 across four panels (×1.5: 0.0031). The H1 axis is under-dosed on Arc-like subsets (×2 → ×3 +0.0005 … +0.0026) but fails leave-top-3-out. Stakes ≤ 0.002. The 09-29 read (457dr vs 453dr) tests A/B/C's panel, not D/E/F's. | Extend `dayof438.py`'s q grid to 1.75, re-run on efp, and wire it into the day (27–35 min, ≤ 8.8 GB). |
| 4 | **LODO key (415/460)** | **robust on A/B/C (board key_ok 09-27), licensed for October; moderately fragile in its October size** | Board 446bj2 − 446b = **+0.0022** overall (bin key_ok): mse −0.0145 raw (above the predicted band), pds +0.0002 (5–20× below the predicted band). It transfers as an mse lever. October ≈ +0.0014 … +0.0026 by 438's panel ratios. efr's share rests on one fold of three. | `DEF_DAY_LODO=1`. Read 460dr vs 453dr (the October path's own fit, already built) on its mse. |
| 5 | **Re-level ρ* 0.9793 (447/448)** | **robust, with one known trigger** | ~0.001 overall per context per 1 % past the cliff. The −1.04 % margin covers A/B/C's worst proxy error (C, −0.88 %). The regret of 0.9793 vs the best margin is ≤ 0.0006 per context for sd ≤ 3 %, and ≤ 0.00013 for a 4× slope error. The trigger is 481 F2 (NTC constructs ≠ 46 × 400): −0.017 to −0.09 overall. | Count cells per `ntc_id` in each D/E/F control file on Oct 22, or apply F2's patch. |
| 6 | **q 1.75 (417/433/446)** | **robust** | The calibrated PLOLO retrodicts both later board steps (+0.0044 vs +0.0046; +0.0021 vs +0.0023). On October-like panels q 1.5 → 1.75 is positive on 8/9 fold-panels, and q 2 is a wash (−0.0005 … +0.0004). The regret of 1.75 is ≤ 0.0004 anywhere. | none |

**The headline:** the same-line mode is booked as "likely the largest October lever" (+0.09 at 35 % coverage). For the
matches that give it coverage, it is a re-weighting of a donor we already have. The +0.09 is not a planning number for
October. Its defaults (w 2, `pt` on every covered target, no weak-target gate under `noprog`) are the ones the data
here argue against.

Scope and rules kept: no builds, no `def_day.sh`, no `vcc`, no `cell-eval2`, no git, no edits outside this directory. One
memory deviation, reported under Files.

---

## 1. Same-line rows at w 2 + per-target prior (mandate decision 1; 450, 451, `auto` since 451b)

**Decision and evidence.** A context whose identity call matches a registered public bank blends that bank's per-target
row into our rows: unit(ours) + w·unit(bank), at our norm. It also takes the row's own signs as the block prior (`pt`).
`auto` sets w 2 for banks with own-target knockdown in ≥ 0.9 of rows (Orion HCT116 0.911, K562-GWPS 0.990), else 0.5.
Evidence: 450 §9/§10 on h1-test, Pan 2025 hESC as the bank, 30 of 85 targets covered: +0.0934 at w 2 + pt.
`scripts/def_sameline.py` has the registry and `STRONG_BANKS`.

### 1a. Known answer and H1 robustness (`sameline487.py`, `sameline487b.py`; 112 MB)
From the per-target CSVs and mse legs, paired by seed with 435's SLOPE, every row of 450 §9/§10 reproduces: full w2
+0.0934 ± 0.0002, w1 +0.0837, w4 +0.0933 ± 0.0005, rows-only w2 +0.0790, flip 0.2 w0.5/w1/w2 +0.0323/+0.0316/+0.0268,
flip 0.3 +0.0166/+0.0086/+0.0007. Uncovered targets contribute exactly 0 on every member.

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

The H1 price survives target resampling and the essential-free subset. Four things the record does not state:
1. **H1's covered targets are strong responders.** n_real median 986 for the 30 covered targets; 13 % below 50. Of
   h1-test's 67 non-essential targets, 54 % are below 50 and 39 % below 10. Only 4 covered targets are weak (SLC5A6,
   STUB1, DPF2, TOR1A; per-target gains ×85 = +0.07, −0.15, +0.17, +0.22). Per-target gain rises with strength
   (Spearman with n_real +0.46). A genome-wide October bank covers the weak half too, where H1 has almost no evidence.
   The coverage curve's +0.20 at 100 % coverage (resampling covered targets) is therefore an upper figure.
2. **w 2 is flat in w and on an edge in bank quality.** At Pan quality w 2 ≈ w 4 and beats w 1 by +0.010. At ≈0.61 sign
   quality w 2 loses 0.016 to w 0.5. `auto` classes banks by own-target knockdown, which does not measure a bank's
   agreement with another lab's truth.
3. **The +0.0934 is one context.** One matched context of three moves the overall by a third of it.
4. **It has no H1 → board calibration.** It is in H1-scaled units. No board point exists or can exist before October
   (A/B/C's same-line banks cover 0 of its essential-free panel). `docs/H1-BOARD-CALIBRATION.md` has no row for a bulk
   direction change. The one board test of a cancer-line same-line bank is 446bja (Nadig Jurkat's programme on A, a
   different object from the rows). It carried 0.07× of its H1-priced fid.

### 1b. The October case: the high-coverage banks are already our donors (`lift487.py`, `whitemix487.py`)
451's coverage table on an October-like panel: K562 300/300 and HCT116 300/300, both "strong" (w 2). Every
new-information bank covers ≤ 83/300 (DLD-1, a sparse DE-lfc table, weak class) or 2–9 (HepG2, A549, MCF7). K562-GWPS and
Orion HCT116 are two of our three production donors, and HEK293T (Orion's other line) is the third.

Statistic: direction lift = mean cosine of the blend to an independent same-line truth ÷ mean cosine of "ours"
(targets centred; ratio invariant to independent truth noise; target bootstrap). Calibrated on H1 with the same
construction against the emitter prices 450f measured on h1-test.

| case | "ours" | targets | lift w 0.5 | lift w 2 [95 %] | h1-test price w 2 / w 0.5 |
|---|---|---:|---:|---:|---|
| H1, Pan = new information | whitened | 125 (62 non-ess.) | ×2.82 | **×4.37** [3.91, 4.96]; non-ess. ×5.42 | +0.0934 / +0.0531 (rows only) |
| H1, Pan flipped p 0.2 | whitened | | ×2.06 | ×2.82 | +0.0268 / +0.0323 |
| H1, Pan flipped p 0.3 | whitened | | ×1.69 | ×2.06; non-ess. ×2.50 | **+0.0007** (mse +0.020 raw worse) / +0.0166 |
| **HEK293, Orion HEK293T already in the mix** (truth XuCao, 408 cells/pert) | whitened | 171 (65 non-ess.) | ×1.14 | **×1.16** [1.06, 1.26]; non-ess. ×1.24 [1.03, 1.46] | — |
| **K562, GWPS already in the mix** (truth VIPerturb, Satija) | whitened | 800 non-ess. | ×1.21 | **×1.29** [1.14, 1.46] | — |
| HEK293 as NEW information (ours = K562 + HCT116) | whitened | 171 | ×1.44 | **×1.62** [1.43, 1.83]; non-ess. ×1.78 | — |
| K562 as NEW information (ours = HCT116 + HEK293T) | whitened | 800 non-ess. | ×1.69 | **×2.08** [1.67, 2.70] | — |
| H1, Pan (new information) | unwhitened 3-donor mix | 125 | ×1.87 | ×2.50 [2.26, 2.80] | |
| K562, GWPS if it were new information | unwhitened | 3,000 | ×1.60 | ×1.92 | |
| K562, GWPS in the mix | unwhitened | 2,517 non-ess. | ×1.10 | ×1.09 [1.03, 1.16] | |
| HEK293, Orion HEK293T if new information | unwhitened | 171 | ×1.22 | ×1.23 | |
| HEK293, Orion HEK293T in the mix | unwhitened | 171 | ×1.01 | ×0.95 [0.90, 0.99] | |

- "Whitened" is 481's base-only proxy of the builder's composition: 438's `Donors438.pre` (three donors, normrel cap,
  PM50), the targets' own columns zeroed, whitened with `e2f65ecf` at γ 0.65 with the on-axis median preserve, × 0.982.
  Known answer PASS: the A/B/C panel's pre-whitening median row norm 4.2401 equals `481/work/whiten481.json`.
  Consistency with the harness that priced +0.0934: on h1-test's 30 covered targets 454b measured the harness's own
  composition (cap15) at cosine 0.068 → 0.197 → 0.307 for w 0 / 0.5 / 2, a lift of ×2.9 / ×4.5. The whitened proxy
  reads ×2.82 / ×4.37 on 125 targets. The bank row is the raw delta, as `def_sameline` stores it. Uncentred cosines give the same ordering (in-mix ×1.24 / ×1.58,
  Pan ×4.44, flip 0.3 ×2.09).
- By truth strength (halves by the truth row's centred energy), the in-mix lift at w 2 is flat: HEK293 ×1.16 (top) /
  ×1.20 (bottom), K562 ×1.36 / ×1.16. Uncentred, K562 reads higher (×1.58; top half ×1.66). GWPS and VIPerturb share a
  K562 mean response that the whitened composition lacks. On H1 the bank's mean ran against the truth's (450 §7, cos
  −0.12), so the centred numbers are the ones to plan on.
- With the production dump itself (465's efr panel, 277 VIPerturb-covered targets, `lift487b.py`), VIPerturb has almost
  no signal on those non-essential targets: w 2 ×0.85 [0.54, 1.15]. The check is underpowered but shows no lift.
- Whitening strips some same-line signal (the in-mix lift rises from ×0.95–1.09 to ×1.16–1.29). It stays far below
  the weakest calibrated bank, and it is flat from w 0.5 to w 2. So w 2 buys no direction over w 0.5 and carries w 2's
  downside: at flip 0.3 the h1-test price of w 2 is already ≈ 0, with mse worse than base.
- At w 0.5 the h1-test price tracks the lift excess (lift − 1) roughly proportionally: flip 0.2 is 0.58 of clean's excess
  → 0.61 of clean's price; flip 0.3 is 0.38 → 0.31. An in-mix bank's excess at w 0.5 is 0.08–0.12 of Pan's. That is
  ~+0.004 per context at 35 % coverage, ≤ +0.015 at full coverage if weak targets gained like strong ones (point 1 says
  they gain less). So ≤ +0.005 overall, in H1 terms, for one HCT116-like context at w 0.5.
- **Even as new information, cancer-line banks lift far less than Pan:** whitened ×2.08 (K562) and ×1.62 (HEK293) at w 2,
  against Pan's ×4.37. Both sit at or below the flip-0.3 calibration point (×2.06), where the h1-test emitter prices
  w 0.5 at +0.0166 and w 2 at +0.0007 (w 2 − w 0.5 = −0.016 [−0.029, −0.005]). Their direction lift is also flat from
  w 1 to w 2 (×1.97 → ×2.08; ×1.59 → ×1.62). So the w 2 dose, calibrated on an ESC bank, is not supported for any
  cancer-line bank, new information or not. The ESC line-specific response is the strong case, as 450 §11 found for
  the programme.
- **What 450's K562/VIPerturb replicate supports:** a same-line bank beats other-line banks as new information (×1.9 on
  2,517 non-essential targets here; 450 §4's 2.4× was essential targets against one other line). It does not support
  Pan's size for a cancer line, and it does not cover the October case, where K562's bank is a donor. 450 §6 measured
  the in-mix case (+15 %) on essential targets only, and that number was never carried into the price.

**Registry coverage on an Arc-style panel** (the 2025 H1 panel's 222 non-essential targets, Arc's only other panel; bank
perturbation lists, `def_sameline.REGISTRY`):

| matched line | bank(s) | in our mix? | `auto` w | Arc-style coverage | A/B/C coverage |
|---|---|---|---:|---:|---:|
| HCT116 | Orion | **yes** | 2 | **0.96** | 290/300 |
| K562 | GWPS (+ VIPerturb, K562-ess.) | **yes** | 0.5 (VIPerturb supplies 22/222 and 30/300 rows, a weak bank) | **1.00** | 299/300 |
| DLD-1 | low-MOI DE-lfc table (sparse) | no | 0.5 | 0.27 | 67/300 |
| HepG2 / Jurkat | Nadig (essential screens) | no | 2 | 0.09 | 0 |
| A549 / HT29 / MCF7 / BxPC3 / HAP1 | Jiang (+ Myllia) | no | 0.5 | 0.09–0.13 | 9–22/300 |
| Calu-3 / THP-1 | Sunshine / Wessels | no | 2 / 0.5 | 0.02–0.03 | 1–6/300 |

Every high-coverage match is an in-mix bank. `auto` already gives K562 w 0.5, because VIPerturb (weak class) fills the
targets GWPS lacks. So the w question is HCT116's alone; `pt` applies to both. The new-information banks cover ≤ 27 %
(a sparse, weak table) or ≤ 13 %. Even at Pan's quality, a 9 % strong bank scales H1's +0.0934 (at 35 %) to ~+0.024
per context, ≤ +0.008 overall. The cancer-line new-information lifts measured above (whitened ×1.6–2.1 vs Pan's ×4.37)
cut that further. No registry scenario reaches the +0.03 overall that the H1 price suggests for one matched context.

### 1c. `pt` on an in-mix bank (`pt487.py`, 0.39 GB)
Per target: the truth's top-100 genes by |delta| (own gene out). Compared: sign accuracy of the bank row's own sign
(`pt`), CLEAN4's pinned programme sign (391's `psign`, 84 % down), and all-down.

| truth vs in-mix bank | stratum | n | pt | CLEAN4 | all-down | pt − CLEAN4 [95 %] |
|---|---|---:|---:|---:|---:|---|
| VIPerturb K562 vs K562-GWPS | non-essential | 2,090 | 0.515 | 0.532 | 0.668 | **−0.017** [−0.021, −0.014] |
| | non-essential, weaker half | 1,045 | 0.510 | 0.528 | 0.651 | −0.017 [−0.022, −0.013] |
| | essential | 410 | 0.561 | 0.533 | 0.623 | +0.028 [+0.018, +0.039] |
| XuCao HEK293 vs Orion HEK293T | non-essential | 78 | 0.529 | 0.491 | 0.503 | **+0.038** [+0.009, +0.070] |
| | essential | 108 | 0.566 | 0.553 | 0.612 | +0.013 [−0.015, +0.042] |

The sign of pt − CLEAN4 flips between the two replicates, which fits 450h: a programme is a lab signature, and CLEAN4
is Weissman's. On non-essential targets everything is within ±0.04 of chance. H1's +0.0144 for `pt` came from 30 strong
targets and a bank with 0.764 sign accuracy.

**`pt` bypasses the weak-target protection** (builder, `make_submission_pm_cd4.py` ~l. 5520–5555).
- The prior gate runs only when the programme default is on (`_prog`).
- Even then it gives CLEAN4 only to weak targets *without* a same-line prior (`weak_kept_pt`).
- Under the October default `pt:noprog` there is no gate at all.
- So an HCT116 context (Orion covers ~290/300) puts every covered target on the bank's own signs, the ~40 % design-weak
  targets included. 467/470 found that weak targets want CLEAN4.
- The nearest H1 analogue, 467's `c4off33` (weak third on line-own noisy signs), costs −0.0032 overall-equiv per context,
  all reach.

### 1d. Verdict and recommendation
**Robust on H1, fragile to unsupported for October.** The most likely high-coverage match, HCT116 through Orion, gets
w 2 + `pt` on ~290 targets. Its bank is already a third of our composition, and w 2 + `pt` is the setting the evidence
here argues against. Own-target knockdown, the quantity `auto` reads, cannot tell a new-information bank from an in-mix
one, and it cannot tell an ESC bank's quality (×4.37 lift) from a cancer bank's (×1.6–2.1 as new information).

**Suggested rule change (the coordinator decides). It needs no code change:**
- **w 0.5 for every cancer-line bank.** Pass a numeric w to `def_day.sh` (arg 5 = `0.5`; "a number forces one global w").
  Every October match is a cancer line, and the measured cancer-line lifts (×1.2–2.1 whitened) sit at or below the bank
  quality where w 0.5 beat w 2 on h1-test. In practice this changes only HCT116 (Orion, w 2 → 0.5): K562 already gets
  0.5, and the other strong banks (Nadig, Sunshine) cover ≤ 9 % of an Arc-style panel.
- **Drop `pt`.** Pass arg 6 = `noprog` without the programme source, or arg 6 = `""` with `DEF_DAY_PROGRAMME_SOURCE`
  (`def_day.sh` builds `VCC_SAMELINE=<npz>:<sha>:<w>${SLF:+:$SLF}`, and empty flags pass its `noprog` check).
  - The bank's own signs are ±0.02–0.04 against CLEAN4, and the sign of the difference flips between the two cancer
    replicates.
  - Without `pt`, the prior gate also covers the matched context's weak targets.
- Rehearse once before Oct 22: a 452b-style forced-HCT116 `def_day.sh` on efr with `0.5` and `noprog`. The rehearsals so
  far (452, 452b, 465, 469) all ran with `pt`.
- Plan October as if the same-line mode is worth ≤ +0.005 overall, not +0.03.
- **Finer-grained alternative, patch ready, not applied:** `patches/in_mix_sameline.diff` adds
  `def_sameline.py --in-mix-weak`, which applies w 0.5 and no per-target signs to in-mix banks only. Known answers:
  - flag off rebuilds 465's `sameline.npz` byte for byte (`b4eb2a6d…`);
  - flag on with a non-in-mix match changes provenance only;
  - forced HCT116 changes exactly `w_D` (2 → 0.5) and `signs_D` (zeroed).
  - A one-line builder follow-up (skip `pt` entries whose signs are all zero) keeps the prior gate's CLEAN4 hand-off for
    those targets when the programme is on. See `patches/README.md`.

**Remaining check** (after 04:30Z, ~3 min and 7.6 GB per run, serial; `inmix487.py` is written and writes only here):
- the h1-test emitter price of a same-line bank below the flip-0.3 point, where the in-mix banks sit;
- Pan's rows sign-flipped at p 0.4 / 0.45, through 450c's blend re-implemented in `inmix487.py`, at w 0.5 / 2, `pt` on and
  off. Known answer first: `inmix487.py 2 0.3 pt` must reproduce 450's `h1c_w2_clean4_f0.3s0` exactly.
- Caveat: flipped rows are independent noise, while an in-mix bank shares noise with ours, so this is a proxy.

---

## 2. Arc's programme and the prior gate (mandate decision 6; 458d/458g/468/472; the 09-30 pair)

**Decision and evidence.** `DEF_DAY_PROGRAMME_SOURCE=h1_2025_full` + `DEF_DAY_PRIOR_GATE=q33`, decided by 458d / 458g
on 09-30 (`468-prior-gate/PREREGISTRATION.md`). H1: 467; the target bootstrap is `467/bootstrap467.py`.

`prog487.py` (0.1 GB): block-prior arms share their emission draws, so any gate rule is an exact per-target mix of 435's
CLEAN4 run and 462's ungated Arc-programme run. **Known answer PASS:** 467's gated33 equals ungated on strength ≥ 2 and
CLEAN4 on strength ≤ 1, target by target, member by member, all 3 seeds (0 mismatches).

Per strength class (467's n_gt_0.2), programme − CLEAN4 per target, h1-test's 67 non-essential targets:

| class | n | fid | reach | nmae | overall-equiv per target |
|---|---:|---:|---:|---:|---:|
| 0 | 3 | +0.026 | −0.000 | −0.019 | +0.020 |
| 1 | 23 | +0.007 | **−0.150** | −0.022 | **−0.019** |
| ≥ 2 | 41 | +0.017 | +0.003 | −0.013 | +0.012 |

Expected per-context overall-equiv (H1 terms) on each panel's class shares (0 / 1 / ≥ 2), class-stratified bootstrap:

| panel (shares; source) | ungated | **le1** (≤ 1 → CLEAN4) | eq0 (= q33 once > 33 % have strength ≤ 1) | le1 − eq0 |
|---|---|---|---|---|
| A/B/C (0.07 / 0.21 / 0.72; 453 builder) | +0.0062 [−0.001, +0.013] | **+0.0087** [+0.005, +0.013] | +0.0048 | +0.0039 [−0.002, +0.010] |
| Arc-style (0.14 / 0.26 / 0.60; 481 knn, 2025 H1 non-ess.) | +0.0050 [−0.004, +0.013] | **+0.0073** [+0.004, +0.011] | +0.0024 [−0.006, +0.010] | +0.0049 [−0.002, +0.013] |
| efr (0.23 / 0.31 / 0.46; 465 builder) | +0.0041 [−0.006, +0.013] | **+0.0056** [+0.003, +0.008] | −0.0003 [−0.010, +0.009] | +0.0059 [−0.003, +0.015] |

- le1 is the only rule whose H1 value has its CI above 0 on every panel. It degrades gently with the weak share:
  0.84× A/B/C's on an Arc-style panel, 0.64× on efr. The weak targets carry ~0 and only the strong share shrinks.
- F4's lower delivered dose is a weak-target property on A/B/C (469: median λ 0.953 once they are on CLEAN4), so le1
  removes most of it.
- q33 is a different rule on the likely October panel. It becomes eq0 and gives up 0.005 (Arc-style) to 0.006 (efr)
  against le1 in H1 terms. The CIs span 0: the class-1 reach loss rests on 4 targets.
- **What each 09-30 outcome licenses.** The exact split makes 458g − 446b the programme on the strong two-thirds and
  458d − 458g the programme on the weak third.
  - Gate transfers (458g − 458d ≥ +0.0015) → programme + **le1**, not q33. D/E/F value ≈ (458g − 446b) ×
    (strength ≥ 2 share on the day ÷ 0.723). Read the share from the build's `de_prior_gate` audit or
    `481/transfer_check481.py`.
  - Null (±0.0005) → the prereg ships the ungated programme. D/E/F has 1.3–1.5× A/B/C's strength-1 share, so any
    per-target cost the A/B/C read hides grows with it. le1 is the lower-variance choice here too (CI above 0 on every
    panel), unless 458d − 458g is clearly positive.
  - Harm (≤ −0.0005) → no gate; the ungated programme only if 458d itself read positive.
  - 458d null/harm but 458g − 446b fid ≥ +0.0015 and overall ≥ 0 → programme + le1, scaled as above.
- A matched same-line context bypasses this gate for its covered targets (§1c).
- **Sensitivity to the one read.** Prior-only arms share their emission draws (468: 0 of 900 blocks equal to neither
  parent), so 458g − 458d carries no re-roll noise. D/E/F will be a different realization, but the programme's effect
  moves little across realizations: H1's per-arm seed sd is 0.0004–0.0007 overall-equiv (467). What does not carry is
  the panel (above) and the fid transfer (0.07–0.67× across DE sign priors, `docs/H1-BOARD-CALIBRATION.md`), and the
  read measures the transfer directly.

**Verdict: fragile as q33, robust as le1.** The board read decides the programme itself. How it transfers to D/E/F
depends on what the gate means on a weaker panel, and le1 pins that. **Single settling check:** 472's KA-le1 (le1
rebuilds 458g byte for byte) **PASSED at 22:46Z** (`472/work/chain472.log`, sha `9929218c…`; KA-1 PASS 23:16Z). What
remains is `DEF_DAY_PRIOR_GATE=le1` in DEF-RUNBOOK's decision line, whichever bin 458g lands in.

---

## 3. The October stack H1 ×2 / CD4 ×1 (mandate decision 2; 438, `dayof438.py`)

Evidence: 438's levelled PLOLO surfaces (3 cancer folds) on abc / real / efp / efr (`438/work/summary438.md`,
`robust438.md`). **Known answer PASS:** the per-fold values in `summary438.json` average to 438's recorded CD4 ×1.5 step
(abc +0.00170, real −0.00160, efp −0.00111, efr −0.00236) and to its "×1: max regret 0.0017, mean 0.0006". Regret of each
CD4 dose at H1 ×2 against each panel's best dose (levelled k = 1; calibrated ~0.4–0.65×):

| CD4 at H1 ×2 | abc | real | efp | efr | max regret |
|---|---:|---:|---:|---:|---:|
| ×0.75 | 0.0039 | 0.0001 | 0.0004 | 0 | 0.0039 |
| **×1** | 0.0017 | 0 | 0 | 0.0007 | **0.0017** |
| ×1.25 | 0.0006 | 0.0005 | 0.0004 | 0.0018 | 0.0018 |
| ×1.5 (A/B/C's) | 0 | 0.0016 | 0.0011 | 0.0031 | 0.0031 |

- **The CD4 optimum is panel-specific.** ×1.5 wins on A/B/C (3/3 folds). ×1 or lower wins on every October-like draw
  (×1.5 is 0/3 on efp and efr). An Arc-style panel sits between A/B/C and efr in design strength (481 F3), so ×1 (or
  ×1.25, within 0.0001 on max regret) is the minimax dose.
- The unknown is whether A/B/C's preference comes from Arc's selection process (then D/E/F repeats it) or from its
  particular targets. 438's mechanism: CD4 rows point 0.89× as well as the base on A/B/C, 0.52–0.79× on the draws.
- **The H1 axis is under-dosed on essential-screen-free subsets.** 438's Arc-like column prices H1 ×2 → ×3 at +0.0026
  (real), +0.0015 (efp), +0.0005 (efr). The step fails leave-top-3-out on efp and efr. 471 changes the H1 union
  itself (h1_2025_full); these prices are for the old H1 pair.
- The 09-29 read (457dr CD4 ×1.5 vs 453dr ×1, October path on A/B/C) tests A/B/C's panel, where 438 expects ×1.5 to win
  (+0.0017 at k = 1). A win validates the instrument on the CD4 axis. It does not license ×1.5 on D/E/F.

**Verdict: robust in expectation; the optimum is panel-specific on both axes; stakes ≤ ~0.002.** The day-of re-pricer
exists (`dayof438.py`, 27–35 min, ≤ 8.8 GB, KA-exact on efp) but is not called by `def_day.sh`, and its q grid stops at
1.5 (481 item 9). **Single settling check:** extend its grid to 1.75 and re-run it on efp. If it reproduces the chain's
efp surface, run it on Oct 22 between setup and the dump. Otherwise ship ×2 / ×1 as planned.

---

## 4. The LODO key (mandate decision 5; 415, 460; board read of 446bj2)

Evidence: 415's PLOLO (V1 beats 410's key on 3/3 cancer folds and kolf2; 7 of 8 target-bootstrap intervals exclude 0).
438 on the stack, levelled k = 1: abc +0.0027 [+0.0016, +0.0041], real +0.0032, efp +0.0031, efr +0.0017, all 3/3 or
2/3 folds. Leave-top-3-out gives +0.0021 / +0.0021 / +0.0021 / +0.0010, and the Arc-like subset +0.0041 / +0.0054 /
+0.0017. The day-of fit (460) reproduces A/B/C's `8c82bef3` to 2e-14, and October rules move w by at most 0.0067.

- **Transfer to D/E/F:** October-like panels read 0.63× (efr) to 1.17× (real) of A/B/C's PLOLO value. A board read of +X
  on A/B/C therefore implies ~0.6–1.2 X on D/E/F. efr is lowest because more of its targets have shared energy ≤ 0
  (floored: 14.6 → 27.3 %, 481).
- **Noise on the read:** 446bj2 and 446b differ in the bulk key only, and the ballast pins the DE channel (410b vs 395: DE
  members within 0.0007 raw). pds/mse reproduce exactly, so the key_ok bin (overall ≥ +0.001, pds ≥ 0) is not exposed to
  the 0.0045 ordering sd.

### 4b. The board read (446bj2, entry `YvVru3DAF89Jiyr34olR`, captured 00:27Z 09-27; `read487.py`)
Known answer for the reader: run on 433k's record against 446b, it returns 446b's recorded step exactly (pds +0.0032,
mse −0.0076 raw, nmae +0.0013, overall +0.0023).

| member (raw Δ vs 446b) | board | prereg (433 PREREGISTRATION-0926) |
|---|---:|---|
| **overall** | **+0.00223** (0.2722766 vs 0.2700428) | +0.001 … +0.003 ✓ |
| pds | +0.00021 | +0.001 … +0.004 ✗ (5–20× below) |
| mse | **−0.01447** | −0.004 … −0.012 ✗ (above the top of the band) |
| reach / nmae / fid / jac | −0.0020 / +0.00006 / −0.00002 / +0.00002 | within the ballast's drift ✓ |
| bulk members only (pds + mse, scaled / 6) | +0.00262 | |

- **Bin: key_ok** (overall ≥ +0.001 with pds ≥ 0). By the preregistered rule, `DEF_DAY_LODO=1` for October.
- **It transferred as an mse lever, not a pds lever.** The overall landed in its band by compensation: mse 1.2–3.6× its
  prediction, pds 0.05–0.2×. This matches the record for per-target scaling: 410/410b's key moved mse and mis-read pds
  (`docs/H1-BOARD-CALIBRATION.md` rules 1 and 5). A PLOLO composite for a per-target key should be read for its mse part
  only.
- **The bulk part is exact.** pds/mse reproduce to 16 figures, so +0.0026 is this panel's value with no ordering noise.
  The DE members moved −0.0004 overall (reach −0.002 raw), inside the ballast drift the prereg allowed.
- **October size:** 438's PLOLO puts October-like panels at 0.62–1.16× A/B/C's step, i.e. +0.0014 … +0.0026 if the ratio
  carries. It is fragile on the weakest panel: efr's +0.0017 rests on one fold (K562 +0.0005, HCT116 +0.0003, HEK293T
  +0.0043), while efp is uniform (+0.0027 … +0.0033). 481: the key's floored share rises from 33 % to 39 % on efr.

**Verdict: robust on A/B/C, licensed for October; moderately fragile in its October size.** **Single settling check:** the
October-path pair already built, 460dr (LODO) vs 453dr (both `def_day.sh` on the relabelled A/B/C bundle, config diff
exactly `bulk_ctx_scale` + `bulk_load_sha256`). It prices the day-of fit (`scripts/def_lodo_key.py`) through the day's
re-level on the board. Read its mse, and expect ≈ 446bj2's +0.0022.

---

## 5. The re-level target ρ* 0.9793 and the NTC-proxy cliff (mandate decision 4; 447/448; DEF-RUNBOOK §5.-1)

`relevel487.py`, analytic. Slopes are ntclevel417's: raw mse −0.31 per unit scale below the cliff, +1.78 above. One raw
mse unit in one context = 1/(0.9502 × 18) overall. **Known answer PASS:** ρ*(sd 0.01) = 0.97934, a margin of −1.044 %
in scale.

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

- Random proxy error δ ~ N(bias, sd): shipping 0.9793 instead of the best margin costs ≤ 0.0001 per context for sd ≤ 1 %
  with no bias, ≤ 0.0003 at sd 2 %, and ≤ 0.0006 at sd 3 % with a C-type bias (−0.88 %).
- If the above-cliff slope on D/E/F were 2× or 4× A/B/C's, the optimal ρ* would be 0.972 / 0.966, and 0.9793 costs
  0.00003 / 0.00013 per context.
- The secant loop's ±1 % ρ tolerance (±0.5 % in scale) moves any row above by ≤ 0.0005 per context.

**Verdict: robust, with one structural trigger already found.** The loss is asymmetric but small at plausible errors
(~0.001 per context per 1 % past the cliff). The margin covers A/B/C's worst measured proxy error, and no slope or sd
misspecification tried here costs more than 0.0006 per context. **Single settling check:** on Oct 22, confirm cells per
`ntc_id` = cells_per_pert in each D/E/F control file (one line of pandas), or apply 481's F2 patch (A/B/C known answer
bit for bit).

---

## 6. Sharpening q 1.75 (mandate decision 3; 417/433/446)

Evidence: board steps 417a (q 1.25, +0.0070), 433k (q 1.5 levelled, +0.0046), 446b (q 1.75 levelled, +0.0023). All are
bulk-only (pds/mse reproduce to 16 s.f.), so the A/B/C curve is measured without ordering noise. 438 PLOLO on
October-like panels with 417a's calibration (k 0.65 mse, k_pds 0.42). Per-fold increments re-derived from
`438/work/summary438.json` (levelled, calibrated):

| panel | q 1.5 → 1.75, per fold (K562 / HCT116 / HEK293T) | mean | q 1.75 → 2, per fold | mean |
|---|---|---:|---|---:|
| abc | +0.0015 / +0.0022 / +0.0027 | +0.0021 | −0.0003 / +0.0005 / +0.0008 | +0.0003 |
| real | +0.0003 / +0.0010 / +0.0010 | +0.0007 | −0.0007 / −0.0001 / −0.0007 | −0.0005 |
| efp | +0.0013 / −0.0000 / +0.0008 | +0.0007 | +0.0002 / −0.0008 / +0.0001 | −0.0002 |
| efr | +0.0009 / +0.0013 / +0.0017 | +0.0013 | +0.0003 / +0.0004 / +0.0005 | +0.0004 |

- **Known answer:** the calibrated instrument retrodicts two board steps it was not fitted on. abc 1.25 → 1.5 +0.0044 vs
  board 433k +0.0046; 1.5 → 1.75 +0.0021 vs board 446b +0.0023.
- On October-like panels q 1.75 beats 1.5 on 8 of 9 fold-panels (the ninth is −0.00002), and q 2 is a wash. The optimum
  is 1.75–2 everywhere, and the regret of 1.75 is ≤ 0.0004 on any panel measured.
- Not in the PLOLO price: the ballast's nmae cost, which grows with q and with the re-level factor (446b: +0.0013 raw
  nmae inside its +0.0023; October factors 1.21–1.29 vs A/B/C's 1.10–1.17). The A/B/C board number includes it; it is small.

**Verdict: robust.** No check needed.

---

## Files
- Scripts: `sameline487.py`, `sameline487b.py` (decision 1 on H1); `lift487.py`, `lift487b.py`, `whitemix487.py`,
  `pt487.py` (in-mix lift and `pt`); `inmix487.py` (the emitter check, written for after 04:30Z); `prog487.py`
  (programme/gate by strength class); `relevel487.py` (cliff).
- Patch (not applied): `patches/in_mix_sameline.diff` + `patches/README.md`.
- Outputs in `work/`: `sameline487.json`, `sameline487b.json`, `lift487.json`, `lift487b.json`, `whitemix487.json`,
  `pt487.json`, `prog487.json`, `relevel487.json`; `RESULTS.interim-2230Z.md` (the interim record).
- Memory: every process ≤ 4 GB peak except `lift487.py`, 5.58 GB for ~8 s (mmap pages of full-width donor rows plus a
  float64 copy). free+inactive read 29.2 GiB just before it, under the 33 GiB gate. It was a deviation
  from the mandate's rule. Later scripts read rows one at a time (`lift487b.py` 0.33 GB, `pt487.py` 0.39 GB,
  `whitemix487.py` 3.2–3.8 GB, run while a coordinator build was past its start gate and free+inactive read 36.7 GiB).
- Writes: only this directory (checked with `find -newer` against the resume time; the `scripts/__pycache__` files
  written tonight predate this session's imports and come from other sessions' runs).

---

## Coordinator's disposition (2026-09-27 07:55Z; committed from the peer's working tree at its request)
- **§1, same-line rows: ACCEPTED, October default changed.** The day runs `def_day.sh` args 5–6 = **`0.5`** and
  **`noprog`**, or **`0.5`** and **`""`** when `DEF_DAY_PROGRAMME_SOURCE` is on. That means w 0.5 for every matched bank
  (every October match is a cancer line) and no per-target signs. Plan October as if the mode is worth ≤ +0.005 overall.
  The finer per-bank patch (`patches/in_mix_sameline.diff`) is not applied; the numeric w does the same job without a code
  change.
  - Rehearsal: experiment 510 (queued after 497e): efr with the forced call (D=K562, E=HCT116, F=A549) at `0.5` / `""`,
    together with every other candidate switch.
- **§2, prior gate: ACCEPTED as an amendment to 468's October table, written before any 09-30 entry exists.** With any
  programme source, the gate is `le1`, and a NULL 458g read (|Δ| < 0.0005) no longer drops it. Only a HARM read (≤ −0.0005)
  closes it: a board harm on A/B/C outweighs the H1 panels.
- §3–§6 (stack mix, LODO, re-level ρ*, q 1.75): recorded as robust. 481 F2's ntclevel patch is already in (474, a4b6b4cc).
  The dayof438 q-grid extension is noted, not queued.
