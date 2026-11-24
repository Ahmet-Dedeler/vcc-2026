# 503 — Manuscript claims audit

**70 of 73 claim rows are SUPPORTED. Three are OVERSTATED (C-01, C-59, C-61). No claim row is WRONG, STALE, or UNSOURCED.** Six further quantitative sentences in `PAPER.md` are wrong or stale even where the claim row is fine. The five corrections from 2026-09-26 (bit-reproducibility, kolf2 coverage, same-line r ≈ −0.02, control-noise fraction, gate benefit) are already in the manuscript or were never in it; they are not still live errors.

Re-run: `uv run python experiments/503-cloud-claims-audit/audit_claims.py` (65 checks, 0 missing; log `work/audit_log.txt`). Evidence below is CONFIRMED unless marked SUSPECTED. The paper was not edited.

## Non-SUPPORTED claims

| claim | verdict | what is wrong | replacement |
|---|---|---|---|
| **C-01** | OVERSTATED | +0.243, +0.054 and +0.0199 are three real numbers on three different surfaces (H1, 85 targets; CD4, 297 targets; the cross-line map, 300 targets). They are not "the same panel" and not "the same cells". The map median is **+0.0199**, which the paper rounds to **≈ +0.03**. | "Three different measurements, on different cells: an H1 Flex split-half median residual cosine of +0.243 (85 targets, 200 cells/half), across-donor agreement in primary CD4 T cells of +0.054 (297 targets), and a cross-line raw-cosine median of +0.020 (300 panel targets). They are not one quantity." |
| **C-59** | OVERSTATED | The H1 mapping 0.950 / 0.882 / 0.844 / 0.788 at a = 0.56 / 0.62 / 0.65 / 0.70 matches 434f. "Entries at 0.82–0.84 therefore imply gate sign accuracy near 0.65" does not. FINDINGS §1 (2026-09-26) supersedes the implication: the map assumes iid signs and constant magnitudes, so those accuracies are **upper bounds**. | "Under iid sign flips and L1-optimal magnitudes, H1 nmae is 0.950 / 0.882 / 0.844 / 0.788 at accuracy 0.56 / 0.62 / 0.65 / 0.70. Board nmae of 0.82–0.84 sits on the 0.65 point of that curve. If magnitudes track confidence, the same nmae is reached at a lower average accuracy, so 0.65 is an upper bound on what those entries need, not a measurement of their signs." |
| **C-61** | OVERSTATED | mse −0.0789 vs −0.0797 and nmae +0.0544 → −0.0006 match the board. "All four DE members within 0.0006" does not. nmae is 0.936458 − 0.937101 = **−0.000643**. The board write-up says **within 0.0007**. | "Ballast keeps the mse gain (−0.0789 vs −0.0797 raw) and removes the weld's nmae cost (+0.0544 → −0.0006). The other three DE members move by at most 0.00009. All four sit within 0.0007 of the reference; the largest move is nmae −0.000643." |

## Paper sentences that the claim table does not carry

| where | verdict | source | replacement |
|---|---|---|---|
| Round 10 header, PAPER L36: a quarter-weight curve "keeps **63 %** of the gain available" | WRONG | `450/BOARD-RESULT-446bja.md` L57: the consensus "holds only **37 %** of the w 0 → 1 gain". The body (L510) and C-63 already say 37 %. 63 % is the complement, stated as if it were the retained fraction. | "Dilution does not explain the null: on H1 a quarter-weight consensus keeps only 37 % of the same-line programme's gain, and the board moved fid by +0.0002." |
| §3(j) L493: six-member score **+0.071 … +0.084**, "mse −0.17, nmae −0.07, reach +0.05 … +0.07" | OVERSTATED | `450/RESULTS.md` L183–185: rows at w 0.5 / 1 / 2 are **+0.0531 / +0.0710 / +0.0790**. No row is +0.084. w 1 + per-target prior is +0.0837 (L191); the adopted arm is **+0.0934 ± 0.0002** (L193, 3 seeds). The member figures mix w 1 and w 2. | "Rows at weight 0.5 / 1 / 2 raise the six-member score by +0.053 / +0.071 / +0.079 (mse −0.13 / −0.17 / −0.19 raw). With a per-target sign prior at weight 2 the score is +0.093 (3 seeds, seed sd 0.0002): mse −0.19, nmae −0.081, reach +0.076, pds +0.010." |
| §3(a) L127: "the measured value is **0.330** at 200/200 cells" | OVERSTATED | `273/work/tables273.md` L39, same row: median **+0.2429**, mean **+0.3309**. The sentence never says mean, next to a median of +0.243. | "The retracted 0.481 at 400 cells was a linear extrapolation. At 200/200 cells the measured split-half residual cosine has median +0.243 and mean +0.331." |
| §4.6 L776: "nothing here is measured on the challenge's held contexts A/B/C" | STALE | §3(j) and §3(k) are scored board entries on A/B/C (446bja, 410/410b, 433a). The sentence was true of the 6 September draft. | "The reliability, conservation and construct results are not measured on A/B/C. The programme test (§3(j)) and the two-functional decomposition (§3(k)) are: they are scored board entries on those three contexts." |
| Figure index L815: F-10 is "**8 banks, 3 labs**" | STALE | Body L517–519 and `450/RESULTS.md` L355–378: fifteen banks, four labs; the ≥1,000-perturbation table is the controlled subset of that. `figs/MANIFEST.md` F-10 no longer says eight banks. | "F-10: programme agreement on matched perturbation sets (pairs sharing ≥1,000 perturbations; fifteen banks, four labs)." |
| §3(j) L620: the leader "trains on a **100M-cell**" proprietary atlas | STALE | `docs/FINDINGS.md` L349–351: Illumina's Billion Cell Atlas, **~150M cells by January 2026**, 1B targeted. Tahoe-100M is a different, public, chemical atlas. The sentence has no experiment citation. | "A leader appears to train on a proprietary many-line perturbation atlas (public descriptions: on the order of 10^8 cells across 200+ lines, with a larger target). That is the resource that would supply the oracle row above." |
| §3(j) L569 heading: "A design-strength gate **removes** that cost" | OVERSTATED | The next sentences, and C-71, already give the interval that spans 0. The heading does not. | "A design-strength gate removes the point-estimate reach cost; its own benefit is not established (target-bootstrap interval spans 0)." |

## Already corrected, and not still wrong in the manuscript

| 2026-09-26 correction | manuscript status | evidence |
|---|---|---|
| Day pipeline "not bit-reproducible" | Not a paper claim. 452's own RESULTS L79 corrects it: the hash differed because the dump path is inside the file. C-58 and C-62's "bit-identical" pds/mse are true of the board JSON (16 digits). | CONFIRMED |
| kolf2 October coverage "0" | Not in PAPER or CLAIMS. C-19's "kolf2 273" is a PLOLO fold size. FINDINGS L68 has 181–251/300. | CONFIRMED |
| Same line across labs r ≈ −0.02 | Withdrawn in the source (`450/RESULTS.md` L361 is v2; L391 withdraws it) and in the paper (L534–535, "superseded") and in C-64, which quotes the matched-set v3 numbers. | CONFIRMED |
| Control-noise fraction | C-66 and PAPER L46/L560–562 state the corrected decomposition: 56 % signal, 33 % control-sample noise, 12 % target noise (`462/RESULTS.md` L63), not the first reading that a large control pool makes the programme almost clean. | CONFIRMED |
| Gate benefit quoted with seed sd only | C-71 and PAPER L573–575 include the target bootstrap +0.0043 [−0.0034, +0.0132]. The seed sd 0.0002 is labelled as emission noise. The heading above still overstates; the numbers do not. | CONFIRMED |

## Full table

Unit of a board member is raw `cell-eval2` unless the row says scaled or overall. "Findings L" is `docs/FINDINGS-2026-09.md`. Claims whose headline number sits in that note were re-read there; experiment 275 had already checked that note against the JSON. This audit did not recompute the bootstraps.

| claim | verdict | source line | evidence |
|---|---|---|---|
| C-01 | OVERSTATED | `276/RESULTS.md` L102 median C_raw **+0.0199** (300 targets); `273/work/tables273.md` L39 median **+0.2429** (85 targets); Findings L303 CD4 **+0.054** (297 targets). Paper L61 and L87 say same cells / same 300-target panel. | CONFIRMED |
| C-02 | SUPPORTED | Findings L47–49: +0.047 / +0.028 / +0.028, 5.6 / 32 / 34 %, n 269 / 290 / 294, cells/half 81 / 70 / 95. | CONFIRMED |
| C-03 | SUPPORTED | `273/work/tables273.md` L39: median +0.2429, p10 +0.0673, p90 +0.7479, share ≤ 0 = 0, n 85. The unlabeled 0.330 in the paper is the mean on this row (see above). | CONFIRMED |
| C-04 | SUPPORTED | Findings L49–50: batch +0.028, random cells +0.071, n 294. | CONFIRMED |
| C-05 | SUPPORTED | Findings L59–60: Spearman with pds +0.58 / +0.23 / +0.17 / +0.62; with cells +0.16 / +0.09 / +0.12 / −0.05. | CONFIRMED |
| C-06 | SUPPORTED | Findings L425: +0.0078 [−0.0068, +0.0119] on h1-test; +0.0012 / −0.0017 on K562. The interval covers zero, and the claim shows it. | CONFIRMED |
| C-07 | SUPPORTED | Findings depth-law paragraph: measured − Spearman–Brown from n = 25 is −0.006 [−0.030, +0.033] on h1-test. Paper L133 quotes predicted +0.2455 vs measured +0.2392. | CONFIRMED |
| C-08 | SUPPORTED | `285/RESULTS.md` L214: 748 perturbations, guide cosine +0.0006 [−0.0002, +0.0014], pds 0.509. Cell-half +0.0063 and the +0.0058 difference are the paired contrast in the same section. | CONFIRMED |
| C-09 | SUPPORTED | `285/RESULTS.md` L252: Cas9 guide ρ 0.81 / 0.76. The ladder (batch 0.87–1.12, study 0.79, line 0.35–0.45, donor 0.02–0.07) is the 285 ρ table the claim cites. Ratio rows without stored intervals are labelled as such in the claim. | CONFIRMED |
| C-10 | SUPPORTED | `285/RESULTS.md` L262: n* = 89 (Jurkat), 76 (HepG2), 328 (h1-test). | CONFIRMED |
| C-11 | SUPPORTED | Paper L151 terciles +0.044 / +0.140 / +0.506 and panel +0.0095 vs essential +0.0605 match the 285 strata the claim cites (`285/RESULTS.md` L217 has essential +0.0605 at n = 25). | CONFIRMED |
| C-12 | SUPPORTED | Findings L303: median +0.054 [+0.043, +0.062], 28.7 % ≤ 0. | CONFIRMED |
| C-13 | SUPPORTED | Findings L305: Δ +0.0033 [+0.0004, +0.0072] against cell-sampling +0.099 [+0.092, +0.106]. Three-donor +0.0067 is the matched-budget contrast in the same finding. The claim says the donor interval does not clear zero in Rest / Stim8hr. | CONFIRMED |
| C-14 | SUPPORTED | Retraction. Findings L40–41 withdraws the depth ceiling. The paper does not quote it as a result. | CONFIRMED |
| C-15 | SUPPORTED | Findings L251: ρ +0.273 [+0.158, +0.388]. | CONFIRMED |
| C-16 | SUPPORTED | Findings L251: ρ +0.435 [+0.340, +0.522]; 50/50 random splits positive. | CONFIRMED |
| C-17 | SUPPORTED | Findings L253: PLOLO ρ +0.387 [+0.326, +0.444]; Flex +0.239 is in the same sentence. | CONFIRMED |
| C-18 | SUPPORTED | Off-map +0.274 and +0.210, with the Broad-16 interval touching zero, are the 291 surfaces the claim cites. The claim records the touch-zero interval. | CONFIRMED (291 surfaces as carried into the claim; not re-bootstrapped) |
| C-19 | SUPPORTED | Findings L258: partial +0.26 / +0.25 / +0.23 vs raw +0.41 / +0.39 / +0.37; weakest tertile +0.03 [−0.18, +0.24]. "About 40 %" is that note's summary of the drop. kolf2 273 is the fold size, not an October coverage. | CONFIRMED |
| C-20 | SUPPORTED | The claim is that C1 is a robustness check. Ledger relabel is what the row reports; half-B ρ +0.361 [+0.257, +0.458] is the passed-then-relabelled number. | CONFIRMED |
| C-21 | SUPPORTED | Findings L278: Δmedian −0.01047 [−0.01269, −0.00844], AUC 0.564 [0.548, 0.581]. | CONFIRMED |
| C-22 | SUPPORTED | Family rank and the core-machinery to identity axis are the 279/284 results the claim cites. No single headline number was re-extracted beyond the axis the paper quotes from `fam284.json`. | CONFIRMED (scope matches the claim's "not ZNF-specific") |
| C-23 | SUPPORTED | Findings / 279: AUC 0.517 → 0.551 → 0.614 and 74 % of the deficit surviving reliability correction (`279/RESULTS.md`: AUC 0.561 → 0.567; "74 % of the deficit survives"). The claim states the terciles are not expression-matched. | CONFIRMED |
| C-24 | SUPPORTED | Findings L360 carries ρ +0.805 [+0.679, +0.883] over 78 families under the strict fold. | CONFIRMED |
| C-25 | SUPPORTED | Features − reliability −0.035 [−0.082, +0.033] is the 284 gene-level result the claim cites. Interval covers zero, and the claim says so. | CONFIRMED |
| C-26 | SUPPORTED | Findings L360 carries remainder agreement ρ +0.328 [+0.241, +0.373], ceiling +0.382, permutation null +0.044 ± 0.018. | CONFIRMED |
| C-27 | SUPPORTED | Findings L381: +0.317 [+0.290, +0.391] vs annotation −0.057 [−0.138, −0.021]. | CONFIRMED |
| C-28 | SUPPORTED | The claim is the sign inversion (proliferation +0.118, PC energy +0.115), reported as a refutation. Those are the 290 per-feature intervals. | CONFIRMED |
| C-29 | SUPPORTED | Findings narrows the ablation: hallmark-only +0.158 vs random +0.099 ± 0.019; spectral-only +0.242; full +0.317. The claim includes that narrowing. "Costs only +0.023" is the in-map ablation and is labelled as such. | CONFIRMED |
| C-30 | SUPPORTED | gen_pc20 +0.295 and +0.186, and the geometry-minus-annotation intervals that cover zero, are the failed C3/C4 conditions. The claim says the headline does not travel. | CONFIRMED |
| C-31 | SUPPORTED | Explicitly not claimed. The DNA-binder exception is dropped. | CONFIRMED |
| C-32 | SUPPORTED | Findings L177: shared 0.091 (source JSON 0.090681, which the claim also quotes). Line-specific 0.115 and noise 0.794 match the claim's three-way split. | CONFIRMED |
| C-33 | SUPPORTED | Predicted-residual cosine +0.0110 not separated from the wrong-recipient control, and the full response losing, are the 289 LOLO result. The claim's reconciled wording ("no detected improvement for the two frozen models") matches the confirmation status. | CONFIRMED |
| C-34 | SUPPORTED | HPV-E6 ordering, permutation P = 0.018, and gating +0.053 [+0.015, +0.089] are the 289 p53 result the claim cites. Mechanism is credited to the source paper. | CONFIRMED |
| C-35 | SUPPORTED | C2a fail (1 of 4), C2b pass, detectable component +0.0347 [+0.0229, +0.0468]. The claim does not say the size reproduced. | CONFIRMED |
| C-36 | SUPPORTED | Preregistered ρ +0.262 (wrong sign) and deconfounded −0.144 [−0.403, +0.126] are the 282 result, reported as not confirmed. | CONFIRMED |
| C-37 | SUPPORTED | Explicitly not claimed. Absolute difference +0.0080 [−0.0026, +0.0186] is why the ratio was dropped. | CONFIRMED |
| C-38 | SUPPORTED | The crossed point estimates (+0.0424 / +0.0116 / +0.0141) are the 282 small-n design. The claim does not attach an interval the source lacked. | CONFIRMED |
| C-39 | SUPPORTED | 76.9–78.7 % of directions and 87.5–91.5 % of energy are the 272 spectra. In-sample anti-correlation and the held-line read are in the same analysis. | CONFIRMED |
| C-40 | SUPPORTED | Findings L122: Wiener −0.034 [−0.044, −0.024]. Hard cut −0.069 and soft floor −0.0022 are the same table. Flex intervals cover zero and the claim shows them. | CONFIRMED |
| C-41 | SUPPORTED | 78.1–79.8 % of directions below the independent noise level, 0/12 cells clearing, best −0.0080 [−0.0163, +0.0001]. `286/RESULTS.md` §6 is the cited source; Findings records the same 0/12. | CONFIRMED |
| C-42 | SUPPORTED | The claim is the caveat (the knee moves with the decoy count: slope 9.5 vs 6.0), not a discovery. | CONFIRMED |
| C-43 | SUPPORTED | 17.5 % anti-aligned, 66 % marginal, ρ +0.875, cross-line Spearman +0.76 / +0.70 / +0.74. Findings L162 has the Spearman. Shares are the 262 decomposition. | CONFIRMED |
| C-44 | SUPPORTED | Findings L81: cross-validated −0.018 [−0.027, −0.009]. In-sample +0.056 and retained fraction −0.32 are the same table. | CONFIRMED |
| C-45 | SUPPORTED | Median 5 significant genes, 85/272 at a 10-gene bar. Findings credits Replogle's 30.5 % on the same data. | CONFIRMED |
| C-46 | SUPPORTED | Findings L460: **25 of 300** targets shared with the 2025 H1 library. Ceiling cosines +0.213 / +0.572 / +0.690 are the 293 measurement the claim cites. | CONFIRMED |
| C-47 | SUPPORTED | `295/RESULTS.md` L3 and L66: ρ +0.5016, file-block [+0.1525, +0.5659], gene [+0.3709, +0.6083]. The claim also carries the 298 jackknife (C5b marginal without Adamson). | CONFIRMED |
| C-48 | SUPPORTED | Explicitly not a scientific claim. Calibration constants stay out of the results sections. | CONFIRMED |
| C-49 | SUPPORTED | Findings L460 and paper L155: ρ_guide −0.003 [−0.085, +0.181] (11) vs +0.687 [+0.197, +1.067] (5); 38.2 % vs 16.6 %, p = 0.0015. | CONFIRMED |
| C-50 | SUPPORTED | The row is a rate on the 2025 panel (18.7 % below \|lfc\| 1.0) and the addendum refuses the 2026 / off-target extrapolation. That refusal matches the paper's Round 5 header. | CONFIRMED |
| C-51 | SUPPORTED | Jurkat +0.458 vs +0.096, K562 +0.340 vs +0.001, HepG2 not separated, Broad-16 no split, CD4 flat. Those are the 296 strata. | CONFIRMED |
| C-52 | SUPPORTED | `296/RESULTS.md` L153: null-corrected K562 ≈ HepG2 ≈ CD4 ≈ 0.03–0.06, Jurkat 0.08, Broad-16 0.79. (b) − (x) through zero except CD4, as the claim says. | CONFIRMED |
| C-53 | SUPPORTED | Failing-construct split-half above zero in all five libraries, and the 301/307/309 failing-vs-working intervals that cover zero and flip sign, are what the claim reports. It does not turn that into a shared vector signature. | CONFIRMED |
| C-54 | SUPPORTED | Study cosine +0.168 [+0.120, +0.215] vs construct +0.011 [+0.004, +0.018] on 112 perturbations is the 298 raw-unit inversion. | CONFIRMED |
| C-55 | SUPPORTED | Findings L574: S +0.786 [+0.387, +1.166]. The claim states the gene-sampling null, the known-answer abort, and the KNOWN/EXTENDED audit. | CONFIRMED |
| C-56 | SUPPORTED | Findings L574: S +0.71293 [+0.32419, +1.07351]. The claim states that C9 is not an independent replication and that the headline p resamples genes. | CONFIRMED |
| C-57 | SUPPORTED | The row lists what is not claimed. The paper's §4.4 matches that list, including the withdrawn wet-lab recommendation. | CONFIRMED |
| C-58 | SUPPORTED | `426/RESULTS.md` L1 and L38: nmae −0.00962 ± 0.00016 (5 seeds); coarse −0.00953 ± 0.00021 (`426/RESULTS.md` L114). Board: `433/BOARD-RESULT-433a.md` L11–17 and the two entry JSON files, pds 0.8589520624303233 and mse 0.6798392234982936 identical, Δnmae −0.016786. The claim's ±0.0002 matches the coarse seed sd. | CONFIRMED |
| C-59 | OVERSTATED | `434/RESULTS.md` L95–101 has the mapping and states the implication. `docs/FINDINGS.md` L257–258 supersedes the implication with the upper-bound caveat. The paper L425 still says "they run". | CONFIRMED |
| C-60 | SUPPORTED | `434/RESULTS.md` L40–42: panel-like full features −0.0066 vs plain MDE −0.0111. L80–83: uncalled share of \|l\| 0.709, sign accuracy 0.537. | CONFIRMED |
| C-61 | OVERSTATED | `410/BOARD-RESULT-410.md` L8–9: mse −0.079684, nmae +0.054375. `410/BOARD-RESULT-410b.md` L8–10 and L16: mse −0.078904, nmae −0.000643, "within 0.0007". | CONFIRMED |
| C-62 | SUPPORTED | Board JSON: 433a and 417e pds and mse identical at 16 digits (see C-58). fid, reach and jac are not identical; the claim does not say they are. | CONFIRMED |
| C-63 | SUPPORTED | `450/BOARD-RESULT-446bja.md` L11–16: fid +0.000219, overall −0.0010661, pds and mse 0. L57: 37 % of the w 0 → 1 gain. The Round 10 header's 63 % is a separate paper error (above); this row is the 37 % figure. | CONFIRMED |
| C-64 | SUPPORTED | `450/RESULTS.md` L375–378: +0.806, +0.205, +0.544, +0.175 on pairs sharing ≥1,000 perturbations. L384–392: line component "about +0.2"; v2's −0.020 withdrawn. The claim discloses point estimates and no interval; the source has none. "Lab adds ~0.37" is 0.544 − 0.175. Paper L39 "proves" is correlational wording; the numbers are the correlations. | CONFIRMED |
| C-65 | SUPPORTED | `450/RESULTS.md` L317–326: same-line cross-lab cosine 0.057 (mean of 0.066 and 0.048), other-line mean 0.033; within Weissman 0.216 vs 0.129. Both ratios are about 1.7. Common set 1,129 targets, 6,083 genes (L312). Point estimates, as the claim says. | CONFIRMED |
| C-66 | SUPPORTED | `462/RESULTS.md` L32–33: Broad 0.801 vs 0.121. L59–61: H1 harness 0.793 vs 0.386. L68: full cells shared 0.919, disjoint 0.687. L63: 33 % control-sample noise. One seeded split, as the claim says. | CONFIRMED |
| C-67 | SUPPORTED | `462/RESULTS.md` L34–38: cross-line 0.043 vs within-line 0.121, "about a third". **132 ordered pairs is 12×11, not a stored count** (SUSPECTED as a label; the r values are CONFIRMED). | CONFIRMED on the correlations; SUSPECTED on "132" |
| C-68 | SUPPORTED | `462/RESULTS.md` L88: shared 0.6184 vs disjoint 0.6183. L76–77: fold changes match the scorer to float32 accumulation (median 1.6e-6, max 1.2e-5), which is the precision the claim names. | CONFIRMED |
| C-69 | SUPPORTED | `462/RESULTS.md` L152–154: cross-line 0.453, ratio 0.53. L160–168: different-lab different-line 0.376 (9 pairs) vs same-lab different-line 0.431 (3 pairs). | CONFIRMED |
| C-70 | SUPPORTED | `466/RESULTS.md` L19–28: non-essential same-lab 0.182 / 0.165 and cross-lab 0.059 / 0.094 and 0.050 / 0.076, summarised there as 0.17–0.18 vs 0.05–0.09; essential cross-lab 0.41–0.52 vs same-lab 0.35–0.43. One seeded split, disclosed. | CONFIRMED |
| C-71 | SUPPORTED | `467/RESULTS.md` L7–11 and L34: fid +0.0161, nmae −0.0102, reach +0.0015, overall-equiv +0.0110 [+0.0067, +0.0158]; gate benefit +0.0043 [−0.0034, +0.0132], P(≤0) about 0.15. The claim uses the target bootstrap, not the seed sd, as the uncertainty. The section heading still says the gate "removes" the cost (paper sentence above). | CONFIRMED |
| C-72 | SUPPORTED | `470/RESULTS.md` L57 and `470/PREREGISTRATION.md` L40–45: 0.652 (5,192), 0.633 (508), 0.666 (793), 0.580 (8,728), 0.624 (8,102), 0.605 (147). Weak-third down fraction 0.631 and strong 0.554 are `467/RESULTS.md` L66–67. Point estimates; the source stores no interval. Descriptive, not a causal rate. | CONFIRMED |
| C-73 | SUPPORTED | `496/RESULTS.md` L20: nmae −0.00146, pds and mse exactly 0, pooled overall-equiv +0.00018 [−0.0017, +0.0020], P(≤0) 0.41. Seed sd 0.00002 is the nmae seed scatter (L17–19: −0.00147 / −0.00143 / −0.00147). The claim includes the overall interval. The §3(i) paragraph quotes only the seed sd; it should also quote [−0.0017, +0.0020]. | CONFIRMED |

## Causal wording that the numbers do not support

Kept only where the sentence asserts a mechanism the design did not isolate. Correlations that are labelled as correlations are not listed.

- C-01 / abstract: "we separate them and find them far apart on the same cells" treats three studies as one experiment.
- C-59 / §3(i): "so they run gate sign accuracy near 0.65" turns a calibration curve into a measurement of rivals' signs.
- Round 10 header: "the programme proves to be mostly a lab and protocol signature" is a set of pairwise correlations (C-64).
- §3(j) heading: "a design-strength gate removes that cost" is a point estimate whose target interval includes zero (C-71).

## What this audit did not do

It did not recompute bootstrap intervals from the saved matrices, and it did not re-score any h5ad. Where a claim's interval is the one printed in `docs/FINDINGS-2026-09.md` or in the experiment's RESULTS, the check is that the manuscript copied that printed number, including its unit and its n. A well-supported negative: the five known 26 September corrections are not still sitting in the claim table as if they were current.
