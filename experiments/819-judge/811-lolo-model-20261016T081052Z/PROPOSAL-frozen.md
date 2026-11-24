# 811 PROPOSAL — M811: a many-line consensus, de-shared and reliability-weighted, chosen by nested folds

Written 2026-10-16 ~08:20Z by the 811 worker, after the preregistered baseline (RESULTS.md §2) and before any M811 number exists.
For the judge (PHASEB "What you owe"). Nothing below has been priced.

## 1. What the baseline showed (the reason for this proposal)
- **More public lines carry more transferable direction, on fresh folds.** On BULK-W the all-lines consensus beats the
  sources-only consensus at the same λ on 17/18 folds: +0.0018 / +0.0050 / +0.0105 overall at λ 0.10 / 0.25 / 0.50, with
  literal pds also positive (+0.0012 / +0.0040 / +0.0076). Preregistered prediction 2 (NULL) is falsified. On A/B/C the
  same contrast leans positive (pdsl +0.0011 / +0.0031 / +0.0058, 3/3 folds at 0.25 and 0.50), but the literal-pds intervals
  span 0. On ABC/kolf2 it is strongly positive (+0.0036 … +0.0224).
- **Against the incumbent, the baseline is a W-only effect.** W gains +0.0039 … +0.0093 (18/18 folds). A/B/C is NULL at
  λ 0.10 (+0.0005 [−0.0007, +0.0017]) and HARM at 0.50 (−0.0074). The pattern matches 706: mse improves everywhere
  (r +2 … +7 % on A/B/C) and literal pds pays (−0.003 … −0.034).
- **Mechanism of the pds cost (align811, raw delta space, no pricing).** The all-lines consensus carries 2–4× more
  across-target shared energy than the truth: 0.021–0.046 against the truth's 0.005–0.017. Sources-only carries 0.008–0.015.
  So the extra lines bring panel-common structure, which costs target discrimination (the 343–382 "coherence is the price"
  mechanism). One line is a specific culprit. KOLF2 is a z-score bank with a 0.245 shared share; dropping it from U_all
  raises U's raw per-target cosine with every A/B/C cancer truth (+0.009 / +0.009 / +0.015).

## 2. The model (transfer class; every input computable on Oct 22 from public lines, D/E/F's line never in any bank)
For a panel's target t and the non-held line groups g (RESULTS §1; on D/E/F every public line is non-held):
1. **Per-line centring (factor C).** Each line's unit rows L_g(t) have that line's panel mean (mean over the panel's
   targets that the line carries) subtracted, then are re-normalised. This removes lab blobs, z-score offsets and generic
   responses before averaging. Mechanism: the shared-share excess in §1. Computable day-of over D/E/F's own panel.
2. **Reliability weights (factor W).** w_g ∝ max(a_g, 0) / (1 − a_g), where a_g is the median over the panel's targets of
   cos(L_g(t), mean of the OTHER non-held lines' rows) (inner leave-one-line-out among the non-held lines only). Lines
   that agree with the rest weigh more. The held line is never touched; on D/E/F this uses only public lines.
3. **Blend.** U(t) = normalise(Σ_g w_g L_g(t)); the composition row is blended norm-preserving at dose λ (the baseline's
   formula), then the incumbent's key, sharpening, render and ballast run unchanged.

Grid: C ∈ {off, on} × W ∈ {equal, reliability} × λ ∈ {0.10, 0.25}, so **8 configurations**.

## 3. Protocol (nested, so the price is not a best-of-N read)
- **Select** on W-tune = L1–L3 × 3 held lines (9 folds): the configuration with the highest best-of-8-shifted lower bound
  (literal-pds overall vs the incumbent).
- **Confirm** that ONE configuration, once, on W-conf = L4–L6 (9 folds) and on BULK-ABC (3 cancer folds; kolf2 apart),
  plus the paired contrast against the baseline at the same λ.
- λ 0.25 configurations move ballast K by up to +5.9 % (baseline meta), so the confirmed winner goes through DE-H1
  (803's de803.py) before any dossier if its |dK_rel| exceeds 5 %.
- The sealed LZ is the lab lead's to open, once, for a final candidate only.

## 4. Predictions and bins (before any number)
- **P1 (mechanism).** C on cuts the A/B/C literal-pds cost at λ 0.25 from −0.0127 to ≥ −0.005 and keeps ≥ 60 % of the r gain
  (+5.1 %).
- **P2.** W (reliability) ≥ equal by ≤ +0.0005 overall. Small, because the cancer partners already dominate the agreement.
- **P3 (the bar).** The confirmed winner on BULK-ABC: overall ≥ +0.0010 (literal pds) with the 2.5 % bound > 0, and W-conf
  within ±30 % of W-tune. Bins as DESIGN §6/§8: SHIP-like / NULL / HARM.
- **Falsifier.** If centring removes ≥ 70 % of the r gain, then the many-line gain is mostly a shared (generic-response)
  component, not target-specific direction. That would be scorer mechanics on mse, not transfer. I would report it as such
  and close the many-line family for October. Separately, if the winner is NULL or HARM on BULK-ABC, the family stays a
  W-panel effect, like 706, and is closed for October.

## 5. Novelty against FINDINGS §2 (what it is NOT)
- Not a fifth donor in the whitened stack or a residual at a dose (218/225/229/248 "sources saturate"). The consensus is a
  raw direction blended into the composition rows, and the baseline measured +0.005 (λ 0.25) over sources-only on 17/18
  fresh folds, which those closures do not cover.
- Not basal-similarity weighting, inverse-variance pooling or MODZ (Aug closures). Those re-weighted the three production
  donors inside the whitened stack. W here weighs 6–15 lines by their response agreement among the non-held lines.
- Not 341's cross-line ridge (a learned linear map; own-cos +0.003, retrieval −0.013 … −0.048), and not SVD truncation. No
  map is learned. Per-line centring is the step aimed at the "coherence is the price" closure (343–382).
- Not 706/707 (sources only, top z_src tertile: 810's lane). M811 is all-target and many-line. The two can interact, since
  both blend toward a donor mean, so the lab lead should not stack them without a joint price.
- The mandate's low-rank target × gene model with line-specific gene scalings is deferred, for a measured reason. The lever
  the baseline found is variance and shared-component control across lines, which M811 captures with two interpretable
  switches. A low-rank add-on would be priced afterwards, under the same nested protocol, as one more try.

## 6. Cost, tries, day-of path, dependency
- EC2: one base811-style run with 8 arms on the 9 W-tune folds plus aggregation, then 1 arm on 9 + 4 folds. About $0.5 in
  all (the baseline cost $0.45 including two failed diagnostic launches). DE-H1 if needed: about $0.4 on 803's engine.
- Tries: 8 (W-tune), 1 (W-conf), 1 (ABC). The baseline's 6 configurations stay on record in the same family.
- Day-of: D/E/F's 300 targets give rows from the same banks (CD4 covers ~283 by Arc's coverage rule). Centring and weights
  are computed over D/E/F's own panel, and the blend goes in before key and sharpening. Nothing is fitted to A/B/C.
- **Data (813, not this lane): two near-full-coverage lines are missing on fresh panels.** These are TeloHAEC genome-wide
  (298–300 per panel) and the KOLF2 full rows (244–252; z-scores, so per-line centring applies). Today, on L folds, CD4 is
  the only full non-cancer line. M811 would be re-priced with them as one more try once 813 wires them.
