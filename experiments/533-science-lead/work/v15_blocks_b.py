"""533: Draft v15 text blocks, part B (§3(i), §3(j)). See v15_blocks_a.py for the conventions."""

S3I_ADD = """

*A second delivery confirmed on the board (446bj2z, 2026-09-28; C-74).* An exact-rank pass moves counts between carrier
cells of down-regulated "mass-holder" genes, inside the controls' intervals, so that their arithmetic means change while
every per-gene total, every library, every ballast row and every Wilcoxon rank sum stays fixed.
- On H1 it lowered nmae by −0.00183 ± 0.00007 on top of the coarse MDE (3 seeds), with the other five members bit-identical.
- On the board it moved nmae −0.002179 raw, inside its preregistered −0.0013 … −0.0029. pds, mse, fid, reach and jac were
  bit-identical to the last digit, so the call set, the signs and the ranks were held exactly as designed.
- Transfer: 1.19× the H1 price (433a's MDE: 1.76×). Pure magnitude changes to this member have transferred at 1.2–1.8× their
  H1 price, with the same sign.
- A delivery-defect change transferred below its range. The zero-gene repair at a floor of −2 instead of −5 read −0.0014
  raw against a preregistered −0.002 … −0.008 (446br2). Such a change scales with the defect's prevalence in the scored
  contexts, not with H1's.

*Called genes already sit at H1's L1 optimum (513f; C-75).* Scaling the delivered fold changes of called genes by s, with
calls, signs and order fixed, raises H1's nmae by +0.0157 at s 0.5, +0.0103 at s 1.5 and +0.0356 at s 2. So s = 1 is the
optimum (the curve bottoms near 1.1). That is what the closed form predicts when 40 % of called signs are wrong: called gate
genes are 0.600 accurate (uncalled 0.537), and each unit of magnitude buys only (2a − 1). Magnitude changes have transferred
above their H1 price, which suggests the contexts' called-gene sign accuracy is higher than H1's and their optimum s above
1. That inference is untested.

*The add's direction from a same-lab programme (517; C-76, C-77; Figure F-13).* The 496 hook takes its sign from the
cross-line consensus. 517 took it instead from the benchmark producer's own 2025 H1 screen. The programme was estimated as the
mean z-score over that screen's non-essential perturbations, with h1-test's own targets excluded. It is applied where
|P_down − 0.5| ≥ 0.05, on design-strong targets; the weak third keeps the consensus. Through the production emitter on H1:
- nmae −0.01832 raw (seeds −0.01830 / −0.01826 / −0.01841), 12× the consensus hook's −0.00147;
- fid +0.0002 and jac −0.0001 (no call leak); pds and mse exact; reach −0.0039 (−0.0015 / −0.0026 / −0.0076 by seed);
- overall-equivalent +0.0044 (seed range +0.0037 … +0.0049; seed-0 target bootstrap [+0.0031, +0.0066]).

This is the programme's own line, so it is an upper reference. For the challenge's lines the question is how much of the
programme's information survives the change of line. Its bar is lower than a block prior's:
- On the genes the design does not call, the hook replaces the design's own delivered sign. That sign is right 0.532 of the
  time on H1's uncalled gate genes, where the programme is right 0.685 (P_down, τ 0.05; 0.702 for the mean-z estimate). A
  block prior must instead beat the cross-line consensus on the high-expression block.
- A weaker programme can be emulated by flipping a random fraction r of its signs. On nmae, the member it moves, the hook
  beats the consensus hook until r ≈ 0.42 (mean z) … 0.38 (P_down). Over all six members it stays positive to r ≈ 0.28 …
  0.37, against r ≈ 0.17 (every target) … 0.24 (design-strong targets) for the block prior. The six-member reading carries
  a small-gate reach term from r 0.2 that H1 cannot price for sign priors (Figure F-13).
- The estimator (mean z over non-essential rows rather than P_down) was chosen on these instruments and on the cross-line
  analogue of §3(j), so its H1 price is in-sample.
- The board test is preregistered: arm 446bj2zmz, primary Δnmae raw vs 446bj2z; ≤ −0.0030 transfers, −0.0030 … −0.0008
  partial, > −0.0008 fails; pds and mse must be bit-identical. It has not been read."""

F18_SOURCES_OLD = "same-line effect through a benchmark's own scorer. Sources: experiments 440, 441, 450.**"
F18_SOURCES_NEW = ("same-line effect through a benchmark's own scorer. Sources: experiments 440, 441, 450, 462, 466, 467, "
                   "516, 517. The same laboratory's programme from another line is tested on the challenge's own lines by "
                   "the 09-30 entries: [[PENDING 09-30]].**")

S3J_NEG_ADD = """

*Target-specific public sources do not supply the missing signs (516; C-79).* If the recipient's programme is out of reach,
can other screens' rows for the same target recover the uncalled genes' signs? No construction passed the bar (+0.02 in
pooled pool-sign accuracy on H1 and on a cancer fold, intervals excluding 0).
- The best construction is a Stouffer consensus of six same-target public rows (K562, two Orion lines, three CD4 conditions),
  noise-standardised per bank and applied only to genes the design does not call. It raises H1's pooled pool-sign accuracy
  by +0.0157 [+0.0069, +0.0240] (0.5611 → 0.5768; 85 targets, 65,045 significant gene–target pairs). Through the analytic
  emitter that is nmae −0.0067 [−0.0110, −0.0024] raw, with a reach cost of −0.030 [−0.067, −0.002] from four small-gate
  targets.
- On three cancer lines held out cross-study (no source from the held line or its study), it reads +0.0033 [+0.0004,
  +0.0062] (HCT116), −0.0002 [−0.0027, +0.0023] (HEK293T) and +0.0010 [−0.0030, +0.0053] (K562) against the production
  analogue.
- The null is not power: re-thinning H1's truth to 100–250 cells, or to the folds' depth, leaves the H1 gain at +0.015 …
  +0.018.
- Mechanism: about +0.017 of the H1 gain is the carriers' whitening. Unwhitened donor signs are better on H1's weak,
  uncalled genes, and the same unwhitened-carrier family inverted nmae on the board (+0.0089 raw, 417c). The remaining
  +0.007 … +0.014 is the consensus, which the cancer folds reproduce at 0.000 ± 0.005.
- Ceilings: an oracle choice of the best single public source per target reaches +0.030 pooled on H1. Supervised stacking
  over seven sources adds ≤ +0.002, and weights learned on the cancer folds invert H1's.
- So "sign sources are at the public-data limit" now covers target-specific combinations too, at these constructions.

*Nor do recipient-side or technical features (517; C-80).*
- The recipient's control identity axis ("genes higher in this line than in others go down") helps on H1, and strongly on
  HEK293T (+0.105 over a cross-line consensus). It reverses on K562 (−0.013), RPE1 (−0.011) and HepG2 (−0.007), and is near
  zero on Jurkat (−0.003) and HCT116 (+0.003). It is a property of some lines, not a law, and nothing in a line's controls
  says which kind of line it is.
- Flex library-size effects are real: a per-gene elasticity of CPM to library size is reliable (split-half 0.97 on H1), and
  strong perturbations shrink their cells' libraries (Spearman −0.44 with response strength). But they are about two orders
  of magnitude below gate-gene effects. They sign H1's gate genes at 0.523 … 0.527, and lower every programme they are added
  to (Arc 0.640 → 0.618; the consensus 0.573 → 0.563).
- Per-gene excess dropout in the controls carries no sign: gate genes are 0.43–0.58 down across expression × dropout cells,
  with no trend. It predicts which genes are called instead: in the lowest expression quintile the call rate rises 0.033 →
  0.067 → 0.109 across dropout terciles. That bears on the call set (§3(l)), not on signs.

*A same-lab programme from another line helps where the design is weakest (517; C-78).* Public screens from two
laboratories (Xaira's HCT116 and HEK293T [ORION]; Weissman's K562 [REPLOGLE22] and three essential screens [REPLOGLE22],
[NADIG25]) reproduce the challenge's situation
without H1. Hold out one line, build the "design" from another laboratory's rows, and take the programme from the held line's
same-lab partner.
- On uncalled genes (the hook's territory) the same-lab programme beats the other-lab design in all three genome-wide cases:
  +0.036 / +0.081 / +0.019 in sign accuracy over every gene, and +0.055 / +0.112 / +0.038 on its most confident 63 %
  (held-out HCT116 / HEK293T / K562). On H1's own line the figure is +0.153.
- HEK293T's programme lowers HCT116's block accuracy (0.488, against 0.523 for the design and 0.557 for another laboratory)
  while raising its uncalled genes' (0.561 against 0.506). A null block prior and a positive hook from one programme: the
  §3(i) prediction, observed.
- On Arc-sized programmes (275 rows, about a third essential), the essential rows carry nothing for another line's
  non-essential genes (0.500 / 0.480). The mean z over the non-essential rows beats the P_down rule (0.548 / 0.573 against
  0.528 / 0.544). That is why 517 estimates the hook's programme that way.
- How much a same-lab partner helps depends on the pair. HCT116's programme predicts HEK293T at 0.595 against 0.541 for the
  best other-lab line, while HEK293T's only ties for HCT116 (0.538 against 0.542); other-lab lines can even anti-predict
  (Jurkat → HEK293T 0.448). Only a board read can give the size for the challenge's producer.

*The same laboratory's programme from another line: the challenge's own test (2026-09-30; C-89).* F18's central open
question is whether a programme estimated in one line by one laboratory carries to other lines that laboratory screened with
the same chemistry. The benchmark's producer released its 2025 CRISPRi screen of H1 on 10x Flex, and screened the 2026 lines
on the same platform.
- The arm: Arc's H1 programme from all public H1 cells (about 1,071 cells per perturbation against 38,176 controls;
  disjoint-control reliability 0.69, C-66; the A/B/C panel's 25 shared targets excluded) replaced the consensus as every
  context's block prior (458d, one variable against 446b: the configuration differs in exactly `sameline`).
- The gate: a second entry gave the programme only to design-strong targets (458g; 214 / 217 / 217 of 300 per context), and
  the rest kept the consensus.
- An exact decomposition: the emitter's random draws do not depend on the prior. Every one of 900 458g blocks equals the
  458d block or the 446b block, and 0 equal neither (468). So 458g − 446b is the programme on the strong two thirds and
  458d − 458g the programme on the weak third, with no re-roll noise between them.
- Preregistered reading (fid raw, 458d − 446b): ≥ +0.0015 carries; +0.0005 … +0.0015 partial; |Δ| < 0.0005 null;
  ≤ −0.0010 harm. Falsifiers: pds and mse unchanged to 1e-6 and jac within ±0.0010. A known call leak (build metas: median 1,
  maximum 1,655 unsuppressible calls per target, against 446b's 0 and 86) was pre-diagnosed to read at −0.001 … −0.004 raw
  jac.
- For scale: on H1's own line the same prior bought fid +0.0198 raw (3 seeds), and sign priors have transferred 0.07–0.67×
  of their H1 fid to the board.

**[[PENDING 09-30 — 458d vs 446b (programme on every target) and 458g vs 458d (the gate): insert the read here. The text for
each preregistered outcome is in `experiments/533-science-lead/RESULTS.md` §4.]]**"""

WHY_OLD = ("unseen line is set by whether anyone has perturbed that line, not by who did it or how. A leader in this challenge\n"
           "appears to train on a proprietary many-line perturbation atlas (public descriptions: on the order of 10^8 cells across\n"
           "200+ lines, with a larger target), which is exactly the resource that would supply the\n"
           "oracle row above.")
WHY_NEW = ("unseen line is set by whether anyone has perturbed that line, not by who did it or how. **[[PENDING 09-30: this\n"
           "sentence stands if 458d reads null or harm. If it reads positive, the clause \"not by who did it or how\" is\n"
           "contradicted: a same-laboratory, same-chemistry screen of another line would then carry a line's programme.]]** A\n"
           "leader in this challenge appears to train on a proprietary many-line perturbation atlas (public descriptions: on the\n"
           "order of 10^8 cells across 200+ lines, with a larger target), which is exactly the resource that would supply the\n"
           "oracle row above. No public source names or screens the final three lines (§3(h)), so for the benchmark's final test a\n"
           "line's programme can come, from public data, only from screens of other lines.")

ANCHORS_B = [
    ("realised |z| order, so an analytic price with the order frozen overstates it.", "after", S3I_ADD),
    (F18_SOURCES_OLD, "replace", F18_SOURCES_NEW),
    ("- The benefit needs the same line and a Perturb-seq-grade measurement.", "after", S3J_NEG_ADD),
    (WHY_OLD, "replace", WHY_NEW),
]
