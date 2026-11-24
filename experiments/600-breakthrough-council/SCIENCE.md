# What the council learned that is scientifically useful (2026-10-02, draft for the manuscript owner)

Five findings, each with the claim, the evidence and its files, how new it is, and its limits. Written by [4d0207] from
the council's own records. These are measurements, not yet board reads: the board arms land from 10-06 / 10-08.

## 1. A shared control pool inflates "programme" oracles in Perturb-seq evaluation (methodology)
- **Claim.** Every target's log fold change is taken against ONE control pool, so the pool's sampling error is common to
  all targets. A leave-one-target-out (LOTO) "context programme" therefore partly predicts each target's realised signs
  from the control pool's noise rather than from biology.
- **Size.** With 400 perturbed vs 18,400 control cells, knowing the control mean's error exactly is worth
  0.5 + arctan(√(400/18400))/π = **0.547 sign accuracy on a gene with no true effect**.
- **Evidence (635, d4cd5f).** h1-test; each gene's controls are split into independent halves, and the truth is taken
  against half B, the oracle against half A. On the non-DE high-CPM genes (fid's "filler"), the LOTO oracle drops from
  0.543–0.553 (sharing controls) to **0.525**, against CLEAN4's 0.519. On reference-DE genes it barely moves (0.625 →
  **0.618**), because real effects dwarf the shared noise. 7772dc raised the issue (7772dc.md §4), and 491's halves
  check (+0.0055 → +0.0024) had already hinted at it.
- **Replicated as a RELIABILITY artifact (462, 2026-09-26; same mechanism, two datasets).**
  - Broad-16 (Cas9 KO, ~500 controls per line): a line's "programme" split-half reliability is 0.801 when both halves
    use the SAME control, and 0.121 with DISJOINT control halves. About 85 % of the apparent programme is the control's
    noise. Honest cross-line sharing is 0.043.
  - Arc's H1 (18,399 controls): 0.793 shared vs 0.386 disjoint on expressed genes. Even a large control pool leaves
    about a third of a weak common programme as noise.
  - Target split-halves cannot detect this; control (or cell) halves can.
- **Why it matters beyond us.** Any benchmark or model-selection procedure that scores against a shared control pool
  overstates how much a "generic response" is learnable for weak or no-effect genes. That includes
  sign-precision-type metrics like vcc2026's fid. A clean oracle needs disjoint control halves.
- **Novelty (literature check 12:35Z): partly known.**
  - AmbiModBench (arXiv 2609.32527, Sep 2026) cites concurrent work (Nicol et al. 2026) that "control-population reuse
    inflates conventional scores". It calibrates with disjoint control partitions.
  - Systema (Nat Biotech 2025) shows models mostly capture systematic perturbed-vs-control variation.
  - What we add: the specific inflation of LEAVE-ONE-TARGET-OUT programme oracles and of SIGN-precision metrics, a
    closed form for its size on null genes (0.547 at 400 vs 18,400 cells), and the split showing it is mostly on non-DE
    genes (filler) and negligible on DE genes.
- **Limit:** measured on one Flex surface (h1-test).

## 2. A line's control-derived identity predicts the SIGN of its DE genes, zero-shot (biology; TEMPERED 13:45Z)
- **Claim (as revised by 7772dc, 4cd40535).** A line's identity I = log1p(ctrl CPM) − mean of other lines' log1p(ctrl
  CPM), computed from its own controls, predicts the SIGN of its called DE genes beyond public priors. Genes the line
  over-expresses tend to be called down. Sign scores from z(CLEAN4) + λ·z(I) beat the best cross-line consensus alone,
  leave one line out. It replicates in CRISPRi (4 surfaces) and Cas9 KO (12 lines, another lab and chemistry). It is not
  a shared-control artifact or an expression-level artifact.
- **Not claimed: a whole-transcriptome "identity erosion" drift across modalities.** A continuous drift toward a common
  state shows only in pooled drug screens (Tahoe), where ambient RNA inflates it. It is absent in single-line KO.
- **Evidence.**
  - 639 (d4cd5f): leave-one-line-out gate-sign gain vs CLEAN4 at λ 1 is h1 +0.016, hct116 +0.015 and hek293t +0.024
    (CIs exclude 0); k562 is a low-power null. It survives the control-halves check (0.5884 vs 0.5887 full).
  - 622 (7772dc): only Flex references work; 3′-only references flip it negative on hek293t. I is a line property
    (Spearman 0.9994–0.9997 across disjoint NTC-guide halves).
  - 622g: independent replication on the Broad-16 (another lab, Cas9 KO, 3′, same-lab references): mean +0.009 over
    12 lines at λ 0.5, 0 lines significantly negative. halves622g: the KO gate-sign gain survives control halves
    (+0.0093 vs +0.0091 full). It is identity, not expression level: with I made orthogonal to expression, +0.0052 of
    the +0.0081 remains.
  - **Where the gain sits is NOT settled.**
    - 639d (d4cd5f; HCT116 and HEK293T, the two Orion lines) puts ≥ 0.83 of the |l|-weighted net on the CPM BLOCK.
    - 622m (7772dc; Broad-16) puts it on the UNCALLED genes: block-restricted λ 0.5 mean −0.0022 (7/12 lines negative,
      min −0.0265, 2 CIs < 0) vs uncalled +0.0204 (0/12 negative).
    - The Broad-16 negatives are mostly epithelial (6/9), including cervix and head & neck look-alikes, which is what
      contexts B (HeLa) and C (CAL33) are. There are no blood lines in the set.
    - So the block localisation rests on the two Orion lines only.
    - **Settled by 609 (4d0207; KAs exact against both analyses): the Orion block result is an artifact.**
      - On HCT116 / HEK293T the block's gate genes are 72–74 % DOWN, from the depth skew of Orion truth.
      - An expression-only rider earns as much block gain as identity or more (+0.047 / +0.084 vs +0.040 / +0.033).
      - With expression partialled out, HEK293T's block gain is −0.006.
      - On the clean surfaces the gain sits on UNCALLED genes: h1 share 0.69, block +0.006; Broad-16 under production's
        sign rule share 0.92, block −0.001.
      - The all-gate-gene sign result survives the identity⊥expression control (h1 block +0.008 / uncalled +0.013;
        Broad-16 all +0.005). The block localisation does not.
      - Methodological lesson: any Orion-fold claim about HIGHLY EXPRESSED genes needs an expression-only control.
    - **The deliverable part is the LOW-identity tail going UP (622n, 7772dc).** On uncalled gate genes, the bottom 10 %
      of I (genes the line expresses BELOW other lines) are truly UP 0.61 / 0.71 of the time on H1's real emitter
      (f 1 / 0.35, against today's delivered 0.52), and 0.78 on the Broad-16 (12/12 lines; CLEAN4 0.62 on the same
      pairs).
      - The top-I tail goes DOWN on H1 (0.58–0.60) but NOT on the KO lines (0.30), so there is no identity DOWN rule.
      - Read as biology: perturbation pulls a line's expression back toward the cross-line mean, and the UP half of
        that regression is the robust one across modalities.
      - Delivery: 630's rank-neutral UP pass is the route that gives these genes real magnitude (630 × identity
        compounding, −0.0105 vs −0.0080 raw nmae at f 1).
  - Delivered value (612, 21dd4c). On H1's real emitter (3 seeds), the key-gated identity hook ties the plain CLEAN4
    hook at board-like regimes (nmae −0.0002 / +0.0007 raw at f 0.5 / 0.35; reach +0.004). No held-out instrument is
    validated for hook-class levers (612 CORRECTION 2), and the fold-based dominance claim is withdrawn. The sign science
    above is unaffected.
  - 626: robust to a sister reference.
  - 628 (7772dc, 7d084ec4): the effect also shows INSIDE the public CRISPRi banks' own screens. A bank's P_down
    correlates with its own line's identity (against 3′ references; own-line and sibling profiles excluded):
    - K562 GWPS r +0.16;
    - K562 essential r +0.28;
    - Nadig Jurkat r +0.13;
    - HepG2 ≈ 0.
    That adds another lab's genome-wide screen to the evidence. The derived lever is negative: de-identifying CLEAN4
    does not transfer better leave-one-line-out (hct116 −0.009, hek293t +0.006, h1 −0.002). The likely reason is that
    "identity" against a reference mean is partly shared across cancer lines.
- **Drugs (627, 7772dc): a drift that does not carry to KO.**
  - Tahoe (50 cancer lines × 1,137 drug-dose treatments; I built leave-plate-out; the drug's cross-line mean removed):
    the line-specific deviation vs −I has Spearman +0.147 (positive in 70 % of 56,777 pairs; wrong-line null +0.013).
    It grows with response strength (+0.085 → +0.231 by quintile).
  - The ambient-RNA fingerprint (genes absent in the line gaining reads) is present but 4.5–6× smaller than the
    identity-gene drop. Roughly ½–⅔ of the Tahoe effect is beyond ambient RNA (semi-quantitative).
  - **broad627: the same statistic on the Broad-16 KO screen (half-split or other controls, 12 lines × 87 targets) is
    −0.019 (null −0.003), with no growth with strength.** The continuous drift does not replicate in single-line KO.
- **Mechanistic reading (hypothesis, narrowed).** Among genes that respond, those carrying the line's
  lineage-specific expression tend to fall. That is a sign bias on responders, not a global regression toward the
  cross-line mean.
- **Novelty.** Partly known. The July 2026 preprint "Perturbation response decomposition enables biologically aligned
  generalization…" (bioRxiv 10.64898/2026.07.24.740459, four CRISPRi Perturb-seq screens) finds the GLOBAL response
  component low-dimensional and inferable from control expression; 441 PS6 saw the cell-line template partly along the
  control identity axis. What is new here:
  - a leave-one-line-out test of SIGN prediction on DE genes across labs, chemistries and modalities;
  - the chemistry sensitivity (Flex-only references);
  - the negative result that the continuous drift seen in pooled drug screens is absent in single-line KO;
  - a deployable zero-shot form (604 v3).
- **Limit.** The effect is modest (+0.01–0.04 sign accuracy). It is a component of the response, not the response, and
  its delivery route (block vs uncalled genes) differs between surfaces.

## 3. A lab's shared programme transfers across its cell lines (biology / experimental design)
- **Claim.** The target-agnostic response pattern a lab's screens share (protocol, timing, selection) is largely
  lab-specific. Learned on 3 lines, it transfers to the lab's other lines.
- **Evidence.** 622h (7772dc), Broad-16: a 3-line same-lab programme gives +0.054 gate-sign accuracy on 12 / 12 held-out
  lines, cross-batch too. A 16-module compression keeps +0.020. This is consistent with 450h's within-lab r 0.54 vs
  cross-lab 0.21.
- **Consequence.** In a multi-line challenge run by one lab, the validation lines carry information about the test lines'
  shared programme. This motivates 625's A/B/C module readout.
- **Limit.** Shown on one lab (Broad-16). Arc's own cross-year transfer (H1 2025 → A/B/C 2026) HARMED at reduced dose
  (458d). Same lab and same campaign may be the condition.

## 4. Mean-based and rank-based DE metrics can be decoupled (benchmark design)
- **Claim.** vcc2026's nmae reads arithmetic-mean fold changes, while calls, fid, jac and reach read Wilcoxon ranks. A
  submission can move nmae with every other member exactly fixed: raise one carrier per gene above the control maximum,
  paid by ballast donations, with per-gene totals and libraries exact.
- **Evidence.** 630 (d4cd5f), real production emitter on H1: Δnmae −0.0066 raw at f 1; fid / jac / reach / pds / mse
  exactly 0; 0 call flips; max ΔAUC 1.0e-4. The board shows many teams' pure-nmae steps (610: −0.017 … −0.070 raw), mostly
  catch-up from worse bases.
- **Why it matters.** It is a lesson for metric suites. Mixing mean- and rank-based statistics lets a submission
  optimise one without improving the biology the other measures. A benchmark could read lfc from the same
  rank-robust statistic, e.g. a Hodges–Lehmann shift, to close it.
- **Limit.** This is scorer mechanics. It does not make the model better at biology, and we say so in every record.

## 5. Technical skews in evaluation truth can masquerade as gene-class biology (methodology)
- **Claim.** Where an evaluation surface's truth carries a technical skew, such as a depth mismatch between perturbed and
  control cells, an analysis that localises a predictor's gain to highly expressed genes can report the skew rather
  than biology. An expression-only control separates the two.
- **Evidence.**
  - `docs/HARNESS-DEPTH-CONFOUND.md`: a harness whose perturbed cells were 28 % shallower than its controls called 91 %
    of DE genes DOWN. Re-thinning inverts the skew.
  - 609: on two Orion-derived folds, the CPM block's gate genes are 72–74 % DOWN (0.51–0.56 elsewhere).
    - An identity-based sign rider looked localised to that block (+0.033 … +0.040).
    - An expression-only rider earned as much or more (+0.047 / +0.084).
    - With expression partialled out, one fold's block gain vanished (−0.006).
    - On a clean same-chemistry surface (H1) and a 12-line KO screen, the gain sat on low-expression uncalled genes.
  - 609b: the skew has TWO tails. On the same folds the LOW-expression uncalled genes skew UP (truth 0.83 / 0.68), and
    an expression-only "bottom decile → UP" rule earns +0.26 / +0.15. That fits CPM normalisation of shallower
    perturbed cells inflating low counts. On clean H1 the same expression-only control is flat (0.580 = CLEAN4), while
    identity's bottom decile gives +0.082, so H1's low-identity UP signal is identity, not expression.
- **Why it matters.** Benchmarks and model analyses that stratify by expression ("the model helps on highly expressed
  genes") should check the truth's sign balance by stratum and run an expression-only control first.
- **Limit.** Two folds from one source dataset. Part of the skew is our own harness construction (thinning / batch
  composition), not the public data.

## Also worth a line in the manuscript (negatives with mechanisms)
- 629 (7772dc): Arc's H1 programme does not transfer to other lines, and identity erosion is NOT why.
  - corr(programme, H1 identity) is only +0.16.
  - De-identified, detrended or low-|I|-restricted versions stay 0.06–0.17 below CLEAN4 on HCT116 / HEK293T / K562,
    and below CLEAN4 on 0/12 Broad lines.
  - On other lines it is near a coin flip (0.51–0.54): uninformative, not reversed. It is worst on genes H1
    expresses below the cancer lines.
  - The identity principle explains a target line's OWN deviations, not why another line's whole programme fails
    to carry over.
- 601: the compositional (CPM closure) filler rule is real (Spearman 0.63) but worth +0.006 at oracle.
- 634 / 638 / 621: public gene-level, responder-weighted and supervised sign sources all lose to CLEAN4 leave-one-line-out.
  The held-out surface's base rate P(down | DE) shifts, and no cross-surface learner can know it.
- 621c: reach is order-limited given sign knowledge (oracle +0.03 overall), but no public signal knows which calls are
  right.
- 623 / 624: identity erosion is not a per-target bulk component. A Testh-style ridge magnitude is real but small
  (+1.3 % aggregate cosine).
