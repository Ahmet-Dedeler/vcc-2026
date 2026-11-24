# 320 — RESULTS: the novelty audit F16 never got

Reading and citation only. **No compute** — no `uv run python`, no `vcc`, no `cell-eval2`, no build, nothing
written under `submissions/` or `data/`, nothing written outside `experiments/320-f16-novelty/`. Every source
with its URL, provenance code and the exact sentence relied on: `work/refs320.md`. The audit's own known-answer
rows (reproductions of 301/307/309's published numbers against the files that hold them, 6/6) are
`work/refs320.md` §G.

---

## THE LINE

**320 (F16 novelty): NEW 0 / EXTENDED 1 / KNOWN 4 / CONTRADICTED 0.**

The phenomenon, the rule, the rule's poor predictive precision, the partner-guide control, and the
aggregate-beats-per-item statistical logic are **all published**, three of them in a paper that used the
**same parent library** two years before F16. What is not published is the specific design — construct as the
unit, with a null that resamples candidate *sets* matched on a *measured* library-wide off-target propensity —
and that is a recombination of two published methods, not a new one. **No scientific claim in F16 survives as
NEW.** The one thing in its neighbourhood I could not find reported anywhere is the Horlbeck supplementary-table
escape defect, and that is data hygiene, not science.

---

## 1. The five items

| # | item | verdict | the prior report it repeats / extends | citation | one-sentence justification |
|---|---|---|---|---|---|
| 1 | dCas9/CRISPRi seed-mediated off-target **transcriptional effects exist** | **KNOWN** | Off-target activity in CRISPRi is pervasive and mostly explained by PAM-proximal seed complementarity, shown with screens **+ RNA-seq + dCas9 ChIP-seq occupancy**, in the **Horlbeck/Weissman genome-wide CRISPRi library** — i.e. our constructs' parent library; independently systematised across three genome-wide Perturb-seq datasets | Rohatgi *et al.*, *Cell Genomics* **4**(11):100693 (2024), doi:10.1016/j.xgen.2024.100693; Hartman *et al.*, bioRxiv 10.64898/2026.03.27.714658v2 (2026) | F16 measures the transcriptional shadow of a phenomenon already demonstrated at the level of **binding** (98–1,141 reproducible dCas9 peaks per guide, 10–33 % promoter-proximal) **and** expression, in the same library family, two years earlier. |
| 2 | a **12-nt seed + NGG match in a promoter window predicts repression** of that gene | **KNOWN** | The rule is Rohatgi's ("10-12 nucleotides… at the 3′-end", PAM-adjacent, promoters "+/- 2000 bp centered around the TSS") and Hartman's (10–12 PAM-proximal bases, **≥12 bp** for high-confidence calls, ±1–2 kb of TSS); **both evaluated it predictively and both report low precision** — 21 of ~562 seed-matching sgRNAs at the GSDMD promoter behaved like on-target guides (≈ 3.7 %), and 47.2–59.8 % of guides carry a ≥12 bp promoter match against 0.83 % (0.2–1.8 %) with a nominated event (≈ 1.4–1.8 %) | Rohatgi 2024; Hartman 2026. Mechanistic basis for the seed: Jost *et al.*, *Nat Biotechnol* **38**:355–364 (2020) **[S]** | The rule, its window, its 12-nt threshold **and its measured predictive value** are published; F16 inherits the rule via 301 and changes only the readout (an aggregate shift rather than per-event nomination). |
| 3 | **construct as the experimental unit**, with a joint null that resamples candidate **SETS** matched on measured off-target hit probability, promoter GC, TSS count/density, mappability and basal expression | **EXTENDED** | Two published components, recombined: (a) the aggregate seed-shift test — "the distribution of log2(fold changes) for transcripts with seed matches… shifted to the left… compared to a background set", one-sided KS, **gene as the unit, no covariate matching**; (b) covariate-matched null loci as a method — propensity-score / stratified matching, **one** matched null set, no repeated-draw empirical p. The unit-of-resampling framework (competitive vs self-contained; gene- vs subject-sampling) is older still | (a) Cazares *et al.*, *Bioinformatics* **40**(1):btae011 (2024), crediting Grimson 2007 / Janas 2018 / McGeary 2019 / Schlegel 2021-22; lineage Birmingham *et al.*, *Nat Methods* **3**:199–204 (2006) **[X, abstract only]**; (b) Davis *et al.*, *Bioinformatics* **39**(5):btad197 (2023); (c) Goeman & Bühlmann, *Bioinformatics* **23**(8):980–987 (2007) | Every ingredient is published and I found no source combining them — so this is a recombination applied to a new substrate (a *measured* 209,070-spacer propensity as the matching covariate, the construct as the unit), which is EXTENDED, not NEW; and see §3(d): the headline p is a **gene-sampling** p, not a construct-level one. |
| 4 | the **partner-construct specificity test** — the same candidates read in the cells of a different construct targeting the same gene | **KNOWN** | The founding observation of the seed off-target field: "few genes regulated in common by different siRNAs to the same target gene"; "the vast majority of the transcript expression patterns were siRNA-specific rather than target-specific"; and their own prescription that "incorporating multiple siRNA duplexes… will increase the confidence with which an observed phenotype… can be linked to target gene silencing". Used as **exactly this control in CRISPRi Perturb-seq**: "the TCR inactivation phenotype was observed exclusively for the APPL2 and LRBA guides… and not for the remaining four guides… **despite repression of their intended targets**"; and as the *discovery device* in Rohatgi ("certain guides displayed a high fold change that deviated from the other guides designed to target the same gene") | Jackson *et al.*, *Nat Biotechnol* **21**(6):635–637 (2003), doi:10.1038/nbt831 **[read in full]**; Hartman 2026; Rohatgi 2024 | The control is standard practice, is 23 years old, and is already in use on CRISPRi Perturb-seq data for this exact purpose; ours is its systematic paired application across all 34 constructs — a quantification, not an idea. |
| 5 | the **aggregate/per-pair distinction** (301/305 find nothing per pair; F16 finds an aggregate) | **KNOWN** | The textbook motivation for aggregation testing: per-item tests are underpowered *and* carry the multiplicity burden, so aggregation both accumulates signal and reduces the number of tests — with the standard caveats that a set-level rejection is evidence for the **set** and not any member, that results are contingent on the pooling definition, and that burden-style statistics assume same-direction effects. 305's specific pathology is not power but **resolution**: a permutation/rank test "gives no basis to distinguish between two gene sets that both hit the floor p-value of 1/(M+1)" | Goeman & Bühlmann 2007 (gene-set version); the rare-variant aggregation literature, e.g. Lee, Abecasis, Boehnke & Lin, *Am J Hum Genet* **95**:5–23 (2014) **[S]** | The tension is not a phenomenon needing explanation — it is the expected behaviour of the design change 307 deliberately made, and the standard treatment (aggregate, then read the result as evidence for the set) is exactly what 307/309 did. |

---

## 2. What, if anything, survives as genuinely new

**Nothing biological.** Every biological content-word in F16's headline — seed, promoter, repression, follows
the protospacer not the target — names something a published paper already reports, and two of the three
sources did it on the same parent library (Rohatgi 2024) or on genome-wide CRISPRi Perturb-seq of the same
kind (Hartman 2026). F16 is a **re-measurement of a known phenomenon on a new sample**: the challenge's own
2025 library. That is worth having, and it is worth saying that it is what it is.

Three things do survive, none of them a discovery:

1. **The design, as a recombination (item 3, EXTENDED).** The construct-level unit with a joint null over
   covariate-matched candidate *sets*, and in particular the use of an **empirically measured** off-target hit
   probability — all 209,070 hCRISPRi-v2 protospacers rescanned through the identical rule — as the matching
   covariate rather than a modelled off-target score. `matchRanges` publishes the matching; the RNAi seed
   literature publishes the aggregate shift; nobody I found publishes them together with the perturbation as
   the unit. This is the only column where our version is ahead of the literature, and it is a methods point.
2. **A defect in a published resource.** The `="SYMBOL"` Excel escape in Horlbeck 2016's supplementary tables,
   which silently stopped an own-target exclusion from firing and moved the per-spacer L = 12 hit rate from
   **0.107 to 0.597**. The *class* of defect is documented and even has a named origin — that escaping style is
   what the Escape Excel tool emits (Welsh *et al.*, *PLoS One* **12**:e0185207, 2017), and downstream parsers
   reading the raw text see the wrapper as part of the value. **This instance, in this table, with a measured
   downstream consequence, I could not find reported anywhere.** It is a data-hygiene note, and F16 already
   scopes it as one.
3. **A scope fact, not a novelty claim.** The challenge's own library has never had its specificity assessed.
   Arc's data notes describe a dual-guide construct, protospacers taken from Replogle 2022, and knockdown QC
   ("83% of cells with >80% knockdown") — and **"off-target effects, specificity filtering, and off-target
   impact assessments are not mentioned"**. The parent library's only off-target control is a *design-time
   alignment filter* with no measured rate ("96.6% of sgRNAs incorporated passed at this threshold"; Horlbeck
   2016 performed no measurement of off-target transcriptional effects). So F16 is the first measurement of
   seed off-target burden on this library. "First on this dataset" is EXTENDED in 287's rubric, not NEW.

**One place F16 is too hard on itself.** §7's point 2 calls the non-monotonic seed-length ordering (L\* = 11
carrying the largest aggregate, L\* = 12 a smaller one) "the single most awkward fact for a simple seed-binding
reading". The source that established the 12-nt threshold hedges it — "complementarity that extends as far as
12 bp of seed sequence **can sometimes** be crucial" — and reports directly that "there must be factors in
addition to 9-nt seed sequence complementarity that determine the binding stability of a sgRNA to genomic
DNA", having found sgRNAs with seed matches close to the TSS **and** in open chromatin that showed no effect.
S6's ordering is therefore *consistent* with the published account, not awkward against it. The awkwardness is
against a strict monotone-seed model that the literature does not hold.

**What would settle the one open tie.** Item 3's EXTENDED rests on a negative search: I could not find the
construct-as-unit matched-set design published, in this field or in RNAi/gene-set methodology. What would
settle it is a targeted methods-literature sweep for perturbation-as-unit competitive tests with matched nulls
(or simply asking the Hartman or Rohatgi groups whether they know of one). It is not worth an experiment.

---

## 3. Wording the audit shows is too strong

Each item quotes the current text, then gives the weaker sentence that is true.

**(a) The partner-construct control is claimed as ours. It is not.**

> `PAPER.md` §3(g), the F16 marker: "What is ours is a preregistered aggregate test with the *construct* as the
> experimental unit on this particular library, a selection-matched joint null with an attainable floor three
> orders of magnitude below α, **and a partner-construct control.**"
>
> `PAPER.md` §4.1: "what we have not found done before on this library is the construct-as-unit aggregate test
> with a selection-matched joint null **and a partner-construct specificity control**"

**True, weaker:** "What is ours is a preregistered aggregate test with the construct as the experimental unit
on this library and a selection-matched joint null over candidate sets. The partner-construct control is not
ours: reading the same candidate genes in the cells of a different guide against the same target is the
standard specificity control of the seed off-target field — it is Jackson *et al.*'s 2003 founding observation
that different siRNAs against one transcript give *siRNA-specific rather than target-specific* signatures, and
it is used on CRISPRi Perturb-seq for this exact purpose by Hartman *et al.* 2026. What is ours is only its
systematic paired application to all 34 constructs."

**(b) "The design's best idea" claims authorship of a standard control.**

> `docs/FINDINGS-2026-09.md` §F16: "**The specificity test is the design's best idea:**"
> `PAPER.md` §3(g): "*The partner-construct test is the design's load-bearing element.*"

**True, weaker:** keep "load-bearing element" — that describes its role in *this* design and is accurate.
Replace "the design's best idea" with: "The specificity test carries most of the design's weight, and it is a
standard control rather than an idea of ours."

**(c) The headline presents a 2003 result as a co-equal finding.**

> `docs/FINDINGS-2026-09.md` §F16 heading: "…and **the effect follows the PROTOSPACER rather than the intended
> target**"

**True, weaker:** "…and, as the seed off-target literature has held since Jackson *et al.* 2003 and as Rohatgi
*et al.* 2024 report for CRISPRi, the effect follows the protospacer rather than the intended target." The
observation is the field's premise, not our result; ours is that it holds in this library at this effect size.

**(d) The headline p is a gene-sampling p; the resolution claim is about draws, not units.**

> `docs/FINDINGS-2026-09.md` §F16: "Attainable floor p = 1.0e-5, **three orders below α**"
> `PAPER.md` §3(g): "With **B = 100,000** joint draws the smallest attainable p is **1.0 × 10⁻⁵**, three orders
> of magnitude below α = 0.05."

Both sentences are literally correct about Monte-Carlo resolution, and 309 already states that the repair cost
one independent check. The audit adds what neither says: **the joint null resamples genes, not constructs.**
In Goeman & Bühlmann's terms it is a *competitive, gene-sampling* null — "gene sampling uses genes as
independent sampling units" — and they show such nulls are anti-conservative under residual dependence between
the resampled items ("the use of the hypergeometric test is anti-conservative; it may greatly understate the
true P-values if the genes in the gene set are not independent", with simulated rejection rates "up to 50,000
times the nominal level"). The covariate matching removes the confounders it measures; it does not make genes
exchangeable.

**True, sentence to add wherever the 2.0 × 10⁻⁵ is quoted:** "That p is a gene-sampling p — the resampled unit
is the candidate gene, not the construct — and its three-orders-below-α resolution is a property of the
100,000 draws, not of the 34 units. The construct-level readings are the construct bootstrap interval
**[+0.32419, +1.07351]**, the paired partner-construct Wilcoxon **p = 0.00459**, and, in C8 only, the
wrong-recipient construct permutation at its floor **p = 2.0 × 10⁻⁴**; under C9's per-construct centring that
permutation is the same object as KA-4R and is no longer independent, so C9's construct-level inference rests
on the bootstrap interval and the partner Wilcoxon. Those, not 2.0 × 10⁻⁵, are the numbers a construct-level
null supports." (This weakens nothing about the conclusion — the interval excludes zero and the paired test is
significant — it corrects which number carries the weight.)

**(e) The occupancy assay is offered as if the mechanism were unmeasured.**

> `docs/FINDINGS-2026-09.md` §F16: "The decisive next experiment … is a **sequence-mutant control** … or a
> **dCas9 occupancy assay** ($5–10 k)." (Same framing at `docs/READ-FIRST.md` §6 and `PAPER.md` §3(g).)

**True, weaker:** "The decisive next experiment is the sequence-mutant control — it is the only design here
that separates seed-mediated binding from a construct artefact that travels with spacer identity. A dCas9
occupancy assay would add construct-specific confirmation rather than the mechanism: dCas9 ChIP-seq on
seed-off-target guides in this library family has already been run and reported (Rohatgi *et al.* 2024: four
guides, 98–1,141 reproducible peaks each, 10–33 % promoter-proximal, with 12 bp of seed complementarity
sometimes crucial for strong occupancy)."

**(f) A comparability note the propensity paragraph needs.**

> `PAPER.md` §3(g): "matched jointly on the *measured* probability that a real library protospacer hits that
> gene off-target under the same rule — obtained by rescanning all **209,070** hCRISPRi-v2 protospacers, not
> modelled"

Nothing here is wrong, but the number invites a false comparison. **Sentence to add:** "This rate is not
comparable with published seed-match rates: our per-spacer L = 12 promoter-hit probability is 0.132353 for the
68 real 2025 spacers and 0.160879 for the pool, against the 47.2–59.8 % of guides Hartman *et al.* report with
a ≥12 bp match within 1 kb of a TSS, because our rule is confined to the CRISPRi active window (TSS offset
−50 … +300) and applies a 10 kb distal filter. Ours is the stricter rule, by design."

**(g) The provenance placeholder is now spent.**

> `PAPER.md` §4.1: "**One result postdates both audits and therefore carries no provenance verdict: F16
> (§3(g)).** We state its status rather than assert its novelty."
> `PAPER.md` §3(g): "Provenance: this result postdates both the 287 audit and its 297 extension and therefore
> carries no audit verdict; under §1.4's rule we say so rather than claiming novelty."

**Replacement:** "Audited by 320: **KNOWN 4, EXTENDED 1, NEW 0.** The phenomenon (Rohatgi *et al.* 2024;
Hartman *et al.* 2026), the seed-plus-promoter rule and its measured predictive precision, the partner-guide
specificity control (Jackson *et al.* 2003; Hartman *et al.* 2026) and the aggregate-over-per-item statistical
logic are all published. What is ours is the design — construct as the unit, with a joint null over candidate
sets matched on a measured library-wide off-target propensity — and a re-measurement of the phenomenon on the
challenge's own library, whose specificity its producers did not assess."

---

## 4. Limits of this audit

Two and a half hours of reading, and the CRISPR off-target literature is large. Concretely:

- **Two sources were unreachable.** The journal version of Rohatgi 2024 (`sciencedirect.com` HTTP 403;
  `cell.com` 403 as 287 also recorded) — every quote attributed to it is from the **bioRxiv preprint**, which I
  read in full, and the two versions could differ. Birmingham 2006 (`nature.com` 303 → IdP) — I quote only its
  EuropePMC abstract and do **not** quote its statistical design; the "seed-match vs background set" idiom is
  attributed instead to SeedMatchR, which I opened.
- **Three claims rest on search-index summaries [S]** and are labelled as such in `work/refs320.md`: Jost 2020's
  seed-mismatch result (corroborating only), the CFD/Elevation framing (used only to say what the off-target
  *score* literature is not), and the rare-variant aggregation background for item 5. Item 5's verdict does not
  depend on the [S] source — Goeman & Bühlmann, which I opened, carries it.
- **Hartman 2026 is a preprint** (v2, 2026-03-30) and its author list came from a search result, not the paper.
- **I searched the literature, I did not exhaust it.** In particular item 3's EXTENDED is a negative search
  result; a construct-as-unit matched-set design could exist in a methods venue I did not query.
- **No re-analysis of our data and no new experiment is proposed**, per the mandate. The five verdicts are
  graded against what each paper's **data shows**, except item 2's "tested predictively", which is graded
  against what Rohatgi and Hartman **report as measurements** (21/562; 47.2–59.8 % vs 0.83 %) rather than
  against their headline claims.

## 5. Files

| file | what |
|---|---|
| `RESULTS.md` | this: the five-row table, what survives, the too-strong wordings, the tally |
| `work/refs320.md` | every source with URL, provenance code [F]/[P]/[S]/[X], access date and the exact sentence relied on; §G the six reading known-answer rows |

Nothing was written outside `experiments/320-f16-novelty/`. No claim, ledger entry or preregistration was
created — this audit grades existing claims C8 and C9, it does not add one.
