# Research programme, 2026-09-05 → (the long game)

## September 9 coordinator update

Current synthesis and next-step decision rules live in
`docs/SCIENTIFIC-REVIEW-2026-09-09.md`. This programme is an experiment history,
not an instruction to launch more workers. The resumed task uses **no
subagents** and **no automatic midnight submissions**, per Ahmet.

**323 completed in the subsequent solo takeover.** The current-family H1
regression panel is rebuilt with the proper target-excluding operator and a
six-metric-validated actual-count evaluator. H1-local anchors and paired
target intervals give 4/4 correct historical point signs but 0/4 resolved
contrasts: not qualified for standalone small-arm selection. The large
target-shuffle control passes. No model promotion or submission. Primary
record: `experiments/323-current-family-validation/RESULTS.md`.
**324 completed:** that source-only feasibility inventory reproduces the
historical coverage but finds no literal two-accession/two-context rectangle
on the panel, and only six targets with same-context/same-method repeats.
Independent publication, samples, controls and constructs are unresolved.
No response fit launched. Next: official-metadata provenance for the ten
same-method accession pairs, especially KHSRP's possible publication-level
crossing. See `experiments/324-residual-label-feasibility/RESULTS.md`.

Round 7's best entry is **319-cap, 0.2186597394**; 319d3 is built but unsubmitted.
317/318 reject their tested H1 constructions, not entire modeling directions.
322 also re-extracted the upstream KOLF atlas: the on-disk panel subset is not
the source's total breadth. F16's novelty/wording is governed by 320/321:
re-measurement, gene-sampling p, no new wet-lab recommendation.

## Current Round 5 interpretation — September 7, 2026

This section and `docs/TAKEOVER-2026-09-07.md` supersede the historical
end-of-Round-4 status below.

- **304 complete and checked:** F15's observed numbers reproduce, but
  orthogonality, a latent-mirror bound, an ORF depth law and a matching-specific
  geometry effect are not established. Its coordinator review also limits
  overstatements in the audit itself. C6 stays VOID and R4 stays spent.
- **305 complete and checked:** two primary sequence-linked repression
  observations, zero full-family discoveries. The post-hoc tiny enrichment
  P value is not accepted as confirmation. Neither failure-specific
  off-targeting nor its absence has been established.
- **303 complete and independently verified:** C7a passes (median −.02633,
  four-cell Holm p=.001999), C7b/c/d fail. R5 was read once after commit
  and is now spent. This is a conditional cross-study/cross-line
  observed-response result, not a universal inverse or a reliability estimate.
  See `experiments/303-second-overexpression-library/COORDINATOR-RESULT.md`.
- **Competition:** no result in 303–305 is a priced production candidate.
  cd4svd5 and the build configuration remain unchanged. 300 rejects the
  tested grid, not all conceivable modelling; its raw-PDS and overall-score
  uncertainty units must not be mixed.

The next priorities are a selection-matched, construct-level
follow-up to 305's sequence-linked candidates and a competition-focused
data/model feasibility proposal. That proposal must identify what is
genuinely different from 224/269/286/289, specify independent-context and
capacity learning curves, and retain the wrong-recipient and shuffled-label
controls. No new experiment is claimed launched merely by listing it here.

Written by the lead scientist on 2026-09-05 13:30Z after nine bounded measurements closed every knob-level
hypothesis around the incumbent (`docs/READ-FIRST.md` v3.7) and after two external reviews (GPT, 2026-09-05)
argued for widening the programme. The funder's brief: **beat the challenge AND produce scientifically
beneficial results; think long-term; parallel experiments for days are fine.**

## The organising question

> **Which perturbation responses are conserved across cell states, which depend on the recipient's state, and
> when can a model know that its transfer prediction is unreliable?**

Everything below serves that question. A better challenge model is one consequence of answering it; a
defensible answer is a contribution on its own.

## Two evaluation tracks, kept apart and both labelled

| track | what counts as evidence | ground truth |
|---|---|---|
| **Competition** | the six `vcc2026` members on the real scorer, production rendering, preregistered slot decisions (`docs/DEF-RUNBOOK.md`, `AGENTS.md` §6) | the board, for the challenge score ONLY |
| **Scientific** | perturbation-specific prediction on INDEPENDENT data; reproducibility across guides / batches / studies; context dependence; uncertainty calibration; response heterogeneity | the held experiments themselves, each "real evidence about its own population" |

A board improvement is evidence of success on that benchmark, not evidence that the model understands the biology;
a scientific result transfers to the challenge only through a measured transfer, never by assumption
(`docs/RESEARCH-METHOD.md` §1, scoped 2026-09-05).

## Rules that bind both tracks (in addition to `AGENTS.md`)

1. **Name the experimental unit.** Study, cell line, biological donor, batch / gemgroup, guide, cell — every split,
   every "replicate" and every σ says which unit it resamples. A per-fold Δ divided by 162's 0.001494 is not a σ
   (READ-FIRST §3).
2. **Reserve confirmation data agents cannot keep inspecting.** Declared reserves: (a) the `orion_hek293t` PLOLO fold's
   split-half B truths (267 builds them; scored only through `scripts/confirm.py`-style logged access, once per
   preregistered claim); (b) four Broad-16 lines named by 268 at its start and never used for fitting or selection.
   Every access is appended to `docs/CONFIRMATION-LEDGER.md` with the claim it served.
3. **Refit preprocessing inside the training split** (operators, centring, gene filters, programme bases). A basis fit
   on all lines and then "evaluated" leave-one-line-out is not held out.
4. **Strong simple baselines in every table**: target-only, recipient-only, equal-weight mean, the incumbent row; and
   the control pair that separates learning from fitting: correct recipient vs WRONG recipient vs shuffled context labels.
5. **Multiplicity and correlated targets**: report the number of cells screened; use paired target bootstraps and
   family-held splits; a per-cell p ignores the screen that selected the cell (265).
6. **A failed feasibility test closes its specified construction, not a field.** Write the closure at that scope.
7. **Scaling is measured on three separate axes**: independent contexts / studies, cells and replicates per condition,
   model capacity. Buy independent information before capacity; let learning curves decide.
8. **Operational**: heavy (≥ 15 GiB) jobs serial on this Mac under the memory gate; genuinely parallel heavy runs go to
   isolated cloud workers with spending caps, immutable inputs and auto-termination (`scripts/cloud_score.sh` pattern);
   no production-builder edit while another experiment reads it; competition safeguards stay durable and are not
   allowed to consume the programme.
9. **Depth-condition every pass rule on a multi-line or multi-study set** (added 2026-09-06 14:30Z after 287): the
   Broad-16 source paper had published its UMI–power law (r = 0.85, elbow ≈ 5 M UMIs per perturbation) before C2a was
   frozen on an unconditioned threshold, and the reserved quarter was the shallow quarter. Stratify or covariate-adjust
   on UMIs per perturbation, and read the source paper's own reliability statistics before freezing a claim on its data.
10. **Audit novelty before claiming it**: every finding carries a provenance line (NEW / EXTENDED / KNOWN /
   CONTRADICTED, with the prior report) before it enters a write-up (E287's table is the template).

## The four workstreams (feasibility stage first, then days for the ones that earn it)

| # | workstream | first deliverable (feasibility, ≤ 4 h) | what earns stage 2 |
|---|---|---|---|
| **267** | **Reference reliability** — the replicate-validated donor-selection advantage | split each held line's cells per target into independent halves (by guide where possible, else batch, else cells, each with its own controls); choose the donor on half A, score on half B; report how much of 263's +0.046 oracle survives, per fold, per 262 group, per reliability tercile | any reproducible advantage → a replicate-aware response model (activity + uncertainty on a modest number of programmes, combined across sources), which must beat the scalar-reliability and shrinkage arms already tested |
| **268** | **Multi-context biology** on Broad-16 (90 targets × 16 lines on disk) | the reproducible context-residual map: shared / line-specific / noise variance components with split-half reproducibility; whether the line-specific part is low-rank (programme-level) and whether it relates to basal state under leave-one-line-out, with shuffled-line and wrong-recipient controls; declare the four reserved lines | a reproducible, programme-level, basal-state-related line-specific component → stage 2 predicts it (LOLO), then tests the prediction on our donors' held lines; 224's conditioned heads are the baseline to beat |
| **269** | **Selective prediction** — a leakage-safe utility router | labels: per-target utility of a FROZEN action set (keep the incumbent row / switch to donor j norm-matched / the deflation on-off) from 262–264's arrays on four PLOLO contexts; features: target biology (DepMap co-essentiality, pathway membership on disk), response reliability, recipient basal programmes; splits: leave-one-context-out AND held target families, reported separately; baselines: fixed actions, shuffled features, the oracle | learned decisions beating fixed actions and shuffled-feature controls on genuinely held data in BOTH settings, with paired bootstraps — not a few winning targets |
| **270** | **Data reconnaissance** — GPSAdb 2.0, PerturbFate, and the multi-context sets on disk (Tahoe, Nadig, scPerturb) | metadata-only audit: coverage of the 300 targets across lines / studies / perturbation types, provenance and terms; one frozen target-level feature (cross-study sign stability) designed and, if metadata allows, computed and tested against 262's per-target `pds` and 267's reliability | a feature that separates conserved from state-dependent responses on held targets → it enters 269's feature set and 268's stage 2 |

The two competition-facing bridges are 267 (does the headroom exist in reproducible biology?) and 269 (can a model
know when not to intervene?). 268 and 270 are the scientific core and feed both.

## What is explicitly NOT on the programme

Finer knob grids around `cd4svd5` (softer deflation, union doses, k ladders, majority heuristics — all closed at their
constructions, READ-FIRST §5); re-running 224 with a larger network; treating an unreleased model (X-Cell, Pisces) as a
dependency; spending a slot on rank or on a candidate that clears only a point forecast.

## Cadence

Each workstream: a feasibility mandate (≤ 4 h, one subagent, measurement only), a written result with the closure or
the go at ITS scope, then a stage-2 mandate with days of budget for the ones that earn it. Slots are spent only on a
preregistered candidate priced on the coherent-pessimistic rule; the board remains the only judge of the competition
track. Progress is logged in `docs/LOG.md`; the live state stays in `docs/READ-FIRST.md`.


## Synthesis after rounds 1–3 (2026-09-06 09:00Z; 267–282, all closed or landed; `docs/FINDINGS-2026-09.md` F1–F12)

**Established (units named, intervals in the findings note):** a single 150–200-cell perturbation response on this
non-essential panel barely resolves a direction (F1); the donor-selection oracle was reference noise (F2); the whitening
operator learns its rows' sampling noise and the retrieval metric needs those directions (F3); `pds` is a threshold at
≈ 0.02 transfer cosine (F4); Broad-16's line-specific component is real, high-rank, not basal-state-predictable (F5);
K562's real DE sets are too sparse for a local DE surface (F6); the board calibrations (F7); cross-line conservation is
a reproducible target property that predicts transfer under LOLO, **robust on the split half** (F8; C1 relabelled 13:40Z — the half-B cells were inside 276's own evaluation, so it is not an untouched confirmation); the
family axis from core machinery to cell identity, **confirmed on 17,106 held genes** (F9); across four human donors a
Flex CD4 response reproduces at +0.054 — the "Flex-400 regime" was an H1 property, but NO cell-budget ceiling is established (F10, corrected 13:40Z: 280 thinned counts, not cells, and its extrapolation mixed cohorts); within one person the across-state loss of agreement is +0.0080 [−0.0026, +0.0186] — not distinguishable from zero — and small pseudobulks are donor × guide construct effects (F11, corrected);
conservation is family-predictable (ρ +0.81) but not gene-predictable beyond reliability (F12; the family split leaked non-first memberships and "not a gene property" was never tested against noise — 284 reruns it).

**Closed at their constructions (competition track):** every knob around `cd4svd5` (257–264), donor selection and
routers (263, 265, 269), noise-aware donor estimation (271), the noise-corrected operator (272), the union's directions
without their magnitude cost (274), class-conditional rendering (277). No priced candidate exists; the slots stay
unspent by rule. The board's own factors (k_pds 1.22–2.48 realised vs PLOLO) say its references resolve MORE than
PLOLO's noisy truths — somewhere between the CD4 regime (+0.05) and H1's (+0.24) — which is the one unknown that
decides how much per-target direction is worth there, and it can only be measured with a board slot spent on a
direction arm whose local reading is unambiguous. None exists today.

**What the programme now says the next investments are:**
1. ~~A resolving reference needs a different cell type, not more cells or donors~~ — **WITHDRAWN 13:40Z** (third external
   review, verified): 280 never measured independent cells, so "more cells cannot help" is not established; **285**
   measures the cell axis, the depth axis and the guide / donor / batch / study terms on ONE matched design, in
   response agreement and in the scorer's statistic. The somatic Perturb-seq sets that would add contexts (Pisces:
   HepG2 / Jurkat) remain unreleased; they are not the only investment. **285 answered it (15:52Z, F14): cells and depth
   are one axis and Spearman–Brown holds; above 76–330 cells per side the binding term is the GUIDE CONSTRUCT (ρ 0.09–0.41
   for CRISPRi), which no amount of measurement of the same construct removes; a second construct, or the challenge's
   own library, is the investment — not more cells or donors of one construct.**
2. **Construct effects bound every one-construct reference** (F11, F1): a second construct per target in any new
   data is worth more than more cells of the first.
3. **The scientific deliverable is the pair of maps and the reliability laws** (F1, F8–F12; `E276/work/conservation276.csv`,
   `E279/work/znf279.csv`, `E280/work`, `E282/work/state282.csv`): a write-up of "which perturbation responses are
   conserved, which depend on cell state, and what a reference can resolve" is supportable now, with C1 and F9's held
   test as its confirmations; the remaining preregistered confirmation to run is F5 on the reserved Broad-16 lines.
4. **Competition track:** hold; re-price only if a new mechanism arrives with a local reading that the reliability
   regime cannot explain away; keep the safeguards durable; spend nothing on rank.
5. **Modelling capacity is NOT exhausted** (review, 13:40Z): what failed are scalar reliability weights and DIAGONAL
   per-gene shrinkage (271), re-weighting the FIXED fitted eigenbasis (272) and fixed programme bases (249–261). A
   full-matrix estimator whose directions come from CROSS-PRODUCTS of independent views (halves, guides, donors), with
   a full-rank residual path, is a different object → **286** (held views → held contexts → the Flex reference; the
   production rendering only after a cleared gate). Negatives close their constructions, not the field (rule 6).
   **286 measured it (15:35Z): DOES NOT CLEAR — 0/12 cells, best −0.0207 pessimistic vs +0.0089; the mechanism
   survives re-estimating the basis from independent views (80 % of the incumbent's whitened prediction energy sits
   where the unbiased signal covariance says there is no reproducible signal; the shared response is itself
   high-rank, 4–11 % of a row's covariance). Open at their scope: target-specific use of Ŝ; a per-target gate for the
   anti-aligned 17.5 % (needs a signal 263 / 264 did not find); the one CD4 cosine cell where the learnt directions
   beat their controls without losing.**


## Round 4 (launched 2026-09-06 13:42Z on Opus, after the third external review)

| # | question | first deliverable | what earns stage 2 |
|---|---|---|---|
| **284** | does 281's family-level result survive a split that excludes EVERY membership of a held family, with component (target-cluster) bootstraps — and is the gene-level remainder a reproducible gene property or noise? | S1 all-membership split and S2 component split next to 281's; the remainder's cross-set Spearman over disjoint source sets with a frozen rule | a reproducible remainder → the map is the necessary object and a per-gene prior is worth building; none → the family axis is the scientific claim |
| **285** | which disagreement between two measurements of one perturbation dominates at the challenge's budget — cells, UMI depth, guide, donor, batch, or study — on matched targets and strata, in agreement AND `pds`? | the benchmark table over h1-test, K562-GWPS (+ the essential screen as the same-line study), Orion, Nadig, Broad-16 and the CD4 bank; the measured cell-ladder depth law; one bounded CD4 cell-level extraction if range-addressable | the dominant removable term names the next measurement (cells vs depth vs a second construct) and bounds every reference claim |
| **286** | can shared response directions learned from cross-products of independent views, with a full-rank residual path, estimate a target's response better than the equal-weight mean — on held donors, held contexts, and the Flex reference? | Stage A (CD4 held donors / held state), Stage B (PLOLO, three cancer folds, priced), Stage C (273's instrument); shuffled-pairing, wrong-source and random-covariance controls; 271/272's constructions as the baselines to beat | a Stage B arm that clears the pessimistic gate 3/3 with controls below half → a preregistered build mandate |
| **287** | which of F1–F12 / C1 / C2 are new, extended, known or contradicted by the source papers and the transfer literature? | the 24-row audit table with citations (`E287/NOVELTY-AUDIT.md`) | landed 14:07Z: NEW 4 / EXTENDED 11 / KNOWN 7 / CONTRADICTED 2 → FINDINGS v1.2, rules 9–10, 288, 289 |
| **289** | do the Broad-16 paper's own state features (p53 status, HPV-E6, lineage, proliferation / p53 programmes, arm CN) predict the line-specific residual where two expression PCs did not; do the paper's mechanistic predictions reproduce on our residual object? | 268's LOLO rerun on the 12 working lines with frozen state features, the three controls, P1 / P2 as preregistered per-target tests, depth as a covariate | a state arm that beats both controls on the residual AND gains on the full response → F5's stage 2 (a state-conditioned model) |

**Round 4 results as they land.** 284 (14:15Z): the family axis holds under folds that share no family with training (+0.805), the ZNF deficit is 57 % predicted (not 84 %), and the gene-level remainder is REPRODUCIBLE across disjoint source sets (ρ +0.328 vs ceiling +0.382) — F12's "mostly measurement" withdrawn; the map is the necessary object (FINDINGS v1.3). 289 (14:36Z): the paper's own state features do not predict the line-specific residual transcriptome-wide (0/25 cells; worse in deep lines), but the paper's P1 reproduces on our data and the LOLO finds it blind (UBE3A / HPV), and p53-axis sign gating beats chance — C-1 reconciled: mechanism per target, no aggregate gain at 11 lines (FINDINGS v1.4). What a stage 2 would need is more LINES, not more features or capacity (rule 7). 290 (15:04Z): 284's reproducible remainder IS carried by the responses' own geometry measured in other lines (ρ +0.317, at the ceiling; annotation −0.057; 40/40 cells) and the preregistered mechanism inverted — broad, spectrally central, proliferation-loaded responses are MORE conserved, programme identity adds little; DNA binders are the exception (F13, FINDINGS v1.5). → 291 (15:30Z; C3 / C4 frozen first): the SIGN confirms on both outside surfaces (`gen_pc20` +0.30 on Broad-16, +0.19 on the CD4-vs-K562 statistic, surviving magnitude), the fitted model does not beat annotation off-map (both claims FAIL by rule; "geometry beats annotation" was a property of 290's remainder target), and the map's own score transfers off-map (+0.27 on 1,773 genes) — F8 extended, F13 narrowed to its sign (FINDINGS v1.6). The conservation track's write-up is now supportable: a reproducible, off-map-predictive map; a family axis that survives strict holds; a reproducible per-gene remainder; a geometry sign confirmed on two outside surfaces. 286 (15:35Z): the cross-product estimator does not clear (0/12) and F3's mechanism survives an unbiased basis — the competition track stays on hold with a measured reason (FINDINGS v1.7). 285 (15:52Z): the matched benchmark — one removable axis (total UMIs, Spearman–Brown), the construct the dominant irreducible term, batch free, study ρ 0.79, line 0.35–0.45, donor ≈ guide in CD4 (F14, FINDINGS v1.8). Round 4's data-generating runs are complete; 292 declared reserve R3 (15:55Z: eight untouched scPerturb files, 174 held-map genes pooled, once-per-claim rule); the library fact is settled from the data (25/300 overlap; the ceiling binds every entry); 293 measured the challenge-regime construct ceiling (16:30Z: ρ_guide +0.075 as measured, +0.687 knockdown-matched — the term is guide efficacy; zero-shot ceiling cos 0.21–0.57, ≈ 5× above the field; ~1 truth in 5 a vector response); 295 spent R3 (16:15Z: C5a PASS +0.502 [+0.153, +0.566], C5b PASS +0.371 [+0.064, +0.415] — the programme's first confirmation on untouched data); 296 tests whether F14's public-library construct term is also knockdown failure; the write-up (294) landed 16:45Z (manuscript draft, 50-row claims ledger, 8 figures) and 297's audit extension graded the post-audit results (NEW 3 / EXTENDED 7 / KNOWN 5): the round's genuinely new content is 284's reproducible remainder, the study-beats-construct inversion and the geometry-vs-identity ablation, with the designs and the preregistered off-map confirmations as extensions; F13's sign and F14's construct term are published (Wang et al. 2026; the CD4 source). The draft is published as an artifact (https://claude.ai/code/artifact/b7545115-c967-4dd0-b670-7329a788f142; regenerate with `experiments/294-writeup/build_html.py`). 296 (17:10Z) settled F14's wording: efficacy explains the construct term in immortalised-line CRISPRi libraries, not in primary T cells; null-corrected, the CRISPRi construct terms are ≈ 0.03–0.08 and only Cas9's survives; failing constructs are reproducible at scale but not a shared vector. Round 4 is complete. The fourth review's competition pilot (299, 20:00Z) — reliability-selected training rows for the whitener — does not clear (−0.0675 `pds`, below every control; the reversed tilt gains +0.0009 3/3 but misses the gate on price and retention); zero slots on 09-06. What the three operator experiments now say together: the rank metric is paid by the donor rows' spread over many directions, and any composition rule that concentrates the covariance on reproducible rows loses — the one composition axis thought to help (062's essential-gene prune) is NOT established on either instrument once read with intervals (300). 300 (22:55Z) closed the reversed-tilt thread as a mechanism: real, monotone, controls passing, and structurally unable to clear (price and amplitude rise together); the whole effect is a third of the ordering-noise SD. **State of the competition track after Round 4:** every construction around the incumbent — doses, sources, renderings, selectors, routers, reliability weights, the operator's basis, its covariance, its training rows — is closed at its scope with a measured mechanism; the measured ceilings (F14, 293) sit ≈ 5× above the incumbent's median transfer cosine, so the headroom is in per-target transfer direction on the marginal band, and no on-disk signal has been found that predicts it (263 / 264 / 269). The next investments are new data (a second construct per target; the challenge's own library in a public screen; Pisces / X-Atlas somatic lines when released) or a new modelling idea with a stated mechanism — not another pass over the same constructions. 301 (23:45Z) closed the write-up's T7 as a measured negative (no seed off-target enrichment among failing constructs) and produced two assets: the 2025 library's spacer table (`E301/work/vcc2025_spacers.csv`) and the resolved JoungZhang TF-ORF atlas (reserve candidate R4: 1.15 M hESC cells, 1,647 TFs, ORF overexpression — a different modality from the map's knockdowns; declare it for one preregistered claim only). 302 (03:50Z) declared and spent R4 on claim C6 (does overexpression reverse knockdown in hESC?): VOID by the preregistration's KA clause (the atlas's count columns are inconsistent), descriptively orthogonal at the atlas's depth with a measured bound; F14's depth law and construct term extend to ORF overexpression; F13's geometry predicts the overlap's size, not its sign (F15). **All four reserves are spent; the next confirmation needs data not yet on disk.**

## Round 6 (launched 2026-09-07 15:30–16:00Z on Opus, after Codex's Round 5 takeover was verified)

Round 5 (303–305) was executed and reviewed by a Codex coordinator (`docs/TAKEOVER-2026-09-07.md`, HEAD `89d38f7`):
C7a PASS on reserve R5 (a bulk-depth second overexpression library shows a **modest** target-matched negative
alignment with knockdown, median cosine −0.02633 [−0.05203, −0.00949], Holm p 0.001999; C7b/c/d FAIL — conservation
and geometry do not predict the reversal), 304's audit **withdrew** F15's stronger readings (no latent-cosine bound,
no established ORF depth law, no demonstrated target-specific geometry effect), and 305 produced **no accepted
multiplicity-controlled discovery** (0/11 primary pairs pass Holm; the design cannot attain the required FWER at its
resolution). All five reserves R1–R5 are now spent. Draft v8 carries the revision.

**The state Round 6 starts from.** Competition: every construction around the incumbent is closed at its scope
(READ-FIRST §5); `cd4svd5` 0.2160944, rank 3/762; leader 0.238932, gap +0.0228. Scientific: the write-up's T7 is
unresolved, and two repression leads survive as observations needing a properly designed test.

| # | question | why it is not a closed family | what earns stage 2 |
|---|---|---|---|
| **306** | what perturbation data and what published method exists **today** that did not on 2026-08-30, and what is publicly knowable about the methods now scoring 0.22–0.24? | 270's reconnaissance is 2 days old and the data freeze is 8; X-Atlas/Pisces was "coming soon" on 09-05; the programme's own conclusion is that the next gain needs **new data** | a released dataset that adds LINES (F5), a SECOND CONSTRUCT per target (F14), or the challenge's own library in a public screen (293) → a costed acquisition |
| **307** | do the frozen T7 seed candidates show aggregate repression, with the **construct** as the experimental unit, against a joint null that preserves seed/PAM, promoter architecture, GC, mappability **and** expression? | 305's per-pair design cannot attain its own FWER (floor 11/158 ≈ 0.0696) and its post-hoc Poisson-binomial was rejected by the coordinator as a rescue; the aggregate test was never run | a bounded aggregate estimate either way; a positive would still not establish off-target binding |
| **308** | does weighting each donor row by the **on-target knockdown its own construct achieved** beat equal weights — and beat reliability weighting? | F14 makes the construct the dominant irreducible term; 296 makes that term guide efficacy; **293 shows a failing construct still yields a reproducible transcriptome**, so every reliability/norm/agreement-based scheme tried (016, 097a, 262/263, 299) is blind to it by construction. 097 route (d) closed |kd| as a *prediction target* (the scorer drops the own coordinate), not as a *weight* | 3/3 PLOLO folds, every control at or below half, the `nmae` leg intact, priced ≥ +0.0089 → a preregistered build mandate |

**Ideas checked and NOT launched** (each already closed, recorded here so Round 7 does not re-propose them):
conservation/geometry as a per-target *lever* (276 limit 4: diagnostic only — it says where the consensus points
away, not how to move it; 263 closed the arbiter); a utility router over {keep, switch donor, undeflate} (269, 0/3
settings, at its shuffled controls); shrink-toward-the-panel-mean (097: loses to shrink-toward-zero 12/12 —
`cos(our mean, the reference's mean)` +0.016/−0.066, capture is cos²); Arc's pretrained STATE checkpoints (182:
NO-GO, cannot be conditioned on a single one of our targets).

**Round 6 results (all landed 2026-09-07 16:00–16:45Z).**

**306 — no acquirable data exists, and the frontier is not using data we lack.** Between the 2026-08-30 freeze and
today **no perturbation dataset touching the 300 targets, the six challenge lines, or the challenge's own library
became downloadable**: X-Atlas/Pisces is still an empty repo (`usedStorage 0`, unchanged 165 days), and GEO,
ArrayExpress, Zenodo and scPerturb return nothing qualifying. Arc has published no data; Oct 22 / Nov 5 confirmed.
Of 770 board entries, **neither team ≥ 0.220 describes a method**, the best describing entry is rank 42 (0.17543),
and **every rival who describes one names only sources already on this disk**, scoring 0.136–0.150 against our
0.2161 — at least four have independently rebuilt the incumbent's own recipe. A rival's public claim of "21/300"
2025-context coverage is not reproducible; our 25/300 reproduces exactly. **With 225's measured saturation of the
public-source channel (+5 % of board cosine), the leader's advantage is method, not data.** Two conditional
acquisition triggers recorded; nothing downloaded.

**308 — the competition arm DOES NOT CLEAR, and the mechanism is refuted at its construction.** 0/16 arms, 51
cells, **64/64 known-answer rows**. The premise held: |kd| is build-time observable on 226/290/294 of 300 targets,
absent from the incumbent's base by construction, not a norm proxy (ρ +0.03), only weakly a reliability proxy
(+0.04…+0.19), with a large acting set (106/300 disagreeing targets) — so this is not a power failure. Best arm
`pds` +0.00077 [−0.00114, +0.00283], priced −0.00145 pessimistic vs +0.0089, **indistinguishable from
entropy-matched random re-weighting of the same rows**, and behind 097(a)'s `normrel^−1.5`. **The measured reason:
on the disagreeing set the higher-|kd| donor's row is NOT the one closer to the held truth** — dropping the failing
donor moves more targets down than up in all six cells by 1.5–2.4×; the variance cost of averaging two rows instead
of three exceeds the bias benefit of dropping a failed construct.

**307 + 309 — the round's scientific result, and the discipline that makes it credible (F16).** C8 measured the
aggregate candidate-repression statistic with the **construct** as the unit over all 34 targeting constructs
against a joint null that resamples candidate *sets* matched on off-target hit probability, promoter GC, TSS count
and density, mappability and basal expression (attainable floor 1.0e-5, against 305's 0.0696). It met both pass
conditions and **aborted itself** on 2 of 7 known-answer rows rather than loosen a threshold — the discipline 305
failed. **309 repaired exactly those two defects and re-tested at the same α: C9 PASS, `S` +0.71293
[+0.32419, +1.07351], p 1.99998e-05, 9/9 rows**, with the bias repair moving the estimate **down 9.5 %** exactly as
declared in advance. 26/34 constructs positive; a real-spacer null agrees; and the specificity test puts the signal
on the **protospacer**, not the intended target (partner-construct difference +0.745 [+0.229, +1.258], Wilcoxon
p 0.00459). Both experiments state the limits in their abstracts: this is **not** evidence of dCas9-KRAB off-target
binding, C8 and C9 are two tests of one hypothesis on one dataset rather than a replication, and the decisive next
step is wet-lab (a sequence-mutant control ~$3–6 k, or an occupancy assay ~$5–10 k).

**Where the two tracks now stand.** *Competition:* blocked, with one more measured closure and no acquirable data;
both 09-07 slots left unspent. The bar is unchanged and quantified in READ-FIRST §1.1 — **+25 % relative transfer
cosine to clear the gate, +48.6 % to match the leader**, against a best-ever arm of +5.16 %. *Scientific:* Round 6
produced the programme's first PASS on a genuinely new question since C5, and it is the finding most likely to
interest someone outside this competition.

## Round 7 — the board loop (2026-09-07 15:30Z → 09-08, after the slot policy was corrected)

**The policy change that started it.** Slots before Oct 22 are prize-irrelevant AND expire nightly, so the
+0.0089 gate — which exists to decide what to BUILD — must never decide whether to spend one. Four had been left
to rot. From 09-07 the loop is: preregister the reading in bins *before* the entry exists, submit the best
available arm, capture the record immediately, read it against the bins.

| # | question | result |
|---|---|---|
| **312 / 312b** | does per-target donor weighting (097 a6 / 308 `normrel1.5`) pay on the real scorer, and where is its dose optimum? | **+0.0016245** at exp −1.5 (0.2177189, our first gain since 09-05); **the curve turns over**: exp −3.0 gives 0.2164511 with `pds` still rising (+0.00202) while `mse`, `nmae` and `reach` all fall — 143's pattern on the real scorer. **The exponent axis is CLOSED at ≈ −1.5.** |
| **313** | is the shipped arm the object 308 measured? | **NO.** PerturbMatch is applied at the uniform share and outside the norm match, so delivered rows run 0.806–2.155 (whole matrix +6.0 %) where 308's arm is exactly shape-only. Max \|diff\| 2.700 against a 3.39e-07 float32 floor. **This explains the `nmae` sign error**: an amplitude-changing object loses `nmae`. |
| **316** | does restoring exact shape-only recover it? | **The preregistered IMPORTANT NEGATIVE fired: 0.2169733, −0.00075.** But **both arms beat the incumbent** — the weighting alone is **+0.00088** and the accidental distortion adds **+0.00075**. `nmae` recovered +0.00339 exactly as predicted, and was outweighed by `mse` +0.00399 and `pds` +0.00225 the distortion had been buying. |
| **317** | can we trade `mse` for `nmae` + `reach`, as the board's other shape suggests? | **NO, 28.6–40.8 : 1 against.** Emission axis empty (best +0.000131 = 7 % of gate, 1.4× a re-seed); amplitude axis cannot clear its gate **even counting `mse` free** (+0.002374 best, +0.005190 oracle). Confirmed on the board's own history: over **1,961 consecutive same-team pairs**, teams that BOUGHT `mse` improved **83.2 %** of the time, teams that SOLD it **33.9 %**. **§4.1c re-established on far stronger evidence.** 094's call-set saturation confirmed on the real six-member scorer for the first time. |
| **318** | is PER-TARGET magnitude redistribution a separate lever from the uniform axis 8.1d closed? | running |

**Four instrument-domain traps, all in one day, all now written down.** (1) A raw-vs-scaled **units error** made a
+0.00184-accurate forecast look like a 4.7× underpricing. (2) 197's `mse`→cosine inversion was applied **across
emission families** where it is not valid — the leader inverts to 0.01546 while their raw `pds_cosine` (0.86832)
is BETTER than ours (0.85706); when the two disagree, `pds_cosine` wins. (3) "Shape-only to 1e-16" held only at an
**intermediate** matrix. (4) h1-test's `mse` optimum sits at amplitude ≈ 0.85 **because that surface is scored with
the H1 union off** (oracle-leaked), so it cannot locate the board's optimum — 187/211 locate it at 0.982, where we
ship. **Three of the four were mine; GPT 6 Astra caught the first, I caught the second, and experiment 317 caught
the fourth before I could act on it.**

**What the competition track now knows that it did not on 09-07 morning.** The per-target weighting mechanism is
**real and small** (+0.00088 shape-only). Its dose axis is closed. The emission and uniform-amplitude axes are
closed on the real six-member instrument. **Buying `mse` is what wins on this board** — measured, not inferred.
The frontier is not using data we lack (306) and the current leader's +0.00725 jump was `fid` alone, closing a
deficit we do not have. **The open question is whether the +0.00075 the accidental distortion earned is a
deliberately reachable per-target magnitude rule (318), or an accident of one profile.**
