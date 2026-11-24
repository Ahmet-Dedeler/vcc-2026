# 320 — The novelty audit F16 never got: is the seed-match repression finding NEW, EXTENDED, KNOWN or CONTRADICTED?

You are a research subagent. Work in `experiments/320-f16-novelty/`. **Edit nothing outside it.** You MAY use web
search and web fetch; record every source (URL, DOI, date accessed) and the exact sentence you rely on in
`work/refs320.md`. No `vcc`, no `cell-eval2`, no build, nothing under `submissions/` or `data/`; `uv run python`
only if a check needs a metadata field (state it and its known-answer row). **Time box: 2.5 hours.**

## Why
Programme **rule 10**: no claim enters the write-up without its NEW / EXTENDED / KNOWN / CONTRADICTED verdict and
the prior report it extends or repeats. **F16 has none.** It landed after 287's audit and after 297's and 304's
extensions, and `experiments/294-writeup/PAPER.md` §4.1 currently says so explicitly — that it postdates both
audits, that the phenomenon itself is published, and that what is ours is the construct-as-unit test. **That
placeholder is not a verdict. Your job is to supply the real one, adversarially, preferring KNOWN wherever a paper
already reports it.**

## The claim to audit, in its exact wording
`docs/FINDINGS-2026-09.md` §F16, and the two experiments behind it:
`experiments/307-offtarget-enrichment/RESULTS.md` (claim C8, `UNRESOLVABLE` — self-aborted on 2 of 7 known-answer
rows) and `experiments/309-c9-calibration-repair/RESULTS.md` (claim **C9 PASS**). In one line: **in the 2025 VCC
CRISPRi library, genes whose promoter carries a match to a construct's protospacer seed are more repressed in that
construct's own cells than sequence-, promoter- and expression-matched genes the rule did not select — construct as
the experimental unit, n = 34, `S` = +0.71293 [+0.32419, +1.07351], joint-null p = 1.99998e-05 — and the effect
follows the PROTOSPACER, not the intended target** (the same candidate genes read at wrong-recipient level in the
cells of the partner construct that silences the same gene with different spacers: difference +0.74500
[+0.22896, +1.25766], paired Wilcoxon p = 0.00459).

## Audit these five items separately — they have different novelty
1. **dCas9/CRISPRi seed-mediated off-target transcriptional effects exist.** Almost certainly **KNOWN**; find the
   primary reports and say how strong the prior evidence is (ChIP/occupancy vs expression; magnitude; how many
   sites).
2. **A 12-nt seed + NGG match inside a promoter window predicts repression of that gene.** Whose rule is this, and
   has anyone tested it as a *predictive* rule rather than annotated it post hoc?
3. **The construct as the experimental unit, with a joint null that resamples candidate SETS** matched on the
   measured off-target hit probability, promoter GC, TSS count/density, mappability and basal expression. Is this
   design published anywhere, in this field or an adjacent one?
4. **The partner-construct specificity test** — the same candidate genes read in the cells of a *different*
   construct targeting the *same* gene. Has that control been used before to separate spacer-carried effects from
   intended-target biology?
5. **The aggregate/per-pair distinction** — 301/305 found no failing-vs-working enrichment per pair, F16 finds an
   aggregate effect. Is that tension itself a known phenomenon (power vs aggregation) with a standard treatment?

Check at least: Jost/Weissman and the CRISPRi mismatch/specificity literature; Horlbeck et al. (hCRISPRi-v2 —
**and note our own finding of an Excel-escaping defect in its supplementary gene labels**); the Replogle 2022 and
Wang 2026 reliability papers already in `REFERENCES.md`; anything the 2025/2026 VCC data notes say about their own
library's specificity; and the off-target prediction literature (seed rules, CFD/Elevation-style scores).

## Rules of evidence
- **Prefer KNOWN.** A phenomenon reported qualitatively still makes our version EXTENDED at best.
- Distinguish what a paper **claims** from what its data **shows**, and say which you are grading against.
- Where a source is paywalled or unreachable, say so; never infer contents from a title.
- **Do not re-analyse our data and do not propose new experiments** beyond naming what would settle a tie.

## Deliverable
`RESULTS.md`: a five-row table (item | verdict | the prior report it extends/repeats | citation | one-sentence
justification), a paragraph on **what, if anything, survives as genuinely new**, any wording in F16 or
`PAPER.md` §3(g)/§4.1 that the audit shows is **too strong** (quote it and give the weaker true sentence), and one
line — **320 (F16 novelty): NEW n / EXTENDED n / KNOWN n / CONTRADICTED n**. Plus `work/refs320.md`. Then stop.
