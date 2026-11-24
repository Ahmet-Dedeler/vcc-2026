# 304 — Complete the interrupted F15 audit

The user asked Codex to take over on 2026-09-07. The original `MANDATE.md` is an
unfinished draft, retained as history. This complete mandate governs this run.

## Scope and operational limits

You own only `experiments/304-f15-audit/`. Read `docs/READ-FIRST.md`, `AGENTS.md`,
`docs/RESEARCH-PROGRAMME-2026-09.md`, `docs/FINDINGS-2026-09.md` F15 and the
results and preregistration of experiment 302. The coordinator owns all shared
documents and the confirmation ledger. Do not edit them, commit, publish, submit,
build predictions, or start other agents.

Use the strongest configured model, public primary sources and read-only source
code inspection. Do not read any reserved expression values: no Nakatake X or
control rows, no Joung atlas X, and no other reserve. Existing 302 summaries,
saved sums and analysis code may be inspected for an audit, not repurposed as
an independent confirmation. Use `apply_patch` for authored files and `uv run
python` for analyses. Stay below 6 GiB RSS; no large downloads. Time box 90 minutes.

303 is preparing an independent second-library test but is NOT authorized to read
reserve expression yet. 305 tests the existing T7 expression candidates separately.
Do not duplicate their computations. Your result does not depend on either result.

## Deliverable

Write `RESULTS.md`, `NOVELTY-AUDIT.md` and `work/refs304.md`. Each audit row needs
the exact local claim and location, primary source with DOI/URL and publication
date, what that source actually demonstrates, a verdict (NEW / EXTENDED / KNOWN /
CONTRADICTED / UNRESOLVED), and proposed defensible wording. Distinguish an
unaudited claim from evidence of novelty. Document inaccessible sources rather
than inventing their contents. Preserve every contrary result.

Audit the following in order:

1. **Population and provenance.** Verify the original cell lines, assay/platform,
   induction times, reference conditions and preprocessing for Joung 2023,
   Nakatake 2020 and Pan/GSE295214. "hESC" is a cell class, not a cell-line ID.
   The Nakatake resource also has a 2021 Scientific Data paper; identify which
   original measurements supplied the SCMG database's Nakatake2020 rows.
   Follow source code and metadata, never a filename alone. Send a concise
   provenance note early if it invalidates the planned same-line interpretation.
2. **F15's orthogonality / lack of mirroring.** C6 is VOID by its frozen KA-E
   clause. Can a non-significant median be called biological orthogonality?
   Compare with primary literature on gain-vs-loss-of-function transcriptomes,
   including Tian CRISPRi/a and the source papers. Separate observed noisy
   alignment from the latent biological alignment and equivalence testing.
3. **Attenuation / ceiling.** 302 substitutes 273's 0.243 into the KD side without
   measuring this bank. Audit sqrt(rho_OE*rho_KD), half/full-depth conventions,
   correlated errors, residualization and using medians of reliabilities to
   rescale a median cosine and its interval. Is the claimed exclusion of a
   strong mirror justified? Provide an assumption-labelled alternative.
4. **The depth law.** Determine whether 302 actually downsampled the SAME TFs
   across depths or compared different TFs in cell-count strata. A cross-target
   association does not by itself establish a within-target learning curve.
   Audit the fitted constant and extrapolated "3,000 cells" requirement against
   prior Spearman–Brown work and the available data.
5. **The isoform-construct term.** Distinct splice isoforms may be different
   interventions; near-zero agreement is not automatically an irreducible
   construct ceiling. Check study designs, centering and shared-cell statements.
6. **Geometry predicts |alignment|.** Audit the descriptive +0.38 result and its
   partials. Check selection, multiplicity, shared KD objects / mathematical
   coupling and whether the association exceeds a target-shuffled alignment
   baseline. "Not a measurement artefact" requires more than those partials.
7. **Narrow, related wording check.** 301's P=0.679 seed enrichment test does not
   establish absence of promoter off-targets. Give a concise correction for the
   sentence "the term is knockdown failure, not off-targeting"; 305 owns the
   expression test, so do not perform it here.

Source search starting points (not assertions of content):
Joung et al. 2023, Cell, doi:10.1016/j.cell.2023.01.026;
Nakatake et al. 2021, Scientific Data, doi:10.1038/s41597-021-00806-1;
Nakatake et al. 2020, Cell Reports 31:107655;
Pan's SCMG_scripts repository and GSE295214. Cite primary text or source code,
not search snippets. Read the full relevant methods before accepting a claim.

Do not turn a literature review into a universal negative: close only claims
the evidence directly contradicts. End with the three most consequential
corrections and the paths you wrote.
