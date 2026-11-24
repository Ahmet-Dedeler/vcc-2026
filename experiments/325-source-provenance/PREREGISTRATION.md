# 325 — Provenance behind the ten apparent same-method accession repeats

September 9, 2026. Codex is performing this bounded metadata pass while the
Claude Code subscription connection awaits user approval. **No model call,
response measurement, matrix download, reserve access, or submission.**

Question: do 324's ten same-target/context/method accession pairs identify
independent samples, controls, constructs and publications? In particular,
could KHSRP form a two-publication/two-context crossing once the separate
ENCODE experiment accessions are grouped by their actual provenance?

Input: the frozen
`324-residual-label-feasibility/work/panel__id_fallback__same_method__pairs.csv`.
Use every pair. Do not select pairs based on response values.

Before source retrieval, commit this protocol, the retrieval/parser code,
tests, and their input hashes. Three official metadata pages were located
through search to establish the endpoints (GSE145984, ENCSR561CBC, and the
GSE123957/GSE123958 pair); no response values were read.

## Allowed retrieval

Official GEO SOFT **self metadata**, not expression tables, for the named GSE
records; official ENCODE experiment JSON for the named ENCSR records.
Follow publication, biosample, library, genetic-modification and sample
metadata references only when needed to resolve units. Never follow a count,
FASTQ, BAM, DE, supplementary-expression or signed-download link.

Sample-level priority is KHSRP (GSE145984 and its two ENCODE experiments),
then the MTA1 GSE123957/GSE123958 pair. Bound follow-up to **80 metadata
objects** total, **30 MiB**, **45 minutes**, no authentication or paid service.
Missing/unavailable fields remain unknown; preserve HTTP errors.

## Interpretation fixed before retrieval

- Distinct accessions alone prove neither independent publications nor
  independent biological samples, controls, constructs, or analysis pipelines.
- Same publication is a provenance grouping, not proof that experiments are
  duplicates. Report experimental independence separately.
- Shared sample IDs, identical data-file identifiers, or explicit reuse links
  are evidence against independent measurement. Disjoint accession sample
  lists alone are not proof that all underlying material is independent.
- Same guide/construct is a construct repeat, not an independent construct.
  Different sequence/construct IDs must be verified, not inferred from names.
- A publication-level crossing may exist even when 324's literal
  two-accession rectangle count is zero. Report this correction if supported.
- Even verified metadata does not establish downstream-response reliability.
  It may earn a separately registered development pilot, never a score or
  confirmation claim. All six competition metrics remain unmeasured.

## Known answers and deliverables

Parsers must reproduce synthetic GEO/ENCODE IDs, preserve repeated metadata
fields, distinguish a shared publication from shared samples, and fail on a
returned accession mismatch. Every successful real retrieval must bind back
to its requested accession and retain URL, UTC time, byte count and SHA-256.

Deliver `work/series-summary.json`, exact source records, a per-pair provenance
table and `RESULTS.md`. Do not modify 324's inputs, the confirmation ledger,
the manuscript or production code.
