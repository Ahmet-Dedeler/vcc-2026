# 325 — Accession repeats do not yet supply independent residual labels

Completed September 9, 2026 by Codex. **COMPLETE_METADATA_ONLY_NO_LABEL_PILOT.**
Nineteen official records (13 GEO series, six ENCODE experiments) were fetched
at 09:35:48–09:36:05 UTC; all returned accessions and receipt hashes verified.
Total source payload **5,330,269 bytes**. No sample follow-up was needed to
resolve the priority KHSRP/MTA1 questions. No expression, count, DE, FASTQ or
other response payload was requested; no reserve or paid service accessed.

## Decision

The promising KHSRP publication-level crossing is **not independent
biological perturbation material**: GEO explicitly describes reuse of ENCODE
KD/control samples. MTA1's pair belongs to one study and changes synchronization
state. SRSF5 is a useful **within-study, distinct-construct** contrast, not an
independent-publication repeat. The remaining pairs are unresolved at the
specific sample/control/construct level.

**Zero pairs are verified here as fully independent same-state biological
repeats.** This is a verification count, not a claim that the unresolved pairs
are duplicates or that GPSAdb's annotations are wrong. The metadata does not
earn a cross-study residual-label training pilot.

## All ten pairs

| target/context | accessions | finding and remaining limit |
|---|---|---|
| EZH2 / K562 | GSE162623 / GSE193486 | Different PMIDs. First explicitly describes EZH2 KO; second describes silencer-region KOs and drugs. Exact GPSAdb group-to-sample mapping and modality unresolved. |
| EZH2 / HAP1 | GSE133183 / GSE285461 | Different PMIDs; first includes EZH2KO, second mixed regulatory-element/CRISPR assays. Exact groups and constructs unresolved. |
| EZH2 / MCF7 | GSE183458 / GSE76271 | Different PMIDs and shEZH2 designs. Exact constructs, controls, samples and material independence unresolved. |
| EZH2 / NHDF | GSE109064 / GSE98768 | Both describe shEZH2; days 4/8 versus stable fibroblast lines. Name-only context, donor/construct identity, and second publication unresolved. |
| KDM5A / HEK293T | GSE179112 / GSE277160 | Different parent series; latter is a broad ZFP knockdown study. Exact KDM5A samples and constructs unresolved. |
| KHSRP / K562 | ENCSR561CBC / GSE145984 | Explicit ENCODE material reuse in GEO; not independent biological perturbation replication. Exact GSM-to-ENCBS mapping not established. |
| KHSRP / HepG2 | ENCSR850CKU / GSE145984 | Same explicit reuse statement, same limit. |
| MTA1 / HCT116 | GSE123957 / GSE123958 | Both SubSeries of GSE123959; nonsynchronized versus synchronized. Shared study, different cell-cycle conditions. |
| SRSF5 / HepG2 | ENCSR447UCG / ENCSR781YNI | Shared ENCODE lab/publications, **different constructs** and sample/control accessions. Underlying material independence is not proven by disjoint IDs. |
| TAF15 / K562 | ENCSR031RZS / ENCSR611ZAL | Shared ENCODE lab/publications and **same construct**, despite different modification, sample and control accessions. |

These are series/experiment metadata judgments. A series containing multiple
modalities cannot adjudicate an individual GPSAdb group without its sample map.
No response-value selection or downstream-reliability claim was made.

## Priority provenance, exactly

GSE145984 (PMID 33287884) says its KHSRP KD/control samples were originally
generated for the ENCODE project by Brenton Graveley's lab, followed by total
RNA sequencing for circular-RNA analysis. The two ENCODE experiments both use
**TRCN0000013254** and cite PMID **32728246** among their references.
Grouping by publications therefore creates a two-context/two-publication
crossing, correcting 324's purely literal-accession framing, **but the
publication units are not independent biological samples**. This could be a
technical/library-preparation comparison; exact material mapping remains unknown.

SRSF5:

- ENCSR447UCG: TRCN0000000141, control ENCSR491FOC,
  biosamples ENCBS551QFT / ENCBS582REX.
- ENCSR781YNI: TRCN0000284799, control ENCSR042QTH,
  biosamples ENCBS762HAX / ENCBS311FWL.

TAF15:

- ENCSR031RZS: ENCGM221YUK, control ENCSR942UNX.
- ENCSR611ZAL: ENCGM895GEB, control ENCSR344XID.
- Both genetic modifications name **TRCN0000020143**.
  Distinct modification accessions are not distinct construct sequences.

## Evidence and verification

- `work/sources/`: raw official metadata plus URL/time/bytes/SHA-256 receipts.
- `work/evidence.json`: narrow, source-bound fields and all ten judgments.
  SHA-256 `b7634403ffc12eced2f62120ebe9991b5cb64fc284b9ab759d0cbccc8ca51178`.
- `work/pairs.csv`: compact table; `work/series-summary.json`: complete parsed
  retrieval output. `report325.py` independently reparses and hashes every source.
- Retrieval protocol/code/input hashes were committed at **40e4c547**, before
  source retrieval. Five parser tests pass.
- The first report assertion omitted the word "the" from the source's literal
  reuse phrase and failed before writing outputs. Corrected to the exact phrase;
  no source or scientific interpretation changed.

The 45-minute bound was respected. No broad metadata expansion or response
pilot followed. This finishes 324's priority provenance question without
claiming that all possible external independent labels have been exhausted.
