# 597 — Does 590's 0.65 floor hold on real Flex for the exact failure pair (HEK293 vs HAP1)? (preregistration)

Written 2026-10-01 ~15:50Z by the session "Explore contribution areas [57577c]", **before the data file is downloaded**
(nothing from it has been read). Approved by the coordinator; the download waits for Ahmet's explicit OK. Light only:
no builds, no refits, every process ≤ 4 GB (the RDS is read with 424's streaming parser `rds_stream.py`, not
`readRDS`), outside the coordinator's windows (10-04 12:00Z – 10-05 18:00Z).

## Why
590 found that the day-of identity call fails on lines absent from CCLE: both HEK293 profiles were called HAP1_ENGINEERED
(a registry line) by all eight methods. Their lowest top-1 correlation (`min_top1_score`) was low:
- 0.554 / 0.606 in 3′ chemistry;
- 0.506 / 0.561 after a simulated Flex bias;
- 0.616 in a full-axis Orion context file (595).

The floor (0.65) rests on those values and on real-Flex correct calls: A/B/C 0.772 / 0.719 / 0.733, and VIPerturb's
GEM-X Flex v2 K562 at 0.765 (590 addendum). No real-Flex HEK293 profile has been scored. The final lines are measured
on Flex.

## Data (to be downloaded only after Ahmet's OK; kept in the session scratchpad, never committed)
- File: `multimodal_cell_line_mixing_pilot.rds`, Zenodo record 18460279 (Satija lab, VIPerturb-seq), 1,579,169,440
  bytes, CC BY 4.0.
- Contents (per 596): a pool of three monoclonal KRAB-dCas9-MeCP2 lines (K562, HAP1, HEK293) on 10x Flex v1. A 100-sgRNA
  mTOR library plus non-targeting guides; hashing (HTO) separates the lines.
- Fixation 1 % formaldehyde for 10 min, not 10x's standard protocol.

## Procedure (fixed now)
1. Parse the RDS with `experiments/424-viperturb-k562/rds_stream.py`; record its structure (assays, layers, meta.data
   columns).
2. Select singlet cells whose hashing call names one line AND whose guide assignment is a single non-targeting guide.
   - Use the object's own per-cell columns, as named in its meta.data; the exact column names are recorded in RESULTS
     before any pseudobulk.
   - If the object offers no per-cell non-targeting label, use every singlet of the line and say so.
3. Per line, sum the raw `counts` of the RNA assay (never a normalized layer; if only normalized values exist, stop).
   Drop guide / custom-barcode / ADT / HTO features.
4. Place the sum on `data/challenge/gene_names.csv` (NaN where the probe set lacks a gene) and score it with 590's
   matcher unchanged: `match590.counts_query` → `def_identity`'s 8 methods with the same three anchors and Borda
   consensus. No refit.
5. **Known answer (must pass before anything is read):** the method weights computed here equal 590's recorded
   `method_quality` (`experiments/590-identity-validation/work/match590.json`) exactly. A/B/C re-scored in the same run
   reproduce 590's calls, unanimity and min_top1 (0.772 / 0.719 / 0.733).
6. Report per line: cells used, call, unanimity, `min_top1_score`, consensus top 5.

## Reading (fixed now; the coordinator's rule, made exact)
- **HEK293 `min_top1_score` < 0.65 → the floor is CONFIRMED on real Flex.** If its call is unanimous HAP1, this is the
  exact failure, blocked.
- **HEK293 ≥ 0.65 AND a unanimous call of a registry line → the floor must be RAISED.**
  - Proposed new floor = the midpoint between HEK293's value and the lowest real-Flex correct call among A/B/C (0.719),
    VIPerturb K562 (0.765), and this file's K562 / HAP1 if they are called correctly.
  - It is proposed only if that midpoint sits ≥ 0.02 below the lowest correct real-Flex call. Otherwise report that a
    similarity floor cannot separate them on Flex, and propose a different guard (e.g. no same-line rows on a
    non-anchored line without a human check). Nothing is applied either way; the coordinator decides.
- **HEK293 ≥ 0.65 but split or non-registry → no switch would happen.** Report the value and the neighbour; the floor is
  neither confirmed nor broken by this profile.
- **HAP1 (a registry line) below 0.65, or not called HAP1 → a real-Flex cost of the floor** (a true HAP1 context would
  lose its rows). Report it beside the HEK293 reading.
- K562 is the in-CCLE positive control inside the same object. Below 0.65 or a wrong call means the object's chemistry
  or processing degrades every profile; then the HEK293 reading is reported as uninterpretable rather than as a pass.

## Caveats stated in advance
- Flex v1 with custom barcode probes and a non-standard fixation.
- CRISPRi machinery is present in every cell.
- A pilot's per-line cell counts may be small; the count goes in RESULTS.
- HEK293 (not HEK293T) is the line here; 590's two profiles were HEK293T (Orion) and HEK293 (Xu).

## Columns (recorded 2026-10-01 17:36Z, after the structure inspection and BEFORE any pseudobulk or score)
`work/inspect597.json`. The file parsed with a clean end of stream (two byte-code items skipped by `rds_stream597.py`,
which subclasses 424's reader).
- Assays: RNA (layers counts / data / scale.data), guides, hto, adt; 11,113 cells.
- **Line:** `hash.ID` ∈ {HEK293 4,297, K562 3,437, HAP1 3,379}. Only the three line labels occur, so the object is
  already restricted to hashing singlets.
- **Non-targeting:** `Gene_assignment == "NO.TARGET"` (1,462 cells in total).
- **Counts:** the RNA assay's `counts` layer.
