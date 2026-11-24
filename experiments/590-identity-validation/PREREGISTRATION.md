# 590 — Does the day-of identity call switch same-line rows on only when it should? (preregistration)

Written 2026-10-01 ~13:40Z by the session "Explore contribution areas [57577c]" (numbers 590–599, agreed with the
coordinator and the 550–559 peer), **before any validation query was scored**. Analysis only: no build, no `vcc`, no
`cell-eval2`, every process ≤ 4 GB. Class (AGENTS §11.1c): **scorer mechanics / day-of procedure** — it changes nothing
fitted; it measures a gate.

## Why
On Oct 22, `scripts/def_identity.py` matches each new context's control pseudobulk to CCLE (041's arithmetic: 4 feature
sets × Pearson/Spearman = 8 methods, anchor-weighted Borda consensus). `scripts/def_sameline.py` then switches same-line
rows on for a context when the consensus call is a REGISTRY line (12 lines) AND the top-1 is identical under all 8 methods.
A false unanimous call injects another line's rows at w 0.5 (290 rows for HCT116, 272 for K562). The matcher has been
checked on 6 profiles only (A/B/C + the K562/HCT116/Jurkat anchors). This measures it on every labelled line profile on
disk.

## Instrument
Queries: control pseudobulks (log2 CPM over the official 18,533-gene axis ∩ CCLE, as for a D/E/F context), one per line
and dataset:
- McFarland 2020 MIX-seq (10x 3'): up to 209 CCLE lines, control cells (DMSO `control`, `sgLACZ`, `sgOR2J2`), ≥ 50 cells;
- sci-Plex 2/3/4 (A549, MCF7, K562 vehicle), Lotfollahi A549, Sunshine Calu-3 (NTC), Papalexi / Wessels THP-1,
  Datlinger Jurkat, K562 sets (Adamson, Dixit, Gasperini, Norman, Schraivogel, Xie), Xu HEK293, Myllia A549 if readable;
- donor control profiles already on disk (`data/donors/*/b_ctrl.npy`): k562_gwps, k562_essential, viperturb_k562,
  nadig_jurkat, orion_hct116, orion_hek293t, rpe1;
- non-CCLE profiles (abstention cases): Arc H1 2025 (Flex), KOLF2.1J iPSC, Joung hESC, Tian iPSC / neurons, Frangieh
  patient melanoma, Zhao patient glioblastoma, Shifrut primary T cells;
- Flex rows: A/B/C (Jurkat / HeLa / CAL33 by 041 + 356), Arc H1 2025.
The matcher is 041's code imported from `def_identity.py` unchanged, with the same three anchors, so the method weights
are exactly the Oct 22 weights. **Known answer:** the harness run with only A/B/C + anchors must reproduce
`def_identity.py`'s method_quality, consensus calls and unanimity on `data/challenge` exactly, or nothing is read.

## Primary quantities
1. **FP_reg**: over queries whose true line is NOT registry line R, the fraction that receive a unanimous call of R (the
   event that switches the wrong line's rows on), pooled over R. Reported overall, for "sister" queries (same ONCOTREE
   lineage as some registry line), and per registry line.
2. **TP_reg**: over queries whose true line IS a registry line present in the CCLE matrix, the fraction with a unanimous
   correct call (the event that switches the right rows on).
3. Top-1 accuracy, P(unanimous), P(correct | unanimous), and the consensus margin (top-1 − top-2 Borda score) for correct
   vs wrong calls, overall and by chemistry and by cell count.
4. Abstention: what the non-CCLE profiles are called, and whether any call is unanimous.

## Decision rule (fixed now)
- **FP_reg = 0 on this set, and no unanimous wrong call among sister queries** → the current rule stands; report the
  margins for the record.
- **Otherwise** → propose an additional gate (a minimum consensus margin and/or a minimum per-method Pearson gap, fitted
  to give FP_reg = 0 here) and report its TP_reg cost; unapplied, the coordinator decides. A gate that also rejects any of
  A/B/C's (Flex) calls is not proposed.
- Chemistry caveat (peer's note, 2026-10-01): MIX-seq and most sets are 3' chemistry; D/E/F are Flex. A Flex-like read is
  added: (a) A/B/C's own margins placed in the 3' margin distribution; (b) the 3' queries re-scored after adding a per-gene
  Flex − 3' log2 bias estimated from same-line pairs (Arc Jurkat A vs Nadig / Datlinger Jurkat; Arc H1 2025 vs Joung
  hESC), pairs never used as validation queries for the bias they fit. (b) is reported beside the plain read; the decision
  uses the worse of the two.
