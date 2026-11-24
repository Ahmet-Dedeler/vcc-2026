# 451 — the October same-line mode (plan, 2026-09-25 03:00Z, coordinator)

## Why
450 (H1, analytic): a same-line bank from another lab and chemistry (Pan 2025, H1 hESC, 10x 5′) carries the recipient's
response. Against h1-test truth:
- block-sign programme 0.625 vs CLEAN4 0.581 (oracle 0.653);
- per-target gate signs 0.764 vs our 0.551;
- per-target centred cosine 0.260 vs 0.069.

DEF-RUNBOOK §3 row 220 left the same-line case open ("must be finished before Oct 22 or not attempted"). The six lines
come from different tissues, and blood is taken by Jurkat, so the likely October tissues are colon, liver, lung, breast,
pancreas, kidney. We hold banks for many of their canonical lines.

## Bank registry (CCLE sample prefix → banks, primary first)
| CCLE line | banks (perturbations) | kind |
|---|---|---|
| HCT116 | orion_hct116 (17,766) | genome-wide |
| K562 | k562_gwps (9,786), viperturb_k562 (17,488), k562_essential (2,020) | genome-wide |
| DLD1 | dld1_lowmoi (4,715) | broad |
| JURKAT | nadig_jurkat (2,269) | essential: programme only on essential-free panels |
| HEPG2 | nadig_hepg2 | essential: programme only |
| A549 | myllia_a549 (175), jiang_A549 (218) | targeted |
| HT29 / MCF7 / BXPC3 | jiang_HT29 / jiang_MCF7 / jiang_BXPC3 (218 each) | targeted, cytokine-stimulated |
| CALU3 | sunshine_calu3 (197) | targeted |
| THP1 | wessels_thp1 (185) | targeted |
HEK293T, RPE1 and H1 are not in CCLE; the identity stage cannot name them, so they are not registered.

## Mechanism (two parts, both per context, both only when the identity call is unanimous)
1. **Programme prior (DE).** The context's default block sign is the bank's P_down (fraction of its perturbations with
   delta < 0 over non-zero entries; panel targets excluded; uncovered genes down). This is 396's `psign_A` rule, so on
   A/B/C with A = Jurkat it must reproduce `clean4q80ja` (446bja) byte for byte. That is the known-answer row.
2. **Per-target rows (bulk and DE).** For every panel target the bank covers, the design row and the bulk row become
   `unit(ours) + w·unit(bank row)`, rescaled to our row's own norm (magnitude kept, direction moved). Panel columns are
   untouched. Optionally the covered targets' DE prior is the bank row's own sign. w and the per-target prior are chosen
   from 450b/450c's H1 six-member prices.

## Pieces
- `scripts/def_sameline.py`: identity JSON + registry → `<work>/sameline.npz` (per context: banks, covered targets,
  unit rows on the official axis, sign rows, programme). A driver stage after `identity`.
- Builder knob `VCC_SAMELINE="<npz>:<sha256>:<w>[:pt]"`. Unset must be bit-identical (prove by rebuilding 446b).
- Known answers:
  - unset reproduces 446b;
  - A = nadig_jurkat with w = 0 reproduces 446bja;
  - a rehearsal on a synthetic bundle whose panel includes bank-covered targets exercises the per-target path.
- DEF-RUNBOOK: the stage, the rule ("unanimous identity + registry hit → same-line mode"), and the fallback (no hit:
  nothing changes).

## Validation limits
The board can only test part 1, on context A (446bja): Nadig Jurkat covers none of the essential-free panel. Part 2's
only instrument is H1 (Pan → h1-test, cross-lab, cross-chemistry, 30 covered targets).
