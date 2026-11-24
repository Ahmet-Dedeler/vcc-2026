# 445 — DLD-1 in the bulk under ballast: merged into the H1 union it dilutes; as its own residual it adds ~+0.8 % cos. SMALL (2026-09-24, coordinator)

## Question
229 killed every builder-expressible DLD-1 construction on nmae (welded channels). Ballast (410b) unwelded the bulk, and
418/420 re-priced the composition knobs but never an extra residual source. Is 223's one live source (DLD-1 low-MOI
CRISPRi, 59–66 of each fold's targets) worth adding to today's bulk (stack|K3 + q 1.5)?

## Instrument
`comp445.py` on 420's PLOLO harness (realised KMAX 16, β-emulated pds). **Known answer: the H1-pair union at g 0.56
reproduces comp420's `stack|K3+sh1.5` with |d| = 0 on all four folds.** Two constructions: (a) DLD-1 merged into the
single H1 union (builder-expressible today: `VCC_UNION_SOURCES=gse295214_h1,h1_2025,dld1_lowmoi`); (b) DLD-1 as its own
post-whiten residual, carrying-row sized at its own dose (223's construction; needs a builder slot). Controls: DLD-1's
rows permuted over targets.

| arm (vs stack; rel cos / Δ pds@.5) | k562_gwps | orion_hct116 | orion_hek293t | kolf2 |
|---|---|---|---|---|
| (a) merged, g 0.42 | −2.30 % / +0.0031 | −3.44 % / −0.0009 | −1.47 % / +0.0016 | −12.8 % / −0.0088 |
| (a) merged, g 0.56 | −1.59 % / +0.0026 | −2.44 % / −0.0028 | −1.23 % / −0.0005 | −6.3 % / −0.0079 |
| (a) merged, g 0.70 | −1.82 % / +0.0019 | −2.29 % / −0.0052 | −1.68 % / −0.0033 | −1.0 % / −0.0073 |
| CTRL merged, permuted | −4.92 % / −0.0028 | −5.00 % / −0.0110 | −4.07 % / −0.0040 | −7.3 % / −0.0042 |
| **(b) own residual, g 0.2** | **+0.87 % / +0.0016** | **+0.83 % / +0.0003** | **+0.72 % / +0.0011** | +0.11 % / −0.0014 |
| (b) own residual, g 0.35 | +0.85 % / +0.0027 | +0.80 % / −0.0003 | +0.70 % / +0.0011 | −0.32 % / −0.0027 |
| (b) own residual, g 0.5 | +0.37 % / +0.0034 | +0.13 % / −0.0013 | +0.13 % / +0.0002 | −1.24 % / −0.0040 |
| CTRL own residual, permuted, g 0.35 | −1.30 % / +0.0008 | −1.07 % / −0.0052 | −0.87 % / −0.0018 | −0.77 % / −0.0010 |

## Verdict
- **(a) is closed:** a per-gene mean with the H1 pair dilutes the pair's residual on the ~20 shared targets; worse than
  the stack at every dose on every fold's cosine.
- **(b) is real and small:** DLD-1 carries target information (it beats its permuted control by ~2 % cos on every cancer
  fold), worth about +0.8 % relative cos on today's bulk at g 0.2 (223 measured +1.4 % on the older, weaker bulk), with
  pds slightly positive on cancer folds and slightly negative on kolf2. At composition-step transfer (k ≈ 1–1.3) that is
  **~+0.001–0.0015 overall**. It needs a builder slot for a second residual (a generic `VCC_EXTRA_RESIDUAL=<source>:<g>`)
  and a new stack dump. Queued as a low-priority rider; its October argument is breadth (DLD-1 covers ~73 of a disjoint
  panel on 225's pessimistic draw, one of the few sources that gains coverage).
