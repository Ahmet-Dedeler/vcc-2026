# 643 — Does a construct's seed off-target signature travel with the construct into another cell line? (preregistered)

Written 2026-10-02 before any K562 expression value is read. Session [537757]. Analysis only.

## Why now
Arc's 2026-10-01 post: "We used the same sgRNA sequences in every cell line". Arc's 2025 H1 library reused Replogle 2022's
K562 day-8 (GWPS) dual-guide pairs for 297/317 targeting constructs (301). F16 (307/309) found, in H1, that genes whose
promoters match a construct's protospacer seed are repressed in that construct's own cells and NOT in the cells of the
partner construct of the same target (own +0.713 vs partner +0.043, construct-level z). If that signature is a property of
the protospacer and not of H1, then (a) K562 GWPS rows carry the construct-specific part of Arc's truth wherever Arc reuses
Replogle's pairs, and (b) construct effects measured anywhere are shared by all six 2026 lines.

## Test
- Units: the 34 constructs of 309's 17 two-construct targets (`309/work/enrichment309.csv`), each with its seed-rule
  candidate genes (sequence-only selection, no expression used).
- Data: K562 GWPS cells (scPerturb `ReplogleWeissman2022_K562_gwps.h5ad`, `guide_id` = the dual-guide pair), the same
  protospacer pairs mapped through `301/work/vcc2025_spacers.csv` (sgID set equality). Controls: non-targeting cells.
- Statistic per construct c with partner c' and candidates S_c: own_c = mean over S_c of Δ(c) and partner_c = mean over S_c
  of Δ(c'), Δ = mean log1p(CP50k) in the construct's cells minus controls. D_c = own_c − partner_c.
- PRIMARY: the median of D_c over constructs present with both arms in K562 (and candidates expressed, control mean
  CP50k ≥ 1). A negative median with a paired Wilcoxon p < 0.05 (one-sided, own < partner) = REPLICATES in K562.
- Secondary: the per-construct correlation of D_c(K562) with 309's H1 construct contrast (delta_own − partner, same sign
  convention); the on-target knockdown of each construct in K562 vs H1.
- Null/known answer: the same statistic on 10 random gene sets per construct, matched on control expression decile
  (expect median ≈ 0). If the matched-null median is itself below −0.01, the test is UNRESOLVABLE.

## What it would change
- REPLICATES: construct effects are line-transferable. K562 GWPS rows are construct-matched to Arc's truth if Arc reused
  Replogle's pairs in 2026 (the 2025 precedent). This would explain h1-test's K562 preference (097 a3 argmax share 0.60;
  175 sign accuracy K562 0.567 vs Orion 0.53-0.54) as a library property, and a donor-weight read on the A/B/C board would
  then be a property of the shared library that carries to D/E/F, unlike other board-tuned knobs. No recipe change by
  itself: 178 priced K562-heavy weights at -0.0012 after the operator on H1.
- DOES NOT REPLICATE: F16's signature is H1-specific (chromatin or expression context). The "same-atlas relative"
  reading of 097 a3 stands.
