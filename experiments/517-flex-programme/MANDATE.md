# 517 — Target-independent DE sign information: the lab/Flex programme, ready for the 09-30 read (worker mandate)

Written 2026-09-28 13:15Z by the coordinator. You are a research worker (Claude Opus 5.5, max effort). Read `AGENTS.md`,
then `docs/FINDINGS.md` §0–§2 (especially "Programme priors are lab signatures" and the 462–470 notes), then the
experiments named below. Your directory: `experiments/517-flex-programme/` (write ONLY there).

## Objective
On 2026-09-30 (~00:30Z) the board reads 458d (Arc's H1 programme from all public H1 cells as every context's DE block
prior, vs 446b) and 458g (458d + the design-strength gate q33, vs 458d); `experiments/468-prior-gate/read468.py` applies
the preregistered bins. That read decides whether a same-lab, same-chemistry (10x Flex) programme carries across cell
lines. Your job is to make the next steps ready BEFORE the read, for every outcome, and to find any other
target-independent sign/magnitude information we are not using.

1. **Positive or partial read → exploitation arms, built on paper now.**
   - The uncalled-gene sign hook (496, `VCC_DE_MDE_SIGN=clean4:0.05`, board read 10-05 as 496a) hard-wires CLEAN4 as its
     source. Specify (and price on the instruments below) the same hook sourced from Arc's programme, and the gated form.
   - A programme-strength dose (458e is the 0.5 blend with CLEAN4, already built) and where the q80 protection should
     sit when the prior is Arc's.
   - The October artifact: Arc's programme computed without any D/E/F rows, on an arbitrary panel's gene axis (the
     day-of path is `DEF_DAY_PROGRAMME_SOURCE`, see `scripts/def_sameline.py --programme-source`); confirm what exists.
2. **Null or harm read → the best alternative.** Which target-independent source is next best for A/B/C and D/E/F
   (CLEAN4 variants, other labs' programmes weighted by measured cross-line agreement, own-control axes (491), a
   chemistry-technical component)? Rank them with numbers.
3. **Flex technical components.** In 10x Flex (probe-based) data, part of a "programme" may be technical: library-size /
   composition shifts between perturbed and NTC cells (CPM compositional effects), probe-level detection, per-gene
   dropout. Measure, in the public H1 2025 Flex data and in the challenge's own control cells (46 NTC guides × 400 cells
   per context), whether such components exist and whether they predict H1's DE signs beyond the biological programme.
   Testh's public descriptions (in board snapshots) mention a "library-size tilt" in their DE emission; find out whether a
   compositional effect is a real, transferable signal.

## Territory (disjoint from the other workers)
- YOURS: per-gene, target-independent information: programmes, lab/chemistry signatures, technical composition effects,
  control-derived axes, the hook's source and dose.
- NOT yours: target-specific donor rows and their combination (worker 516), the bulk channel (518), production code,
  builds, slots (the coordinator).

## Instruments
- H1 is Arc's own line: a programme estimated from H1 cells and scored on h1-test is circular unless you split controls
  and targets properly — follow 462's rule (judge a source with CONTROL halves, never target halves) and 441/462e/f's
  fid-like and gate-block rows. The Broad-16 (one lab, 12 lines, matched targets; `data/external/broad16line`) and the
  Weissman lines (K562, RPE1, HepG2, Jurkat) and Orion (HCT116, HEK293T) give cross-line, within-lab structure.
- Emitter-level pricing on H1: `experiments/496-uncalled-sign-hook/` (run496.py, the hook), `experiments/468-prior-gate/`,
  `experiments/462*`, `experiments/467*`, `experiments/470*`. Reproduce one of their recorded numbers as your known answer
  before reading anything new (e.g. 496's hook nmae −0.00146 on seed 0, or 462's fid-like 0.5545 / CLEAN4 0.5244).
- Analytic DE members: `experiments/513-reach-walk-anatomy/walk513.py` reproduces the production emitter's H1 fid, jac,
  nmae, reach exactly; reuse it to price sign substitutions.

## Prior art — read, do not repeat under a new name
441, 446bja (Jurkat programme: board null), 450h–450l (lab signatures), 458/458d/458e/458g (PREREGISTRATIONs and
RESULTS), 462 (single-control noise), 466, 467 (strength gate), 468 (gate), 470 (the weak third wants CLEAN4), 472, 478
(uncalled-gene delivery through the coarse MDE), 491 (own-control programme: small), 496 (the hook with CLEAN4: weak).

## Deliverables (hard deadline for §1–2: 2026-09-30 00:00Z)
- `RESULTS.md` with known answers, numbers, and for each 09-30 outcome (positive / partial / null / harm) the ranked
  next arms: exact env knobs, the builder/emitter patch each needs (as a patch file in your dir, not applied), the
  expected config diff vs its base, H1 prices with target-bootstrap CIs, and a board-bin draft the coordinator can turn
  into a preregistration.
- §3 (technical components) can follow after the deadline.

## Rules
- Memory: a production build queue runs until ~02:00Z 2026-09-29. Until then keep EVERY process ≤ 4 GB RSS and check
  `vm_stat` free+inactive ≥ 33 GiB before any step expected above 2 GB; after that ≤ 10 GB per process, one heavy step
  at a time. Report peak RSS (`/usr/bin/time -l`).
- No builds, no `vcc`, no `cell-eval2`, no git, no edits outside your directory. Never start long-lived shell commands
  whose arguments name production scripts (the build chains watch process lists).
- Known-answer row in every screen; target-bootstrap CIs; negative results are results. Finish with RESULTS.md and a
  final message summarising the numbers.
