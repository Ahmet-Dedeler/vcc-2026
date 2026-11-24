# 259 — The union WITHOUT the source that carries its `nmae` cost, on top of `cd4svd5`, in the two constructions 258 did not test. Measurement only, 2–3 hours. DO NOT BUILD, DO NOT SUBMIT.

You are a research subagent. Work in `experiments/259-union-without-igvf/`. **Edit nothing outside it.** No build, no
`vcc`, no `cell-eval2`, nothing under `submissions/`. `uv run python`. One PLOLO fold at a time, peak RSS under ~10 GB;
a sibling subagent (260) may run a light inventory alongside. **Time box: stop and write up after ~3 hours of wall
clock even if the grid is incomplete; a partial table with its retrodictions is a valid result.**

## Why
258 (`experiments/258-union-dose-on-cd4svd5/RESULTS.md`) priced the six-source union on top of `cd4svd5` and found that
its `nmae` cost is carried by **`igvf_cardiomyocyte`** (the H1 pair + igvf alone charges +0.0144 raw `nmae` at g = 0.28
for no `pds`), while the pair + `dld1_lowmoi` is `nmae`-cheap (+0.0014 raw, e +0.20) and the only subset with `pds`
positive on 3/3 folds (+0.00119, 0.8 LOO-σ, leave-top-3-out −66 %). **But 258's `stack5` dropped Myllia, not IGVF: the
four-source union without the blamed source — H1 pair + `dld1_lowmoi` + `jiang_5line_noK562` + `myllia_a549` — was
never measured on this base.** 225 (`experiments/225-source-sweep/RESULTS.md` §"LOO −igvf") read the analogous
leave-IGVF-out stack at +0.00336 = 1.77× the gate — on the OLD base (`cd4f034`, no H1 pair in the incumbent, no
deflation) and the old rate card; that number is history, not today's prediction. And 258's dose sweep lowered ONE
shared dose for the whole union, which under-doses the H1 pair the incumbent already carries at 0.28 (cosine −0.8 …
−4.8 % below g = 0.28) — so "the union fails at lower doses" never isolated whether the added sources are useful at a
dose of their own.

## READ FIRST
`docs/READ-FIRST.md` (v2.7), 258's `RESULTS.md` §0–§6 and `price258.py` (REUSE its instrument: the `apply_frozen_cd4_residual`
route on the `cd4svd5` base, the two functionals, the production-frame composition, its fold caches under
`experiments/258-union-dose-on-cd4svd5/work/plolo/` — byte-copy them, do not recompute), 257's `RESULTS.md` §0 item 9
and §1.3 (constants vs the captured entries), `experiments/248-second-residual/RESULTS.md` (the builder's SECOND union
residual: `VCC_UNION2_SOURCES` + `VCC_UNION2_G`, sized by `solve_carrying_row_fraction`, its mixing order and why its
DLD-1 rider was not cleared: 1.14× gate, leave-top-3-out 25 %), 225 §LOO, `AGENTS.md` §8.1b, §11.1 rule 10, §14.

## The two constructions (three cancer folds; `kolf2` apart; per-fold rows; LOO σ)
A. **Single union replacing the pair**, as 229/258 built it: `VCC_UNION_SOURCES=gse295214_h1,h1_2025,dld1_lowmoi,jiang_5line_noK562,myllia_a549`
   at g ∈ {0.20, 0.28, 0.35}, on `cd4svd5` (k = 5, f 0.34).
B. **The incumbent's pair kept at 0.28, plus a SECOND residual** `VCC_UNION2_SOURCES=dld1_lowmoi,jiang_5line_noK562,myllia_a549`
   at g₂ ∈ {0.15, 0.28, 0.35}, through 248's `union2` path in the production mixing order (one merged add, one
   preserve — verify against 248's own known-answer rows; do not invent the order).
Plus the controls that make the reading honest: (i) the same B with IGVF added back (`dld1,jiang,myllia,igvf`) at the
best g₂ — the direct test of the attribution; (ii) B with `dld1_lowmoi` alone at the best g₂ (248's rider, re-priced);
(iii) for the best cell of A and of B: leave-top-3-out `pds` retention, the per-target `pds` share of the top 3, and
the energy-matched random-row control if the construction admits one cheaply (else say so).

## Known-answer rows before any candidate number is read
R1 258's `cd4svd5`, `stack6 g=0.35`, `stack3 g=0.35` and `pair+igvf g=0.28` rows reproduce from YOUR runner at
0.0e+00 (same caches, same functionals). R2 the recipe retrodicts BR249 (`cd4svd5` vs `h1union` +0.0043 ± 0.0008) and
BR237's ARM6 by member with the factor set 258 §0 found (state it). R3 248's second-residual known-answer row
reproduces through your union2 path.

## The price — COHERENT PESSIMISTIC scenarios, not one "conservative" coefficient set
257 §0 item 9 and 258 §0 showed a single "conservative" coefficient set can flatter a candidate (a 0.85 `pds` factor
SHRINKS a `pds` LOSS; a 3.76 cosine factor HELPS a cosine GAIN). Price every cell three ways and report all three:
- **PESSIMISTIC** = for each leg, the end of its measured band that LOWERS the total, chosen per the sign of that
  leg for that cell: `pds` k ∈ {0.70, 2.48} (243's class value / ARM6-realised A-B-C mean), cosine k_cos ∈ {0.94,
  3.76} (new-source / re-shape), DE-three rate ∈ {0.0425, 0.24}, `nmae` identity × {1.0, 1.5} on a loss and ×0 credit
  on a gain, amplitude leg always charged.
- **POINT** = ARM6-realised for the union legs (k_pds 2.48, k_cos 0.96, DE-three 0.0425, identity ×1.0).
- **OPTIMISTIC** = the other ends. A cell is CLEARED only if PESSIMISTIC ≥ +0.0019 over `cd4svd5` AND leave-top-3-out
  retention ≥ 50 % AND `pds` positive on 3/3 folds at ≥ 2 LOO-σ. Beating ARM6 or `h1union` is irrelevant: the
  comparison is `cd4svd5`.

## Deliverable
`RESULTS.md`: R1–R3 first; the A and B tables; the IGVF-back and DLD-1-alone controls; the best cell in full (six
members at the three scenarios, retention, amplitude, the exact env the builder would need — every flag exists:
`VCC_UNION_SOURCES` / `VCC_H1_UNION_G` / `VCC_UNION2_SOURCES` / `VCC_UNION2_G` on top of `VCC_CD4_DEFLATE_K=5`); one
line **CLEARS / DOES NOT CLEAR**, and, if it clears, `PREREGISTRATION-DRAFT.md` in 249's form with five falsifiers. A
clean negative — "the attribution to IGVF is right but removing it does not buy a cleared arm" — is a valid,
useful result. No model name, no description, no submission. Then stop.
