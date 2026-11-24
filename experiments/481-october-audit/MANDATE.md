# 481 — October transfer audit: which A/B/C-tuned settings change meaning on a new panel?

You are a research teammate (Claude Opus 5.5, max effort) in the Virtual Cell Challenge 2026 repo. Read `AGENTS.md`, then
`docs/FINDINGS.md` §0–§2, `docs/DEF-RUNBOOK.md` and `docs/REVIEW-2026-09-26-second-opinion.md` §2 before anything else.
You own `experiments/481-october-audit/` and write nothing outside it. Budget: up to ~9 hours. Write interim findings
into `RESULTS.md` as you go (sessions die without warning); finish with a complete `RESULTS.md` and a short final
message.

## Why this matters

Only the final entry on D/E/F (three new cancer lines, a new 300-target panel, released Oct 22, no leaderboard) counts.
Every knob in the pipeline was tuned or validated on A/B/C. A setting defined relative to a panel's own distribution
(a quantile, a rank, a count, a coverage fraction) can silently become a different rule on a new panel while every
check still passes. A peer review found one on 2026-09-26: `VCC_DE_PRIOR_GATE=q33` cuts at `np.quantile` of an integer
strength; on A/B/C the cut is 2 (weak = strength {0,1}), on the efr synthetic panel it is 1 (weak = {0} only), so the
switch the board validates is not the switch October would run (fixed as 472's absolute `le1`). Find the rest of this
class before Oct 22.

## Scope

The October path: `scripts/def_day.sh`, `scripts/def_pipeline.py`, `scripts/make_submission_pm_cd4.py` (the builder;
READ ONLY, it is pinned and patched by running chains), `scripts/def_sameline.py`, `scripts/def_lodo_key.py`,
`experiments/448-def-rehearsal-q175/relevel_iter.py` (the re-level secant loop), and the preflight/census/audit/gate
scripts they call. Every `VCC_*` knob and every hard-coded constant on that path, including (not only): sharpening q,
H1/CD4 stack doses and `dayof438.py`'s day-of re-pricing, the re-level target ρ and bounds, cliff checks, the census
floors, identity unanimity, same-line `auto` weights and bank classes, programme coverage rules, the LODO key clip, the
MDE `c,cap,min_count`, the zero-repair floor, ballast K/KMAX/LB, q80 protection, call-set size, CD4 sizing, the prior
gate, `pertscale` clips, and anything read from A/B/C-specific artifacts.

Also look for fail-open paths: a check that passes vacuously on an unusual panel, a default that silently applies when
an artifact is missing, a context-count or context-label assumption.

## Method

For each item: how it is defined; whether it is absolute or relative to a panel/context distribution; its realised value
on the relabelled A/B/C bundle vs the synthetic October panels; whether its meaning shifts; the expected score impact if
it does; a proposed fix with a known-answer test (e.g. "on A/B/C the fix selects the same set, so it rebuilds 446b /
458g byte for byte").

Materials:
- Bundles: `experiments/180-def-readiness/work/def-bundle` (A/B/C relabelled D/E/F: A→D, B→E, C→F, verified by
  content), `experiments/238-essential-free-panel/work/ef-bundle-r` (efr, realistic) and `ef-bundle` (efp,
  pessimistic).
- Day-of runs with every artifact and meta: `experiments/453-abc-dress-rehearsal/work/run` (A/B/C relabelled),
  `experiments/465-both-switches-efr/work/run` (efr, both switches), `experiments/452-def-day-rehearsal/work/efr*`,
  `experiments/448-def-rehearsal-q175/work/efr`, `experiments/469-gated-october-path/work/run` (gated, A/B/C).
- Stack dumps (`dump_*.npz`: bulk, carriers, perts, genes) sit in those run dirs. The reviewer reproduced the builder's
  gate strength from `carriers` plus each context's control gate (CPM > 10 and b_ctrl > 0) to within ~5 %.
- Known answer: your method must rediscover the prior-gate shift (A/B/C threshold 2, efr threshold 1) before you trust
  it on anything else.
- There is no efp stack dump. If an item needs one, say so rather than building it.

## Rules

- Analysis only: no builds, no `def_day.sh` runs, no `vcc`, no `cell-eval2`, no git, no edits to `scripts/`, `src/` or
  any production file. Put proposed patches as separate files in your experiment directory. Peak RSS ≤ 10 GB; check
  `memory_pressure` before any load over 4 GB (the coordinator runs ~20 GB builds).
- Never touch the coordinator's running chains (470b, 471, 472; 473 tomorrow), `submissions/`, or any `.pid` file.
- Rank findings by expected October impact; separate confirmed shifts (measured) from risks (reasoned).
- Other teammates: 480 studies the nmae/magnitude question; 482 (Grok) hunts literature and new datasets. Do not
  duplicate them.
