# 524 — the October recipe dress rehearsal on the board (preregistration, 2026-09-29 01:10Z)

Written before any 09-29 entry exists (457dr was scoring, unpublished), before the reads that decide the pending October
switches (09-30, 10-02, 10-03, 10-05), and before the arm is built. Slot: 10-06 slot 1 (`chain1006.sh`).

## The arm
**524dr** = `def_day.sh` on 180's relabelled A/B/C bundle with EXACTLY `scripts/october_recipe.py command --allow-pending`
as resolved at ≥ 10-05 02:00Z (`recipe_dr524.sh` saves it as `work/recipe524.json` with a copy of
`docs/october-decisions.json`), 453's day operator cloned (66610be7), relabelled to A/B/C, prepped, gated by `gate524.py`
(config diff vs 453dr ⊆ the keys the resolved switches explain; every ON switch visible in the meta; postbuild, cliff,
zero audit, novelty). Reference: **453dr = 0.2663613** (entry `q3VuTImh4FRGKXH2waCp`), the plain October path.

Switches as of this writing: DECIDED on — `lodo` 1, `de_mag` dn_k250_zs, q 1.75, w 0.5 + `noprog` (vs 453dr's w auto + `pt`:
no effect on relabelled A/B/C, where the matched Jurkat bank covers 0 targets as rows), `cd4_stack` x1. PENDING —
`programme_source`/`prior_gate` (09-30), `sharpen_pt` (10-02), `lodo_hi_pct` (10-03), `mde_sign` (10-05). A switch still
PENDING at resolve time runs OFF and is listed in `recipe524.json`.

## Why
The final entry runs this path on D/E/F. 453dr validated the plain path (bulk −0.0031 vs production, DE pinned); every
switch decided since was read on PRODUCTION A/B/C arms, never on the October path. This read checks that the day script
carries the adopted switches, and that they add up there as they did in production. A failure found on 10-06 leaves 16 days
to fix it; one found on Oct 22 leaves none.

## Prediction (fixed rule; the numbers come from reads that all precede the 524dr entry)
Σ = the sum, over the switches ON in `recipe524.json`, of each switch's own board read on A/B/C (the entry that decided it
minus that entry's reference), per member and overall:
- `lodo` = 446bj2 − 446b: **+0.0022338** overall (mse −0.0145 raw, pds +0.0002, DE bit-identical);
- `de_mag` = 446bj2z − 446bj2: **+0.0005494** overall (nmae −0.002179 raw, the other five bit-identical);
- `programme_source` (+ `prior_gate`) = the adopted arm − 446b on 09-30 (458d or 458g per read468.py);
- `sharpen_pt` = the adopted 479 arm − 446bj2 (10-02); `lodo_hi_pct` = 495 − 446bj2 (10-03); `mde_sign` = 496a − 446bj2
  (10-05).
With only the two DECIDED switches, Σ = +0.0027832.

## Bins (overall Δ = 524dr − 453dr)
- **|Δ − Σ| ≤ 0.0020: PASS.** The October path carries the recipe additively. Nothing changes.
- **Δ < Σ − 0.0020: LAG.** Before any other October work, find the lagging member: bulk members (pds/mse) point at the
  LODO key or sharpening on the October path's re-level (compare ρ, factors and the dump with 460dr/497dr); DE members point
  at a DE switch that did not reach the emitter (meta keys, the def_day log). A switch whose October-path contribution reads
  ≤ 0 comes out of the recipe unless the defect is found and fixed with a rehearsal.
- **Δ > Σ + 0.0020: super-additive.** Record it; no action.
- **Falsifiers (integration defects, diagnose before reading Δ):** a DE member moving by more than its switches' Σ ± (nmae
  0.002, fid 0.003, reach 0.012, jac 0.001 raw); pds or mse moving while every ON switch is DE-only; a member moving
  that no ON switch touches on A/B/C by more than the tolerances above.

The tolerance 0.0020 is ~2× the bulk spread the October re-level showed against production (453dr −0.0031 in a
−0.0045 … −0.0005 band); DE switches that share their emission draws (468) are deterministic across paths.
