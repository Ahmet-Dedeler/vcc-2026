# 248 — A SECOND union residual in the builder (229's "CUM" construction), inert by default, then price and build the DLD-1 rider

You are a build subagent. Work in `experiments/248-second-residual/`. **You may edit exactly one file
outside it: `scripts/make_submission_pm_cd4.py`** (247 has finished with it — read its RESULTS and its
commits first; do not undo them). You are authorised to build packages into `submissions/` under the
`248-` prefix for THIS mandate only, after every gate below. **NEVER `vcc submit`.** `uv run python`.
Builds are ~20 GiB and `vcc prep` ~25 GiB: strictly serial, never alongside another build/prep
(`ps aux | grep -E "make_submission|vcc prep|def_pipeline"`), and gate on `memory_pressure` ≥ 40% free
AND `vm_stat` free+inactive ≥ 30 GiB. Record peak RSS. Commit the builder with its known-answer hashes.

## Why
229 (`experiments/229-stack4-build/RESULTS.md` §6) found the best arm on the table — `dld1_lowmoi` as a
SECOND residual on top of the shipped H1 union (`CUM1` +0.00228, `CUM4` +0.00367 rate-carded; 225 §7) —
**is not expressible in the single-union builder**: forcing DLD-1 into the H1 union (`ARM3`) reads
+0.00148, and the six-source union (`ARM6`) is what tonight's slot tests. A rider on top of the shipped
stack needs a builder that carries a second residual sized on ITS OWN carrying rows, merged into the one
add and the one preserve (218 §7.2's `merged()` with three terms). 208's rider rule: price against the
CURRENT incumbent, never the one it was measured on.

## READ FIRST
`docs/READ-FIRST.md`; `experiments/247-oct22-guards/RESULTS.md`; 229 §5–§6 and `build229.sh`; 222
(`build222.sh`, `stage222.py`, `repro222.py`, `profile222.py`, `metadiff222.py`, `probe222.py` — reuse
them); 218 §7 (`compose218.py` `merged()`); 225 §7 (the CUM rows, and INBOX-CONCENTRATION: DLD-1's Δ`pds`
is concentrated in its top-3 targets); 243 §4.4 (the class rule: a first source into empty rows k_cos
0.94–1.00, k_pds 0.70) and §5.2 (measure the board cosine off the package, never the floor); 216
`six216.py` (the exact `nmae` identity — it has mis-signed the board twice; report it as a magnitude
bound, and separately); `AGENTS.md` §11.1 rule 10, §14.

## A — the builder change
Add `VCC_UNION2_SOURCES` (registry names) and `VCC_UNION2_G` (dose), both unset by default. When set:
CD4 (its own rule), union 1 (`VCC_UNION_SOURCES`/`VCC_H1_UNION_G`, carrying-row rule) and union 2
(carrying-row rule) are each solved against the same base denominator, merged into ONE residual, added in
ONE step, followed by production's ONE global preserve. The meta records an `effects.union2` block with
the same fields as `effects.h1_union` (sources + four-file shas, carrying rows, `g_per_unit`, solved
fraction, `shipped_median_rule_would_solve`).
**Known-answer rows (abort on failure):**
1. Flags unset → the staging effects reproduce 222's `stage-on028.json` `effects_array_sha256`
   **`abeb0749…`** exactly, and — memory permitting — the A/B/C package rebuilds byte-for-byte
   (`.h5ad` `8df4dadd…`, 12,648,807,800 B). If 247's flag is on by default for A/B/C, use 247's
   recorded reference instead and say so.
2. Order invariance: `VCC_UNION2_SOURCES` = the H1 pair with union 1 unset reproduces the same staged
   effects as the shipped configuration (union 1 = H1 pair, union 2 unset), bit-for-bit.
3. `union2 = dld1_lowmoi @ g2 0.35` on top of the shipped stack reproduces 229/225's `CUM1` PLOLO row
   through 162's instrument on the builder's own staged effects (rel `cos_agg` +1.405%, Δ`pds` +0.00256
   at 0.35; three cancer folds) — 222's `repro222.py` pattern.

## B — price the rider, both instruments, all six members
`CUM1` (dld1 @ 0.35) and `CUM4` (225's four as union 2, if the registry carries them) on top of the shipped
stack, three cancer folds with `kolf2` separately: `cos_agg`, macro, `pds` (must be ≥ 0), `q_w`;
falsifiers: target-shuffle at every dose, row-matched CD4 control, non-carrying rows, **leave-top-3-
targets-out on cosine and `pds` separately** (225 found DLD-1's `pds` gain concentrated). Then two prices:
(i) 243's class rule with the arm's own delivered amplitude (stream the staged package's `A_c`/`B_c` as
`totals243.py` does) for the direction legs plus the DE-four propagation; (ii) 216's exact identity
for `nmae`, reported separately with its two-for-two mis-sign record. Give the composite under each.

## C — build, conditionally
If (i) clears **≥ 1.5× the shape-only gate (+0.0019)** with `pds` ≥ 0 on 3/3 cancer folds and the
leave-top-3-out retention ≥ 50% on `pds`, build `submissions/248-cum1-dld1-g035.{h5ad,vcc}` through
222's full gate set (rule-10 metadiff against `222-h1union-g028.h5ad.meta.json` with exactly the
pre-declared differences; profile split by carrying row, scale-corrected; non-vacuity; challenge `n_pred`
±2%; `vcc prep --dry-run`), record sha256 and bytes, and stop. **The lead decides whether it is
submitted, after tonight's board reading of ARM6 (237).** If it does not clear, do not build; say why.

## Deliverable
`RESULTS.md`: the three known-answer rows, the pricing table under both instruments with every
falsifier's number, the package (or the reason there is none), the commit hash of the builder change
with its inertness proof, and one line: **RIDER BUILDABLE / NOT**. Then stop.
