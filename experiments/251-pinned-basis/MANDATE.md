# 251 — The PINNED deflation basis in the builder: October ships A/B/C's `U_5`/`V_5`, never a day-of re-fit. Inert by default.

You are a build subagent. Work in `experiments/251-pinned-basis/`. **The only file you may edit outside it
is `scripts/make_submission_pm_cd4.py`** (changed today by 222, 247, 249, 248 — read their RESULTS; do not undo
them; 248's byte copy `experiments/248-second-residual/work/make_submission_pm_cd4.248.py` is git HEAD's
builder). You may write a pinned-basis ARTIFACT under `data/donors/cd4_deflation_basis/` (rule 11: a path and
a sha) and build known-answer packages into `experiments/251-pinned-basis/work/` only. **NEVER `vcc submit`.**
`uv run python`. Heavy stages (≥ 15 GiB) strictly serial, alone, gated on `memory_pressure` ≥ 40% free AND
`vm_stat` free+inactive ≥ 30 GiB; record peak RSS. Commit the builder with its known-answer hashes.

## Why
250 (`experiments/250-deflation-october/RESULTS.md`, read all of §4 and §6) found that 246/249's CD4
deflation transfers to disjoint 300-target panels in its `pds` mechanism, but a **day-of re-fit** of the
projector on the new panel's own rows removes panel-specific signal too: its cosine cost is 4–14× A/B/C's and
its price is ≤ 0 on both low-coverage draws. **The PINNED A/B/C basis (249's `U_5`, sha `f06e0262…`) beats
the re-fit on 4/4 disjoint panels** (`pds` +0.0153 / +0.0037 / +0.0073 / +0.0107 at −0.1 … −2.0% cosine,
price +0.0061 / +0.0009 / +0.0014 / +0.0030, never negative), and any foreign basis behaves the same: a basis
fitted on other rows removes the shared T-cell hub and nothing else. 250's runbook row: *"ship the deflation on
D/E/F as the PINNED A/B/C basis, k = 5 — never a day-of re-fit; pinning needs a builder path distinct from
`VCC_CD4_DEFLATE_K` (rule 10), a 249-style byte-identity known answer on A/B/C and a 199-bundle rehearsal, all
before Oct 15."* This mandate is that path. Whether it ships is decided by the 2026-09-05 board point on
`cd4svd5`; the capability must exist either way.

## READ FIRST
`docs/READ-FIRST.md`; 250 RESULTS §4.3–§4.4, §6 and its `cross250.py`/`basis250.py` (how the pinned basis is
applied: gene-space `V_5` from 249's staged A/B/C residual, applied AFTER centring on the target panel's own
carried rows, support re-zeroed); 249 RESULTS §1 and `deflate_cd4_residual` in the builder; 246 `lib246.py`;
247 §3 (carrying-row sizing off A/B/C — the pinned deflation composes with it); 233's per-panel caches
(`experiments/233-october-retention/work/{folds,ops,cd4}`) and 238's `ef-bundle`; `AGENTS.md` §11.1 rule 10, §14.

## The change
1. **Export the basis as a pinned artifact:** from 249's staged A/B/C CD4 residual (the same object whose
   deflation audit recorded `basis_u_k_sha256 f06e0262…`), write `data/donors/cd4_deflation_basis/V5_abc.npz`
   holding the gene-space orthonormal basis `V_5` (18,533 × 5, or the support-restricted equivalent 250 used),
   the row-space `U_5`, the centring convention, and provenance (source residual sha, k, date). Record its
   sha256 in the builder's pin table (the `UNION_SOURCE_REGISTRY`/`WHITENER_REGISTRY` pattern) — never a
   description.
2. **New flag `VCC_CD4_DEFLATE_BASIS=<registry name>`**, mutually exclusive with `VCC_CD4_DEFLATE_K` (refuse
   both; refuse an unregistered name; refuse if the file's sha does not match the pin). When set: centre the
   CD4 residual over ITS carried rows, project out the pinned `V_5` in gene space, re-zero outside the support,
   then size as usual (247's rule off A/B/C). Meta: `effects.cd4.deflation` with `mode: pinned`, the basis name
   and sha, `energy_removed_frac`, and **the fraction of the day's centred CD4 energy captured by the pinned
   basis** (250's diagnostic: 0.10–0.14 on disjoint panels, 0.198 on A/B/C).
3. Wire nothing into `def_pipeline.py`'s defaults (the runbook decision is the lead's, after the board point);
   but make `def_postbuild.py` accept and echo the new keys if 247's shipped-arms check would otherwise refuse
   them — report exactly what you had to touch there, if anything.

## Known-answer rows (abort on failure)
1. Flags unset: staging `effects_array_sha256` `abeb0749…` reproduces; the A/B/C package rebuilds byte-for-byte
   (`8df4dadd…`, 12,648,807,800 B) — 248 §2.2 did this an hour ago through HEAD's builder; if its verified
   copy is gone, rebuild once.
2. **On A/B/C the pinned basis IS the re-fit basis**, so `VCC_CD4_DEFLATE_BASIS=V5_abc` must reproduce
   `VCC_CD4_DEFLATE_K=5`'s staged effects bit-for-bit (249's staged `on` effects sha — read it from
   `experiments/249-cd4-deflate/work/stage/`) and the same `effects.cd4.deflation` energy numbers.
3. **On 233's `real` and `efr` panels** through the builder's staging path (233's per-panel operator and CD4
   re-extraction, 247's carrying-row sizing), the pinned path reproduces 250's pinned rows through 162's
   instrument — `pds` +0.0153 (`real`) / +0.0107 (`efr`) at −0.1 … −2.0% cosine — to the digit 250 reports,
   and the captured-energy diagnostic reads 0.10–0.14.
4. **The rehearsal row:** run `scripts/def_pipeline.py` on 238's `ef-bundle` with `VCC_CD4_DEFLATE_BASIS=V5_abc`
   exported (nothing else), stages through `postbuild` and `vcc prep --dry-run --contexts D,E,F`; the meta
   carries the pinned deflation block, `def_postbuild` passes, and the diagnostic prints. (244 measured 29.7
   min for the full sequence; reuse 244's fitted operator for this bundle if the driver's skip branch accepts
   it, and say so.)

## Deliverable
`RESULTS.md`: the artifact's path and sha, the diff summary, the four known-answer rows with verbatim lines,
the captured-energy diagnostic on every panel you touched, the commit hash, and one line: **PINNED PATH READY
/ NOT**. Propose the `DEF-RUNBOOK` §3 row text and the census trigger (250 §6: CD4 carried rows < 250/300 →
re-fit reads ≤ 0; the pinned path stays ≥ 0). Do not edit the runbook. Then stop.
