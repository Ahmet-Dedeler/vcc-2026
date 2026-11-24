# 199 — Oct 22 made mechanical: the whole D/E/F transition pre-staged and rehearsed end to end

Date: **2026-09-03** (relaunched 20:48 UTC+3 after the memory pause; the first instance's four
staged scripts and its 16.8 GB CD4 source download were kept and verified). Readiness only. **No
`vcc submit`, no `cell-eval2`, nothing written to `submissions/`, no tuning.** Peak RSS: every
analysis stage ≤ 2.7 GiB; the operator fits 8.1 GiB at the eigensolve (process high-water 12.9 GiB,
the excess being the file-backed basis memmap); the two builder runs 19.80 / 18.82
GiB, which is the three-donor floor 180 §11 documented and which the pause note explicitly deferred
rather than forbade. Both ran strictly alone, after a free-memory check the driver now enforces.

**Deliverables:** `docs/DEF-RUNBOOK.md` (the numbered day), six new scripts under `scripts/`
(`def_preflight.py`, `def_identity.py`, `def_census.py`, `def_postbuild.py`,
`def_pin_panel_artifacts.py`, `def_pipeline.py`), and this record. No production file's behaviour
changed; the builder was edited only transiently by the rehearsal's re-pin and byte-restored (§7 —
with one incident the coordinator must act on).

---

## Verdict in eight lines

1. **The path is rehearsed end to end and it takes 29 min of wall clock** from a new bundle
   to a prepped, falsified artifact: 33 s of checks and re-extraction, 7.3 min operator fit,
   9.2 s effects run, 21 min build, 13.9 s `vcc prep`. Add ~15 min of upload.
   **A same-day rebuild-and-submit is ≈ 1 h 15 min before scoring starts.**
2. **180's rehearsal had a hole that would have lost the round silently.** Two "ours, never moves"
   frozen inputs are functions of the panel: the whitening operator (axis = union − panel; the A/B/C
   file passes every builder check on a new panel and mis-pairs **98.8%** of its coordinates) and the
   CD4 banks (panel-restricted; `cd4_complete_targets` → 0). Both are now re-made on the day
   (7.3 min and 2 s) and re-pinned by a scripted, diffed source edit (`def_pin_panel_artifacts.py`).
3. **Every stage carries a known-answer row and all of them reproduce**: 041's A/B/C identity call
   (8/8 variants, 24 anchor ranks), AGENTS.md §2's control-file numbers, the 269/248/276 census, the
   cloud CD4 banks **byte-for-byte on 36/36 files** from the pinned source, and the production
   operator **`e2f65ecf…` reproduced exactly** by a refit (bitwise-deterministic on this machine).
4. **`panel_support` is predicted from set intersections before any build**, and the ~19 GiB
   effects run confirmed the prediction exactly on the synthetic panel (141 / 167).
5. **The blind decision table is written** (runbook §3): global `VCC_ALPHA 0.982` ships, the
   `knee × (1 − 1.25%)` alternative is recorded with a numeric trigger; gamma 0.65, lam 0, the operator
   by sha, `rc_callset` (board 0.1918638 today), 192's call set, CD4 f 0.1, no mass fix, the a8/a12/a14
   flags, equal donor weights, the seeds, and a rule for identity-selected donors.
6. **Coverage drops have pre-decided responses** (runbook §5), sized on a pessimistic synthetic panel
   (K562 150/300, CD4 174/300, 7 single-carrier targets) and a realistic one (300/300 everywhere).
7. **Rule 10 holds: the A/B/C build is bit-identical under everything 199 added** — 50/50 recorded meta hashes recompute identically, 9/9 production files byte-identical to `ab8aa3e`, 0 production references to the six new scripts, 65/65 builder tests, `scripts/frozen_inputs.json` unchanged (`a654d66b…`).
8. **One incident:** a sweep commit (`04101c5`, 21:16) captured the rehearsal's transient builder
   pin. HEAD's builder therefore pins the SYNTHETIC CD4 banks; the working tree is restored and
   must be committed (§7).

---

## 1. State on arrival and what was already there

The first instance (paused 12:05 UTC by the coordinator for memory contention) had written
`scripts/def_identity.py`, `def_census.py`, `def_preflight.py`, `def_postbuild.py`, the synthetic
bundle builder, a rule-10 checker, and had completed a 16,786,240,107-byte download of the CD4
atlas source (`GWCD4i.DE_stats.h5ad`) into `work/`. That size equals `source_size_bytes` in the
cloud compact's manifest; the file opens as valid HDF5 (33,983 × 10,282, six layers) and hashes to
`c355f535ff32cf7ba1edc49cf9c6039fe84f2c9ebe4d005515cba75790cfbb62` (31 s). It is now the pinned
source for CD4 re-extraction (runbook §1).

---

## 2. Item 2b — the CCLE matcher, pre-staged

`scripts/def_identity.py` is 041's arithmetic verbatim with the ingest replaced by a one-off parsed
cache (`data/external/ccle_broad_2025/ccle_cache.npz`, 127 MB, sha `1dafee29…`; mirror TSVs pinned
with their sha256 in the runbook). Gene axis read from the bundle's `gene_names.csv` at match time.

| run | wall | peak RSS | call | unanimous 8/8 | anchors K562 / HCT116 / Jurkat (hvg2k, hvg5k Pearson) |
|---|---:|---:|---|---|---|
| `--verify` on A/B/C | 15.7 s | 2.7 GiB | A = JURKAT, B = HELA, C = CAL33 | yes | 1,2 / 1,1 / 2,1 — identical to 041's `results.json` (24/24 ranks) **PASS** |
| synthetic D/E/F (pessimistic) | 14.7 s | 2.7 GiB | D = JURKAT, E = HELA, F = CAL33 | yes | same |
| synthetic D/E/F (realistic) | 10.4 s | 2.7 GiB | same | yes | same |

Minutes, not hours: 041's original spent its wall clock in `pd.read_csv` on the 313 MB TSV.

---

## 3. Item 3 — the census, and what the two synthetic panels say

`scripts/def_census.py` (now with `--cd4-root`) reports per-target carrier counts over the three
production donors, predicts `panel_support` from set intersections, and lists the single-carrier
targets. Known-answer row on A/B/C: 269 / 248 / 276 **PASS** (0.02 s).

| | A/B/C | realistic draw (K562 ∩ Orion ∩ CD4, seed 1990) | pessimistic draw (any donor, seed 199) |
|---|---:|---:|---:|
| K562-gwps / HCT116 / HEK293T | 269 / 290–294 | 300 / 300 / 300 | **150** / 291 / 293 |
| carriers 3 / 2 / 1 / 0 | — | 300 / 0 / 0 / 0 | 141 / 152 / **7** / 0 |
| CD4 Rest / Stim8 / Stim48, union | 283 / 287 / 284, 291 | 291 / 296 / 293, 300 | 173 / 173 / 168, **174** |
| `panel_support` pm / cd4 (predicted) | 248 / 276 | 273 / 286 | 141 / 167 |
| `panel_support` observed by the effects run | 248 / 276 (191, 202) | — | **141 / 167** ✓ |
| panel ∩ essential list | 0 | 34 | 30 |

The pessimistic panel is the stress test: it exercises the 1-carrier path (all-zero agreement key,
alpha = 1 on those targets), 047's CD4 gates (they fire at 173 < 275 and had to be relaxed
explicitly), and the operator's `dropped_in_panel` (30 rows). The realistic panel is what §2.1's
selection criterion predicts Oct 22 will look like. Responses per outcome: runbook §5.

---

## 4. Item 4 — the preflight

`scripts/def_preflight.py`: 66 checks against the §6.1 contract, every value read from the bundle
(manifest keys, per-context counts, the header row on both CSVs, gene order vs each control file's
`var` index, `X` encoding/dtypes/shape/indptr, integrality and non-negativity of every stored value
streamed in 8M blocks, library-size range, NTC count, the density-cap arithmetic), plus what the
frozen table currently pins and whether the two re-pins are needed. 4–5 s, 0.3 GiB. Known-answer
row = AGENTS.md §2's measured numbers, **PASS**. Negative control: the first instance's
`def-bundle-badaxis` (two symbols transposed) fails the order check.

---

## 5. The rehearsal — timings, stage by stage (pessimistic bundle, `work/run1/`)

| stage | wall | peak RSS | evidence |
|---|---:|---:|---|
| preflight | 4.5 s | 0.32 GiB | 66/66 |
| identity | 14.7 s | 2.73 GiB | D/E/F = JURKAT/HELA/CAL33, unanimous |
| cd4-extract (+ source sha256) | 2.2 s (+9 s) | 0.26 GiB | 173/173/168, union 174 |
| census | 0.4 s | 0.10 GiB | predicts 141 / 167 |
| whitener fit (`--max-rss-gib 14`) | **437.6 s** | 8.08 GiB at eigh; 12.88 GiB high-water | sha `68dd20dc…`, 17,942 cols |
| repin-frozen | 2.1 s | 0.23 GiB | 3 contexts out, 3 in, `panel_support → null` |
| pin-artifacts | 0.5 s | 0.05 GiB | 17 lines: 1 registry entry + 15 CD4 hashes |
| repin-panel-support (effects run) | **9.2 s** | **19.80 GiB** | observed 141 / 167 = predicted |
| build (`--emitter rc_hybrid`) | **1287.1 s** | **18.82 GiB** | 13,210,848,728 bytes, nnz 1,647,463,577 (34.7% of cap), sha `28bbe934…` |
| postbuild | 0.2 s | — | 14/14 PASS vs 191's meta (first pass 13/14: the intended `expect_contexts` key, now allowed) |
| `vcc prep --dry-run --contexts D,E,F` | **13.9 s** | 21.69 GiB | ✓ dry run; 360,000 cells; D/E/F = 120,000 each; targets verified; counts-preserved |
| build (`--emitter rc_callset`, the incumbent) | **1255.8 s** | **19.89 GiB** | 12,623,641,976 bytes, nnz 1,574,950,396 (33.2%), sha `9d8cf570…`; postbuild **14/14 vs 202's meta** (diff = `contexts`, `expect_contexts`, `whitener_id`); prep 11.8 s at 23.84 GiB, clean |
| **total, bundle → prepped artifact** | **1,748 s = 29 min** | | |

Known-answer operator refit on the A/B/C bundle: **409.2 s, `e2f65ecf…` MATCH** (process high-water
10.95 GiB). The first attempt of both fits at a 10 GiB guard aborted at the covariance step at
10.02 GiB because the guard counts the resident 4.2 GB basis memmap; 14 GiB is the documented
setting (the anonymous footprint is ~8 GiB).

---

## 6. Items 5/6/7 — the decision table and the checklist

Runbook §3 (decision table), §5 (coverage responses), §6 (submission discipline with the
time budget: a22's 3.5 GiB upload took 15m31s at 3.8 MiB/s; one re-submission ≈ 1 h before scoring
starts; nothing new after Nov 5 20:00 UTC). Two rows deserve a note here:

- **`VCC_ALPHA`.** The coordinator's decision is global 0.982 (board-validated; abc-knee falsified the
  per-context knee at Δ`mse` −0.00226). 187's D/E/F advice — ship `knee × (1 − 1.25%)`, +0.0009
  expected — is recorded as the alternative with a numeric trigger (a context whose blind knee alpha
  exceeds 1.13, i.e. sits below ρ ≈ 0.80 under 0.982, where 180 §8 prices the forfeit at ≥ 0.018
  scaled `mse`). Below that trigger the two rules differ by less than the proxy error and 0.982 ships.
- **Identity-selected donors.** The one time a CCLE identity selected a donor, it shipped as a
  board-validated residual (041 → 062 → a15). The pre-decided rule is that identity may add a donor
  only through a mechanism that already has a board point; a new mechanism cannot be priced on D/E/F.

---

## 7. Rule 10 — bit-identity of the A/B/C build, and the sweep-commit incident

`rule10_199.py` (incumbent now `202-rc_callset.h5ad.meta.json`, the newest meta written by the
HEAD-of-`ab8aa3e` builder): **PASS.** [1] all 50 `hashes.implementation` + `hashes.inputs` fields of the incumbent's meta recompute from disk with 0 differences — a build started now writes exactly that hash block; [2] all 9 production files are byte-identical to the launch commit `ab8aa3e` (the builder at `6cc91501…`, which is the hash 202's meta records); [3] `git grep` finds 0 references to any of the six new scripts from `scripts/` or `src/` other than from each other; [4] `python -m pytest scripts/test_make_submission_pm_cd4.py`: 65 passed (62.8 s, including the test that resolves the pre-`--emitter` builder out of git and compares the emitted `.h5ad` sha256 through the real writer); [5] `scripts/frozen_inputs.json` sha256 `a654d66b…`, unchanged. Record: `work/rule10.json`.

**Incident.** The rehearsal's `pin-artifacts` stage edits the builder's source (by design: the
operator and CD4 hashes live there) and the driver byte-restores it after `prep`. Between the
first pin (21:2x) and the restore, the coordinator's sweep commit `04101c5` (21:16:20 +03) captured
the pinned builder — `scripts/make_submission_pm_cd4.py` now differs from `ab8aa3e` by 29 lines at
HEAD: a `union_pruned_def199` registry entry (harmless, additive) and **the 15 CD4 hashes pointing
at the synthetic banks** (an A/B/C build at HEAD aborts loudly at "frozen model input hash drift for
cd4/…"; it cannot mis-build). The working tree has been restored to the `ab8aa3e` bytes
(`6cc91501…`), so `git status` shows the builder as modified relative to HEAD (`git diff --stat HEAD --
scripts/make_submission_pm_cd4.py`: 12 insertions, 17 deletions = the pin removed; `git diff ab8aa3e`
on it is empty). **The coordinator should commit that restore before anyone builds A/B/C at HEAD.** The pin script's `--restore` now restores from the byte copy saved at
pin time (`--record`), not from git HEAD, precisely because HEAD can move during a run.

---

## 8. What was NOT done, and why

- **No production behaviour change.** The right long-term fix is to move the operator and CD4
  pins into `scripts/frozen_inputs.json` next to the bundle half (180's own reasoning applies to
  them). That is a builder change with a full A/B/C bit-identity proof; it was outside this
  mandate, so the day uses a scripted, diffed source edit instead. Recommended before Oct 22.
- **The heavy stages exceed the 10 GB guidance** (18.82 GiB). Stated, not worked around: the
  floor is the three-donor model load (180 §11); the driver enforces "alone, and ≥ 26 GiB free".
- **The whitener was not re-fit through the driver** (it ran in a chained background job with the
  same command); its 437.6 s is carried into the total by hand.
- **No slot spent, no `vcc submit`, no `cell-eval2`.** The rehearsal artifacts are in `work/run1/`
  and are not submissions.
- **The realistic bundle ran only the cheap stages** (23 s); its heavy stages would measure the
  same clock.

---

## 9. Files

| file | what |
|---|---|
| `docs/DEF-RUNBOOK.md` | the day, numbered; decision table; coverage responses; checklist; pins |
| `scripts/def_pipeline.py` | the driver: 11 stages, timed, fail-closed, memory-guarded, never submits |
| `scripts/def_pin_panel_artifacts.py` | re-pins the operator (additive registry entry) and CD4 hashes in the builder; `--restore --record` |
| `scripts/def_identity.py` / `def_preflight.py` / `def_census.py` / `def_postbuild.py` | stages (first instance; census gained `--cd4-root`) |
| `defbundle199.py`, `defbundle199r.py` | the pessimistic and realistic synthetic bundles |
| `rule10_199.py`, `work/rule10.json` | the bit-identity proof |
| `work/run1/` | rehearsal logs, `timings.json`, meta, postbuild, pin record + pre-pin builder copy |
| `work/run-r/` | realistic bundle, cheap stages |
| `work/whitener_union_pruned_def199.npz` (+`.json`), `work/whitener_union_pruned_abc_ka.npz` | the re-fit and the known-answer refit |
| `work/cd4_compact_abc_ka/` | the byte-identical re-extraction of the shipped banks |
| `work/GWCD4i.DE_stats.h5ad` | the pinned CD4 source, 16.8 GB, `c355f535…` |
