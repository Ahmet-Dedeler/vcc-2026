# 256 — The CD4 deflation is WIRED into the October driver: the re-fit `VCC_CD4_DEFLATE_K=5` on A/B/C, the PINNED basis `VCC_CD4_DEFLATE_BASIS=V5_abc` on every other panel; `def_postbuild.py` asserts it by value and by panel with `cd4svd5` as the incumbent. **Known-answer row 1 (240's, re-run): the driver's DEFAULTS rebuild the scored `cd4svd5` h5ad byte-for-byte** (`aba27758…`, 12,648,807,800 B), postbuild 34/34 inside the driver. The `.vcc` re-prepped under `vcc` 0.2.0 is content-identical to the scored package and NOT byte-identical — the CLI upgrade changed the container, which is a fact about the novelty guard, fixed here. **WIRED and REHEARSED — verdict line at the bottom (§6).**

**Date 2026-09-05 (01:35–01:50 UTC).** Deterministic. **No `vcc submit`, no model name, no description, no slot.** Two
files changed outside this directory, `scripts/def_pipeline.py` and `scripts/def_postbuild.py` (PRE256 byte copies in
`work/`), plus `scripts/check_submission_novelty.py` (§5) and `submissions/SUBMITTED_LEDGER.json` (backfilled h5ad
digests). The 20 GiB build ran alone under the gate (`memory_pressure` 77 % free, free+inactive 24.9 GiB ≥ 23, zero
other builders, zero `vcc prep`; the binding gate is 20 GiB, AGENTS.md); peak RSS 18.30 GiB, wall 1383 s; prep 18.9 s
at 10.79 GiB (0.2.0's halved prep memory, as its release mail said). Nothing under `submissions/` was written; the
rebuilt h5ad and its container were deleted after hashing (byte / content duplicates), the meta and logs are kept.

| row | required | measured | verdict |
|---|---|---|---|
| **KA 1** 240's row 1 vs the NEW incumbent: driver defaults on `data/challenge`, every `VCC_*` unset | `.h5ad` `aba2775836df3bcf9619a520809fc7a71f325c0a81a881e6f9eb16c970f65faa` on 12,648,807,800 B | **`aba2775836df3bcf…` on 12,648,807,800 B — MATCH**; `[arms]` block: `cd4_deflate_k` asked 5 / carries 5, `cd4_deflate_basis` asked absent / absent, the five 240 rows OK; postbuild 34/34 inside the driver (`work/ka1/logs/postbuild.log`) | **PASS** |
| **KA 1b** `vcc prep` of that h5ad vs the scored package | `6e8b41c6…` (the scored `.vcc`) | **`7f8336fb…` on 2,824,192,000 B vs 2,824,683,520 B — bytes DIFFER, content IDENTICAL**: both containers' `pred.h5ad.zst` decompressed and compared dataset by dataset — `X/data`, `X/indices` byte-identical; `X/indptr` values equal (written `int64` by 0.2.0, `int32` by the old CLI); every `obs`/`var` string and code array equal by value; the new container carries a `meta.json` member `{"cli_version": "0.2.0", …}` the scored one lacks. `vcc` was upgraded to 0.2.0 at 22:35 local on 09-04; the scored package was prepped before that | **PASS on content; the byte row is a CLI-version fact (§4)** |
| **KA 2** `def_postbuild.py` on the 249 meta against itself | 34/34 | **34/34** (`work/pb249.json`): `config.cd4_deflate_k == 5` and basis-absent rows PASS, `moved=[]` | **PASS** |
| **KA 3** the guard demonstrated FIRING: the 222 (`h1union`) meta as the build | fails on the deflation row and the config diff, nothing else | **31/33**: FAIL `config.cd4_deflate_k == 5 on A/B/C: '<absent>'`; FAIL `config differences … moved=['cd4_deflate_k']` (`work/pb222.json`) | **PASS** |
| **KA 4** 251's D/E/F rehearsal meta (`lantern-under-the-stairs`, ef-bundle) against the 249 incumbent | every deflation row PASS; the re-fit → pin swap allowed in the rule-10 diff | **38/41**: `cd4_deflate_basis == 'V5_abc' on a non-A/B/C panel` PASS, `cd4_deflate_k ABSENT while a pinned basis is on` PASS, config diff allowed set now includes `cd4_deflate_basis` and `cd4_deflate_k` PASS; the 3 failures are the operator-registry rows (`union_pruned_def` is pinned per rehearsal and restored, 251) | **PASS (expected failures only)** |
| **KA 5** driver `--dry-run` defaults | A/B/C: `VCC_CD4_DEFLATE_K '5'`, no basis; ef-bundle: `VCC_CD4_DEFLATE_BASIS 'V5_abc'`, `VCC_CD4_SIZING 'carrying'`, no k | as required, both panels | **PASS** |

## 1. The change
`scripts/def_pipeline.py`: `ABC_ENV = {"VCC_CD4_DEFLATE_K": "5"}` merged into the day's defaults on A/B/C only;
`DEF_ENV` gains `"VCC_CD4_DEFLATE_BASIS": "V5_abc"` for every other panel; `day_env` composes them; `ARM_BOARD_DELTA`
gains both variables at **+0.0042634 (`cd4svd5`, entry `Xv0ayvuoHpsKYWv2iu42`)**; `assert_build_carries_arms` re-reads
`config.cd4_deflate_k` and `config.cd4_deflate_basis` off the finished artifact (asked value, or asserted absence);
`--incumbent-meta` defaults to `submissions/249-cd4svd5-f034.h5ad.meta.json`; the A/B/C print block names the basis as
unset there. The builder refuses both variables at once (rule 10), which is why the two tables exist.

`scripts/def_postbuild.py`: `SHIPPED_META`/`SHIPPED_BOARD` → `cd4svd5`; `SHIPPED_CD4_DEFLATE_K = 5`,
`SHIPPED_CD4_DEFLATE_BASIS = "V5_abc"`; on A/B/C `config.cd4_deflate_k == 5` is required and a basis is refused; on a new
panel `config.cd4_deflate_basis == "V5_abc"` is required and, when the reference carries `cd4_deflate_k == 5`, that key
is allowed to disappear in the rule-10 diff (the intended swap of the same k); a reference without `cd4_deflate_k == 5`
is reported stale. `--allow-arm-drift VCC_CD4_DEFLATE_K` / `VCC_CD4_DEFLATE_BASIS` waives the row, as for the other arms.

## 2. What is now true for October
The runbook's row (RB §3.1, 250's trigger, fired by BR249) is executed by defaults: an operator who exports nothing
builds D/E/F with the pinned A/B/C basis, and the two guards (driver `[arms]`, `def_postbuild`) refuse a D/E/F build that
dropped it or carried the re-fit instead. 240's row 1 holds for the new incumbent. The `DEF-RUNBOOK.md` table row and
§0.0 still describe the flag as "unset until the point is read" and must be rewritten SHIPS (a docs edit, §6).

## 3. KA 6 — the D/E/F rehearsal with NOTHING exported (done 02:17–02:40Z, `rehearsal256.sh`, chained behind 258)
251's `rehearsal251.sh` with its one `export` removed and every `VCC_*` unset: the driver on 238's ef-bundle, 244's
operator, `--cd4-thresholds relaxed`, stages preflight → … → build → postbuild (rc 0, 22.4 min, peak 19.42 GiB), then
prep `--dry-run --contexts D,E,F` + `--restore-builder-after` (rc 0, 10.96 GiB). **`[arms]`: `cd4_deflate_basis` asked
`'V5_abc'` / carries `'V5_abc'`, `cd4_deflate_k` asked absent / absent, `cd4_sizing` `carrying`, `cd4_target_g`
0.3443655278649641, the four 240 rows OK. `def_postbuild` 41/41** (the operator-registry rows pass here because the
run pins `union_pruned_def` for its own duration). The meta's deflation block: `mode pinned`, `basis V5_abc`,
`captured_frac_of_centred` **0.1334** (250's band 0.10–0.14). Builder restored to `3615d20f…` (= intended), `git
status scripts/` empty. Artifact `work/rehearsal/umbrella-in-the-attic.h5ad` (12 GB, kept as 251 kept its own; delete
freely), meta and logs committed. **The October path now carries the deflation by default, end to end.**

## 4. The `.vcc` byte row, and why it matters beyond this experiment
A byte-identical h5ad now preps to a **different `.vcc` digest** than before 22:35 local 09-04. Two consequences:
- **Board provenance by `.vcc` sha is CLI-version-dependent.** The pairs on record: `h1union` `8df4dadd…` → `f13ce40b…`
  (old CLI), `cd4svd5` `aba27758…` → `6e8b41c6…` (old CLI) → `7f8336fb…` (0.2.0, content-identical). Record the h5ad
  digest as the identity of a prediction; the `.vcc` digest identifies an upload.
- **The novelty guard had a hole.** `scripts/check_submission_novelty.py` keyed the ledger on the `.vcc` digest only, so
  re-prepping a scored h5ad under 0.2.0 would have printed `OK: novel package` and let a slot be spent on a known
  answer — the exact mistake the guard exists for (2026-08-31). §5.

## 5. The novelty guard, closed
`check` now also derives the package's prediction identity — `output.sha256` from the sibling `<stem>.h5ad.meta.json`
(the builder's own record; no 12 GB hashing), or the sibling `<stem>.h5ad` hashed when there is no meta — and refuses
when that digest is in any ledger record's `h5ad_sha256`. `record` stores `h5ad_sha256` and `vcc_cli_version`. Every
scored record whose meta exists is backfilled (`work/backfill-ledger.log`). Demonstrated: a scratch container with a
novel `.vcc` digest and a sibling meta naming `aba27758…` is REFUSED; the same with a novel h5ad digest passes.

## 6. Verdict
**WIRED and REHEARSED.** `scripts/def_pipeline.py` and `scripts/def_postbuild.py` carry the deflation by panel; the driver's defaults
reproduce the incumbent h5ad byte-for-byte (`aba27758…`) and build D/E/F with the pinned basis with nothing exported (KA 6, 41/41); the six known-answer rows hold; the `.vcc` byte difference is
the 0.2.0 container, content-identical, and the novelty guard now keys on the prediction as well as the upload.
