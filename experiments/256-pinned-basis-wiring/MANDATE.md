# 256 — Wire the CD4 deflation into the October driver: the re-fit on A/B/C, the PINNED basis elsewhere. Known-answer row 1 re-run against the new incumbent.

Written 2026-09-05 01:45 UTC by the coordinator, after the two board points of the night (BR237, BR249) and
`docs/READ-FIRST.md` v2.2. Deterministic; no `vcc submit`, no model name, no description, no slot.

## Why now
`cd4svd5` (`VCC_CD4_DEFLATE_K=5`) is the incumbent at 0.2160944 (entry `Xv0ayvuoHpsKYWv2iu42`), and
`docs/DEF-RUNBOOK.md` §3.1's pre-committed trigger fired (board raw `pds` 0.83572 → 0.85064): the deflation
SHIPS on D/E/F as the PINNED A/B/C basis `V5_abc` (251), never the day-of re-fit (250). Until this experiment
the driver was blind to both variables (R251 §7): its defaults rebuilt `h1union`, and `def_postbuild.py`
asserted `h1union` as the incumbent.

## The change (one commit, `scripts/def_pipeline.py` + `scripts/def_postbuild.py`; PRE256 byte copies in `work/`)
- `def_pipeline.py`: `ABC_ENV = {VCC_CD4_DEFLATE_K: 5}` merged into the defaults on A/B/C only;
  `DEF_ENV` gains `VCC_CD4_DEFLATE_BASIS: V5_abc` for every other panel; `ARM_BOARD_DELTA` gains both;
  `assert_build_carries_arms` gains a `cd4_deflate_k` and a `cd4_deflate_basis` row (asked vs artifact,
  absence asserted when not asked); `--incumbent-meta` defaults to `submissions/249-cd4svd5-f034.h5ad.meta.json`.
- `def_postbuild.py`: `SHIPPED_META`/`SHIPPED_BOARD` → `cd4svd5`; by panel: A/B/C requires
  `config.cd4_deflate_k == 5` and no basis; a new panel requires `config.cd4_deflate_basis == V5_abc` and
  allows the incumbent's `cd4_deflate_k` to be replaced by it in the rule-10 diff (the intended swap, same k);
  a reference without `cd4_deflate_k == 5` is stale.

## Known-answer rows (all must pass; a failure reverts the two files to their PRE256 bytes)
1. **Row 1 of 240, re-run**: the driver's DEFAULTS on `data/challenge` (`--stages build,postbuild`, both deflation
   variables unset in the environment) reproduce the incumbent **byte-for-byte**: `.h5ad` sha256
   `aba2775836df3bcf9619a520809fc7a71f325c0a81a881e6f9eb16c970f65faa` on 12,648,807,800 B; `vcc prep` of it
   `6e8b41c6244b32dbe2a9765c7d2f2293b2032c49e59b455982b95bc6bb074533` (the scored package). `run256ka1.sh`.
2. `def_postbuild.py` on the 249 meta against itself: PASS (34/34) — done 01:38Z, `work/pb249.json`.
3. `def_postbuild.py` on the 222 (`h1union`) meta as the build: FAILS on exactly `cd4_deflate_k` and the config
   diff — the guard demonstrated firing — done 01:38Z, `work/pb222.json`.
4. `def_postbuild.py` on 251's D/E/F rehearsal meta (`lantern-under-the-stairs`, ef-bundle) against the 249
   incumbent: every deflation row PASS, the k → basis swap allowed; the only failures are the three
   operator-registry rows, expected outside a rehearsal (251 pins `union_pruned_def` per run and restores the
   builder) — done 01:40Z, `work/pb251def.json`.
5. Driver `--dry-run` prints `VCC_CD4_DEFLATE_K '5'` and no basis on A/B/C; `VCC_CD4_DEFLATE_BASIS 'V5_abc'` and
   no k on the ef-bundle — done 01:40Z.

## Gate and discipline
Memory gate before the build and before the prep: `memory_pressure` ≥ 40 % free AND free+inactive ≥ 23 GiB
(the binding gate is 20 GiB, AGENTS.md; 251 used 30 with a freer machine) AND zero other builders AND zero
`vcc prep`; the driver's own `--min-free-gib 23`. Strictly serial. Nothing under `submissions/` is written;
the rebuilt artifact lives in `work/ka1/` and is deleted on MATCH (a byte duplicate), kept on mismatch.
Never `git add -A`. The armed-slot scripts are not touched (both slots are spent).
