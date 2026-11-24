# 551 — C10b, R6's last access (running record)

Owner: session "Explore contribution areas" (2026-09-30). Preregistration: `PREREGISTRATION.md`. Ledger: `docs/CONFIRMATION-LEDGER.md`,
R6 / C10b. Code: `lib551.py` (lib534 + per-target rows), `designs551.py`, `full551.py`, `comp551.py`, `price551.py`,
`aggregate551.py`, `gate551.py`, `drive551.py`. Nothing here touches R6 before the configuration F is recorded in the ledger.

## Status
- 2026-09-30 21:23Z: C10b timing amendment in the ledger (after the 10-04 read; da543a7c).
- 2026-09-30 ~22:00Z: code and preregistration written; the efp known answers (`drive551.py ka`) wait for the gate
  (549's build, then ≥ 27 GiB free+inactive).
- 2026-10-01 00:32Z: **efp known answers PASS, |d| = 0.0 everywhere** (`drive551.py ka`; logs `work/logs/full551-efp.log`,
  `comp551-efp.log`; peak RSS 5.3 / 7.6 GiB). Two false starts first, both fixed and committed: 062's 9 GiB guard when one
  process computed three donor banks' statistics (now one bank per process, each byte-identical to 534's cache), and my own
  KA comparing `key410` without `equal_nan` (5 targets carry no 410 key; identical otherwise).

| known answer (efp, 3 cancer folds) | result |
|---|---|
| KA-FOLD: fold caches == 438's | bit for bit, 3/3 |
| KA-438 / KA-493 (M4) / KA-546 (M5) / KA-522 (M6, + hi, AB, K_mean) | |d| 0.0, 3/3 |
| KA-546o: P 96.647 at q 1.75, + P2g, + S_P2 (design and realised 16) | |d| 0.0, 3/3 |
| KA-493h: V1_5_3 + P2g_carried, + S_P2 (design and realised 16) | |d| 0.0, 3/3 |
| KA-534: every 534 design and levelled twin == 534's efp comp534 | |d| 0.0, 3/3 |
| full551: KA-COMP (vs ef_b1 meta) | worst rel 8.9e-06 (534's own: same) |
| full551: KA-KEY (410 key == c3e429dd), KA-CLAIMS, KA-E438, KA-534FULL (energies, ρ, level factors, keys) | 0.0 |
| full551: KA-PT0 (per-target path at constant q == scalar path, M5 and M6) | bit for bit |

The machinery is ready for the access (10-04 05:00–12:00Z): `record_config551.py` with the coordinator's message, then
`drive551.py access` (operators re-fitted with `--expect-sha256`, KA-REPRO against C10a before any new number).


## R7a driver (added 2026-10-02; not part of C10b, changes no C10b file)
`r7a551.py` runs this instrument on 552's reserve R7a (T1–T4) after C10b, with fresh operator fits (no C10a values exist
for them); `aggregate_r7a551.py` aggregates the claim (552's `CLAIM-R7-DRAFT.md`). Both live here because guard551 confines
the instrument's outputs to this directory (T outputs under `work/T1..T4`, the driver's state and logs under `work/r7a/`).
Checks on 10-02: `plan` (552's four bundles equal draws.json by sha256; the instrument's 30 repo files, the import closure
over every import statement including function-level ones (lib551 imports ntclevel417 inside a function), are
byte-identical to `e7be21ca`; the package versions are recorded and `access` refuses if they change; uv.lock and the
venv are unchanged since 08-21; no T output exists); KA-CD4 on efp PASS (`work/r7a/ka-cd4.json`: 38 files byte-identical to
238's efp banks, 4 json records equal but for `peak_rss_gib`); the aggregation's smoke test on efp's KA outputs
(`work/r7a/test/`, 200 draws: runs; efp is tuned-on, so its numbers are not a reading). Heavy steps reuse gate551 with
R7a's own window and C10b's floors (module globals patched in-process; gate551's file is unchanged).

## 10-02: a tercilec-vs-tercile row in C10b? Declined, with the reason (the coordinator's question; 652 audit F4)
The audit asks for tercilec's held-out bulk read before October adopts it. **C10b's claim row is that read.** designs551's
P2g is the carried-only form (`lib551.q_p2g`: cut points from the carried targets' key weights, non-carried MID, every
floor target LOW; 493h `P2g_carried`). That is 598's `tercilec`: the same rules, read side by side on 10-02, and 598's A/B/C
known answer moves exactly 493h's 4 targets. With F = M6+P2g (or M5+P2g), C10b prices tercilec(2.5, 2.0, 1.75) against the
uniform q 1.75 on R6's 12 folds, sealed and preregistered, with its permuted control. No new row is needed.

A tercilec-vs-plain-tercile row is declined for three reasons:
- **It would choose between forms on the sealed reserve.** R6's rules say C10b confirms, it does not choose, and "only the
  configuration's own steps are priced" (PREREGISTRATION §2).
- **The instrument has no plain-tercile design.** Adding one means a code change to the instrument, a fresh efp known answer
  (~1–2 h heavy) and a re-pin of R7a's driver, all inside the 36 h before the window.
- **The row would not faithfully price the builder's plain form.** That form depends on how the builder resolves ties: the
  non-carried targets sit at v = r/c = 1 ± 1e-14, and a tie block that straddles the lower cut is split by float noise
  (493g's emptied MID). The instrument's exact key weights cannot reproduce that.

The form comparison already exists off the reserve: 493h on 438's October-like panels (efp: carried +0.84 % rel cos, 3/3
folds; plain +0.48 %, 2/3), and on A/B/C the two forms differ on 4 of 300 targets.
