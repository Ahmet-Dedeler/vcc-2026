# 460dr BOARD RESULT — 0.2692306 (+0.0028694 vs 453dr): the October path's day-of LODO key transfers; DEF_DAY_LODO=1 stays

Entry `LrkhsIXWMjEsvRXSGToa`, submitted 2026-10-01 00:01Z by `464/arm1001.sh` slot 1 (name from the queue, description
empty), uploaded 2.6 GiB in 11 m 50 s, published 00:27Z, captured 00:27:29Z, rank 10. Package
`submissions/460dr-october-path-lodo-abc-q1.75.vcc` (`0a911b57…`). Counterpart 453dr `q3VuTImh4FRGKXH2waCp` 0.2663613.
Bins: `PREREGISTRATION.md` (2026-09-26 10:22Z). Reader: `464/read1001.py lodo` (written 09-30 11:11Z). Read 2026-10-01T00:29Z.

460dr = 453dr (the October path: `def_day.sh` on 180's relabelled A/B/C bundle) + `DEF_DAY_LODO=1` (415's LODO-V1 key
re-fitted on the day). Config diff vs 453dr exactly `bulk_load`, `bulk_load_sha256`, `bulk_ctx_scale` (gate460).

| member | Δ raw vs 453dr | Δ scaled | Δ overall |
|---|---:|---:|---:|
| pds | +0.001044 | +0.002324 | +0.000387 |
| mse | −0.013060 | +0.013736 | +0.002289 |
| nmae | +0.000125 | −0.000220 | −0.000037 |
| fid | +0.000003 | +0.000012 | +0.000002 |
| reach | +0.001223 | +0.001350 | +0.000225 |
| jac | +0.000005 | +0.000015 | +0.000002 |
| **overall** | | | **+0.0028694** |

Scaled members: pds 0.81433, mse 0.36313, nmae 0.17145, fid 0.00547, reach 0.23262, jac 0.02838.

## Against the bins
- The A/B/C key's own board effect (446bj2 − 446b): pds +0.000208, mse −0.014473 raw. Ratios here: **mse 0.90×** (inside
  0.7–1.3), **pds 5.02×** (outside; the denominator is tiny). **BIN: same sign, outside the band: it transfers at a
  different size; keep it ON (overall ≥ +0.001); ratios recorded.** `DEF_DAY_LODO=1` (already DECIDED) stands.
- DE pinned within ballast noise (nmae |Δ| 0.000125 ≤ 0.0015, reach |Δ| 0.0012 ≤ 0.010). In 537's preregistered noise check
  (`537/PREREGISTRATION-noise.md`, 3σ bands |Δfid| ≤ 0.0001, |Δjac| ≤ 0.00005, |Δreach| ≤ 0.0040 raw) this bulk-only pair is
  inside every band.

## What it adds
- On the October path the day-of LODO key is worth +0.0029 overall on A/B/C (vs +0.0022 for the board key on the A/B/C
  path): the mse gain reproduces (0.90×) and pds gains more. 460dr − 446bj2 = −0.0030460: the October path at CD4 ×1 sits
  one CD4 half below production, as 457dr measured (+0.0031 for ×1.5).
- Class (AGENTS §11.1c): transfer. Held-out evidence on record: 438 §3.5 (efp/real/efr) and 534's sealed R6 (+0.00099
  pooled, retention ~0.30 — the weakest October switch off-board).
