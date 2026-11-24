# 392 BOARD RESULT — CPM-block signs from the 4-bank cross-line prior through the decoupled DE channel: 0.2226221, Δ +0.00251 vs 319p2 — NEW BEST; fid in its top bin

Entry `TJTOj9qrbiobJBpqLSa3`, submitted 2026-09-13 07:46:36Z (upload 11m11s), published 08:10Z, rank 17 at publication
(13th of 974 on the 2026-09-16 03:14Z snapshot). Package `submissions/392-cap15-partial50-designclean4.vcc` sha256
`20b987a4…`, model name "the stapler went on sabbatical". Full record `data/board-entries/TJTOj9qrbiobJBpqLSa3.json`.
Reference 319p2 `lqWNXX5YaYArdglu6Stn` 0.2201159. Bins frozen in `PREREGISTRATION.md` (commit 175d8454) before the build.

| member | raw 319p2 → 392 | Δ raw | Δ scaled | Δ overall points |
|---|---|---:|---:|---:|
| pds | 0.854738 → 0.854738 | **0 (exact)** | 0 | 0 |
| mse | 0.809898 → 0.809898 | **0 (exact)** | 0 | 0 |
| nmae (lower better) | 0.934784 → 0.938837 | +0.004053 (worse) | −0.006981 | −0.001163 |
| fid | 0.507933 → 0.516157 | **+0.008224** | **+0.026214** | **+0.004369** |
| reach | 0.277798 → 0.273878 | −0.003920 | −0.004182 | −0.000697 |
| jac | 0.040903 → 0.040898 | −0.000006 | −0.000015 | −0.000002 |
| **overall** | 0.2201159 → **0.2226221** | | | **+0.0025062** |

## Reading against the preregistration
- **Falsifiers: PASS.** `score_pds` and `score_mse` identical to 319p2 to every printed digit (the DE channel left the bulk
  untouched on the real A/B/C build, as on H1); jac −0.000006 raw, inside ±0.0005.
- **Primary, Δ score_fid = +0.0262 → top bin (≥ +0.020):** the prior's per-gene signs carry A/B/C information and the
  channel delivers it. Predicted +0.041 (H1 1:1); realised 64 % of it — consistent with the addendum's dose diagnostic
  (A/B/C feasible dose median 0.90 vs H1 1.00) and with H1 being neither the strongest nor the weakest line for the prior.
- **Overall Δ +0.0025 → bin "+0.0005 … +0.004: partial transfer; incumbent because fid is in its top bin."** 392 is the
  incumbent and the best-scored package.
- **Secondary members did not follow H1.** nmae worsened +0.0041 raw (H1 predicted +0.0024, unresolved) and reach fell
  −0.0039 raw (H1 predicted +0.0175, unresolved). Both are the members that read the *strongest* block genes: nmae on the
  reference's own significant set, reach on the front of our order. A block gene with a large |E0| is the one most likely
  to be truly DE and to sit early in the reach walk, and the prior's precision there (0.52) is barely above the delivered
  sign's — flipping it at full magnitude costs 2|lfc| when wrong. That is exactly what 394's q80 protection removes.

## Consequences
- 392 is the incumbent (0.2226221). Every further DE-channel arm is registered against it.
- Next one-variable step: `VCC_DE_SIGN_PRIOR=clean4q80` (394: strongest 20 % of block calls keep the delivered sign),
  ported in commit 7f1cddd6, emitter re-pinned 6ff0b866… (clean4 path byte-identical: the diff is the five `protect_q`
  lines only). Registered as experiment 395.
- Housekeeping: nobody was at the machine 2026-09-13 08:11Z → 2026-09-16 05:00Z; the 09-14 and 09-15 slots (four) expired
  unspent. Board moved meanwhile: Illumina 0.2979, cqawesome 0.2741, CellWM 0.2473, nailong king 0.2459, chill 0.2426.
