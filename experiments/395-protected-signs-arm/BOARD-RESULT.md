# 395 BOARD RESULT — q80-protected block flips on 392: 0.2239363, Δ +0.00131 vs 392 — NEW BEST, incumbent; reach carried it

Entry `RzcGA26xZEigSuXgUp5t`, submitted 2026-09-16 05:50:22Z (upload 17m40s at 2.4 MiB/s), published 06:21Z, rank 16 at
publication. Package `submissions/395-cap15-partial50-designclean4q80.vcc` sha256 `19435f09…`, model name "the thermostat has a
union rep". Full record `data/board-entries/RzcGA26xZEigSuXgUp5t.json`. Reference 392 `TJTOj9qrbiobJBpqLSa3` 0.2226221. Bins
frozen in `PREREGISTRATION.md` (commit d86fd5b3) before the build.

| member | raw 392 → 395 | Δ raw | Δ scaled | Δ overall points |
|---|---|---:|---:|---:|
| pds | 0.854738 → 0.854738 | **0 (exact)** | 0 | 0 |
| mse | 0.809898 → 0.809898 | **0 (exact)** | 0 | 0 |
| nmae (lower better) | 0.938837 → 0.937101 | −0.001736 (better) | +0.002794 | +0.000466 |
| fid | 0.516157 → 0.514228 | −0.001929 | −0.006138 | −0.001023 |
| reach | 0.273878 → 0.283909 | **+0.010031** | **+0.011247** | **+0.001875** |
| jac | 0.040897 → 0.040891 | −0.000007 | −0.000019 | −0.000003 |
| **overall** | 0.2226221 → **0.2239363** | | | **+0.0013142** |

## Reading against the preregistration
- **Falsifiers: PASS** (pds and mse identical to 392 to every digit; jac −0.000007 raw).
- **Primary (nmae/reach recovery):** reach +0.0100 raw clears its +0.003 margin; nmae −0.0017 raw just misses its −0.002 margin →
  bin "exactly one recovers: partial; read which member; q50 only if reach carried it". **Reach carried it** (+0.0019 overall),
  nmae added +0.0005, fid gave back −0.0010 (the protected strongest-20 % block calls include correct prior flips that 392 had
  delivered; 394 forecast fid −0.0008 raw, board −0.0019).
- **Overall +0.00131 → bin "+0.0003 … +0.0015: small gain; incumbent if Δfid raw ≥ −0.002"** → −0.0019 ≥ −0.002: **395 is the
  incumbent**, and the dose step (q50: protect the strongest half) is authorised because reach carried the gain (390 §6 priced q50
  on H1 at fid +0.0081, reach +0.025, nmae −0.0001 vs cap15).
- Prediction check: predicted +0.0011 overall (fid −0.0004, reach +0.0011, nmae +0.0004); realised +0.0013 with the same signs on
  every member — the H1 direct contrast (394) transferred.

## What it proved
The strong-|E0| block calls are where the marginal tilt's flips cost reach on A/B/C: protecting them buys reach at the price of a
third of the fid gain. Reach is the member that rewards the FRONT of the order, so the protection axis is a reach lever. Next
constructions on this channel start from 395 (399 today: per-target neighbourhood tilt, which raises the precision of the flips
themselves; then the q50 dose or a per-target protection quantile).
