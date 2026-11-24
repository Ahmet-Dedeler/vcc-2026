# 539 — How the leaders climbed (board snapshots, 2026-09-15 → 09-30)

Coordinator, 2026-09-30. `traj539.py` reads only our own leaderboard snapshots (`data/leaderboard/*.json`, every distinct
entry per team) and decomposes each team's climb since 09-15 into the six scaled members (overall-equivalent = scaled Δ / 6).
Output: `traj539.json`. No network, no new board data.

| team | 09-15 → last | climb | entries | median step | pds | mse | nmae | fid | reach | jac |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Illumina AI | 0.2763 → 0.4192 | +0.1429 | 26 | +0.0076 | +0.0099 | **+0.0662** | +0.0225 | +0.0232 | +0.0134 | +0.0077 |
| cqawesome | 0.2708 → 0.3276 | +0.0568 | 26 | +0.0020 | +0.0013 | +0.0173 | +0.0158 | +0.0173 | +0.0037 | +0.0014 |
| Testh | 0.2151 → 0.3061 | +0.0910 | 17 | +0.0046 | +0.0062 | +0.0228 | **+0.0287** | +0.0083 | +0.0240 | +0.0009 |
| nailong king | 0.2435 → 0.2927 | +0.0492 | 17 | +0.0005 | +0.0045 | +0.0208 | +0.0174 | −0.0004 | +0.0039 | +0.0031 |
| us (displayed rows) | 0.2224 → 0.2651 | +0.0427 | 11 | +0.0024 | +0.0046 | +0.0296 | +0.0108 | −0.0021 | −0.0002 | −0.0001 |

(Our row is the displayed = latest entry, 458g on 09-30; our best is 446bj2z 0.2728260.)

## What it says
- **Every leader climbed in many small steps,** 17–26 entries in two weeks with median steps of +0.0005 … +0.0076. None
  jumped from a single new data source. Illumina's first entry on record (09-11, `iCell`) was already 0.287; its lead grew
  from ~0 to +0.146 over 09-15 → 09-30.
- **The largest share of each climb is mse** (Illumina +0.066 of +0.143) with small pds moves: magnitude allocation at
  roughly fixed direction. cqawesome's pds has been flat at 0.824–0.833 since 09-13 while its mse rose 0.35 → 0.45.
  Testh's climb is mostly DE (nmae +0.029, reach +0.024).
- **Why this matters for the final:** a parameter found by stepping on the A/B/C board (per-context levels, magnitude
  allocations, call-set sizes) cannot be stepped on the dark D/E/F board. Our October path uses no A/B/C board information
  beyond global choices (453dr: the path costs nothing but its CD4 dose on A/B/C; 534 Part B: the bulk recipe keeps 0.70 of
  its A/B/C gain on sealed October-like panels). The leaders' D/E/F retention is unknown; ours is measured.
- Limits: members are averages over three contexts; a step's cause is not observable; teams may hold stronger entries than
  the displayed ones.
