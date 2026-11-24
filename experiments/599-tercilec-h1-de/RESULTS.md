# 599 — On H1's real emitter, tercilec costs nothing on the DE side of a thin-coverage panel. Read: **R1, clean**.

- **The read:** tercilec − tercile, post-zs, 5 seeds, f 0.35 gives DE overall-eq **+0.00002** (standard error 0.00006)
  and Δnmae raw **+0.00002**.
- **The cause of the concern:** tercilec raises the ballast by +0.96 mean K (+2.0 on each non-carried target).
- **Bulk:** pds +0.0029 (5/5 seeds), mse +0.0006 (5/5 seeds, slightly worse); approximate bulk overall-eq **+0.0010**.
- **The K model fails here.** 492g's model (+0.0007 nmae) is refuted on this instrument; the null holds.

2026-10-02, session "Explore contribution areas" ([36a81a]), at the coordinator's request (652 F4).
- Registered in `PREREGISTRATION.md`: committed ec25947a before any arm run. The KA amendment (94fcff50) also came
  before any arm run, and so did the extension's K registration (aceebfe8).
- Every heavy step went through gate602: prep 6.0 GB; runs 5–7 GB, ~6 min each.
- No builds, no `vcc`, no `cell-eval2` CLI; the DE members come from lib521's exact functions.

**Class (AGENTS 11.1c):**
- This is a held-out read of the DE side of 598's transfer-rule fix. It is H1, not A/B/C, and nothing is fitted here.
- Tries: one rule pair, 5 seeds, plus one anchor arm.

## Setup (registered)
- **Harness.** run496's 446bj2z construction: production emitter 0f56209d, CallSet(400, 2200, 1.9), CLEAN4 q80, coarse
  MDE, zs. It is driven by monkeypatch (`run599.py`) and nothing is copied.
- **Bulk.** 410's key on the 85 H1 targets (all carried as fitted).
  - Coverage-thinned (seed 599) to **43 non-carried** (w = 1), as on efp.
  - Clip 0.5 / 6, energy-preserving.
  - Sharpened per target by the builder's own `profile_from_effect(..., sharpen_q=q_t, protect=panel)`.
  - q_t comes from 598's patched `sharpen_q_by_tercile`.
- **Arms.**
  - t = `tercile:2.5,2.0,1.75`, tiers 21 / 0 / 64. Both cut points land on the w = 1 block (edges 0.59459 = 0.59459),
    so MID is empty and all 43 non-carried targets go HIGH: efp's failure mode, reproduced.
  - tc = `tercilec:2.5,2.0,1.75`, tiers 14 / 57 / 14, with the 43 non-carried targets in MID.
  - s175 = scalar q 1.75, the anchor.
- **Ballast K** is fixed before emission and was computed in prep.
  - tc − t: ΔK +0.96. The 43 non-carried targets rise by about +2.0 each; the carried ones move −0.15.
  - The K-16 share goes 0.28 → 0.35.
  - **All 13 arm runs used exactly the registered K.**

## Known answers
- **KA-bulk:** mode ka's bulk == run496.ps_bulk, bit for bit. **PASS**.
- **KA-K:** ka's emulated K == 435's recorded ballast_K for every target, seeds 0–2. **PASS**.
- **KA run:** mode ka, seed 0, reproduces 435's base run. Every stats array (q, p, lfc, sums, jk, targets) is bit-equal,
  pre- and post-zs.
  - The pre-zs raw nmae differed in the last bit (…616 vs …618).
  - The scorer is not bit-deterministic on identical inputs: 435's own arrays score …617 in three fresh processes
    (`hashseed599.py`).
  - The KA was amended **before any arm run**: arrays exact, metrics to 1e-12. **PASS**.

## The read: tc − t (post-zs)
| regime | fid | jac | reach | nmae | **DE overall-eq** | per seed |
|---|---:|---:|---:|---:|---:|---|
| f 1 | −0.00003 | −0.00001 | +0.00304 | −0.00004 | +0.00056 | +0.00049, +0.00246, −0.00031, +0.00052, −0.00036 |
| f 0.5 | −0.00007 | −0.00000 | −0.00040 | +0.00005 | −0.00012 | +0.00022, −0.00006, −0.00073, −0.00001, −0.00005 |
| **f 0.35 (primary)** | −0.00003 | −0.00002 | +0.00026 | **+0.00002** | **+0.00002** | −0.00015, +0.00006, +0.00015, +0.00014, −0.00010 |
| cmp | −0.00003 | −0.00001 | +0.00065 | +0.00007 | +0.00008 | +0.00026, +0.00075, −0.00039, +0.00010, −0.00030 |

- **Bulk (the runs' own raw members):**
  - pds +0.00289 (per seed +0.0021, +0.0039, +0.0028, +0.0028, +0.0028);
  - mse +0.00056 (+0.00056, +0.00072, +0.00046, +0.00056, +0.00049);
  - approximate bulk overall-eq **+0.00097**, using the BRIEF's mid anchors.
- **Pre-zs** gives the same picture: f 0.35 +0.00001 (Δnmae +0.00005); bulk identical.
- **The registered 3-seed read** (seeds 0–2) was already R1: f 0.35 +0.00002 (−0.00015, +0.00006, +0.00015), Δnmae
  +0.00007. Its mean sat within 2 standard errors of 0, so the registered extension (t and tc, seeds 3–4) was run, and
  the read above is the 5-seed one.

## Reading (the preregistered rule)
- **R1 (clean): DE overall ≥ −0.0003 and Δnmae ≤ +0.0010. Both hold, with a wide margin.**
  - The DE side does not block adoption. tercilec's net on a thin panel is its bulk gain less ≈ 0.
  - On this H1 panel the bulk gain is itself positive (pds-led).
- **H_K is refuted on the production emitter.** It predicted Δnmae +0.0003 … +0.0010 from 492g's board slope of 0.00073
  per unit of mean K. The measured value is +0.00002 at f 0.35 and −0.00004 at f 1.
- **H_0 holds.** On H1's real emitter, the ballast count itself is not a DE cost at this scale, consistent with 416 (K
  forced: |Δ| ≤ 1e-5) and 490 (the footprint, −0.0001).
  - So 492g's board slope must come from something that travels WITH K on the board (the re-level, the bulk shape), not
    from K itself.
  - S_P2's board failure (+0.0018 raw nmae) is not reproduced by tercilec's K rise here.

## Anchors: each per-target rule against the scalar q 1.75 (3 seeds, post-zs)
| | ΔK | f 1 | f 0.5 | f 0.35 | cmp | pds | mse | bulk ≈ |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| t − s175 | +0.64 | −0.00094 | +0.00015 | −0.00030 | −0.00121 | +0.0011 | −0.0030 | +0.0009 |
| tc − s175 | +1.60 | −0.00007 | −0.00004 | −0.00028 | −0.00100 | +0.0040 | −0.0024 | +0.0019 |

The DE entries are overall-eq. Their driver is reach: t − s175 −0.0018 (f 0.35) to −0.0068 (cmp); tc − s175 −0.0014 to
−0.0054.
- On this thin panel, per-target sharpening itself, in either rule, has a small **reach** cost against the scalar q
  (cmp about −0.001 overall-eq) and a larger bulk gain.
- tercilec is no worse than tercile on that reach cost (better at f 1), and has about twice the bulk gain.
- This is the DE side of `sharpen_pt` itself, not of 598. It is in line with P2g's board read (+0.0010, DE in noise).

## Caveats
- **One H1 panel** (85 targets) with a **random** coverage thinning. efp's non-carried targets are not random.
- **The key.** 410's key at an October-like clip, not the day-of LODO key.
- **No re-level.** The production re-level is absent here, so mse reads the un-levelled render. mse is slightly worse
  under tc on 5/5 seeds (+0.0006); pds carries the bulk gain.
- **Approximate conversion.** The bulk overall-eq uses the BRIEF's mid anchors.
- **Instrument status.** 611 validated the instrument for delivery levers; reach is provisional there. The f 1 / cmp
  reach entries are the noisiest rows (per-seed spread ±0.001).

## Files
- **Code:** `PREREGISTRATION.md` (with its KA amendment), `common599.py`, `prep599.py`, `prepx599.py`, `run599.py`,
  `ka599.py`, `hashseed599.py`, `score599.py`, `drive599.sh`, `drivex599.sh`.
- **Design and registration:** `work/prep599.{json,npz}`, `work/prepx599.npz`, `work/key599.npz`.
- **Scores:** `work/score599.json` (5 seeds; t and tc with seeds 3–4) and `work/score599.3seed.json` (the registered
  3-seed read).
- **Per-run records:** `work/h1_<mode>/metrics-*.json`, `work/h1_<mode>/run599_s*.json`.
- **Logs:** `work/logs/`.
- The stats arrays (`work/h1_*/stats-*.npz`, about 0.9 GB) stay local.
