# 467 — Arc's programme gated to design-strong targets (H1, production emitter)

**Result.** On H1, giving Arc's clean H1 programme only to the design-strong two thirds of targets, and CLEAN4 to the
weakest third, keeps most of the programme's fid/nmae gain and removes its reach cost. Overall-equivalent: +0.0110
gated vs +0.0067 ungated. Weak targets want CLEAN4: not Arc's programme, and not the design's own signs.

**Uncertainty (added 20:13Z 09-26 after a peer review; `bootstrap467.py`).**
- Over h1-test's 85 targets, the gated programme vs CLEAN4 is robust: +0.0110, 95 % [+0.0067, +0.0158], P(≤ 0) 0.000.
- **The gate's own benefit (gated − ungated) is not established on H1:** +0.0043, 95 % [−0.0034, +0.0132],
  P(≤ 0) 0.14–0.17 on every seed. It rests on four targets. Without CNN3, COL6A1, ZNF664 and AASS it reverses: ungated
  +0.0152, gated +0.0115.
- The seed sd quoted earlier (0.0002) measures emission re-rolls only, and understates the relevant uncertainty about
  20×. The board pair 458g vs 458d (~250 weak target-contexts, exact split) is the real test.

## Config
- Harness: 435's h1-test run (`run435.py`, `RUN426_EMITTER=prod`, `--bulk ps --mde-c 0.1 --mde-cap 0.1
  --mde-min-count 0.3 --mde-coarse`), seeds 0–2. Deltas are raw, vs 435's CLEAN4 (q80) run at the same seed.
- Programmes:
  - Arc = P_down over `data/donors/h1_2025_full` (all public VCC 2025 H1 cells, 458d's bank) with h1-test's own targets
    removed;
  - CLEAN4 = the incumbent `prior_signs_full`.
- Strength per target = the count of genes with |key| > τ, where `key` is the emitter's gated design row.
  - `record` mode measured τ 0.05 / 0.1 / 0.2 and the L2 norm. It chose n_gt_0.2 by Spearman with h1-test's n_real:
    0.48 (others 0.13 / 0.17 / 0.20). The measure was fixed before any gated number existed.
- Gate: a target is strong if its strength ≥ the q-quantile over the 85 targets.
- Known answer: `record` (CLEAN4 through the wrapper) equals 435 seed 0 bit for bit (`work/ka467.txt` PASS).
- Peak RSS 7.1–7.6 GB, about 3 min per run. Logs in `work/logs/`, summary in `work/summary467.json`
  (`summary467.py`).

## Numbers (raw Δ vs CLEAN4, mean of 3 seeds; pds and mse Δ exactly 0.0 in every arm)
| arm | strong targets get | weak targets get | fid | nmae | reach | jac | overall-equiv |
|---|---|---|---|---|---|---|---|
| ungated (462l) | Arc | Arc | +0.0198 | −0.0148 | −0.0388 | −0.0026 | +0.0067 |
| **gated33** | Arc (57) | CLEAN4 (28) | +0.0161 | −0.0102 | +0.0015 | −0.0022 | **+0.0110** |
| gated50 | Arc (45) | CLEAN4 (40) | +0.0126 | −0.0076 | +0.0029 | −0.0012 | +0.0091 |
| gated67 | Arc (31) | CLEAN4 (54) | +0.0067 | −0.0041 | +0.0017 | −0.0002 | +0.0051 |
| arcoff33 | Arc (57) | design signs (28) | +0.0139 | −0.0124 | −0.0099 | −0.0034 | +0.0078 |
| c4off33 | CLEAN4 (57) | design signs (28) | −0.0021 | −0.0022 | −0.0114 | −0.0012 | −0.0032 |
| c4off50 | CLEAN4 (45) | design signs (40) | −0.0028 | −0.0017 | −0.0203 | −0.0024 | −0.0059 |

- "Design signs" means the weak block keeps the design's own signs: `protect_q` 0, so only the block's minimum-|E0| gene is
  flippable. The coarse MDE and the shrink floor are unchanged.
- overall-equiv = Σ Δraw / span / 6, with H1-like spans: pds 0.45, mse −0.955, nmae −0.60, fid 0.30, reach 0.90,
  jac 0.37.
- Paired (same seed) gated − ungated: q33 +0.0043 (seed sd 0.0002; target-bootstrap 95 % [−0.0034, +0.0132], see the
  uncertainty note above), q50 +0.0024, q67 −0.0016.
- Seed sd of one arm's overall-equiv: 0.0004–0.0007. reach carries most of it (sd 0.002–0.003 raw).

## Mechanism
- The q33 gate marks all four of 462 §11's reach-loss targets as weak: CNN3, COL6A1, ZNF664 and AASS each have strength 1
  (n_real 2, 2, 8, 17).
- Weak set at q33: 28 targets, median n_real 9, against 342 for the strong set. 64 % of weak targets have n_real ≤ 20
  (strong: 30 %). By n_real tercile, the weak set holds 15 / 10 / 3 targets.
- For H1's weak targets' reach, the order is CLEAN4 > design signs > Arc's programme:
  - with Arc on strong targets: +0.0015, −0.0099, −0.0388;
  - with CLEAN4 on strong targets: 0 (base), then −0.0114 for design signs.
- A consensus prior is right more often than the line's own programme on the few DE genes of weak targets. Those genes
  run against the line's common programme (462 §11). The design's own signs on weak targets are noisy.
- H1's strength is coarse (integers; 5 targets at 0, 23 at 1, 12 at 2). Any q from 0.05 to 0.30 gives the same 5-target
  weak set, and q 0.33–0.40 give the same 28. q < 0.33 was not run.
- **Per target (17:35Z):** on strong targets the ungated and gated arms are identical (the same prior): reach +0.0022,
  fid +0.0239 mean vs CLEAN4. On weak targets Arc gains fid +0.0113 mean but loses reach −0.127 mean (sum −3.31).
  - Nearly all of the loss is 4 targets: COL6A1 −1.00, CNN3 −1.00, ZNF664 −0.71, AASS −0.69. No weak target gains reach.
- **Why those four:** their reference DE genes are almost all DOWN (target excluded): AASS 15/16, ZNF664 7/7, CNN3 1/1,
  COL6A1 1/1 (ZNF32 5/6).
  - Across the weak third, DE genes are down 0.631 of the time per target, 15 of 26 majority-down. The strong two
    thirds are near balance (0.554 per target, 0.522 pooled).
  - CLEAN4 is 84 % down on the full axis and keeps those genes down. Arc's programme is 55 % down and flips some.
  - So a weak CRISPRi response is mostly down, while a strong one needs a real programme. 470's `wdown` tests the
    extreme.

## Caveats
- The measure (1 of 4) and q (1 of 3) were chosen on H1, so the H1 gain is in-sample. The board is the out-of-sample test.
- `docs/H1-BOARD-CALIBRATION.md`:
  - DE sign priors carried 0.07–0.67× of their H1 fid to the board;
  - a sign prior's reach is not priceable on H1;
  - so the size of the reach recovery on A/B/C is unknown.
- On A/B/C the design rows are the donor-transfer model's. Their strength distribution differs from H1's, and the
  builder knob uses each context's own quantile.
- A better strength measure would need information the panel lacks. On H1, K562-essential membership predicts n_real
  about as well as design strength (Spearman 0.49 vs 0.48; 0.59 combined; essential median n_real 2,068 vs 17). But the
  2026 panel is essential-free (436), and none of the four reach-loss targets is essential.

## What it changed
- Builder knob `VCC_DE_PRIOR_GATE=q<NN>` (468; weak targets fall back to CLEAN4).
- Arm 458g = 458d + q33, for 09-30 s2 beside 458d (`468-prior-gate/PREREGISTRATION.md`).

## Panel-like subset: H1's non-essential targets (18:43Z 09-26)
The 2026 panel is essential-free, while h1-test has 18 K562-essential targets of 85. The same per-target metrics were
restricted to the 67 non-essential targets. Raw Δ vs CLEAN4 (3 seeds); DE-only overall-equiv as above.

| arm | subset | fid | reach | nmae | jac | overall-equiv |
|---|---|---:|---:|---:|---:|---:|
| ungated | non-essential (67) | +0.0138 | −0.0497 | −0.0149 | −0.0024 | **+0.0015** |
| gated33 | non-essential (67) | +0.0103 | +0.0020 | −0.0092 | −0.0018 | **+0.0078** |
| gated50 | non-essential (67) | +0.0082 | +0.0038 | −0.0084 | −0.0015 | +0.0069 |
| gated67 | non-essential (67) | +0.0035 | +0.0022 | −0.0043 | −0.0002 | +0.0034 |
| ungated | essential (18) | +0.0418 | −0.0002 | −0.0145 | −0.0033 | +0.0257 |
| gated33 | essential (18) | +0.0375 | −0.0002 | −0.0126 | −0.0036 | +0.0227 |

- On an essential-free panel, the ungated programme nets almost nothing on H1 (+0.0015): the reach cost eats the fid
  gain. The gate keeps +0.0078.
- Most of the programme's fid gain on the full H1 panel comes from essential targets: +0.042 per target vs +0.014.
- So 458d (ungated) may read near null on the board even if the programme carries, and 458g is the arm that isolates
  the useful part.
