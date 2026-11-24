# 470 — which block prior the design-weak targets want (H1, production emitter, 467's q33 gate)

**Result (20:33Z 09-26):** nothing beats 467's gated33, CLEAN4 on the weak third and Arc on the strong two thirds.
- All-down on the weak third ties it.
- Arc programmes estimated from the bank's weak or strong half lose on weak targets.
- Arc's strong-half programme on the strong side is slightly ahead (+0.0009 overall-equiv on all 85, +0.0018 on the 67
  non-essential targets), but inside target-sampling noise. By the preregistered rule it is not adopted.
- The builder's gate (CLEAN4 for weak targets) stands, and no new knob is needed.

## Config
- `run470.py` with the 435/467 harness, `RUN426_EMITTER=prod`, 3 seeds per arm, and 467's strength record and q33 (weak
  = strength ≤ 1, 28/85).
- Known answer: `wc4` seed 0 = 467's gated33 seed 0, stats bit-equal (`work/ka470.txt` PASS).
- Priors on the full axis (down fraction; agreement with CLEAN4):

  | prior | down fraction | agreement with CLEAN4 |
  |---|---:|---:|
  | CLEAN4 | 0.837 | — |
  | all-down | 1 | 0.837 |
  | Arc all | 0.549 | 0.595 |
  | Arc weak half | 0.521 | 0.556 |
  | Arc strong half | 0.556 | 0.602 |

  The halves are split by response size: |delta| > 0.2 on the 4,000 most expressed genes, median 17.
- 467's gated33 is the reference. Peak RSS 7.1–7.6 GB, ~3 min per run.

## Numbers (raw Δ vs CLEAN4 everywhere, mean of 3 seeds; pds and mse Δ exactly 0 in every arm)
| arm | strong 2/3 get | weak 1/3 get | fid | reach | nmae | jac | overall-equiv (all 85) | non-essential (67) |
|---|---|---|---:|---:|---:|---:|---:|---:|
| gated33 (467; = `wc4`) | Arc | CLEAN4 | +0.0161 | +0.0015 | −0.0102 | −0.0022 | **+0.0110** | **+0.0078** |
| `wdown` | Arc | all down | +0.0157 | −0.0003 | −0.0106 | −0.0022 | +0.0107 | +0.0072 |
| `warcs` | Arc | Arc strong half | +0.0200 | −0.0257 | −0.0133 | −0.0023 | +0.0090 | +0.0045 |
| `warcw` | Arc | Arc weak half | +0.0174 | −0.0374 | −0.0144 | −0.0034 | +0.0052 | −0.0004 |
| `sarcs` | Arc strong half | CLEAN4 | +0.0150 | +0.0088 | −0.0093 | −0.0012 | +0.0120 | +0.0097 |

Target bootstrap of each arm − gated33 (3,000 draws per seed; seeds agree):
- `sarcs`: +0.0009, 95 % [−0.0022, +0.0045], P(≤ 0) 0.27–0.32. On the non-essential subset +0.0018,
  [−0.0009, +0.0062], P(≤ 0) 0.12–0.16.
- `wdown`: −0.0004 on average, [−0.0038, +0.0016], P(≤ 0) 0.44–0.81.
- `warcs`: −0.0020, [−0.0108, +0.0052].
- `warcw`: −0.0058, [−0.0149, +0.0015], P(≤ 0) 0.93–0.94. Non-essential −0.0082, P(≤ 0) 0.95.

## Against the preregistration
- `wdown`: predicted a tie (reach within ±0.005, overall +0.000 … +0.002). Got −0.0004, reach −0.0018 vs gated33. As
  predicted. CLEAN4 is 84 % down on the full axis, so the two priors are nearly the same on weak targets.
- `warcw`, `warcs`: predicted reach 0.02–0.04 below gated33. Got −0.039 and −0.027. As predicted.
  - A per-gene programme from this line misleads its weak targets, whichever half of the bank it is estimated from.
  - The weak half is worse: its signs are 52 % down, the least down-leaning of the priors.
- `sarcs`: predicted fid within ±0.002 of gated33. Got −0.0011, with reach +0.0073 and jac +0.0010 better.
  - The rule (≥ gated33 + 0.001 overall) is missed at +0.0009, and the target interval straddles 0. Not adopted.
  - It is the best strong-side programme seen on H1. 470b (Arc's non-essential vs essential perturbations) is the
    follow-up question.

## Side analysis: weak CRISPRi responses lean down in every line (`downfrac470.py`)
Proxy DE sets (|delta| > 0.3 on each bank's 2,000 most expressed genes) give a mean down fraction for weak perturbations
(1–20 genes) of 0.58–0.67:
- K562 0.652, RPE1 0.633, Jurkat 0.666, HCT116 0.580, HEK293T 0.624, H1 0.605.

Strong responses also lean down in most lines (K562 is the exception at 0.52). Manuscript C-72.

## What it changed
- Nothing in the builder: CLEAN4 stays the weak targets' prior. The gate stays as built, and its October form becomes
  absolute (`le1`, 472).
- The strong side's programme is an open question for 470b.

## 470b (22:13Z 09-26): the strong side's programme by essentiality — no gain
Weak third on CLEAN4 as in gated33; the strong two thirds get Arc's programme from the 155 non-essential (`sarcne`) or
the 60 essential (`sarcess`) non-h1-test perturbations. 3 seeds each.

| arm | fid | reach | nmae | jac | overall-equiv, all 85 | non-essential 67 |
|---|---:|---:|---:|---:|---:|---:|
| gated33 (Arc, all perturbations) | +0.0161 | +0.0015 | −0.0102 | −0.0022 | +0.0110 | +0.0078 |
| `sarcne` (Arc, non-essential) | +0.0142 | +0.0053 | −0.0095 | −0.0026 | +0.0103 | +0.0080 |
| `sarcess` (Arc, essential) | +0.0137 | +0.0063 | −0.0085 | −0.0024 | +0.0100 | +0.0072 |
| `sarcs` (Arc, strong half; 470) | +0.0150 | +0.0088 | −0.0093 | −0.0012 | +0.0120 | +0.0097 |

- **Against the preregistration:**
  - `sarcne` was predicted within ±0.002 of gated33 on the non-essential subset. Got +0.0002. The rule's
    ≥ +0.001 is not met: not adopted.
  - `sarcess` was predicted fid lower by 0.002–0.006 on the non-essential subset. Got −0.0019 (the edge of the range).
- Splitting Arc's programme by essentiality adds nothing for an essential-free panel. Both halves trade a little fid
  for a little reach, as `sarcs` does, and none separates from the full programme beyond target sampling.
- **Consequence:** the October programme stays the full `h1_2025_full` programme with the gate.
