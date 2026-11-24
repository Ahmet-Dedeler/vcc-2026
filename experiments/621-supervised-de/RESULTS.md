# 621 — Supervised DE direction + membership model (chill-style), leave-one-surface-out: NEGATIVE as a new sign source

> **AUDIT 616 (7dac7a98, [21dd4c], 10-02): the thinned Orion gate tables' DE truth was ~98 % selection artefact.** Here the hct116 / hek293t columns (Orion gate tables) are VOID. On h1 and k562 the HGB direction head still loses to CLEAN4, so the negative stands.

2026-10-02, council session [7772dc]. Portfolio item 2: the last public sign-source candidate before the pivot rule.
Spec from [21dd4c] (chill's day-1 design), instrument [d4cd5f]'s 636 (gates built with the scorer's rules by 480).
Analysis only: no build, no `vcc`, no `cell-eval2`. `build621.py` ran under `gate602.sh` and peaked at 5.26 GB RSS. I
had declared 4 GB, so this was 1.3 GB over; the gate's 20 GiB floor held. The other scripts stayed under 1 GB.

## What was built
- **Features per (target, gene)** on each held-out surface, from banks 636 allows there only (no own line, no
  same-study sibling, no H1 bank on h1):
  - per-target source evidence z_d = δ_d / σ_d(g), where σ is the robust SD over ≤ 3,000 bank rows;
  - a Stouffer z per study family (Weissman / Satija / Xaira / Arc);
  - the family-balanced consensus Zbal, with agreement, max |Z| and the number of families;
  - the family-balanced gene programme and gene responsiveness;
  - the held-out line's own control log-CPM, and target strength.
- **Heads** (HistGradientBoosting): direction P(down | gate) and membership P(gate | tested).
- **Training** uses the other surfaces with every held-out target removed. The Orion sibling is dropped when Orion is
  held out. Each training surface gets equal weight.

## Results (gate sign accuracy, 636's gates)
**1. The HGB direction head fails everywhere.** Paired against CLEAN4 (LOLO) on shared pairs:

| surface | h1 | hct116 | hek293t | k562 |
|---|---:|---:|---:|---:|
| Δ vs CLEAN4 | −0.013 | −0.099 | −0.205 | +0.001 |

Cause: P(down | gate) shifts across surfaces. Orion's gates carry its known technical down-blob; Arc's H1 is
balanced. A learner trained on other surfaces cannot know the held-out base rate.

**2. The membership head** reads AUROC 0.68 (h1), 0.79 (hct116), 0.74 (hek293t), 0.74 (k562), with precision@n_real
0.05–0.18. These rows are in `work/score621.json` for 636's comparison with the call-propensity block.

**3. The robust signal is per-target evidence strength** (`switch621.py`, CLEAN4-covered gate pairs):
- At |Zbal| ≥ 2, the family-balanced consensus sign is right 0.81–1.00 on all four surfaces. CLEAN4 is right 0.60–0.73
  on the same pairs.
- At |Zbal| < 1 (76–80 % of gate pairs on h1 and Orion), the consensus is near a coin flip.
- Switch rule (consensus where |Zbal| ≥ c, CLEAN4 elsewhere) vs CLEAN4:

  | c | h1 | hct116 | hek293t | k562 |
  |---|---:|---:|---:|---:|
  | 1 | +0.030 | −0.004 | −0.023 | +0.125 |
  | 2 | +0.010 | +0.012 | −0.001 | +0.105 |

**4. Against what we SHIP** (`vsprod621.py`, `byclass621.py`). This is h1 only, the one surface with a production
emitter run (435 via 521's lib), on identical gate pairs.
- All targets: production delivered 0.578, switch c1 0.604 (+0.025 [+0.008, +0.043], target bootstrap).
- Non-essential targets (the October-relevant subset, 14,516 pairs): production 0.563, switch c1 0.576
  (+0.012 [−0.007, +0.038]). Not decisive.
- By delivery class (non-essential, c1):

  | class | production | switch | Δ |
  |---|---:|---:|---:|
  | head (top-400 \|key\|) | 0.727 | 0.702 | −0.024 |
  | CPM block | 0.565 | 0.569 | +0.005 |
  | complement (uncalled) | 0.538 | 0.564 | +0.026 |

  - The switch must never touch the head.
  - The complement gain is mostly CLEAN4 replacing production's near-random uncalled signs. That is 478/496's
    uncalled-gene hook, already built and priced weak (446bj2zmz / 446bj2zcph).
  - The genuinely new part, a confident consensus overriding production's sign, is 298 block + 125 complement
    overrides at 58 % / 62 % right. That is ≈ +0.005 gate sign accuracy, about +0.001 overall.

## Verdict
**621 misses the council's bar as a NEW sign source** (≥ +0.01 gate sign over the incumbent on most held-out lines).
Learned heads fail on label shift, and the robust part (per-target evidence) is what production's head already uses.
Per the pivot rule in SYNTHESIS v4, public sign sources close for this round. Class: transfer (negative).

**By-products:**
1. The uncalled-gene CLEAN4 hook is re-confirmed: complement gate signs 0.557 vs production's 0.538 at c2 on
   non-essential H1 targets.
2. Per-target donor evidence above |Z| ≈ 2 is right 81–100 % of the time on every surface. Any future DE design should
   key its confident set on noise-normalised, family-balanced z, not on raw |key| alone. This was not tested against
   q80's choice of genes; that would be a small follow-up.
3. Membership AUROC 0.68–0.79 from public features is available for 636's comparison with call propensity.

## 621b — learned membership vs production's call set (`member621.py`, `work/member621.json`)
h1, non-essential targets. Pooled AUROC for gate membership over tested pairs, and per-target φ (the fraction of calls
that are reference-DE) at K calls:

| ranking | AUROC | φ@400 | φ@1,000 | φ@2,600 | recall@2,600 |
|---|---:|---:|---:|---:|---:|
| production composition (head, block, rest) | 0.609 | 0.061 | 0.086 | 0.069 | 0.530 |
| control log-CPM alone | 0.611 | 0.116 | 0.091 | 0.068 | 0.440 |
| \|Zbal\| | 0.527 | 0.062 | 0.054 | 0.047 | 0.383 |
| 621 HGB membership | **0.641** | 0.116 | 0.089 | 0.067 | 0.460 |

The learned head ranks best overall, but at our call-set size (2,600) every ranking reaches φ ≈ 0.07. **Not a call-set
lever** (closes 3(a)'s learned-membership variant). Production's 400-gene head is chosen for sign, not membership: its
φ is half that of the CPM ranking.

## 621c — direction-confidence REACH ORDER, analytic ceiling (`reach621.py`, `work/reach621.json`)
Sweep candidate #2. Production's delivered signs and calls are pinned (435's emitter on h1-test, via 521's lib), and
only the scorer's order inside the called set changes, as if an emitter realised it perfectly. Δ reach vs today's
realised order, with target-bootstrap 95 % CIs:

| order inside the called set | f 1, all | f 1, non-ess | f 0.35, all | f 0.35, non-ess |
|---|---:|---:|---:|---:|
| \|Zbal\| desc | −0.012 [−0.026, +0.000] | −0.015 | −0.001 | −0.005 |
| agreement-confidence (+\|Z\| if consensus agrees with our sign, −\|Z\| if not) | +0.007 [−0.005, +0.021] | +0.006 | +0.009 [−0.005, +0.024] | +0.003 |
| same, applied to the uncalled tail too | +0.006 | +0.005 | +0.001 | −0.009 |
| ORACLE (right-signed calls first) | +0.176 | +0.154 | +0.198 | +0.165 |

- **Closed without emitter runs.** The best realisable order is worth ≤ +0.009 raw reach (≈ +0.001 … +0.002 overall)
  even if delivered perfectly, and every CI spans 0. That is below the ±0.0005 board read only marginally, and below any
  slot's worth.
- The oracle row shows that reach IS order-limited given perfect sign knowledge (+0.15 … +0.20 raw ≈ +0.03 overall).
  No public signal knows which of our calls are right well enough to capture it. This is consistent with 513's
  "sign-limited" reading: order only pays through sign knowledge.
