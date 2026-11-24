# 528 — DRAFT board bins for a DE-only arm: the call-propensity block (not armed; for the coordinator)

Written 2026-09-29 by the 528 worker from the H1 real-emitter runs (RESULTS §4), before any such arm is built and before any
09-30 entry exists. The coordinator fixes the slot day and may re-price against the incumbent of that day (§ Interactions).

## The arm (one variable vs 446bj2z, entry `be90smZUDcIV7j5MgQ03`, 0.2728260)
- **446bj2zcp** = 446bj2z + `patch_emitter528.py` (emitter) + `patch_builder528.py` (builder) +
  `VCC_RC_CALLSET_BLOCK=prop:experiments/528-callprop-block/work/callblock528_abc.npz:04fb24ecf81582f0af49b4e44b4dfb7b35def66425b20afd6e834c0d907bf0dd`
  (swap 1, steer `kept`: the CallSet's 2,200-gene CPM block is chosen per context by log10 control CPM + 22 · excess
  dropout (median-rescaled to H1's dropout scale), from that context's own controls; only the block's CPM part keeps the
  CLEAN4 steering). Swaps (no head): A 546, B 572, C 684 of 2,200.
- Known answers before the arm: KA-1 (unset) and KA-2 (the inert artifact `callblock528_inert_abc.npz`, e62ee765…) must
  both rebuild 446bj2z's output `ec3a5ae6…` byte for byte (RESULTS §5).
- Expected config diff vs 446bj2z: exactly `rc_callset_block` (+ `de_sign_emitter_sha256`, 499 F2's no-output-effect pin).
- The bulk is untouched (pass 1 keeps every gene total; the emitter asserts it): pds and mse must equal 446bj2z's to
  ≤ 1e-6 raw. On H1 they are equal to the last bit on every run.

## Why (H1, the production emitter, 3 emission seeds, vs 435 = 446bj2z's DE construction; RESULTS §4)
| regime | overall [95 % target CI] | Δfid | Δjac | Δreach | Δnmae (raw) | fid+jac overall | board: 0.72 × fid+jac |
|---|---|---:|---:|---:|---:|---:|---:|
| f 1.0 | **+0.0051** [+0.0033, +0.0070] (seeds +0.0052 / +0.0049 / +0.0053) | −0.0013 | +0.0091 | +0.0096 | +0.0004 | +0.0034 | **+0.0025** |
| f 0.5 | +0.0048 [+0.0017, +0.0075] | −0.0018 | +0.0089 | +0.0085 | −0.0004 | +0.0031 | +0.0022 |
| f 0.35 | +0.0045 [+0.0021, +0.0072] | −0.0018 | +0.0083 | +0.0089 | −0.0001 | +0.0028 | +0.0020 |
| cmp | +0.0053 [+0.0033, +0.0073] | −0.0018 | +0.0092 | +0.0086 | −0.0019 | +0.0032 | +0.0023 |

The mechanism is gene-level and measured, not fitted to H1's targets: at matched expression, genes with more zeros in
the controls than Poisson sampling predicts (cell-state programme genes: G2/M, immediate-early/ISR, interferon, hypoxia)
respond more (×1.4 mean |lfc|) and are called 1.3–3× as often; the block key was fitted on h1-trainval (126 other
perturbations) only; an expression-matched shuffle of the dropout feature LOSES jac (−0.0052 raw on H1, as emulated); on
public cancer truths the dropout feature beats its shuffle on 3/3 folds and the block beats the CPM block on 2/3
(K562 +0.071, HEK293T +0.039, HCT116 −0.012 capture; RESULTS §3).

## Noise model (as 523)
The swap changes the emitter's pass-2 placement draws, so reach carries the ordering channel's re-roll spread (sd ≈ 0.0045
overall between digest-identical arms). fid and jac are functions of the designed call set and the deterministic
sign/magnitude design; their re-roll spread is small (H1 seed sd of the fid + jac overall: ≤ 0.0002). **Primary = Δ(fid +
jac) in overall units vs 446bj2z** (scaled Δfid + Δjac, / 6).

## Predictions (board vs 446bj2z; fid/jac at 192's 0.72×; nmae ×1; reach unpriced, three columns)
- Δ(fid + jac) overall: **+0.0020 … +0.0025** (all four regimes).
- Δjac raw ≈ +0.0060 … +0.0066 (0.72 × +0.0083 … +0.0092); Δfid raw ≈ −0.0009 … −0.0013.
- Δnmae raw: −0.0019 … +0.0004 on H1 (≈ 0: the swapped-in genes keep today's unsteered treatment).
- Overall with reach at ×1 / ×0 / ×−0.76: +0.0037 … +0.0044 / +0.0020 … +0.0028 / +0.0007 … +0.0016.

## Bins (primary = Δ(fid + jac) overall vs 446bj2z)
- **≥ +0.0012 (≥ half the prediction): POSITIVE.** Adopt the block for A/B/C DE arms, and for October apply
  `patch_defday528.py` and run `DEF_DAY_CALLSET_BLOCK=prop` (the day's keys from D/E/F's own controls; nothing refitted).
- **+0.0004 … +0.0012: PARTIAL.** The dropout-call link transfers weakly to A/B/C. Keep it as an October rider only if the
  same entry's Δjac raw is ≥ +0.002 and Δnmae ≤ +0.001 raw; otherwise keep the CPM block.
- **−0.0004 … +0.0004: NULL.** A/B/C's truth does not reward the swap (the HCT116-fold pattern: high-CPM, low-dropout genes
  dominate the calls). Keep the CPM block; record it; no October rider.
- **≤ −0.0004: NEGATIVE.** As NULL, and the fold result for HCT116 (not K562/HEK293T) is the better model of A/B/C.
- **Falsifiers (diagnose before reading):** pds or mse moving by more than 1e-6 raw (a bulk leak: stop); Δjac ≤ 0 while
  Δfid is at or above its prediction (the call set did not move: a construction or artifact error); |Δnmae| > 0.003 raw
  (the swapped genes' delivered lfc moved more than the design implies).
- Reach and nmae are reported, not binned.

## If 446bj2zmz (527's hook) is the incumbent when this is armed: the variant 446bj2zmzcp
Same knob on 446bj2zmz's env (KA-1 then rebuilds 446bj2zmz's output). On H1 (3 seeds, vs `cpm+mzg` = 446bj2zmz's DE
construction; RESULTS §4.5) the pair is worth MORE, because the hook's programme signs reach exactly the genes the new
block calls (under steer `kept` they sit outside the steered block, i.e. in the hook's complement):

| regime | overall [95 % CI] (seeds) | Δfid | Δjac | Δreach | Δnmae | fid+jac overall | board 0.72 × fid+jac |
|---|---|---:|---:|---:|---:|---:|---:|
| f 1.0 | **+0.0085** [+0.0064, +0.0108] (+0.0083 / +0.0079 / +0.0093) | +0.0032 | +0.0089 | +0.0120 | −0.0015 | +0.0058 | **+0.0042** |
| f 0.5 | +0.0074 [+0.0042, +0.0105] | +0.0031 | +0.0088 | +0.0068 | −0.0016 | +0.0057 | +0.0041 |
| f 0.35 | +0.0075 [+0.0050, +0.0104] | +0.0031 | +0.0083 | +0.0088 | −0.0014 | +0.0054 | +0.0039 |
| cmp | +0.0082 [+0.0060, +0.0104] | +0.0031 | +0.0091 | +0.0095 | −0.0022 | +0.0058 | +0.0042 |

Bins for this variant (primary Δ(fid + jac) overall vs 446bj2zmz; predicted +0.0039 … +0.0042): POSITIVE ≥ +0.0020;
PARTIAL +0.0007 … +0.0020; NULL −0.0007 … +0.0007; NEGATIVE ≤ −0.0007; the same falsifiers (pds/mse ≤ 1e-6; Δjac > 0;
|Δnmae| ≤ 0.003 raw). Do NOT arm the coupled form (steer `key`) on the hook: on H1 it collapses to +0.0012 (nmae +0.009).

## Do not stack with 523's K 1,600 as is
At `CallSet(400, 1200)` the kept swap makes ~377 extra (unprescribed) calls per target on H1 (low-to-mid CPM complement
genes; base 66) and prices −0.0012 … −0.0027 at f ≤ 0.5 / cmp over 2 seeds (+0.0054 at f 1; RESULTS §4.5). If 523's pair reads K 1,600 positive, this
knob needs its own re-design before it rides on the smaller call set.
