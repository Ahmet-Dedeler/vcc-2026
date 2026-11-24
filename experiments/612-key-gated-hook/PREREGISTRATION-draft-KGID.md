# 612 — [RETIRED 15:25Z: replaced by idall3, the ungated v3 hook, per the coordinator's rule after idall3 passed its registered H1 test; see PREREGISTRATION-draft-IDALL3.md] DRAFT board bins for the KEY-GATED IDENTITY hook (10-08+), read against 496a's plain-hook effect (for the coordinator; not armed)

Written 2026-10-02 ~12:30Z by [21dd4c] (council 600), BEFORE any 10-05 or later entry exists and before the arm is built.
The H1 real-emitter rows are filled from `work/score612.json` (3 seeds). **Amendment A1** (committed 3efebfd7, 12:52Z,
before ds612c existed) added a retrodiction gate for the fold instrument. **ds612c FAILED it (12:53Z), so A1's FAIL
branch is the decision rule** (see "A1 result" below): a TIE goes to the plain CLEAN4 hook, and KGID needs the WIN bin. Patch: 612 (`VCC_DE_MDE_SIGN_KEYGATE`), applied 10-06 after 598. Identity artifact: 604
([4d0207]).

## The arm
- **446…kgid** = the best entry on 10-07 (the "base"), plus:
  - **[A2, 13:30Z]** `VCC_DE_MDE_SIGN=file:experiments/604-identity-prior/work/ident604v3_flex_lam0.5_ABC.npz:6cd84dc258cb609d08d93ec4c6c6ec7c7134ff8e3093eccde5ac15bbd2ac618f:0`
    (604 v3: score = z(P_C4 − ½) + 0.5·z(I) on CLEAN4-COVERED genes only, side = sign(score), Flex references;
    `pdown_<c>` finite on 6,346 genes per context; KAs pass per [4d0207]). **τ 0** = the v3 sign on every covered gene.
  - (Superseded by A2: the v1 artifact `ident604_flex_lam0.5_ABC.npz` carried the identity-only extension, 10,962
    genes. That is the extension 639d measured wrong on UNCALLED genes of cancer lines, HEK293T 0.27–0.34 and HCT116
    0.49, and the hook acts only on uncalled genes. Found by [d4cd5f]. The H1 rows priced it in idsign612's 622 form,
    WITH the extension; the v3-form re-run `idkg3`, seeds 0–2, replaces them below.)
  - `VCC_DE_MDE_SIGN_KEYGATE=0.5` (612: the hook's sign only on each target's weak-key half; the strong half keeps its
    delivered key sign).
  - If the base already carries a hook, the arm swaps its source and adds the gate. The config diff then names exactly
    those keys.
- **Expected config diff vs the base:** `de_mde_sign` (+ `de_mde_sign_keygate`), plus `de_sign_emitter_sha256` only if the
  emitter moved. pds and mse equal the base's to ≤ 1e-6 raw: a DE-only, prior-class arm.
- **Not** the identity BLOCK prior (VCC_SAMELINE with programme on). That is a separate switch with its own read.

## Why: evidence before the read
1. **[SUPERSEDED by A1: an unvalidated instrument, kept for the record]** **Held-out folds, "delivery-faithful"**
   (612/kgh612c, kgh612d: 480's `base + s·add` shape; τ 0 for the identity and key-gated forms). This model applies the
   hook to every gene, including production's 2,600 CALLED genes, which the hook never touches. With the call set
   masked (ds612c), kgID vs no hook falls from +0.021 to +0.0013 (HCT116) and the K562 "tail" from −0.010 to ≈ −0.0015. Δnmae vs no hook, positive = better:
   - plain CLEAN4 (τ 0.05): h1 −0.0006, HCT116 +0.0149, HEK293T +0.0148, K562 −0.0102;
   - key-gated CLEAN4 (τ 0): +0.0010, +0.0179, +0.0163, +0.0016;
   - **key-gated identity (τ 0): +0.0057, +0.0210, +0.0203, +0.0003.**
   - Key-gated identity − key-gated CLEAN4, at matched coverage: +0.0048 / +0.0031 / +0.0040 / −0.0013 (ns).
   - The plain hook's tail is −0.010 on K562, the line whose CLEAN4 lacks its own banks. Both key-gated forms avoid it.
2. **H1, real production emitter** (612/run612 via 537's harness; 521's regimes; score612; same base, bulk sums bit-equal).
   Raw nmae Δ (negative = the first form better), mean over seeds (per-seed sd), then overall-equivalent:
   | contrast | f 1 | f 0.5 | f 0.35 | cmp |
   |---|---:|---:|---:|---:|
   | idkg − c4h (KGID − plain hook), 3 seeds | −0.00315 (0.00026) | −0.00019 (0.00018) | **+0.00067** (0.00019) | −0.00365 (0.00026) |
   | overall-eq, same contrast | −0.00136 | +0.00080 | +0.00058 | +0.00109 |
   | c4kg − c4h (the key gate on CLEAN4), 2 seeds | +0.00034 | +0.00126 | +0.00108 | +0.00011 |
   | idkg − c4kg (identity vs CLEAN4 on the weak half), 2 seeds | −0.00333 | −0.00136 | −0.00033 | −0.00361 |
   | idkg − base (KGID vs no hook), 3 seeds | −0.00461 | −0.00238 | −0.00146 | −0.00437 |
   - fid and jac move by at most 0.0001 in every row, so this is a delivery-class change.
   - reach, idkg − c4h: −0.0121 / +0.0037 / +0.0042 / +0.0006.
   - c4kg s2 is pending. Its two seeds agree to a per-seed sd of ≤ 0.00004.
   Diagnosis (612/diag612b), two separate effects:
   - (a) the identity source works: weak-key uncalled genes outside CLEAN4's coverage go from 0.507 to 0.541 sign
     accuracy, −0.003 nmae at f 1;
   - (b) the key gate gives up the plain hook's SHRINKAGE: on strong-key CLEAN4-covered genes the opposite add shrinks
     |lfc| where two sources disagree, worth +0.0011 nmae against idkg at f 0.35, although the key's sign there is the
     better one (0.657 vs 0.646).
   H1 is tuned-on and CLEAN4-friendly (its fold row: plain ≈ key-gated ≈ 0). The folds are the held-out read. The October
   switch needs both (AGENTS 11.1c.2).
3. **Mechanism (612 diag612).** The hook sets only the direction of a small MDE add on a key-signed base. It flips 2–6 % of
   signs and otherwise shrinks |lfc| where the source disagrees with the key. The gate stops it touching strong-key genes,
   whose key sign is the more accurate one. The identity source raises sign accuracy on the weak half (639/622: gate-gene
   sign +0.009 … +0.042 over CLEAN4 leave-one-line-out).

## Primary statistic, noise and predictions
- **Primary:** Δ(fid + nmae), overall units, arm − base. The DE noise for a hook-class knob is 0.3–0.6× a re-roll
  (537 §2.4): σ ≈ 0.0001 overall, 2σ ≈ 0.0002. reach is reported, never primary.
- **Board prediction (611 multipliers: nmae × 1.2–1.6 for delivery levers):**
  - KGID vs no hook (idkg − base) at f 0.35: −0.00146 raw × 1.2–1.6 = −0.0018 … −0.0023 raw nmae ≈ +0.0005 … +0.0006
    overall.
  - KGID vs the plain hook (the 496a comparison):
    - 611's nmae ratios hold at f 1 / 0.5 / 0.35 alike, but this contrast turns sign across them.
    - At the board-like f 0.35–0.5, H1 says +0.0007 … −0.0002 raw, i.e. +0.0011 … −0.0003 raw on the board.
    - In the primary Δ(fid + nmae) that is −0.0003 … +0.0001 overall: **the TIE bin is the most likely read**.
  - The folds' expectation (kgh612c, Orion mean +0.021 raw vs no hook: −0.004 … −0.008 raw on A/B/C) is SUSPENDED
    until A1's retrodiction passes. The fold model over-states KGID − plain on h1 (−0.0064 / −0.0034 / −0.0041 at
    f 1 / 0.5 / 0.35, ds612b, vs the real emitter's −0.0032 / −0.0002 / +0.0007).
- **Reach band (611 rule for delivery arms):** ±(H1 |Δreach| at f 0.35 + 0.004) raw = ±(0.0041 + 0.004) = **±0.0081 raw**
  for the board's Δreach (arm − base). H1 idkg − base Δreach at f 0.35 is +0.0041 (3 seeds).

## The comparison: 496a's plain-hook effect (10-05)
Δ_496a = Δ(fid + nmae) overall of 496a − 446bj2. 496a is the plain CLEAN4 hook on 446bj2's construction. Both arms are
prior-class, so each reads exactly against its own base.

## Bins (written before either read)
| outcome (Δ_arm and Δ_496a in overall units, primary statistic) | decision for October's hook |
|---|---|
| Δ_arm ≥ Δ_496a + 0.0002 **and** Δreach within the band | **KGID**: identity source + key gate (`DEF_DAY_MDE_SIGN` from the day-of identity artifact + KEYGATE 0.5) |
| \|Δ_arm − Δ_496a\| < 0.0002 (tie) **and** Δ_arm ≥ 0 | **KGID**, on MIN-MAX REGRET alone: equal on A/B/C, and the folds put the plain hook's worst case at −0.010 nmae (K562) against the identity hook's −0.0013 (ns). It is an expected-value-vs-tail choice, not a dominance (H1 seed 0 has the plain hook ahead at f 0.35) |
| Δ_arm < Δ_496a − 0.0002 **and** Δ_496a > 0 | **plain hook on A/B/C entries only**. October still takes KGID unless the board gap exceeds 0.0010: the held-out folds favour KGID by +0.0055 … +0.0105 raw nmae on every line, and A/B/C are three lines, not D/E/F. Write the reason down |
| Δ_arm < −0.0002 (harm vs its own base) | **no identity hook**. Re-check the artifact (604's KA λ 0 == CLEAN4) before closing the family |
| pds or mse ≠ the base by > 1e-6, or \|Δjac\| > 0.001 raw | **falsifier**: a delivery leak, stop and audit |

## Amendment A1 (2026-10-02 ~12:55Z, at the coordinator's request; written and committed BEFORE any ds612c row exists)

Why: the fold model behind "Why" item 1 (kgh612c) applies the hook to every gene. Production adds only on UNCALLED
genes (`rc_emit_k.py`: `use_s = ~head & ~blk & (asg != 0)`). On h1 that model gets the KGID − plain contrast wrong
(ds612b: −0.0064 / −0.0034 / −0.0041 at f 1 / 0.5 / 0.35, against the real emitter's −0.0032 / −0.0002 / +0.0007). So
its Orion and K562 rows, which carry the min-max-regret basis (the K562 tail), are not validated (AGENTS 11.1 rule 3).

**ds612c, specified before it runs:**
- The fold model of ds612b, with production's call set masked out of every hook form:
  - head = each target's top 400 genes by |donor consensus key|;
  - block = the top 2,200 of the rest by the surface's control expression (lib521's N_HEAD 400 / N_CPM 2,200);
  - on head ∪ block, every form takes the no-hook prediction, so those cells cancel in every contrast.
- Folds: h1, HCT116, HEK293T, K562. Regimes: 521's pools at f 1 / 0.5 / 0.35 (ds612b's `pool_T`).
- Forms:
  - nohook;
  - plain (CLEAN4 τ 0.05) = c4h;
  - kgC4_05 (the gate on CLEAN4 τ 0.05) = production c4kg;
  - kgID = idkg;
  - IDw_C4s.
- Known answer: with an empty call set, ds612c must reproduce ds612b's rows exactly, or ds612c is void.

**(1) Retrodiction criterion: ds612c passes only if its h1 cell meets ALL of these.**
- (a) kgID − plain at f 1 is negative, within ×0.5–2 of the real emitter's −0.00315 raw: in [−0.0063, −0.0016].
- (b) kgID − plain at f 0.35 is positive (point estimate > 0; real +0.00067, 3 seeds, about 6 SE from 0).
- (c) kgID − plain at f 0.5 is "−/0": point estimate ≤ +0.0002 (the real mean −0.00019 plus 2 per-seed sd), and between
  ds612c's own f 1 and f 0.35 values (monotone in f, as the real emitter is).
- (d) The mechanism: kgC4_05 − plain is positive at f 0.5 and at f 0.35, i.e. the gate costs nmae as it does on the
  real emitter (+0.00126 / +0.00108). The held-out decision rests on exactly this strong-half shrink (it IS the K562
  tail), so a model that is right on (a)–(c) for another reason does not pass. f 1 (+0.00034) is not required.
- The bands are fixed here; they are not re-set after the rows exist.

**(2) Decision rules after ds612c.** They replace the TIE and BELOW rows of the bins above. WIN, HARM and the
falsifier rows do not change.

If ds612c FAILS (no validated held-out instrument):
- The min-max-regret basis is unsupported.
- TIE (|Δ_arm − Δ_496a| < 0.0002) → the **plain CLEAN4 hook**: the incumbent, already board-read, simpler on the day.
- BELOW (Δ_arm < Δ_496a − 0.0002) → the plain hook, for A/B/C AND October.
- KGID needs the WIN bin: Δ_arm ≥ Δ_496a + 0.0002 with Δreach inside ±0.0081 raw.

If ds612c PASSES, its Orion and K562 rows replace kgh612c's in "Why" item 1, and:
- TIE → the form with the smaller **max regret over the 4 folds × {f 0.5, f 0.35}**. Here regret(form, fold, f) = the
  better of {plain, KGID}'s Δnmae vs no hook − this form's Δnmae, which is ≥ 0. Equal max regret → plain (the
  incumbent).
- BELOW → the plain hook on A/B/C entries. October takes KGID only if ds612c has KGID ahead of plain on all 4 folds at
  both f 0.5 and f 0.35 AND the board gap is ≤ 0.0010. Otherwise October takes plain.

## Amendment A2 (2026-10-02 13:30Z, before any board read): the artifact is 604 v3, and the H1 rows are re-priced
- The arm named the v1 artifact by mistake (see "The arm"). The fix changes no bin and no rule.
- The H1 real-emitter rows in "Why" item 2 (idkg) used the identity-only extension. `idkg3` (612/idsign612v3.py: 6,161
  signed genes on H1 = CLEAN4's coverage; the 4,616 extension genes dropped; equal to idkg's signs elsewhere) re-prices
  KGID in v3's form, seeds 0–2 (chain612c).
- The board prediction and the expected TIE bin are re-derived from the idkg3 rows when they land. Until then, read
  the idkg numbers as an upper bound on the extension's effect on H1, where the extension's signs were ~0.55 accurate.
  On cancer lines 639d measured them as harmful.

## Amendment A3 (2026-10-02 ~13:32Z): a second hook instrument (613d), tested out of sample on idkg3 BEFORE idkg3 is scored
Committed while idkg3's seeds were still running (chain612c). No idkg3 contrast has been computed or read.
- **The instrument:** 613/block613d.py, committed before its first run (2c254e95).
  - It applies 426's add rule to UNCALLED gate genes on 521's pools. A hook flip moves |error| by ±2·add, and a gene
    the key left unsigned by ±add. There are no base magnitudes, which is what sank ds612c.
  - Known-answer rows on H1's real emitter (computed after the real values were known), raw nmae at f 1 / 0.5 / 0.35:
    | contrast | 613d | real emitter |
    |---|---|---|
    | c4h − base | −0.0011 / −0.0025 / −0.0019 | −0.0015 / −0.0022 / −0.0021 |
    | c4kg − c4h (613d2) | +0.0001 / +0.0014 / +0.0012 | +0.0003 / +0.0012 / +0.0011 |
  - So its largest error on a known contrast is 0.0004.
- **The out-of-sample test, fixed now.**
  - 613d's h1 KGID (v3) − c4h is −0.00202 / −0.00034 / +0.00006 at f 1 / 0.5 / 0.35.
  - The real emitter's idkg3 − c4h (3 seeds, score612) must lie within ±0.0006 of each value (the largest known error
    0.0004, plus 2 SE of a 3-seed mean, about 0.0002):
    - f 1 ∈ [−0.0026, −0.0014];
    - f 0.5 ∈ [−0.0009, +0.0003];
    - f 0.35 ∈ [−0.0005, +0.0007].
  - PASS iff all three hold.
- **If A3 PASSES:**
  - 613d becomes KGID's held-out instrument. Its cancer-fold rows replace ds612c's.
  - KGID − c4h, raw nmae at f 1 / 0.5 / 0.35:
    - HCT116 +0.0002 / −0.0009 / −0.0005;
    - HEK293T −0.0012 / −0.0015 / −0.0012;
    - K562 −0.0007 / −0.0020 / −0.0018.
    - So KGID leads on all 3 cancer folds at f 0.5 and f 0.35.
  - A1's PASS-branch rules then apply with these rows:
    - TIE → the smaller max regret over {h1 + the 3 cancer folds} × {f 0.5, f 0.35}, with equal → plain;
    - BELOW → plain on A/B/C, and October takes KGID only if it leads on all folds at both regimes and the board gap is
      ≤ 0.0010.
- **The PASS branch is VOID if 613d's own 10-05 test closes it for hooks.**
  - 613d registered, before 496a's read, a prediction of +0.00017 raw nmae for 496a − 446bj2, i.e. ≈ 0.
  - If the board instead moves < −0.0005 raw, the cancer-fold rows do not represent A/B/C for hooks, and A1's FAIL
    branch returns.
- **If A3 FAILS:** A1's FAIL branch stays: TIE → plain, and KGID needs WIN.

## A3 result (2026-10-02 14:52Z): PASS (`a3verdict612.py`, `work/a3verdict612.out`)
- The real emitter's idkg3 − c4h raw nmae (3 seeds), against 613d's registered values:
  | | f 1 | f 0.5 | f 0.35 |
  |---|---|---|---|
  | real emitter | −0.00173 | −0.00021 | +0.00005 |
  | 613d (registered) | −0.00202 | −0.00034 | +0.00006 |
  - Every regime is within 0.0003, inside the ±0.0006 bands.
- **613d is KGID's held-out instrument**, unless 496a's 10-05 read closes it for hooks (Δnmae < −0.0005).
- The full idkg3 rows on H1:
  - vs c4h: reach +0.0017 / +0.0029 / +0.0001; overall +0.0009 / +0.0006 / −0.0001.
  - vs base: nmae −0.0032 / −0.0024 / −0.0021; overall +0.0011 / +0.0009 / +0.0006.
  - v3 vs v1 (idkg3 − idkg): reach at f 1 is **+0.0138**. Dropping the identity-only extension recovered the −0.012
    reach loss that idkg showed at f 1.
- **Decision rules now in force (A1's PASS branch, with 613d's rows):**
  - TIE → the smaller max regret over {h1, HCT116, HEK293T, K562} × {f 0.5, f 0.35}: KGID's 0.00006 (h1 at f 0.35)
    against the plain hook's 0.00204 (K562 at f 0.5). **So a TIE goes to KGID.**
  - BELOW → the plain hook on A/B/C. October would need KGID to lead on all folds at both regimes, and it trails by
    0.00006 on h1 at f 0.35. So October also takes the plain hook.
  - WIN → KGID.
- **Coordinator rule (NEXT-SLOT-DECISION, 13:50Z):** if idall3 (the ungated v3 hook) passes ITS registered H1 test, it
  REPLACES KGID in the queue and KGID retires. That test is running (chain612d2).

## A1 result (ds612c, 2026-10-02 12:53Z; `work/ds612c.out`, `work/ds612c.json`)
- Known answer: an empty call set reproduces ds612b on 4 folds × 3 regimes × 3 forms exactly. PASS.
- Retrodiction on h1 (ds612c vs the real emitter), raw nmae, negative = the first form better:
  | contrast | f 1 | f 0.5 | f 0.35 |
  |---|---:|---:|---:|
  | kgID − plain, ds612c | −0.00504 | −0.00188 | −0.00268 |
  | kgID − plain, real emitter (3 seeds) | −0.00315 | −0.00019 | +0.00067 |
  | kgC4_05 − plain, ds612c | +0.00048 | +0.00023 | +0.00051 |
  | kgC4_05 − plain, real emitter (2 seeds) | +0.00034 | +0.00126 | +0.00108 |
  - (a) PASS;
  - (b) **FAIL**: the f 0.35 sign is wrong;
  - (c) **FAIL**: not monotone in f;
  - (d) PASS: the gate's cost has the right sign, at 20–50 % of the real size.
- **Verdict: FAIL.** No fold instrument is validated for hook-class decisions. **The decision rules in force are A1's
  FAIL branch:**
  - TIE → plain CLEAN4 hook;
  - BELOW → plain hook (A/B/C and October);
  - KGID needs WIN: Δ_arm ≥ Δ_496a + 0.0002 in Δ(fid + nmae), with Δreach inside ±0.0081 raw.
- Expected read (H1 f 0.35–0.5): −0.0003 … +0.0001 overall in Δ(fid + nmae), i.e. most likely a TIE, so the plain hook
  stays. P(WIN) is low: the H1 point predictions sit below the +0.0002 WIN line. The reach gain (+0.0037 … +0.0042 raw,
  about +0.0008 overall) is not in the primary statistic.
- For the record only (an unvalidated instrument; not used for any decision), ds612c's held-out rows, Δnmae vs no hook at
  f 0.35:
  - plain: HCT116 −0.0000, HEK293T +0.0008, K562 −0.0013;
  - KGID: −0.0000, +0.0016, +0.0006.
  - Most of the old fold gain sat on weak-key genes of the CPM block. Production steers those with the block prior, not
    with the hook.

## Try count (11.1c rule 4)
- Gate: on/off (1 knob). Sign sources compared: plain CLEAN4, key-gated CLEAN4, key-gated identity (3).
- q: 0.5 pre-specified; an 8-value scan after (qscan612: flat 0.5–0.8), NOT used to select.
- λ: {0.5, 1} on 622's folds; 0.5 chosen on nmae (forking path, noted by 7772dc).
- References: flex / mixed / 3′-only (622); flex chosen.
- Earlier zero-base models (dosecurve610: 3 sources × 6 doses; kgh612): superseded and not used for this decision.
- Board reads for this switch: **this one** (plus 496a's plain-hook read as the comparison).
