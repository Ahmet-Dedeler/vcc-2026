# 613 — The block-prior instrument retrodicts 458d (PASS), and identity's null-block reading restricts 639's arm to TRANSFERS

2026-10-02 13:11–13:20Z, council session [21dd4c], handed over by [d4cd5f]. Analysis only: no build, no `vcc`, no
`cell-eval2`. Peak RSS 5.0 GB; 21 s. The rules come from `experiments/639-identity-axis-lolo/PREREGISTRATION-BLOCK-DRAFT.md`
(commit 187ee529). `block613.py` was committed (2ac94af2) BEFORE its first run. Outputs: `work/block613.out`,
`work/block613.json`. The descriptive refinement (not used by any rule) is in `block613b.py` / `work/block613b.out`.

## Config
- Instrument: 636's LOLO surfaces (h1, HCT116, HEK293T, K562). CLEAN4 with any bank forbidden on a surface dropped
  (K562: Nadig Jurkat + HepG2 only).
- Block = decomp639's call set minus the head: lib521's real call set on h1; ds612c's construction on the folds.
- Gate genes: q < 0.05. Null genes: q ≥ 0.05 with lfc ≠ 0.
- Only paired cells count (both the source and CLEAN4 carry a sign).
- Primary statistic: the pooled paired sign-accuracy difference vs CLEAN4, which is 639d's "gain".
- KA: identity λ 0.5's gate-block gains reproduce `decomp639.json` exactly on all four surfaces. PASS.

## Results (pooled Δ vs CLEAN4, [95 % target bootstrap])
| surface | identity, gate block | identity, NULL block | 458 H1 programme, gate block | 458, NULL block |
|---|---:|---:|---:|---:|
| h1 (458 = its own line's programme) | +0.0057 [+0.0004, +0.0111] | +0.0014 [−0.0003, +0.0029] | +0.0874 | +0.0279 |
| HCT116 | +0.0399 [+0.0310, +0.0496] | **−0.0137** [−0.0149, −0.0125] | **−0.1695** | **−0.0263** |
| HEK293T | +0.0333 [+0.0279, +0.0391] | **−0.0081** [−0.0097, −0.0066] | **−0.1941** | **−0.0378** |
| K562 (CLEAN4 without K562's banks) | −0.0289 [−0.0580, −0.0010] | +0.0001 [−0.0010, +0.0012] | **−0.0853** | **−0.0082** |

- **(A) 458d retrodiction: PASS.**
  - 458's programme is below CLEAN4 on BOTH readings on all 3 cancer folds; the rule needed 2 of 3.
  - The instrument sees what moved the board in 458d (fid −0.0055, reach −0.0267).
  - It also gets the size about right. At full dose on the applied cells (613b), the fid estimate is −0.004 … −0.034
    raw, i.e. about −0.002 … −0.015 at 458d's 0.43 dose, which brackets the board's −0.0055.
  - The same programme is +0.087 / +0.028 on its own line (H1): a line's mean-response programme does not transfer.
- **(B) Null-block reading: identity is below CLEAN4 on 2 of 3 cancer folds (HCT116, HEK293T), so per the prereg the arm
  ships ONLY on a board TRANSFERS read (≥ +0.0015 in Δ(fid + reach + nmae)).**
  - The identity gains on gate genes (+0.033 … +0.040), but the block's ~1,400–1,500 null genes per target outnumber its
    40–140 gate genes.
- **H1 control halves (635's construction; not used by a rule):** null-block identity − CLEAN4 is +0.0013 on the full
  pool, against +0.0012 from half A versus the truth from half B. The shared-control agreement is negligible here.
- **613b (descriptive): fid exposure on the builder's APPLIED cells** (`clean4q80`: block genes with |E0| ≤ the
  target's block 80th percentile, ~1,760 of 2,200), Δfid ≈ Δcorrect / 2,600:
  - identity: h1 +0.0009, HCT116 **−0.0060**, HEK293T −0.0004, K562 −0.0002 raw (HCT116 −0.0033 overall);
  - 458: −0.019 / −0.034 / −0.004 raw on the cancer folds.
  - This is an upper bound at full dose. v3's confidence map and the builder's shrink/τ leave some applied cells
    unflipped.

## One honest sentence
The restricted-LOLO instrument earns its validation: it retrodicts the only board pair of a block-prior source change, in
sign and roughly in size. Read with it, the identity prior is a two-sided bet. It is right more often than CLEAN4 on
gate genes, but wrong more often on the far more numerous null block genes on two of three cancer lines. So its fid
exposure (≈ −0.003 overall on HCT116-like lines) can cancel the gate gain, and only a clear board TRANSFERS read should
ship it.

## 613c: block sign-source screen for null-block precision (fid's filler); `block613c.py` committed before running (64a19596)
- **Retrodictions: R1 (458d) PASS. R2 (392) PASS.** CLEAN4 beats the key's own sign on the null block on all 3 cancer
  folds (+0.036 / +0.049 / +0.004; gate block +0.147 / +0.155 / −0.048). That matches 392's board fid +0.0082 raw. The
  instrument now retrodicts two board points of block sign changes.
- **No lead** (the bar was > +0.005 on the null block on all 3 cancer folds, with the gate block ≥ −0.01). CLEAN4 is the
  best tested public block sign source for fid's filler on held-out cancer lines. The pack's +0.014 fid edge is not a
  better public block sign source of this kind.
- The one mixed signal, descriptive only:
  - Nadig Jurkat ALONE beats 4-bank CLEAN4 on both Orion lines (null +0.012 / +0.033, gate +0.007 / +0.052).
  - It loses on H1 (−0.0065 / −0.017) and on K562 (−0.0035, where CLEAN4 is already Jurkat + HepG2).
  - The K562 banks alone are much worse off-line (k562_gwps null −0.036 / −0.053, gate −0.089 / −0.147 on Orion). So
    they act like 458's H1 programme: a line-specific programme that dilutes CLEAN4 away from its own line.
  - Not actionable from these folds: the two Orion lines share a platform and agree for that reason too, and H1 and K562
    disagree. A down-weighted K562 CLEAN4 would need its own prereg and the Broad-16 tail.
- Try count: 9 candidates × 4 surfaces (2 of them retrodiction rows).

## 613d / 613d2: a hook-class instrument (uncalled genes, 426's add rule), with a registered 496a prediction
- `block613d.py` was committed (2c254e95) before its first run and before 496a's 10-05 board read.
- The model: a hook flip on an uncalled gate gene moves |error| by ±2·add, and a gene the key left unsigned moves it by
  ±add (add = fold480.mde_add, add < |lfc|). The prediction is the mean over 521-pool targets of Σ Δ|err| / Σ|T|. There
  are no base magnitudes, which is where ds612c failed.
- **Retrodiction on H1's real emitter (raw nmae, f 1 / 0.5 / 0.35):**
  | contrast | 613d | real emitter | |
  |---|---|---|---|
  | c4h − base | −0.0011 / −0.0025 / −0.0019 | −0.0015 / −0.0022 / −0.0021 | pre-registered gate, PASS |
  | c4kg − c4h (613d2) | +0.0001 / +0.0014 / +0.0012 | +0.0003 / +0.0012 / +0.0011 | known answer |
  - Its largest error on a known contrast is 0.0004.
- **Registered prediction for 496a − 446bj2 (board, 10-05): Δnmae +0.00017 raw, i.e. ≈ 0 to slightly worse.**
  - Cancer-fold mean at f 0.35: HCT116 +0.00001, HEK293T −0.00076, K562 +0.00127. H1 itself predicts −0.0019 here and
    −0.0021 on the real emitter.
  - Read rule:
    - a board move > +0.0005 confirms the instrument for hooks on held-out lines;
    - a move < −0.0005 closes it for hooks;
    - |Δ| ≤ 0.0005 is uninformative.
- **Depth sensitivity (613d3, descriptive; the registered prediction is unchanged).**
  - With the add computed at A/B/C's depth (library 20,109) instead of each fold's own 11–14k, the cancer-fold mean at
    f 0.35 is +0.00029 (registered: +0.00017).
  - HCT116 −0.0002 and HEK293T −0.0005 favour the hook; K562 +0.0016 opposes it.
  - So the "≈ 0" prediction is a balance between two CLEAN4-friendly lines and the CLEAN4-weakened K562 fold.
- **KGID (v3) − c4h by this instrument:**
  - h1 −0.0020 / −0.0003 / +0.0001;
  - HCT116 +0.0002 / −0.0009 / −0.0005;
  - HEK293T −0.0012 / −0.0015 / −0.0012;
  - K562 −0.0007 / −0.0020 / −0.0018.
  - The h1 row is tested out of sample by idkg3 under KGID prereg amendment A3 (2c24e322, ±0.0006).

## 613e: held-out premise check for 630 (rank-neutral UP magnitude on uncalled UP genes; `block613e.py`, committed first)
Uncalled gate genes whose delivered sign is UP (the key's sign; with c4h's beside it):
| surface | real-UP share, f 1 | f 0.35 | cells per target, f 1 / f 0.35 | median \|T\|/add (share > 2·add) |
|---|---|---|---|---|
| h1 | 0.535 (c4h 0.544) | 0.576 (0.582) | 189 / 56 | 3.4 (0.94) |
| HCT116 | 0.699 (0.735) | 0.701 (0.722) | 8 / 1 | 5.6 (1.00) |
| HEK293T | 0.730 (0.766) | 0.786 (0.828) | 38 / 7 | 4.7 (1.00) |
| K562 | 0.726 (0.719) | 0.889 (0.833) | 1 / 0 | 7.4 (1.00) |
- **The premise holds off H1.** UP-signed uncalled gate genes are right 70–83 % of the time on cancer lines, with ample
  headroom.
- ~~The footprint is small off H1: 1–38 such cells per target, against H1's 189.~~ **RETRACTED (13:45Z, [d4cd5f]'s
  point; 613e2).** Those counts required add > 0. The Orion folds' low depth (library 12–14k against ~20k on H1 and
  A/B/C) zeroes the add on many genes, and cell counts conflate pool size.
- **613e2, in nmae's own units:** the share of each target's Σ|T| on uncalled UP-signed gate genes:
  | surface | f 1 | f 0.35 |
  |---|---|---|
  | H1 | 0.25 | 0.23 |
  | HCT116 | 0.185 | 0.10 |
  | HEK293T | 0.184 | 0.11 |
  | K562 | 0.134 | 0.08 |
  - That is a 1.3–1.9× smaller footprint at f 1 and ~2–3× at f 0.35, not 5–25×.
  - Part of the gap is the folds' head capturing more of the gate: on each Orion fold the consensus key includes the
    other Orion line, a same-platform donor.
  - So 630's expectation for A/B/C is H1-like to ~2× smaller, and its board read settles it.

## 613f / 613g: CLEAN4's bank composition, read on 12 independent Broad lines plus 636's surfaces (both committed before running)
- **613f (single banks and subsets): CLEAN4's equal-weight 4 banks beat every composition.** Δ vs CLEAN4 on the
  block, 12 working Broad lines (the 4 reserved lines never read):
  | composition | null block: wins / 12, mean | gate block: wins / 12, mean |
  |---|---|---|
  | k562_gwps alone | 2/12, −0.0061 | 1/12, −0.0560 |
  | k562_essential alone | 0/12, −0.0061 | 0/12, −0.0319 |
  | nadig_jurkat alone | 1/12, −0.0083 | 0/12, −0.0287 |
  | nadig_hepg2 alone | 4/12, −0.0020 | 0/12, −0.0270 |
  | Jurkat + HepG2 | 4/12, −0.0011 | 4/12, −0.0048 |
  | CLEAN4 without k562_gwps | 4/12, −0.0005 | 2/12, −0.0021 |
  - 613c's Orion-only signal (Jurkat alone better) does not generalise: it wins 1 of 12 lines here.
- **613g (more banks): the ensemble hypothesis fails its bar.**
  - Adding every legal public bank (+ h1_2025 + the two Orion banks, 7 on Broad): null 6/12 (−0.0002), gate −0.0074,
    and HEK293T null −0.007 / gate −0.040. H1's programme hurts again, as in 458d.
  - Adding only the Orion banks: null 8/12 (+0.0013), gate −0.0067, K562 gate −0.020.
- **Reading.**
  - Single-line programmes do not transfer.
  - An equal average of 4 diverse CRISPRi lines does better than any subset, and adding more lines does not help it.
  - **CLEAN4 stays as October's block prior.** That is a protective negative against a tempting Orion-only change.
  - The bars are 613f's: ≥ 9/12 Broad wins on the null block, gate ≥ 0, no 636 surface below −0.005. Try count: 8
    compositions.

## 613d4 (descriptive): KGID − c4h = the gate (c4kg − c4h) + the weak-half source (KGID − c4kg), by 613d's add rule
| surface | f 1 | f 0.5 | f 0.35 |
|---|---|---|---|
| h1 | −0.0020 = gate +0.0001 + source −0.0022 | −0.0003 = +0.0014 − 0.0017 | +0.0001 = +0.0012 − 0.0012 |
| HCT116 | +0.0002 = +0.0010 − 0.0008 | −0.0009 = +0.0000 − 0.0010 | −0.0005 = +0.0001 − 0.0006 |
| HEK293T | −0.0012 = +0.0008 − 0.0019 | −0.0015 = +0.0002 − 0.0017 | −0.0012 = −0.0002 − 0.0011 |
| K562 | −0.0007 = −0.0006 − 0.0001 | −0.0020 = −0.0016 − 0.0004 | −0.0018 = −0.0013 − 0.0004 |
- The weak-half SOURCE helps on every surface. That is v3's identity sign, plus its coverage at τ 0 against c4h's
  τ 0.05.
- The key GATE costs on H1 and HCT116 and pays only on the CLEAN4-weakened K562 fold.
- So if A3 passes, KGID's held-out case rests on its weak-half source, and the gate is roughly neutral off K562. Splitting
  identity from τ-coverage would take one more arm (CLEAN4 τ 0 gated), which is not run.

## 613d5 / 613d6 (descriptive): ungated hook forms by 613d's add rule (raw nmae vs c4h; f 1 / 0.5 / 0.35)
| surface | ungated v3 identity (idall3's form) | ungated CLEAN4 τ 0 (c4all0) | KGID (v3) |
|---|---|---|---|
| h1 | −0.0016 / −0.0023 / −0.0023 | −0.0008 / −0.0008 / −0.0011 | −0.0020 / −0.0003 / +0.0001 |
| HCT116 | −0.0016 / −0.0011 / −0.0020 | −0.0019 / −0.0027 / −0.0032 | +0.0002 / −0.0009 / −0.0005 |
| HEK293T | −0.0039 / −0.0026 / −0.0011 | −0.0010 / −0.0009 / +0.0006 | −0.0012 / −0.0015 / −0.0012 |
| K562 (weakened CLEAN4) | +0.0005 / +0.0016 / +0.0024 | +0.0012 / +0.0016 / +0.0023 | −0.0007 / −0.0020 / −0.0018 |
- Ungated v3 is the best or second-best form on H1, HCT116 and HEK293T. On the K562 stress fold it is the worst.
- Its H1 row is the registered out-of-sample test (612/PREDICTION-idall3.md, 633f3a39; chain612d).
- c4all0 (CLEAN4 at τ 0, ungated) is an UNREGISTERED forking path, recorded only. Identity adds over it on H1 and
  HEK293T, not on HCT116.

## 613d7: a third known-answer row for 613d, recorded BEFORE the A3 / idall3 verdicts (14:15Z)
KGID in its v1 form (identity alone fills CLEAN4's gaps; key-gated) − c4h on h1, against the real emitter's idkg − c4h:
| | f 1 | f 0.5 | f 0.35 |
|---|---|---|---|
| 613d | −0.0039 | +0.0002 | +0.0018 |
| real emitter (3 seeds) | −0.0032 | −0.0002 | +0.0007 |
| \|error\| | 0.0008 | 0.0004 | 0.0011 |
- Signs match at f 1 and f 0.35, and f 0.5 is near zero on both. But the errors reach 0.0011, against ≤ 0.0004 on the
  two covered-gene contrasts (c4h − base, c4kg − c4h).
- 613d is less accurate on forms that sign CLEAN4-UNCOVERED genes (the v1 extension, where the key is often unsigned).
- idkg3 and idall3 sign covered genes only, so their ±0.0006 bands (A3, PREDICTION-idall3) stand as registered. A
  narrow miss should still be read with this row in mind.

## 613d OUT-OF-SAMPLE PASS #1 (14:52Z): idkg3 on H1's real emitter (KGID prereg A3)
- 613d predicted, before the run: KGID (v3) − c4h −0.00202 / −0.00034 / +0.00006.
- The real emitter gave −0.00173 / −0.00021 / +0.00005 (3 seeds). That is within 0.0003 at every regime, against
  ±0.0006 bands.
- With two covered-gene known-answer rows and one out-of-sample pass, 613d is the first validated held-out instrument
  for hook-class levers. Its cancer-fold rows now carry KGID's held-out reading, unless 496a's 10-05 read closes it.

## 613d OUT-OF-SAMPLE PASS #2 (15:21Z): idall3 (the ungated v3 identity hook) on H1's real emitter
- 613d predicted, before the mode existed (633f3a39): idall3 − c4h −0.00160 / −0.00232 / −0.00225.
- The real emitter gave −0.00186 / −0.00218 / −0.00223 (3 seeds). That is within 0.0003 again.
- A form chosen on 613d's 5-form screen held up out of sample on the real emitter. The board's 496a read (10-05) is
  613d's cross-surface exam.
