# 612: board bins for the UNGATED v3 identity hook (`idall3`), read against 496a's plain-hook effect (FINAL, 15:25Z)

Drafted 2026-10-02 ~13:55Z by [21dd4c], on the coordinator's rules (NEXT-SLOT-DECISION, 10-02 ~13:50Z).
- **Status: FINAL.** idall3 PASSED its pre-registered H1 test at 15:21Z (`work/idall3verdict612.out`), within 0.0003 at
  every regime. Committed BEFORE 10-05's 496a read. Per the coordinator's rule, idall3 replaces KGID (KGID retires).
- (It existed only if idall3 PASSED its pre-registered H1 real-emitter test, 612/PREDICTION-idall3.md, 633f3a39.)
  On a FAIL this draft is void, and the hook family stays as it is.
- On a PASS it REPLACES KGID in the queue, and KGID retires.
- It must be committed in final form BEFORE 10-05's 496a read, so that nothing in it is tuned on that read.
- ⟨H1⟩ rows are filled from `work/score612.json` (idall3, 3 seeds) when chain612d lands. The bins below do not move
  when they do.

## The arm
- **base** = the best entry at slot time, plus
  `VCC_DE_MDE_SIGN=file:experiments/604-identity-prior/work/ident604v3_flex_lam0.5_ABC.npz:6cd84dc258cb609d08d93ec4c6c6ec7c7134ff8e3093eccde5ac15bbd2ac618f:0`
  (604 v3: z(P_C4 − ½) + 0.5·z(I) on CLEAN4-covered genes, τ 0). **No** `VCC_DE_MDE_SIGN_KEYGATE`.
- The builder already parses this form (`file:<npz>:<sha256>:<tau>`). October serves it as
  `DEF_DAY_MDE_SIGN=hookbank:identity:0` (604's patch, 10-06), so no new builder code is needed.
- If the base already carries a hook, the arm swaps the source.
- Expected config diff vs the base: exactly `de_mde_sign` (+ `de_sign_emitter_sha256` only if the emitter moved). pds and
  mse are bit-identical: a DE-only, prior-class arm.

## Evidence before the read (11.1c.2: two readings)
1. **H1, real production emitter.**
   - idall3 − c4h, 3 seeds (per-seed sd ≤ 0.0001), raw nmae **−0.00186 / −0.00218 / −0.00223** at f 1 / 0.5 / 0.35.
     613d had registered −0.00160 / −0.00232 / −0.00225: PASS.
   - Other members vs c4h: reach +0.0045 / +0.0029 / +0.0002; fid and jac |Δ| ≤ 0.0001. Overall-eq +0.0014 / +0.0011 /
     +0.0006.
   - idall3 − base: nmae −0.0033 / −0.0044 / −0.0044; overall +0.0015 / +0.0015 / +0.0012.
   - idall3 − idkg3 (KGID v3): nmae −0.0001 / −0.0020 / −0.0023; overall +0.0005 / +0.0005 / +0.0007. The ungated form
     beats the key-gated one at every regime on H1.
2. **Held-out folds** (613d's add-rule instrument; validated by known-answer rows, plus this H1 test), vs c4h at f 1 / 0.5
   / 0.35:
   - HCT116 −0.0016 / −0.0011 / −0.0020;
   - HEK293T −0.0039 / −0.0026 / −0.0011;
   - **K562 stress tail +0.0005 / +0.0016 / +0.0024.**
     - K562's CLEAN4 lacks its own two banks. INSTRUMENTS.md reading rule 1 treats it as a weakened-CLEAN4 stress test,
       not a typical D/E/F line.
     - It is still a real risk if a D/E/F line's biology sits far from all four CLEAN4 lines.

3. **Broad-12, a different lab and modality (Cas9 KO, 3′)** ([7772dc] 622o, 74042c4a; 613d's `pred_dnmae`, verified
   equal; reserved lines never read). With a production-like DONOR key (K562 GWPS's response to the same target),
   idall3 − c4h is:
   - worse on 8 / 5 / 7 of 12 lines at f 1 / 0.5 / 0.35, mean +0.0001 … +0.0005 raw;
   - by as much as K562's +0.0024 on 1–2 of 12 lines (T3M5 +0.0030 at f 1; ONE58 +0.0027 at f 0.35);
   - about neutral once balanced for KO's UP lean (worse on 4–6/12, mean ≈ 0).
   - Mechanism: identity's part helps on 9–12/12 lines (idall3 − c4all ≈ −0.001). The τ 0 coverage expansion hurts on
     8–10/12 (c4all − c4h ≈ +0.001). That is the reverse of HCT116 in 613d6.
   - **So the held-out evidence for idall3 over the plain hook is mixed by modality:** favourable on the CRISPRi folds,
     neutral to slightly negative on KO lines, with a K562-magnitude loss on ~10–15 % of lines. A/B/C, and so D/E/F,
     are CRISPRi Flex from one lab, so the CRISPRi folds and the board read weigh most. This is still why a TIE goes to
     the plain hook.
   - (7772dc's candidate, c4h's coverage plus identity only on |z(I)| ≥ 1.5 tails, is consistent but tiny on Broad:
     −0.0001 … −0.0002 raw. Not pursued.)
   - **2026-10-03 ~21:45Z note (bins unchanged):**
     - (a) This prereg never quoted 622n's "Broad 0.78 / 12/12" low-I → UP tail. 7772dc found that tail holds on 8/12 lines,
       +0.05 over each line's own base rate. The 622o hook rows above are base-rate-balanced and offset-free (7772dc, 657),
       so they stand.
     - (b) "Favourable on the CRISPRi folds" above referred to 613d's Orion rows, now VOID (616). The CRISPRi evidence is
       now:
       - H1, the registered real-emitter PASS;
       - 618, the held-out H1 panel with the faithful key: idall3 − c4h nmae −0.0024 / −0.0006 / −0.0007 at f 1 / 0.5 /
         0.35 on the real emitter. 661 registers a board TIE vs Δ_496a from it (PRED-661-idall3.md, 00ad7165).

**Addendum (2026-10-02 ~16:20Z, before 10-05; bins unchanged): additivity with 630** ([d4cd5f], 630 RESULTS, 3 seeds on
the idall3d dumps).
- 630's rank-neutral UP pass is worth, overall-eq at f 1 / 0.5 / 0.35 / cmp:
  - on the idall3 base: +0.0027 / +0.0029 / +0.0028 / +0.0035;
  - on the plain-hook base: +0.0023 / +0.0030 / +0.0029 / +0.0034;
  - 0 calls changed.
- So the two compose additively, within 0.0004 at every regime. If the base at slot time already carries 630 (applied
  10-06), the arm's Δ against its own base stays comparable to Δ_496a (read on 446bj2, no 630). The primary statistic is
  unaffected.

## Primary statistic, noise, prediction
- **Primary:** Δ(fid + nmae) overall, arm − base, read against Δ_496a = Δ(fid + nmae) of 496a − 446bj2 (the plain hook's
  own effect).
- Noise for a prior-class hook arm: σ ≈ 0.0001 overall, so 2σ ≈ 0.0002 (537 §2.4).
- **Prediction:** board Δnmae vs the plain hook ≈ H1's f 0.35 value (−0.00223) × 1.2–1.6 (611) = −0.0027 … −0.0036 raw
  ≈ **+0.0007 … +0.0010 overall** in Δ(fid + nmae) relative to Δ_496a. fid ≈ 0. That is well clear of the WIN line
  (+0.0002).
- **Reach band:** ±(\|H1 Δreach at f 0.35\| + 0.004) raw = ±(0.0002 + 0.004) = **±0.0042 raw**.

## Bins (written before either read)
| outcome | A/B/C entries | October |
|---|---|---|
| Δ_arm ≥ Δ_496a + 0.0002, with Δreach in its band (WIN) | idall3 | idall3, if 613d survives 496a's read (see below). Otherwise idall3 only with a written reason (11.1c.2) |
| \|Δ_arm − Δ_496a\| < 0.0002 (TIE) | the plain CLEAN4 hook | the plain CLEAN4 hook: the incumbent, already board-read. idall3 carries the K562 stress tail and a best-of-N selection |
| Δ_arm < Δ_496a − 0.0002 (BELOW) | the plain hook | the plain hook |
| Δ_arm < −0.0002 vs its own base (HARM) | no identity hook | no identity hook. Re-check the artifact's KA (λ 0 == CLEAN4) before closing the family |
| pds or mse ≠ the base by > 1e-6, or \|Δjac\| > 0.001 raw | **falsifier**: a delivery leak, stop and audit | — |

**October section (coordinator rule 2).**
- 613d registered a 496a prediction before its read: Δnmae +0.00017 raw (≈ 0).
- If 496a's board Δnmae is < −0.0005, 613d is closed for hooks. idall3's fold rows then carry no weight, and its October
  case rests on H1 + this board read only. Adopting it then requires a written reason in october-decisions.json.
- If 496a reads > +0.0005, 613d's fold rows stand as idall3's held-out reading.
- If |Δ| ≤ 0.0005, the read is uninformative about 613d. Its rows are then reported, not relied on.
- **Stress-tail probability for a D/E/F line:** ~10–15 % of lines lose K562-magnitude nmae (+0.0024 … +0.0030 raw) under
  idall3 vs the plain hook (Broad-12, KO, donor key). Even in a WIN, write that in October's decision record.
- **2026-10-02 19:55Z (bins unchanged):** K562 is INCONCLUSIVE for idall3, neither "against" nor "favourable". The
  coordinator's ruling is in docs/NEXT-SLOT-DECISION.md (10-02 19:49Z). It weighs 614's real-emitter row (idall3 − c4h
  −0.00026, descriptive) against 613d's analytic row (+0.0005 … +0.0024); both read a low-power truth. Also from 616:
  the 496a/613d discrimination rule above is withdrawn (NEXT-SLOT-DECISION 19:10Z). The Orion fold rows are void.

## Try count (11.1c.4)
- Hook forms priced on H1's real emitter today: c4h, c4kg, idkg (v1), idkg3, idall3, plus the c4kgall KA and the idkgd
  dump, which are duplicates.
- Hook forms on 613d: c4h, c4kg, KGID v3, ungated v3, ungated CLEAN4 τ 0 (the last unregistered, recorded only).
- On the void fold models: kgh612b/c/d (4 forms), qscan612 (8 q values), ds612 (13 variants).
- That is ≥ 8 hook forms looked at. idall3 was CHOSEN from 613d's 5-form screen (best-of-5 on that instrument). Its H1
  real-emitter run is a single pre-registered test of that choice, not a best-of-N on H1.
- Board reads for this switch: this one, plus 496a's plain-hook read as the comparison.
