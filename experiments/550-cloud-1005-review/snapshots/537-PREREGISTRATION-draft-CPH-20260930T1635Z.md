# 537 — DRAFT board bins for 446bj2zcph: the call-propensity block WITH the CLEAN4 uncalled-gene hook (not armed; for the coordinator)

Written 2026-09-30 by the 537 worker from real production-emitter H1 runs (RESULTS §3.2, P1), before any 10-01 entry exists
and before the arm is built. The coordinator fixes the slot day (the 10-07 rule; `rule1007_draft537.py` is a draft of it) and
adopts or amends these bins before that day's first entry.

## The arm
- **446bj2zcph** = 446bj2z (entry `be90smZUDcIV7j5MgQ03`, 0.2728260) + `VCC_RC_CALLSET_BLOCK=prop:experiments/528-callprop-block/work/callblock528_abc.npz:04fb24ecf81582f0af49b4e44b4dfb7b35def66425b20afd6e834c0d907bf0dd`
  (528's block, swap 1, steer `kept`: exactly 446bj2zcp's knob) + `VCC_DE_MDE_SIGN=clean4:0.05` (exactly 496a's knob).
- Both knobs are in the committed builder 4daae237 / emitter 0f56209d and compose in one emitter call (the builder passes
  `block_key/block_steer` and `lfc_target["mde_add_sign"]` together; `make_submission_pm_cd4.py` near its `rc_hybrid_block_k`
  call). No new code, so no new known answer is needed beyond the ones each knob already passed (535's KA535a/b/c; 527's KA-1 /
  KA-496).
- Build: `run_arm_generic.sh 446bj2zcph 1.75 <446bj2z meta> <name from the queue> <stem> 4daae237… 0f56209d…
  de_mde_sign,de_sign_emitter_sha256,rc_callset_block` with 446bj2z's env (rebuild535.sh's block) plus the two knobs.
- **Expected config diff vs 446bj2z: exactly `rc_callset_block`, `de_mde_sign`, `de_sign_emitter_sha256`; vs 446bj2zcp: exactly
  `de_mde_sign`.** pds and mse must equal 446bj2z's (and 446bj2zcp's) to ≤ 1e-6 raw.

## Why (H1, the production emitter, 3 emission seeds; `price537.py`)
| comparison | f 1 | f 0.5 | f 0.35 | cmp |
|---|---|---|---|---|
| hook on the block (446bj2zcph − 446bj2zcp) | **+0.0016** [+0.0002, +0.0033] | +0.0011 [+0.0003, +0.0018] | **+0.0014** [+0.0008, +0.0020] | +0.0006 [−0.0008, +0.0018] |
| members of that increment (fid / jac / reach / nmae raw) | +0.0006 / +0.0002 / +0.0035 / −0.0018 | +0.0010 / +0.0002 / −0.0013 / −0.0026 | +0.0010 / +0.0003 / +0.0000 / −0.0027 | +0.0010 / +0.0003 / −0.0004 / −0.0001 |
| block + hook vs 446bj2z | +0.0067 [+0.0047, +0.0088] | +0.0059 | +0.0059 | +0.0059 |
| block alone (446bj2zcp's construction) | +0.0051 | +0.0048 | +0.0045 | +0.0053 |
| hook alone on 446bj2z | +0.0002 | +0.0004 | +0.0007 | −0.0002 |

- **Mechanism:** under `kept` the ~616 swapped-in genes per target are called but unsteered (design sign, 0.508 accuracy); they
  lie in the hook's territory (not head, not the steered block), so the hook gives the CLEAN4-covered ones CLEAN4's sign with an
  MDE-sized magnitude. fid recovers +0.0006 … +0.0010 raw of the block's −0.0013 … −0.0018 cost; nmae falls −0.0018 … −0.0027.
- **Deliverability:** unsuppressible calls 8.3 per target with the hook vs 20.4 for the block alone (base 7.0).
- **A/B/C reach of the mechanism:** CLEAN4 τ 0.05 covers 29 / 19 / 33 % of A / B / C's swapped-in genes vs 19 % on H1
  (`coverage537.py`); on the blob-free K562 fold CLEAN4 (without its K562 banks) is right on 0.553 of the covered swapped-in
  gate genes (H1 0.582; `foldsign537.py`).
- **At A/B/C-like CLEAN4 quality** (the hook's source attenuated to r 0.17, 2 seeds, RESULTS §3.2 P8): the increment is
  +0.0014 / +0.0007 / +0.0008 / +0.0005 (f 1 / 0.5 / 0.35 / cmp; fid +0.0005, nmae −0.0012 … −0.0014 raw) — ~55–65 % of the
  unattenuated value, CIs above 0 at f ≤ 0.5. This is the central estimate.

- **At zero sign information** (P9: CLEAN4's signs permuted among its covered genes, 3 seeds): the hook on the block is
  −0.0006 / +0.0001 / +0.0002 / +0.0001 (f 1 / 0.5 / 0.35 / cmp; it still stops the block's leaks: jac +0.0003, fid +0.0003
  raw); its break-even sign accuracy is ≤ 0.52, and CLEAN4's accuracy over the hook's territory is 0.59–0.67 on the three
  held-out public lines (H1 0.60). For October (D/E/F) this read is evidence about CLEAN4's sign information on lines of the
  challenge's kind: POSITIVE says it is there on A/B/C too; a NULL says A/B/C sit near break-even — the bins below then drop
  the hook, and P9 prices that choice at ≈ 0 either way at f ≤ 0.5 (the leak-stopping part is sign-free), so the NULL branch
  loses little. Only Δfid < −0.0003 raw (anti-informative signs) says the hook would hurt on lines like A/B/C.

## Noise (RESULTS §2)
Both knobs keep pass 1 and change pass 2 (RESULTS §2.6), so the read carries at most the re-roll sd of today's law (B1:
reach 0.0013 raw, 0.00024 overall); H1's paired sd for exactly this pair (hook on the block, 3 seeds, board-scaled) is 0.00011
overall; fid ≈ 0.00003, jac ≈ 0.00001, nmae ≈ 0.0001 raw. **Primary = Δ(fid + nmae) in overall units vs 446bj2zcp**
(σ ≈ 0.00005 overall), with jac and reach reported.

## Predictions (board vs 446bj2zcp)
- fid: +0.0006 … +0.0010 raw × the sign-prior class's fid transfer (0.07 … 0.67×, docs/H1-BOARD-CALIBRATION.md rule 3) =
  +0.00004 … +0.0007 raw.
- nmae: −0.0018 … −0.0027 raw × T, T = 496a's board Δnmae / −0.00146 (10-05). Before 496a reads: T 0.35 … 1.76 → −0.0006 …
  −0.0048 raw.
- jac +0.0001 … +0.0002 raw; reach unpriced (H1 −0.0013 … +0.0035 by regime).
- **Δ(fid + nmae) overall +0.0002 … +0.0017; overall +0.0003 … +0.0017** (reach at ×0); central (P8, A/B/C-like CLEAN4
  quality, nmae at T ≈ 1): **≈ +0.0006 … +0.0008**. Versus 446bj2z: 446bj2zcp's own read + this.

## Bins (primary Δ(fid + nmae) overall vs 446bj2zcp)
- **≥ +0.0004: POSITIVE.** October runs the block WITH the hook (`DEF_DAY_CALLSET_BLOCK=prop` + `DEF_DAY_MDE_SIGN=clean4:0.05`).
- **+0.0001 … +0.0004: PARTIAL.** The hook rides on the block only if 496a TRANSFERS (the mechanism confirmed twice).
- **< +0.0001: NULL / NEGATIVE.** No hook on the block; October keeps whatever 446bj2zcp and 496a decided separately.
- **Falsifiers (diagnose before reading):** pds or mse ≠ 446bj2zcp's by more than 1e-6 raw (stop); |Δjac| > 0.001 raw vs
  446bj2zcp (a call leak: compare `call_set_unsuppressible` in the metas); Δfid < −0.0003 raw (the hook's signs reached the
  called genes the wrong way round: check that `de_mde_sign`'s audit shows the swapped-in genes covered).

## If 446bj2zcp read NULL/NEG on 10-01 with Δjac raw ≥ +0.002 and Δfid < 0
The composition transferred and the unsteered signs cost fid: this arm is the repair. Same bins; the read vs 446bj2z decides
whether the block ships at all (POSITIVE vs 446bj2z ≥ +0.0012 on Δ(fid + jac), 535's bin, with the hook on).
