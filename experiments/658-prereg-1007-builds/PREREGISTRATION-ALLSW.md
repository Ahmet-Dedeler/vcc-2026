# 658a — ALLSW, the all-switches production arm that doubles as 637's factorial validation (PREREGISTRATION)

Council session [7772dc], assigned by the coordinator 10-04 06:18Z (NEXT-SLOT-DECISION: "658 assigned"). Written 10-04 from
06:40Z and committed before ALLSW is built, before 496a's entry (10-05 s1) exists, and before tercilec is decided (10-06 step 9).
Every number below comes from captured board entries (`data/board-entries/`) through 531's `readdef.delta`.
**Reader:** `read658.py allsw` (tests: `test658.py`, 23/23).

> **AMENDMENT 1 (2026-10-04 06:48Z; ALLSW not yet built or read; 496a, 524dr and tercilec not yet read or decided).**
> - R1 now uses 637's **v2 leak band** (centre −0.0002, band +0.0002 / −0.0002 / −0.0006). [d4cd5f] re-centred it in 0b484be2,
>   on the 9 pairs in §7; the coordinator asked for this change before the read.
> - R1 PASS iff |ΔDE(F) − ΔDE(Y) − (−0.0002)| ≤ 0.0006. With X's measured leak (+0.0000573) that is I ∈ [−0.0008573,
>   +0.0003427]. v1's window was [−0.0009573, +0.0002427].
> - Under v2, R1's PASS window contains R2's whole ADDITIVE bin. §6's "R1 FAIL by the centre" sliver no longer exists.
> - R2, the falsifiers, the prediction and every other rule are unchanged.
> - `read658.py` pins v2 (`ALLSW_LEAK_BAND`), and each verdict records its band. `test658.py` is 25/25.

## 1. The arm
- **ALLSW = 498p2g's production build, unchanged, plus the DE switches decided by build time.** Concretely, it is chain498.sh's
  round-2 `498p2g` build call (`498-combined-bulk/chain498.sh`, `build()`) with:
  - `REF_META` = `submissions/498p2g-cap15-partial50-clean4q80-lodo6-ballast-sharpenp2g-stacklodo05hi6load-ctxlevelr-mde010coarse-zr5.h5ad.meta.json`;
  - Q = `tercile:2.5,2.0,1.75`;
  - `VCC_BULK_LOAD=experiments/495-lodo-cap6/work/dump_h1x2_cd4x15_lodo05hi6.npz:d58750b83915d40d488125926614c83793bd6d7260ceef1192b4b10c19eeb45c`;
  - `VCC_BULK_CTX_SCALE=A=1.0687,B=1.1223,C=1.1259` (498p2g's recorded factors; **no re-level**: DE switches cannot move ρ);
  - `VCC_DE_MDE=0.1,0.1,0.3,coarse`, `VCC_DE_ZERO_REPAIR=-5`;
  - chain498's unsets kept (`VCC_SAMELINE VCC_DE_PRIOR_GATE VCC_UNION_SOURCES VCC_RC_AUDIT_DUMP`);
  - **`VCC_DE_MAG=dn_k250_zs` added** (zs, DECIDED 09-28);
  - **`VCC_DE_MDE_SIGN=clean4:0.05` added only if `mde_sign` is DECIDED** from 496a's 10-05 read by its own bins (TRANSFERS).
    PARTIAL or FAILS leaves it out.
- **EXPECT (run_arm_generic's exact config diff vs 498p2g's meta):** `de_mag`. If the hook is on, also `de_mde_sign`.
  `de_sign_emitter_sha256` joins EXPECT only if the post-10-06 emitter's sha differs from 498p2g's `edb5b8dc…`. That is benign
  without the hook (no output effect, 499 F2; the all-unset KA proves it) and required with it.
- **Builder / emitter:** the post-10-06 tree after APPLY DAY step 8's commit and every step-6 KA PASS. BUILDER_SHA and
  EMITTER_SHA are taken from that commit.
- **Queue:** 630 → 524dr2 → **ALLSW** → DNMAG → idall3 (coordinator, 10-04). Built 10-07 after 630 and 524dr2, serial, under
  the memory gate.
- **Board:** `-m` from `docs/MODEL-NAMES.md`'s queue (`next_model_name.py --claim allsw`), `-d` empty.

### 1.1 tercilec stays OUT of ALLSW, even if step 9 decides it (a deviation from the assignment text, with the reason)
- tercilec is a **bulk** switch, the form of `sharpen_pt`: its cut points come from the carried targets only. 599 measured it
  on H1: pds +0.0029, mse +0.0006 raw, DE ≈ 0.
- **On A/B/C it re-tiers exactly 4 targets** (598's KA moves 493h's 4). So an A/B/C board read of it is ≈ 0 by construction.
  Its value is on thin-coverage panels, and its held-out read is C10b.
- **Inside ALLSW it would cost the validation.** The bulk would no longer equal 498p2g's, so F1 and 637's bulk criterion fail
  by design.
- **It would also add noise.** A bulk change shifts the file's single emission stream (the ballast K; 537 v4 §2.6), so the DE
  re-rolls (sd 0.00025) and R2's exact interaction read is lost.
- **Its October-path rehearsal is 524dr2**, which carries the resolved recipe, tercilec included (658b).
- **If the coordinator puts tercilec in anyway:** F1 and R1 are void. The bulk Δ vs 498p2g reads tercilec on A/B/C (expected
  |Δbulk| ≤ 0.0003 overall, else a defect), and the DE side is read only with a 2-sd re-roll tolerance (0.0005).

## 2. Pre-slot gates (all must pass; the coordinator's chain runs them)
1. run_arm_generic's gates: `[ALLSW] PASS`, CLIFF OK, zero-gene audit 0, `vcc prep` dry run then real, novelty rc 0.
2. **Bulk identity before the slot** (F1 checked locally):
   - Restore 498p2g's content with `scripts/restore_prediction_h5ad.py submissions/498p2g-….vcc <tmp>.h5ad`.
   - Run `uv run python experiments/658-prereg-1007-builds/bulkid658.py <ALLSW>.h5ad <tmp>.h5ad --out <work>/bulkid_allsw.json`.
   - It must print BULK IDENTICAL, i.e. 0 differing (context, target, gene) sums. Anything else means ALLSW's bulk is not
     498p2g's: diagnose and do not submit.
   - Delete the restored copy afterwards (12 GB).
   - This works because the DE switches move counts between cells only (zs, hook), never between genes or blocks.

## 3. Prediction (fixed now)
| read | zs only (hook not decided) | zs + hook |
|---|---|---|
| ALLSW − 498p2g | **+0.0005493** overall: nmae **−0.002179** raw, the other five members **0 exactly** | + (496a − 446bj2) member by member |
| ALLSW − 446bj2 | **+0.0031917** → ALLSW **0.2754683** (rank-of-record irrelevant) | + (496a − 446bj2) |

- X = 498p2g − 446bj2 = +0.0026424. Its bulk part is +0.0025851; its DE part is the measured re-roll leak, **+0.0000573**.
- Y = 446bj2z − 446bj2 = +0.0005493: nmae only, the other five bit-identical. With the hook on, add 496a − 446bj2.
- Expectation: **ADDITIVE.** zs is a post-hoc magnitude pass. H1 (496z) found zs and the hook additive.

## 4. Falsifiers (read first; any one stops the read: diagnose)
- **F1:** |Δpds|, |Δmse| raw (ALLSW − 498p2g) ≤ 1e-6. DE switches never move the bulk: 0 / 2.4e-8 on every board pair.
- **F2:**
  - zs only: |Δfid|, |Δreach|, |Δjac| raw (ALLSW − 498p2g) ≤ 1e-6, since zs is call-neutral (bit-identical on 446bj2z).
  - Hook on: within 496a's own leak rule of 496a's deltas: fid 0.001, jac 0.001, reach 0.01 raw.
- **F3:** the meta's config diff vs 498p2g is exactly EXPECT (run_arm_generic enforces this before prep).

## 5. Reads and bins
**R1: 637's validation, exactly as registered** (`readdef_factorial637.validate`; PROPOSAL §1):
- FACT = ALLSW, X_ALONE = 498p2g, BASE = 446bj2.
- Y_ALONE = 446bj2z (zs only), or the synthetic entry 446bj2 + (446bj2z − 446bj2) + (496a − 446bj2) when the hook is on.
- PASS iff:
  - |Δpds(F) − Δpds(X)| ≤ 1e-6 and |Δmse(F) − Δmse(X)| ≤ 1e-6;
  - |ΔDE(F) − ΔDE(Y) − (−0.0002)| ≤ 0.0006 overall. **[Amendment 1: v2 centre −0.0002; was v1's −0.0003.]**
- PASS → the D/E/F feedback plan may pair one bulk and one DE switch per entry (637 §1). FAIL → one switch per D/E/F entry.

**R2: the exact interaction** (new here; sharper than R1, because ALLSW shares 498p2g's bulk and emission draws):
- I = ΔDE(ALLSW − 498p2g) − Σ ΔDE(single reads vs 446bj2), in overall units. Noise ≈ 1e-5 (prior-only arms on shared draws).
- Bins are set in units of the same-law DE re-roll sd (0.00025 overall), the noise of any D/E/F DE read that changes the bulk.

| bin | overall | nmae raw (zs only: the only moving member) | fid / jac raw (hook) | reach raw (hook) |
|---|---|---|---|---|
| **ADDITIVE** | \|I\| ≤ 0.00025 | \|I_nmae\| ≤ 0.0010 | ≤ 0.0001 | ≤ 0.0013 |
| **WEAK INTERACTION** | 0.00025 < \|I\| ≤ 0.0005 | 0.0010 … 0.0020 | 0.0001 … 0.0002 | 0.0013 … 0.0026 |
| **INTERACTION** | \|I\| > 0.0005 | > 0.0020 | > 0.0002 | > 0.0026 |

The raw columns use A/B/C's scale: 1 raw nmae = 0.2521 overall (446bj2z's −0.002179 raw = +0.0005493), and reach's re-roll
rms is 0.00132 raw (537 §2.5). The overall column decides; the member columns locate the interaction.

## 6. What follows (fixed now)
- **R1** decides 637's D/E/F pairing as 637 registered it. R2 does not override R1.
- **R2 ADDITIVE:** the October Σ-based predictions (524dr2, the D/E/F plan) stay additive.
- **R2 WEAK:** the factorial DE bins on D/E/F widen by |I|. Recorded in FINAL-PHASE-FEEDBACK-PLAN.
- **R2 INTERACTION:** DE switches are read as single-switch D/E/F entries whatever R1 says. Every Σ prediction that stacks a DE
  switch on the cap6+P2g bulk carries I.
- **If ΔDE(ALLSW − 498p2g) ≤ 0:** zs (and the hook) give nothing on the current bulk. Their October values are re-opened by the
  coordinator with a written reason; a decided switch is not reversed on less than this.
- **[Superseded by Amendment 1: under v2 the sliver below is gone.]** **R1 and R2 can disagree only in a sliver.** With X's measured leak (+0.0000573), R1 passes iff I ∈ [−0.0009573, +0.0002427]
  (`read658.py` prints this window). R1 is miscentred because 637's centre (−0.0003) is 0.00036 below this X's leak. So R1
  FAIL with R2 ADDITIVE is possible only for I ∈ (+0.0002427, +0.00025]. That case is read as **R1 FAIL by the centre**: the
  pairing waits for 637's owner to re-centre the leak band (§7) and is not adopted from this read.

## 7. Registered before the read: 637's leak band no longer matches the board
637's premise was "a bulk switch leaks into DE by a consistently negative −0.00020 … −0.00039". With the newer pairs (all
bulk-only, `readdef.delta`, DE part in overall units):

| pair | DE leak |
|---|---:|
| 495 − 446bj2 (cap 6) | +0.0000034 |
| 498p2g − 495 (P2g on cap 6) | +0.0000538 |
| 479p2g − 446bj2 (P2g) | −0.0001986 |
| 479sp2 − 446bj2 (S_P2) | −0.0003655 |
| 446bj2 − 446b (LODO key) | −0.0003882 |
| 471c − 446b (union swap) | −0.0004120 |
| 460dr − 453dr (October-path LODO) | +0.0001928 |
| 433k − 433h, 446b − 433k (637's table) | −0.000213, −0.000241 |

- n 9: mean **−0.000174**, sd 0.000213, rms 0.000266, range −0.000412 … +0.000193.
- That is the same-law re-roll (537 v4 §2.6: a bulk step shifts the single emission stream), not a fixed negative offset.
- 637's band (−0.0001, −0.0003, −0.0005) misses the 3 newest pairs.
- **Recommendation to 637's owner ([d4cd5f]), DONE in 0b484be2 (v2 default; v1 kept for the registered past reads):** before Oct 22, re-centre LEAK_CENTRE to −0.0002 and widen the band to
  (+0.0002, −0.0002, −0.0006). R1 here was to be read as registered (−0.0003); Amendment 1 moved it to v2 before any read.

## 8. Class, tries, noise
- **Class:**
  - R1 / R2: scorer mechanics. They test the D/E/F reading design and whether the October switches stack; no model change.
  - ALLSW − 498p2g is also zs's second board read (and the hook's, if on) on the current bulk. That makes it a transfer
    confirmation for de_mag / mde_sign on a second base.
- **Tries:** 1 arm, no tuning. Nothing in it is fitted to A/B/C beyond the already-decided switches.
- **Noise:**
  - Bulk members reproduce to 16 significant figures.
  - ALLSW − 498p2g shares draws (≈ 1e-5).
  - ALLSW − 446bj2 carries 498p2g's measured re-roll, which is known exactly and so is not noise here.

## 9. Read command
```bash
uv run python experiments/658-prereg-1007-builds/read658.py allsw <ALLSW entry id> [--hook <496a entry id>] \
  --out experiments/658-prereg-1007-builds/work/read_allsw.json
```
Exit codes: 0 = read, 2 = members hidden, 3 = falsifier. Board record: `experiments/658-prereg-1007-builds/BOARD-RESULT-ALLSW.md`
and a BOARD-LEDGER row.
