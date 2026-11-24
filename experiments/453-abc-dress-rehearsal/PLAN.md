# 453 — the October day, run blind on A/B/C and scored on the board (plan, 2026-09-26 00:12Z)

**Question.** What does the exact October path (`scripts/def_day.sh`) cost when it runs blind against the
tuned A/B/C incumbent? 452/452b rehearsed its mechanics on synthetic bundles, which cannot be scored. On the
real A/B/C bundle it produces a board-scorable entry. That entry is the only end-to-end known answer we can get
for the final entry before the board goes dark.

**What the blind path changes against the tuned incumbent (446b and its successors).** Each item is a day-of
component that never had its own board point on A/B/C:
1. The day's whitener refit (`union_pruned_def`), the day's 410 key, the CD4 extraction and census, all run
   inside `def_day.sh`. 452 showed a refit is float noise (ρ agrees to 4 dp).
2. DEF_ENV's CD4 sizing: carrying rows at g 0.344 × m (`VCC_CD4_TARGET_G`), not ABC_ENV's all rows at f 0.34.
3. The October stack, H1 ×2 / CD4 ×1 (438's lowest-regret dose), not A/B/C's CD4 ×1.5. 438 prices the CD4 ×1.5
   half at +0.0018 on abc (k = 1 levelled).
4. The NTC-proxy re-level to ρ* 0.9793 per context (448's secant loop), not the matched ρ of the tuned chain
   (446b: 0.9908 / 0.9827 / 0.9873).
5. The same-line stage. On A/B/C it matches A = Jurkat (Nadig): 0 panel rows covered, so under `noprog` it is a
   no-op. If 446bja reads positive, `noprog` is dropped, and A gets 446bja's programme.

**How (revised 00:20Z).** The driver branches on the labels: an A/B/C bundle takes the legacy A/B/C path (ABC_ENV
CD4 sizing, pinned key, `--verify` stages). So the rehearsal runs `scripts/def_day.sh` unmodified, on its default D/E/F
path, over 180's `def-bundle`. That bundle is the A/B/C controls and panel relabelled A->D, B->E, C->F: all three
control files and both CSVs were verified identical to `data/challenge` at 00:19Z. Then:
- `relabel453.py` APFS-clones the final artifact and renames only obs['context']'s categories back (D->A, E->B,
  F->C), checking X, indptr and codes against the source. Cells are named `cell_{i}`, so nothing else carries a
  label.
- `vcc prep` (dry run, then real) against `data/challenge`'s CSVs, then the novelty check.
- `chain453.sh` (pid 47467, armed 00:21Z) does all of it. It waits for chain0926's END and for `work/args.env` (Q,
  FLAGS), which the coordinator writes after the 09-26 reads. It never submits.
- A short-lived `DEF_DAY_CONTEXTS` edit to `def_day.sh` was reverted (00:20Z). The script is byte-identical to HEAD
  `b517b6d2`, and the chain refuses to run otherwise.

**Reading (the preregistration is written before the entry exists).** Δ = 453 − the tuned counterpart carrying the
same riders. The prediction is the sum of the 438 CD4 term, the cliff-model level term and ~0 refit noise.
- Inside the predicted band: the October pipeline is validated end to end.
- Outside it: a day-of component is mis-priced. Attribute it by rebuilding with one component swapped
  (a local build, no slot).

Slot: 09-28 (after the 09-27 riders read). Name: a fresh joke. Description: empty.
