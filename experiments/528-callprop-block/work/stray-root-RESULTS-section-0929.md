# Stray section found appended to the ROOT RESULTS.md (uncommitted, written ~2026-09-29 by worker 528's run;
# superseded by experiments/528-callprop-block/RESULTS.md). Moved here by the coordinator on 2026-10-01T13:04Z when the root file was restored.

## 4. Real production-emitter H1 runs (IN PROGRESS; seed 0 so far)
Harness: `run528.py` = 435's (run496, RUN426_EMITTER=prod: CallSet(400, 2200, 1.9), CLEAN4 q80 prior, coarse MDE, + the zs
pass = 446bj2z's DE construction on H1) with 528's patched emitter; KA cpm s0 == 435 bit for bit (§1). Each arm vs 435's
own run of the same seed. Units: overall (192's INV / 6); board columns = global521.board_cols [reach ×1 | ×0 | ×−0.76].

**Seed-0 reads (pre-zs; the zs pass changes nmae by ≤ 0.00006 and nothing else):**
- `lcap` (coupled: the swapped-in genes are called AND take the block's CLEAN4 steering): jac **+0.0090** raw (96 % of the
  emulated +0.0094), fid −0.0016, nmae **+0.0040** (worse), reach +0.0033 → overall +0.0027 (f 1) / +0.0020 (0.5) / +0.0020
  (0.35) / +0.0035 (cmp). Deliverability is not the problem: swapped-in genes are called at 0.971 (0.966 at < 1.5
  counts/cell, 0.973 at 1.5–3) against 0.988 for the kept block and 0.984 for today's block; unsuppressible 4.8 / target
  (base 7.0). The nmae cost is the block's treatment on these genes: their median delivered |lfc| falls 0.121 → 0.085
  (the kf shrink) against a truth median of 0.314 (they are the responsive genes), +0.0035 of the +0.0040.
  fid: the called swapped-in genes' signs are right 0.509 (CLEAN4) where the swapped-out were 0.516.
- `shuf` (the falsifier): jac **−0.0052** (emulated −0.0052), overall −0.0059 (f 1): the gain is the dropout information.
- `lcap-dec` (block_steer "cpm": call the propensity block, keep steering the CPM block): nmae pinned (±0.00000), jac
  +0.0092, but ~30 % of the steered-not-called CPM genes cannot be held below BH (181 leaked calls / target, uniform across
  target strength): fid −0.0027 and jac +0.0069…0.0077 at f ≤ 0.5. Overall +0.0064 (f 1) / +0.0034 / +0.0027 / +0.0035.
- A third form, `lcap-kept` (block_steer "kept": steer only the called block's CPM part; swapped-in genes are called
  unsteered, swapped-out genes are neither called nor steered: no leak, no shrink on the responsive genes), is queued.
