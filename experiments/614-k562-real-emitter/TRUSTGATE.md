# 614 — Trust gate for the K562 real-emitter DE harness (fixed BEFORE any harness run or lever price)

Written 2026-10-02 ~16:45Z by [21dd4c]. The K562 truth was still streaming, and no harness run, reference or score
existed. Coordinator's assignment; spec: `experiments/631-nmae-headroom/SPEC-614-K562-harness.md` ([d4cd5f]).
AGENTS 11.1 rule 3: the harness prices a lever class only after it retrodicts that class's board pairs. If this gate
FAILS, every 614 number is descriptive only.

## The harness, as it will be built (stages in PLAN.md)
- **Truth:** K562 GWPS (scPerturb ReplogleWeissman2022) through 480's `k562harness480.py`.
  - 18,400 random NTC cells, plus up to 250 cells for each A/B/C panel target with ≥ 100 cells.
  - Then a seeded random 100 of those targets (seed 614), fixed before any run.
  - Native depth (~11k median UMI; it cannot be thinned up to A/B/C's ~20k). 8,248 genes.
- **Reference:** evaluate323's `reference()` logic on that truth, without H1's frozen-axis and 400-cell asserts (the
  same control blocks, Wilcoxon DE and stored statistics).
- **Scorer:** evaluate323's `score()`, i.e. cell-eval2's own DE member functions (vcc2026 preset) on prepared frames.
  This is the exact evaluator's member code, not a re-implementation. The regimes f 0.5 / 0.35 come from 521's
  pools on the K562 reference's own p.
- **Emitter:** run496's driver with the production emitter (RUN426_EMITTER=prod), through a wrapper (`run614.py`) that
  swaps only:
  - the reference directory;
  - the design key (below);
  - CLEAN4's banks: Nadig Jurkat + HepG2, dropping both K562 banks.
- **KA for the wrapper:** pointed at H1's own reference and key, run614 must reproduce run496's H1 stats byte for byte.
- **Design key, first version = PROXY:** fold480's donor consensus from the two Orion lines only (no K562 donor; no
  PerturbMatch K562 replacement).
  - Its scale is matched to the production cap15 key ON H1 ONLY: one scalar s = median |cap15| / median |3-donor
    consensus| over H1's nonzero key cells. No K562 truth is used.
  - If the gate FAILS on the proxy, the faithful key (run323.stage() with Orion-only donors, CD4 kept, and the whitener
    pruned of K562's targets) is built and the SAME gate is re-run once, with this file unchanged.

## Amendment A1 (2026-10-02 16:45Z, before any VALID harness run; bands unchanged)
- **The problem:** every K562 run crashed in the production emitter ("rc_hybrid_k did not meet the per-cell totals").
  - The emitter balances each cell's exact library with sub-gate genes (< 5 CPM, outside the DE universe).
  - The GWPS release keeps only 8,248 expressed genes, so the K562 truth had 0 sub-gate genes; H1 has 7,301.
  - A gate614 "FAIL" written from those crashed runs is VOID (`work/gate614.VOID-crashed-runs.json`): no G1/G2/G3 data
    existed.
- **The fix (`sub614.py`):** the truth gets an H1-shaped sub-gate reservoir.
  - 7,301 synthetic sub-gate genes carrying H1's control sub-gate share (0.58 % of each cell's library), each gene at
    its H1 counterpart's share.
  - Drawn identically in control and perturbed cells, seed 614.
  - They shift every cell's CPM normaliser by the same ~0.58 %, as H1's real sub-gate genes do.
  - **Measured (sub614.json, 16:42Z, before any valid run finished):** 20 of the 7,301 sit just above the 5-CPM gate
    (max 5.68). K562's lower depth raises a low-count gene's per-cell CPM mean, so those 20 enter the DE universe.
    They are null by construction (identical in control and perturbed cells, key 0), so their effect on the gate members
    is negligible. Recorded, not refit.
- The reference is rebuilt on it (`work/ref_k562s`), and so is the G3 lcap key (`key614cp.py`, re-pointed). The design
  key maps by gene name, so the synthetic genes get 0 effect.
- **Reading rule A1-R** (the coordinator, 16:45Z; written before the first valid run): capacity-bound levers (630's UP
  pass, 631 / v7, anything limited by sub-gate donor room) get two things on 614.
  - (a) Every such row is labelled **"H1-reservoir delivery"**: its delivery is H1-shaped by construction.
  - (b) Before any such row is read, the same lever is re-run with the reservoir at **×0.5 and ×2** (library share
    0.29 % / 1.16 %; `sub614.py` with s0 scaled, same seed). If the lever's conclusion (bin or sign) flips across that
    range, 614 decides only its SIGN, not its size or capacity.
  - The trust gate's pairs (zs, coarse MDE, call set) and the sign / design levers (hook forms) are read on the ×1
    reservoir as specified.
- **Everything else in this file is unchanged:** the gate pairs, bands, seeds and f 1.

## The gate: the board's known DE pairs, 3 seeds each, at f 1 (the K562 truth's own BH set)
| pair | harness construction | board raw Δ (611) | required on K562 |
|---|---|---|---|
| **G1 zs** (446bj2z − 446bj2) | the same run, pre-zs vs `--zs dn_k250_zs` (prodem, coarse MDE, zero repair) | nmae −0.00218 | Δnmae < 0, and board / K562 ratio in [1/3, 3], i.e. K562 Δnmae in **[−0.0065, −0.00073]** |
| **G2 coarse MDE** (433a − 417e) | prodem without MDE vs prodem + MDE 0.1 / 0.1 / 0.3 coarse (no zs, no zero repair) | nmae −0.01679 | Δnmae < 0, ratio in [1/3, 3], i.e. **[−0.0504, −0.0056]** |
| **G3 call set** (446bj2zcp − 446bj2z) | 528's call-propensity block ("kept") vs the CPM block, on the zs base | fid −0.00358, jac +0.00399 | **Δfid < 0 AND Δjac > 0** (signs) |

- **PASS iff G1, G2 and G3 all hold** (mean of seeds 0–2).
- Reported beside it, with no rule attached:
  - each pair at f 0.5 / 0.35;
  - G2's reach (board +0.0024; H1 got this sign right only at f 0.35);
  - G3's reach and nmae (611: H1 never prices reach for call-set changes);
  - per-seed values;
  - the same pairs on H1@20k (611's existing numbers) and, if built, on **H1@11k** (h1-test binomially thinned to
    K562's median depth; the coordinator's control). H1@11k separates "cancer line vs H1" from "11k vs 20k depth" in
    every reading.
- Depth note: the MDE add depends on counts per cell, so K562's 11k depth gives larger adds than A/B/C's 20k. The
  ×3 ratio band allows for that. Lever prices after a PASS are also reported d20-rescaled (633's convention).

## After a PASS (not before): register, then price
Each lever's prediction is committed before its K562 run:
- 630 v3 `b_k100_q20r3p10` (K 1.0, vs 2–3 if bracketed);
- 631 / v7 when [d4cd5f] delivers;
- the hook forms on the same base: the plain CLEAN4 hook (c4h) vs idall3.
For the hook, CLEAN4-minus-K562 makes K562 the weakened-CLEAN4 stress case (INSTRUMENTS.md rule 1). Read the hook's
K562 row as the stress tail, not a typical line.

## Gate result and the fallback rule (written 2026-10-02 ~17:55Z, BEFORE the H1@11k number exists)
- **The gate FAILED** (`work/gate614.json`, 4db74d02): G1 PASS; G2 FAIL (nmae −0.00075); G3 FAIL (fid +0.0001).
- **Fallback (coordinator, conditional waiver).** H1@11k's coarse-MDE Δnmae (`h1t614.py`, f 1, mean of seeds 0–2) decides
  the faithful-key re-run.
  - **If it is ≤ −0.0056** (≥ 1/3 of the board's −0.0168), depth alone does not shrink the lever, so the failure points
    at the K562 line or the proxy key. **The faithful key is then built** (run323.stage(), Orion-only donors, CD4 kept,
    the whitener pruned of K562's targets), and this gate is re-run once, unchanged.
  - **Otherwise the re-run is WAIVED**: depth (and the power it costs) explains the failure, and a better key cannot
    restore it.
  - The waiver or the re-run is recorded below with the H1@11k number, before anything else is read off 614.

## The H1@11k number and the fallback rule's outcome (recorded 2026-10-02 20:20Z, before anything else is read off 614)
`h1t614.py` (`work/h1t614.json`): h1-test thinned to K562's median library, 3 seeds, f 1, raw Δ.

| pair | H1@11k | H1@20k (611) | board |
|---|---|---|---|
| G1 zs Δnmae | **−0.00337** | −0.00183 | −0.00218 |
| G2 coarse MDE Δnmae | **−0.00648** (reach +0.0009, fid +0.0010, jac −0.0003) | −0.00966 | −0.01679 |

- **Rule outcome:** H1@11k's coarse-MDE Δnmae is −0.00648 ≤ −0.0056, so depth alone does not shrink the lever below 1/3
  of the board's. Per the rule written above (~17:55Z, before this number existed), **the faithful-key re-run is
  TRIGGERED**: run323's stage / operator / compose with Orion-only donors (CD4 kept, the whitener pruned of K562's
  targets), then this gate once more, unchanged.
- **Known before this number, and relevant to whether the re-run can pass:**
  - 614's own diagnosis: the K562 truth has a median of 5 DE genes per target, and nmae scores 34 / 100 targets;
  - 615 (Orion HEK raw) stopped on the same power diagnostic.
  - A better key cannot add DE genes to the truth, so a faithful-key re-run is likely to fail G2 again on power.
- **Execution:** run618.py (run323 parameterised) makes the faithful key cheap to build. The coordinator decides whether
  to run it as the rule says, or waive it with a written reason given the power evidence. Until they rule, it is queued
  behind 618's gate.

## The faithful-key re-run: scope and construction, fixed BEFORE it runs (coordinator, 20:25Z: run it, no waiver)
**Scope (coordinator note 1).** H1@11k changed depth only: it kept H1's 400 cells per target and H1's response sizes. So it
rules out depth and nothing else.
- Still open as causes of 614's FAIL: cells per target (≤ 250 vs 400); K562's line biology (its knockdown and response
  sizes, i.e. the truth's power, median 5 DE genes per target); and the key.
- **The re-run tests ONLY the key.**
  - A PASS says the proxy key caused the FAIL, and gives a validated K562 harness.
  - A FAIL closes 614. It does not separate cells per target from line biology.
- **The faithful key** (`fkey614.py`, run323's stage / operator / compose through run618's validated, KA'd code, with
  exactly these intended differences):
  - (1) targets = the K562 reference's 100;
  - (2) donors = orion_hct116 + orion_hek293t ONLY (K562 is the truth);
  - (3) no PerturbMatch, which is K562 data; the builder refuses its K562-slot replacement without a K562 donor anyway;
  - (4) CD4 extracted for these 100 targets with 047's extractor (kept, as the rule says);
  - (5) the operator refit on the 2-donor union, pruned of the 100 targets;
  - (6) arms without PM blending only (cap15 is the one used);
  - (7) 323's H1-only KAs (historical pre-whitening, the 319 A/B/C hashes, 77 complete CD4 targets) do not apply.
    618's ka panel already proved this code path.
- **The gate:** the same 12 runs as the first gate (G1 / G2 / G3 × seeds 0–2, `--bulk base`, the same flags). Only the
  key changes: `run614.py k562f / k562fcp`. The same bands, read by `gate614.py` on the `*f` dirs (`gate614f.json`).
- **Priority (coordinator note 2):** behind 618's KA / ref_tv and 655 / 656 in gate602. Inside the 10-04 05:00Z – 10-05
  18:00Z window it waits; there is no deadline.
- **Amendment (2026-10-03 ~15:00Z, before any re-run output existed):** the first `fkey614.py stage` died on the builder's
  `donor_weights()` guard (3-donor production order only). Production's weights are uniform when VCC_DONOR_WEIGHTS is unset
  (asserted), so **difference (8): donor weights = 1.0 for each of the 2 Orion donors.** The inputs were re-frozen after
  the fix; the old manifest is kept as `inputs.json.failed-stage-1452Z`. No output had been written.

## Faithful-key re-run verdict (2026-10-03 17:39Z): **FAIL**, so 614 is CLOSED (`gate614f.py`, `work/gate614f.json`)
- G1 zs: −0.00299, **PASS**.
- G2 coarse MDE: −0.00186, **FAIL** (band ≤ −0.0056). It was −0.00075 with the proxy key: 2.5× better, still ~1/9 of the
  board.
- G3 cp: fid +0.00147, **FAIL** (needs < 0); jac +0.00034.
- **Reading, per the scope note written before the run:** the re-run tested only the key, and it fails, so 614 closes.
  The key is not what fails the K562 harness.
- A FAIL cannot separate cells per target from line biology. The truth-power diagnosis (median 5 DE genes per target) is
  the standing explanation.
