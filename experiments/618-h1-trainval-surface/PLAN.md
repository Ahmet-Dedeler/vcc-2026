# 618 — a second, unworn H1 real-emitter surface: h1-trainval with the FAITHFUL production key

[21dd4c], approved by the coordinator 2026-10-02 ~19:25Z ("Go 618 (faithful)").

## Why
- 614 (K562) and 615 (Orion HEK raw) fail on truth power, and the thinned Orion folds are void (616).
- That leaves h1-test (85 targets) as the only real-emitter DE surface, and it is the most-tuned surface in the repo
  (AGENTS 11.1c rule 5).
- h1-trainval has 126 OTHER targets × 400 cells: the same 2025 Arc H1 chemistry, ~20k UMI, disjoint from h1-test's
  panel (0 shared targets).

## Conditions (coordinator, binding)
1. **KA first.** The faithful-key pipeline (run323's stage / operator / compose, generalised as `run618.py`) must
   reproduce 323's `cap15` on h1-test's targets **bit for bit** before it builds trainval's key.
2. **Controls.** Check whether h1-full's 18,399 controls are the same pool as h1-test's.
   - **Measured (19:2xZ):** h1-test 18,400 controls (unique names); h1-full 18,399 rows but only 15,571 unique names
     (barcodes collide across batches); trainval 6,000 controls, all inside h1-full's names, 2,840 inside h1-test's.
   - **Construction:** 618's truth uses h1-test's 18,400 controls (528's `truth528` precedent: trainval's cells against
     h1-test's controls, KA'd exact).
   - So **the two surfaces share their control pool, and any panel-draw noise estimate from h1-test vs 618 is a LOWER
     bound** (shared control noise cancels in the comparison).
3. **Trainval's prior use, exactly** (grep of every experiment that loads `h1-trainval`):
   - 004 (gene-weight fits, an August closed family);
   - 006 (cosine variants scored on it);
   - **528: lcap β 22 selected by trainval capture; the `ltv` / `ptv` call-propensity keys fitted on trainval's 126
     perturbations**;
   - 537's cp block inherits β 22.
   - **So call-set / lcap rows (the gate's G3 and any cp lever) are WORN on 618 and are marked so. Delivery and hook
     rows are unworn.**
4. **The trust gate** is 614's (G1 zs / G2 coarse MDE / G3 cp vs the board's pairs), with bands committed before any
   run, plus one row re-checking 613d's registered H1 passes (idkg3, idall3) on the fresh draw. That is the first real
   out-of-panel test of the hook instrument.
5. **Heavy stages** run before 10-04 05:00Z or after 10-05 18:00Z, through gate602.
   - **Leakage guard:** `data/donors/h1_2025` never enters a donor, prior or bank. 517's Arc-hook / H1-programme forms
     are not priced here (own line).

## Stages
- **A, the faithful key** (`run618.py <ka|tv> <freeze|cd4|stage|operator|compose>`, a parameterised copy of run323):
  - `ka` = h1-test's 85 targets (323's own axes). `tv` = trainval's 126.
  - CD4: `ka` uses 273's `cd4_compact_h1` (323's input). `tv` gets a fresh extraction with 047's `extract_cd4.py`
    (the D/E/F pipeline's cd4-extract stage; 0.34 GiB) on the trainval panel.
    - KA for the extraction path: the same extractor on h1-test's panel must equal 273's `cd4_compact_h1` array for
      array.
  - The operator is `build_union_whitener.py --basis-prune essential` on the mode's bundle (pruned of its own panel;
    ~14 GB, gate602).
  - **The KA is passed iff `ka` reproduces `323/work/effects.npz`'s cap15 bit for bit** (plus 323's own stage KAs:
    historical H1 pre-whitening and the 319 A/B/C hashes).
- **B, truth and reference:** `truth618.py` (trainval's 126 × 400 cells + h1-test's 18,400 controls; 614 format), then
  `614/ref614.py` → `work/ref_tv`, then the power diagnostic (`615/diag615.py` form) before any gate run.
- **C, wrapper and gate:**
  - `run618e.py`, a copy of run614 with REFDIR = `work/ref_tv`, `native_cap15` = the tv cap15, the four CLEAN4 banks,
    and CLASSES3 from the reference;
  - its own h1ka KA (== 435);
  - then TRUSTGATE.md (614's bands, committed first), the G1 / G2 / G3 runs, and the 613d re-check.
