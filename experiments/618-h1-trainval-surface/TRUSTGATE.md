# 618 — Trust gate for the h1-trainval real-emitter surface (fixed BEFORE any 618 emitter run)

Written 2026-10-02 ~19:35Z by [21dd4c].
- **State at writing:** stage A (the faithful key) was at its KA panel's stage step, and stage B (truth / reference) was
  queued. No 618 emitter run, reference or key existed.
- **Coordinator's condition 4:** 614's trust gate, bands committed before any run, plus one row re-checking 613d's
  registered H1 passes on the fresh draw.
- AGENTS 11.1 rule 3: if this gate FAILS, every 618 number is descriptive only.

## The surface
- **Truth:** h1-trainval's 126 targets × 400 cells (disjoint from h1-test's 85) against h1-test's own 18,400 controls,
  on the 18,080-gene H1 axis (`truth618.py`).
- **Reference:** 614's `ref614.py` → `work/ref_tv`.
- **Key = the FAITHFUL production key:** run323's stage / operator / compose on this panel (`run618.py tv`):
  - donors k562_gwps + orion_hct116 + orion_hek293t (the donor tables, clean per 616);
  - CD4 extracted for this panel;
  - the operator refit and pruned of these 126 targets;
  - h1_2025 never a donor.
  - **It is used only if `ka618.py` passed:** the same pipeline on h1-test's panel reproduces 323's cap15, operator and
    CD4 bit for bit.
- **Emitter:** run496's driver, production emitter, through `run618e.py` (a run614-style wrapper):
  - REFDIR = `work/ref_tv`, `native_cap15` = the tv cap15, the four CLEAN4 banks.
  - CLASSES3 = the reference's control-library quartiles. These are h1-test's controls, so they must equal H1's
    constants (asserted).
  - Its own h1ka KA (== 435) runs before any tv row is read.
- **Shared controls:** the controls equal h1-test's, so every 618 − h1-test comparison cancels control noise. Any
  panel-draw spread read from the pair is a LOWER bound.

## The gate: 614's pairs and bands, UNCHANGED (3 seeds, f 1 = the trainval truth's own BH set)
| pair | construction | board raw Δ (611) | required |
|---|---|---|---|
| **G1 zs** (446bj2z − 446bj2) | the same run, pre-zs vs `--zs dn_k250_zs` (prodem, coarse MDE, zero repair −5) | nmae −0.00218 | Δnmae in **[−0.0065, −0.00073]** |
| **G2 coarse MDE** (433a − 417e) | prodem without MDE vs prodem + MDE 0.1 / 0.1 / 0.3 coarse (no zs, no zero repair) | nmae −0.01679 | Δnmae in **[−0.0504, −0.0056]** |
| **G3 call set** (446bj2zcp − 446bj2z) | 528's call-propensity block ("kept") vs the CPM block, on the zs base (zero repair −5) | fid −0.00358, jac +0.00399 | **Δfid < 0 AND Δjac > 0** |

- **PASS iff G1, G2 and G3 all hold.**
- **G3 is WORN on this surface:** 528's lcap β 22 was selected by h1-trainval's capture. G3's result is reported, but a
  G3 PASS here is weaker evidence than on an unworn surface. The gate's verdict is also reported with G3 excluded
  ("delivery-only PASS"), and that is the verdict for delivery and hook levers.
- Reported with no rule attached: f 0.5 / 0.35, reach, per-seed values, and the zero-total counts (`zerocount654.py`).
- **Power diagnostic before the gate:** median n_real, targets with n_real ≥ 10, against h1-test (89; 58 / 85). If it is
  K562-like (median ≤ 10 or < 45 % of targets ≥ 10), that is said before the gate runs.

## The 613d re-check on the fresh draw (registered now)
- **Instrument:** 613d's add-rule instrument (`613/block613d.py` + `block613d5.py`), run on the trainval surface.
  - Inputs: 480's gate rules on `tv-truth.h5`; the tv cap15 key; 521's pools from `ref_tv`.
  - It runs through 616's sandbox, so the code is unchanged and only its H1 inputs are redirected.
- It predicts, at f 1 / 0.5 / 0.35:
  - **c4h − no hook**, **idall3 − c4h** (and idkg3 − c4h);
  - for the real emitter's `run618e tvc4h / tvidall3` (3 seeds) against the zs base.
- **Rule:** PASS iff |613d − real| ≤ 0.0006 raw nmae (the tolerance of A3's registered test) for each contrast at every
  regime with ≥ 20 scored targets.
- **Reading:**
  - PASS → 613d's hook validity is not specific to h1-test's panel.
  - FAIL → its two h1-test passes do not generalise across panels, and hook prices must come from the real emitter on
    each surface.
- This row is run only after the gate's delivery-only verdict is PASS. If the gate fails, it is descriptive.

## Verdict (23:32Z): **PASS**, all three pairs; the delivery-only verdict PASSES
- G1 −0.00163 and G2 −0.01473 (both inside their bands). G3: fid −0.00408 < 0 and jac +0.00872 > 0 (worn).
- `work/gate618.json`; RESULTS.md has the table.

## 613d re-check verdict (12:03Z): **FAIL** by the registered rule
- c4h − base and idkg3 − c4h pass at every regime.
- idall3 − c4h misses at f 1 (0.00082) and f 0.5 (0.00067). 613d under-states idall3's gain on this panel.
- Registered reading: 613d's h1-test passes do not generalise to every hook form across panels, so hook prices come
  from the real emitter on each surface. RESULTS.md has both tables.
