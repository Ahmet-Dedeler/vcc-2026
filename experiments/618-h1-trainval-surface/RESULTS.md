# 618 — a second, unworn H1 real-emitter surface (h1-trainval, faithful production key)

[21dd4c], from 2026-10-02 ~19:20Z. The design and the coordinator's conditions are in PLAN.md; the gate (614's bands)
is in TRUSTGATE.md; the stack predictions are in PRED-618-stacks.md.

## Stage A — the faithful key: KA PASS (condition 1), 22:00Z (`ka618.py`, `work/ka618.json`)
`run618.py ka` (run323's stage / operator / compose, parameterised by panel) on h1-test's 85 targets reproduces 323:
- `effects.npz`: **bit-equal on every array** (cap15, plain15, plain30, uniform, parity15, shuffled_cap15, permutation,
  genes, targets);
- the operator `whitener_h1_pruned.npz`: **bit-equal** (evals, evecs);
- the CD4 extraction (047's `extract_cd4.py` on h1-test's panel): **equal to 273's `cd4_compact_h1`, every .npy / .txt**;
- 323's own stage KAs pass under the current builder: historical H1 pre-whitening ≤ 1e-10, and 319's A/B/C pre / cand
  hashes.

The same code builds trainval's key (`run618.py tv`: CD4 extracted for its panel, operator refit and pruned of its 126
targets).

## Stage B — truth and power (before any gate run)
- **Truth** (`truth618.py`): 126 targets × 400 cells (disjoint from h1-test's 85) plus h1-test's own 18,400 controls,
  on the 18,080-gene H1 axis.
  - 440M nonzeros; library median 19,964 perturbed / 20,190 control.
- **480's gate rules on it** (`gate-h1_tv.npz`): median **149** DE calls per target, **94 / 126** targets with ≥ 10,
  against h1-test's 90 and 59 / 85. H1-class power, unlike 614 / 615.

## Ops notes
- The big gate tickets first headed gate602's FIFO queue (19:30–19:51Z) and blocked smaller jobs. They were withdrawn
  and moved to a polite wait outside the queue, but that wait could not be met: the compressor holds ~16 GiB, so
  free+inactive stays ~21 GiB even when quiet.
- 4d0207 shipped gate602 v5 (bounded backfill) and v6 (opt-in quiet admission).
- Stages A + B now run as one serial chain (`chain618ab4.sh`), so quiet jobs never overlap. The ka operator
  quiet-started at 21:54:48Z.

## Stage A (tv) and Stage B — built 22:00–22:10Z
- **Trainval faithful key** (`run618.py tv`, inputs committed c9bf7af8): stage 7.5 GB, operator 13.63 GB, compose 4.5 GB
  (peaks from the gate log).
  - The same validated code as the KA panel. CD4 was extracted for the 126 targets, and the operator was refit and
    pruned of them.
- **Reference** `work/ref_tv` (614's `ref614.py`, peak 12.73 GB; declared 17 from the extrapolation, so the estimate ran
  cold).
- **Power diagnostic, read before any gate run** (`diag618.json`): **NOT K562-like.**

| truth | median n_real | q90 | targets with n_real ≥ 10 | DE calls | DOWN share |
|---|---|---|---|---|---|
| **h1-trainval (618)** | **149** | 2,614 | **94 / 126 (75 %)** | 110,502 | 0.522 |
| h1-test (614's ref_h1ka) | 89 | 2,691 | 58 / 85 (68 %) | 65,130 | 0.526 |
| K562 (614) | 5 | 64 | 34 / 100 | 2,213 | 0.435 |

So the gate runs (TRUSTGATE: 614's bands, plus the delivery-only verdict because G3 is worn here).

## The trust gate: **PASS** on all three pairs; the delivery-only verdict PASSES (`gate618.py`, `work/gate618.json`, 23:32Z)
Wrapper KA first (`ka618e.py`): run618e h1ka == 435's H1 stats, bit for bit. Then 614's pairs and bands, unchanged:
3 seeds, f 1, raw Δ, `--bulk ps`.

| pair | 618 h1-trainval | band / rule | board (611) | h1-test (611) | board / 618 |
|---|---|---|---|---|---|
| **G1 zs** Δnmae | **−0.00163** | [−0.0065, −0.00073] | −0.00218 | −0.00183 | 1.34 |
| **G2 coarse MDE** Δnmae | **−0.01473** (reach +0.0005, fid −0.0001, jac −0.0001) | [−0.0504, −0.0056] | −0.01679 (reach +0.0024) | −0.00966 | **1.14** |
| **G3 call set** (WORN here) | **Δfid −0.00408, Δjac +0.00872** (nmae +0.0011, reach −0.0017) | fid < 0 AND jac > 0 | fid −0.00358, jac +0.00399 | — | fid 0.88, jac 0.46 |

- **618 is the council's second validated real-emitter DE surface.** It is the first one on an unworn panel: 126 targets,
  disjoint from h1-test's.
- **Better calibrated than h1-test for delivery.**
  - On the coarse MDE, 618 reads 88 % of the board (board / 618 = 1.14); h1-test reads 58 % (1.74).
  - On zs, 618 is 75 % of the board (1.34); h1-test 84 % (1.19).
  - The trainval panel's response distribution is closer to the board's than h1-test's for add-based delivery.
- **G3 (worn, because 528's lcap β was selected on trainval):** fid matches the board (−0.0041 vs −0.0036), and jac is
  2.2×.
- **Zero totals** (`zerocount654.py`): 0 in every zero-repair run. The no-repair G2 arms (the board pair's construction)
  have exactly 1 gated zero-total cell per run, in 1 target.
- **Next (chain618c → chain618d):** the hook forms (c4h / idall3 / idkg3 × 3 seeds); the registered 613d re-check;
  the two magnitude stacks; then 614's faithful-key re-run.

## The registered 613d re-check on the fresh panel: **FAIL** by its rule (`recheck618.py`, `work/recheck618.json`, 12:03Z)
- **Rule:** |613d − real| ≤ 0.0006 raw nmae for every contrast at every regime with ≥ 20 scored targets (94 / 85 / 82).
- **Comparison:** 613d's add-rule instrument on the trainval surface (block613d / block613d5 unchanged, through 616's
  sandbox map) against the real emitter (`run618e`, 3 seeds, scored by `score618.py`).
- **First attempt:** it failed at 01:53Z on 612's exact bulk-sums guard. Zero repair's +1 floor reacts to the hook on one
  target (0–4 cells, |Δ| 1); the guard is now bounded and recorded (d15286b6).

| contrast | 613d f 1 / 0.5 / 0.35 | real emitter f 1 / 0.5 / 0.35 | verdict |
|---|---|---|---|
| c4h − no hook | −0.00218 / −0.00079 / −0.00149 | −0.00227 / −0.00074 / −0.00139 | **PASS** (errors ≤ 0.0001) |
| idkg3 − c4h | −0.00082 / −0.00030 / −0.00053 | −0.00116 / −0.00076 / −0.00073 | **PASS** (≤ 0.00046) |
| idall3 − c4h | −0.00162 / +0.00011 / −0.00029 | −0.00244 / −0.00056 / −0.00074 | **FAIL** at f 1 (0.00082) and f 0.5 (0.00067); f 0.35 within (0.00045) |

- **Reading (the registered FAIL branch):** 613d's h1-test passes do not fully generalise across panels. It reproduces
  the plain hook almost exactly and idkg3 within tolerance. It UNDER-states idall3's advantage on this panel, i.e. it is
  conservative for idall3.
- **Hook prices therefore come from the real emitter on each surface,** not from 613d.

## Hook forms on the real emitter, h1-trainval (the held-out H1 read for idall3's October case; `work/score618.json`)
3 seeds, raw Δ. "ov" = 521's overall-equivalent INCLUDING reach. Under 619's amended rule, reach is unpriced for delivery
levers, so read nmae.

| contrast | nmae f 1 / 0.5 / 0.35 / cmp | reach f 1 / 0.5 / 0.35 / cmp | ov f 1 / 0.5 / 0.35 (incl. reach) |
|---|---|---|---|
| c4h − no hook | −0.00227 / −0.00074 / −0.00139 / −0.00042 | −0.0004 / +0.0005 / +0.0006 / +0.0009 | +0.0005 / +0.0003 / +0.0005 |
| idall3 − c4h | **−0.00244 / −0.00056 / −0.00074** / −0.00084 | −0.0001 / +0.0009 / **−0.0048** / −0.0045 | +0.0006 / +0.0003 / −0.0008 |
| idkg3 − c4h | −0.00116 / −0.00076 / −0.00073 / −0.00130 | +0.0040 / +0.0037 / −0.0006 / −0.0026 | +0.0012 / +0.0009 / −0.0000 |
| idall3 − no hook | −0.00471 / −0.00130 / −0.00213 | −0.0005 / +0.0014 / −0.0042 | +0.0011 / +0.0005 / −0.0003 |

- **nmae:** both identity forms beat the plain hook in every regime, and they tie at f 0.35 (−0.00074 vs −0.00073).
  - At 618's nmae × 1.0–1.2 (619), idall3 − c4h → board ≈ −0.0007 … −0.0009 raw at f 0.35, i.e. **+0.0002 … +0.00025
    overall from nmae**.
- **reach (unpriced, reported):**
  - idall3 drops reach at f 0.35 / cmp (−0.0048 / −0.0045, 3-seed sd of the overall 0.00006) but raises it at f 0.5.
  - idkg3 raises reach at f 1 / 0.5.
  - 619 showed reach mis-signed for a delivery pair on this surface, so this is not a price. But it is the one place where
    the two identity forms differ, and the board read of idall3 is what settles it.

## The two magnitude stacks (PRED-618-stacks.md, registered 14b5fda0): **PICK (b) 630 + DN** (`stacks618b.py`, `work/stacks618b.json`)
- **Runs:** 3 seeds, zero repair −5. Plain and ×1.25 dumps (`run618e tv stk_base / stk_x125 … --dump 126`).
  - KA: stk_base == the gate's G1 runs, 6 / 6 stats files bit-equal.
  - d4cd5f's rnup630 (UP pass; + DN for (b)) ran unchanged through `rnup618.py` (evaluate323's reference dir → ref_tv).
    Its own zs KA matched, and the UP pass moved 0 calls.
- **Pricing:** price630's own functions vs plain (`stacks618b.py`).
  - KA: (b) and 630-alone reproduce price618 exactly (max |diff| 0.0).
  - Why a second pass was needed: price630 refused ×1.25 and (a) on its "per-gene sums pinned" assert. **×1.25 is NOT
    bulk-neutral in this harness:** its per-gene block totals differ from plain in ~1.35M (target, gene) cells over all
    126 targets, by up to ~5e4 counts. 630 and DN do pin the sums.

| arm, overall-eq vs plain (3 seeds) | f 1 | f 0.5 | f 0.35 | cmp | **MIN(0.5, 0.35)** | nmae-only MIN |
|---|---|---|---|---|---|---|
| (a) ×1.25 + 630 | +0.00363 | +0.00360 | +0.00374 | +0.00522 | **+0.00360** | +0.00362 |
| **(b) 630 + DN** | +0.00306 | +0.00416 | +0.00404 | +0.00423 | **+0.00404** | +0.00404 |
| 630 alone (×1) | +0.00277 | +0.00328 | +0.00269 | +0.00296 | +0.00269 | |
| ×1.25 alone | +0.00085 | +0.00030 | +0.00092 | +0.00205 | +0.00030 | +0.00032 |

- **Decision (the coordinator's rule, mechanical):** MIN(a) − MIN(b) = **−0.00044**.
  - That is inside the 0.0005 band, which goes to (b), and (b) is higher anyway → **(b) 630 + DN**.
  - The nmae-only reading (619's amendment) gives the same pick (−0.00042).
- **The registered H1 prediction does not transfer.** H1 said (a) by +0.0014; 618 says (b) by 0.0004.
  - ×1.25's DE value is much smaller here: +0.0003 … +0.0009 at f 0.5 / 0.35, against H1's +0.0011 … +0.0019.
  - 630 + DN is about as large as on H1 (+0.0040 vs H1's +0.0030 … +0.0037).
- **Caveat for October:** ×1.25 also moves the bulk in this harness, and that pds / mse effect is NOT in these DE values.
  (b) is bulk-neutral by construction, so the pick does not depend on the unpriced bulk term unless ×1.25's bulk effect
  were large and positive.
- **Zero totals:** 0 in the plain-construction runs; 1 gated zero-total cell per ×1.25-based run (1 target), repaired by
  the −5 floor.
- **The bulk term, priced afterwards** (7772dc, 654 §3b, 6d563b29): ×1.25's pds / mse term, from 435's slopes on H1's
  harness bulk, is **+0.0005** (3 / 3 seeds positive). It is indicative only, because production re-levels the bulk per
  context.
  - Adding it to (a)'s 618 MIN gives ≈ +0.0041 vs (b)'s +0.00404. That is still inside the 0.0005 tie band, so **the pick
    stays (b)**.
  - The pick is robust to the unpriced bulk term at its H1 size.
