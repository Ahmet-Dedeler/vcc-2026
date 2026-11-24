# 654 — DE carrier amplitude on H1's real emitter: ×0.8 hurts, and under production's zero repair ×1.25 HELPS (+0.0019 overall-eq at f 0.35, 3/3 seeds)

2026-10-02 16:59–19:57Z (commit times), council session [7772dc]. Approved by the coordinator (range 654–655; low priority);
[d4cd5f] confirmed no overlap with v7.
- PREREGISTRATION b0bbfa20, before any run.
- Amendment A1 b1590abd, after seed 0's as-registered read and before any zero-repair run.

Class: **scorer mechanics** (the DE carriers' design amplitude; the bulk is fixed). Nothing is fitted to A/B/C, but the dose is
read on H1, so a held-out read is owed before October (11.1c.2).

## 0. Verdict
1. **Lowering the amplitude hurts.** ×0.8 vs ×1 is negative in every regime on every seed: f 0.35 −0.0028 overall-eq
   (−0.0036 / −0.0028 / −0.0020), nmae +0.0069 raw, reach −0.0042. **P2 (my bet that ×0.8 is better) is LOST; bin CLOSED
   below ×1.**
2. **The as-registered harness was not production-faithful, and it mattered for ×1.25.**
   - 537's command (which 654 copied) runs without production's zero repair (`VCC_DE_ZERO_REPAIR=-5`).
   - ×1.25 drives 10 gated genes per panel to a zero block total (delivered lfc ≈ −35), 3 of them inside the f 0.35 pools.
     Each costs ~+31 of error where production floors it near log2 −5.
   - As registered, ×1.25 reads −0.0033 (nmae +0.0147), all from those collapses.
3. **Production-faithful (A1, `--zero-repair -5`): ×1.25 BEATS ×1 at every board-like regime on all three seeds.**

| ×1.25 − ×1, overall-eq | f 1 | f 0.5 | f 0.35 | cmp |
|---|---:|---:|---:|---:|
| mean (seeds) | +0.0001 (+0.0013 / −0.0011 / +0.0003) | **+0.0011** (+0.0015 / +0.0006 / +0.0013) | **+0.0019** (+0.0011 / +0.0027 / +0.0019) | **+0.0012** (+0.0012 / +0.0012 / +0.0013) |
| members at f 0.35 (raw) | fid −0.00002 | jac −0.00006 | reach +0.0043 | nmae **−0.0042** |

4. **P1, retrodiction of 410 §7's "×1.25 is worse" on today's construction: FAILS.**
   - As registered, ×1.25 is worse (overall and nmae), but reach rises instead of falling.
   - Production-faithful, ×1.25 is better.
   - 410 §7 was measured 09-22, before zero repair existed (428, 09-23). Its "DE amplitude is not the lever" is overturned
     for today's production up to ×1.25.
5. **P3 ([d4cd5f]: the kf pins head DOWN) FAILS.** Amplitude reaches head DOWN almost proportionally: delivered |lfc| −16 %
   at ×0.8, +22 % at ×1.25 (f 0.35, zero repair). [d4cd5f] has corrected 631's mechanism note.

## 1. Where it acts (631's regions, the base's delivered direction, lib521.nmae_rows' conventions; zero repair; 3 seeds; f 0.35)
| region × dir | cells/target | nmae share ×0.8 / ×1 / ×1.25 | mean delivered \|lfc\| ×0.8 / ×1 / ×1.25 |
|---|---:|---|---|
| head DOWN | 17 | 0.0702 / 0.0666 / **0.0640** | 0.346 / 0.413 / 0.505 |
| head UP | 17 | 0.0268 / 0.0264 / 0.0266 | 0.297 / 0.345 / 0.399 |
| block DOWN | 143 | 0.2214 / 0.2208 / 0.2205 | 0.056 / 0.060 / 0.062 |
| block UP | 53 | 0.0662 / 0.0655 / 0.0652 | 0.096 / 0.107 / 0.114 |
| uncalled DOWN | 121 | 0.3272 / 0.3275 / 0.3277 | 0.117 / 0.127 / 0.135 |
| uncalled UP | 111 | 0.2402 / 0.2385 / **0.2372** | 0.132 / 0.142 / 0.145 |

- Region shares sum exactly to members()' nmae (asserted).
- The gain is head DOWN (−0.0026) plus uncalled UP (−0.0013), the two regions 631 found under-delivered.
- Without zero repair, ×1.25's head-DOWN share was 0.0829 (+0.0163 over ×1; 3 seeds). On seed 0 the loss sits in 9 cells
  pushed below lfc −3, 3 of them to a zero total (summed error +92). Every other head-DOWN cell improved (summed error
  −19.8).

## 2. What else moves
- **Calls per target:** 2,540 / 2,571 / 2,597 (×0.8 / ×1 / ×1.25). ×1.25 is a small CALL-SET change, not a pure delivery
  lever ([d4cd5f]), so 611's call-set rule prices it: fid ×2 (cost), jac ×0.5 (gain), reach not trusted on H1.
- **The bulk is not exactly pinned** (the prereg's premise was wrong). The ballast count K moves with the carriers' amplitude,
  so per-gene block totals change. Library totals change by up to 155 % per target, but the composition pds/mse read moves
  by only mean |Δ log1p CP10k| 2.3e-4. Harness raw pds +0.0007 … +0.0025 and mse +0.0006 at ×1.25 (f 1).
  - The prereg's rule therefore voids a clean "DE-only" reading. The deltas above are "carrier amplitude including its
    ballast-K drift", which is what the production knob does too (§3).
- **Zero-total cells:** ×1.25 + zero repair leaves 6 per panel. All are own-target rows (never repaired, never scored), 0 in
  any pool. Base has 1 (own target), so base_zr == base bit for bit on all 85 targets (A1's structural KA, 3/3 seeds).

## 3. The production knob, and what it would take for October
- **The knob already exists: `VCC_ALPHA`.** Production and def_day load the bulk (`VCC_BULK_LOAD`: the loaded bulk "rides the
  ballast while the carriers keep this build's design", builder 1139–1143), so α scales only the DE carriers and the
  ballast count. That is exactly this experiment's arm. ×1.25 corresponds to α 0.982 → 1.2275.
- **Board-scale price (611's rule, call-set change):** nmae −0.0042 raw × 1.2–1.6 gives ≈ +0.0014 … +0.0018 overall. fid ×2
  and jac ×0.5 cost ≈ −0.00003. Reach is unpriced (H1 says +0.004 raw). **Expected ≈ +0.0013 … +0.0018 overall on A/B/C.**
- **Before October (11.1c.2), in order:**
  1. 655's dose curve: ×1.5 / ×2.0 with zero repair. Running; its preregistered condition is met, 3/3 seeds.
  2. [d4cd5f]'s stacking test of the DOWN pass (UPDN) on ×1.25, from the dumps in `work/blocks/` (seq654c).
  3. A builder KA: `VCC_ALPHA` × s moves the carriers exactly as `eff_de` × s does here.
  4. A held-out cancer read with zero repair: 615's Orion raw rows or 614's lever rows. 631f found cancer surfaces want MORE
     DOWN magnitude than H1, so H1's dose is plausibly conservative.
  5. An A/B/C board read with its own prereg.
- **Tries:** 3 amplitudes (×0.8, ×1, ×1.25) × 2 constructions on H1. Dose points are 655's.

## 3b. Note added 10-03: "overall-eq" here is DE-only, and the bulk term is small
- Every overall-eq number in this file is lib521.overall_delta, i.e. fid / jac / reach / nmae. The harness's pds / mse drift
  (§2) is not in it.
- Priced with 435's slopes, ×1.25 zr adds +0.0005 (3/3 seeds positive) and ×0.8 adds +0.0001 (noisy). See 655's held-out
  section for the table, 618's read (×1.25 alone +0.0003 … +0.0009) and the October pick: (b) ×1 + 630 + DN; ×1.25 is the
  fallback.

## 4. Known answers and records
- **KA:**
  - base s0 == 537's `h1_cpm` stats bit for bit (q, p, lfc, sums; pre-zs and zs);
  - base s1 == 630's stored stats bit for bit;
  - base_zr == base on 85/85 targets, seeds 0–2 (A1's structural KA).
- **Harness:** `run654_496.py` = run496 + one design branch (`x080`); `run654.py` = run537 importing it. Diffs are in `work/`.
- **Analysis fix mid-run (de785273, before the final read):** decomp631's estimator keeps targets with < 10 pool genes and
  does not sum to lib521's nmae. `analyze654.py` now uses lib521.nmae_rows' conventions, with the sum asserted.
- **Outputs:** `work/analyze654.{out,json}` (as registered), `work/analyze654zr.{out,json}` (production-faithful),
  `work/seq654{,b,c}.log`. The stats npz files stay local (`work/h1_cpm/`, ~36 MB each). `zerocount654.py` is the shared
  zero-total counter (also used by [21dd4c]'s 615).
- **Gate:** every run went through gate602, one ticket per run. Peak RSS was 5.8–7.6 GB.

## 5. Lessons
- **A harness that omits a production repair step misprices exactly the levers that trigger the repair.** 537's command lacks
  `--zero-repair -5`. Any DOWN-magnitude lever that can create a zero total must run with it. [d4cd5f] checked DN / UPDN: no
  zero totals, prices stand. [21dd4c] checked 614: its lever rows already pass it. Count lfc < −30 cells
  (`zerocount654.py`) before trusting a DOWN-magnitude price.
- **Region shares must use the scorer-faithful estimator.** Otherwise they do not add up to the member being explained.
