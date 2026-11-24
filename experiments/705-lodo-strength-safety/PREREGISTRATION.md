# 705: 522's θ 0.25 LODO strength rider, an A/B/C SAFETY read before it ships dark (preregistration)

Written 2026-10-15 ~05:05Z by the science session [bf003b], on the coordinator [400b6e]'s GO (10-17 s1). At this time the
rider dump exists (below), but no 705 build, package or entry does.

## The arm (one variable vs ALLSW, entry `4WNmsNfdQhYx6bfib0dR`, 0.2760851)
- **705lstr** = ALLSW's exact env, with `VCC_BULK_LOAD` = `work/dump_h1x2_cd4x15_lodo05hi6_s025.npz`
  (sha256 ee8ff715da5c76f614dba5c108a820b9e0a2f4bd19ae8e47ba86eadb59e750e5). Then re-levelled as 498's and 702's arms:
  r1 at ALLSW's factors → ntclevel417 → r2 at factors × EV (rho within 0.5 % of 0.9793; at most one r3).
- **The dump is exact** (`make_rider_dump705.py`, provenance `<dump>.prov705.json`). The builder writes the dump right
  after `apply_pert_scale`, so bulk = carriers × (w·c). The rider dump is carriers × (w′·c′), with w′ = 522's
  floor-preserving rider on 495's cap-6 key: w′ = max(0.5, w × exp(0.25 z_src)) above the floor, and w at the floor
  (artifact 392f1b93…, from `def_strength_key522.py`).
  Known answers, all PASS before the dump was written:
  - KA-sha: dump d58750b8 = ALLSW's VCC_BULK_LOAD, key dc374281, rider 392f1b93;
  - KA-cap6: carriers × (w₆·c₆) reproduces the cap-6 dump's bulk np.array_equal;
  - KA-rider: every key array but w equals dc374281's; def_strength_key522 `--ka` (KA-z / KA-inert / KA-arm) PASS today.
- What moves: per-target bulk magnitude only. Non-panel energy is preserved (78.39671 both). Row multipliers are 0.838 …
  1.075 (median 0.950): 54 rows up, 246 down. The carriers (the DE design) are untouched. The P2g tercile edges re-derive
  from v = ‖bulk‖/‖carriers‖, as they will on the day.
- Expected config diff vs ALLSW: r1 `bulk_load`, `bulk_load_sha256`, `bulk_sharpen_pt`; r2 adds `bulk_ctx_scale`
  (+ `de_sign_emitter_sha256` if the emitter pin moved). Built by `build705.sh` on the committed builder: no builder edit,
  after K-UP's build, from ~18:00Z 10-15.

## Class and the transfer evidence (AGENTS 11.1c)
**Transfer.** A per-target amplitude allocation from cross-line source strength. z_src is the donors' robust strength rank:
nothing from A/B/C enters, and it is computable for any panel from the donor banks plus `pert_counts.csv`.

Held-out reads (522; target bootstrap 2,000 draws, rider − percentile cap):

| surface | Δ direction | 95 % interval |
|---|---|---|
| efp (October-like) | +1.36 % | [+0.71, +2.28] |
| efr (October-like) | +1.81 % | [+0.96, +2.63] |
| real (October-like) | +0.71 % | [+0.36, +1.07] |
| A/B/C PLOLO folds | +0.07 % | [−0.35, +0.48] |

- 9/9 October-like folds > 0, with either 479 partner (S_P2 / P2g).
- Arc-like re-draws (essential-like removed, the strength mix matched to A/B/C): +1.58 / +1.22 / +0.72 points, 8/9 folds.
- The 702 pricing rules (FINDINGS §0):
  1. Every CI lower end is > 0 on the October-like panels.
  2. No seen-context surface is involved: the folds are leave-one-line-out on public lines.
  3. The change moves energy between targets, not between genes, so the instrument's gene axis sees it as it sees the base.
  4. The LODO family has board points (495, 498p2g; 495's instrument read was 15–40 % high but right in sign).
- October value ≈ +0.0004 … +0.0009 overall (522 §7). On A/B/C it is ≈ +0.0001, which is why this entry is a SAFETY
  read, not a gain read.

## Primary statistic, prediction, bins (board, vs ALLSW's own entry)
- **Primary: Δ(pds + mse) in overall units** = (scaled Δpds + scaled Δmse)/6. The direction members reproduce exactly
  (rule 13), so the bulk read is exact.
- **Prediction:** primary ≈ +0.0001, in −0.0003 … +0.0004 (the A/B/C fold interval through 518's identity). DE members
  move by re-roll only (σ ≈ 0.00025 overall; |Δ| ≲ 0.0005).
- **HARM, primary ≤ −0.0005: VETO.** The rider stays off for October. Record what the October-like instruments missed,
  as 702's post-mortem did.
- **Otherwise, primary > −0.0005: the rider is cleared for October** (`lodo_strength` θ 0.25 via DEF_DAY_LODO_STRENGTH_THETA).
  It ships only after [17adab] applies the unapplied pass-through patch (522 §8) with its known answers and one rehearsal.
  If that cannot happen before the freeze, it ships off (522: "ship the percentile cap alone").
- **Falsifiers (STOP and diagnose before reading):**
  - DE members (fid + jac + reach + nmae) move by more than 0.0015 overall in total (a bulk-only arm moves them only by
    re-roll);
  - a config diff beyond the expected keys;
  - a final rho outside 0.5 % of 0.9793.

## Tries
- Board, on this switch: 1, a safety read.
- Held-out: 522's θ grid (0.25 / 0.5 / 1 / 2 and the bounded-key forms) on 4 surfaces. θ 0.25 was chosen by minimax
  leave-one-panel-out with A/B/C included; the more aggressive θ that the October-like panels prefer were not chosen.

## What it does not test
The October gain itself, which is unreadable on A/B/C (+0.07 %). This entry can only veto.
