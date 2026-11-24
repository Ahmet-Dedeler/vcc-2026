# PRED-661-idall3: the ungated v3 identity hook (idall3) on the best base, as its prereg's arm — 618-registered prediction

[21dd4c], 2026-10-03 ~19:05Z, committed BEFORE idall3's board read. Coordinator's 661 item (3). The arm, the primary
statistic and the bins are those of `612/PREREGISTRATION-draft-IDALL3.md` (FINAL, unchanged). This file only adds 618's
prediction beside H1's.

## The statistic
- Δ(fid + nmae) overall of arm − base, read against Δ_496a (the plain CLEAN4 hook's own effect).
- That is the board analog of idall3 − c4h. WIN means ≥ Δ_496a + 0.0002 with Δreach inside ±0.0042 raw; TIE means
  within ±0.0002.

## Source (618, real production emitter, 3 seeds; `618/work/score618.json`, "idall3-c4h")
| 618 raw Δ, idall3 − c4h | f 1 | f 0.5 | **f 0.35** | cmp |
|---|---|---|---|---|
| nmae | −0.00244 | −0.00056 | **−0.00074** | −0.00084 |
| fid | −0.00010 | −0.00010 | **−0.00020** | −0.00020 |
| jac | +0.00005 | −0.00004 | **−0.00008** | −0.00009 |
| reach (unpriced, 619) | −0.00007 | +0.00094 | **−0.00476** | −0.00453 |
| Δ(fid + nmae) overall | +0.00061 | +0.00010 | **+0.00009** | +0.00012 |

H1 (h1-test, the prereg's evidence) had idall3 − c4h nmae −0.00186 / −0.00218 / −0.00223, i.e. about 3× 618's at
f 0.35.

## Board prediction (619: nmae × 1.0–1.2; the fid of a call-set side effect × 0.8–2)
- **At f 0.35 (the board-like regime):**
  - nmae −0.00074 × 1.0–1.2 → +0.00020 … +0.00024;
  - fid −0.00020 × 0.8–2 → −0.00009 … −0.00022;
  - **Δ(fid + nmae) − Δ_496a ∈ [−0.00002, +0.00015], point +0.00007** (×1.1, ×1.4).
- **At f 0.5:** [+0.00004, +0.00014].
- **At f 1:** [+0.00056, +0.00075] (a WIN).
- **Registered:** **TIE** (|Δ_arm − Δ_496a| < 0.0002) at the board-like regimes. 618 does NOT predict the WIN that H1
  predicted (+0.0007 … +0.0010). A WIN would require the board to behave like f 1.
- **Reach risk, stated:** 618's reach for idall3 − c4h at f 0.35 is −0.0048, just outside the prereg's ±0.0042 band.
  - 619 made reach unpriced for delivery levers on 618 (the MDE pair was mis-signed), so this is a risk flag, not a
    price.
  - But a board WIN on fid + nmae could still fail the WIN bin's reach condition.
- **Members:**
  - **fixed:** pds and mse (a DE-only, prior-class arm; bit-identical by construction);
  - **~fixed:** jac (|Δ| ≤ 0.0001 on 618);
  - **moving:** nmae (gain), fid (a small cost from hook-induced call changes), reach (unpriced, risk).
- **What the read then tells us:**
  - a TIE confirms 618 over H1 for this lever;
  - a WIN with reach in band says h1-test was right and 618 under-reads identity hooks;
  - a WIN with reach out of band is the risk above, realised.

## I-662 ALTERNATIVE (added 2026-10-04 ~17:55Z, BEFORE idall3's read; the TIE registered above stands)
- **Why:** K20's read put the board's level for call-set changes at f ≈ 0.61 (662).
- **618 idall3 − c4h at f\* 0.61** (`alt662.json`; known answer vs score618 exact): nmae −0.00095, fid −0.00013, jac −0.00001,
  reach −0.00011.
  - Band f 0.60–0.62: nmae −0.00073 … −0.00095, fid −0.00010 … −0.00013.
- **Multipliers:** nmae × 1.05–1.30 (618 at f\*); fid × 0.8–2.7 (I-662's call-set fid range, both panels).
- **I-662 prediction:** Δ(fid + nmae) − Δ_496a = **+0.00020 (+0.00006 … +0.00028)**.
  - **Borderline TIE / WIN:** the point sits on the WIN bar.
  - The reach risk flagged above is gone at f\* (−0.0001, well inside ±0.0042).
- **Reading:**
  - A clean TIE (< +0.0002) with reach in band fits both lines.
  - A WIN of +0.0002 … +0.0003 with reach in band favours I-662's level.
  - Reach outside ±0.0042 favours the f 0.35 line.
