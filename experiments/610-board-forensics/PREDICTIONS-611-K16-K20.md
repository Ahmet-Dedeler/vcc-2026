# 611's registered board predictions for K16 (10-06) and K20 (10-07): the first out-of-sample tests of its call-set multipliers

Written 2026-10-02 ~13:45Z by [21dd4c], BEFORE either entry exists. Arms:
- K16 = 446bj2zk16 (`VCC_RC_CALLSET_CPM=1200`: 400 head + 1,200 CPM block);
- K20 = 446bj2zk20 (`VCC_RC_CALLSET_CPM=1600`).
Both read against 446bj2z. The slot decisions are 523's / 537's (their bins stand). This card tests the INSTRUMENT,
611's multipliers.

## Inputs: H1 real emitter, 3 seeds (537 `price537.json`, P3 = K 1,600 alone, P2 = K 2,000 alone), raw member Δ
| arm | regime | fid | jac | reach | nmae |
|---|---|---:|---:|---:|---:|
| K16 | f 1 | −0.0185 | −0.0065 | +0.0013 | −0.0001 |
| K16 | f 0.5 | +0.0005 | −0.0001 | −0.0046 | +0.0001 |
| K16 | **f 0.35** | **+0.0034** | **+0.0019** | −0.0025 | −0.0011 |
| K20 | f 1 | −0.0100 | −0.0029 | −0.0020 | +0.0003 |
| K20 | f 0.5 | +0.0009 | +0.0013 | −0.0045 | −0.0003 |
| K20 | **f 0.35** | **+0.0035** | **+0.0024** | −0.0083 | −0.0022 |

## Registered predictions (611: f 0.35; fid ×~2 and jac ×~0.5 for call-set changes; nmae ×1.2–1.6; reach not priced)
| arm | Δfid raw | Δjac raw | Δnmae raw | reach band (raw) | Δ(fid + jac) overall |
|---|---|---|---|---|---|
| K16 | **+0.0069** (×1–2.8: +0.0034 … +0.0097) | **+0.0010** (×0.43–1: +0.0008 … +0.0019) | −0.0013 … −0.0017 | ±0.0065 | **+0.0043** (×1: +0.0028) |
| K20 | **+0.0070** (+0.0035 … +0.0098) | **+0.0012** (+0.0010 … +0.0024) | −0.0026 … −0.0035 | ±0.0123 | **+0.0044** (×1: +0.0030) |

## How the reads grade the instrument
- **Regime.** fid's sign flips between H1's f 1 (−0.010 … −0.019, the yield penalty when n_pred < n_real) and f 0.35.
  - A board Δfid > 0 confirms that the board's n_real sits in the f 0.35–0.5 range (192 §4.2).
  - A board Δfid < 0 means the board is f 1-like for call-set SIZE. Then every f 0.35 call-set price is wrong-signed.
- **Multipliers.** Board Δfid / H1 f 0.35 Δfid:
  - ≈ 2 → the ×2 holds beyond the one composition pair (cp);
  - ≈ 1 → size changes do not inherit cp's ×2.
  - Likewise for jac, against ×0.5.
- **Reach** is read inside its band only. 611 never prices reach for call-set changes.
- After each read, `desk611.py` takes the pair in one line. The multiplier table in RESULTS-611.md is then re-fitted
  on 3 call-set pairs instead of 1.

## GRADE (2026-10-04, after K20 read; appended by [21dd4c] at the coordinator's request): K20 **FAIL, wrong sign**
- **The read:** K20 = entry `kJfe2dP9LEJawDEf0pMm`, 10-04 s2, vs 446bj2z.
  - Board raw Δ: fid **−0.00169**, jac +0.000185, nmae −0.000059, reach +0.00060.
  - Δ(fid + jac) **−0.00086** overall. pds / mse exactly 0.
  - Record: `523-callset-size/BOARD-RESULT-446bj2zk20.md`, 59543bfe.
- **Against the card:** it said Δfid +0.0070 (+0.0035 … +0.0098), Δjac +0.0012, Δ(fid + jac) +0.0044.
  - fid has the WRONG SIGN. jac is ~6× too large. The card's nmae (−0.0026 … −0.0035) missed the board's ≈ 0 by ~40×.
  - By the card's own regime rule, the board is not f 0.35-like for call-set size, so every f 0.35 call-set price is void.
- **The multiplier question did not get as far as being asked:** the regime failed first.
- **K16's prediction above is WITHDRAWN.** It is superseded by `experiments/662-callset-regime/PRED-662-K16.md`,
  registered before any 10-06 entry.
- **Which level the board matches (662, `fit662.py`): f ≈ 0.61 (0.60–0.62), not f 1.**
  - f 1 matches K20's fid sign only: it over-reads it 6×, mis-signs its jac, and needs fid multipliers 17× apart for cp and
    K20.
  - At f 0.61 one fid multiplier (×2.69) fits cp and K20, their jac ratios agree (0.43 / 0.49), and every priced sign of the
    four board DE pairs is right.
- **Replacing 611's "re-fit on 3 pairs":** `desk611.py` is not re-run. 662's I-662 (fitted on cp + K20) is the re-fit, and
  K16 is its first out-of-sample test.
  - I-662 predicts K16 Δfid −0.0116 [−0.0130, −0.0107], Δ(fid + jac) −0.0067 [−0.0075, −0.0062].
  - The f 1 reading predicts −0.0031 / −0.0030.
