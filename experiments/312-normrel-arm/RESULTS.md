# 312 / 312b — Per-target donor weighting on the real board: the dose curve, and the arm that produced our best entry

> **September 9 review:** this is the historical plain-dose result; 319-cap is
> now best. The 0.806–2.155 / +6% values below came from a synthetic fixture,
> not the shipped panel. Experiment 318 §1 supplies the delivered values and
> shows that 316 is norm-matched before whitening, not after emission.
> Therefore “the distortion was worth having” is not a pure magnitude
> attribution. Published scores here remain valid.

**Written by the coordinator** (this arm was built and submitted directly, not by a research subagent, so this file
is the record). **Two slots of 2026-09-07, one preregistered dose curve, three board points on ONE shipped
implementation.**

## The line
**The arm pays at exp −1.5 (+0.0016245, our first gain since 09-05) and the curve TURNS OVER by exp −3.0. The
exponent axis is closed at its peak ≈ −1.5.**

| exp | entry | model name | overall | vs `cd4svd5` |
|---|---|---|---:|---:|
| 0 | `Xv0ayvuoHpsKYWv2iu42` | `cd4svd5` (incumbent) | 0.2160944 | — |
| **−1.5** | **`Q519I2pLWRCW7PAEx5ul`** | `the cat sat on the modem` | **0.2177189** | **+0.0016245** |
| −3.0 | `1n5LFk6PSEHq15UHENIE` | `a goose in the server room` | 0.2164511 | +0.0003567 |

All three records captured under `data/board-entries/`. Preregistrations: `PREREGISTRATION.md` (hashed 644a43d4…,
frozen before the first entry) and `PREREGISTRATION-dose2.md` (73037e9b…). Dose 2 landed in its preregistered
**Δ₂ ≤ −0.0008** bin ("the curve has turned over"); dose 1 in its **+0.0010 ≤ Δ < +0.0040** bin.

| member | exp 0 | exp −1.5 | exp −3.0 | Δ(−1.5→−3) |
|---|---:|---:|---:|---:|
| `pds` | 0.77809 | 0.79039 | 0.79241 | **+0.00202** |
| `mse` | 0.17477 | 0.18013 | 0.17634 | −0.00379 |
| `nmae` | 0.11047 | 0.10277 | 0.09980 | −0.00296 |
| `reach` | 0.22139 | 0.22026 | 0.21791 | −0.00235 |

**`pds` rises monotonically with |exp| while every other member falls past −1.5** — 143's pattern reproduced on the
real scorer. Note w ∝ r^−3 is the **square** of the r^−1.5 weight ratios, not twice the dose.

## What was built
`VCC_NORMREL_EXP` (`scripts/make_submission_pm_cd4.py`, `apply_patch.py`): per-target donor weights
w_j(t) ∝ (‖row_j(t)‖ / that donor's own panel median)^exp, renormalised over the donors carrying t, the row then
rescaled to the uniform-weight row's norm. **Default unset = the original expression, so a normal build is
bit-identical.** 097 route (a) claimed the builder already took a per-target weight dict — **that was wrong**
(`weights` is `Mapping[str, float]`, one scalar per donor) — but the wiring is still small because
`build_v9_pre_whiten_effects` already loops per target, so `MultiDonorModel` needed no change.

**Known-answer rows (`test312.py`, 4/4 PASS, re-verified against the committed builder):** exponent 0 == unset at
**0.000e+00** (so an unset build cannot change); row norms preserved at **2.05e-16**; the donor whose row is 5× its
own panel median is down-weighted 0.452 → 0.059. KA-3b first read 1.36e-07; a float64 control proved the logic
**exact** (error 0.0) and the residual purely float32 representation, so the norm-match was moved to float64 rather
than the threshold loosened.

## The two things this arm turned out NOT to be — both established after the fact
1. **It is not the object 308 measured** (313). Production applies the PerturbMatch correction at the **uniform**
   share and **outside** the norm match, because the patch touched only the pre-whiten function.
2. **It is not shape-only as delivered** (313, 316). Delivered per-target row norms run **0.806–2.155** (whole
   matrix +6.0 %). The 2.05e-16 identity above holds at the **pre-whitening matrix only**.

**And the distortion was worth having.** 316 built the exactly-shape-only version: **0.2169733**, i.e. the
weighting alone is **+0.00088** and the accidental per-target norm distortion adds **+0.00075** on top. `nmae`
recovered +0.00339 under parity exactly as 313 predicted, and was outweighed by `mse` +0.00399 and `pds` +0.00225
that the distortion had been buying — consistent with 317's board-wide measurement that **buying `mse` is what wins
here** (83.2 % of teams that bought it improved). Whether that per-target magnitude effect is a deliberately
reachable rule is 318's question.

## Files
`apply_patch.py` (idempotent), `apply_cap_patch.py` (Astra's `max(1,r)^−p` variant, **staged unapplied**),
`test312.py`, `run312.sh` / `run312b.sh` (gated two-invocation builds), `PREREGISTRATION*.md`,
`work/build*/` (metas, driver logs), `work/submit312*.log`.
