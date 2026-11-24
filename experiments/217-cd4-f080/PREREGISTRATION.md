# PREREGISTRATION — `cd4f080`, written 2026-09-04 00:40 UTC, BEFORE any build

Coordinator (lead). Committed before the builder runs. Not editable after the board result.

## 1. The arm — ONE variable against the current incumbent

```
VCC_CD4_TARGET_F : 0.34  ->  0.80
everything else  : cd4f034, unchanged (emitter rc_callset, VCC_ALPHA 0.982, operator e2f65ecf)
```

Incumbent = `cd4f034`, entry `QSe3JlNUYZ2BzvlcsL14`, board **0.2049195830**, rank 7.

**`f = 0.80` is not an extrapolation. It is 207's own measured PLOLO cosine maximiser** (+31.5%
against `f = 0.10`), which 207 declined to ship only to protect `pds` under a guard we now know was
calibrated on an understated instrument.

## 2. Why now, and why this size of step

`cd4f034` delivered **+0.0130557 against a +0.00437 prediction — 3.0×** — and the post-mortem
(`AGENTS.md` §8.4) locates the error in two systematic conventions, not in luck:

1. **A direction gain propagates to five members, not one.** The four we carried at ~zero delivered
   **60%** of the realised gain. `mse`-only pricing understates by **2.51×**.
2. **PLOLO understates the board's cosine response by 1.80×** (212: +13.83% at `f = 0.34`; board
   +24.9%, `cos` 0.07228 → **0.09027**).

So 207's sweep chose 0.34 with an instrument that understated the payoff by ~4.5× **and** understated
the `pds` cost by 2.0×. The payoff error is much the larger of the two.

## 3. Predictions

Scaling 207's PLOLO curve (`f = 0.34` → 0.80 is **2.28×** the cosine gain over the `f = 0.10` base)
by the board factor:

| quantity | prediction |
|---|---|
| board cosine | **+57%** over `f = 0.10`, i.e. `cos ≈ 0.113` (`cqawesome` 0.13211, `jw_13288` 0.13117) |
| the five positive members | +0.0144 × 2.28 = **+0.0328** over `rc_callset` |
| `pds` cost | −0.001345 × 2 to × 4 = **−0.0027 to −0.0054** over `rc_callset` |
| **net over `rc_callset`** | **+0.0274 to +0.0301** |
| **net over `cd4f034`** | **+0.0143 to +0.0170** |
| **predicted board overall** | **0.219 – 0.222** |

The leader is **0.22082**. This arm is preregistered as the one that could take first place, and it
is stated that way so the result is read against an ambitious claim rather than a safe one.

**Gate.** Under §8.4's corrected rule the `mse`-only price is ~+0.0131 and the six-member price
~+0.0143 to +0.0170, against a +0.0089 amplitude-changing floor. **Clears by 1.6–1.9×** even at the
worse floor. `VCC_ALPHA` stays **0.982** — 207 §6.2 measured α\* = 0.982 exactly and 211 re-solved
it at 0.9808/0.9879, all inside noise; **do not re-solve, that would make it two variables.**

## 4. Falsifiers — pre-committed

1. **raw `mse` does not improve on `cd4f034`'s 0.83635.** It is the load-bearing member; if a 2.35×
   dose increase does not move it, the dose-response has turned over between 0.34 and 0.80 and the
   whole ladder stops.
2. **raw `mse` outside [0.775, 0.815].** Point prediction ≈ 0.795.
3. **`pds` falls by more than 0.012 raw** (3.3× the realised −0.00365 at `f = 0.34`). 207's own guard
   put `f = 0.40` at 2.43σ, so a large `pds` cost here is the expected failure mode and this is the
   number that decides whether the ladder continues upward or stops.
4. **`reach` or `fid` falls.** Both rose at `f = 0.34` (+0.0264 / +0.0091 scaled). §8.4's mechanism
   says they track sign accuracy, so a fall would mean cosine did not actually rise and the `mse`
   move was amplitude.
5. **Config diff is not exactly 2 fields** against `210-cd4f034.h5ad.meta.json` (`cd4_target_f` and
   its deterministic `effects.cd4.fraction` consequence), with `operator.sha256`
   `e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123`, `alpha` 0.982,
   `fully_supported_targets` 276, `after_all_median_norm` unmoved.

## 5. What the slot buys even if it fails

The **third point on the board's own CD4 dose-response curve** (0.10, 0.34, 0.80). Two points give a
slope; three give the curvature and the turnover, which is what decides whether experiment 216's
extrapolation to `f ≈ 2.6` is worth any further slots at all. Either outcome is decisive.

## 6. Risk, stated plainly
This is a **2.35× step on the one parameter that has ever moved our direction channel**, taken on a
single board point's calibration. The conservative alternative was `f = 0.60`. I am taking 0.80
because it is 207's *measured* maximiser rather than an extrapolation past it, because the board has
just told us the instrument that chose 0.34 understated the payoff by ~4.5×, and because with ~96
slots remaining before Oct 22 the cost of one bad point is small against the cost of under-sampling
a lever this size. **If falsifier 3 fires, the next arm is `f = 0.60`, not `f = 1.2`.**

## 7. Discipline
Joke model name, **empty description** (`AGENTS.md` §14). This is the second of the two 2026-09-04
UTC slots.
