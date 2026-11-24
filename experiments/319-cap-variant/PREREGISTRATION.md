# 319 — Preregistration, frozen BEFORE the entry exists

**The arm.** `VCC_NORMREL_EXP=-1.5 VCC_NORMREL_MODE=cap`, package `submissions/319-normrel15-cap.vcc`.
Config differs from the winning entry `Q519I2pLWRCW7PAEx5ul` (**0.21771888483638538**) in **exactly one key**,
`normrel_mode` — verified. **Reference for Δ is that entry.** Second slot of 2026-09-08.

## The hypothesis being separated
`r^exp` does two things at once: it **suppresses** donor rows above their own panel median AND **up-weights** rows
below it. `cap` clamps `r ← max(1, r)`, so only oversized rows are touched. The dose curve motivates it: past
−1.5 `pds` keeps rising (+0.00202 at −3.0) while `mse` −0.00379, `nmae` −0.00296 and `reach` −0.00235 all fall —
143's precedent ("`pds` up, signed credit down"), which is what amplifying small noisy rows would look like.

**Measured, so this is a real arm and not a near-incumbent:** the cap touches **49.94 %** of the real panel's
(target, donor) cells and retains **52.5 %** of the shipped arm's median weight move away from uniform. Against
the shipped arm the delivered rows have median cosine **0.999731** but a minimum of **0.065347** — most targets
barely move and a few move a great deal, which is what a rule that only acts above the median should do.

## Bins, fixed now (Δ vs 0.21771888)
| bin | reading |
|---|---|
| **Δ ≥ +0.0005** | **suppression alone beats the combined rule** — the up-weighting half was actively harmful. The dose axis reopens *for the cap form only*, and the next arm pushes `cap` to a larger \|exp\| where `plain` turned over. |
| **−0.0005 … +0.0005** | the two halves are not separable at this dose. The family is closed at ≈ +0.0016 and the weighting story is finished; no further variant earns a slot. |
| **Δ ≤ −0.0005** | **the up-weighting half was CONTRIBUTING**, not damaging. The combined rule is right as shipped, the 143 pattern at high dose has some other cause, and the family is closed. |

**Declared secondary, no verdict:** `pds` and `nmae` specifically — if the cap's story is right, it should keep
more of the shipped arm's `pds` gain while giving back less `nmae`. And whether Δ ≈ half the shipped arm's
+0.0016245, which is what a linear response to a 52.5 %-sized perturbation would give (that would say the rule is
simply dose, with no mechanism distinction at all).

## Stated in advance
- Expected magnitude is **small either way**; this does not approach the **+0.0296** gap to first and will not
  move us off 4th. It is spent because the slot is prize-irrelevant, expires tonight, and this is the last
  untested hypothesis in the family.
- The board shows the LATEST entry: our displayed row is currently 0.2169733 (the parity arm). **If Δ < 0 it falls
  further. Not to be re-submitted to restore rank.**
- `Q519I2pLWRCW7PAEx5ul`, `1n5LFk6PSEHq15UHENIE` and `2K5ppxkowMSj2ugoZqYb` are already captured.
- **This is one variable.** 318 separately closed the per-target *magnitude* axis (NO-GO, beaten by shuffles of its
  own multipliers), so a magnitude explanation for any movement here is already excluded.
