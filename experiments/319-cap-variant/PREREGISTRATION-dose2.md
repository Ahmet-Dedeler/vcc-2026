# 319d3 — Preregistration for the cap form's SECOND dose, frozen before the entry exists

**The arm.** `VCC_NORMREL_EXP=-3.0 VCC_NORMREL_MODE=cap`, package `submissions/319d3-cap30.vcc`. Config differs
from the winning entry `xuDelEik7zGWMXdNu0R4` (**0.21865973944346428**, exp −1.5 cap) in **exactly one key**,
`normrel_exp` −1.5 → −3.0 — verified against that package's own meta. **Reference for Δ is that entry.**
First slot of 2026-09-09.

## Why this is the decisive dose point
The `plain` rule peaked at −1.5 (0.2177189) and **collapsed to 0.2164511 by −3.0** — `pds` still rising (+0.00202)
while `mse`, `nmae`, `reach` all fell. 319 showed the collapse's likely cause: `plain` also **up-weights**
below-median rows, and removing that (`cap`) gains **+0.00094** at −1.5. **If the up-weighting half caused the
turnover, the cap form should hold up at −3.0 where plain fell apart.**

| | exp −1.5 | exp −3.0 |
|---|---:|---:|
| `plain` | 0.2177189 | 0.2164511 (**−0.00127**) |
| `cap` | **0.2186597** | *this entry* |

## Bins, fixed now (Δ vs 0.21865974)
| bin | reading |
|---|---|
| **Δ ≥ +0.0005** | the cap form's dose axis is genuinely open where `plain`'s was not. Push further (−4.5) and re-price the family; this is the first arm in the programme to gain twice on the same axis. |
| **−0.0005 … +0.0005** | flat: the cap removed `plain`'s turnover but adds nothing beyond −1.5. **The family closes at +0.0026**; no further dose earns a slot. |
| **Δ ≤ −0.0005**, but **> plain's −0.00127** | the cap *softens* the turnover without removing it — informative, and the family closes at −1.5. |
| **Δ ≤ plain's −0.00127** | the turnover is not caused by the up-weighting half at all; 319's +0.00094 has some other explanation and the mechanism story needs re-opening before any further arm. |

**Declared secondary, no verdict:** `pds` (rose +0.00202 for plain at this dose — does cap keep it?) and `nmae`
(fell −0.00296 for plain — does cap protect it, as it did at −1.5 with +0.00475?).

## Stated in advance
- **A ceiling is expected somewhere.** Pushed far enough, `cap` converges on a hard gate at the median, and 308
  measured hard gates losing badly (`gate:rel` −0.00863). **−3.0 may already be past the cap's own optimum**, and
  that is a real possibility rather than a hedge.
- The board shows the LATEST entry. Our displayed row is currently our best (0.2186597); **if Δ < 0 it falls.
  Not to be re-submitted to restore rank.**
- One variable. 318 closed the per-target *magnitude* axis, so a magnitude explanation is already excluded.
