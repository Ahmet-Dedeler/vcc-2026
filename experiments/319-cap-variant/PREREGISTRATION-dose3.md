# 319d45 — Preregistration for the cap form's THIRD dose (−4.5), frozen before any dose-2 or dose-3 entry exists

Written 2026-09-09 22:41Z, while `319d3-cap30.vcc` (dose −3.0) is uploading as the first slot of the day
and before its score exists. The build below is **speculative**: it starts now because build + prep take
~60 min and the day's second slot expires at 00:00Z. Whether the package is ever uploaded is decided by
the rule in §2, not by whoever is awake at 23:45Z.

**The arm.** `VCC_NORMREL_EXP=-4.5 VCC_NORMREL_MODE=cap`, `VCC_NORMREL_PM` unset. One key away from
319d3 (`normrel_exp` −3.0 → −4.5) and from 319 (−1.5 → −4.5). Builder and driver must hash identically
to the 319d3 run (`1318cd3c…`, `3fb0c4bc…`); the run script aborts otherwise. Package
`submissions/319d45-cap45.vcc`.

## 1. Notation
- S1 = 0.21865973944346428 (cap −1.5, entry `xuDelEik7zGWMXdNu0R4`)
- S3 = the cap −3.0 score, unknown at freeze time
- S45 = this entry's score

## 2. Upload rule, fixed now
Upload **only if S3 − S1 ≥ +0.0005**, i.e. dose 2 lands in its own preregistered bin 1 ("the cap form's
dose axis is genuinely open; push further (−4.5)"). In every other dose-2 bin the frozen dose-2 reading
already says *no further dose earns a slot*; the package is kept, unsubmitted, and the slot expires. An
expiring slot is not a reason to override a frozen reading.

## 3. Bins for Δ = S45 − S3, fixed now
| bin | reading |
|---|---|
| **Δ ≥ +0.0005** | the cap axis is still open past −4.5. Reconcile with 308's hard-gate loss before any −6.0: cap at large dose converges on a median hard gate, so a gain here means the *soft* gate matters, not the limit. |
| **−0.0005 … +0.0005** | plateau from −3.0 to −4.5. Family closes at max(S3, S45). Ship −3.0 (fewer suppressed rows) unless the secondaries below favor −4.5 on `nmae`. |
| **Δ ≤ −0.0005** | turnover located between −3.0 and −4.5; the cap optimum is ≈ −3.0. Family closes there. |

**Secondaries, no verdict:** `pds`, `nmae`, raw `mse` (326 measured cap30's raw-MSE cost at +0.0152 on H1;
the expectation is that the cost grows with dose and the board's unclamped `mse` member pays it).

## 4. Stated in advance
- The board shows the LATEST entry. If this uploads and loses, the displayed row falls below both 319 and
  319d3. That is costless and is **not** to be corrected by re-submitting.
- If prep overruns 00:00Z and S3 was in bin 1, the upload proceeds as the next day's first slot; the bins
  above still apply.
- Nothing here changes production defaults. `VCC_NORMREL_EXP` remains a build-time switch.

## 5. Pre-outcome amendment, 2026-09-09 22:46Z — dose is −4.0, not −4.5
The first launch (`work/build_d45/`, 22:38Z) stopped at `repin-frozen`: the byte-identical builder rejects
`VCC_NORMREL_EXP out of range [-4,4]`. The builder is not edited (that would break the one-variable claim), so
dose 3 is **−4.0**, the largest magnitude the builder accepts. Script `run319d4.sh`, package
`submissions/319d4-cap40.vcc`, model name "a raccoon audited the pantry". No dose-2 or dose-3 score exists at
this time. §2 upload rule and §3 bins are unchanged; every "−4.5" above reads −4.0.
