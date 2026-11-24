# 319 — Build and verify Astra's cap variant `VCC_NORMREL_MODE=cap`: suppress OVERSIZED donor rows without up-weighting tiny ones. One variable against the winning arm. STOP before submitting.

You are a research subagent. Work in `experiments/319-cap-variant/`. You **are** authorised to edit
`scripts/make_submission_pm_cd4.py` and `scripts/def_pipeline.py` — those two only, plus your own directory.
**NEVER run `vcc submit`.** You MAY run `vcc prep`. **One** full build. `uv run python`. Heavy build serial under
the memory gate. **Time box: 4 hours.**

## Why — this is the last untested hypothesis in the family, and the dose curve motivates it
`VCC_NORMREL_EXP` sets w_j(t) ∝ r^exp with r = ‖row_j(t)‖ / that donor's own panel median. **That rule does two
things at once: it SUPPRESSES donor rows above the median AND UP-WEIGHTS rows below it.** Board evidence for the
family, all on one implementation:

| arm | overall | vs `cd4svd5` |
|---|---:|---:|
| `cd4svd5` incumbent | 0.2160944 | — |
| exp −1.5 (shipped) | **0.2177189** | **+0.0016245** |
| exp −3.0 | 0.2164511 | +0.0003567 |
| exp −1.5 + `VCC_NORMREL_PM=blend` (parity) | 0.2169733 | +0.0008789 |

**The curve turns over past −1.5 with `pds` still RISING (+0.00202) while `mse`, `nmae` and `reach` all fall** —
143's precedent is exactly "`pds` up, signed credit down across every family". **The hypothesis: the damage at high
dose comes from the UP-WEIGHTING half (amplifying small, noisy donor rows), and suppression alone would let the
dose go further without the collateral.** 318 separately closed the per-target *magnitude* axis (NO-GO, not
separable from its own shuffles), so this is a **weighting** hypothesis, not a magnitude one.

## The change — already written, apply it
`experiments/312-normrel-arm/apply_cap_patch.py` is staged and idempotent: it adds `VCC_NORMREL_MODE` ∈
{`plain`, `cap`}, default **`plain`** (today's behaviour), and under `cap` clamps the ratio with
`ratio = max(1.0, ratio)` before the exponent, so a donor at or below its own panel median keeps weight 1 and only
oversized rows are down-weighted. It also records `normrel_mode` in the meta config under the same
"written only when set" convention. **Read it, apply it, and verify it does what this paragraph says.** Register
`VCC_NORMREL_MODE` in `def_pipeline.py`'s `ARM_BOARD_DELTA` exactly as 312/316 registered theirs.

## Known-answer rows — ALL must pass before you build
1. **`VCC_NORMREL_MODE` unset (or `plain`) with exp −1.5 reproduces the CURRENT builder's arm BITWISE** — compare
   against a byte copy of the builder taken before your patch, loaded as an independent module (316's `test316.py`
   shows the pattern; reuse it).
2. **`VCC_NORMREL_MODE=cap` with exp 0 reproduces `cd4svd5` bitwise.**
3. **The cap actually caps:** on 313's fixture, every donor with r ≤ 1 has weight exactly equal to its uniform
   weight under `cap`, and every donor with r > 1 is down-weighted; report the share of (target, donor) cells the
   cap touches on the fixture **and, if cheaply available, on the real panel** — if the cap touches almost nothing,
   say so, because then the arm is nearly the incumbent and that is the result.
**If any row fails, STOP and report. Do not build.**

## Then
One full build on the real A/B/C bundle with **`VCC_NORMREL_EXP=-1.5 VCC_NORMREL_MODE=cap`** — the SAME exponent
as the winning arm, so the cap is the ONLY variable against entry `Q519I2pLWRCW7PAEx5ul`. Follow
`experiments/312-normrel-arm/run312.sh`'s gated two-invocation pattern (expect one postbuild config-drift FAIL
naming the new key). Verify `config.normrel_exp` **and** `config.normrel_mode` landed in the meta, and that the
config differs from `submissions/312-normrel15.h5ad.meta.json` in **exactly one key**. Then `vcc prep` to
`submissions/319-normrel15-cap.vcc`, copy the meta beside it, run `scripts/check_submission_novelty.py check`.
**STOP THERE** — the coordinator preregisters the reading and submits.

## Rules
**Do not quote any number as "delivered" that was measured on 313's synthetic fixture** — that error was made
twice in the last day (Correction 5: the fixture's 0.806–2.155 / +6.0 % against the real panel's 0.9321–1.2349 /
+0.85 %). If you report a magnitude or norm, **name the stage and the surface it was measured on.** Programme
rules 1–9 bind.

## Deliverable
`RESULTS.md`: the diff; the three known-answer rows with measured numbers; the share of cells the cap touches; the
build's postbuild line; the one-key config diff; prep and novelty output; and one line —
**319 (`normrel15-cap`): BUILT AND VERIFIED / BLOCKED (KA row N failed: …)**. Then stop.
