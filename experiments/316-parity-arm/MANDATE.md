# 316 — Build the parity-correct arm `normrel15-pm` (`VCC_NORMREL_PM=blend`) to 313's specification, prove it with three known-answer rows, and STOP before submitting.

You are a research subagent. Work in `experiments/316-parity-arm/`. You **are** authorised to edit
`scripts/make_submission_pm_cd4.py` and `scripts/def_pipeline.py` — those two files only, plus your own directory.
**NEVER run `vcc submit`.** You MAY run `vcc prep`. You may run ONE full build. Nothing under `submissions/` except
the prep output named below. `uv run python`. **Time box: 4 hours.** Heavy build serial under the memory gate.

## Why
Entry `Q519I2pLWRCW7PAEx5ul` (`VCC_NORMREL_EXP=-1.5`) scored **0.21771888483638538**, Δ **+0.0016245** over
`cd4svd5` — a real gain. But **313 proved the shipped object is not the one 308 measured**: PerturbMatch is applied
at the uniform donor share and OUTSIDE the norm match, so the delivered rows are amplitude-changing (per-target
delivered norm **0.8062–2.1549** on scored columns, whole-matrix **+6.00 %**) where 308's arm is exactly shape-only
(ratio 1.000000). We preregistered shape-only and predicted an `nmae` **benefit** of −0.00095; the board returned
**−0.00771**. The dose curve then turned over at −1.5 (312b: exp −3.0 gives 0.2164511, `pds` still rising but
`mse`/`nmae`/`reach` all falling), so **the exponent axis is closed and the remaining lever is the SHAPE of the
object, not its size.** Restoring true shape-only is the best-motivated arm on the table: if it recovers most of
the −0.00771 `nmae` leg at the same weighting, it lands near **+0.0029**.

## READ FIRST
`experiments/313-normrel-parity/RESULTS.md` — **especially §8, which is your specification**; its fixture
`work/fixture313.npz` and `work/parity313.json`. Then `experiments/312-normrel-arm/` (`apply_patch.py`,
`test312.py`, `PREREGISTRATION.md`), `docs/READ-FIRST.md` §2 and §7, `AGENTS.md` §8.1d and §11.2b.

## The change — implement 313 §8 exactly, do not improvise
`VCC_NORMREL_PM` = `blend`; **default unset = today's shipped behaviour**. Setting it without `VCC_NORMREL_EXP`
must be a **hard error**. In `_normrel_pre_whiten_effects`, per target `t` carried by both K562-GWPS and
PerturbMatch, compute in ORD_A's order:

```
corr(t)  = alpha * effect_scale * PM_BLEND * (pm_standardized(t) - incumbent_k562(t))   # on pm_target_columns
base_pm  = effect(t, uniform)         + (1/|have|) * corr(t)     # blended UNIFORM row = the norm reference
cand_pm  = effect(t, normrel weights) + w_k562(t)  * corr(t)     # blended WEIGHTED row
row(t)   = cand_pm * ||base_pm|| / ||cand_pm||                   # norm match LAST -> exactly shape-only
```

Two supporting edits, both mechanical: (i) `build_production_effect_matrix` must **not** apply
`apply_perturbmatch_k562_replacement` a second time to targets this path already blended — pass the handled-target
set through and keep emitting the same `pm_audit` dict so the postbuild drift check is unchanged; (ii) register
`VCC_NORMREL_PM` in `def_pipeline.py`'s `ARM_BOARD_DELTA` exactly as 312 registered `VCC_NORMREL_EXP`.

**Do NOT touch `apply_perturbmatch_k562_replacement`'s `k562_share`** — that function is also the incumbent's
PerturbMatch path and changing it would move `cd4svd5` and every other arm. **Freeze the two second variables**:
D5 — take the weights over 308's support (the panel-zeroed operator columns, `wcols`), and say so; D7 — keep
production's broad PerturbMatch column policy. **Do not sweep either in this arm.**

## Known-answer rows — ALL THREE must pass before you build, and they are the point of this mandate
1. `VCC_NORMREL_PM` **unset** reproduces `submissions/312-normrel15` **byte for byte** (compare the built h5ad's
   `output.sha256` against `submissions/312-normrel15.h5ad.meta.json`).
2. `VCC_NORMREL_PM=blend` with `VCC_NORMREL_EXP=0` reproduces **`cd4svd5` byte for byte**
   (`submissions/249-cd4svd5-f034.h5ad.meta.json`).
3. On 313's saved fixture, the arm's delivered row norm equals the uniform row's for **every** target to ≤ 1e-6
   relative — i.e. it is exactly shape-only, the property the shipped arm lost.
Rows 1 and 2 are cheap to check on the fixture first; only run the full build once the fixture versions pass.
**If any row fails, STOP and report — do not build, do not prep.**

## Then
One full build on the real A/B/C bundle with `VCC_NORMREL_EXP=-1.5 VCC_NORMREL_PM=blend`, following
`experiments/312-normrel-arm/run312.sh`'s gated two-invocation pattern (expect the same single postbuild
config-drift FAIL on the new key; that is the guard noticing a deliberate arm). Verify `config.normrel_exp` **and**
`config.normrel_mode`/`normrel_pm` landed in the meta. Then `vcc prep` to
`submissions/316-normrel15-pm.vcc`, copy the meta beside it, and run
`scripts/check_submission_novelty.py check`. **STOP THERE.** The coordinator preregisters the board reading and
submits; you do not.

## Deliverable
`RESULTS.md`: the diff you made (function by function); the three known-answer rows with their measured numbers;
the fixture comparison against 313's `parity313.json`; the build's postbuild line; the prep and novelty output; and
one line — **316 (parity-correct `normrel15-pm`): BUILT AND VERIFIED / BLOCKED (KA row N failed: …)**. Then stop.
