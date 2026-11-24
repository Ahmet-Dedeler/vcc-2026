# 210 — BUILD `cd4f034` and verify every falsifier. DO NOT SUBMIT.

You are a build subagent. The arm is preregistered at
`experiments/210-cd4-build/PREREGISTRATION.md` — **read it first; it is binding and must not be
edited.** Build it, then try hard to break it. **NEVER run `vcc submit`** — the coordinator submits.

Work in `experiments/210-cd4-build/`, plus `submissions/` for the package (authorised for this
mandate only). **Do not modify `src/` or `scripts/` — this arm needs no code change**, only an
environment variable. If you believe a code change is required, stop and say so.
`uv run python`. **Peak RSS ≤ 22 GB, builds strictly SERIAL** (191 measured 19.83 GiB, 202 measured
16.95 GiB). Two research agents are running; check free memory before the emission stage and wait if
it is under 25 GB.

## READ FIRST
1. `experiments/210-cd4-build/PREREGISTRATION.md` — binding, including its six falsifiers.
2. `experiments/207-cd4-dose/RESULTS.md` §8 — the operating point, the exact invocation, and the
   pre-build meta assertions. **§8.2's table of required meta values is your checklist.**
3. `experiments/202-callset-rider/build202.sh` — **the script you are copying**, and
   `experiments/202-callset-rider/RESULTS.md` for `digest202.py` / `probe202.py`.
4. `experiments/202-callset-rider/BOARD-RESULT.md` — the incumbent's board row.
5. `AGENTS.md` §6.1 (format contract), §11.1 **rule 10**, §14.

## The change

Copy `build202.sh`. Change **one line**: `export VCC_CD4_TARGET_F=0.1` → `0.34`. Keep
`--emitter rc_callset`, `VCC_ALPHA=0.982`, and every other flag byte-identical.

**Set `f`, never `fraction`.** 207 §4 footnote 1: the builder solves `f = 0.34` to
`fraction ≈ 0.322` on the challenge panel (and 0.359 on PLOLO). `VCC_CD4_FRACTION` is the legacy
path and `effective_f_at_cd4_fraction` is a diagnostic of it, not of this arm.

**Do NOT re-solve `VCC_ALPHA`.** 207 §6.2 measured α\* = **0.982 exactly** at this dose — the
re-solve is a literal no-op, and moving α would make this two variables.

## Known-answer row (mandatory, abort on failure)
Rebuild the **incumbent** through your script first: `build210.sh rc_callset` with
`VCC_CD4_TARGET_F=0.1` must reproduce `submissions/202-rc_callset.h5ad`'s sha256
**`2a894aa543db7f574ebc024dbb3ebd5e36f923827c4e9a01f263092fe49692e0`** on **12,648,807,800 bytes**,
verified with `shasum` on the files, not read from a log. If it does not, your environment differs
from 202's and nothing after it is readable.

## Falsifiers — report each PASS/FAIL with its measured number
1. **The six pre-committed falsifiers in the preregistration §5**, as far as they are checkable
   before the board: #6 (config diff exactly 2 fields, `operator.sha256`,
   `fully_supported_targets` = 276, `alpha` = 0.982, `after_all_median_norm` unmoved) is fully
   checkable now and is the one that decides whether this is one variable. Check it against
   `submissions/202-rc_callset.h5ad.meta.json`, field by field, and print every difference rather
   than only the expected ones.
2. **Every meta value in 207 §8.2's table**, including
   `effects.cd4.residual_frame.target_f_requested` = 0.34, `fraction_source` = `VCC_CD4_TARGET_F`,
   and `solved_fraction` ≈ 0.32.
3. **The delivered profile actually moved, and by the predicted amount.** This arm is *not*
   pseudobulk-preserving — unlike the last two, the per-gene panel totals **must** change. So:
   report the per-gene panel-total digest against `202-rc_callset.h5ad` (expected to **differ** on
   900/900), and then the quantities 207 predicted — **delivered median row norm ≈ 1.0009×**,
   **delivered `A` within 1.0%**. A digest that matches on any block means the residual did not
   apply and the arm is void.
4. **Non-vacuity, explicitly.** 200 and 207 both found the CD4 residual can be *identically zero*
   while every check still passes, because the sizing rule takes a median. Assert the delivered
   effect at `f = 0.34` differs from the `f = 0.10` artifact by a stated norm, and report that norm.
   **This is the falsifier most likely to catch a silent failure.**
5. **Challenge `n_pred` unchanged** within ±1% of `rc_callset`'s **2,599**, via `probe202.py`
   against `data/challenge/context_{A,B,C}.h5ad`. The call set must not move; if it does, `fid`/
   `jac`/`reach` are in play and preregistration falsifier 5 is at risk.
6. **`vcc prep --dry-run` clean**, nnz against the 4.75e9 cap (`rc_callset` sat at 33.22%) and no
   cell above the 1e6 counts cap.

## Deliverable
`RESULTS.md`: the known-answer row, every falsifier with its number, the full config diff verbatim,
the package path and sha256, the measured non-vacuity norm, and a one-line
**BUILDABLE / NOT BUILDABLE** verdict. Then stop.

**Do not submit. Do not choose a model name or write a description.**
