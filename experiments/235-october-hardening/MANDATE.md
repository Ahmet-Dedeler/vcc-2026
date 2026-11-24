# 235 — Make the Oct-22 silent failure LOUD, and prove the October sources are sizeable. 233's items 2b and 2c.

You are a build subagent. Work in `experiments/235-october-hardening/`. **You MAY modify
`scripts/make_submission_pm_cd4.py` and `scripts/build_union_whitener.py`.** You may **NOT** run
`vcc submit`, choose a model name, or write a description. `uv run python`.

## MEMORY RULE — corrected, because my earlier one was wrong
233 showed my previous instruction ("hold if `vm.swapusage free < 2 GB`") is **not a pressure
reading**: macOS swap is a *dynamically sized file*, its free tail sat at 0.8–1.7 GB all session
while `memory_pressure` reported 86% free, and obeying it literally would stall indefinitely.
**Use 233's gate instead:** proceed when
`memory_pressure ≥ 40% free` **AND** `vm_stat` free+inactive **≥ 20 GiB** **AND** swap free
**≥ 400 MiB**. Run heavy jobs **strictly serially**. Report peak RSS with `/usr/bin/time -l`.
A job that exceeds memory here is killed **without a traceback** (219, 226).

## Why this experiment exists
233 measured the Oct-22 transition and named three items. **One is already done** — 229 registered
the four October-only sources with paths and shas (`0204168`). **Two are not, and both are cheap and
pay only in October.**

### Item A — the operator's silent failure is still silent 🔴
`DEF-RUNBOOK` §0.1: the whitening operator is **a function of the panel** and fails **silently** on a
new one — **98.8% of coordinates paired with the wrong gene while every check passes.** 233 priced
skipping the refit at **−0.021 overall on `pds` alone**, independently confirming 212's **−0.0429**.

**The only check that exists is `n_cols`, and it is 17,942 on every panel** (`make_submission_pm_cd4.py:1248`,
`build_union_whitener.py:363`). It cannot detect the failure it is there to detect.

**Fix:** write an **`axis_sha256`** — a hash of the operator's ordered column/gene axis — into the
operator's metadata at fit time, carry it in the build's `inputs` block, and **assert it matches the
panel the build is running against, raising loudly on mismatch.** Follow rule 11: a path and a sha,
never a description.

### Item B — prove the carrying-row rule covers the October sources
233: the four October-only sources price at **+5.5% to +17.7% relative cosine each** on a disjoint
panel (dose-monotone, clean shuffle control) — **but the all-rows median is `0.0000` for all four on
both panels**, so under the old sizing rule they deliver **identically zero**. 222 put the
carrying-row rule in for the H1 union. **Prove it covers any source selected through the registry**,
including these four, and add a test that fails if it ever does not.

## Known-answer rows — all must pass, and they govern the change
1. **Flag unset** → `.h5ad` hashes to `7036e89f8a03c82d69e79a5c6690b275adae52dbb3ff72bb6ab23e7ca622bf7a`.
2. **H1 pair at `g = 0.28`** → `.vcc` hashes to
   `f13ce40bca908af7cb85b8e84887af473aae4588410c238529139b7d1f245b1e` (2,824,785,920 bytes).
3. **Re-fitting the operator for the A/B/C bundle reproduces `e2f65ecf…`** — and its **new
   `axis_sha256` must be stable across two independent refits.**
4. **`solve_cd4_fraction` untouched**; `effects.cd4.fraction` **0.32207201067578145** unmoved.
5. **The assertion must FIRE.** Deliberately run the A/B/C operator against a synthetic disjoint
   panel (`DEF-RUNBOOK` §5's bundles) and show it **raises**. **An assertion never observed to fire
   is not a guard** — this row is the point of the experiment.

## Scope discipline
- **Do not change any shipped number.** This is a guard plus a test. If a byte of any existing
  build changes, you have done it wrong.
- **Do not price arms, do not tune doses, do not touch `solve_cd4_fraction`.**
- **Commit both files when your known-answer rows pass**, staging those paths explicitly — never
  `git add -A`, five agents share this repo — with the hashes in the message (`AGENTS.md`'s
  build-agent rule).

## Deliverable
`RESULTS.md`: the five known-answer rows (row 5 quoted verbatim — the raise), the diff summary, the
`axis_sha256` for the A/B/C operator, the carrying-row proof for all four October sources, and a
one-line **HARDENED / NOT HARDENED** verdict. Then update `docs/DEF-RUNBOOK.md` §2 so the Oct-22
driver checks the axis hash at the step where it re-fits.
