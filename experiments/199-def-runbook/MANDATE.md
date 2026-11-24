# 199 — Make Oct 22 mechanical. Pre-stage and rehearse the whole D/E/F transition, end to end.

You are a build/readiness subagent. Work in `experiments/199-def-runbook/`; you may also add
**new** files under `scripts/` (never modify an existing script's behaviour without proving
byte-identity). No `vcc submit` ever. No `cell-eval2`. Peak RSS ≤ 10 GB. `uv run python`.

## Why this exists and why now

**Oct 22: the final panel (contexts D/E/F, a NEW 300 targets) is released and the leaderboard goes
dark. Nov 5 23:59 UTC: hard deadline. Only the LAST entry counts, and every submission in that
window is blind.** Nothing below can be done under time pressure on the day, and everything below is
currently either undone or exists only as prose in `docs/CURRENT_STATE.md`'s transition plan.

Items 1 and 2 are done (experiment 180: the builder reads contexts from the manifest, the
frozen-input table lives in `scripts/frozen_inputs.json` with `--repin-frozen-inputs`, and a
synthetic D/E/F bundle was rehearsed through `vcc prep --dry-run` producing byte-identical
`X/data`, `X/indices`, `X/indptr` to a23's). **Items 2b through 7 are open. You own them.**

## READ FIRST
`AGENTS.md` §2, §2.1 (how the panel was built), §3 (timeline), §6 (submission mechanics and the
slot banner), §11.1 rules 10 and 11. Then `docs/CURRENT_STATE.md`'s **"D/E/F TRANSITION PLAN"**
section in full, `experiments/180-def-readiness/` (what is already done, and its §10 runbook),
`experiments/041-ccle-match/` (the matcher), `experiments/166-donor-acquisition/` (`census.py`).

## The items

### 2b — the CCLE matcher, pre-staged to run in minutes
041 identified the validation contexts from control profiles alone against the DepMap 25Q3 CCLE
mirror — **A = JURKAT, B = HELA, C = CAL33** — the same top match under all eight feature/rank
variants, calibrated first on three known controls (K562 and HCT116 at rank 1). That identification
is what motivated the **CD4 residual now shipping in a15**, so line identity selects auxiliary
donors and has already produced a shipped component. §2.1 adds that all six lines are **cancer**
lines, which narrows the search space.

Deliverable: a single script that takes the D/E/F control `.h5ad` files and emits ranked identities
with the calibration rows re-run, in **minutes, not hours**. Verify it by reproducing 041's A/B/C
call exactly. Pin the CCLE mirror locally now — do not depend on a download working on Oct 22.

### 3 — donor coverage census on the new panel
K562-gwps covers **269/300** of the current panel and Orion 290–294. §2.1 says targets are selected
for *"adequate representation in existing public datasets"*, so coverage should hold — **verify,
don't assume.** Wire 166's `census.py` to run on a new `pert_counts.csv` in one command, and report
what the consensus looks like at 3, 2 and 1 covering donors so a coverage drop has a pre-decided
response rather than an improvised one.

### 4 — gene axis and manifest verification
The writer asserts the gene axis and `cells_per_pert`. Make a one-command preflight that checks the
new `gene_names.csv` order, `pert_counts.csv`, the manifest's contexts and counts, and the control
files' shapes and dtypes, against the contract in §6.1. **Remember both CSVs have a header row even
though `vcc prep --genes` documents "headerless"**, and that gene order must be read from the file,
never assumed.

### 5/6 — the blind decision rule, written down before it is needed
`n_real` **cannot be measured on D/E/F** and the leaderboard is dark. So every final-round choice
must be pre-decided from A/B/C evidence. Assemble, from the experiment record, a single decision
table: for each shippable flag, its A/B/C-validated setting, the instrument that validated it, and
the rule for what to do if the D/E/F controls look different. Specifically carry:
- **`VCC_ALPHA`: global 0.982, NOT the per-context knee** — abc-knee's board Δ`mse` was −0.00226
  against a local +0.005. 187 adds the expected-value rule: under 1% cliff uncertainty the optimum
  is **knee × (1 − 1.25%)**, worth +0.0009 over a global 0.982, and the raw knee is expected to gain
  ≈ 0. Record both and mark which ships.
- the emitter, the call-set policy, `WHITEN_GAMMA` 0.65 (170: a six-member local optimum), the
  operator **by path and sha256** (rule 11: eight operators exist on disk with near-identical names
  and max eigenvalues 8.036 / 6.257 / 5.895).

### 7 — the submission discipline, as a checklist
Submit the best-supported package **early** in the window, then replace it only for something
strictly better on the local instruments with a preregistered prediction. **`vcc prep` hardcodes
`REQUIRED_CONTEXTS = ("A","B","C")` — pass `--contexts D,E,F`** (180 §10). Include the upload-time
budget: a 3.5 GiB package took **15m31s** at ~3.8 MiB/s, so a same-day rebuild-and-submit is ~1h of
wall clock before scoring even starts. Put that in the checklist; it is the kind of thing that loses
a round.

## Rehearsal — the deliverable that matters
Build a synthetic D/E/F bundle (A/B/C renamed, new manifest, a permuted 300-target list drawn from
the donors' coverage so it is *not* the current panel) and run the **entire** path end to end:
preflight → CCLE matcher → census → build → `vcc prep --dry-run --contexts D,E,F`. Time each stage.
Report the total wall clock, because that number decides how the Oct 22 day is planned.

**Rule 10 is binding: the A/B/C build must stay bit-identical under every change you make.** Prove
it against the current builder resolved out of git, and report the diff of the written meta.

## Deliverable
`RESULTS.md` plus `docs/DEF-RUNBOOK.md` — a numbered, executable runbook a tired person can follow
on Oct 22 without rediscovering anything, with the decision table, the timings, and every path and
sha256 pinned. **No production behaviour changes; no submission.**

## Territory
192–198 own the science. You own **readiness**. Do not tune anything.
