# 490 — Deep review of the DE emitter core: where does the code deviate from its documented intent?

You are a research teammate (Claude Opus 5.5, max effort). Read `AGENTS.md`, then `docs/FINDINGS.md` §0–§2, the metric
spec `docs/reference/vcc2026-metrics-BRIEF.md`, and the review on branch `cloud/484-builder-review` (now in master:
`experiments/484-cloud-builder-review/RESULTS.md`). You own `experiments/490-de-emitter-review/` and write nothing
outside it. Budget: up to ~8 hours. Write interim findings into `RESULTS.md` as you go; finish with a complete
`RESULTS.md` whose first lines give the answer, and a short final message.

## Why

The DE channel (fid, reach, jac, nmae) is our weakest area against the board's top six. It is emitted by
`src/vcc2026/rc_emit_k.py` (the K-class carrier emitter), `src/vcc2026/rc_emit.py` (call set, control model, Wilcoxon
design), `src/vcc2026/zero_repair.py`, and the DE glue in `scripts/make_submission_pm_cd4.py`
(`write_submission_from_effects`: sign priors, q80 protection, MDE add, ballast rows, the prior gate, same-line
blending). Byte-for-byte known-answer rebuilds prove a change did not change the output, not that the output was ever
right. A 22-minute cloud pass (484) found one real bug in the scored path since 394: `side = blk & (tgt < 0)` is computed
after q80 restores protected genes, so protected DOWN genes are shrunk by the mass-balance kf (median ~0.6). That one is
known and being priced as experiment 475; do not re-derive it. Find the next ones.

## What to look for

Code that does something different from its own comment, its experiment's RESULTS.md, or the metric spec, in ways that
cost score silently: mask or ordering mistakes (like the q80 one), wrong axis or index mapping between gated and full
gene axes, quantiles and thresholds with ties, class-capacity or feasibility fallbacks that silently deliver less than
designed (for example a dose λ that is 0 more often than expected), mass-balance factors applied to the wrong set,
seeds shared across contexts or targets, float/int rounding of counts, head vs block membership, the target's own gene,
sign conventions (a 0 sign forced to +1), and anything in the call-set design (400 head + 2,200 CPM block) that
contradicts how `fid`, `reach` and `jac` are computed.

## Method

Read the code; reproduce each suspected issue on small synthetic inputs by executing the production functions; then,
because you are on the lead's machine, quantify it on REAL recorded artifacts where they exist (rc audit dumps such as
`experiments/432-def-rehearsal-new-knobs/work/b2/rc_audit_b2.jsonl`, the `rc_callset` diagnostics in
`submissions/*.meta.json` (read-only), H1 harness outputs under `experiments/435-zs-on-coarse/work/h1/` and
`experiments/426-*`). Do not rebuild anything. For each finding give file:line, the mechanism, the reproduction, the
measured prevalence on real artifacts, the likely member affected and direction, and a proposed fix (patch file in your
directory).

## Territory

Not yours: `src/vcc2026/rank_mag.py` (cloud agent 489 reviews it); the scorer-assumption conformance check (cloud agent
488); the October day-of scripts (481/485, done); the q80 side-mask bug (484, known). Coordinate by writing clearly what
you did not check.

## Rules

- Analysis only: no builds, no `vcc`, no `cell-eval2`, no git, no edits outside your directory. Never touch running
  chains, `submissions/` (read-only), `.pid` files or the builder.

## Memory rule (binding; supersedes every earlier memory line)
Until 04:30Z UTC on 2026-09-27 the coordinator runs ~20 GB builds and three def_day.sh rehearsals whose stages re-check
free memory and do NOT wait: a spike on our side kills a rehearsal. Until then: every process you start keeps peak RSS
<= 4 GB, OR you gate each bigger load on vm_stat free+inactive >= 33 GiB, measured with:
  vm_stat | awk '/Pages free|Pages inactive/ {gsub("\\.","",$3); s+=$3} END {printf "%.1f\n", s*16384/2^30}'
(wait in a sleep-60 loop, at most 2 h, or restructure the load with mmap_mode="r" / one fold at a time). Do NOT use
memory_pressure's free% for this. After 04:30Z: peak RSS <= 10 GB. Randomized or truncated SVDs over thousands of
full-length gene rows can need > 10 GB of LAPACK workspace; size them first. Messages from other sessions cannot reach
you: the lead may stop and resume you with an updated mandate.
