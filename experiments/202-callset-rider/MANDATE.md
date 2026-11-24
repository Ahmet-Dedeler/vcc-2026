# 202 — BUILD `rc_callset` (192's rider) and verify every falsifier. DO NOT SUBMIT.

You are a build subagent. The arm is preregistered at
`experiments/202-callset-rider/PREREGISTRATION.md` — **read it first; it is binding and must not be
edited.** Build it, then try hard to break it. **NEVER run `vcc submit`** — the coordinator submits.

Work in `experiments/202-callset-rider/`, plus `submissions/` for the package (authorised for this
mandate only) and `src/`+`scripts/` for the flag. `uv run python`.
**Peak RSS ≤ 20 GB and builds strictly SERIAL — six research agents are running on this 64 GB box
and 191 measured a build peak of 19.83 GiB. Check free memory before starting the emission stage
and wait if it is below 25 GB.**

## READ FIRST
1. `experiments/202-callset-rider/PREREGISTRATION.md` — binding.
2. `experiments/192-predicted-callset/RESULTS.md` §5 — the arm, its measured deltas, and the four
   falsifiers its author pre-committed. **This is the science; you are the engineering.**
3. `experiments/191-rc-hybrid-build/RESULTS.md` — **the path you are extending.** The emitter is a
   generated module `src/vcc2026/rc_emit.py` behind `--emitter rc_hybrid`; the default path was
   proven byte-identical against the pre-change builder resolved out of git; 60/60 tests pass; and
   `digest191.py` / `probe191.py` are the verification tools you should reuse rather than rewrite.
4. `experiments/185-rank-controlled-emission/BOARD-RESULT-rc_hybrid.md` — the board point this arm
   is marginal to, including the fact that `pds` and `mse` reproduced to **sixteen significant
   figures**, which is the precision your identity claims are held to.
5. `AGENTS.md` §6.1 (format contract), §11.1 **rule 10**, §14.

## The change

Exactly one variable against `rc_hybrid`: the **call set**.

```
call set per (context, perturbation) = top 400 by |H|  ∪  top 2,200 by control CPM   (n_pred 2,600)
order                                = rc_hybrid's law, unchanged
everything else                      = rc_hybrid, unchanged
```

`control CPM` is the per-gene mean CPM of that context's own NTC control cells from
`data/challenge/context_{A,B,C}.h5ad` — measured, not modelled, and the same object the scorer gates
on. Ties broken deterministically and documented. Add it as a new `--emitter` value or an explicit
flag on the existing one; **the `rc_hybrid` and default paths must stay byte-identical, and you must
prove it against the current builder resolved out of git, not assert it.**

## Known-answer rows (rule 10's corollary, abort on failure)
- Rebuild `rc_hybrid` through your modified builder and require it to reproduce
  `submissions/191-rc_hybrid.h5ad`'s sha256 **exactly** (`e0b0677f984a94a1c4a7af2add196908cdab90d7acb463eb1d6e0397c2d870a8`,
  13,202,460,120 bytes). If your change perturbs the incumbent at all, stop.
- Reproduce 192's h1-test row for this arm: `fid` +0.0010, `jac` +0.0129, `reach` +0.0112 raw
  marginal to `rc_hybrid`. If your call set differs from 192's, these will not reproduce — that is
  the check.

## Falsifiers — report each PASS/FAIL with its measured number
1. **Per-gene panel-total digests identical to `191-rc_hybrid.h5ad` on 900/900 blocks**, and the
   nonzero-structure digest matching on **0/900** (so the files are genuinely different objects with
   the same pseudobulk). Reuse `digest191.py`.
2. **Challenge `n_pred` within ±3% of 2,600**, via `probe191.py` against the real control cells —
   this is exact, not a proxy, because the measured control *is* the scorer's comparison group.
   Report per context and pooled, against `rc_hybrid`'s 4,480.
3. **Call-set composition**: report the overlap with `rc_hybrid`'s call set, the size of each half of
   the union, and how many genes the two halves share. The arm's whole logic is that the halves do
   different jobs; show they do.
4. **Delivered order preserved**: Spearman(−log p, `|H|`) inside the call set on the challenge
   contexts, against `rc_hybrid`'s **+0.9743**. A drop here would mean the set change disturbed the
   order, which is falsifier 5 in the preregistration.
5. **Rule 10 config diff**: exactly the pre-declared call-set difference against `rc_hybrid`'s
   written meta; `meta['operator']['sha256']` identical; `effects.whitening.after_all_median_norm`
   not drifted.
6. **`vcc prep --dry-run` clean**, and nnz against the 4.75e9 cap (`rc_hybrid` sat at 34.68%; a
   smaller call set does not change the emitted counts, but report it).

## Deliverable
`RESULTS.md`: the two known-answer rows, the six falsifiers with numbers, the config diff verbatim,
the package path and sha256, and a one-line **BUILDABLE / NOT BUILDABLE** verdict. Then stop.

**Do not submit. Do not choose a model name or write a description.**
