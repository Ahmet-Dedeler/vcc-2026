# 201 — BUILD `rc_hist_a` and verify every falsifier. DO NOT SUBMIT.

You are a build subagent. The coordinator has preregistered this arm at
`experiments/201-rc-hist-a/PREREGISTRATION.md` — **read it first; it is binding and must not be
edited.** Produce the artifact and try hard to break it. **NEVER run `vcc submit`.**

Work in `experiments/201-rc-hist-a/`, plus `submissions/` for the package (authorised for this
mandate only) and `src/`+`scripts/` only if a flag genuinely must be added. `uv run python`.
Peak RSS ≤ 14 GB — eight research agents are running on this 64 GB box. Builds strictly serial.

## READ FIRST
1. `experiments/201-rc-hist-a/PREREGISTRATION.md` — binding.
2. `experiments/191-rc-hybrid-build/RESULTS.md` — **the path you are reusing.** It added the emitter
   as a generated module `src/vcc2026/rc_emit.py` behind `--emitter rc_hybrid`, proved the default
   path byte-identical against the pre-change builder resolved out of git, and passed 60/60 tests.
3. `experiments/185-rank-controlled-emission/RESULTS.md` §3a (the law), §3b (the arm table), §7
   (the spec, and the `abs` vs `signed` note).
4. `experiments/185-rank-controlled-emission/BOARD-RESULT-rc_hybrid.md` — the board point.
5. `AGENTS.md` §6.1 (format contract), §11.1 **rule 10**, §14.

## The change
`rc_hist_a` differs from `rc_hybrid` in exactly one respect: the `|z|` multiset is re-assigned by
descending `|H|` across the **whole gene axis**, rather than within a15's BH call set and within its
complement separately. Everything else — a15's config, the per-gene totals, the `T = S // 400` weld,
the `abs` sign policy, the sub-gate tail law, greedy placement with sub-gate repair — is unchanged.
Add it as a second `--emitter` value; **the default path must stay byte-identical and you must prove
it, not assert it.**

## Known-answer rows (rule 10's corollary, abort on failure)
Your build must reproduce 185's h1-test `rc_hist_a` row: `reach` **0.23278**, `n_pred` **6040**,
`jac` **0.05932**, `fid` **0.50808**, `nmae` **0.96670**, `P(1)` **0.8780**, and the twin at
`reach` 0.23272. 191 showed the production emitter reproduces 185's numbers exactly through this
path, so anything else means your re-key differs from 185's.

## Falsifiers — report each PASS/FAIL with its number
1. **Pseudobulk identity against `submissions/191-rc_hybrid.h5ad`**: per-gene panel-total digests
   must match on **900/900** blocks, and the nonzero-structure digest must match on **0/900**. 191's
   `digest191.py` does exactly this — reuse it.
2. **`nmae`**: recompute the scored lfc from both artifacts against the challenge control means;
   within **0.0005 raw** of `rc_hybrid`.
3. **Challenge `n_pred`** via 191's `probe191.py` against `data/challenge/context_{A,B,C}.h5ad`
   (the measured control *is* the scorer's comparison group, so this is exact): report it against
   `rc_hybrid`'s 4,480 and a15's 4,474, and report the **call-set overlap** with `rc_hybrid` — that
   overlap is the variable this arm exists to price, so quantify it.
4. **Delivered order**: Spearman(−log p, `|H|`) across the full axis on the challenge contexts,
   ≥ +0.95 (`rc_hybrid` read +0.9743 inside its call set).
5. **Rule 10 config diff**: exactly the one pre-declared emitter difference against `rc_hybrid`'s
   meta; operator sha256 identical; `effects.whitening.after_all_median_norm` not drifted.
6. **`vcc prep --dry-run`** clean, and report nnz against the 4.75e9 cap (`rc_hybrid` sat at 34.68%).

## Deliverable
`RESULTS.md`: the known-answer row, the six falsifiers with numbers, the config diff verbatim, the
package path and sha256, and a one-line **BUILDABLE / NOT BUILDABLE** verdict. Then stop.

**Do not submit. Do not choose a model name or write a description.** The coordinator does both.
