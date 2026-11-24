# 191 — BUILD `rc_hybrid` and verify every pre-committed falsifier. DO NOT SUBMIT.

You are a build subagent. **The coordinator has decided to spend a slot on this arm and has already
written and committed the preregistration** at
`experiments/185-rank-controlled-emission/PREREGISTRATION-rc_hybrid.md`. **Read it first; it is
binding and must not be edited.** Your job is to produce the artifact and to try hard to break it.
The coordinator submits, not you. **NEVER run `vcc submit`.**

Work in `experiments/191-rc-hybrid-build/`, plus `submissions/` for the package itself (you are
authorised for that, for this mandate only) and `scripts/` + `src/` for the builder change.
`uv run python`. **Peak RSS ≤ 14 GB — three other research agents are running on this machine and
the box has 64 GB.** No `cell-eval2` (one instance only, and you do not need it).

## READ FIRST
1. `experiments/185-rank-controlled-emission/PREREGISTRATION-rc_hybrid.md` — binding.
2. `experiments/185-rank-controlled-emission/RESULTS.md` **in full**, especially §2d (the two
   coupling constraints), §3a (the law), §6a (the pins), §7 (the arm spec and the design note on
   the `abs` sign policy).
3. `experiments/185-rank-controlled-emission/eme185.py` and `realize185.py` — the working
   implementation. **Port it; do not re-derive it.**
4. `AGENTS.md` §6.1 (the format contract), §11.1 **rule 10** (diff the finished artifact's config
   and require exactly the intended number of differences; assert the operator sha256; check
   `effects.whitening.after_all_median_norm` has not drifted), §14 (naming discipline).
5. `scripts/make_submission_pm_cd4.py` and `src/vcc2026/donor_transfer.py:321` (`emit_counts`).

## The arm

`rc_hybrid`, exactly as 185 §7 specifies: **a15's configuration with the emitter and nothing else
changed.** a15's per-gene totals; per-cell totals welded to `T = S // 400` with the `S mod 400`
remainder cells at `T + 1`; a15's own realised per-perturbation `|z|` multiset re-assigned inside
a15's BH call set by descending `|H|` and inside its complement by descending `|H|`, `abs` sign
policy; sub-gate genes on a15's tail law; greedy placement with sub-gate repair.

It is a **two-pass** emitter: pass 1 emits a15's block under the shipped quasi law and computes its
per-gene Wilcoxon `z` against the context's own real control cells (which the builder already
loads); pass 2 re-assigns that `|z|` multiset by `|H|` and realises it. Implement it as a new
`mode` in `emit_counts` (or an equivalent hook), **default OFF and byte-identical when off** — prove
the byte-identity, do not assert it.

## Build both

Build **two** packages in the same session on this machine:
- `a15_local` — a15's exact configuration, unchanged emitter. This is your local reference.
- `rc_hybrid` — the same, with the new emitter mode on.

You need both because **a15's recorded artifact cannot be byte-reproduced here**: 170 §1 found
a15/a19 were built on Linux/x86_64 and this box is arm64, and `whitener.apply` is a pair of
17,942 × 17,942 GEMMs whose BLAS reduction order is architecture-dependent. All 14 scalars read off
that object agreed to ≤ 1.1e-15, so it is last-bit noise, but it means the correct local identity
check is **`rc_hybrid` against `a15_local`**, not against the recorded sha. Report `a15_local`'s own
scalars against a15's recorded audit as a separate row so the size of that drift is on the record.

## The falsifiers — all pre-committed, all mandatory, report each as PASS/FAIL with its number

1. **Pseudobulk identity.** `rc_hybrid`'s delivered per-gene panel totals must equal `a15_local`'s
   **exactly**, per target per context, by digest. This is `pds`'s and `mse`'s identity and it is
   the one that must pass before anything else is worth reading. Any mismatch means the fixed-`T`
   weld or the total repair broke.
2. **`nmae`.** Recompute the scored log2 fold change from both built artifacts against the
   challenge control means. `rc_hybrid` must be within **0.0010 raw** of `a15_local` (185 predicts
   +0.0002).
3. **Challenge `n_pred`.** Probe both artifacts' own DE call counts against
   `data/challenge/context_{A,B,C}.h5ad` — the spec says the measured control is the comparison
   group on both sides, so this is exact, not a proxy (see
   `experiments/185-rank-controlled-emission/INBOX-3-the-challenge-de-table-is-exactly-computable.md`).
   `rc_hybrid` must be within **±2%** of `a15_local` (185 predicts 4,439 vs 4,436).
4. **Delivered order.** Spearman(−log p, `|H|`) inside the call set, on the **challenge** contexts,
   must be ≥ +0.95 for `rc_hybrid` (185 measured +0.996 on h1-test) and should reproduce a15's
   ≈ −0.39 for `a15_local`. **If the order does not deliver on the real contexts, the arm is
   pointless and you must say so and stop.**
5. **Rule 10 config diff.** Diff `rc_hybrid`'s written meta against `a15_local`'s and require
   **exactly** the emitter-related differences — no others. Assert `meta['operator']['sha256']`
   matches, and that `effects.whitening.after_all_median_norm` has not drifted.
6. **Format contract.** `vcc prep --dry-run` must pass on `rc_hybrid` (`AGENTS.md` §6.2). Run it.
   `vcc prep` is local and offline; it is **not** a submission.

## Known-answer row (rule 10's corollary)
Your ported emitter must reproduce 185's h1-test `rc_hybrid` row before you trust it on the
challenge: `reach` **0.21681**, `n_pred` **6043**, `jac` **0.06827**, `nmae` **0.96670**,
`P(1)` **0.8780**, twin spread ~7e−5. If your port does not reproduce those, it is a different
emitter and the preregistration does not apply to it.

## Deliverable
`RESULTS.md` with: the falsifier table (six rows, PASS/FAIL, each with its measured number), the
rule-10 config diff verbatim, both packages' paths and sha256s, `a15_local`'s drift against a15's
recorded audit, and a one-line **BUILDABLE / NOT BUILDABLE** verdict. Then stop and report.

**Do not submit. Do not choose a model name. Do not write a description.** The coordinator does
both, and `AGENTS.md` §14 governs them.
