# 229 — Generalise the source registry and BUILD 225's four-source stack for the 2026-09-05 slot

You are a build subagent. Work in `experiments/229-stack4-build/`. **You MAY modify
`scripts/make_submission_pm_cd4.py` and you MAY build a submission package.** You may **NOT** run
`vcc submit`, choose a model name, or write a description — the coordinator spends the slot.
Peak RSS ≤ 22 GB; `vcc prep` alone peaks at **22.33 GiB** and must **never** run alongside a build.
Five other agents are live (221, 226, 227, 228, and 225 finishing) — **check before any heavy run**,
and never run `cell-eval2` concurrently. `uv run python`.

## READ FIRST, in this order
1. **`docs/NEXT-BUILD-SPEC.md`** — the whole job is specified there. Read it before anything.
2. **`experiments/225-source-sweep/RESULTS.md`** — the arm you are building and its evidence.
3. **`experiments/222-h1-union-build/RESULTS.md`** — the pattern to copy. 222 did exactly this job
   for the H1 union: default-off flag, byte-identical default path, three known-answer rows, five
   falsifiers. Its arm scored **+0.0069114, rank 3 of 656**. **Copy its discipline exactly.**
4. `AGENTS.md` §8.4, §11.1 (rule 13 and the digest test), §14.

## The arm

225 swept the space (83 assets + GEO + BioStudies + CZI + Zenodo + HuggingFace + IGVF + Arc;
**nothing panel-covering exists that we do not hold**) and priced 15 candidates through a harness
that reproduces 223's five headline rows to the last digit and 218's baselines at **Δ = 0.0e+00**.
Its recommendation:

> **`dld1_lowmoi ∪ jiang_5line_noK562 ∪ igvf_cardiomyocyte ∪ myllia_a549` at `g ≈ 0.35`** —
> **2× the shape-only gate** on a leave-one-fold-out dose, **3/3 folds positive**, F1 and F3 passed
> at every dose, against DLD-1 alone at 1.07×. Cumulative PLOLO price **+0.00421** over the shipped
> `h1union` arm.

**Take 225's tables and its dose. Do not re-derive them, and do not re-tune `g` against your own
build** — that would be fitting the dose to the instrument that scores it. If you believe `g` is
wrong, say so in `RESULTS.md` and build the pre-registered one anyway.

## The code change — this is the only production edit authorised

`VCC_H1_UNION_G`'s machinery is already source-agnostic (standard compact-donor layout, carrying-row
sizing, merge with the CD4 residual into one add, single global preserve rescale). **The only
H1-specific thing is the hardcoded `H1_UNION_SOURCES` tuple** (`scripts/make_submission_pm_cd4.py`
≈ line 636). Make the source list **registry-driven**: a name→path table with **sha256 per entry**
(rule 11: a path and a sha, never a description), plus an env var naming which entries are active.

**Three known-answer rows govern this change and all three must pass:**
1. **Flag unset** → the `.h5ad` must hash to `7036e89f8a03c82d69e79a5c6690b275adae52dbb3ff72bb6ab23e7ca622bf7a` (`cd4f034`).
2. **H1 pair selected at `g = 0.28`** → the `.vcc` must hash to
   `f13ce40bca908af7cb85b8e84887af473aae4588410c238529139b7d1f245b1e` (2,824,785,920 bytes), the
   shipped `h1union`.
3. **`solve_cd4_fraction` untouched** — CD4 sits at 93–98% coverage, its own solve is deliberately
   outside this flag (216/217's territory). 222 proved the fraction did not move by one bit;
   **prove it again** (`effects.cd4.fraction` **0.32207201067578145**, `target_f_requested` 0.34,
   `f_per_unit_fraction` 1.0556645369046553, `fully_supported_targets` 276).

**A build that cannot reproduce rows 1 and 2 does not ship.** Diff the config against `h1union`'s
and show that every difference is the pre-declared source-list change plus its deterministic
consequences — that is 222's falsifier 1 and it is the one that matters most.

## Falsifiers on the new package (222's five, adapted)
1. **Config diff** is exactly the source-list change and its consequences; `operator.sha256`
   unchanged, γ 0.65, λ 0.0, `VCC_ALPHA` 0.982, `VCC_CD4_TARGET_F` 0.34, emitter `rc_callset`.
2. **Digest test** — the per-gene panel-total digest, rule 13's single-valued test. Report it. A
   profile-CHANGING emission law costs **−0.0146 scaled**; profile-class costs 0.
3. **Non-vacuity** — every delivered block must move; report the minimum relative displacement.
4. **Carrying/non-carrying separation** — carrying rows displaced materially more than
   non-carrying (222 measured 3.3×).
5. **`vcc prep --dry-run` clean**, and pre-prep the `.vcc` so the coordinator can submit instantly.

## Also report, because the coordinator needs it to size the arm
- The **six-member** PLOLO price at the built dose, three cancer folds, **`kolf2` separately** — it
  is pluripotent, gains 3–5×, and inflated a prediction on 2026-09-04 badly enough to force a band
  reissue mid-flight.
- **Leave-top-3-out retention on cosine and `pds` separately.** On `h1union` the board realised
  ~**100%** of the cosine prediction and only **70%** of the `pds` one.
- **Do not convert PLOLO to board points with a scalar** — the factor is **not identified** (1.384
  on `cd4f034`, 1.000 on `h1union`, matched `cos_agg`). Quote PLOLO units; name the fold set,
  functional, gate and tag.
- **October**: each source's projected coverage under `DEF-RUNBOOK` §5's two draws. The panel is
  **replaced 2026-10-22 with 300 disjoint targets** and only the last entry after that counts.

## Stop conditions
- **If the source-list change alters CD4's solved fraction by even one bit, STOP and report.** A
  two-variable arm is worth less than an unspent slot.
- **If known-answer row 1 or 2 fails, STOP.** Do not "explain" a hash mismatch.

## Deliverable
`RESULTS.md` with a one-line **BUILDABLE / NOT BUILDABLE** verdict, the package path and sha256, the
three known-answer rows, the five falsifiers, the six-member price, retention, and the October
table. **Do not submit. Do not choose a model name or write a description.**
