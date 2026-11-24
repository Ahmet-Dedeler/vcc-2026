# 200 — Leave-one-component-out on the SHIPPED pipeline. Every part was validated when it was added; none has been re-validated since.

You are a research subagent. Work ONLY in `experiments/200-production-ablation/`.
No `vcc`, no `cell-eval2`, no submission build, no writes outside your directory.
Peak RSS ≤ 12 GB. `uv run python`. Long compute is fine.

## READ FIRST
`AGENTS.md` §4, §4.2, §8.1, §11.1 rules 1–13 with the audit-181 amendments. Then:
- `scripts/make_submission_pm_cd4.py` — **the shipped builder. Read its flag documentation; it is
  extensive and it records why each default is what it is.**
- `experiments/162-retrodicting-validator/` (PLOLO — the only `pds` instrument validated against a
  board point, 6/6 signs, LOO RMS 0.001494 scaled).
- `experiments/188-magnitude-calibration/RESULTS.md` Phase 3 — **`debench.Surface.pds()` must never
  size a `pds` arm.** Standing rule. And its DE-pool rule: build the pool by **significance**
  (`|lfc|·√cpm` at the reference's BH depth, median ~89 genes/row), never by magnitude, or `nmae`
  deltas invert in sign.
- `docs/BOARD-LEDGER.md` — every scored arm, all six members.

## The gap

a15 is a stack of components, each added by its own experiment and each validated **at the time it
was added, against the then-current baseline, on the then-current instruments**. In order:
the union whitening operator and `WHITEN_GAMMA = 0.65` (017/019/020, gamma re-confirmed as a
six-member local optimum by 170), the pruned basis (091), PerturbMatch K562 replacement (058/061),
the post-whitening CD4 residual (062, and a15's own +0.007256 scaled `pds`), the gene
normalisation and CPM gate (006), equal donor weights (019), `VCC_ALPHA = 0.982` (a14/a11, and 187
confirms it is the joint `mse`+`nmae` optimum), the tiered spread policy
(`--keep-top 8000 --keep-spread 0.25 --tail-spread 1.0`), `--overdispersion 4.0` (a8), and a23's
mass fix.

**Nobody has ever removed one and re-measured.** Three reasons that matters now:

1. **The instruments have changed underneath them.** Several were validated on surfaces we have
   since disqualified — the Orion harnesses carry a target-independent blob (26.4% / 44.1% of truth
   energy), `debench.Surface.pds()` disagrees in *sign* with PLOLO on the 5–10 CPM band, and DE
   pools built by magnitude invert `nmae`. A component validated on one of those has an unvalidated
   verdict.
2. **The objective has changed.** 187 closed calibration, 184 closed line adaptation, 189/190 closed
   the ordering key. A component whose value was measured against `pds` alone may be negative on the
   six-member sum, and vice versa — a13 was gated positive on both harnesses and lost **−0.0203**
   on the board.
3. **Interactions.** Components were added sequentially against a moving baseline, so each one's
   measured value is conditional on the stack below it and on nothing above it. a23's mass fix, for
   instance, was worth **+0.004014 scaled `pds`** on PLOLO and cost **−0.0079 scaled `reach`** on
   the board — a trade nobody re-examined after `reach` became the identified deficit.

## Work

### Phase 0 — known-answer rows (mandatory, abort on failure)
- Reproduce a15's PLOLO row exactly, and 169's `q = 1.00` identity (16/16 bit-identical, worst
  |Δpds| **0.00e+00**), and the a15 → a23 mass-fix price **+0.004014 scaled**.
- Reproduce `claim(a15 as shipped)` = **0.24581** and a15's raw `mse` / `score_mse` =
  **0.88228 / 0.11112**.
- **Rule 11: record a path and a sha256 for every operator and frozen input you load.** Eight
  whitening operators exist on disk with near-identical names and max eigenvalues 8.036 / 6.257 /
  5.895, and an evening was once lost to loading the wrong one.

### Phase 1 — leave-one-out
Turn each component **off, one at a time**, from the full a15 stack, and price all six. Also run
the reverse where it is cheap (add-one-in from a bare consensus) — the two together separate a
component's marginal value from its interactions, which sequential addition cannot.

Instruments, per member: `pds` on **PLOLO** (162's four folds, challenge frame, pre-committed
unweighted mean); `mse` through **096's `claim` route** on the real challenge panel; `nmae`,
`reach`, `fid`, `jac` on **h1-test**, emitted, ≥3 matched seeds, with a **significance-selected
pool**. Donors on h1-test: `k562_gwps` + `orion_hek293t` only — `data/donors/h1_2025` covers 85/85
of h1-test's targets and is pure oracle leakage.

Report each component's six-member sum in **overall board units** using the scoring map
(`dOverall/dRaw`: `pds` +0.3690, `mse` −0.1752, `nmae` −0.2755, `fid` +0.5607, `reach` +0.1882,
`jac` +0.4544), against the composed gate (+0.0019 shape-only, +0.0089 amplitude-changing).

### Phase 2 — the question that pays
**Is any shipped component now NET NEGATIVE?** If one is, removing it is free score and it is the
cheapest arm available. Rank every component by its measured six-member value and say plainly which
are load-bearing, which are noise, and which are costing us.

Pay particular attention to the ones whose trade has shifted:
- **a23's mass fix** (+0.004014 `pds` on PLOLO, **−0.0079 scaled `reach`** on the board) — with
  `reach` now the identified deficit, does it still clear?
- **the tiered spread policy** — 171 §12c measured that it parks the signal block at the
  **minimum-|z|** point while the noise tail leads the reach walk (61–79% of head positions at sign
  accuracy 0.5212 against the keep block's 0.6629), and our board `reach` peaked at **a8**, which
  used a flat 0.25 spread and no tier (`docs/BOARD-LEDGER.md` §3). It buys `mse` through `claim`;
  price both halves.

### Phase 3 — falsifiers
1. **Every ablation must leave the untouched members where they belong.** An emission-only ablation
   must leave the per-gene panel totals digest-identical (191 proved the check works, 900/900
   blocks); a profile-touching one must not.
2. **Rule 10 config diff** on every arm: exactly the intended number of differences in the written
   meta, operator sha matching, `effects.whitening.after_all_median_norm` not drifted.
3. **Retrodict where a board point exists.** a23's mass fix, a19's ridge (−0.01352 on `mse`), a8's
   `--overdispersion 4.0` (+0.021398 scaled `reach`) and a15's CD4 residual (+0.007256 scaled `pds`)
   are all board-measured. Your ablation of each should recover its sign and rough size, and audit
   181's rule 1 says you must show that before your other numbers are read.

## Deliverable
`RESULTS.md`: known-answer rows with shas, the full leave-one-out and add-one-in table in overall
units, the retrodictions, and a ranked verdict naming any component that should be dropped. **"The
stack is coherent and every component earns its place" is a full result** and would be the first
time anyone has been able to say it.

## Territory
192 call set · 193 operator eigenvalue map · 194 donor tables · 195 learned transform · 196 emitted
depth · 197 cosine inversion · 198 closure re-audit · 199 D/E/F readiness. **You own the shipped
stack as a whole.** Where your ablation touches another's territory, price it and hand the
follow-up over rather than tuning it.
