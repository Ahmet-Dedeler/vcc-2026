# 192 — Predict the reference's OWN DE set and call exactly it. `jac` and `fid`, where we rank 169th and 176th of 612.

You are a research subagent. Work ONLY in `experiments/192-predicted-callset/`.
No `vcc`, no `cell-eval2`, no submission build, no writes outside your directory.
Peak RSS ≤ 12 GB. `uv run python`.

## READ FIRST
`AGENTS.md` §4, §8.1 (the two-channel split), §11.1 rules 1–13 with the audit-181 amendments.
`docs/reference/vcc2026-metrics-BRIEF.md` §3–§6 — **read the actual definitions of `fid` and `jac`,
do not work from a summary.** Then:
- `experiments/188-magnitude-calibration/RESULTS.md` — the discovery this mandate is built on.
- `experiments/185-rank-controlled-emission/RESULTS.md` §4 (Part C) and its
  `INBOX-1-per-perturbation-sizing.md` and `INBOX-5-the-fid-jac-payoff.md`.
- `experiments/191-rc-hybrid-build/RESULTS.md` — the production emitter that can realise any call set.
- `experiments/092-call-set-truncation/`, `093-spread-truncation/BOARD-RESULT.md`,
  `164-callset-artefact/`, `171-reach-prefix/`.

## The observation that opens this

**188 found, while investigating something else, that the reference's own BH-significant gene set is
PREDICTABLE from quantities we know exactly.** Its pool-reconstruction check reports that
`|lfc| · √cpm` reproduces the reference's significance-selected pool at **0.7219**, against `|h|` at
0.6406 and `|lfc|` alone at **0.3796**, at a **median pool size of 89** — which independently
matches audit 181's board-implied `n_real` of 90–110.

That is not a curiosity, it is the missing key. The reference's calls are driven by **statistical
power × effect size**, and we hold both inputs: the control CPM is measured exactly from the
challenge's own 18,400 NTC cells per context, and `lfc` is our prediction. **We can estimate which
genes the reference will call.**

**185 measured the call-set frontier with `|H|` as the key and found `|H|` "no better than random
for `jac`".** That is consistent with 188 — `|h|` scores 0.6406 and `|lfc|` 0.3796, and neither is
the power-weighted product. **185 measured the right frontier with the wrong key, and 188 found the
right key hours later for an unrelated reason. Nobody has combined them.** That is your job.

## Why this is worth a full run

`fid = k / max(n_pred, n_real)` and `jac = |R ∩ R̂| / |R ∪ R̂|`. We currently call ~4,474 genes per
(context, perturbation) against a reference set of ~89–200. Solving our own board row: our 4,474
calls contain ≈123 of the reference's genes — **61% recall at 3.1% precision** — and `fid` 0.5038 is
simply the sign accuracy of all 4,474. At `n_pred = n_real` the yield cap pins at 1 and `fid`
becomes the accuracy of the **head**, and `jac` becomes `p/(2−p)` in the head's precision.

We rank **169th (`jac`) and 176th (`fid`) of 612**, and the field's best are `jac` 0.0374 and `fid`
0.5192. **But the ceiling here is not the field's best** — `jac`'s replicate anchor is 0.3972, i.e.
a real repeat experiment agrees with itself on DE membership only ~40% of the time, and nobody on
the board is within 10× of it.

## Work

### Phase 0 — known-answer rows (mandatory, abort on failure)
- Reproduce 188's pool-reconstruction numbers: `|lfc|·√cpm` **0.7219**, `|h|` **0.6406**, `|lfc|`
  **0.3796**, median pool size **89**. State exactly what statistic "reproduces at 0.7219" is
  (rank correlation? AUC? recall at n?) — 188's own definition, not one you invent.
- h1-test: production `reach` **0.18283**, `n_pred` **6035**, `jac` **0.06824**, `fid` **0.50701**
  (185 §3b). Donors `k562_gwps` + `orion_hek293t` only — `data/donors/h1_2025` covers 85/85 of
  h1-test's targets and is pure leakage.

### Phase 1 — build the membership predictor properly
Predict, per (perturbation, gene), the probability that the **reference** calls it significant.
Inputs must all be available at build time: our predicted `lfc`, the control CPM (exact), our
predicted effect in `h` space, the control cells' per-gene detection rate and dispersion, and the
reference's own power law. Start from 188's `|lfc|·√cpm` and ask whether a fitted form beats it —
fit **leave-one-target-out** on h1-test and, separately, on the Orion surfaces, and report both.

Report the **enrichment curve**: `p(n)` = fraction of your top `n` that the reference calls, and
the fold-enrichment over `n_real/n_tested`, at n = 25, 50, 89, 150, 300, 600, 1200, and at the true
per-perturbation `n_real_p`. **This curve is the whole result** — everything below is arithmetic on it.

### Phase 2 — the joint (`fid`, `jac`, `reach`) optimum over the call set
With the call set free (191 proved the production emitter realises any prescribed set at
bit-identical pseudobulk), compute all three members as a function of `n_pred` and of the key.
Remember:
- `n_real` is a **distribution**, not a scalar: audit 181 puts the board median at 90–110 with a
  heavy tail, and 188's independent pool reconstruction says median **89**. Report against a swept
  level, preserving the shape — **scale each perturbation's pool proportionally, never
  `Surface(n_real_target=N)`**, which deletes the small-`n_real` stratum carrying 98% of our reach.
- The penalty is **one-sided**: flat above `n_real_p`, linear below (093/v14 lost −0.0355 board
  `fid` exactly this way). So the expected-value optimum sits **above** the predicted `n_real_p` by
  roughly the predictor's error sd. Size that margin explicitly; it is the design.
- `reach` is a no-op in the call set (171, proved) but shares the emitter, so report it as a check,
  not as a target. 189 and 190 closed the ordering-key family; do not reopen it.

### Phase 3 — falsifiers, all mandatory
1. **Shuffled-key control** at every `n`. The enrichment must beat it.
2. **Retrodict 093/v14**: `--keep-top 850` gave board `fid` **−0.0355**, `jac` **+0.0043**,
   `reach` +0.0137. Your model must reproduce the *sign and rough size* of all three given v14's
   realised `n_pred`. Audit-181 rule 1: this family has a board point and it must be retrodicted
   before your numbers are read.
3. **Retrodict a21**: `n_pred` 4,420 → 1,141 with board `fid` 0.5034 → 0.4945 and `jac` −0.017.
4. **Leakage audit**: state how many of the 300 panel targets are in any fitting pool, and show the
   result survives deleting them.

### Phase 4 — price it
Six-member ledger in overall units. `pds` and `nmae` are pinned by the emitter (191 proved the
per-gene panel totals are digest-identical on 900/900 blocks); `mse` through 096's `claim` route.
Rule 13 charges emission-law arms −0.0146 scaled `reach` — **a board point for the
profile-preserving class is in flight right now** (`rc_hybrid`, submitted 2026-09-03 11:20 UTC,
preregistration in `experiments/185-rank-controlled-emission/`). Report your ledger **both ways**,
with and without the toll, and say which board outcome would make the arm shippable.

## Deliverable
`RESULTS.md`: known-answer rows, the enrichment curve, the three-member frontier, the sizing margin,
all four falsifiers, and the two-way ledger. A negative with the enrichment curve in it is a full
result — that curve is the number this project has never had.

## Territory
193 owns the whitening operator. 194 owns the donor tables. 195 audits. You own the **call set**.
