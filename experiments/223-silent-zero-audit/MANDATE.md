# 223 — Twenty-one partial-coverage donor experiments may have been measured at a dose of EXACTLY ZERO

You are an audit + measurement subagent, in the tradition of 181, 198 and 219. Work ONLY in
`experiments/223-silent-zero-audit/`. Read-only elsewhere. No `vcc`, no builds, no `cell-eval2`.
Peak RSS ≤ 12 GB. `uv run python`.

## READ FIRST
1. `experiments/218-independent-donors/RESULTS.md` **§7.1** — the mechanism. It is short and it is
   the whole basis of this mandate. **Verify it in the source before doing anything else.**
2. `experiments/215-leader-decomposition/RESULTS.md` §4.2 and its Phase-3 verdict — why this matters
   strategically.
3. `AGENTS.md` §8.4 (as corrected by 219: the honest multiplier is **2.51×**, and §8.2's single-arm
   bar is **+11.6–13.6%** amplitude / **+2.5–2.9%** shape-only, not +31%).
4. `scripts/make_submission_pm_cd4.py` — `solve_cd4_fraction`, around line 1258.

## The bug, and why it invalidates a class of NO-GOs

`solve_cd4_fraction` sizes an added residual by the ratio of its **median row norm over every panel
row** to the whitened base's. **A source covering ≤ 50% of the 300-gene panel has a median row norm
of exactly 0, so `f_per_unit = 0` and the solved fraction is exactly 0.**

Production *raises* on a fully-zero residual, so a CD4 build fails loudly. **A partially-covering
source is not fully zero — it is, say, 21% covered with a zero median — and nothing anywhere catches
it.** 218 reproduced this on a third surface: its `telohaec` and `gse295214_h1` arms returned
`fraction 0.00000` with endpoints **bit-identical to the baseline at all six doses on all four
folds**. 200 found the same class on CD4; 207 verified it on h1-test.

**So any experiment that evaluated a sub-50%-coverage source through this path measured a dose of
nothing, and its NO-GO is void — it tested the null, not the arm.**

218 already showed what the fix is worth on two such sources: sized on **carrying rows**
(`median(‖added_i‖ : ‖added_i‖ > 0)`), the Pan-2026 H1 + `h1_2025` union (77/300) prices at
**+0.00919 / +0.01054**, **5.5× the gate**, with `pds` going **UP** — which no CD4 dose arm has ever
achieved. Experiment 222 is building it now.

## And the strategic reason this is urgent

215 decomposed the leaders' trajectories and found their direction **arrived in MANY SMALL STEPS
inside the readable window** — `cqawesome` in **5 separately resolved steps**, largest only 28% of
the total; `Haystack` in 3. **Nobody shows a single step carrying their direction.** Our own history
is the opposite: one big step (the whitening operator) and then nothing.

**A stack of many partial-coverage sources, each worth a few percent of cosine, is exactly what
"many small steps" looks like** — and it is precisely the thing our sizing rule has been silently
discarding. **That is a hypothesis, not a finding, and your job is to test the half of it we can
test: how much is sitting in our own closed experiments.**

## Work

### Phase 0 — verify the mechanism (mandatory, abort if it fails)
Read `solve_cd4_fraction` and confirm the median is taken over **all** panel rows. Then construct a
minimal reproduction: a synthetic residual covering 21% of the panel, and show its solved fraction
is **exactly 0** and the resulting effect is **bit-identical** to the baseline. If it is not, this
mandate is void — say so and stop.

### Phase 1 — classify all 21 candidate experiments
These evaluated a partial-coverage source:
`055-kolf-direct`, `057-kolf-union-whiten`, `059-feng-ipsc`, `079-kolf-local-centering`,
`085-kolf-pm-gate`, `105-telohaec-residual`, `107-jurkat-high-moi`, `115-gse295214-h1`,
`116-dld1-lowmoi`, `119-altair-a549`, `120-myllia-a549`, `121-a549-cross-source-concordance`,
`126-hcasmc-panel-confirmation`, `129-gse295214-transfer`, `130-telohaec-panel-confirmation`,
`131-hues8-hepatocyte-metadata`, `135-telohaec-atlas-confirmation`, `137-jurkat-tangent-cone`,
`138-dld1-replicate-sign`, `139-jurkat-donor-simplex`, `176-cmap-donors`.

For each, determine from its own code and logs: **(a)** its source's panel coverage, computed not
eyeballed; **(b)** whether it went through `solve_cd4_fraction` or an equivalent median-sized path;
**(c)** whether its reported effect was **bit-identical to its baseline** (the signature of the
silent zero — 218 found `fraction 0.00000` recorded in its own artifacts); **(d)** its verdict and
margin.

Classify: **VOID** (measured at zero — the NO-GO tested nothing), **VALID** (a different path, or
coverage > 50%), **UNKNOWN** (cannot tell from what is on disk — say so rather than guessing).

**`105-telohaec-residual` is the one to check first** — its name says *residual*, and 218 measured
`telohaec` returning `fraction 0.00000`.

### Phase 2 — re-measure the VOID ones, on carrying rows
For every VOID source with a built donor table, re-run 218's construction — un-whitened
post-whitening residual, **sized on carrying rows**, swept over dose — on 162's four PLOLO folds.
Report **Δcos (primary, ±0.20% relative per 204)**, **Δ`pds`**, and `q_w`, per source and for the
**union of everything that survives**.

**Price under §8.4's corrected rule and never carry a member at zero.** Compare against the
corrected bar (**+2.5–2.9%** shape-only), not the retired +31%.

### Phase 3 — falsifiers, all mandatory
1. **Shuffled-target control per source, at every dose.** 218's was strictly negative at every dose;
   204's and 209's arms both died because a permutation control tied or beat them. **A source whose
   shuffle ties it is dead regardless of coverage.**
2. **Non-carrying rows must not move** (218 measured |Δ| ≤ 0.018 pp). If they do, the arm is
   row-reweighting, not adding information.
3. **The row-matched CD4 control**: same rows, same sizing, CD4's own data. 218 measured it
   reproducing **1/3.9** of the union's gain *and costing* `pds`. Any source that does not beat its
   own row-matched CD4 control is adding nothing new.
4. **Leak audit.** Several of these sources are H1 or iPSC and `data/donors/h1_2025` covers 85/85 of
   h1-test's targets. State the donor↔surface overlap for every arm; PLOLO's challenge frame is the
   instrument, not h1-test.
5. **Do not double-count 222.** Pan-2026 H1 and `h1_2025` are already being built; report the
   *marginal* value of each additional source **on top of** that union.

## Deliverable
`RESULTS.md`: the Phase-0 reproduction, the 21-row classification table with coverage and verdict
per experiment, the re-measured VOID sources with all four falsifiers, and a ranked shortlist of
what is worth adding to 222's union — each with its marginal price and what it would cost to build.

**"All 21 were VALID and nothing was measured at zero" is a full result** and would close this in an
afternoon. **"N of them were void" is a much bigger one**, and it would mean this project has been
rejecting donor sources on a null for weeks.

## Territory
216 CD4 dose/form (closed it at `f = 0.34`) · 220 ordering keys · 221 cloud scorer ·
222 is building the Pan/h1_2025 union. **You own the silent-zero audit and the partial-coverage
sources.** Do not touch CD4's dose or the builder — 222 owns the sizing-rule change.
