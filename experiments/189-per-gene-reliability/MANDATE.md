# 189 — The per-gene cross-line reliability prior: the one confidence signal 172 excluded by design

You are a research subagent. Work ONLY in `experiments/189-per-gene-reliability/`.
No `vcc`, no `cell-eval2`, no submission build, no writes outside your directory.
Peak RSS ≤ 10 GB (four other agents are running). `uv run python`.

## READ FIRST
`AGENTS.md` §4, §11.1 rules 1–13 with the audit-181 amendments. Then, and you must know each of
these well enough not to repeat it:
- `experiments/172-sign-confidence-key/RESULTS.md` — learned confidence keys, **NO-GO**.
- `experiments/171-reach-prefix/RESULTS.md` §12c — the p-order is expression-dominated.
- `experiments/154-confidence-key-screen/RESULTS.md` — "reach is model-limited, not order-limited".
- `experiments/161-crossline-residual/RESULTS.md` — the residual is not predictable **from context**.
- `experiments/187-joint-amplitude/RESULTS.md` §5 — `nmae` moves only through the **|t|-weighted**
  sign accuracy `q_w` (|p|-weighted 0.62 on h1 against count accuracy 0.55).
- `docs/BOARD-LEDGER.md` — our `reach` trajectory and the shrunk-cohort control.

## The gap in the record

172 learned a per-(target, gene) sign-confidence key from the three production donors and closed the
family. Its feature set was, in its own words, *"four measurable scalars … **no gene identity and no
target identity**"* — it deliberately mimicked `AgingLab`'s public description. **So a prior attached
to the GENE has never been fitted here**, and 154's screen predates the `|H|`-ordering diagnosis.

That prior is estimable at a sample size nothing else in this project enjoys. `orion_hct116` and
`orion_hek293t` share ~17,900 perturbation targets; with `k562_gwps` the three-way intersection is
**9,073 targets**. For each gene that is ~9,000 paired cross-line observations of whether the
response sign agrees — an extremely well-determined per-gene quantity, fitted entirely on donor
lines and therefore **leakage-safe for the challenge and for D/E/F by construction**.

## Pilot — reproduce these first, they are your known-answer rows

Run by the coordinator, 2026-09-03 02:55 UTC. Deterministic. Donor `delta.npy` matrices, common gene
axis (7,976 genes), three-way target intersection (9,073 targets), **each matrix centred across
targets** (removing the shared programme, as production does), a per-gene sign-agreement statistic
gated to entries above each line's own per-gene median |effect|, genes with >50 supported pairs:

| quantity | value |
|---|---|
| common genes / three-way targets | **7,976 / 9,073** |
| genes with support | **7,916** |
| per-gene `r`(HCT116, HEK293T) | mean **0.5379**, sd **0.0467**, p5 0.499, p50 0.530, p95 0.594 |
| spearman( r(HCT,HEK), r(HCT,K562) ) | **+0.3833** |
| spearman( r(HCT,HEK), r(HEK,K562) ) | **+0.3359** |
| spearman( r(HCT,K562), r(HEK,K562) ) | **+0.4763** |
| top vs bottom decile of r(HCT,HEK), read on **held-out** r(HCT,K562) | 0.5515 vs 0.5113, spread **+0.0402** |

And the confounder test, which is the interesting half:

| quantity | value |
|---|---|
| spearman( r_HE, control expression `b_ctrl` ) | **+0.5215** |
| spearman( r_HE, mean \|effect\| ) | +0.1125 |
| spearman( r_HE, across-target sd ) | +0.1578 |
| **after partialling all three out (rank space):** spearman( resid r_HE, resid r_HK ) | **+0.2133** |
| same, resid r_EK | **+0.2251** |
| held-out decile spread of the **residual** prior | **+0.0173** |
| held-out decile spread of **expression alone** | **+0.0505** |

**Two findings, and the second is the surprising one.** (1) A per-gene reliability prior exists,
generalises to a held-out line, and survives partialling out expression, magnitude and spread —
worth +0.017 held-out sign accuracy between deciles on top of everything observable. (2) **Control
expression is by itself the strongest per-gene reliability signal on the board of candidates,
worth +0.0505 between deciles.**

(2) sits in direct tension with 171, which measured that our emitted p-order is **expression-set**
(`rho −0.391` against `|effect|` inside the call set) and treated that as the defect to remove. If
expression genuinely predicts which genes transfer, then part of what 171 called a defect is signal,
and the right key is a **composite** rather than either one alone. **Resolving that tension is the
first thing you should do**, because it decides whether 185's delivery work should target `|H|` or
something else. Note 172 reports that no ad-hoc composite (× or ÷ expression, centred, z-scored,
squared) beat `|H|` on `reach`; a *fitted* prior is a different object, and you must beat `|H|`'s
**0.23533** on the same instrument to claim anything.

## Work

### Phase 0 — reproduce the pilot table above. If you cannot, abort and say so.
Also reproduce `|H|`'s h1-test reach **0.23533** and production's **0.18283** (172/171) before
quoting any reach number.

### Phase 1 — build the prior properly
Estimate `r_g` leave-one-LINE-out so the prior used to score a line never saw it. Consider adding
`nadig_jurkat`, `rpe1` and `k562_essential` as extra *fitting* lines — 161 used a five-line core
(K562/RPE1/Jurkat/HCT116/HEK293T, 1,883 shared targets, none on the challenge panel) and its
`inventory.py`/`common.py` already assemble it; reuse that rather than rebuilding. Report how the
prior's held-out predictive power scales with the number of fitting lines.

Shrink `r_g` properly — it is a binomial proportion per gene with known `n`, so use an empirical-Bayes
posterior mean rather than the raw fraction, and say what the shrinkage does to the decile spread.

### Phase 2 — does it buy anything on the members?
Three uses, priced separately. **You own the KEY and the CALL SET composition; 185 owns the
mechanism that delivers an order; 188 owns per-gene magnitude.** Do not build an emitter.

1. **As an ordering key.** `reach(key)` on h1-test with `key = f(|H_g|, r_g, expression_g)`, against
   the `|H|` incumbent **0.23533** and the H-oracle **0.62119**. 172's calibration applies: *reach
   reads only the front, so precision at the top is the objective, not rank correlation* — report
   top-1% and top-5% sign accuracy over the pool alongside reach.
2. **As a call-set filter.** Sign precision of the top-`n` genes as a function of `n` and of the
   prior, which is what `fid = k/max(n_pred, n_real)` reads. Coordinate with
   `experiments/185-rank-controlled-emission/INBOX-1-per-perturbation-sizing.md` — read it, it holds
   the board arithmetic bounding `n_real`, and do not duplicate its per-perturbation sizing work.
3. **As a weight on `q_w`.** 187 §5 shows `nmae` moves only through the |t|-weighted sign accuracy.
   Measure whether down-weighting low-`r_g` genes raises `q_w` on h1-test, and by how much. Do not
   change magnitudes to do it — report the achievable `q_w` and hand the magnitude question to 188.

### Phase 3 — the honest falsifiers, all of which must be run
- **Permutation null:** shuffle `r_g` across genes and re-run every endpoint. The prior must beat it.
- **Expression-only control:** the whole pipeline with `r_g` replaced by control expression. If
  expression alone does as well, the finding is "use expression", which is still a result but a
  different and cheaper one — report it as the headline in that case.
- **The h1-test surface's own trap:** `data/donors/h1_2025` covers 85/85 of h1-test's targets and is
  pure oracle leakage. Use `k562_gwps` + `orion_hek293t` only on h1-test, and audit donor↔surface
  target overlap before every screen.
- **Panel disjointness:** state how many of the 300 challenge panel targets are in your fitting pool.
  161's core5 pool has zero by construction; keep that property or justify losing it.

## Deliverable
`RESULTS.md`: the reproduced pilot, the prior's construction and its LOLO scaling, the three
priced uses against their incumbents, all three falsifiers, and a verdict. If the answer is
"expression alone, and 171's defect is half signal", say that plainly — it would redirect 185.

## Territory
184 owns Tahoe. 185 owns the emission law and how an order is delivered. 186 owns intel.
188 owns per-gene magnitude. **You own the per-gene PRIOR and which order to deliver.**
