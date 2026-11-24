# 198 — Which NO-GOs were closed on an instrument we now know is WRONG? A systematic re-audit of 190 closures.

You are an adversarial audit subagent. Work ONLY in `experiments/198-closure-reaudit/`.
Read-only on everything else. No `vcc`, no `cell-eval2`, no builds, no submission.
Peak RSS ≤ 8 GB. `uv run python`. This is mostly reading and arithmetic on numbers already on disk.

## Why this exists

**Twice in the last twelve hours a closed result turned out to have been closed on the wrong
measurement**, and both times the correction was worth more than the original experiment:

- **188** found that 169's `nmae` verdict — a 56× disagreement between two harnesses that had
  blocked a decision for a day — was **DE-pool composition**, not instrument bias. The reference's
  BH pool is *significance*-selected; `|lfc|·√cpm` reproduces it at **0.7219** (`|h|` 0.6406,
  `|lfc|` 0.3796) at median size **89**. Priced on a correctly-built pool, a third surface agreed
  with the one we trusted. Deep or magnitude-selected pools **invert the sign** of the `dU` term.
- **188 Phase 3** found that `debench.Surface.pds()` and PLOLO disagree in *sign* on the same arm,
  and the cause is the **space**, localised to the 5–10 CPM band (11.0% of genes, 21.4% of the
  `|lfc|` norm). Converting both sides of the same table into the scorer's `h` space turned
  −0.02339 into **+0.00504**.

`AGENTS.md` rule 3 says *"validate the validator"*. This experiment applies it retroactively.

## The task

**Sweep every `experiments/*/RESULTS.md` and identify every closure whose verdict rests on a
measurement that a later result has invalidated.** Then rank the reopenable ones by how much the
correction could plausibly be worth.

### The invalidating findings, each of which retro-condemns a class of measurement

| finding | source | what it invalidates |
|---|---|---|
| The reference's DE pool is **significance**-selected, median 89; magnitude-selected pools invert `nmae`'s sign | 188 | **any `nmae` or DE-pool number built by taking the top-N genes by magnitude** |
| `debench.Surface.pds()` must never size a `pds` arm; the disagreement is the space, in the 5–10 CPM band | 188 Phase 3 | **any `pds` verdict taken from an h1-test DE table rather than PLOLO** |
| The Orion thinned harnesses carry a target-independent blob (26.4% / 44.1% of truth energy) that is the long-gene/nuclear-fraction technical axis | `AGENTS.md` harness warning, 159, 175/178 | **any absolute `pds`, transfer cosine, `fid`, `reach`, `jac` or `n_real` read off an Orion harness** |
| `n_real` is a **designed** distribution spanning orders of magnitude (the panel was binned by DE count), median ~89–110, not a scalar | `AGENTS.md` §2.1, audit 181, 188 | **any conclusion that assumed a uniform or single-valued `n_real`, and every use of `Surface(n_real_target=N)`, which deletes the small-`n_real` stratum carrying 98% of our reach** |
| 093's `n_real ≈ 2,670` is refuted (board reach would read 0.0034) | audit 181 | **anything downstream of that inference** |
| h1-test is **not** validated for `reach` on emission-law changes (a21/a22 invert it); it **is** for spread-policy changes | audit 181 | reach verdicts on emission-law arms |
| Rule 13's `−0.0146` reach toll is an **emission-law** constant; 169 falsified it for profile-class arms | audit 181, 169 | **profile-class arms vetoed by a toll that does not apply to them** |
| Cross-line adaptation is closed in the strong form on 50 lines — the line-specific residual is not a basal-conditioned gain | 184 | conclusions that kept that family open, and conversely any arm rejected *because* it did not adapt |
| Global amplitude is at the joint `mse`+`nmae` optimum; `nmae` is ~93% direction | 187 | arms whose verdict rested on an amplitude argument |
| The ordering key family is closed (`|H|` is best; the 0.90 purity bar caps it) but the **delivery** is real and worth 95% of the counterfactual | 189, 190, 185 | reach verdicts that conflated key quality with delivery |

### Method

1. **Enumerate.** Build a table of every experiment with a NO-GO/terminal verdict, the endpoint its
   verdict rested on, and the surface that endpoint was measured on. Machine-readable, from the
   files — do not work from memory or from `CURRENT_STATE`'s summaries.
2. **Flag.** Mark each against the table above. Report counts per class.
3. **Rank the reopenable.** For each flagged closure, estimate what the corrected measurement would
   plausibly read and what it would be worth in **overall board units**, using the scoring map in
   `docs/CURRENT_STATE.md` (`dOverall/dRaw`: `pds` +0.3690, `mse` −0.1752, `nmae` −0.2755,
   `fid` +0.5607, `reach` +0.1882, `jac` +0.4544). Be explicit about which are cheap to re-measure.
4. **Adversarial pass on the LIVE work.** Six experiments are running (192 call set, 193 operator
   eigenvalue map, 194 donor re-estimation, 195 learned retrieval transform, 196 emitted depth,
   197 cosine inversion) and one arm is on the board (`rc_hybrid`, preregistration in
   `experiments/185-rank-controlled-emission/`). **Read their mandates and try to break their
   premises.** Any premise you can falsify from data already on disk is worth more than a ranked
   list.
5. **Two specific candidates to adjudicate**, named because the coordinator flagged them and did
   not resolve them:
   - **125** (hierarchical panel posterior) was rejected for losing "signed credit" while gaining
     **+0.0122 top-20% sign precision**. Is "signed credit" any of the six scored members? If it is
     not, 125 was closed on an endpoint the competition does not read, and 125's own text says
     *"a future experiment could test a direction-only confidence/ranking layer while retaining
     122's effect magnitudes"*.
   - **AgingLab's `GeroAI_v9 → v11`** gained **+0.1002 raw board `reach`** with `pds` going *down*
     0.008 — the largest ordering-shaped move on the board — while **172** tried their published
     feature family and got nothing. Reconcile, or state precisely what evidence would.

## Rules
- **Quote the number and the file for every claim.** A closure you cannot cite is not a finding.
- Where a closure is *not* affected, say so — the value of this audit is as much in confirming the
  map as in reopening it.
- Do not reopen a family on taste. The bar is a *specific* measurement that the invalidating finding
  shows was wrong, plus a plausible corrected value.

## Deliverable
`RESULTS.md`: the full enumeration table, the flag counts per class, a ranked list of reopenable
closures with expected value in overall units and re-measurement cost, the adversarial pass on the
six live mandates, and the two named adjudications. A finding of "every closure survives" is a full
result and would be the strongest statement anyone has made about this project's map.
