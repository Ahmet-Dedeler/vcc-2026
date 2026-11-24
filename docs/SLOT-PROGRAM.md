# The slot program — the board is a measurement instrument with ~94 uses left

> ## 🔄 REVISION 2026-09-04 01:20 UTC — the program worked, and item 2 is now closed. Read this first.
>
> **Three slots spent, +0.0308331: a15 `0.1740865` → `cd4f034` `0.2049196`, rank 24 → 7 of 647.**
> Gap to `cqawesome` **0.04383 → 0.01590**. Every one of the three was a designed measurement that
> also happened to score, which is what §1 argued for.
>
> | slot | arm | Δ | what it measured |
> |---|---|---:|---|
> | 09-03 | `rc_hybrid` | **+0.0093325** | **item 2 — CLOSED.** Rule 13's −0.0146 `reach` toll does **not** apply to profile-preserving arms; it was charged and the arm realised **+0.054940**. And the h1-test → board `reach` factor for designed-order arms is **1.290**, not 0.4277 — h1-test *understates* by 3× (item 4, first real answer) |
> | 09-03 | `rc_callset` | **+0.0084448** | items 1 and 3 — `jac` went **positive for the first time**, `n_pred` 4,480 → 2,599, and `pds`/`mse` reproduced to **sixteen significant figures** for the second time |
> | 09-04 | `cd4f034` | **+0.0130557** | **a new item the program did not anticipate** — the first arm ever to move the **direction** channel on the board (`cos` 0.07228 → **0.09027**), and the post-mortem that produced `docs/FINDINGS.md` §8.4 |
>
> ### The program's own lesson, and it outranks the items below
> `cd4f034` predicted +0.00437 and delivered **+0.0130557**, because **we priced it on `mse` alone
> and carried `nmae`/`fid`/`reach`/`jac` at zero — and those four delivered 60% of the gain.**
> §8.4 (as corrected by 219) puts the honest multiplier at **2.51×**, and 219 additionally showed
> **§8.2's "+31% cosine" single-arm bar was itself `mse`-only and 2.0–2.3× too high** — corrected to
> **+11.6–13.6%** (amplitude floor) and **+2.5–2.9%** (shape-only).
>
> **Consequence for every future slot: price all six members, always.** The rule-12 exemption
> ("unmeasured members cannot move") is for `spread`/`phi`/`keep_top`/call-set arms where the
> pseudobulk is bit-identical. **A direction arm is the opposite case and the exemption must never
> be claimed for one.**
>
> ### The revised priority order
> | # | what | status |
> |---|---|---|
> | **1** | **the CD4 dose ladder** — the only lever that has ever moved the direction channel | `f = 0.34` **done** (+0.0131); `f = 0.80` **building**; branch tree pre-committed in `experiments/217-cd4-f080/DECISION-TREE.md`. Note 219: the propagation factor **falls with dose** (`pds` turns over at `f ≈ 0.20`), so every extrapolation must use the **local** factor |
> | **2** | **a fifth independent donor line** | re-priced by 219 from +0.0062 to **+0.0140–0.0163** — from below the old bar to **1.6–1.8× above** the corrected one. Experiment 218 is searching under the right criterion (166's was *"genome-wide covering the panel"*, which **would have missed the CD4 atlas**) |
> | **3** | old items 1 and 3 (`n_real`, `fid`/`jac` vs `n_pred`) | **largely answered without spending slots on them**: 192 measured the board's `n_real` distribution from the existing `fid` series (median 31–45, mean 270–383), and `rc_callset` mapped the call-set response. Do not spend a dedicated ladder on these |
> | **4** | the per-member transfer factors | now partly answered (`reach` ×1.290; `jac` ×0.44 from 181; PerturbMatch ×1.533 from 214), and **experiment 221's cloud `cell-eval2` worker will produce the full six-member table** from the four board-anchored artifacts we already hold — **no slots required** |
>
> **The standing rule is unchanged and was honoured on all three:** one variable per slot, a
> preregistration committed before the build, and falsifiers named in advance. All three arms fired
> at least one falsifier, every time in the favourable direction, and each firing taught more than
> the score did.



Written 2026-09-03. **Read with `AGENTS.md` §6's banner (never spend a slot on rank) and §11.1
rule 13 (the composed gate).** This does not relax either. It says what the *other* legitimate use
of a slot is, and designs it.

## 1. Why a program rather than one arm at a time

Three of the six members — `fid`, `jac`, and `reach`'s normalisation — are **unreadable on every
local harness** (rule 9: `n_real` is 68–92 on the repaired Orion harness, 222 on h1-test, ~89–110
on the challenge; the challenge is in none of those regimes). Audit 181 additionally removed
h1-test's validation for `reach` on emission-law changes. **So for half the score, the board is the
only instrument we have.**

Against that: **2 slots/day, 49 days to Oct 22, ≈96 slots remaining**, of which we have used 25 in
six weeks — about **one every other slot**. And the scorer is **deterministic**: a byte-identical
resubmission reproduced every member to float representation, so a board delta carries no
run-to-run variance and a single point is a real measurement, not a sample.

**A slot spent measuring something no harness can read is a slot spent on first place.** A slot
spent re-scoring something we already know is not, and the banner still forbids it.

## 2. What only the board can tell us, in priority order

| # | quantity | why it matters | how one slot measures it |
|---|---|---|---|
| 1 | **the `n_real` distribution** | `fid = k/max(n_pred, n_real)`; the knee in `fid(n_pred)` **is** `n_real`. It has already cost four board predictions (rule 9) and 093's estimate of 2,670 was refuted by audit 181 | an emission-only arm at a prescribed `n_pred` well below a15's 4,474. Two or three points bracket the knee |
| 2 | **the emission-law `reach` toll for PROFILE-PRESERVING arms** | rule 13 charges −0.0146 scaled on the strength of a21+a22, both of which *also* moved the profile. It currently vetoes an entire arm class | `rc_hybrid` — **in flight, submitted 2026-09-03 11:20 UTC** |
| 3 | **`fid` and `jac` versus `n_pred`** | the two members where we rank 176th and 169th of 612, and the yield cap makes the response one-sided and therefore cheap to map | the same ladder as (1); every point reads both |
| 4 | **the h1-test → board transfer factor, per member per arm class** | every local prediction passes through it. Only two independent estimates exist: 164's **0.4277** and the a8/a15 pair's **0.390** (`docs/BOARD-LEDGER.md` §4), and the second is confounded | any arm with a clean h1-test prediction; accumulates for free from every other probe |

## 3. The design rules

1. **One variable per slot** (rule 5). v3 moved two and its −0.018 was unattributable for a day.
2. **Prefer emission-only arms for measurement.** They hold the per-gene panel totals
   bit-identical, so `pds` moves by **exactly zero** — board-verified twice, on v14 and a8 — and
   `nmae` to five decimals. 191 proved the production emitter can realise any prescribed call set
   and any prescribed order under that constraint (digests identical on 900/900 blocks). **Every
   board movement is then attributable to the DE members alone**, which is what makes a single
   point a measurement rather than a data point.
3. **Preregister before building**, with the two hypotheses the slot discriminates and the
   falsifiers, and commit it. The preregistration for `rc_hybrid` is the template.
4. **Predict every member**, or name the ones you cannot and the mechanism that pins them (rule 12).
5. **Retrodict the class's existing board points first** (audit 181 rule 1). For the call-set family
   those are 093/v14 (`fid` −0.0355, `jac` +0.0043, `reach` +0.0137 at `--keep-top 850`) and a21
   (`n_pred` 4,420 → 1,141, `fid` 0.5034 → 0.4945, `jac` −0.017).
6. **Joke model name, empty description** (§14), on every one of them.

## 4. What this is not

It is **not** licence to submit small improvements. The composed gate still governs *shipping*:
predicted > +0.0019 shape-only, > +0.0089 amplitude-changing. A measurement slot is justified by
the information, is preregistered as a measurement, and is expected to be roughly break-even on the
displayed row — which is informational only, since no leaderboard entry is prize-eligible and only
the last entry after Oct 22 counts.

**And it stops the moment the answer is known.** Once the `n_real` distribution is pinned, further
ladder points buy nothing and the banner applies again in full.

## 5. Standing sequence

- **Now:** `rc_hybrid` (row 2), in flight.
- **Next, pending experiment 192's enrichment curve:** the first `n_real` ladder point, at the
  `n_pred` 192's curve says is closest to the median. 192 owns choosing it; do not guess it.
- **Then:** two more ladder points bracketing the knee, one above and one below.
- **Throughout:** every arm's h1-test prediction goes in `docs/BOARD-LEDGER.md` alongside its board
  result, so row 4 accumulates without spending a slot on it.

**Before Oct 22 the ladder must be finished**, because `n_real` cannot be measured on D/E/F at all
(transition-plan item 5) and the call-set policy for the final round has to be decided from A/B/C
evidence. That is a deadline, not a preference.

---

# 6. REVISION 2026-09-04 — the source ladder. Branch A fired; this supersedes §5's sequence.

§5's `n_real` ladder is **complete** — 192 closed `fid` to the call set and 220 bounded the whole
ordering channel at **+0.0077**, less than one direction arm. The ordering programme is done.
**(Updated 2026-09-04 by 231: that bound is itself now stale — our own `h1union` row moved `reach`
to world rank 7 of 656, so the remaining prize is **+0.002689**, below the **0.004502** ordering
noise. `docs/FINDINGS.md` §8.1c: do not open another ordering arm.)**

## 6.1 What changed

`h1union` scored **+0.0069114** (board 0.2118310, **rank 3 of 656**), **all six members up**, on the
single variable of adding a fourth donor estimate. That is the first board evidence that **stacking
sources works**, and it fired `DECISION-TREE-AMENDED.md`'s branch A: *rule and sources validated.*

## 6.2 The budget nobody has been counting

**~48 days × 2 slots = ~96 slots before Oct 22, and NOT ONE of them is prize-eligible.** Only the
last entry after Oct 22 counts. So every pre-Oct-22 slot is **purely an instrument**, and the
constraint on spending them is not risk — it is our ability to *learn* something from each one.

**We have been rationing a resource that expires worthless.** Three slots in two days is
under-spending against a 96-slot budget.

**Before the first ladder build, read `docs/NEXT-BUILD-SPEC.md`** — the `VCC_H1_UNION_G` machinery
is fully source-agnostic except for one hardcoded source tuple, which must become registry-driven
before any new source can ship. The two known-answer hashes that govern that change are in §3 there.

## 6.3 The programme: one source per slot, greedily, measured against the shipped stack

1. **Each slot adds exactly one source** to the shipped stack and measures its **marginal** board
   value. One variable per slot, as always.
2. **Order by 225's ranking**, which is by *projected pessimistic-draw October coverage*
   (`INBOX-OCTOBER.md`), then by leave-top-3-out retention (`INBOX-CONCENTRATION.md`), then by
   headline marginal. **Not by today's panel overlap** — that is fitted to a panel about to be
   replaced.
3. **Stop adding a source the moment two consecutive additions land below the shape-only gate**
   (+2.5–2.9% relative cosine, §8.4). That is saturation, and it is the number 225's Phase 4 is
   built to find; the board is the arbiter, not PLOLO.
4. **Reserve the last four pre-Oct-22 slots** for the re-rehearsal artifact (`DEF-RUNBOOK` §5.-1)
   and for whatever the transition breaks. Do not arrive at Oct 22 having spent everything.

## 6.4 Why this beats holding slots back

Every source we do not board-test before Oct 22 is a source we must ship **blind** into D/E/F on the
one entry that counts. The PLOLO→board factor is **not identified** (`AGENTS.md`: 1.384 vs 1.000 on
two arms in the same family), so local pricing alone cannot substitute for a board point. **Each
slot spent now converts an unidentified factor into a measured one on that specific source.**

## 6.5 Unchanged

Joke model name, **empty description**, on every one (§14). Never re-submit a scored package. Never
spend a slot to improve the displayed row. Preregister every arm before it builds, and **reissue any
band before the result lands, never after** — 2026-09-04 is the worked example in both directions:
the `kolf2` correction was right, and the transfer-factor correction applied on top of it was wrong.

---

# 7. REVISION 2026-09-04 05:25 — §6's source ladder is CLOSED. What the ~90 remaining slots are actually for.

**§6.3's "one source per slot, greedily" is obsolete.** 225 swept the space exhaustively (GEO,
BioStudies, CZI, Zenodo, HuggingFace, IGVF, Arc — **nothing panel-covering exists that we do not
hold**) and measured **saturation at 1.84× DLD-1 at five sources, 1.77× at thirteen**. There is no
ladder to climb. Write the source programme down to **one arm**, not a sequence.

## 7.1 The arms that actually exist right now

| arm | owner | status | expected |
|---|---|---|---|
| the expressible source union | **229** | building; the literal 4-source arm **regresses** — a single union cannot express 225's construction. Fallback ladder: `dld1+igvf+myllia` (1.38×), then `dld1` alone (1.07×, +0.00228, and its October coverage **improves** 67→73) | one slot |
| the per-gene reliability gain | **228** | pricing 224's diagonal (**+8.9%** panel cosine, Δ`pds` **+0.0118**, gate cleared 3.1×) on all six members through PLOLO | one slot **if** it survives the scorer read-back |
| the operator's fitting rows | **226** | 50% subsample collapses (rank-deficient); the full-rank segment is the decisive one | unknown |
| **Hamilton (208's rider)** | **230** | **+0.002203, clears its floor at 1.16×** — *unless* rule 13's toll applies. 230 is testing step-vs-scaling | one slot **if** the toll scales |

## 7.2 The honest position on the budget

**We have ~90 pre-Oct-22 slots and four candidate arms, at most two of which are ready.** §6.2 said
we were under-spending a resource that expires worthless; that is still true, but **the binding
constraint has moved from willingness to inventory.** Slots are not the scarce thing any more —
**hypotheses are.**

**Say this plainly rather than manufacturing arms to fill slots.** A slot spent on a
one-variable arm we cannot interpret is worse than an unspent one, because it consumes the
build/verify capacity that a real arm would need.

## 7.3 So the programme is

1. **Spend both 2026-09-05 slots** on whichever two of §7.1 are ready and single-variable.
2. **After that, spend a slot only when an experiment produces an arm that clears its gate.** Do not
   schedule slots; let the science schedule them.
3. **Re-run 219's audit trick instead of inventing arms.** It found **three re-opened closures and
   was forbidden from measuring any of them** — one of those (Hamilton) is now a live +0.0022
   candidate found 30 days after it was closed. **214's and 200's re-openings are still
   unharvested.** Auditing our own closures has a better hit rate right now than new hypotheses.
4. **Reserve the last four slots before Oct 22** for the re-rehearsal artifact (`DEF-RUNBOOK` §5.-1)
   and whatever the transition breaks.
5. **The October asymmetry outranks marginal score.** `DEF-RUNBOOK` §5.-2: `kolf2` and `feng_ipsc`
   project to **exactly zero** on a disjoint panel; our shipped H1 union degrades **64→34** and
   **20→3**; four sources we hold are worth **zero today and ~30 targets each in October** and
   **cannot be priced until the panel changes**. **An arm that is panel-independent by construction
   — 224's gain is fitted on non-panel targets — is worth more than its headline number says.**
