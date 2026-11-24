# 239 — KNOWLEDGE AUDIT: the facts we measured and are not using

**Date:** 2026-09-04. **Fact audit, not a re-pricing.** Read-only over the repo. **No `vcc`, no
`cell-eval2`, no build, no submission, no production change, nothing re-measured, no slot
proposed.** Nothing written outside `experiments/239-knowledge-audit/`. Peak RSS **< 0.1 GB** —
this was `grep`, `sed` and reading; the one script written (`work/headlines.txt` extraction) is a
shell loop over `head`. All Python would have run under `uv run python`; none was needed.

**Scope read:** `AGENTS.md` (all 1,699 lines), `docs/LOG.md` (via targeted retrieval),
`docs/CURRENT_STATE.md`, `docs/BOARD-LEDGER.md`, `docs/DEF-RUNBOOK.md`, `docs/SLOT-PROGRAM.md`,
`docs/STRATEGY-2026-09-04.md`, `docs/NEXT-BUILD-SPEC.md`, the three prior audits
(**181**, **198**, **219**) in full, and the headline block of **all 216 `experiments/*/RESULTS.md`**
(`work/headlines.txt`, 6,448 lines, 216 blocks). Primary sources — `scripts/def_pipeline.py`,
`scripts/def_postbuild.py`, `scripts/build_union_whitener.py`,
`data/donors/whitener_union_pruned.json`, `experiments/202-callset-rider/build202.sh`,
`experiments/222-h1-union-build/build222.sh` — were opened for every finding that rests on one.

**Prior audits read so as not to repeat them.** 181 audited *claims closed on one day*. 198 audited
*closures against instruments known to be wrong*. 219 audited *closures against a corrected
pricing multiplier*. All three re-priced. **This one asks a different question: which measurements
in the record contradict, or are absent from, what we are doing right now.** Where a prior audit
already found something, I say so and do not re-report it as mine.

---

## Headline, in five lines

1. **🔴 The Oct-22 driver cannot build our current entry, and it silently overrides an attempt to
   fix it.** `scripts/def_pipeline.py:58` hardcodes `VCC_CD4_TARGET_F = "0.1"` and passes it as an
   env **override**, never sets `VCC_H1_UNION_G`, and diffs the result against `rc_callset`'s meta.
   Run as documented on the only prize-eligible entry, it drops **both** board-validated direction
   arms — **−0.0199671 overall** — and `def_postbuild.py` **passes**.
2. **🔴 `docs/STRATEGY-2026-09-04.md` §4b records a candidate as *"never been tried"* that
   experiment 091 measured at **+0.002409 (t = +19.85)** and that **production has shipped since
   a15**.** The shipped operator's own metadata says `"basis_prune": "essential"`. Same failure mode
   as the essential-gene case, in the document written today.
3. **🟠 The board's `n_real` has two board-derived measurements 2–3× apart and neither is in
   `AGENTS.md`.** 181 says median **90–110** (from board `reach`); 192 says **31–45** (from board
   `fid` coverage, corroborated by `jac`). `AGENTS.md` §2.1 and §8.1f both carry **89**, which is
   h1-test's. 198 §7 item 4 asked for this two audits ago.
4. **🟠 219's own #1-ranked recommendation was never harvested**, and `DEF-RUNBOOK` §3.1 still
   quotes the number 219 flagged as unsafe.
5. **The record is otherwise healthier than I expected.** Six leads I chased as candidate orphans
   turned out to have been closed by the right experiment (§4). **The two red findings are both
   in documents written in the last twelve hours, not in the deep record** — which is itself the
   finding: the failure mode is *summary drift at the head*, not lost measurements in the tail.

---

## 1. The ranked table

Ranked by **decision relevance** — *"would this change something we are doing or about to do?"* —
not by how interesting it is. Everything below row 8 is recorded for completeness and changes
nothing.

| # | cat | the fact | the record says | what we are doing | changes |
|---:|---|---|---|---|---|
| **1** | 🔴 CONTRA | The Oct-22 driver builds `rc_callset` at `f = 0.10` with no H1 union | `docs/BOARD-LEDGER.md:92-94` + `docs/CURRENT_STATE.md:14-17`: `cd4f034` **+0.0130557**, `h1union` **+0.0069114**, shipped, board 0.2118310 | `scripts/def_pipeline.py:58` `INCUMBENT_ENV = {"VCC_ALPHA":"0.982","VCC_CD4_TARGET_F":"0.1"}`; `:152` `--incumbent-meta` default `202-rc_callset...`; `docs/DEF-RUNBOOK.md:112,64,158` | **−0.0199671 overall on the only entry that counts**, silently. §2.1 |
| **2** | 🔴 CONTRA | "Fit the operator on non-essential rows only" is shipped, not untested | `experiments/091-basis-prune/RESULTS.md:14` **+0.002409, se 0.000121, t = +19.85**, composition not size; `data/donors/whitener_union_pruned.json` `"basis_prune":"essential"`, 4,393 rows dropped | `docs/STRATEGY-2026-09-04.md:176` *"Fitting the operator on **non-essential rows only** has never been tried"*, queued conditional on 238 | Kills a queued experiment; corrects §4b's premise. §2.2 |
| **3** | 🟠 DISAGREE | Board `n_real` median: **90–110** vs **31–45**, one-parameter family, two board members | `181/RESULTS.md:503` (board `reach` ⇒ f ≈ 1.05, median 90–110) vs `192/RESULTS.md:329-352` (board `fid` coverage ⇒ f ∈ [0.35,0.50], median 31–45, mean 270–383; `jac` agrees at f ≈ 0.35) | `AGENTS.md:74` *"median of 89"*; `AGENTS.md:829` *"89 of 17,806 scored genes = 0.50% of the axis"* — both h1-test numbers used as board numbers | Sets `fid`'s denominator and the October call-set decision; §8.1f's headline share is 2–3× too generous. §2.3 |
| **4** | 🟠 ORPHAN | 219's #1 recommendation is unharvested and the number it flagged is still quoted | `219/RESULTS.md:433` row 1, cost **LOW**, confidence **high**; `:404` *"200's 'six load-bearing, four inside their gate' ordering is not safe to quote until it is re-run with `VCC_MASS_FIX` unset"* | `docs/DEF-RUNBOOK.md:176` still quotes 200's PerturbMatch row (−0.001615) as the live evidence for an October ship decision; `SLOT-PROGRAM.md:206` *"214's and 200's re-openings are still unharvested"* — 214's was harvested (231), 200's was not | The October component ranking rests on ten `pds` legs measured on a renderer only a23 ever shipped (a23 scored −0.00073) |
| **5** | 🟡 CONTRA | `AGENTS.md` §8.2 advertises a reopening that was closed the next day | `experiments/203-reopen-125/RESULTS.md:1` *"Experiment 125 reopened, adjudicated, and **closed for the right reason**"*, composed gate **−0.001393** vs +0.0089, `nmae` −0.002468; `DEF-RUNBOOK.md:182` lists "the hierarchical estimator (203)" as closed | `AGENTS.md:1005-1009` *"Its measured +0.0122 top-20% sign precision is a direction gain that **was never adjudicated**"* | Sends the next agent at closed work, from the section that "reorders every priority" |
| **6** | 🟡 STALE | `--overdispersion 4.0`'s justification is a board point 200 showed does not apply | `200/RESULTS.md` headline: *"In a8's flat-0.25 context it is +0.00393 raw… In a15's tiered context the **same knob** is **−0.02610 raw** (t = −7.49). The board point that justified it was measured in a context we no longer ship"* | `DEF-RUNBOOK.md:160` evidence column: *"a8 board (`reach` peak 0.1890 raw)"*, action *"never"* | Documentation only: §8.1c caps the whole `reach` channel at **+0.002689 < 0.004502 noise**, so nothing is actionable. Fix the citation, not the flag |
| **7** | 🟡 CORRECT | 230 understates the toll's realised reach | `230/RESULTS.md` *"the charged row never contained a scored arm, so the toll's only realised effect on this project was the veto on 208's rider"* | `177/RESULTS.md:256`, `179/RESULTS.md:42,345`, `189/RESULTS.md:37-39,519` were **each** vetoed with `−0.00243` as the binding term | **No verdict reverses** — 177/179 die on a21's realised board deltas (181 §1), 189 at +0.00083 vs a +0.0019 floor. But "the toll never mattered" is not safe to reuse |
| **8** | 🟢 ORPHAN(null) | 198's own #1 — *"the single highest value-per-hour item in the audit"* — was never run, and is now dead | `198/RESULTS.md:156,576` re-score 150/151 on h1-test; expected `[−0.0004,+0.0013]`, carrying `jac` **+0.00245** *"where we rank 169th of 612"* | `BOARD-LEDGER.md:116` `jac` is now **world 1 of 646**; `AGENTS.md` §8.1c closes the ordering channel at +0.002689 | Nothing. **Record it as dead** so it is not re-opened as "highest value per hour" |
| 9 | 🟢 DISAGREE | `d reach / d log mse`: **+0.0127** vs **+0.03468** (2.7×) vs partial corr **+0.543** | `186` (quoted at `AGENTS.md:938`) vs `198/RESULTS.md:335` refit on the full 613-entry snapshot vs `BOARD-LEDGER.md:65` | Used to argue "amplitude explains almost nothing" for `reach` | Conclusion survives at either slope (198's residual sd 0.0324 ≈ 1σ). Recorded, not actionable |
| 10 | 🟢 DISAGREE | PLOLO→board `k_cos` **1.384** (`cd4f034`) vs **1.000** (`h1union`), 7.84× spread over four pairs | `AGENTS.md` §8.1c / `231/RESULTS.md:4` | **221 owns it** — its six-member cloud calibration is the instrument | Not mine; do not re-derive |
| 11 | 🟢 DISAGREE | Delivered-`pds` γ maximiser **0.50–0.55** (170) vs **0.80** (212) | `170` §2 pt 5 / `212` §1 | γ closed three ways, the binding one a **ratio** argument that is instrument-independent | Nothing |

---

## 2. The top three, in detail

### 2.1 🔴 The Oct-22 driver is pinned to the 2026-09-03 incumbent, and it overrides you

**This is the finding. It is a code defect on the one path that produces the prize-eligible entry,
and it is exactly the class of silent failure 199 found on the operator axis and 235 hardened.**

**What the record says we ship.** Two board-validated single-variable arms landed after the runbook
was written:

| entry | arm | the one variable | Δ overall | source |
|---|---|---|---:|---|
| `QSe3JlNUYZ2BzvlcsL14` | `cd4f034` | `VCC_CD4_TARGET_F` 0.10 → **0.34** | **+0.0130557** | `CURRENT_STATE.md:16` |
| `iUziqE1gewfcPKzYUDaB` | `h1union` | `VCC_H1_UNION_G` = **0.28** | **+0.0069114** | `CURRENT_STATE.md:17` |
| | | | **+0.0199671** | board 0.1918638 → **0.2118310**, rank 3 of 656 |

**What the Oct-22 path actually does.** `docs/DEF-RUNBOOK.md:90` gives one command,
`scripts/def_pipeline.py`, and §2's preamble is *"Follow it top to bottom. Do not improvise on the
day."* Inside that driver:

```python
# scripts/def_pipeline.py:56-60
# The incumbent's build flags, byte-for-byte from experiments/191-rc-hybrid-build/build191.sh
# (reused unchanged by 201 and 202). The env pair is exported for the effects run too.
INCUMBENT_ENV = {"VCC_ALPHA": "0.982", "VCC_CD4_TARGET_F": "0.1"}   # the ONLY env the incumbent sets
INCUMBENT_FLAGS = ["--min-cpm", "10", "--tail-spread", "1.0", "--keep-top", "8000",
                   "--keep-spread", "0.25", "--overdispersion", "4.0"]
```

The comment names its own provenance: **`build191.sh`**, the 2026-09-03 incumbent, two shipped arms
ago. And the env is applied as an **override**, not a default:

```python
# scripts/def_pipeline.py, Runner.run
proc = subprocess.run(full, cwd=ROOT, env={**os.environ, **(env or {})}, ...)
```

`{**os.environ, **env}` puts `INCUMBENT_ENV` **last**, so `VCC_CD4_TARGET_F=0.34 uv run python
scripts/def_pipeline.py …` is silently reset to `0.1` for both the effects stage (`:253`) and the
build (`:257`). `VCC_H1_UNION_G` is not in the dict at all, so it would pass through from the
environment — but only on top of `f = 0.10`, a combination that has never been built or scored
(`h1union` was measured on top of `f = 0.34`; `build222.sh` sets both).

**And the guard cannot catch it.** `def_pipeline.py:152` defaults `--incumbent-meta` to
`submissions/202-rc_callset.h5ad.meta.json`, and `DEF-RUNBOOK.md:113` step 10 passes exactly that to
`def_postbuild.py --allow-config-diff whitener_id`. The artifact the driver builds **is**
`rc_callset`'s config, so the "config differs in exactly `contexts`, `whitener_id`" check
**passes**. The runbook's own §4 known-answer row is `build202.sh rc_hybrid` — also pre-`cd4f034`.

**The three places the runbook still carries the stale config, verbatim:**

- `DEF-RUNBOOK.md:64` — *"| **incumbent package** `rc_callset` (board **0.1918638**, rank 11/643,
  entry `31oB9yWy8CSNBDEkY4iW`, 2026-09-03) | `submissions/202-rc_callset.h5ad` |"*
- `DEF-RUNBOOK.md:112` (step 9, the build) — *"the incumbent's flags byte-for-byte from
  `experiments/202-callset-rider/build202.sh` (`VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.1`, …)"*.
  `build202.sh` confirms: `export VCC_CD4_TARGET_F=0.1`, no `VCC_H1_UNION_G`.
- `DEF-RUNBOOK.md:158` (§3, *"the blind decision table — every shippable flag, pre-decided"*) —
  `VCC_CD4_TARGET_F` **SHIPS `0.1`**. There is **no row for `VCC_H1_UNION_G` at all.**

**The runbook contains its own pre-decided rule and it has already fired.** `DEF-RUNBOOK.md:175`:

> *"**`cd4f034` (entry armed for the 2026-09-04 00:00 UTC slot) decides this row** — if it lands
> positive, `f = 0.34` supersedes `0.1` for D/E/F."*

It landed at **+0.0130557**. The row was never updated.

**Why this is not merely a stale doc.** §5.-1 *does* record that the rehearsal is stale
(`DEF-RUNBOOK.md:287`, *"`VCC_H1_UNION_G` = 0.28 — a code path added by 222 … The rehearsal never
executed it"*). So the *rehearsal's* staleness is on the record. **The decision table's and the
driver's are not** — and those are what execute on the day, blind, with no board to read back.

**What would settle it.** Nothing to measure. It is a three-line correction plus a re-rehearsal
that `DEF-RUNBOOK` §5.-1 already schedules for **no later than 2026-10-15**:
(a) `INCUMBENT_ENV` → `{"VCC_ALPHA": "0.982", "VCC_CD4_TARGET_F": "0.34"}` and add
`VCC_H1_UNION_G: "0.28"` (or make the driver take them as arguments so a future arm cannot be
silently overridden again); (b) `--incumbent-meta` default → the `h1union` meta;
(c) `DEF-RUNBOOK` §1/§2 step 9/§3 rows and the §4 known-answer row → `build222.sh`.
**I did not make these edits — they are outside this experiment's directory and the mandate forbids
it.** The owning agent should carry `build222.sh`'s known-answer row (`VCC_H1_UNION_G` unset ⇒
sha256 `7036e89f…` on 12,648,807,800 bytes) as the proof the change is inert on A/B/C.

---

### 2.2 🔴 §4b's candidate is shipped, not untested — the essential-gene failure, one more time

`docs/STRATEGY-2026-09-04.md` §4b, written today, records a candidate generated by the
essential-gene finding (`:164` *"One candidate the essential-gene finding generates — recorded, not
launched"*). Its premise, quoted in full:

> `:170` **"But the whitening operator is fitted on donor effect vectors that INCLUDE
> essential-target rows"** — and essential knockdowns are exactly the rows dominated by the shared
> death/stress programme that §4.2 says whitening exists to remove.
>
> `:174-176` **"This is untested.** 193 closed the **regulariser**, 212 the **exponent**, 226 the
> **fitting sample size** — and 226 states it held *"same three donors, same prune rule,
> **composition held fixed**"*. **Composition is not sample size.** Fitting the operator on
> **non-essential rows only** has never been tried."

**All three factual claims are wrong, and the primary source is one `cat` away.**

**(a) It was tried.** `experiments/091-basis-prune/RESULTS.md:1` — *"091 — Pruning common-essential
rows from the whitening basis"*, verdict *"**POSITIVE, and it is a composition effect, not a size
effect**"*:

| arm | mean `avg6` |
|---|---:|
| **`BASIS_PRUNE=essential`** | **0.122493** |
| `BASIS_PRUNE=random` | 0.120084 |

> `091:14` **"delta `+0.002409`, se `0.000121`, t `= +19.85`. Per-seed ranges do not overlap."**
> `091:17` *"Both arms remove the **same number of rows** (3,810 …). The only difference is
> *which* rows. So basis size is held constant and the entire effect is composition."*

That is §4b's experiment, run, with §4b's own confound (size vs composition) controlled by design.

**(b) It ships.** `scripts/build_union_whitener.py:254-256` — `--basis-prune {none,essential}`,
*"'essential' drops basis rows whose perturbation is in the common-essential list"*. And the
production operator's own metadata, `data/donors/whitener_union_pruned.json`:

```json
{"n_cols": 17942, "basis": [29394, 17942], "basis_prune": "essential",
 "essential_list_size": 2021, "rows_dropped_total": 4393, "dropped_in_panel": [],
 "sha256": "e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123"}
```

`e2f65ecf…` is the sha `DEF-RUNBOOK.md:60` pins as *"A/B/C production operator (rule 11: this sha,
not the filename)"*, and `DEF-RUNBOOK.md:108` step 5 refits it on Oct 22 with
`--basis-prune essential`. **4,393 essential-target rows are already out of the operator we ship,
and have been since a15.**

**(c) 226's "same prune rule" is the essential prune.** `226:139` reads *"same three donors, same
prune rule, composition held fixed"* — 226 held it fixed **because it is production's rule**, which
is the opposite of the inference §4b draws from the sentence.

**Corroborating measurement nobody has joined up either:** `experiments/042-auxiliary-whitening/RESULTS.md:18`
added the essential libraries to the covariance from the other direction —
`+ K562-essential + RPE1` → **−0.00055**, and `+ all auxiliary, equal source mass` → **−0.00498**.
Essential rows in the operator were measured as harmful in **August**, twice, from both directions.

**What is left, honestly.** 091 §"What it does not yet establish" is explicit: the comparison is
essential-prune vs **random**-prune at equal row count, *not* against the **unpruned** basis on that
harness. So a narrow question survives: *essential-pruned vs unpruned, on today's instrument.* It is
not the question §4b asks, and three measurements bracket it against paying — 091's +0.0024 (the
gain is real and already banked), 090's t = −13.9 for adding unrepresentative rows, and 226's cliff
(dropping rows *reduces* `n`, and production sits at `n/p = 1.638` against a collapse below 1.4).

**What would settle it:** nothing needs measuring to correct §4b. If the narrow variant is ever
wanted, it is 226's instrument with `--basis-prune none`, reproducing `e2f65ecf…` on the pruned fit
as the known-answer row — but the prior is against it and §4b's stated conditional ("if 238 says the
essential distinction is material") should be struck, because the answer does not depend on 238.

**The process point, and it is the reason this audit exists.** `AGENTS.md` §8.1g already wrote the
rule after the last instance: *"Before building on a synthetic panel — or resurrecting any source
previously called dead — **grep the log for it first.**"* This is the same error one section later,
on a different object, in a document written after that rule was landed. **The rule needs to cover
the operator and the builder, not just sources and panels:** before recording any candidate as
"never been tried", grep `experiments/*/RESULTS.md` for the mechanism **and** read the shipped
artifact's own metadata. The metadata is authoritative and it is 569 bytes.

---

### 2.3 🟠 The board's `n_real` has two board measurements 2–3× apart; `AGENTS.md` carries a third

**Same model, same data source, two members, two answers, no adjudication.** Both estimates fit the
*same* one-parameter family — scale h1-test's per-perturbation `n_real` distribution by `f` — to
board rows:

| route | member | fitted `f` | median `n_real` | mean | source |
|---|---|---:|---:|---:|---|
| board `reach` = 0.1742, proportional power matching | `reach` | **≈ 1.05** | **90–110** | ~765 | `181/RESULTS.md:490-503` |
| board `fid(n_pred)` over 5 emission-only arms, coverage term | `fid` | **0.35–0.50** | **31–45** | 270–383 | `192/RESULTS.md:329-352` |
| a15 board `jac` 0.030080 at `n_pred` 4,474 | `jac` | ≈ 0.35 | ~31 | ~270 | `192/RESULTS.md:344-350` |

192's two routes agree with each other to 1.4× in `f` and are explicit about the disagreement:

> `192:351-353` *"Band **f ∈ [0.35, 0.50]**: median `n_real` **31–45**, mean **270–383**. **Below
> audit 181's board-implied median of 90–110** and well below h1-test's 89/765."*

181 is equally explicit about its own instrument's limits (`181:509-513`, *"`n_real` is a
distribution, not the scalar 240"*), and 208 read the conflict and used 192's number:

> `208/RESULTS.md:175-176` *"h1-test's `n_real` has **median 89 but mean 765** … **192 §4.2 puts the
> board's `n_real` at median 31–45 / mean 270–383.** Scaling h1-test's distribution by a scalar
> preserves a tail the board does not have."*

**What `AGENTS.md` carries is neither of them.** `AGENTS.md:74` — *"188 independently reconstructed
the reference pool at a **median of 89** genes"* — is h1-test's own median, used as the challenge's.
And it propagates into the sign channel's closing argument:

> `AGENTS.md:829` (§8.1f) — *"`q_w` lives on the **reference's** DE pool — **89 of 17,806 scored
> genes = 0.50% of the axis.** A key cannot know that pool, so **99.5%** of a pool-blind arm's flips
> land where they can earn no sign credit."*

Under 192's band that reads **31–45 of 17,806 = 0.17–0.25%**, i.e. **99.75–99.83%** wasted.
**The direction of the error is favourable to the conclusion**, so §8.1f's verdict does not reverse
— it strengthens. But the *number in the master document is an h1-test number presented as a board
number*, and **198 §7 item 4 asked for exactly this correction two audits ago**:

> `198/RESULTS.md:584-585` *"**update AGENTS.md §2.1's `n_real` median from 89–110 to 31–45**, which
> 192 measured twice from two board members."*

It was not done, and 219 (the next audit) classified the ordering family as 220's territory and did
not revisit it.

**Which is right.** Neither is established; they are one-parameter fits with honest error bars
(192: *"a one-parameter fit to five points inside an assumed shape family"*, RMS 0.00126 vs a
no-coverage null of 0.00432). But they are not equally weighted for the same question: **`fid`'s
coverage term is a direct function of the `n_real` distribution and of nothing else** (192 §3.2),
while `reach` is dominated by the small-`n_real` stratum and by ordering, which is exactly the
channel that has moved twice since 181 measured it (`rc_hybrid`, `rc_callset`, `h1union` took raw
`reach` 0.1742 → 0.276545). **181's `f ≈ 1.05` was solved from a board `reach` of 0.1742 that no
longer exists.** That alone is a reason to prefer 192's, and it is a re-derivation nobody has done.

**What changes if it is 31–45 rather than 90–110.**

- `fid = k / max(n_pred, n_real)`: at `n_pred` 2,599 (`rc_callset`) the yield cap is saturated under
  both, so **today's call-set policy does not move**. `DEF-RUNBOOK.md:157` says *"ship the board
  winner"* — unaffected.
- `AGENTS.md` §8.1f's 0.50% becomes 0.17–0.25%: the sign channel closes **harder**.
- It is the input to the one October decision that cannot be measured on D/E/F
  (`CURRENT_STATE.md` transition item 5, *"No D/E/F `n_real` can be measured"*), so it should be
  right before Oct 22 rather than after.

**What would settle it — and it is cheap, and it is not mine to run.** Re-run 181 §14's
`nreal181.py` proportional power-matching against the **current** board `reach` (0.276545, `h1union`)
instead of a15's 0.1742, and against the current `n_pred` 2,599 rather than 4,474. If 181's route
then lands inside 192's `f ∈ [0.35, 0.50]`, the disagreement was an artefact of a stale board row and
the two members agree. If it still reads `f ≈ 1`, the two members genuinely disagree and the
`AGENTS.md` line should carry both with their provenance. **Either way `AGENTS.md:74` and `:829`
should say which surface the number comes from.** No new data, no build, one script that already
exists.

---

## 3. What I looked for and did not find

Reported because a negative from a fact audit is a result, and because these are the leads a future
agent would otherwise re-chase. **Six candidate orphans came back closed by the right experiment:**

| lead | why it looked orphaned | why it is not |
|---|---|---|
| 189's `theta` as a **flip** key | §8.1e explicitly leaves it open: *"It remains testable as a flip key"* | **234 tested it**: `theta` Δ`q_w` **+0.00037**, `key189` +0.00022, against a +0.01537 bar (`234:157,159,171`) |
| 189's `theta` in the **call-set** channel | 224 moved it to direction; nobody names the call set | **189 tested it**: *"the gain is entirely reproduced by the shuffled-theta control (z = −0.04)"* (`189` §4) |
| 181's item 4 — *"reopen `pds` with one lever that is NOT an un-whitening move"* | flagged REOPEN, load-bearing on *"the entire October research direction"* | **183 ran it** (`183:4-5` cites it by name): unique-component boost, PLOLO +0.00104 raw `pds`, 3× short |
| 196's per-context **tilt** (198 §6.2: *"196 should price the tilt, not the shortfall"*) | 198 said 196's abort condition had fired | **196 repointed**: its headline is the `floor`→largest-remainder fix that removes **97%** of the tilt; it became 208's rider A and 230's +0.002203 |
| 192's `√cpm` / detection-rate prior, *"computable from the 18,400 NTC cells we get on day one"* | the one panel-independent asset in the stack | **Already shipped and pinned for October**: `DEF-RUNBOOK.md:157`, `top 400 by |H| ∪ top 2,200 by control CPM` |
| 227's un-estimable per-row offset (**+0.0542 raw at oracle**) keyed on the reference's gate, vs 234's *"the missing ingredient is knowledge of the reference's DE pool"* | two channels independently naming the same missing object; a partial estimator of that object exists (192, fold-enrichment 2.37) | **234 already ran the oracle-assisted version**: hand the arm the reference's DE pool and it prices **+0.00057** against a +0.0019 gate — *"3.3× short with the answer handed to it"*. A partial estimator cannot beat the full one |

**And the shape of what I did find is worth stating.** The two red findings are in
`docs/STRATEGY-2026-09-04.md` (written today) and `scripts/def_pipeline.py` + `docs/DEF-RUNBOOK.md`
(pinned 2026-09-03, superseded twice on 2026-09-04). **Neither is a lost measurement in the deep
record — both are heads that stopped tracking a tail that moved.** 198's conclusion holds and this
audit re-confirms it from a different angle: *"the invalidating findings were, in almost every case,
applied by the next experiment. The map is self-healing."* What is **not** self-healing is the
**summary layer** — `AGENTS.md` §8.2, `DEF-RUNBOOK` §1/§2/§3, `STRATEGY` §4b — which is written by
the coordinator, is not owned by any experiment, and has no known-answer row. Three of the eleven
rows above are that one defect.

**The cheapest general fix, offered as an observation rather than a proposal:** the build path
already has the pattern that would have caught row 1 — `build222.sh` carries a known-answer sha and
refuses to proceed if the incumbent moved. `def_pipeline.py` has no equivalent assertion that its
`INCUMBENT_ENV` still describes the current board leader's config, and `DEF-RUNBOOK` §4's
known-answer rows check the *operator* and the *A/B/C build*, never *which arm we are rebuilding*.

---

## 4. Integrity

No `vcc`, no `cell-eval2`, no `cell-eval`, no build, no `vcc prep`, no submission, no write to
`submissions/`, no production code touched, **no file outside
`experiments/239-knowledge-audit/` created or modified** (this file and `work/headlines.txt`).
Nothing re-measured; every number above is transcribed from a primary source with its path and
line. No slot proposed. Two sibling agents (221, 238) were live throughout and neither's territory
was entered — row 10 is explicitly handed to 221, and §2.2's conditional is explicitly *not* handed
to 238 because it does not depend on 238's result.

**Known-answer rows carried** (rule 10 — a fact audit's known answer is a fact it can re-derive
from a primary artifact, not a re-scored number):

| row | expected | measured here | |
|---|---|---|---|
| production operator sha in its own metadata == `DEF-RUNBOOK.md:60`'s pin | `e2f65ecf69b7…1db123` | identical | ✓ |
| production operator's `basis_prune` field | `essential` (091 shipped) | `"basis_prune": "essential"`, 4,393 rows dropped | ✓ |
| `build202.sh` sets the CD4 dose the runbook quotes | `VCC_CD4_TARGET_F=0.1` | `export VCC_CD4_TARGET_F=0.1`, no `VCC_H1_UNION_G` | ✓ |
| `build222.sh` sets the shipped dose | `0.34` + `VCC_H1_UNION_G` | `export VCC_CD4_TARGET_F=0.34`; `VCC_H1_UNION_G=0.28 build222.sh` | ✓ |
| `cd4f034` + `h1union` sum to the board move | 0.1918638 → 0.2118310 | +0.0130557 + 0.0069114 = **+0.0199671** | ✓ |
| `experiments/*/RESULTS.md` on disk, all enumerated | ≥ 196 (219's count, earlier the same day) | **216** blocks in `work/headlines.txt`; 217 files after this one was written | ✓ |
| 150/151's six-member averages reconstruct (198's KA row, re-checked) | −0.005045 / −0.004344 | as recorded in both files | ✓ |
