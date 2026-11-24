# Current state

> **Historical notebook, not the live state.** As of September 9, 2026, use
> `docs/READ-FIRST.md` for current artifacts, scores, permissions and the reading
> map. Dated entries below are retained as evidence, including superseded claims.

> # 📌 This file is 300 KB of dated narrative. Read `docs/READ-FIRST.md` first (2026-09-04).
> Everything below is kept as history. Numbers here are not authoritative; the archive and READ-FIRST are.


> # STATE AS OF 2026-09-05 01:20 UTC — both 09-05 slots spent; **incumbent `cd4svd5` 0.2160944, rank 3 of 688**, gap to first +0.00845
>
> `ARM6` (237 gate test) 0.2081907, Δ −0.0036403: the `nmae` identity VALIDATED at scale (leg −0.00716 vs −0.004783).
> `cd4svd5` (`VCC_CD4_DEFLATE_K=5`) 0.2160944, Δ +0.0042634: `pds` +0.00545 overall (k_pds 1.22), `mse` −0.00080
> (re-shape cosine factor 3.76 on `cos_agg`); October ships the PINNED basis `V5_abc` (runbook trigger fired).
> Full records `data/board-entries/`; readings `experiments/237-gate-validation/BOARD-RESULT.md`,
> `experiments/249-cd4-deflate/BOARD-RESULT.md`; `docs/READ-FIRST.md` v2.2 is the live state.

> # STATE AS OF 2026-09-04 04:20 UTC — four slots spent in two UTC days, **+0.0377445 scored**
>
> **Best and displayed entry: `0.2118310`, `h1union`.** **RANK 3 of 656.** Was **0.1740865 /
> rank 24** at 2026-09-03 00:30 UTC. **Leader: `cqawesome` `0.22431`. Gap: `+0.01248`** — was
> +0.01939, so **36% of it closed in one slot.** Implied cosine **0.09947**, leader 0.13356,
> ratio **1.343×** (was 1.481×).
>
> The leaderboard displays a team's **latest** entry, not its best. **Never spend a slot on rank.**
>
> | # | arm | overall | Δ | the one variable |
> |---|---|---:|---:|---|
> | 1 | `rc_hybrid` | 0.1834190 | **+0.0093325** | the emitted p-ORDER, at bit-identical pseudobulk |
> | 2 | `rc_callset` | 0.1918638 | **+0.0084448** | the CALL SET (`n_pred` 4,480 → 2,599) |
> | 3 | **`cd4f034`** | **0.2049196** | **+0.0130557** | `VCC_CD4_TARGET_F` 0.10 → 0.34 — **first arm ever to move the DIRECTION channel** (`cos` 0.07228 → **0.09018**) |
> | 4 | **`h1union`** | **0.2118310** | **+0.0069114** | `VCC_H1_UNION_G` = 0.28 — a fourth donor estimate via the H1 union. **ALL SIX MEMBERS UP**, including `nmae` (+0.00588), which the exact-identity instrument predicted would fall |
>
> **World rank per member:** `jac` **1st**, `mse` **10th**, `reach` **10th**, `pds` 7th,
> `nmae` 59th, `fid` 129th. **Every member improved with `cd4f034`.**
>
> ### What is actually reachable — the board's own ceiling (2026-09-04 04:33Z, 657 entries)
>
> | member | ours | world best | holder | headroom |
> |---|---:|---:|---|---:|
> | `pds` | +0.74537 | +0.82637 | cqawesome | +0.08099 |
> | `mse` | +0.17957 | +0.27101 | MT Zeng | +0.09144 |
> | `nmae` | +0.11074 | +0.19847 | CellWM | +0.08772 |
> | `fid` | −0.01648 | +0.02273 | The_Way | +0.03920 |
> | `reach` | +0.22294 | +0.23907 | luka | +0.01613 |
> | **`jac`** | **+0.02884** | **+0.02884** | **us** | **0** |
>
> - **ours 0.21183** · **actual #1 0.22431** · **best-of-board composite 0.26441**
> - *"our everything + Aginglab's `nmae` alone"* = **0.22404** — arithmetically within 0.0003 of
>   first, **but 227 showed this is NOT reachable at fixed cosine** (it costs +5–12% relative cosine
>   against a ±0.2% gate). Their `nmae` is a *consequence* of better direction, not a substitute for
>   it. **Treat the whole table as headroom, not as a menu of independently purchasable members.**
> - *"cqawesome's direction + Aginglab's `nmae` + our ordering"* = **0.24993**.
>
> **So: ~0.25 is the realistic stretch and 0.30 is not reachable this season** — it is above the composite of *every team's best member taken
> simultaneously* (0.26441), so no combination of publicly demonstrated capability reaches it.
> Ahmet's stated target was "0.2+ and ideally 0.3"; **0.2 is cleared, 0.3 should be retired as a
> target and replaced with 1st place** — which is **+0.01248 away and, after 227, is a
> transfer-cosine problem end to end.**
>
> ### 🔴 RETRACTED 05:10 — "two teams, two reasons" was MY ERROR. There is ONE axis: cosine.
> At 04:33 I claimed `Aginglab.com` (rank 2) beats us on a **separate** `nmae` axis, citing their
> lower *implied* cosine and `corr(implied cosine, scaled nmae) = −0.394`. **227 falsified it and I
> verified the refutation myself.** **15 of those 55 teams have `score_mse` exactly 0.00000** at raw
> `mse` **1.935–4.00** — far past the clamp — so 197's inversion returns cosines of −0.44 to −1.43
> that were never measured. **Without them the correlation is +0.022.** `AGENTS.md` already said the
> inversion is *"vacuous above raw `mse` 1.02361"*; I used it anyway.
> **VALIDITY GUARD: only invert rows with `score_mse` > 0, and say that you did.**
> Aginglab is reconciled — a **2.4% smaller jackknife claim** (their `jac` −0.033 / `fid` +0.020 pair
> indicates it) makes their row what **+12.3% better transfer cosine** prints. **They are ahead on
> our axis; the instrument could not see it.** `nmae` is closed: reaching their credit costs
> **+5.0%/+12.3% relative cosine** against a ±0.2% gate, emission is closed *structurally* and the
> call set *by definition* (which is why `rc_callset`'s +0.0084 cost `nmae` **exactly zero**).
> **187 was right; my amendment of it was the false step.**
>
> ### The gap to `cqawesome` is cosine and almost nothing else
> Behind on all three direction members (`pds` +0.0886, `mse` +0.0947, `nmae` +0.0163), **ahead on
> all three ordering members** (`fid` −0.0396, `reach` −0.0084, `jac` −0.0351). 197's inversion:
> implied cosine **ours 0.09018 vs cqawesome 0.13356 (1.481×)**, worth **+0.01583** of the
> **+0.01939** gap — **cosine is 82% of the distance to first.** `AGENTS.md` **§4.1c**.
>
> **MT Zeng holds the board's best direction (cosine 0.14138, above cqawesome's) and ranks 8th**,
> because their ordering members throw it away. **A point of cosine is worth more to us than to
> anyone on the board. Do not trade ordering for direction.**
>
> ### 🔴 Do not size an arm through a scalar PLOLO→board factor — it is not identified
> Two arms in the same family, matched `cos_agg` functional: **`cd4f034` 1.384** (216 §4.4),
> **`h1union` 1.000** (222, 3 cancer folds — predicted +10.3065%, board delivered **+10.3045%**).
> **221 has both points and owns the reconciliation.** Until it lands, quote arms in PLOLO units
> with the **fold set and functional named** (`cos_agg` vs macro, which folds, which gate). Both of
> the last two mis-predictions — 216's 3.0× under-call and my own under-call on `h1union` — were
> functional-mismatch errors, not modelling errors.
>
> ### Slots
> **2026-09-04: both spent — `cd4f034` (0.2049196) and `h1union` (0.2118310).**
> `cd4f080` was built and **cancelled unspent** — 216 priced it at −0.0064 before the build
> finished. See `experiments/217-cd4-f080/CANCELLED.md`. Next slots: **2026-09-05 00:00 UTC.**
> ARMED for 2026-09-05 00:00Z (lead, 2026-09-04 11:37Z): `submissions/229-stack6-g035.vcc` under
> `experiments/237-gate-validation/PREREGISTRATION.md` — a gate-validation slot (229 measured it −0.00219 on the
> exact `nmae` identity, the rate card says +0.00268; the board adjudicates the instrument). DLD-1 alone
> (+0.00228, 225) is the alternative if the gate is validated.
>
> ### The three things a new agent most needs to know
> 1. **`AGENTS.md` §8.4 (as corrected by 219): a direction arm priced `mse`-only is understated
>    ≈2.51×.** Four members — `nmae`, `fid`, `reach`, `jac` — move with sign accuracy and must
>    **never** be carried at zero for a direction arm. The `+31%` single-arm cosine bar in §8.2 was
>    itself `mse`-only and is **2.0–2.3× too high**; corrected to **+11.6–13.6%** (amplitude floor)
>    and **+2.5–2.9%** (shape-only).
> 2. **`solve_cd4_fraction` sizes a residual by the median row norm over ALL 300 panel rows.**
>    Below 50% coverage that median is **exactly 0** → the source contributes **nothing**, silently
>    (218). Between 50% and ~90% it silently **over-doses**. **223 resolved this: it is a
>    MEASUREMENT bug, not a shipping bug** — the builder raises; the silent `return 0` is in the
>    harness mirror `experiments/162-retrodicting-validator/plolo.py:363-369`. **0 of 21 experiments
>    were void.** 222 measured the over-dose factor: **1.012–1.017× today, ~1.24× at the pessimistic
>    D/E/F draw** (`DEF-RUNBOOK.md` §5.0).
> 3. **The direction channel is not a wall.** 215 decomposed 707 board arm-pairs: the leaders'
>    direction arrived **inside the readable window, in many small steps** (`cqawesome` in 5,
>    largest only 28%). Our own history is one big step (whitening) and then nothing.
>
> ### Closed, do not re-open without new evidence
> donor lines (166, 205) · cross-line adaptation (184, on 50 Tahoe lines) · learned transforms (195)
> · the operator's eigenvalue map (193) · **its exponent (212** — a *ratio* argument, therefore
> instrument-independent**)** · donor re-estimation (194) · weighting (204) · donor–donor
> correlation (209) · the hierarchical estimator (203) · **global amplitude (211** — `claim` at
> 100.0000% of its cap**)** · **the transport law (213** — max +1.6% against the bar; a gain cannot
> reorder**)** · **the CD4 dose (216** — `f = 0.34` is the argmax over a 31-point grid and all three
> forms**)** · **the public-text route to the leaders' method (2026-09-04** — of 2,306 scored
> entries, 764 carry a description and **not one above implied cosine 0.09 does**; the whole top of
> the direction channel keeps it empty**)** · **the CD4/atlas "better data" reading (216** — 12 byte
> ranges match, `LastModified` 2026-05-28, **zero objects moved** on that bucket during the leaders'
> doubling window**)**.
>
> ### The two live hypotheses, both running
> 215's many-small-steps shape has exactly two readings left, and 216 and 218 killed the others:
> **224** — iterative retraining of a learned transfer model, gated at **+2.9% relative cosine on
> held LINES** with the slope in `k` as the decision variable; **225** — public sources added one at
> a time, an exhaustive census plus pricing, with saturation as the decision variable.

## 2026-09-04 10:47 UTC — 240: the Oct-22 driver is FIXED, and the old one is on the record building a −0.0200 submission that PASSED its own check

**The single most valuable fix today.** 239's audit found it; 240 fixed it; both are verified.

```
The driver's defaults, on the A/B/C bundle, with NOTHING exported, rebuild the shipped
`h1union` package BYTE-FOR-BYTE:
    .h5ad  8df4dadd5a63b76c7ca3e3d57a511e7e52000b79d016d140897c55dd24c3ecf6  12,648,807,800 B
    .vcc   f13ce40bca908af7cb85b8e84887af473aae4588410c238529139b7d1f245b1e   2,824,785,920 B
and `def_postbuild.py` passes 20/20 against the CURRENT incumbent.
```

**And the demonstration that matters — row 3, run deliberately and quoted verbatim:**

> **Before 240 the same command built `rc_callset`'s configuration — `f = 0.10`, no H1 union,
> −0.0199671 overall — and `def_postbuild.py` printed 14/14 PASS on it.**

**That is the failure, executed, with its safety check reporting success.** On the only
prize-eligible entry it would have taken board **0.2118310 → 0.1918638** — **larger than our entire
gap to first place (+0.01248)** — silently.

**All 5 known-answer rows pass. No shipped number moved; `make_submission_pm_cd4.py` is untouched
at 235's state.** Three files changed: the driver, the post-build check, and `DEF-RUNBOOK`.

### Why it happened, recorded so the class is closed
`INCUMBENT_ENV` was pinned by its own comment to **`build191.sh` — the 2026-09-03 incumbent, two
shipped arms ago** — and applied as `{**os.environ, **INCUMBENT_ENV}`, i.e. as an **override**, so
even exporting the right value was silently reset. `VCC_H1_UNION_G` was absent entirely. The
post-build check diffed against **`rc_callset`'s** meta — the incumbent from *before* both arms — so
it confirmed the stale configuration rather than catching it.

**This is the third silent-failure class found and closed in one day**, after 235's operator-axis
mismatch (**−0.021 on `pds` alone**) and 223's silent-zero sizing. **All three shared one shape: a
check that passes because it is checking the wrong reference.**

## 2026-09-04 09:18 UTC — 235: the Oct-22 silent failure is HARDENED, and it independently CONFIRMS the essential-gene problem

**All 5 known-answer rows pass. The diff is 181 insertions and ZERO deletions — no existing line of
either file was modified or removed, and no shipped byte moved.** Both files committed with hashes.

**ITEM A — the operator's panel mismatch now RAISES.** Observed firing **four ways**, in ~1 second —
and **the pre-235 builder is shown accepting the same wrong panel with no error at all**, which is
the demonstration that matters: the guard is proved against the failure it exists to catch.
A/B/C `axis_sha256` = **`edf6ba886286d7837fa6adaef152372cee481331a81987789151862d5c0a400f`**.
Previously the only check was `n_cols` = **17,942 on every panel**, and 233 priced skipping the
refit at **−0.021 overall on `pds` alone** (212: −0.0429).

**ITEM B — the carrying-row rule covers every source the registry can select.** **37 of 37
(source, panel) cells size ABOVE ZERO across three panels, and the shipped all-rows rule returns
EXACTLY 0.0 in all 37.** A selected source that carries nothing now **raises** instead of
contributing a silent zero.

### 🔴 And it confirms 8.1g independently, mid-run
235 picked up `AGENTS.md` §8.1g while running and built a **third panel — disjoint AND
essential-free**. On it:

> **`k562_essential`: 30 targets on 199's bundle, 0 here.**

**That is direct confirmation that the October asset evaporates on an essential-free panel**,
arrived at by a different agent on a different instrument from the one I asked to test it. **238 is
still the full measurement** — all four sources, the whole §5.-2 ranking — but the direction is no
longer in doubt.

**Net: the transition is hardened and the October asset is provisionally worth ≈0.** Those are
opposite-signed, and both are better known now than on Oct 22.

## 2026-09-04 07:53 UTC — 233 and 234: October retention MEASURED, and the sign channel CLOSED

### 233 — we retain 65% of our score on a randomly drawn disjoint panel, and it is a LOWER BOUND
Reusing `DEF-RUNBOOK` §5's two synthetic bundles verbatim, which hold cell lines, gene axis and
control cells fixed and vary **only which 300 genes are the panel**: **10 operator re-fits, 3 CD4
re-extractions, 81 six-member rows.**

| channel | retained |
|---|---|
| aggregate transfer cosine | **93%** |
| retrieval (`pds`) | **64%** |
| ordering | **60%** |
| **projected overall** | **0.2118 → 0.137 (65%)**; **0.119 (56%)** at the pessimistic coverage draw |

**The loss is carried by the TARGETS, not by coverage** — coverage explains **78% of the cosine
loss but only 18% of the `pds` loss**. **And that is why it is a lower bound, not a forecast:** the
property carrying the loss is the targets' cross-line transferability, which is **exactly what Arc's
selection protocol optimises and a random draw does not** — measured model-free at **1.3–2.0×**.

**Its verdict on allocation is "neither":** keep the direction programme, and spend the marginal day
on three October items. **229 had already done one** (the four October-only sources are registered
with paths and shas). **235 is doing the other two**, and their prices are now measured:
- **the four October-only sources are LIVE at +5.5% to +17.7% relative cosine each**, dose-monotone
  with a clean shuffle control — against **exactly zero today**;
- **their all-rows median is `0.0000` on both panels**, so under the old sizing rule they deliver
  **identically zero** — the carrying-row rule must cover them;
- **skipping the operator refit costs −0.021 overall on `pds` alone** (independently confirming
  212's −0.0429), and the only guard that exists is `n_cols` = **17,942 on every panel**.

### 234 — the sign channel is CLOSED, and my §8.1e target was wrong by 3.2×
**All five known-answer rows pass**; 232's exchange rate reproduces to **<5e-6** through
independently written pricing code. Then two corrections, **both raising the bar** (now §8.1f):
1. 232's rate is a **two-sided fit of an asymmetric response** — gaining pays **+0.002571**/+0.01,
   losing costs **+0.004827**, a **1.88×** asymmetry → the gate is **+0.00739**, not +0.00476.
2. **232's oracle flips only wrong signs; a real key flips both ways** — realising **0.44×** of the
   oracle rate → a real key's bar is **Δ`q_w` = +0.01537**. **Shortfall 5.7×.**

**The mechanism, in one number: `q_w` lives on the reference's DE pool — 89 of 17,806 genes, 0.50%
of the axis. A key cannot know that pool, so 99.5% of a pool-blind arm's flips earn no sign credit
and only cost cosine.** One arm flipped **9,953 coordinates** for Δ`q_w` of **exactly +0.00000** and
**−1.34% cosine**. **Hand it the pool** and it meets 232's published bar and still prices **3.3×
short**. **Every shippable flip arm is NEGATIVE on all six members at every dose.**

**§8.1e's +0.17325 headroom is real and remains unreachable** — the missing ingredient is knowledge
of a pool that is unknowable by construction.

## 2026-09-04 07:45 UTC — COMPETITIVE READ: cqawesome is PLATEAUING and we are ACCELERATING

Implied transfer cosine over each team's last four scored entries (validity guard applied — only
rows with `score_mse` > 0 are inverted):

| team | entries | best | last submission | cosine trend (last 4) | Δ over the 4 |
|---|---:|---:|---|---|---:|
| **cqawesome** | 12 | 0.22431 | 6.1h ago | +0.126 → +0.127 → +0.132 → **+0.134** | **+0.008** |
| Aginglab.com | 5 | 0.21288 | **26.6h ago** | +0.031 → +0.095 → +0.095 → +0.099 | +0.068 |
| Hongjian Zhou | 8 | 0.21064 | 12.1h ago | +0.098 → +0.127 → +0.130 → +0.130 | +0.032 |
| jw_13288 | 7 | 0.20949 | 5.4h ago | +0.080 → +0.058 → +0.131 → +0.132 | +0.052 |
| MT Zeng | 5 | 0.20744 | 5.4h ago | +0.093 → +0.069 → +0.127 → **+0.141** | +0.048 |
| **us** | 13 | **0.21183** | 3.9h ago | +0.072 → +0.072 → +0.090 → **+0.099** | **+0.027** |

**Three reads:**

1. **`cqawesome` is decelerating.** Their last three steps are **+0.001, +0.005, +0.002** — they gained
   **+0.008** cosine over four entries while **we gained +0.027**. If that holds, the gap closes on
   its own in two or three successful arms.
2. **The field's cosine ceiling looks like ~0.13–0.14.** Nobody has exceeded **MT Zeng's 0.141**, and
   the three Prima Mente accounts (206: one operator) all sit at **0.127–0.141**. **We are at 0.099.**
   Matching the leader needs **≈ +0.035 cosine = +0.0128 overall** at §8.2's rate — which is exactly
   the **+0.01248** gap. The arithmetic closes.
3. **`Aginglab.com` has not submitted in 26.6 hours.** The rank-2 team has stalled.

**Why this matters for allocation.** It argues against panic and for patience: **our problem is that
we ran out of arms today, not that we are being outrun.** It also raises the value of `233`'s
October question — if the field's cosine is near its ceiling and everyone must rebuild on a disjoint
panel on Oct 22, **the transition is a larger differentiator than the last +0.003 on this panel.**

## 2026-09-04 07:18 UTC — 232 addendum: whitening costs 1.37 points of SIGN accuracy (never charged), and a harness-validity flag

**INBOX-1's question, answered.** The whitening operator costs **1.37 points of |t|-weighted sign
accuracy**, worth **−0.00547 overall**, and **no operator experiment has ever charged it** — 193
(regulariser), 212 (exponent) and 226 (fitting sample) all scored the operator on DIRECTION metrics.
The sign channel explains **68%** of whitening's whole `nmae` cost; the rest is magnitude.

**But it must be quoted next to what the operator buys.** On the same three folds, production's γ
against γ = 0 gives Δ`pds` **+0.06761 → o_pds +0.02520**. The sign cost is **22% of the `pds`
credit** — real, previously uncharged, and **the operator remains strongly net positive.** Not an
arm; a correction to the ledger.

### An independent validation of the whole §8.1e framework
232's measured `dnmae/dq = −0.2573` implies **`k = 0.1286`** under 188's identity
`nmae = 1 + k(1 − 2q)`. **187 measured the |t|-weighted median per-gene `k` at 0.12 by a completely
different route.** Two independent measurements, **0.1286 vs 0.12**. The sign-exchange machinery is
not a fitted convenience — it retrodicts a number nobody fitted it to.

### 🔴 RETRACTED 07:30 — the flag below does NOT hold. I routed it without checking it.
**232's premise is a misattribution and I amplified it.** Two independent reasons it fails:
**(1)** `AGENTS.md` §4.2's evidence is **macro `pds` on LODO** (0.7378 → 0.7873 → 0.7978, "+0.060,
three-fold consistent, replicated on an independent 4-line LODO (+0.043)") — it makes **no cosine
claim and no board claim**. **(2)** The cosine↔`pds` trade along γ is **212's measured, closed
result**: `d(overall mse*)/d(cos) = +0.3655`, `d(overall pds)/d(cos) = −0.3517`, i.e. *"the
whitening operator converts bulk cosine into retrieval at **96.2% of the board's own exchange
rate**"* — and 212 closed γ on exactly that basis. **232 rediscovered 212 from the other side.**
Sign agrees; the magnitude differs because 232 crossed the whole 0 → 0.65 range (through the rank
cliff 226 found at `n/p ≈ 1.4`) while 212 measured a local slope at production's γ. **Non-linearity
across a cliff is not a validity failure.** Retracted to 221 so it does not spend effort on it.
**What survives is one narrow question**: does any board pair isolate γ = 0 vs γ > 0? §4.2 is LODO
and a19 is γ-high, so possibly none does — worth knowing, not a contradiction.

### The flag 232 raised (SUPERSEDED by the retraction above — kept for the record)
On PLOLO the **un-whitened** consensus has a **HIGHER** aggregate cosine (**0.040855 vs 0.031299,
+30.6%**) and a much lower `pds` (0.71947 vs 0.78708). **`AGENTS.md` §4.2 calls whitening "the
largest measured gain" — established on the BOARD through `mse`, i.e. through cosine.** So **PLOLO
cosine and board `mse`-implied cosine disagree about the sign of γ's effect.**

232 ruled out the easy explanations: hubness does not account for it (`rr` 0.0031 at γ=0 vs 0.0026
at γ=0.65) and macro and aggregate cosine move the same way, so it is not a functional mismatch.

**This matters far beyond γ, which stays closed** (170 measured the vertex end-to-end; a19 lost
−0.0135 on the board). **Every arm this project prices goes through PLOLO cosine.** And the
PLOLO→board factor is *already* unidentified — 1.384 vs 1.000 on two arms, and 231 found a **7.84×**
spread over four pairs with the only same-class pair reading **0/0**. **A structural local-vs-board
disagreement is a candidate explanation for that instability.** **Routed to 221**, which owns the
local→board transfer table.

## 2026-09-04 07:15 UTC — 229: registry LANDED AND COMMITTED; the ARM is NOT BUILDABLE. And a live question about the gate that killed it.

### The registry is in, and it is byte-inert
Committed as **`0204168`**. All three mandate known-answer rows pass, plus INBOX-2's inertness row
and two new ones at **0.0e+00**. **`VCC_UNION_SOURCES` unset reproduces the shipped `h1union`
package bit-for-bit** — I verified that independently: `229-ka-h1union-g028.vcc` hashes to
**`f13ce40b…`** at 2,824,785,920 bytes, identical to what is on the board. It also added
`verify_union_source_hashes()`, which re-hashes every selected bank before a 25-minute build and
**raises** on drift — three of the six banks live under `experiments/*/work/`, a tree other agents
write to. **The infrastructure for shipping any future source now exists.**

### The arm is dead, two ways
**FINDING 1 — 225's construction is not expressible.** 225 priced its stack as a **second residual**
on top of `h1union` (`cd4@0.34 + h1pair@0.28`, **then** `+ stack4@0.35`). The builder carries **one**
residual at **one** dose, so naming those sources **replaces** the H1 pair instead of joining it.

**FINDING 2 — everything expressible dies on the exact `nmae` identity:**

| construction | rate card | **measured** |
|---|---:|---:|
| `ARM6` (all six sources) | +0.00268 | **−0.00219** |
| `ARM4` (225's four names) | −0.00102 | **−0.00589** |
| `ARM3` (H1 pair + `dld1_lowmoi`) | — | **+0.00067** (0.35× gate) |

**The rate card's `nmae` leg is not a measurement** — it is `0.000980 × (Δcos/0.01)`, **positive by
construction**. The measured leg is **the opposite sign and 58× the magnitude.** All three `nmae`
pools agree and all three cross 1.0.

### 🔴 One correction to me, and one open question I am NOT closing
**229 corrected my INBOX-4.** I wrote that `h1union`'s `nmae` "went UP +0.00588". That is the
**board** number and it is right; on **this instrument** `h1union`'s `nmae` went slightly *down*.
What rose by ~+0.0056 in the instrument was **top-89 ordering precision** (0.59752 → 0.60312).
**Both facts are true, and their disagreement IS the mis-signing.**

**The open question: this `nmae` instrument has mis-signed the board twice.** For the CD4 dose step
it gave **−0.000152** where the board delivered **+0.001762** — wrong in sign, 11.6× in magnitude.
229's defence is that the *magnitude* is decisive anyway: both calibration points sit at ~+0.0004–
0.0006 raw and `ARM6` is at **+0.01722**, i.e. **31–44× outside every point the instrument has ever
been checked at**, with an independently damning `e = −0.033`. **I accept that for the arm.**

**But 229 was asked whether the ARM is worth a slot, not whether the GATE is trustworthy** — and
that gate has now killed 228's arm, `ARM6`, and will kill every future one. **Only the board can
validate it.** `ARM6` is the sharpest available test (largest predicted `nmae` move), the displayed
row is cosmetic before Oct 22, and ~90 slots expire worthless. **Decision deferred to ~17h from now,
after 233 and 234 report** — recording the option now so the reasoning is pre-committed, not
retrofitted.

**Two verified packages are on disk and NOT RECOMMENDED:**
`229-stack6-g035.vcc` `7a3c2d3e…`, `229-stack4-g035.vcc` `efa5cf7b…`.

## 2026-09-04 07:08 UTC — 226: the operator is SATURATED in `n`. Its geometry is now closed on all three axes.

**More fitting rows do not help.** Lifting production's row count **+36%** (29,394 → 39,921; same
three donors, same prune rule, composition held fixed) moves PLOLO transfer cosine **+0.57% macro /
+2.13% aggregate** — **+0.49σ / +1.01σ of the operator's own subsample noise** — and **costs `pds`
−0.001616 raw (−0.000602 overall)**.

**The row axis has exactly one real feature and it is a CLIFF, not a slope.** Below `n = p = 17,942`
the operator collapses (cosine **−35% to −53%**, `pds` to 0.586); above **`n/p ≈ 1.4`** it is flat to
the noise floor. **Production sits at 1.638 — past the knee.** 19 eigendecompositions of a 17,942²
covariance, run strictly one at a time.

**So the operator's geometry is closed on all three axes**: the **regulariser** (193 — "the ridge is
optimal"), the **exponent** (212 — a ratio argument, therefore instrument-independent), and now the
**fitting sample** (226). Adding sources buys nothing here either, which matters because it was the
one place the closed source space (225) might still have paid.

**Reported honestly, and worth keeping:** 226's own 11.0 GiB guard **fired once** at 12.30 GiB —
fancy-indexing a memmap pulls every touched page into process RSS. It changed the reader to plain
`read()` (page cache, not RSS) rather than raising the cap. That is the right fix and the right
disclosure.

## 2026-09-04 05:56 UTC — 228: the per-gene gain SURVIVES the read-back and is killed by `nmae`. NO-GO, and it hands us a design rule.

**224's diagonal is real after the scorer read-back — and unshippable.** 334 scored arm×fold rows,
preregistration committed before the first fit.

| quantity | measured |
|---|---:|
| relative aggregate cosine, after read-back | **+4.55%** (224's effect-level number was +8.9% — **156's warning fires at 1.9×, not ∞**) |
| Δ`pds`, 3 cancer folds | **+0.00484**, positive **3/3** |
| **|t|-weighted SIGN accuracy** | **+0.00027** — nothing |
| Δ`nmae` raw | **+0.0373** |
| efficiency `e = ΔM/ΔU` | **0.072** vs break-even **0.5** |
| six-member composite | **−0.00702** (best of ten doses: −0.00128) vs gate **+0.0019** |

**None of it is amplitude** — re-solving α to match the shipped stack's delivered norm (α 0.982 →
0.889) leaves the cosine and `pds` gains unchanged or larger. **All four of 224's nulls reproduce**
through this instrument, including the wrong-recipient signature (wrong **+1.71%** vs right
**−6.19%**). **And it fails under every free choice taken jointly** — un-validated pool, larger
`k_cos`, absolute `jac` transfer, dose picked on the evaluation folds — best composite anywhere
**+0.00177**, still under the gate.

### The rule this gives us — now `AGENTS.md` §8.1d
**A cosine gain earned by reshaping per-gene MAGNITUDES is destroyed by `nmae`; one earned by
improving SIGNS is not.** `nmae`'s matched term is sign-matched `min(|pred|,|truth|)`, so it
**cannot** collect a magnitude rebalance — the unmatched term collects all of it.

- **Signs → shippable.** `h1union` added a donor *estimate* and **`nmae` rose +0.00588** with all six
  members up.
- **Magnitudes → dead.** 224's gain, 158's per-target scale (`nmae` −0.0245 at oracle), 169's
  sharpening, 187's amplitude. **Four independent arms, one wall.**

**Ask this before building any direction arm: does it move signs or magnitudes?**

It also **rescues the retracted §4.1c correlation by a different route** — my board-level
`corr(cosine, nmae) = −0.394` was a clamp artefact and stays retracted, but 228 reproduced the
underlying tension **inside our own pipeline, causally and dose-responsively**. The coupling is
through **magnitude**, not a separate `nmae` axis.

## 2026-09-04 05:52 UTC — 231: Arm A REFUTED on arithmetic, and the ORDERING CHANNEL IS CLOSED

**I opened this lead and it closes negative — but it closes a whole channel with it, which is worth
more than the arm was.**

### Arm A is refuted twice, before any re-measurement was needed
- Its **+0.0088** requires a board raw `reach` of **0.31223** against a **world best of 0.290835** —
  **250% of the entire remaining gap to the best `reach` team on Earth**, from removing one
  retrieval component. Even 214's own ×1.533-deflated **+0.0054** implies **0.29983**, still above
  the field best.
- Its transfer factor — *"which is the whole arm"* (219) — is **NOT IDENTIFIED for its class**. Four
  consecutive single-variable board pairs give **0.390, 1.2927, 1.4346, 3.0563** (**7.84×** spread,
  geometric sd **×/2.343**), and the only pair in Arm A's own class (`cd4f034 − rc_callset`) is
  **0/0**: h1-test moved **exactly zero** (bit-identical effect matrices) while the board moved
  **+0.023416 raw**. **n = 0 usable pairs.** 214's headline borrowed 1.290 from a different class.

### 🔴 The finding that outlives the arm: 220's bound has shrunk 65%, because WE shrank it
220 priced the whole ordering channel off our **`rc_callset`** row (raw `reach` **0.24988**). **The
shipped `h1union` row reads 0.276545 — `reach` world rank 7 of 656.** The remaining prize for
becoming the best `reach` team in the world is now **+0.002689 overall, not +0.0077.**

**Three things close the channel, and they compound:**
1. **+0.002689 is BELOW the 0.004502 ordering noise** (230). The entire remaining prize is **0.6 sd**
   — taking *all* of it would not be distinguishable from a re-roll.
2. **Any amplitude-changing ordering arm is dead at any transfer factor**: floor **+0.0089** against
   a total prize of **+0.002689** — **the floor is 3.31× the whole prize.** Arithmetic, so it holds
   **even at oracle**.
3. **No usable transfer factor exists for the class**, per above.

**`AGENTS.md` §8.1c now says: do not open another ordering arm.** Not keys (220), not the call set
(192), not a renderer rider (230), not a component removal (231). **Every remaining point is in
DIRECTION.**

### Where that leaves the board
Three closures in ninety minutes — 230 (Hamilton, below noise), 231 (Arm A, refuted), and 226
converging on noise. **Both live arms are DIRECTION arms and therefore measurable** (§8.1b: `pds`
and `mse` reproduce to sixteen significant figures): **229's `ARM6`** (+0.00268 PLOLO, 3/3 cosine,
3/3 `pds`, built) and **228's per-gene gain** (pricing now). Those are tomorrow's two slots.

## 2026-09-04 05:45 UTC — 230: rule 13's `reach` toll is MIS-KEYED. Hamilton survives it and is still not shippable.

**27 known-answer rows, all pass.** The question I posed (step or slope?) was **mis-posed**.

- **Rule 13's own calibration points sit in the wrong row.** a21 and a22 deliver per-gene panel
  totals **identical to the incumbent's on 900/900 blocks** — by the table's own digest test they
  belong in the **zero-toll** row. **The row the `−0.0146` is charged to contains ZERO board
  points.** It has never been charged to a scored arm.
- **The digest is anti-correlated with what it predicts**: it fires on `cd4f034` (realised
  **+0.026** `reach`) and not on a21 (realised **−0.018**).
- **The cause is separated at an exactly-zero profile coefficient.** 100% of the `−0.0146` belongs
  to what a21/a22 changed *besides* the profile — the **cell channel** (emitted library CV
  **0.019 → 0.382**, `n_pred` 4,420 → 1,141, recall 0.669 → 0.271).
- **Not a call-set toll either.** `rc_callset` changed the call set by **more** (jaccard 0.2072,
  `n_pred` 4,480 → 2,599) at zero cell-channel change and realised **+0.022590**. **Size does not
  carry the sign; quality does.**
- **A profile change is not, by itself, a `reach` cost** — a23, `h1union` and `cd4f034` average
  **+0.007369**.

### 🔴 The finding that outranks the toll, and it changes what we should ship
**The unexplained spread of board Δ`reach` among DIGEST-IDENTICAL arms is sd 0.027011 scaled =
0.004502 overall — 2.4× the entire slot floor.** No toll of any shape can be estimated to the
precision this gate needs. **So "stop shipping small arms" now binds harder than the toll ever did.**

**Applied to Hamilton itself:** at the correct toll of 0 it reads **+0.002203 against +0.0019 —
1.16×** — and that **0.0003 margin is ~1/15th of the unpriced noise**. **It is a live candidate, not
a cleared one. It should NOT ride a slot on that number.** I opened this lead; the honest close is
that it does not clear.

**Nothing else re-opens** — the charged row never contained a scored arm, so the toll's only
realised effect on this project was the veto on 208's rider. `AGENTS.md` rule 13 is corrected:
**charge no `reach` toll; use the digest to classify an arm, never to price it.**

## 2026-09-04 05:20 UTC — an UNHARVESTED arm found in 219's audit: the Hamilton rider, worth +0.0022, blocked on one untested question

219 re-priced every closure and found **three re-open (208, 214, 200)** — then its own mandate
**forbade it from re-measuring any of them** (*"your job is to say what is worth re-measuring"*).
**Only 208's rider *rule* was harvested (by 223/225); the arm itself never was.**

**208's Hamilton rider, corrected.** 219 §2.2: 208 and 200 both measured on
`VCC_MASS_FIX=minnorm`, a renderer **only a23 ever shipped** — and a23 scored −0.00073 and was
abandoned. On the renderer that actually ships, 214 measured the factor at **1.84×**:

| | value | vs the +0.0019 floor |
|---|---:|---|
| 208's original rule-13 composite | +0.001197 | VETO |
| **corrected to the shipped renderer** | **+0.002203** | **CLEARS at 1.16×** |
| six-member, corrected | **+0.003909** | |
| …with the **−0.00243** emission-law toll charged | **−0.000227** | fails |

**The arm lives or dies on rule 13's toll alone**, and the toll has only ever been read as a
**class constant**. Its two calibration points (a21, a22) moved the profile *catastrophically* —
a21 took `n_pred` **6,034 → 1,913** and `nmae` **−0.064 scaled** — while Hamilton moves the
challenge tilt from **−0.01253 to +0.00030 log2**, a **0.03%** change that *removes* a systematic.

**230 is testing whether the toll is a STEP or SCALES with profile-change magnitude**, fitting
against every scored arm with a known emission law (a21, a22, `rc_hybrid`, `rc_callset`, `cd4f034`,
`h1union`) with an out-of-sample retrodiction. **Either answer is valuable, and one is much larger
than Hamilton**: a21/a22 changed the **call set** as well as the profile, so if the −0.0146 is
really a call-set toll misattributed to the profile, **every arm rule 13 has ever charged re-opens.**

**Also settled from that audit:** 158 — the near-reversal — stays dead, and 219 kills it by a
*different* route than 158's own (165 §7's 12–48× Orion-pair deflation, vs 158's `nmae` argument).
Two independent routes to the same verdict.

## 2026-09-04 05:10 UTC — 227: `nmae` is CLOSED, my "separate axis" premise was a clamp artefact, and 187 was right all along

**I launched 227 on a hypothesis that was wrong. It found the error, and the error is the durable
lesson.**

### What I claimed at 04:33, and why it was wrong
I read `Aginglab.com` (rank 2) as beating us on a **separate `nmae` axis**, citing (a) their lower
*implied* cosine and (b) `corr(implied cosine, scaled nmae) = −0.394` over the 55 teams above 0.15.

**I applied 197's inversion outside its documented validity domain.** `AGENTS.md` already stated the
inversion is *"one-sided, vacuous above raw `mse` 1.02361."* **15 of those 55 teams report
`score_mse` exactly 0.00000, at raw `mse` from 1.935 to 4.00** — 2–4× past the clamp — so the
inversion returns "cosines" of **−0.44 to −1.43** that were never measured. **Drop them and the
correlation is +0.022, not −0.394.** I verified this independently before accepting 227's refutation.
The "catastrophic cosine, great `nmae`" cluster is simply teams who **never α-shrank**, sitting at
`k ≈ 1` where the identity gives `nmae → 2(1−q) ≈ 0.89`.

> ### ⚠️ STANDING GUARD — only invert a row whose `score_mse` is strictly positive, and say you did.
> A clamped row carries **no** cosine information. Every board read in `AGENTS.md` §4.1c depends on
> this, and every rival comparison must state that the guard was applied.

### Aginglab is reconciled, not exceptional
A **2.4% smaller jackknife claim** than ours — a tighter emission, independently indicated by their
`jac` −0.033 / `fid` +0.020 pair — makes their board row **exactly what a submission with +12.3%
better transfer cosine prints.** 197's inversion misreads a 1% claim difference as **5.20%** of our
cosine. **They are ahead on our axis; the instrument could not see it.**

### `nmae` is closed on all three candidates, two of them structurally
- **Unreachable at fixed cosine.** Over four sign-fixing orderings — **two of them oracle** — the
  minimum cost of reaching their `nmae` credit is **+5.0% relative aggregate cosine on h1-test and
  +12.3% on PLOLO's three cancer folds**, against a ±0.2% gate. **Failed by 25–60×.**
- **Emission: closed structurally.** The quasi law's per-gene totals are `floor(N·L·profile)` —
  **independent of `spread`** — and every emitted cell carries the same library total, so the
  arithmetic mean of per-cell CPM **is** the design pseudobulk. Measured: ≤ **2.8e−3** raw over a
  16× spread and 16× dispersion range; **≤3.6e−4** at production's own dispersion.
- **Call set: closed by definition.** `nmae`'s gate is the **reference's** significant set, and our
  predicted lfc is a function of the design profile alone. **This retroactively explains a fact we
  had never accounted for: `rc_callset`'s +0.0084 overall cost `nmae` exactly zero.**
- **The largest channel is not estimable.** A per-row additive offset worth **+0.0542 raw at
  oracle** — larger than the whole gap to Aginglab — correlating **+0.975** with the reference
  gate's median true lfc and **−0.029** with the whole-row shift. Best honest fit on our own
  observables: **R² = −0.66.** It is a property of the reference's gate selection, not of our
  prediction.

### Net
**There is ONE route, not two: transfer cosine.** 187's original verdict — *"`nmae` is a direction
problem, not an amplitude problem"* — stands, and my 04:35 amendment of it was the false step.
`nmae` is handed back to the direction programme. No flag changes; nothing here is an arm.

## 2026-09-04 04:48 UTC — 225: the public-source space is CLOSED, sources SATURATE at 1.8×, and BOTH readings of the leaders' climb are now dead

**225 swept the space exhaustively and its pre-committed decision rule fired on the negative branch.**

### The census is closed
83 assets scanned from their own bytes, plus **GEO** (192 human CRISPR single-cell series,
2025–2026), **BioStudies/ArrayExpress** (191 hits), **CZI CELLxGENE** (388 public collections,
**zero** CRISPR perturbation screens), **Zenodo** (which 166 could not reach), **Hugging Face**, the
**IGVF portal** (98 released files) and the **Arc Virtual Cell Atlas**. **Nothing panel-covering
exists that we do not already hold.** X-Atlas/Pisces is still an empty repo (`usedStorage: 0`,
unchanged since 2026-03-25, re-checked).

### Phase 0: 20/20 gates
223's five headline rows reproduce to the last digit (`dld1_lowmoi@0.35` **+0.00228**, `@0.28`
+0.00224, `telohaec@0.28` −0.00120, `cmap_xpr@0.28` −0.00077, `kolf2@0.17` −0.00174); 218's three
baselines match at **Δ = 0.0e+00** on all three cancer folds; the known answer reproduced
**bit-identically across three separate processes**.

### 🔴 SATURATION — the number that kills the hypothesis
Cumulative price over the shipped `h1union` arm: **1 source +0.00228 → … → 5 = +0.00421, i.e.
1.84× DLD-1 alone.** The **fifth source buys 10% of what the first bought**. **Thirteen sources are
WORSE than five** (+0.00403 vs +0.00421), and **twelve sources without DLD-1 do not clear the gate
at all** (+0.00137 = 0.72×) — removing DLD-1 costs more than removing all four others combined.

225's mandate pre-committed: *"If five sources give 5× DLD-1 the leaders' shape is explained… If
they saturate at 1.5×, the hypothesis is dead."* **Measured 1.84× at five, 1.77× at thirteen. The
public-source ladder cannot produce `cqawesome`'s five separately-resolved steps.**

### 🔴 So BOTH readings of 215 are now dead
- **224** killed *iterative retraining of a learned transfer model*: every capacity above a per-gene
  diagonal loses, and the slope in `k` is negative.
- **225** killed *sources added one at a time*: saturation at 1.84×.

**The leaders' climb is now unexplained by any hypothesis this project has tested.** That is the
honest state. Note that neither closure weakens what each experiment *found* — a shippable per-gene
gain (224 → 228) and a shippable four-source stack (225 → 229).

### One new source, real and far too small
`data/external/igvf-cardiomyocyte/…` carries **87/300** panel targets — more than DLD-1's 67. 081
tested it in Aug 2026 through a *retrieval-gallery* adapter and closed it; **it had never been
through the post-whitening residual arm** (208's rider rule). It passes all three falsifiers and
prices at **+0.00038 = 0.20× the gate**, positive on **1/3** folds — a single-recipient effect, not
a source. Cause is visible: a **significant-hits-only** deposit, **median 6 non-zero genes per panel
row** against DLD-1's dense rows.

### 🔴 A correction: `sunshine_calu3` is DEAD, and I flagged it prematurely
Reading 225's *in-flight logs* I saw `sunshine_calu3` at +4.59% relative cosine with a negative
shuffle control and called it "stronger than DLD-1" **before its falsifier suite finished**. The
completed suite kills it: **it fails F1 — its own target-shuffled control scores higher.** Final
price +0.00148 (0.78× the gate). **Do not read conclusions out of a running agent's logs; wait for
the falsifiers.** That is why the falsifiers exist.

### A leak found and repaired rather than argued away
Jiang 2025's six lines **include K562**, the held line on the `k562_gwps` fold. Dropping that one
column (`jiang_5line_noK562`) **raises** the price +0.00079 → **+0.00096** — the contamination was
costing, not flattering.

### What rides
**`dld1_lowmoi ∪ jiang_5line_noK562 ∪ igvf_cardiomyocyte ∪ myllia_a549` at `g ≈ 0.35`** —
**2× the shape-only gate**, 3/3 folds positive, F1 and F3 passed at every dose, against DLD-1 alone
at 1.07×. **Nearly double 223's recommendation, for zero new data.** **229 is building it.**

## 2026-09-04 04:45 UTC — 224: a trained model DOES beat the donor mean, it is the SMALLEST one, and the mechanism is a statistic we already owned

**224 is the first positive answer to "can any trained model beat the donor mean on held-LINE
transfer" — and it answers it at the bottom of the ladder, not the top.**

### The result
A **per-gene diagonal** — one scalar per gene, ~**5,070 parameters**, fitted by inner
leave-one-training-LINE-out on `(donor consensus → held-out line)` pairs:

- **+5.16%** relative transfer cosine on the 7-line surface, 95% CI **[+4.63, +5.72]**,
  **7 of 7 held lines positive**
- **+8.9%** at production's `k = 3` **on the 300-gene challenge panel**, at **Δ`pds` +0.0118**
- **+9.63%** cosine and **+0.0124** `pds` **on top of production's own whitening** `γ = 0.65` —
  *more* than it adds to the unwhitened consensus, so it is **orthogonal to the operator**
- Pre-registered gate was **+2.9%**. It clears by **3.1×**.

### Every higher-capacity rung LOSES, and the stop rule fired
- feature ridge **replacing** the consensus: **−38.6%** (132 reproduced)
- feature ridge **stacked** on it: **+7.8% on held perturbations, −9.1% on held lines** — 195's wall
  reproduced quantitatively
- a **25M-parameter** learned gene-space operator: **+11.5% cosine at −0.032 `pds`** — it buys cosine
  by **inflating the decoy**
- MLP / factorised bilinear over (target embedding × recipient basal): **+4.4%**, no better than the
  plain diagonal it contains, **every recipient-conditioned parameter driven to zero**
- **the slope in `k` is negative or flat on every rung that clears** — the diagonal peaks at `k = 3`
  (+7.80%) and decays to `k = 6` (+5.16%)

**So capacity is measurably not the binding constraint: a 5,070-parameter diagonal beats a
25M-parameter operator on the endpoint that matters. No GPU was rented; $0 spent.**

### The mechanism is identified, and 189 already measured it in the wrong channel
The fitted gain **is 189's per-gene cross-line reliability prior**, recovered by a completely
different estimator: `r(gain, θ)` = **+0.81 / +0.90** on the two surfaces, `r` with 224's own
reliability statistic **+0.94**, and gains fitted for **different held lines correlate at r = 0.95**.

> **189 measured this statistic in the ORDERING channel, where it bought nothing. Nobody had put it
> in the DIRECTION channel.**

Four nulls pass: gene-shuffle **−7.2%**, moment-matching **−3.5%**, label-derangement **−25.7%**,
wrong-recipient **better than the arm**. Leak audit clean (0 panel coordinates in any gene axis,
`h1_2025` in no surface).

### The panel negative was POWER, not substance
7,312 per-gene gains cannot be fitted from 73–120 panel targets. Fitting on **1,758 non-panel
targets** and evaluating on **232 panel targets** of a held line recovers **+8.3% ± 0.3 over five
seeds**. **That also makes the arm panel-independent by construction — the only candidate we have
that is October-robust by design.**

### What it closes
**215's "iterative retraining of a learned transfer model" reading is NOT supported by any capacity
above a per-gene diagonal.** Recipient conditioning on CRISPRi is dead a second time (7 CRISPRi
lines, wrong-recipient null beating the arm, an MLP zeroing its own conditioning heads) — **184 is
now over-determined.** Basal-similarity donor weighting stays closed (010).

**So the leaders' climb remains unexplained by learned models, and 225's source-stacking reading is
the last one standing.**

### 🔴 It does NOT justify a slot, and 224 says so itself
This is an **effect-level** measurement on donor tables — not run through `profile_from_effect` +
the scorer read-back, so per 156 it measures **an object the board never sees**; not a PLOLO number,
so not priceable through rule 13; and neither the transfer factor nor §8.4's 2.51× multiplier is
measurable on that instrument. **228 is running 224 §10's specified successor**: wire the gain as a
per-gene multiplier, fit LODO on non-panel targets, and price **all six members** through PLOLO —
including the two things 224 could not see, the `VCC_ALPHA` interaction (211 has `claim` at 100.0000%
of its cap, so the gain may be spending amplitude we do not have) and the **DE call set moving**,
which makes `fid`/`jac` priceable rather than zero.

## 2026-09-04 03:45 UTC — BOARD: rank 5 of 656 at 0.20492. Ahmet's 0.2 target is cleared.

`cd4f034` scored **0.2049196** (model `The Ribosome Has Filed A Grievance`, submitted
2026-09-04T00:00:06Z, description empty). Up **+0.0130557** from 0.1918638. **Rank 5 of 656; only
five teams are above 0.20 and only one above 0.22.** Per-member world rank improved on every
member: `jac` **1st**, `mse` **10th** (was 26th), `reach` **10th** (was 29th), `pds` 7th, `nmae`
59th (was 83rd), `fid` 129th (was 194th).

### The gap to 1st is cosine and almost nothing else

| member | ours | cqawesome | gap |
|---|---:|---:|---:|
| `pds` | 0.7378 | 0.8264 | **+0.0886** |
| `mse` | 0.1593 | 0.2540 | **+0.0947** |
| `nmae` | 0.1049 | 0.1211 | +0.0163 |
| `fid` | −0.0201 | −0.0598 | −0.0396 |
| `reach` | **0.2193** | 0.2109 | −0.0084 |
| `jac` | **0.0284** | −0.0067 | −0.0351 |
| **overall** | **0.20492** | **0.22431** | **+0.01939** |

We are **behind on all three direction members and ahead on all three ordering members.** 197's
inversion reads implied transfer cosine **ours 0.09018, cqawesome 0.13356 — a gap of +0.04338
(1.481×)**, worth **+0.01583** at §8.2's exchange rate against an actual gap of +0.01939:
**cosine explains 82% of the distance to first.** See `AGENTS.md` §4.1c.

### MT Zeng holds the board's best direction and ranks 8th

Implied cosine **0.14138** — above cqawesome's — for overall 0.19332, because their ordering
members throw it away. **We rank 5th on a cosine 36% below theirs.** A point of cosine is worth
more to us than to anyone on the board, and trading ordering for direction has a visible price.

### Closed today: the public-text route to the leaders' method

Mined the `description` field of **all 2,306 scored entries across 19 snapshots; 764 are
non-empty. Not one entry with implied cosine above 0.09 has a description.** cqawesome, the
Prima Mente cluster, Garrrbage and we all leave it empty — the whole top of the direction channel
independently converged on §14's discipline. **No public text will reveal their method; do not
spend another agent reading descriptions.**

The census still dates the field: `k562` 181, `replogle` 81, `gwps` 75, `jurkat` 48, `state` 47,
`ridge` 31, `x-atlas` 28, `cd4` 25, `esm` 24, **`marson` 6**, `tahoe` 4, `depmap` 4. The CD4 atlas
206 identified as the underexploited lever is named by **6 of 764** — it was genuinely underused,
which is why taking it paid +0.0131.

### What this leaves

Every non-cosine route is closed: amplitude (211, 187, 216 — three ways), dose (216: `f=0.34` is
the argmax over the whole 31-point grid; the cosine argmax costs **−0.01195** overall and even
there we reach only **0.799×** the leaders), `nmae` (187: "a direction problem, not an amplitude
problem"), `mse` (§4.1c: it *is* the cosine readout), ordering (220, bounded). **The only live
question is transfer cosine, and 224 (learned transfer) and 225 (public sources) are its two
standing hypotheses.**

## 2026-09-04 03:35 UTC — 216 FINAL: the atlas hypothesis is dead at byte level; the dose curve turns over at f*=1.018

**216 is NO-GO on every arm, and its Phase 1 is the most decisive negative we have.** The
hypothesis was that the leaders' doubling came from a refreshed CD4/atlas source. It did not.

- **It is the same object.** 12 byte ranges spanning 0 -> 16.79 GB match the live bucket;
  `LastModified` **2026-05-28**. **Zero objects on that bucket moved inside 2026-08-26 -> 09-03** --
  the exact window in which 215 resolved the leaders' climb. Whatever they did, it was not this.
- The one artefact refreshed 08-24 added **14 columns that were already `.obs` fields of the h5ad we
  already hold**. Panel coverage is identical to the published ceiling (291/276).
- The only materially different object is the pseudobulk (**17,771 vs 10,113** official genes, public
  since May). Its extra genes carry 11.94% of gated truth energy, but measured elasticity is
  **-0.55 to -0.97** -- **cutting that tail improves cosine**. Worth <= +3.2% at the most generous
  law, negative at the measured one.
- **Phase 0 passed 20 instrument gates**, rebuilt the four folds at max|diff| **0.000e+00**, and
  re-extracted the CD4 banks **byte-identical on 38 of 42 files** (the 4 differ only in local paths
  and peak RSS). It also re-derived 215's two leader windows from raw snapshots: `cqawesome`
  **x2.056 (rho 0.94, 12.4x the leak bound)**, `jw_13288` **x2.071 (rho 1.13, 14.9x)**. 215 holds.
- **Phase 2: the dose curve turns over at f* = 1.018.** Consistent with 216's earlier -0.00642
  pricing of f=0.80 and with `nmae`'s turning point.

**What this closes.** The "they got better data" reading is now excluded for the CD4/atlas route
specifically, on byte-level evidence rather than inference. Combined with 218 (independence) that
leaves two live readings of 215's many-small-steps shape: **iterative model retraining (224)** and
**public sources added one at a time (225)**. Both are now running.

## 2026-09-04 03:37 UTC — the October over-dose risk is MEASURED, not estimated

222's `six222.json` prices the sizing trap at both today's coverage and the pessimistic draw, on all
four PLOLO folds:

| fold | carrying | coverage | over-dose now | over-dose at q0.14 |
|---|---|---|---|---|
| k562_gwps | 251 | 0.933 | 1.0171 | **1.2436** |
| orion_hct116 | 267 | 0.921 | 1.0142 | **1.2453** |
| orion_hek293t | 271 | 0.922 | 1.0120 | **1.2381** |
| kolf2 | 253 | 0.927 | 1.0153 | **1.2369** |

Today's exposure is **1.2-1.7% -- negligible**, which is why the current dose work is unaffected. At
the pessimistic D/E/F draw it is **~1.24x across every fold**, tightly clustered. This is the number
`DEF-RUNBOOK.md` Sec 5 needed and did not have.

## 2026-09-04 03:20 UTC — 223 FINAL: the bug is a MEASUREMENT bug not a shipping bug, 0 of 21 void, and DLD-1 is live

**One correction to my mandate that matters more than the audit.** The **builder raises** on a
zero-coverage residual. The silent `return 0` is in the **harness mirror**,
`162/plolo.py:363-369`. **So this is a measurement bug, not a shipping bug** — production has never
shipped a silently-zeroed source; our *instruments* have been reporting zeros as if they were
measurements. That is a materially better position than I assumed, and it re-scopes the October risk
(below).

**0 of 21 VOID** — 15 VALID, 6 NOT-RUN. The reason is chronological: the median-sized residual is a
late construction (099 → 062 → 207 → 218) and these span 055–176. Two details worth keeping:
`105-telohaec-residual` **measured direct transfer and then closed the residual adapter without ever
running it**; and **`059-feng-ipsc` had already implemented the carrying-row rule with a loud zero
guard — thirteen days early.**

### The live source: DLD-1
Nine partial-coverage sources re-run with carrying-row sizing, marginal over 222's union, four PLOLO
folds. **Seven dead.**

| | |
|---|---|
| **`dld1_lowmoi` (67/300)** | +1.40% macro cosine, **Δ`pds` +0.0025 on 3/3 cancer folds**, price **+0.00228** = **1.20× the shape-only gate**, LOFO **+0.00188** |
| replication *(a falsifier I did not ask for)* | **both of DLD-1's independent libraries are separately positive, and the two-split mean beats both** (+0.565 / +1.080 / +1.403%) |
| `telohaec` | negative, and loses to its row-matched CD4 control at 0.43× — **105's NO-GO survives** |
| CMap (122/300, the largest source we hold) | **dead in all three variants** |
| `feng_ipsc` on top of the union | fails the row-matched CD4 control **at every dose** — the union already took its shared direction. This is the second branch of INBOX-1's question, answered |

**Δ`pds`-first ranking (INBOX-1 §1) selects `dld1` alone** — adding `myllia` costs `pds`.

### The third failure mode, quantified — and the October risk is smaller than I feared
Above 50% coverage the rule **over-doses**, exactly as derived: **`feng` by 1.265×**, with the
all-rows median landing at the predicted **18.1st percentile** of its carrying rows. **CD4 is
1.014×**, so **216's and 217's dose work is untouched** and `cd4f034` is unaffected.

The October exposure therefore stands only where D/E/F coverage falls: `docs/DEF-RUNBOOK.md` §5's
pessimistic draw puts CD4 at 174/300 (58%), which is *below* feng's 61% and so worse than 1.265×.
**That is a real risk and it is now sized rather than speculative** — and because the bug is in the
harness mirror rather than the builder, the fix is a runbook step, not a production emergency.

## 2026-09-04 03:05 UTC — 220 closes the ordering-key family and BOUNDS the whole channel; 223 refutes my silent-zero hypothesis and finds a live source anyway

### 220 — the ordering keys are closed, and both of my mandate's premises were wrong

**Premise 1 wrong: the 77x tightening is in the WRONG AXIS.** It reproduces exactly
(`0.00695 -> 0.00009`, 76x) but that is the **emission-seed** sem. 189's `+/-0.00424` and 190's CI
are **paired bootstraps over the 85 targets**, and the *target* se moved only **1.08x**. On the new
emitter the target se is **167x** the seed sem. Worse, for a counterfactual key the two instruments
are **the same object** - the per-target `reach` vectors for `|H|` and 189's key are
**bit-identical** on both DE tables. **189 and 190 were already on the "new" instrument.**

**Premise 2 wrong: 189's structural cap is off by 23x.** The measured transfer is
`d reach / d acc@10 = +0.143`; at `c = 0.30` front accuracy is **0.872 - below 0.90 - and `reach` is
already +0.0379**. The 0.90 bar binds **per perturbation**, not on the mean. 190 inherited the same
error. So the family is **neither** noise-limited **nor** structurally capped.

**It dies for a sharper reason that generalises.** The gate needs **+0.055** of acc@10; 189's key
delivers **+0.010** and 190's is negative - theta is **2.4x flatter exactly where the key reads**
(sd 0.140 inside the top-400 `|H|` genes against 0.340 outside). And **the honest ceiling of the
entire gene-level-key family is NEGATIVE**: in-fold gene oracle **+0.0968**, cross-fold **-0.0114**.

**Concentration kills the remainder:** **60% of 189's effect and 111% of 190's is one
perturbation**; dropping 189's best three targets gives **-0.00191**. And **190's real network is
1.2 sigma WORSE than its own sign-shuffle**, beaten 5-of-6 and 6-of-6.

**The loophole was real and is closed by learnability.** 189's gain *is* where `reach` lives (160%
inside the 31 targets carrying 94.4% of it), and a **perfect** target gate reaches **+0.00205
overall - just over the gate**. But no self-side stratifier correlates above **|r| = 0.10**, the
split-half-honest gate gives **+0.00011**, and `corr(headroom, incumbent reach) = -0.900`.

**What replaces the closure - a specification, not a verdict.** A key clears rule 13 **iff** it is a
per-pair confidence at **AUC ~ 0.639** (corr ~ 0.30) on *"is our sign right"*. Measured: `|H|`
**0.535**, theta **0.511**, expression **0.505**. Anything proposed later is testable in minutes.

### The bound that should drive allocation
**Becoming the world's best `reach` team is worth +0.0077 overall. `cd4f034` alone delivered
+0.0131.** The entire remaining ordering channel is **smaller than one direction arm**. 220 also
notes that **AGENTS.md 8.4's 2.51x multiplier does not apply to ordering arms** - they carry no
member at zero and do not move the cosine.

**And `fid` is bounded too — asked and answered 2026-09-04, because the rank invites the question.**
We are `fid` **rank 129 of 656** (scaled −0.0201; world best `The_Way` +0.0227), a gap of 0.0428
scaled = **+0.00713 overall** — comparable to the whole `reach` prize, and the obvious thing to point
an agent at. **It is closed, because `fid` decomposes into two parts we have already bounded
separately.** `fid` is direction-fidelity *yield*: a sign-accuracy factor, which is pure DIRECTION and
therefore information-limited exactly like `pds`/`mse`/`nmae`; times a yield factor
`min(1, n_pred/n_real)`, which is the **call set** — and 192 swept it at **7 keys × 10 sizes × 4
`n_real`** and never beat a15 by more than **+0.0041 raw = +0.0023 overall**. There is no third
mechanism. Our low `fid` rank is our cosine showing through a different lens, not an unworked lever.

**So the ordering channel is closed on both remaining members**, and every point of `fid` we want
arrives the same way every other point does: through cosine.

*(All 33 delivered arms had `pds`/`mse` pinned, digest 85/85, `nmae` |delta| <= 4e-7, toll 0 - the
profile-preserving identity holding for the fourth independent time.)*

### 223 (interim) - my silent-zero hypothesis is REFUTED, and it found a source anyway

**Phase 0 passes: the mechanism is exactly as 218 described, verified three ways.** A synthetic
residual covering **21%** of a real PLOLO fold solves to `f_per_unit = 0.0` and `fraction = 0.0`
**exactly**, effect matrix **bit-identical** to baseline. A 40% source is silently zero too; 60% and
92% are live.

**Phase 1: 0 of 21 are VOID** - 15 VALID, 6 NOT-RUN. **Not one ever called `solve_cd4_fraction` or
applied a median-over-all-rows rule to a partial-coverage residual**, across 100 Python files, every
median-norm hit read by hand.

**The reason is chronological.** The post-whitening median-sized residual is a **late** construction
(099 -> 062 -> 103 -> 207 -> 210 -> 216 -> 218); the 21 experiments span **055-176** and are older or
differently shaped. **The bug is real and has cost real measurements - on the `cd4`-path experiments
(200, 207, 218), not on these.**

**But re-running the residual arm those experiments never ran finds one live source: DLD-1, at
+0.0025 marginal over 222's union, `pds` UP, all four falsifiers passed.** Replication across
DLD-1's two independent libraries still in flight.

## 2026-09-03 — 184: TAHOE-100M CLOSES CROSS-LINE ADAPTATION IN THE STRONG FORM. The five CRISPRi failures were NOT sample-size limited.

`experiments/184-tahoe-transfer/RESULTS.md`. The most important negative this project has produced,
because it closes a family rather than an arm.

**Acquisition — the atlas was never needed.** The right object is `tahoebio/tahoe-de-rhaister`,
**1.2 GB** against a 120 GB budget, yielding a **1,137 treatment × 50 line × 2,000 gene** pseudobulk
tensor that is **99.96% complete** (every line shares every well, so line and treatment are
unconfounded), plus per-(line, plate) DMSO basals and 7,800 replicate pairs. Rejected as
unnecessary: h5ad 336 GB, cell-level parquet 338 GB, full DESeq2 pseudobulk 92 GB.

**A basal-conditioned transfer operator works spectacularly — on Tahoe.** Line-specific residual
cosine, leave-one-line-out over 50 lines:

| fitting lines `k` | 3 | 5 | 10 | 20 | 40 | 49 |
|---|---:|---:|---:|---:|---:|---:|
| arm | +0.2293 | +0.3159 | +0.3891 | +0.4244 | +0.4557 | **+0.4635** |
| context-free control | +0.0026 | +0.0113 | +0.0221 | +0.0435 | +0.0351 | +0.0316 |
| line-permutation null | +0.0195 | +0.0269 | +0.0316 | +0.0188 | +0.0089 | +0.0389 |
| oracle ceiling | +0.3896 | +0.4217 | +0.4595 | +0.4809 | +0.5148 | +0.5224 |

+0.1467 cosine over baseline on **49 of 50** lines, rising monotonically in `k` exactly as the
1/k argument predicted.

**And it is an AMPLITUDE law, so it cannot touch the thing we need.** Sign accuracy moves
**+0.0037 at best**; the oracle gain itself scores **−0.0003** on sign, because the winning family
is a per-gene gain and a gain is sign-preserving.

**The decisive row.** The identical measurement on our own CRISPRi donors (161's exact five lines):
arm **+0.0665**, context-free control **+0.0667**, wrong-recipient null **+0.0663** — *identical to
three decimals*, against an oracle of +0.1126. That also reproduces 161's 0.0592 / 0.0606 from a
different model class.

**And it is not a noise floor, which is the objection that kept this family alive.** Degrading
Tahoe to CRISPRi reliability (ρ 0.281 → 0.047) preserves the specificity: at ρ = 0.076, `k` = 5, the
arm reads **+0.0676 against a null of +0.0057 — 12×** — where CRISPRi at the same `k` and the same
reliability reads **1.0×**. Model-free confirmation: the correlation between the line-specific
per-gene gain and line-specific log basal is **+0.666 on Tahoe** (still +0.640 after 6× noise)
against **+0.166 on CRISPRi**.

**So chemical perturbations carry a basal-dependent gain and genetic ones do not.** The five
failures (104, 160, 161, 174, 175) were not underpowered — the object they were trying to learn does
not exist in the receiving data. **Cross-line adaptation is closed. Do not reopen it without a
mechanism that is not a basal-conditioned gain.**

Consequence for allocation: with 187 closing calibration, 184 closing line adaptation, 188 closing
the magnitude map, and 189/190 closing the ordering key, **the only two live levers are the
ordering DELIVERY (185/191, in flight on the board) and better estimation of the shared component
itself** — which is what experiments 193, 194, 195 and 196 now own.

## 2026-09-03 11:20 UTC — `rc_hybrid` BUILT, all six falsifiers PASS, and SUBMITTED

`experiments/191-rc-hybrid-build/RESULTS.md`. Preregistration (unedited, committed before the
build): `experiments/185-rank-controlled-emission/PREREGISTRATION-rc_hybrid.md`.

| package | sha256 | bytes | nnz | % of cap |
|---|---|---:|---:|---:|
| `191-a15_local.h5ad` | `26fb3834…1bc85e3` | 19,049,370,200 | 2,377,814,190 | 50.06% |
| `191-rc_hybrid.h5ad` | `e0b0677f…d2c870a8` | 13,202,460,120 | 1,647,149,337 | 34.68% |

**Falsifiers, all PASS.** (1) per-gene panel-total digests identical on **900/900** blocks, nonzero
structure identical on **0/900** — so the two files are genuinely different and the pseudobulk is
genuinely the same. (2) Δ`nmae` **+0.00011 raw** against a 0.0010 bound. (3) challenge `n_pred`
**4,480 vs 4,474 = +0.13%** against ±2%. (4) Spearman(−log p, `|H|`) in the call set
**+0.9743 mean** across A/B/C against a ≥ +0.95 bound; expression correlation +0.279 → +0.043.
(5) 78 effects/operator fields, **0 differences**; the config diff is exactly the five pre-declared
emitter keys; operator sha and `after_all_median_norm` exact. (6) `vcc prep --dry-run` clean.

**Three findings the build produced that no analysis had.**

1. **`a15_local` IS a15 — byte-identical**, sha256 matching a15's recorded `output.sha256` across all
   19,049,370,200 bytes. So 170's arm64/x86_64 drift is real in the effect matrix (worst 6.5e−16;
   both *pre*-whitening digests exact, so it enters only at `whitener.apply`) but
   `floor(n·L·profile)` absorbs it entirely. **Falsifier 1 is therefore against a15 itself, not a
   local proxy** — a much stronger statement than the preregistration claimed.
2. **The incumbent's ordering defect is 2.6× smaller on the challenge than on h1-test**: a15 reads
   Spearman **−0.152** on A/B/C, not h1-test's −0.391, and the delivered span is 82% of h1-test's.
   **Part of the 0.4277 h1-test → board `reach` shrink may be this, not a transfer failure** — which
   is a different and more benign reading of a constant every `reach` prediction passes through.
3. The artifact is **31% sparser** than a15 (34.7% vs 50.1% of the density cap). Breaks nothing;
   noted because it is visible and because it leaves headroom (see experiment 196).

Submitted 2026-09-03 11:20 UTC, joke model name, **empty description**. Board result and the
falsifier adjudication go in `BOARD-RESULT.md`, not here.

## 2026-09-03 05:10 UTC — 188, 189, 190 ALL NO-GO. The ORDERING KEY is closed; only the DELIVERY survives.

Three experiments landed together. Each is a clean negative with a mechanism, and **two of them
falsify claims the coordinator made earlier tonight.** Corrections first.

### 🔴 CORRECTION 1 — the expression pilot is RETRACTED (189)

The 03:00 entry reported a coordinator pilot finding *"control expression is by itself the
strongest per-gene reliability signal available"* (held-out decile spread +0.0505) and concluded
that *"part of what 171 called a defect is signal"*. **189 shows that was a gate artefact.**

Replace the pilot's median-|effect| gate with a **significance-like** gate — the gate the scorer
actually applies — and `spearman(r_HE, expression)` runs **+0.5215 → −0.0696 → −0.3353**, while
the prior's own held-out spread *rises* +0.0402 → +0.1094 → **+0.1422**. On h1-test, inside `|H|`
decile 1, high expression signs **0.0724 worse** (P = 0.019), and every expression tilt loses
`reach`. 190 measured the same thing independently: `|H| × expression` is **−0.00961** below `|H|`
(189: −0.01001, the two agree to 0.0004).

**171 was right. The expression dominance of our p-order is a defect, not half signal, and 185
should keep targeting `|H|`.** The lesson is procedural: a magnitude gate and a significance gate
select different gene populations, and every per-gene statistic in this project must state which
one it used.

### 🔴 CORRECTION 2 — see the 02:40 entry: `nmae`'s `k ≈ 0.75` reading was also mine, also wrong

Two coordinator pilots falsified by subagents in one night, both in the same direction: a quick
screen produced a large, plausible number that did not survive a correctly-specified control.
**Coordinator pilots generate hypotheses; they do not enter the record as findings.**

### 188 — the L1-optimal magnitude map: NO-GO, and `|p|` is not a magnitude variable

The estimand needed correcting first (the L1 optimum is the conditional median of `r·sign(p)`, not
of `|r|`; the `|r|` version is the screen's worst arm at **−0.1497 scaled `pds`**). Then the shape
is the finding: **over a 330× range of `|p|`, median `|r|` moves +9.5% while sign accuracy moves
1.8 pp.** The map the data supports is nearly flat, and flat is catastrophic for `pds`. Every
fitted variant is strictly dominated by `sign(h)|h|^q`. **`|p|` is a sign-confidence variable, not
a magnitude variable** — which is why a calibration cannot rescue 169.

**And 169's 56× `nmae` disagreement is now resolved — against the arm.** It was **DE-pool
composition**, not instrument bias. The reference's BH pool is *significance*-selected, and
`|lfc|·√cpm` reproduces it at **0.7219** (against `|h|` 0.6406 and `|lfc|` 0.3796) at **median size
89** — which independently matches audit 181's board-implied `n_real` of 90–110. Priced on a pool
built that way, **PLOLO — a third independent surface — reads +0.02392 against h1-test's +0.04640:
same sign, within 2×.** The HCT116 branch that said GO was the wrong-pool branch. Deep pools invert
the sign of the `dU` term, and `e = dM/dU < ½` everywhere: **`nmae` opposes concentration
structurally.** Six-member sums for q1.15: −0.0094 (h1-test) / −0.0032 (PLOLO) against +0.0019.

**Two standing instrument rules fall out, both worth more than the arm:**
1. **Build a DE pool by significance, not by magnitude.** Use `|lfc|·√cpm` at the reference's own
   BH depth (median ~89 genes/row). A top-N-by-magnitude pool inverts the sign of `nmae` deltas.
2. **`debench.Surface.pds()` must never size a `pds` arm.** Converting both sides of the *same*
   h1-test table into the scorer's `h` space turns q1.15L's −0.02339 into **+0.00504**; the
   difference is localised to the 5–10 CPM band (11.0% of genes, 21.4% of the `|lfc|` norm). PLOLO
   was right, and the decoy count is *not* the cause.

### 189 — the per-gene cross-line reliability prior: REAL, ORTHOGONAL, NOT SHIPPABLE

The prior is unambiguously real. Built on 5 lines / 10 cross-line pairs with a hierarchical logit,
its LOLO held-out decile spread is **+0.1047 against a permutation null of −0.0003 ± 0.0025 —
z = +41.6**, rising monotonically with fitting lines (+0.0868 → +0.0979 → +0.1047). EB shrinkage
does nothing and the arithmetic says why: median `n_g` = 13,141, binomial SE 0.0064 against a
between-gene sd of 0.0667. It is **orthogonal to `|H|`** (spearman +0.031) and it *does* order
accuracy where `reach` reads: inside `|H|` decile 1 the top-θ tercile signs **0.6694 vs 0.5873**.

**And it still cannot move `reach`, for a structural reason that closes the whole key family:**
`reach` needs **90% purity**, the median pure prefix is **3 genes**, and 0.587 → 0.669 does not
cross 0.90. Correctly applied within `|H|` strata it is worth **≤ +0.0005 raw reach**. Priced:
ordering key +0.00408 ± 0.00424 (λ selected on the evaluation surface, so optimistic); call-set
filter +0.0088 precision at n=50 against a **shuffled control at +0.0091**; `q_w` +0.00026.
Slot gate **−0.0016. VETO.**

Also from 189, and it should govern every future confidence screen: **the in-fold oracle a naive
screen reports is 22.3× inflated.** The honest cross-fold per-gene ceiling on h1-test is
+0.01065 `q_w` = **+0.0019 overall**.

### 190 — the regulatory-prior ordering key: NO-GO, and `tydymy` does not retrodict

Best key over 8 networks × 13 λ: `reach` **0.23670 against `|H|`'s 0.23533** = +0.00026 overall,
CI [−0.0014, +0.0052]. **The sign-shuffled regulon, selected identically, reads 0.23699 — higher
than the real network.** Four of eight networks select λ = 0. The oracle over every pair CollecTRI
covers is **+0.00015 overall, 13× below the ship floor** — 189 derived the same cap independently
from the 0.90 purity bar.

**Mechanism: it is not coverage, it is accuracy.** On the exact pairs each network speaks about its
sign is *less* accurate than ours — CollecTRI **0.5862 vs `|H|`'s 0.6379**; propagated 0.5103 vs
0.5364; co-expression 0.5158 vs 0.5373. Coverage is also thin (41/300 panel targets are CollecTRI
TFs, 0.18% of pool pairs), but that is the second problem, not the first.

**And the board arm this was built on does not reproduce.** At `tydymy`'s published c = 0.25, 190
recovers 3 of 5 moves and gets **both gain members backwards**: `reach` −0.0131 (they +0.0353),
`fid` −0.0004 (they +0.0133). A shuffled network reproduces the cost side at 0.67σ. **M2 is
closed**, and a published competitor config is not evidence until it retrodicts on our own
instruments. (CollecTRI is **GPL-3.0**; OmniPath's academic and commercial queries return
byte-identical files.)

### Where this leaves the ordering channel

**The KEY side is closed.** `|H|` is the best key available; neither a per-gene reliability prior
(z = +41.6 and genuinely orthogonal) nor a regulatory prior beats it materially, and both are
capped by the same structural fact — `reach`'s 0.90 purity bar against a median pure prefix of 3
genes. **The DELIVERY side is open and measured**: 185 takes the delivered order from
Spearman −0.391 to +0.996 and realises 95% of the `|H|` counterfactual.

So the ordering channel is worth **the delivery alone**, ≈ +0.0028 overall if it transfers. That is
smaller than the 03:30 plan assumed, and it is why `rc_hybrid`'s justification is now **primarily
about pricing rule 13's toll for the profile-preserving class** — a standing constant that
currently vetoes an entire arm family — rather than about unlocking future keys. The slot decision
stands; the preregistration is unchanged and was committed before any of this landed.

## 2026-09-03 evening — FIVE MORE CLOSURES. The direction channel is now down to ONE untested lever.

`rc_hybrid` on the board (**0.183419**, rank 20 of 631) plus 197's inversion left **direction** as
the entire remaining gap (§8.2 in `AGENTS.md`). Five experiments then attacked it. All five are
NO-GO, each with a mechanism, and together they close the geometry and the estimator.

### 195 — fitting the transform to the retrieval objective: NO-GO, and the boundary is the cell line
**080 is this idea and is terminal** (its held-panel macro was −0.0000225, CI [−0.001863, +0.001803],
`promotion_gate.passed=false`); **084 measured nothing** — its RESULTS.md reads "implementation and
bounded benchmark pending". 195 re-ran on 70,496 held-line-safe retrieval problems against 080's 600.

**The direct objective picks `γ = 0.65` on its own** — donor-dev argmax 0.65, panel argmax 0.68
beating it by 8e−5. **Nothing to gain inside the family.** Outside it, **21 of 23 configurations lose
the line**: the pre-committed arm reads −0.00554 (CI [−0.01049, −0.00077], P(positive) 0.0118) and
`corr(Δdev, Δpanel) = −0.227`. Overfitting is *not* the mechanism — 1,088 parameters at 25 queries
each generalise fine to unseen *perturbations* (+0.0100) and still lose. The leave-one-training-LINE-out
decomposition is the finding: **+0.00904 across unseen perturbations, −0.00618 across an unseen cell
line**, and −0.00880 / −0.01273 with panel-shaped two-predictor queries, negative in 12/12 cells.
**Three lines is not enough diversity to fit a line-general correction** — the same wall 184 hit on
50 lines from the other side. Permutation: real beats deranged **0/12** on dev and loses **12/12** on
the panel; doing nothing beats what was learned.

### 193 — the operator's eigenvalue map: NO-GO, and two premise corrections
**The mandate's premise was wrong and 198 caught it first: production ships `WHITEN_LAM = 0.0`.
There is no ridge** — only `Whitener.fit`'s 1e−12 relative eigenvalue floor. And the shipped
operator is **29,394 × 17,942** (sha `e2f65ecf`), not §4.2's stated 12,000 × 7,760, which matches no
file on disk. **§4.2's dimensions are now known to be wrong.**

**There is no Marchenko–Pastur bulk to clip.** Condition number is exactly **1.000e+12** on all six
operators — the floor, binding — and the spectrum decays smoothly over **twelve decades**, so the MP
fit collapses (97–99% of eigenvalues declared spikes, 0.00% of the trace bulk) and clipping is the
identity (+0.000002 scaled).

**55 maps × 4 folds, nothing wins, and the best covariance estimator is the worst arm.**
`λ = 0` is the ridge optimum. Ledoit–Wolf linear shrinkage is algebraically inside the ridge family
and reads −0.000144. The **Ledoit–Péché RIE scores +82 to +93% PRIAL on this spectrum's actual
shape — a genuinely better covariance estimate — and reads −0.002552 scaled.** *Improving the
covariance estimate makes retrieval worse*, which is §4.2's "discriminability is not
signal-to-noise" measured a third way. Mechanism: below 1e−6 relative — 12.5% of dimensions and
**99.67% of the operator's gain** — the prediction holds **0.016%** of its energy, so a regulariser
would correct a quantity that is not there. `nmae` independently forbids regularisation: at
`q = 0.5375` and `k = 0.278`, `nmae = 1 + k(1 − 2q)` *improves* with larger `k`.

**Hand-off: basis COMPOSITION, not size, is the live lever.** Lifting the 12,000/donor cap to every
Orion row (+35% rows) costs −0.000305 raw, while the essential prune is **+0.000952** raw,
independently confirming 091.

### 194 — re-estimating the donor tables: NO-GO on eleven estimators, and one real hand-off
Phase 0 rebuilt `k562_gwps` **bit-identically** (sha `0c942228…`). Eleven estimators × four PLOLO
folds: best arm `kd_post@0.25` at **Δpds +0.00086, CI [−0.00040, +0.00210]** — and **its own
permuted-weight falsifier reads a HIGHER p(Δ>0)** (0.947 vs 0.907). Everything at full dose loses,
down to −0.01341 for the combination. **κ² does not move** (0.0951 → 0.0815–0.0978; the knockdown
arms *lower* it). Effective donor count unchanged.

**Two mandate candidates were structurally unavailable, which nobody had checked: 92% of targets in
all three donors carry exactly ONE guide construct**, so guide averaging, median-of-guides and
within-donor MODZ have no axis at all; and Orion's `pass_guide_filter` passes **100.0% of 7,943,468
cells**. `min_cells` is at its optimum (20 → 50 → 100 → 200 gives 0 → −0.00113 → −0.01678 →
**−0.08317**). Two incumbent defects found and both negligible: `build_orion_donor.py`'s `^neg`
regex matches **NEGR1** (154/233 cells pooled into the control; Δpds −0.000004), and the control
pool outweighs the median target 424–1,086× so `log1p` concavity leaves a shared offset (exact
correction ships as `data/donors/194-debias-*`, +0.00033, 2/4 folds — **do not swap it in**).

**The hand-off, and it is now experiment 204: the three donors' achieved repression differs 3.3×**
— median self-gene residual **0.155 / 0.300 / 0.510**, and only **4.4%** of HEK293T targets clear
the challenge's own >80% criterion against **47.6%** for K562. Equal weighting was established in
019 and re-confirmed in 022, **both before that measurement existed.**

### 203 — reopening 125: 198's diagnosis CONFIRMED in full, and the reopening is still a NEGATIVE
The rejection **was** a pure scale artefact — 231 assertions at tolerance 0.0 confirm the
re-implemented estimator equals 125's own, and the homogeneity and invariance claims hold. But the
estimator does not pay on the endpoints the competition reads: **−0.002372 raw `pds`**
(CI [−0.011331, +0.007291]) and **−0.00147 transfer cosine** (CI [−0.002163, −0.000752]) at its own
re-solved amplitude, composed gate **−0.001393** against +0.0089, with `nmae` a further −0.002468.
The direction-only variant 125 itself proposed loses to the shipped `|H|` key by −0.00281 at 20%
depth and −0.01429 at the reference's BH depth. **"The rescale control explains everything and the
estimator adds nothing" is the result** — the family is now closed for the right reason.

### Where the direction channel stands

| lever | status |
|---|---|
| more donor lines | **closed** — 166: no fifth panel-covering line exists; 197: the leader implies ~10.9 effective lines |
| cross-line adaptation | **closed** — 184, on 50 Tahoe lines, in the strong form |
| a learned transform | **closed** — 195, at the cell-line boundary, 12/12 negative |
| the operator's eigenvalue map | **closed** — 193, and the best covariance estimator is the worst arm |
| donor re-estimation | **closed** — 194, eleven estimators, κ² unmoved |
| the panel-level hierarchical estimator | **closed** — 203, for the right reason now |
| **donor mixing weighted by measured repression** | **OPEN — experiment 204, launched** |
| basis composition | open, small (+0.00095 raw measured by 193) |

**204 is the last untested direction lever with a mechanism behind it.** If it fails, the direction
channel is closed on every route anyone has proposed, and the remaining work is the call set
(202, on the board tonight) plus whatever 196 and 200 return.

### Process: four agents died on transient API errors and sat dead for an hour
`API Error: 529 Overloaded` and `500 Internal server error` killed 196, 200, 202 and 203 mid-run.
**`scripts/agent_supervisor.sh` now exists**: it reads each run's final `result` event, relaunches
only when `is_error` **and** the text matches a transient API/overload/timeout pattern, caps at
three attempts, and **never touches a run that wrote `RESULTS.md`** (203 had, and was correctly
skipped). It is armed as a persistent 2-minute monitor alongside a memory guard.

## 2026-09-03 END OF DAY — 0.1740865 → 0.1918638, rank 11 of 643, and the direction channel finally opened

**Two slots, both paid: +0.0177773.** Gap to `cqawesome` **0.04383 → 0.02895**. Fifteen experiments
closed today (192–210). Per-member world rank now: `jac` **1st**, `pds` **7th**, `mse` 26th,
`reach` **29th** (was 192nd), `nmae` 83rd, `fid` 194th.

| entry | arm | overall | Δ | what moved |
|---|---|---:|---:|---|
| `nQUJZhZb3MW2uFJHBX5e` | `rc_hybrid` | 0.1834190 | **+0.0093325** | `reach` +0.054940 scaled; `pds`/`mse` identical to **16 sig figs** |
| `31oB9yWy8CSNBDEkY4iW` | `rc_callset` | **0.1918638** | **+0.0084448** | `jac` +0.027589, `reach` +0.022590; `pds`/`mse` identical again |

### The three findings that reorganise the project

**1. The ORDERING channel converted, and better than any local instrument said.** `rc_hybrid`
realised **3.33×** its own prediction, and the h1-test → board `reach` factor for designed-order arms
is **1.290**, not the 0.4277 every prediction in this project had passed through. Rule 13's
`−0.0146` `reach` toll is **refuted for profile-preserving arms** — it belongs to arms that change
the delivered profile, and both its calibration points (a21, a22) did. `AGENTS.md` §11.1 rule 13
carries the amendment.

**2. The DIRECTION channel opened, on the one component we already ship.** 197 sized the whole
remaining gap as direction (`mse`+`pds`+`nmae` = 135% of it; we now *lead* on `fid` and `jac`) and
put the leaders at cosine ≥0.13 against our 0.0792. 206 then established the leaders are **off**
159's `1/k` curve — implied `k_eff` 10–14 against the **~11 CRISPR-perturbed human lines that exist
in total** — that the four "converging" accounts are **one team, Prima Mente**, and that the
documented exit is a **lineage-matched donor**. 207 measured ours: the CD4 dose was **3.4× too low**,
because 062 chose it by maximising LODO `pds` — the member least sensitive to direction. The cosine
maximiser is **+31.5%**. `cd4f034` is preregistered at `experiments/210-cd4-build/` and building.

**3. The Oct 22 trap.** 199 found the whitening operator is a **function of the panel** and fails
**silently** on a new one — 98.8% of coordinates paired with the wrong gene while every builder check
passes. See `AGENTS.md` §8.3 and `docs/DEF-RUNBOOK.md`. **This would have destroyed the only entry
that counts.**

### The eight closures, each with the number that killed it

| # | verdict |
|---|---|
| **192** | 188's membership key is an **oracle**; with our own `lfc` it reads fold 1.88 against control CPM's 2.05. `fid` is **closed to the call set** (7 keys × 10 sizes × 4 `n_real`, never above a15 by more than +0.0041 raw). `n_real_p` **not predictable** (LOTO r +0.24) so per-perturbation sizing is −0.097 as built. **But its union call set became `rc_callset`, +0.0084 on the board.** |
| **193** | Production ships **`WHITEN_LAM = 0.0`** — no ridge. No MP bulk (condition number exactly 1e12, twelve decades). The **Ledoit–Péché RIE scores +82–93% PRIAL and reads −0.002552 scaled**: improving the covariance estimate makes retrieval *worse*. §4.2's stated operator dimensions are wrong (it is 29,394 × 17,942). |
| **194** | Eleven estimators; best +0.00086 with its **own permuted falsifier reading a higher p(Δ>0)**. κ² unmoved. **92% of targets carry exactly one guide**, so three candidates had no axis. Handed off the 3.3× repression spread. |
| **195** | The direct retrieval objective picks **γ = 0.65 on its own**. Outside the family 21/23 lose, and the split is diagnostic: **+0.009 across unseen perturbations, −0.006 across an unseen cell line**. 080 is this idea and is terminal; **084 measured nothing**. |
| **196** | Depth strongly negative (raw `mse` 0.882 → 1.008 at m=2). Underneath it: `emit_counts` uses **`floor`**, a **−0.01234 log2 size tilt**; Hamilton removes 97% and is worth **+0.001182** with `claim` still saturated. **Free rider — 208 owns it.** |
| **200** | All ten components ablated: **six load-bearing, four inside gate, nothing to drop.** `--overdispersion`'s `reach` value belongs to the **(spread, φ) pair**, not to φ. The mass fix **never cleared**. Found the CD4 median-sizing silent zero. |
| **204** | Equal donor weighting is the **maximum of the transfer cosine over 29 alternatives**; the permutation falsifier **ties** the best arm; and the negative **has power** (model predicted +8.5–16.8%, measured **−1.50%**, CI [−2.38, −0.67]). `a_i` is **not identified** (four estimates span 0.504–1.254). Handed off the Orion↔Orion excess **subspace** → 209. |
| **205** | X-Cell's repo raises `NotImplementedError` in **both** public functions; Pisces is **0 bytes** since 2026-03-25; the preprint's own availability claim is **false**. X-Cell *would* have passed the structural test that killed State (**298/300** panel genes carry a STRING embedding against State's 0/300 over 2,024). A complete Pisces Jurkat drop is a genuine fifth line worth **+0.0062**. |

### In flight

| # | what | status |
|---|---|---|
| **210** | `cd4f034` build — the first direction arm to clear the gate, +0.0044 to +0.0057 predicted | building for the 00:00 UTC slot |
| **208** | Hamilton rounding + the tier-vs-flat trade, adjudicated with `rc_callset`'s `jac` board point | running |
| **209** | deflate the Orion↔Orion excess **subspace** — the last structural idea about our donor set | running |

## 2026-09-04 00:00 UTC — 211 and 212: amplitude and gamma both CLOSED FOR THE LAST TIME, and my premise was wrong

**212 falsified the mandate I wrote it from, and that is the headline.** I argued that §4.2's
measured *"whitening costs 1.5–2.1× on directional accuracy"* and 197's 1.66× gap were the same
number. Measured properly on PLOLO — 4 folds, 300 panel targets, challenge frames — **the operator's
entire cost on transfer cosine is 1.3575×** (per fold 1.124 / 1.481 / 1.616 / 1.182). Abolishing the
operator outright reaches cosine **0.10758, still 1.219× short** of the board's best direction, and
costs **−0.0687 raw `pds` = −0.02561 overall**. **The gap is not inside the operator.**

And the structural reason gamma can never pay: the operator trades cosine for retrieval at **96.2%
of the board's own exchange rate** —

```
d(overall mse*)/d(cos) = +0.3655     (exactly §8.2's rate)
d(overall pds )/d(cos) = −0.3517     (what whitening charges on the gamma axis)
                 sum   = +0.0138     — they cancel to 3.8%
```

— and that ratio does not move with the emission, the dose, or α. The two maximisers **do** differ
as predicted (cosine wants γ = 0, `pds` wants γ = 0.80), and it changes nothing: best on
`rc_callset` **+0.00013**, on `cd4f034` **+0.00042**, most generous anywhere **+0.00321** at γ = 0.40
against a **+0.0089** floor. Re-solving α at each gamma does turn the `mse` leg positive everywhere
below 0.65 (up to +0.00987 at γ = 0) — 170's V was partly an α-mismatch artefact — and `pds` eats
all of it. **Third and final closure.**

**211: the cliff did not move, and 208's number measures the other side of the `min`.** `Ĉ/C_p`
1.125 → 2.940 is the margin before **`claim` falls off its cap** (196's depth knob). The amplitude
cliff is the margin before **`s²B` rises to meet `claim`**, and `claim` is hard-capped because
`PRED_TRACE_CAP_K = 1.0`. **a15 was already at 100.0000% of that cap** (per-context `q_c`
1.000000 / 0.999904 / 1.000000). Cliff α moved **+0.0000% / +0.0048% / +0.0000%**.

The board had already proved it: **context A is the only one on the `claim` branch**, so
`rc_hybrid`'s and `rc_callset`'s Δ`mse` = 0.000000 **to sixteen significant figures** *is* a proof
that `q_A = 1` on all 300 perturbations. Re-solved optimum **α\* = 0.9808** against 187's 0.982 —
worth **+0.000076**. Clearing +0.0089 on amplitude alone needs **α 1.147**; the binding context
permits **−0.12%**. **Global amplitude closed for the second and final time.**

211 also corrected the mandate's own claim that `B` should move: `B` is a **pseudobulk functional**
and 270/270 blocks deliver identical panel totals, so `A`, `B`, `X` and `Kc` are all pinned.
And it measured what 210 did not: `cd4f034`'s `B` tracks its `A` to **0.05%**, opening **+0.60%** of
headroom on the binding context — optimum α **0.9879**, worth +0.00039. **The α that ships does not
depend on whether `cd4f034` lands.**

### Two by-products worth more than either verdict

1. **`cd4f034`'s dose raises PLOLO cosine by +13.83%** — the first *forward* measurement of the CD4
   dose on the cosine channel, independently corroborating 197's +6.8% at `f = 0.10`. That is direct
   support for the arm sitting built and waiting for the 00:00 UTC slot.
2. **§8.3's Oct-22 operator trap now has a price: −0.0429 overall on `pds` alone** — roughly a
   quarter of our entire score, silently, on the only entry that counts.

### Standing rule 211 adds
`claim` is at its ceiling in the one context that reads it, so **no future emission arm should be
re-priced on amplitude grounds.** 208's re-read instruction is correct for the depth knob
(196 §2g's `m = 2, φ = 8`) and **only** for that.

### Where this leaves the direction channel

| route | closed by |
|---|---|
| more donor lines | 166 · 205 |
| cross-line adaptation | 184 (50 Tahoe lines; the winning law is sign-preserving and absent in CRISPRi) |
| a learned transform | 195 (at the cell-line boundary, 12/12 negative) |
| the operator's eigenvalue map | 193 (the best covariance estimator is the worst arm) |
| **the operator's exponent** | **212 (the trade cancels to 3.8%, structurally)** |
| donor re-estimation | 194 (eleven estimators) |
| donor weighting | 204 (equal is the maximum over 29) |
| donor–donor correlation | 209 (the excess is the long-gene axis at rank 2) |
| the hierarchical estimator | 203 |
| **global amplitude** | **211 (second and final)** |
| **the CD4 dose** | **OPEN — `cd4f034` built, +13.83% cosine measured forward** |

**Everything except the CD4 dose is now closed. The leaders hold 1.66× our cosine on public data
(206), and it is not in any of the above.** The one structural assumption nobody has ever tested is
the **transport law itself** — experiment 213.

## 2026-09-04 00:25 UTC — 214: all three arms VETOED, and PerturbMatch is the OPPOSITE of redundant

`experiments/214-cheap-positives/RESULTS.md`. Every known-answer row passed (11 `pds` gates, 7 `mse`
gates, 2 bit-identical DE gates); `art_floor` reproduces `de202.rc_callset.npz` **bit-identically on
85 × 10,779 × 4 fields**.

**The question I cared about is answered with the opposite sign.** I suspected PerturbMatch had gone
net-negative because 207 tripled the CD4 dose and 064 had measured the two as synergistic — i.e.
that they were substitutes and 207's sweep was run against a redundant component. **They are
complements.** PerturbMatch's `pds` value **rises** with the dose (+0.002856 raw at `f = 0.10` →
+0.003797 at 0.34 → +0.004310 at 0.50), the interaction is positive at every dose, and the decisive
corollary: **the `pds`-optimal dose is `f = 0.20` with PerturbMatch and `f = 0.20` without it.**
207's sweep stands.

**And PerturbMatch never flipped.** Its effect-level value reproduces 058's +0.00517 as **+0.004785**
on PLOLO's 300-target panel. Three things changed, none of them the component: the **level**
(rendering costs it 40% of its `pds`), the **renderer** (200 priced it through a23's mass fix, which
has not shipped since a23), and the **accounting** (200 was the first six-member price on a
board-validated surface; 061's positive `nmae` came from the quarantined HCT116 harness).

| arm, marginal to `cd4f034` | `pds` | `reach` | `nmae` | six-member | verdict |
|---|---:|---:|---:|---:|---|
| **A** drop PerturbMatch | −0.001059 | +0.007621 | +0.002909 | +0.0097 (deflated **+0.0060**) | **VETO** — moves delivered amplitude **+4.28%**, so the floor is +0.0089 |
| **B** Hamilton | +0.001534 | +0.000014 | +0.001731 | **+0.0032** | **VETO** — digest differs on 900/900, carries the −0.0146 toll, and fails without it too |
| **A∘B** | +0.000125 | +0.006939 | +0.002843 | +0.0102 | **VETO** |

Composition **measured, not assumed**: `A+B = +0.011300` against `A∘B = +0.010153` — sub-additive by
−0.0011, consistent with 178's law that gains discount when composed and costs do not.

### Three by-products, each worth more than the arms

1. 🟢 **The ordering channel finally has a low-variance instrument.** `reach` under `rc_callset` has
   a three-seed sem of **0.00009** against 200's **0.00695** — a **77× tightening**. Effects that
   were "inside noise" on the old instrument are now resolvable, which **re-opens 189's ordering key
   (+0.00408 ± 0.00424) and 190's (+0.00137) as measurable questions** rather than closed ones.
   Combined with 215's finding that the strict ordering channel is **nearly unused board-wide**
   (9 arms in 707 pairs, two of them ours), this is the cheapest live territory we have.
2. **h1-test overstates any PerturbMatch reading by a measured ×1.533** — two donors, so the K562
   slot carries 1/2 instead of 1/3. **This applies to 200's numbers too.**
3. ⚠️ **200's and 208's prices were both measured on a renderer that has not shipped since a23.**
   Hamilton is **1.84× larger** on the shipped renderer (+0.002069 vs +0.001124); PerturbMatch
   **1.77× smaller** (+0.001065 vs +0.001885). **Every number in 200's ten-component ablation table
   inherits this**, so re-read it before using any of its rows to justify an arm.

### Two corrections to my own mandate
`PM_BLEND = 0.50` is a **hard-coded constant** (`make_submission_pm_cd4.py:572`) with no env
override, and the pinned-support gate would need re-pinning — so "drop PerturbMatch" is a **source
edit**, not one flag. And `cd4f034` was **not yet on the board** when 214 priced against it, so every
number above is marginal to a built-but-unscored incumbent.

## PLAN TO 2026-11-05 — written 2026-09-03 03:30 UTC. Read this before the dated entries below.

### The score has exactly two channels, and they need different work

Tonight's four closures (186, 187, plus the board re-solve and the 189 pilot) let the six members
be split by **what limits them**, which is more useful than splitting them by what they read:

| channel | members | what limits it | who attacks it |
|---|---|---|---|
| **DIRECTION** | `pds`, `mse`, `nmae`, and the *sign* half of `fid`/`reach` | information — 3 donor lines, κ² ≈ 0.34 shared, `1/k` law needs ~10 (159); 166 says no more panel-covering CRISPRi exists | **184** (Tahoe-100M cross-line operator), **188** (per-gene magnitude shape) |
| **ORDERING** | the *order* half of `reach`, `fid`, `jac` | mechanism — entirely ours, and we are −1.92σ below the field's own conditional expectation | **185** (delivery), **189** (per-gene prior), **190** (regulatory prior) |

**187 closed global amplitude** (joint `mse`+`nmae` optimum is a15's 0.982 to four decimals) and
reclassified `nmae` as ≈93% direction. **186 established that the field's entire `reach` spread is
ordering, not direction**, four independent ways. So there is no calibration work left: every
remaining point is either new information or better ordering.

### What each channel is worth, in board points

Gap to first (`cqawesome` 0.21790 vs our 0.17404): **+0.0438**. Empirical ceiling — field-best on
each of the six — **0.25463**.

| item | overall pts | confidence |
|---|---:|---|
| ordering: close our `reach` residual to the field's conditional line | **+0.0077** | high — measured on 41 matched teams |
| ordering: `tydymy`'s regulon arm, gross, with the `nmae` bill removed by delivering it as a key | **+0.0141** | medium — one board arm, one independent ablation |
| ordering: `fid` from a per-perturbation call set (`F_p = k/max(n_pred, n_real)` is per-perturbation and our policy is global) | **unpriced, plausibly large** | low — never tested by anyone |
| direction: Tahoe cross-line operator | **unbounded if it works, 0 if not** | low |
| direction: 169's `pds` lever with the `nmae` cost removed | **+0.0032** | medium-high — PLOLO, 4/4 folds |
| calibration | **≈ 0** | closed by 187 |

**Ordering alone plausibly reaches ≈ 0.19–0.20 (rank 3–5). First place needs ordering AND one
direction win.** That is the honest arithmetic and it should drive every allocation decision.

### The slot plan — we are under-using the only instrument that can price three of six members

2 slots/day × 49 days to Oct 22 = **~98 slots**, plus 28 blind ones in the dark window. We have
used 25 in six weeks. Rule 9 says `fid`, `jac` and `reach`'s normalisation are **unreadable on
every local harness**, and audit 181 says h1-test is not validated for `reach` on emission-law
changes. **The board is the only instrument that can price the channel we have just identified as
the live one.** A slot spent measuring it is a slot spent well; the standing rule is that a slot
must buy information toward first place, and this does.

Two designs worth slots once 185 can hit a prescribed `n_pred`:

1. **The `n_real` ladder.** `F_p = k / max(n_pred, n_real)` has its knee exactly at `n_real`. A
   preregistered 3–4 point ladder in global `n_pred`, emission-only so the pseudobulk stays
   bit-identical and `pds`/`nmae` stay pinned, measures the board's own `n_real` distribution — the
   number that has already cost four board predictions and that no harness can read.
2. **The key A/B.** Two builds differing only in the ordering key, everything else bit-identical.
   Reads the ordering channel's board exchange rate directly, which is the constant every future
   arm in this channel needs.

Both are emission-only, so `pds` moves by exactly zero (board-verified twice: v14 and a8).

### D/E/F endgame — the parts still open

Items 1 and 2 of the transition plan are done (180). Still open, and none can be done under time
pressure on Oct 22:

3. **Re-audit donor coverage** on the new `pert_counts.csv` the day it lands (166's `census.py`).
   K562 covers 269/300 today, Orion 290–294; D/E/F's targets are different.
4. **Verify the gene axis and `cells_per_pert`** from the new manifest before building.
5. **`n_real` is unmeasurable on D/E/F.** Decide the call-set policy from A/B/C evidence in
   October, not November. Design 1 above is how that evidence gets bought.
6. **Instruments are all we will have**: h1-test, PLOLO, 096's claim route. 187 adds one rule —
   the blind per-context knee is expected to gain ≈ 0, and the expected-value optimum under 1%
   cliff uncertainty is **knee × (1 − 1.25%)**, worth +0.0009 over a global 0.982.
7. **Submit the best-supported package early in the dark window**, then replace it only for
   something strictly better on the local instruments with a preregistered prediction. The last
   entry counts and a late unverified change is the one way to lose the round.

### Standing allocation

Direction work is a lottery ticket with a large payoff; ordering work is a near-sure +0.008 to
+0.014. **Run both, but never let a direction experiment block an ordering one** — and when an
ordering arm clears the gate, ship it rather than waiting for the direction bet to land.

## 2026-09-03 03:15 UTC — 186: THE FIELD'S `reach` SPREAD IS ORDERING, NOT DIRECTION. 154 IS OVERTURNED BY THE BOARD.

`experiments/186-intel-0903/RESULTS.md`. 612 teams, every model name and description mined, plus
the PertMind paper, D²R², Lingshu-Cell and the GitHub/arXiv sweep since 2026-08-28.

**Central question — what buys raw `reach` ≥ 0.25 at raw `pds` ≥ 0.84 — answered four ways:**

1. **Sign accuracy is INVERTED across teams.** Inverting `nmae` through `nmae(k, q)`, `Happy`
   (reach 0.1871) and `Huiyan Qian` (0.1857) have strictly better `q` than `Slurm Shadi` (0.2584),
   `Haystack` (0.2494), `Andrew Liu`, `jw_13288`, `Hongjian Zhou` and `cqawesome` — at **every**
   common `k` — and lose 0.06–0.07 raw reach. **15 of 49 high×low pairs sign-invert**; a `q` story
   needs zero. Only 18% of `cqawesome`'s +0.0498 reach jump is the `q` channel.
2. **Six same-team arms moved `reach` with the profile pinned** (|Δpds| ≤ 0.006, |Δlog mse| ≤ 0.12):
   `luka` +0.0374, `Weijing Zhao` +0.0324, `tydymy` +0.0353, `VIVAI` +0.0297, `ggband` +0.0253,
   `BioFarmer` +0.0455 — four with `nmae` getting *worse*.
3. **Amplitude is weak**: over 29 amplitude arms `d reach = +0.0127 · d log(mse)`, residual sd
   **0.0363**, i.e. 2.9× the effect. Six teams gained reach while halving `mse`.
4. **Architecture is worth zero**: three teams swapped STACK→Lingshu-Cell, CellFlow→Lingshu-Cell and
   CPA→STATE; every member moved < 0.006 and all landed at reach ≈ 0.208.

**This overturns 154's "reach is model-limited, not order-limited"** on evidence 154 could not
have had, and it confirms 171/172's delivery diagnosis from outside.

**Our deficit, priced against the field's conditional expectation.** On the 41 coverage-saturated
teams (`fid` ∈ [0.503, 0.512], `pds` > 0.79), `reach ~ pds + nmae + log mse + jac` gives R² 0.582
and our residual is **−0.0411 raw = −0.0077 overall, −1.92σ, 4th-worst of 41.** `corr(reach, jac)`
= −0.569 is the strongest correlate, but trading along that line prices at **+0.0004 overall per
0.005 `jac` — nothing.** The prize is the residual, not the trade.

### The lead: `tydymy` published a board-measured ordering arm, and our own 108 used it wrong

`tydymy` (Complex Disorders Alliance) put their config on the public leaderboard: Orion HCT116 +
HEK293T + Replogle K562, decomposition transfer, `conserved_alpha 1.30`, **"plus a CollecTRI signed
TF to target correction at 25 percent of backbone L2."** The arm that introduced it moved, on the
board: `reach` 0.2196 → 0.2549 (**+0.0353 = +0.0066 overall**), `fid` 0.4789 → 0.4922 (**+0.0133 =
+0.0075 overall**), `jac` +0.0030, `pds` −0.0013, `nmae` 0.9792 → 1.0213 (**−0.0421 = −0.0116**).
Net +0.0021.

**The prior bought the ORDER and the CALLS and paid in the MAGNITUDES.** Corroborated the same week
by **D²R² (arXiv 2608.15288)**, which initialises an ordering policy from a GRN **inferred from
control cells** and ablates: biological-prior ordering beats random, beats uncertainty heuristics,
and degrades every metric when reversed — on Norman19 and **VCC-H1**.

**Our experiment 108 tested a signed regulon as a DIRECTION blend and closed NO-GO. That is the use
that costs `nmae`. As an ordering key it has never been tested here.** And because `reach` and `fid`
read only the order and the call set, delivering the key at fixed pseudobulk collects the +0.0141
gross with `pds`, `mse` and `nmae` provably unchanged — `tydymy`'s −0.0116 bill is never incurred.
**Experiment 190 owns this.**

### Three corrections and one warning from the same sweep

- **PertMind is not why `luka` leads on `reach`.** `luka`/`PertMind_Pro` and `UM`/`BoomLi` uploaded
  the **same file** on 2026-08-28 — identical to 16 significant figures on all six members — and
  `luka`'s roster contains `BoomLi`. Their reach came from a 2026-09-02 ordering change shared
  across ≥ 4 accounts. Its *label rule* is reusable (Apache-2.0): call at `|lfc| ≥ 0.585` and
  `padj ≤ 0.05`; "No" at `|lfc| < 0.20`, `padj ≥ 0.50`, `baseMean ≥ 10`.
- **A direct tension with 177**: `Qigong Li` measured through the **official scorer** that a
  parametric NB with **per-gene dispersion fitted from the context's own NTC cells** beat bootstrap
  resampling on 5 of 6 members (`reach` +0.003). Resolve before any emitter arm ships.
- **`lingshu` = Lingshu-Cell**, arXiv 2603.25240, Alibaba DAMO.
- **`Haystack` (`SE`, `GSA-020`) is the single best team to reverse-engineer**: the board's largest
  ordering surplus, **+2.02σ**, at `pds` 0.8014 — *worse than ours* — with `reach` 0.2494. Twenty-one
  submissions, zero descriptions, no public trace.

## 2026-09-03 03:00 UTC — 189 and 190 launched on the ordering channel

**189 — the per-gene cross-line reliability prior.** 172 closed learned confidence keys, but its
feature set was explicitly *"no gene identity and no target identity"*, so a prior attached to the
**gene** has never been fitted here. Coordinator pilot on the three-way donor intersection (7,976
common genes × 9,073 shared targets, each matrix centred across targets):

| quantity | value |
|---|---|
| per-gene `r`(HCT116, HEK293T) | mean 0.5379, sd 0.0467, p5 0.499, p95 0.594 |
| spearman( r(HCT,HEK), r(HCT,K562) ) — a **held-out** line | **+0.3833** |
| spearman( r(HCT,HEK), r(HEK,K562) ) | +0.3359 |
| decile spread of r(HCT,HEK) read on held-out r(HCT,K562) | **+0.0402** |
| spearman( r_HE, control expression ) | +0.5215 |
| after partialling out expression, \|effect\| and across-target sd | **+0.2133 / +0.2251** |
| held-out decile spread of the **residual** prior | +0.0173 |
| held-out decile spread of **expression alone** | **+0.0505** |

Two findings. A per-gene reliability prior exists, generalises to a held-out line, and survives
every observable confounder. And **control expression is by itself the strongest per-gene
reliability signal available** — which sits in direct tension with 171, where our p-order being
**expression-set** was diagnosed as the defect. Part of what 171 called a defect is signal, and the
right key is a composite. Resolving that decides what 185 should deliver.

**190 — the regulatory-prior ordering key**, on the `tydymy` + D²R² evidence above.

## 2026-09-03 02:40 UTC — 187 CLOSES GLOBAL AMPLITUDE, AND CORRECTS THE `nmae` READING BELOW

**Read this before the 01:00 entry — it falsifies one of that entry's conclusions.**

`experiments/187-joint-amplitude/RESULTS.md`. Four scripts, 82 turns, no build, no `cell-eval2`.

**1. The joint `mse` + `nmae` optimum for a global α is `0.982` to four decimals — exactly where
a15/a23 already sit.** Any increase crosses context A's ρ cliff at α 0.9823, where `mse` loses
**30×** what `nmae` gains per unit `ln α`; any decrease loses on both. **Global amplitude is
CLOSED.**

**2. abc-knee is retrodicted with zero free parameters at the point of prediction**, which is what
the 00:00 entry said had to happen before the knee could be reopened. The C_p-free solve
(`X` 0.019240, `ΣD` 0.229141) reproduces all four of the a12 audit's pre-registered proxy
hypotheses to ≤ 5e-5, and the two board deltas that straddle the cliffs then discriminate them:

| hypothesis | a12→a14 err | abc-knee err (scaled) |
|---|---:|---:|
| NTC proxy as 180 used it | +0.00003 | +0.0073 |
| h1 level bias ×1.0195 | +0.0080 | +0.0073 |
| h1 spread CV 0.052 | −0.0055 | +0.0004 |
| **per-context factor `f = (1, 1, 0.9825)`** | **+0.00003** | **+0.00018** |

The reference's jackknife on **context C** is 1.5–2.0% below the NTC-construct proxy while A and B
sit on it. abc-knee overshot C's true cliff by 0.9% and paid +0.0103 raw there — that is the entire
`mse` miss. The h1-measured level bias and spread are both **falsified on the board**.

**3. Per-context α is worth +0.0011 vs a23 and +0.0004 vs abc-knee — composed gate +0.00075 against
the +0.0089 amplitude floor. VETO by 12×**, and the builder cannot express it anyway. For D/E/F the
blind knee is expected to gain ≈ 0; the expected-value optimum under 1% cliff uncertainty is
**knee × (1 − 1.25%)**, worth +0.0009 over a global 0.982 with the asymmetric downside removed.

### 🔴 CORRECTION to the 01:00 entry — `nmae` is a DIRECTION problem, not an amplitude one

The 01:00 entry read our board `nmae` 0.9432 through the uniform-k identity as **`k ≈ 0.75`** and
concluded *"the leader's entire `nmae` advantage over us is amplitude, not direction, and it is
available today"*. **That is wrong, and 187 shows why.**

The uniform-k identity collapses a per-gene `k` that is **not** uniform: on h1-test the |t|-weighted
per-gene `k` has **median 0.12 and p99 1.11**, and the board's own measured amplitude response
(−0.0416 raw per `ln α`, from two agreeing current-pipeline pairs) is only **75% of the uniform-k
slope**. The exact gene-wise form is

```
nmae = 1 − 2 q_w W_hi − (2 q_lo − 1) K_lo + K_hi
```

which collapses to `q|k−1| + (1−q)(k+1)` only under uniform `k` and gene-independent `q`.

**Consequence:** the collectable `nmae` gain inside the `mse`-feasible amplitude window is
**0.0012 raw**, and abc-knee already collected 0.00116 of it. The gap to the leader is 0.0163 raw
and to Aginglab 0.054. So **≈ 7% of our `nmae` deficit is calibration and ≈ 93% is direction** — the
opposite of the 01:00 reading. The credit that moves `nmae` is `⟨sign(t), p⟩` on the confident genes
(|p|-weighted sign accuracy 0.62 on h1 against count accuracy 0.55), and raising it is 184's
territory, not a scalar's.

A uniform `k` increase is safe on `nmae` to about **+20%** (turnover at +29%); the `mse` cliff binds
an order of magnitude earlier, at **+0% / +3.1% / +5.8%** for A/B/C.

### What this leaves, stated plainly

With global amplitude closed and `nmae` reclassified, the members split cleanly:

| channel | members | status |
|---|---|---|
| **direction** (`q`, cosine) | `pds`, `mse`, `nmae`, and the sign half of `fid`/`reach` | information-limited at 3 donor lines; 184 (Tahoe) and 188 (magnitude shape) are the only live attacks |
| **delivery** (the emitted p-order and call set) | the ordering half of `reach`, `fid`, `jac` | **the one channel that is not direction-bound**, and the one where the board says the defect is ours alone (`docs/BOARD-LEDGER.md`). 185 owns it |

**So the cheap-wins estimate in the 01:00 entry drops from ≈ +0.019 to ≈ +0.006–0.016, all of it in
the delivery channel.** Everything else now requires a better direction.

## 2026-09-03 01:00 UTC — THE BOARD'S OWN CEILING IS 0.2546, WE ARE INSIDE THE PARETO FRONTIER, AND THE GAP IS NOT `pds`

Fresh pull, 612 scored teams (`data/leaderboard/` snapshot). We are **rank 22 at 0.17404**;
leader `cqawesome` **0.21790**. Five experiments launched on this analysis (184–188, below).

### The scoring map, re-solved on 612 entries

Per-unit-raw slopes, fitted on the live board (residual rms in brackets — `pds`, `reach`, `jac`,
`fid` fit to 4 decimals, so these are the scorer's own anchors, not estimates):

| member | dScaled/dRaw | dOverall/dRaw | implied baseline `b` | implied replicate `r` |
|---|---:|---:|---:|---:|
| `pds` | +2.2138 | +0.3690 | 0.5002 | 0.9519 |
| `mse` | −1.051 | −0.1752 | ≈0.988 | ≈0.037 |
| `nmae` | −1.6532 | −0.2755 | 1.0003 | 0.3954 |
| `fid` | +3.3645 | +0.5607 | 0.5124 | 0.8097 |
| `reach` | +1.1292 | +0.1882 | 0.0790 | 0.9646 |
| `jac` | +2.7263 | +0.4544 | 0.0304 | 0.3972 |

### The gap to first, decomposed — `mse` and `reach` are 81% of it

| member | ours raw | leader raw | gap in overall points |
|---|---:|---:|---:|
| `mse` | 0.8745 | 0.7604 | **+0.0200** |
| `reach` | 0.1773 | 0.2592 | **+0.0154** |
| `pds` | 0.8368 | 0.8684 | +0.0117 |
| `nmae` | 0.9432 | 0.9269 | +0.0045 |
| `fid` | 0.5038 | 0.4929 | −0.0061 (we win) |
| `jac` | 0.0301 | 0.0274 | −0.0012 (we win) |
| | | | **+0.0438** |

### Conditional on our own `pds`, we are below the field line on exactly three members

Regress each raw member on raw `pds` over the 206 teams with `pds > 0.70` (`mse` over the 49 with
raw `mse < 1.2`), then read our residual:

| member | ours | field fit at our `pds` | residual | in σ | overall pts to the line |
|---|---:|---:|---:|---:|---:|
| **`reach`** | 0.1773 | 0.2206 | −0.0433 | **−1.04σ** | **+0.0082** |
| `mse` | 0.8745 | 0.8392 | +0.0353 | +0.59σ | +0.0062 |
| `nmae` | 0.9432 | 0.9291 | +0.0141 | +0.43σ | +0.0039 |
| `fid` | 0.5038 | 0.5013 | +0.0025 | +0.06σ | — |
| `jac` | 0.0301 | 0.0298 | +0.0003 | +0.08σ | — |

**+0.0183 sits in `reach` + `mse` + `nmae` with no improvement to the direction at all**, and the
field proves each is attainable at our discrimination.

### The `mse` ↔ `reach` trade is NOT forced, and we are strictly dominated

Partial correlation of scaled `mse` and scaled `reach` given `pds`, among the 51 teams with raw
`pds > 0.80`: **−0.139**. Weak. The Pareto frontier in that plane:

| team | scaled mse | scaled reach | raw pds | overall |
|---|---:|---:|---:|---:|
| Andrew Liu | +0.2407 | +0.1915 | 0.8329 | 0.20834 |
| **cqawesome** | +0.2394 | +0.2034 | 0.8684 | 0.21790 |
| Aginglab.com | +0.1693 | +0.2066 | 0.8199 | 0.21122 |
| Weijing Zhao | +0.0000 | +0.2394 | 0.8194 | 0.17180 |
| **us** | **+0.1193** | **+0.1109** | 0.8368 | 0.17404 |

We are **inside** the frontier — dominated on both axes by two teams at once.

### The decisive comparison is Aginglab, not the leader

`Aginglab.com` beats us by **+0.0378** while having a **worse** `pds` than we do (0.8199 vs 0.8368):

| member | their edge over us, overall pts |
|---|---:|
| `reach` | +0.0160 |
| `nmae` | +0.0149 |
| `mse` | +0.0083 |
| `fid` | +0.0054 |
| `pds` | −0.0062 |
| `jac` | −0.0006 |

**Three of their four gains — `nmae`, `reach`, `fid` — are the sign-accuracy-bound members.** Put
their `nmae` through the 02:15 identity: at `k = 1`, `nmae = 2(1 − q)`, so their raw 0.8892 solves
to **`q` = 0.5554**, and the leader's 0.9269 solves to **`q` = 0.5365 — our own `q` at `k = 1`**.
Our 0.9432 at `q` = 0.537 solves to **`k` ≈ 0.75**.

**So the leader's entire `nmae` advantage over us is amplitude, not direction, and it is available
today.** Aginglab's further edge is a genuinely better `q` (+0.018 absolute). That single arithmetic
line splits the `nmae` gap into a cheap half we own and an expensive half that is 184's territory.

### The empirical ceiling: 0.2546

Take the field-best scaled value of each member and average: **0.25463**. That is what "as good as
the best human on every member" is worth, and no single team is near it (best is 0.2179). Ahmet's
0.25 target maps exactly onto it. Above that requires beating the field everywhere, which needs `q`.

| member | field best (scaled) | who | ours |
|---|---:|---|---:|
| `pds` | +0.8169 | cqawesome | +0.7477 |
| `mse` | +0.2407 | Andrew Liu | +0.1193 |
| `nmae` | +0.1838 | Aginglab.com | +0.0962 |
| `reach` | +0.2394 | Weijing Zhao | +0.1109 |
| `jac` | +0.0242 | Swamy Vijayan | −0.0004 |
| `fid` | +0.0227 | The_Way | −0.0295 |

### The five experiments this launched (2026-09-03 00:40–01:00 UTC)

| # | owns | question |
|---|---|---|
| **184** | Tahoe-100M | 1,100 drugs × ~50 cell lines is the only substrate that can estimate a **cross-line transfer operator**. 159's 1/k law says we need ~10 lines and 166 says no more CRISPRi lines exist — but all five failed line-adaptation experiments had ≤4 lines to fit on. Phase 2 is a leave-one-LINE-out learning curve in `k` inside Tahoe: if the advantage is flat in `k`, the line-specific 70% is closed for good. **Never touched by this project before today.** |
| **185** | the emission law | 171 closed the *spread-schedule* family, but its knob was one scalar inside a fixed parametric law. At fixed per-cell total and fixed per-gene total the pseudobulk members are **pinned** while the rank-reading members are **free** — a rank-sum statistic reads the count histogram's shape and the mean does not. Tests whether the emitted p-order is a designable object, which would deliver the `|H|` key (0.18283 → 0.23533 raw reach). |
| **186** | the outside world | What mechanism gives raw `reach` ≥ 0.25 at raw `pds` ≥ 0.84? Four teams do it. Also: `PertMind` (team `luka`) holds the board's best raw reach 0.2908 and is public — RL over **Tahoe-100M**. |
| **187** | the global α | The joint `mse` × `nmae` optimum, which nobody has solved together. Must retrodict `abc-knee`'s −0.00226 / +0.00208 or admit it is the same broken proxy. |
| **188** | the per-gene magnitude SHAPE | Resurrects 169. Its `pds` gain (**+0.019452** PLOLO, 2.9× a15, the only arm ever to clear Gate 4) died on `nmae` alone, and the identity says why: a power law charges slope **+1** on the 9% of genes already at `|h| > |r|`. `nmae` is an **L1** loss, so its optimal magnitude map is the **conditional median of true given predicted**, fit LOLO — not a monomial. |

## 2026-09-03 00:00 UTC — `abc-knee` published: overall near the prereg, `mse` FALSIFIED

Entry `N3ui4hCkKx9bxExG2pxE`, joke name, empty description. Full write-up:
`experiments/180-def-readiness/BOARD-RESULT.md`.

| | prereg Δ scaled vs a23 | board Δ scaled | |
|---|---:|---:|---|
| pds | 0 ± 0.001 | +0.00012 | holds (scale-invariant) |
| **mse** | **+0.0050** | **−0.00226** | **falsifier fired** |
| nmae | ≈ +0.0002 | +0.00208 | sign right |
| fid / jac | ≈ 0 | +0.0006 / +0.0001 | noise |
| reach | ≈ 0 | +0.00347 | inside seed |
| **overall** | +0.0008 to +0.0012 | **+0.00068** | right size, **wrong member** |

Raw `mse` 0.87235 → 0.87449 (worse). We raised B/C toward a locally computed knee and the real
scorer punished amplitude. The +0.00068 overall is `nmae`+`reach` covering a `mse` miss.

**D/E/F default reverts to a global alpha (a15 / 0.982). Do not ship `VCC_ALPHA=knee` on Oct 22.**
The `C_p` NTC-jackknife proxy — or 096's claim route against the board functional — does not
transfer. Reopen knee only with a proxy that retrodicts this point.

a15 remains best (knee is −0.000043 vs a15). `hy_kt8000` was not submitted.

## 2026-09-03 03:00 — ADVERSARIAL AUDIT (exp 181) OF EVERY 2026-09-02 CLOSURE. Read this before anything below it.

A fresh agent re-derived every number of the day from on-disk artifacts (16 known-answer rows,
all exact) and tried to break each closure. **Every NO-GO holds. The day's one recorded WIN is
refuted by a board point that was on disk the whole time.**

**🔴 CONTRADICTED — "the ordering defect is solved by the emission law" (177) and `hy_kt8000` (179).**
**a21 (2026-09-01, `resample_x 0.14`) IS this mechanism**: keep the top ~2,600 genes by `|effect|`
parametric and resample the rest — `hy_kt≈2600`, exactly the dose where h1-test shows the largest
reach gain (+0.050 scaled predicted). **The board scored it: reach −0.0185, nmae −0.064, fid
−0.029, jac −0.017.** h1-test gets 2 of 4 DE members wrong in sign for this emission law. Both
call-site contracts 177 named were honoured at a21's commit, so a21 is a valid board test of the
law. **Rule 13's −0.0146 reach toll is exactly `mean(a21, a22)` — two exact board measurements of
this very emission law — and it was right.** 157's "validation" of h1-test for reach was a8 vs a15,
a spread-*policy* change inside the quasi law, not an emission-law change; a21/a22 are the only
emission-law board points and both invert h1-test. Corrected `hy_kt8000` expectation:
**−0.009 to −0.004 overall. DO NOT SUBMIT. The `|H|` counterfactual (+0.0525 raw) remains real
and UNCLAIMED — no known delivery mechanism survives the board.** Sixth board measurement that
shrinking the call set has never raised board `fid` (a21: n_pred 4,420 → 1,141, fid 0.5034 → 0.4945).

**Also CONTRADICTED / WEAK:** `corr(raw_reach, raw_pds)` is **+0.701**, not +0.831 (the −1.22σ
residual argument stands); 179's `fid` +0.0029 branch (board series says −0.0005); 093's `n_real`
2,670 (now closed independently: board reach would read 0.0034); the 0.46 DE-delta scaling
(double-counts a shrink); **"production sits at the pds-optimal point; only new information can
help" is WEAK** — every input to it is an un-whitening move, while 170's delivered `pds` peaks at
gamma 0.50–0.55 and 169's sharpening raised `pds` +0.0195 the same day; **"h1-test understates
`pds` costs 6×" is a frame effect, not instrument bias** — read `pds` in the challenge frame (PLOLO);
175's "Orion adds nothing" was measured on the 64% of the pool K562 reaches; **`n_real` is a
distribution, not a scalar** — board reach ⇒ median 90–110, board jac ⇒ heavy tail ~half of h1-test's.

**SOLID (strengthened):** the ρ-cliff *algebra* (board-validated to ±1% on a15's global alpha);
the long-gene technical axis (**gene-length correlation −0.71 / −0.65 over 10,601
genes; K562 control −0.02** — no longer an eyeball); the board partition (exact); 171's BH no-op;
164's tail precision; 176's CMap NO-GO (understated: landmarks are 9.5% of the pool).
**The per-context knee as a D/E/F default is now NO-GO** — abc-knee's board Δ`mse` was −0.00226
against a local +0.005 (2026-09-03 00:00 UTC). The algebra is right; the `C_p` proxy is not.

**Standing rules the audit adds (also in AGENTS.md):**
1. **A screen whose family already has a board point must retrodict it, and the known-answer row
   must be the board's, not the harness's.** 177 and 179 carried eighteen harness-internal gates and
   zero board retrodictions while a21 sat on disk.
2. **Rule 13's reach toll is two constants**: −0.0146 for emission-law arms (a21+a22, exact, correct)
   and ≈ 0 for profile-class arms (169 falsified it there).
3. **Reopen `pds` with one lever that is NOT an un-whitening move** — 170's criterion, norm-preserving
   and raising `X` without `A`, is unexplored outside sign repair.

`experiments/181-audit-0902/RESULTS.md`.

## 2026-09-03 02:15 — the `nmae` landscape in one formula, and the design constraint it imposes

`nmae` is the normalised MAE of log2FC over the reference's DE genes; predicting zero scores 1.000.
For a DE gene with true magnitude `|t|`, our magnitude `|p| = k·|t|`, and sign accuracy `q`:

```
nmae(k, q) = q·|k − 1| + (1 − q)·(k + 1)         (predict-zero = 1.000 for every q)
```

| q | k = 0.5 | k = 0.75 | **k = 1** | k = 1.25 |
|---:|---:|---:|---:|---:|
| 0.537 (ours) | 0.963 | **0.944 ← our board** | 0.926 | 1.176 |
| 0.60 | 0.900 | 0.850 | 0.800 | 1.050 |
| 0.70 | 0.800 | 0.700 | 0.600 | 0.850 |

**The optimum is k = 1 exactly, worth `2(1 − q)`, and the slopes are asymmetric: −(2q − 1) = −0.074
below it, +1 above it.** Undershooting by 25% costs 0.018 raw; overshooting by 25% costs 0.25. This
is why every magnitude-inflating lever has died on `nmae` — sharpening (169), a13, a21's library
unwelding — each pushed some DE-pool genes past k = 1, where the 46% wrong-sign genes are charged
`|p| + |t|`. Our board 0.944 reads as **k ≈ 0.75 at q = 0.537**: we sit on the cheap side of the kink.
The leader's 0.927 is k ≈ 1 at our q, or q ≈ 0.56 at lower k; AgingLab's ~0.89 is q ≈ 0.555 at k = 1.

**Design constraint for the reopened `pds` question (audit 181 item 4):** any norm-preserving
transform that raises `X` must not push DE-pool magnitudes above their true values — magnitude may be
*redistributed* only from non-pool to pool genes and only up to k = 1, which we cannot observe per
gene. A safe move is a small uniform k increase (the per-context knee alpha already does +3–7% on B/C
and 180 measured `nmae` +0.0002); a large one collides with the `mse` knee. **Beyond that, `nmae` —
like `reach` and `fid` — moves only through q.** Three of six members are now shown to be
sign-accuracy-bound by formula, not just by measurement.

## D/E/F TRANSITION PLAN — written 2026-09-02, must be executed BEFORE Oct 22

The final panel (contexts D/E/F, new 300 targets) is released **Oct 22**; the deadline is **Nov 5
23:59 UTC**; the leaderboard is **dark** for the whole window. Every final-round submission is blind
and only the last one counts. Nothing below can be done under time pressure on Oct 22.

> **STATUS 2026-09-03 01:30 — items 1 and 2 DONE by exp 180 (`experiments/180-def-readiness/`).**
> The builder reads contexts from the manifest, the frozen-input table lives in
> `scripts/frozen_inputs.json` with `--repin-frozen-inputs` / `--repin-panel-support`, and a
> synthetic D/E/F bundle was rehearsed end to end through `vcc prep --dry-run` — delivering
> **byte-identical `X/data`/`X/indices`/`X/indptr` to a23's** with only the labels changed. A/B/C is
> bit-identical under the patch (179 effects fields, one difference = the builder's own hash).
> **`VCC_ALPHA=knee` validates locally and FAILED on the board (abc-knee, 23:28 UTC).** Local
> alphas 0.98586 / 1.01275 / 1.04568, budget ratio ≈1, PLOLO `pds` invariant. Board Δ`mse`
> **−0.00226** vs local +0.005. **D/E/F default is a global alpha (a15 / 0.982), not knee.**
> Package submitted as entry `N3ui4hCkKx9bxExG2pxE`. **`vcc prep` hardcodes `REQUIRED_CONTEXTS =
> ("A","B","C")` — pass `--contexts D,E,F`.** Full runbook in 180 §10.

**1. The production builder refuses non-A/B/C contexts.** *(DONE — see status above.)* `scripts/make_submission_pm_cd4.py:2617`
raises unless `spec.contexts == ("A","B","C")`, and the frozen-input hash table (lines ~110–130 and
~1258) hardcodes `context_A/B/C.h5ad` and their sha256s. This is a deliberate validation adapter. It
must be generalised to read contexts from `manifest.json` (the `ChallengeSpec` already does), with a
**re-pin step** for the D/E/F bundle hashes, and **rehearsed on a synthetic D/E/F bundle** (A/B/C
renamed, new manifest) end to end through `vcc prep --dry-run`. Keep the A/B/C build bit-identical
under the change (rule 10).

**2. `VCC_ALPHA` is a single global constant fitted on A/B/C board points, and the `mse` optimum is
context-specific.** *(IMPLEMENTED — then FALSIFIED on the board as `abc-knee`. Keep the code;
do not make it the D/E/F default. See `experiments/180-def-readiness/BOARD-RESULT.md`.)* 170 established `mse` is minimised at the ρ knee `s* = √(claim/B)`, and a15's
0.982 sits on it *for A/B/C*. `claim` depends on the context's `C_p`, which 096 proxies from the
**context's own NTC control cells** (delete-1 jackknife of one construct's 400 cells, stable to ~3%);
`B` is our own prediction's across-perturbation budget. **So the knee is computable per context from
the D/E/F controls alone — no board needed.** To do: implement `VCC_ALPHA=knee` (per-context `s*`,
default off), validate on A/B/C that it reproduces ~0.98–1.01 and the board `mse`, then make it the
D/E/F default. Off-knee in either direction costs `mse` (170's V). **Validated arithmetically (23:45):** at a15's
alpha the contexts sit at ρ = 0.992 / 0.940 / 0.882, i.e. A on the knee, B 3% and C 6.6% below in
amplitude; the global-alpha optimum is where the contexts *straddle* the knee (170's "ctx at ρ=1"
column). Per-context knee alphas are **0.987 / 1.013 / 1.047**, worth ≈ +0.0008 overall on A/B/C
(small because `X` is small; `pds` invariant to a within-context scale). For D/E/F the stakes are
larger: a context landing at ρ 0.70 under a global 0.982 forfeits ≈ 0.03 scaled `mse` that the
knee rule recovers from its own controls.

**2b. IDENTIFY D/E/F WITH THE CCLE MATCHER ON DAY ONE — this is a new transition item and it is
time-critical.** Experiment 041 identified the validation contexts from their control profiles
alone against the DepMap 25Q3 CCLE mirror, with the same top match under all eight
feature/rank variants: **A = JURKAT (T-ALL), B = HELA (cervix), C = CAL33 (head-and-neck
squamous)**. The matcher was calibrated on three known controls first (K562 and HCT116 rank 1;
independent Jurkat rank 1–3). That identification is what motivated the **CD4 residual now
shipping in a15** — primary CD4 T-cells as a biologically matched auxiliary donor for A/Jurkat.

So the line identities are not trivia; they select which auxiliary donors are matched, and the
one time we acted on them it produced a shipped component. **Run `experiments/041-ccle-match/`
on the D/E/F controls the hour they land**, before any build, and re-audit the auxiliary-donor
set against the new lineages. If any of D/E/F is a line with public Perturb-seq of its own
(HepG2, Jurkat, RPE1, A549, MCF7, HT29, HAP1, BxPC3 are all covered by data we already hold or
can fetch), that is in-domain information we do not have for A/B/C.

**Boundary, from 041 and the rules (§7):** public data may **train or condition** a model;
it may not be used to hand-revise output. A learned per-line adapter fit on that line's own
public perturbation data is training and is allowed; a lookup of a measured effect is not.
Note the panel deliberately avoids essential genes, so the essential-gene screens
(`nadig_jurkat`, `rpe1`, `k562_essential`, `nadig_hepg2`) cover **0/300** of the current panel —
they can train an adapter but can never serve as a donor at inference.

**3. Re-audit donor coverage for the new panel.** K562 covers 269–272/300 of the *current* panel;
Orion 290–294. D/E/F's 300 targets are different. Run 166's `census.py` on the new `pert_counts.csv`
the day it lands; if coverage drops materially the equal-weight consensus changes character.

**4. Verify the gene axis** (`gene_names.csv`) and `cells_per_pert` from the new manifest before
building; the writer asserts both.

**5. No D/E/F `n_real` can be measured.** The call-set regime that decides the resample-tail arm
(177) must be judged from A/B/C evidence — either the board reads from an information slot spent
before Oct 22, or the default stays a15's law. Decide this in October, not November.

**6. Local instruments are all we will have.** h1-test (emission-law `reach`, `fid`, `nmae`), PLOLO
(`pds`), 096's claim route (`mse`). Their validation record (162, 165, 168d, 170, 178) is the entire
basis for every blind decision; do not ship anything in the dark that they cannot price.

**7. Submission discipline in the dark:** submit the best-supported package early in the window, then
re-submit only for a package that is strictly better on the local instruments with a preregistered
prediction; the last entry counts, and a late unverified change is the one way to lose the round.

## 2026-09-02 — a23 in flight: the renderer was undoing the geometry

`profile_from_effect` divides the rendered profile by its own sum, and that division is the
identity only if `sum_g clip(expm1(b_ctrl + h_g), 0) == TS`. The control profile and every
reference perturbation satisfy it exactly; **our effect does not** (mass ratio `1.00134` on the
challenge frame), so the renderer applies its own correction — **2.83–3.35% of the row norm on
A/B/C, at row-row cosine 0.303 against an intended row-row cosine of 0.0040.** The contamination
is 75x more mutually-aligned than the signal, and hubness is exactly what `pds` punishes.

`VCC_MASS_FIX=minnorm` imposes the constraint with the minimum-L2 correction instead. All six
members move in our favour on both held Orion lines; aggregate cosine `+13.6% / +56.0%`; sign
accuracy on the reference DE pool `0.5311 -> 0.5417`; of 318 flipped signs, **68.9% land
correct**. Prediction `+0.004` (range `+0.003` to `+0.011`), pre-registered with falsifiers in
`experiments/156-delivered-geometry/A23-PREDICTION.md`. Full record: that experiment's
`RESULTS.md` and the 2026-09-02 entry in `docs/LOG.md`.

**Two consequences larger than the arm:**

1. **Every effect-level screen since experiment 017 scored an object the board never sees.**
   `identity/intended` beats `identity/delivered` by `+0.00324 / +0.00382` raw pds. Harness
   geometry and board geometry have differed all along. Any future geometry screen must render.
2. **Per-target amplitude is a confidence problem, not a magnitude problem.** `pds` is exactly
   invariant to a per-row scale, so it is a free `mse` lever. Oracle per-target scaling is worth
   `+59% / +95%` aggregate cosine — but oracle MAGNITUDE alone (true `||delta_p||`, no per-target
   cosine) is worth `+0.000011`. **18.9% / 27.3% of targets have an optimal scale <= 0**: our
   direction is anti-correlated with the truth for a quarter of targets, and identifying them is
   the whole lever. This refutes the COMPASS-derived "magnitude is the transferable half, so
   monetise it" reading in `docs/COMPETITIVE-INTEL-2026-08-28.md` §2.5.

**Cross-line transfer, measured in the scorer's own space** (120 targets shared by the two Orion
lines): same target in a different line, cosine `+0.1881`; different targets across lines,
`+0.1790`; **target-specific excess `+0.0091`.** Energy split 60.1% shared / 39.9% line-specific;
a perfect context-independent predictor caps at cosine 0.77 to each line. Reproduces the
independent `+0.021` sign-accuracy excess from 2026-09-01. Attenuated by 150-cell noise;
experiment 159 is disattenuating it and converting it into a paradigm ceiling.

**Also corrected:** `AGENTS.md` §4.2 point 1's stated reason for banning panel-fitted operators
("orthogonalises 259 vectors, which improves a ranking by construction") is **wrong on the
metric's algebra** — `k_p` ranks the 300 *measured* effects against one fixed prediction, so
decorrelating predictions cannot improve a ranking by construction. The real hazard is a
rank-deficient covariance. The lever is nonetheless empty for us: our panel matrix carries 0.35%
shared energy and `centroid` pre-render measures exactly `+0.00000`.

**a15 published 2026-08-30 20:14 UTC** under the public model name
`The Fern Filed An Appeal`. It improved a14 by `+0.00366226`; all scaled
members improved except effectively-flat Jaccard. Raw `mse` improved
`0.882775 -> 0.873805`, PDS improved `0.832702 -> 0.835987`, and raw nmae
improved `0.945764 -> 0.944251`. Full audit:
`experiments/096-mse-floor/A15-BOARD-RESULT.md`.

## The structural result — there is ONE problem, not four

`de_wilcoxon_lfc_nmae` is a **direction member, not a magnitude member**, in our regime. Its
denominator and gene gate are reference-only, and our `|lfc_hat|` is 0.14-0.17x the reference's
on its own gate with `|h| <= |r|` on 99.3% of gated genes — so

> `NMAE = 1 - <sign(lfc_ref), lfc_hat>`   (residual < 0.004)

| member | leader | us (a15) | gap |
| --- | ---: | ---: | ---: |
| `nmae` | 0.1838 | 0.0942 | +0.0895 |
| `reach` | 0.2066 | 0.1154 | +0.0912 |
| `fid` | 0.0030 | −0.0304 | +0.0334 |
| `mse` | 0.1693 | 0.1201 | +0.0493 |
| `pds` | 0.7089 | **0.7459** | −0.0369 |
| `jac` | −0.0042 | **−0.0006** | −0.0036 |

**`nmae` + `reach` + `fid` = 96.1% of the net six-member gap, and all three are per-gene
directional accuracy.** `mse` is the main remaining non-directional deficit; our PDS and
Jaccard advantages offset part of those four deficits.

## Current submission state

| arm | change | predicted | status |
| --- | --- | --- | --- |
| **a14** | a12 with `VCC_ALPHA` 1.0108 -> **0.982** — the true PER-CONTEXT cliff | `+0.0036`; observed `+0.00244` | published; previous best |
| **a15** | a14 plus `VCC_CD4_TARGET_F=0.10` — 099's `f`, **not** `CD4_FRACTION` (see below) | `+0.004`; observed **`+0.00366226`** | **published; current best** |
| **a19** | a15 plus whitener ridge, `gamma 0.65->0.80`, `lam 0.0->0.007` | `+0.002` (falsifier at `-0.001`) | **published; NO-GO `-0.01352`** |
| **a15 restore** | byte-identical re-submission of a15, `55f7bd41...` | n/a | restored after a19 |

**a19's board result is the most instructive of the project.** All five DE-side members moved
as predicted — `reach` `+0.00297`, `nmae` `+0.00845`, `fid` `+0.00400`, `jac` `+0.00057`, `pds`
`-0.00118` — and it still lost `-0.01352`, because **`mse` collapsed `0.120057 -> 0.024149`**.
The selecting screen measured `pds`, `reach` and `nmae` on two surfaces and never computed
`mse` at all. A whitener ridge changes the effect vector's shape and `mse` reads the profile
directly, so the single unmeasured member carried the entire result. See AGENTS.md 11.1 rule 12.

**Do not treat the `pds`/`mse`/`nmae` invariance as general.** It is a property of the quasi
emitter's exact-total repair and holds only for `spread` / `phi` / `keep_top` / call-set arms.
Anything touching the operator, the profile, or the emitted per-cell library sizes moves all
six.

a14 rests on a real bug found in our own model: the scorer applies `min(claim, budget)` **inside
each context**, and `min` is concave, so evaluating it at the panel mean is a different function.
Per-context cliffs are at alpha 0.982 / 1.008 / 1.042 and the curve turns up at the LOWEST.
Corrected, the model fits a8/a10/a11 to <=0.0001 and a12 to 0.0012 — but that last figure is
**in-sample**: the model was refit after a12 missed. a14 validated the central `mse`
prediction but not the simple nmae slope. a15 subsequently passed all four preregistered
falsifiers and is the current live best entry.

**a15's knob was mis-specified in this file and is now fixed.** An earlier revision said
"`VCC_CD4_FRACTION` 0.05 -> 0.10", which silently treats production's `CD4_FRACTION` as if it
were experiment 099's `f`. It is not: 062 rescales each condition to `EFFECT_SCALE` before
averaging and never rescales the mean back up, so the two knobs sit on different scales and
099 §5.4 put production's effective `f` at **0.033-0.049** — a *range*, because the
challenge-panel mixture norm had never been measured. 099 §8 asked for that mapping to be
re-solved **inside the builder's own frame** before anything was built, and named this exact
configuration — two codebases with a knob called `fraction` — as what cost this project a 5x
error on `gamma`. Shipping `0.10` flat would have landed `f` near 0.066, two-thirds of the
intended dose, against a cliff edge that 099 places **at** `f = 0.10`.

The builder now measures it instead of assuming it. `solve_cd4_fraction()` reports 099's `f`
for every build (`effects.cd4.residual_frame` in the meta), and `VCC_CD4_TARGET_F` solves the
fraction in closed form — the residual is linear in the fraction, so one measurement pins the
line exactly. Default unset: `CD4_FRACTION` is used untouched and every prior build reproduces
bit-for-bit. **Build a15 with `VCC_CD4_TARGET_F=0.10`, not with a hand-picked fraction**, and
read `effective_f_at_cd4_fraction` out of the meta to learn what the shipped 0.05 was really
worth.

**Built and measured, 2026-08-30.** `a15-cd4f010.vcc` is on disk, checksum-verified against the
worker, 36 tests green through the fail-closed gate. The shipped `CD4_FRACTION = 0.05` is worth
**`f = 0.0528`** — above 099 §5.4's whole estimated 0.033-0.049 range — and `f = 0.10` solves to
`CD4_FRACTION = 0.09473`, **below** §5.4's prescribed 0.10-0.15. Both of that section's bounds
were wrong, in opposite directions from the truth; the honest note is that this file's original
naive doubling to 0.10 would have landed `f = 0.1056`, nearer the target than anything §5.4
recommended, though 5.6% past a cliff edge rather than on it. a15 realises `f = 0.09997`.

## The ordering channel is closed — 2026-08-31

Experiment 154 screened every computable per-gene confidence key against the
scorer's own per-perturbation `reach`. Best buildable key `0.12473` against
`|pred_lfc|` `0.10716` and an oracle of `0.62797` — **3.4% of the headroom**,
median pure head 4 genes of 70 against the oracle's 44.

`reach` is **model-limited, not order-limited**. The emitter can express a good
order; there is no good order to express. Do not build ordering keys, per-gene
`phi` confidence maps, or spread policies aimed at `reach`. Donor concordance
alone is *anti-useful* (`0.05657`) and only helps multiplied by magnitude, and
"more donors" is unavailable — `k562_essential`, `nadig_jurkat` and `rpe1` cover
**zero** panel perturbations. Full record:
`experiments/154-confidence-key-screen/RESULTS.md`.

**Open question with a 10x spread: the challenge's `n_real`.** The `~2,670`
figure in `093-spread-truncation/BOARD-RESULT.md` is **circular** — it inverts
board `fid` using the *harness's* precision (`0.566`), and precision depends on
`n_real`, which A8's own rule says is unreadable on our harnesses. A competing
read puts `n_real ~= 240` from the board's `jac` level, with enrichment `1.335`
against the harness's `1.339`. The two worlds differ by 10x and by the entire
value of the call-set channel: at `~2,500` our call set would be 5.7x worse than
random at hitting the reference's DE set, implying roughly `+0.025` overall
sitting unclaimed; at `~240` the channel is at its ceiling and closed. A
fixed-profile, raise-`n_pred`-only submission separates them, with `pds`/`mse`/
`nmae` invariant by construction.

## Dual-moment emitter (148-152) — cross-context transfer has FAILED, 2026-08-31

One fixed lambda-1 dual-moment emitter was carried across four held surfaces.
It does not transfer, and the reason is now specific.

| exp | held surface | projected avg delta | outcome |
| --- | --- | ---: | --- |
| 148 | H1, 45 held targets | `+0.009707` | positive; authorized cross-context test only |
| 149 | HEK293T | `+0.012919` | **NO-GO** on frozen mechanism gates despite the gain |
| 150 | HCT116, uniform lambda 1 | `-0.005045` | **NO-GO** on score |
| 151 | HCT116, target-gated lambda | `-0.004344` | **NO-GO** on score |
| 152 | HCT116, gene-wise ramp | `-0.003235` | **NO-GO**; family closed |

Per-target PDS was exactly unchanged in every arm — the emitter only ever moves
the DE moment, never retrieval geometry.

**The failure is `reach`, and it is diffuse.** On HCT116 the uniform transform
lost `-0.0451` scaled reach. Experiment 151 halved lambda across the entire
low-consensus half of targets and recovered only `+0.0015` of that. So
cross-context agreement measured **between targets** does not localise the
damage. Meanwhile NMAE and Jaccard improve in every HCT arm. The emitter is
trading directional depth for lfc calibration, which wins on H1/HEK293T and
loses on HCT116.

**The harm is linear in dose.** The three HCT116 arms differ only in how much
relaxation they apply, with mean lambda over all 132 targets:

| exp | placement | mean lambda | delta |
| --- | --- | ---: | ---: |
| 150 | uniform | `1.0000` | `-0.005045` |
| 151 | target-wise | `0.7727` | `-0.004344` |
| 152 | gene-wise | `0.6170` | `-0.003235` |

`delta = -0.004606 * mean_lambda - 0.000539`, `R^2 = 0.9449`, residuals an order
of magnitude below the bootstrap intervals. At zero dose the fit predicts
`-0.000539` — i.e. doing nothing. **No forward dose reaches a positive gate**,
and two independent placement schemes (between targets, within targets) both
failed to beat the dose line by a meaningful margin.

**Do not** revive uniform, target-gated, or gene-gated dual-moment relaxation,
and do not tune its percentiles or coefficients on HCT outcomes. The whole
observed range is explained without a new parameter.

**Post-hoc closure — it was never context-dependence.** A source-only diagnostic
(`analysis/dual-moment-reach-closure.md`) found one law operating in both
contexts: the emitter **replaces** the incumbent's directional ranking with one
of roughly fixed quality instead of refining it. Regressing per-target
`reach_delta` on incumbent reach gives break-even at `+0.0173` on HEK293T and
`+0.0163` on HCT116 — agreement within 6% across two unrelated cell lines. HEK's
incumbent barely reaches (mean `0.013`, 68/132 targets at exactly zero) so
replacement upgrades it; HCT's reaches well (mean `0.101`) so replacement
destroys it. Same law, different input distribution.

The gate that would fix it is unobtainable: per-target incumbent reach does not
transfer across contexts (`spearman = -0.033`;
`P(HCT reach=0 | HEK reach=0) = 0.290` against a base rate of `0.292` — zero
information), and no source feature tested predicts it (`|rho| <= 0.16`). A
greedy per-target oracle using the held HCT outcomes reaches only `+0.007207`,
against experiment 152's `+0.006` gate. **The whole family had at most `+0.0072`
available, and only with perfect foresight.**

**Do not** build further variants that re-weight, re-place, or re-dose this
transform. A successor must *combine* the incumbent and emitter rankings rather
than substitute, and needs a build-time confidence signal that this analysis did
not find. Note experiment 101 already rejected a sparse wrong-sign correction.

## Current research priorities and recent decisions — 2026-08-30

### 146: context-gated sign repair — SOURCE-ONLY NO-GO

The fixed full-feature HGB produced macro signed-credit delta
`+0.0000119081`, but only 4/16 contexts improved and only 4/16 reached 65%
weighted flip precision. Its gain was far below the availability-stratified
permutation-null p95 (`+0.0003156613`), FDR `0.10` sensitivity was negative,
and 63.76% of positive gain came from one context. It did not confidently beat
zero-tail or donor-only HGB controls. Macro raw PDS changed
`-0.0000770895`, with minimum context `-0.0009259259`.

The source package and gate calculation passed independent read-only audit.
No challenge controls, development/held 613/410 outcomes, counts, scorer,
build, or submission were opened or created. Do not promote, integrate, build,
or submit this method. Full closure:
`experiments/146-context-gated-sign-repair/RESULTS.md`.

### 082: ESM2 biological target-kernel residual — NO-GO

The first valid guarded run (`20260830T042050Z`) selected identity before any
held-panel truth was opened. All 60 nonidentity candidates were negative in
all four outer-donor × inner-teacher cells. The least harmful candidate
(rank 32, ridge 1, beta 0.02) had macro raw-PDS delta `-0.005120` and worst-cell
delta `-0.005718`. All 27 recovered artifact hashes verify, and the AWS worker
is stopped. Do not reopen this linear ESM2 residual family without a
materially different mechanism.

### 104: direct-credit adapter — NO-GO

The full STRING/context ridge improved continuous signed credit and PDS on all
three held surfaces, but it reduced H1 top-20% sign precision by `0.0414` and
failed the frozen permutation-null gate on H1 and HCT116. The signal is real
but not safe enough to ship.

### 105: TeloHAEC residual source — NO-GO

The first extraction incorrectly counted cells from sub-threshold guides
toward the 20-cell target minimum. The repaired, hash-sealed extraction retains
2,167 targets and 19,015 targeting guides; every retained target has at least
20 contributing cells. It covers 65/300 VCC perturbations and zeros all 300
panel coordinates exactly. Its frozen direct-source gate then lost on both
held contexts: equal-blend PDS deltas `-0.01825` on HCT116 and `-0.01731` on
HEK293T. The preregistered residual follow-on is closed unless a new mechanism
explains that consistent deficit.

### 106/113: official STATE reproduction and H1 domain check — NO-GO

The Arc source is pinned at commit
`9bbfe78a434a55205e4de834e1ea99f85f7a3add`; the official notebook and
isolated dependency lock are hashed. The fixed step-20,000 `state_sm`
checkpoint completed, froze, and failed decisively. On the 120-target held
screen, STATE PDS was `0.50252 / 0.50028` versus incumbent
`0.80917 / 0.79454`; a 10% direction blend lost `-0.02574` macro PDS.
Experiment 113 then falsified the “context transfer only” explanation on H1:
four independent control banks gave PDS `0.4940–0.5133`, pooled PDS
`0.50448`, versus incumbent `0.86345`. Do not build, submit, blend, or retrain
this checkpoint.

### 110: STATE disagreement arbitration — NO-GO before retraining

The perfect-information HEK293T disagreement-routing ceiling is only
`+0.00686`, below the preregistered `+0.015` gate. The original scaffold also
applies only on a reference-significant mask unavailable at deployment.
The two additional 20,000-step checkpoints were not trained.

### 112: Orion control-label correction — algebraic no-op

The `NEGR1` control-label bug translates every Orion donor row by the same
vector. Donor and whitener row-centering remove that translation exactly, so
the current production prediction cannot change except for floating-point
roundoff. No expression or held truth was opened.

### 109: gauge-invariant guide nuisance — SOURCE-ONLY NO-GO

The repaired extractor independently pools all cells for compact-target
reconstruction while retaining the 20-cell floor for guide effects. On the
complete 1,989,578-cell K562 source it reproduced all 9,787 compact target
effects with maximum absolute error `4.768e-7` against the frozen `2e-5`
gate. The guarded source-only fit then selected exact identity (`eta=0`):
untouched-test cosine delta, normalized-MSE delta, both null p95 values, and
the target-bootstrap lower bound were all `0.0`. The family stopped before
held-context access.

### 114: DepMap co-essentiality — terminal NO-GO

The official DepMap Public 24Q4 source is verified at 428,678,699 bytes,
MD5 `6edf7ade09b9b34199210b559d4745d3`, and SHA-256
`3d8f3ec6dbf2db7ff834b79b508622ec0b226f3518003fe96ecf5a4fcf167e3b`.
Its actual schema is 1,178 cell lines by 17,916 genes, corrected through a
frozen post-download amendment before feature construction. The built
`17,916 x 1,178` float32 feature artifact has SHA-256
`7f42b585f788db95cdc345d44b44926bd76ce0fc5ae94c544fd2f6929b9ddc51`
and passes 20 synthetic tests plus bit-exact sampled-source reconstruction.
The repaired frozen weak-row replacement was then evaluated on its authorized
35-target HEK293T diagnostic. The primary STRING+DepMap arm lost
`-0.12016807` raw PDS; STRING-only and DepMap-only were also strongly
negative. H1 was not opened. Do not reopen leave-target-out co-essentiality
replacement without a materially different mechanism.

### 080: listwise metric adapter — terminal NO-GO

The repaired leakage-safe adapter selected a nonidentity rank-8 model and
passed every source-only pre-gate, with untouched non-panel macro PDS gain
`+0.023376` above its complete-selection null p95 `+0.021217`.

It did not transfer to the held panel. HCT116 was `-0.0006884`, HEK293T was
`+0.0006435`, KOLF was `+0.0006868`, and the HCT116/HEK293T primary macro was
`-0.00002245`. Bootstrap 95% was `[-0.001863,+0.001803]`, below the null p95
`+0.0001197`. `promotion_gate.passed=false`.

All 331 recovered files, the 328-file artifact manifest, source-code hashes,
303 model files, nine query files, three cluster files, selection freeze, and
permutation plan verified. The AWS worker is stopped. Do not integrate,
retune, or run count-level confirmation for this adapter.

### 107: GSE249595 high-MOI Jurkat — terminal NO-GO

Official metadata confirms a genome-scale Jurkat CRISPRi source with 83,401
guide records and 299/300 VCC targets (`ZNF714` absent). It is not a normal
single-guide donor: the paper reports a median of 13 guides per cell, and the reliable
readout is a 375-gene targeted-primer panel. The source model is frozen at
SHA-256
`03759a7b2cafa86d2643a58669235aca7fc38a33599aa403c67439f8979a7e26`
with 17,615 transferable targets and 349 exportable response genes. A strict
27-role artifact contract now covers six Nadig galleries plus H1 and HEK293T,
including last-occurrence duplicate-axis handling and exact a14 settings.
The complete immutable v2 evaluation is in
`experiments/107-jurkat-high-moi/work/strict-evaluation-results-v2/`.
Nested six-gallery PDS improved by `+0.00424563`; the target-cluster bootstrap
lower 95% bound was `+0.00127271`, above the 99-permutation maximum
`+0.00115942`. The effect nevertheless failed its frozen promotion gate:
one held gallery was negative (`-0.000356745`) and directional credit was
negative on both independent held contexts (H1 `-0.000164807`, HEK293T
`-0.000103447`). `promotion.passes=false`. Result-manifest SHA-256:
`8f5b8e748b2613239d3d7f0e9e295441e10f8c6e962865b88286312fbbc2d595`.
Do not tune fractions, masks, or target subsets against this outcome.

### 129: Pan response-neighborhood transfer — SOURCE-ONLY NO-GO

The materially distinct Pan-neighborhood learner used GSE295214 only as a
cross-target geometry; it never added the query target's Pan effect. The
source-only selector nevertheless chose identity. Every nonidentity
neighborhood reduced complete-gallery PDS, with macro deltas from `-0.0110`
to `-0.0209`, even though several candidates improved signed directional
credit. The untouched 55-target source test remains sealed. This result is
the clearest current example of the central tradeoff: adding broad direction
can make each query more similar to decoys faster than to its own truth.

### 130: TeloHAEC panel confirmation — terminally closed

Experiment 130 was a metadata-only fail-closed confirmation protocol for a
passing Experiment 125 candidate. Experiment 125 is an immutable NO-GO, so no
candidate satisfies the unlock condition. `CANDIDATE_BINDING_AMENDMENT.json`
remains absent and no TeloHAEC panel expression value was opened under this
experiment. Any future TeloHAEC confirmation requires a materially new
hypothesis and a new experiment number.

### 131: HUES8 and WTC11 differentiated-cell sources — access-limited

HUES8 ESC, HUES8 definitive endoderm, and WTC11-derived hepatocyte each retain
91/300 VCC targets after QC, with broad 9,519 / 9,584 / 11,138-gene response
axes. WTC11 cardiomyocyte retains 91 targets in the canonical QC summary
(90 with at least 100 cells) and has a 7,260-gene response axis. The useful
canonical H5MU or calibrated gene-level tables require authenticated,
certified Synapse access. A small public cardiomyocyte log2FC table exists,
but it is the same sparse portal family already rejected by Experiment 081
and is not pinned to the canonical H5MU run. Do not download the 16–30 GB
objects until compact signed tables or verified access become available.

### 132: DepMap PCA50/Ridge10 response model — SOURCE-ONLY NO-GO

The exact upstream PCA-50/Ridge-10 family carries real target information but
is far weaker than preserving the observed same-target source mean. On the
613-target frozen selection split, its biological-family macro PDS delta was
`-0.204710`; all four PDS contexts failed and macro signed credit was
`-0.002369`. Identity was selected and the untouched 410-target test remains
sealed. This closes replacement-style target-feature regression; any future
context model must be residual-only and preserve the incumbent response.

### 135: TeloHAEC confirmation of Experiment 122 atlas — terminal NO-GO

The independently reviewed, hash-frozen evaluator opened the previously
unaggregated panel responses once for 65 TeloHAEC targets. It compared the
exact staged-a15 incumbent against the same row after exact Experiment 128
integer-budget transport. Controls used pooled cell pseudobulk, libraries
replayed all 300 targets from one `default_rng(0)`, support did not gate finite
coordinates, and both null families preserved the structural self coordinate.

The atlas improved macro signed credit by `+0.004066`, above both null p95
values (`+0.003255` target and `+0.002738` response; p `0.025` and `0.020`),
but missed the frozen `+0.005` gain and its bootstrap lower 95% was
`-0.000456`. Top-20% magnitude-ranked sign precision worsened by `-0.005245`
(`0.515000 -> 0.509755`), with both null and bootstrap gates failing. Median
transport sign retention was `0.890`, below the `0.900` gate. All integer
invariants passed.

`promotion_passes=false`. TeloHAEC is consumed for this exact atlas. Do not
integrate, rescale, subset, or retune it, and do not build or submit an
a15-plus-panel arm.

### 136: DepMap Ridge residual/margin screen — SOURCE-DEVELOPMENT NO-GO

The Experiment 132 Ridge output was tested only on its already-consumed
613-target selection split as a residual around the exact same-target source
mean. Every nonzero global dose reduced PDS in every scored context. A
deployment-available source-rank-preserving router lost `-0.020302` macro PDS,
while exact nearest-decoy-margin preservation collapsed to identity.

Even an invalid held-outcome oracle choosing one dose per target gained only
`+0.001926` macro PDS, below the `+0.003` gate, and remained negative on the
K562 family. Experiment 132's untouched 410-target split remains sealed.
Do not spend it on this Ridge residual family.

## Closed by measurement — do not re-open without new evidence

Everything calibration: amplitude (both directions — at challenge gate depth `nmae` wants it
DOWN like `mse`, and `mse` charges ~10x per 1% either way), emission `spread`/`phi`/call-set
geometry, the whitening operator (`pds` costs 3.4-8.7x any `nmae` gain), per-gene and
per-perturbation shrinkage (buys `mse`, pays it back ~1:1 across `pds`/`nmae`/`reach`),
per-decile magnitude calibration, the `nmae` normaliser (**closed by reading the spec** —
reference-only), a constant directional bias (**structurally inexpressible**: gated genes are
99.4% of profile mass and renormalisation annihilates it), the on-target channel (the scorer
excludes each row's own target from **every** member), and a low-rank ridge transfer
(20 of 20 permuted maps beat the real fit).

**Donor expansion, closed 2026-08-29.** Three of the seven prepared donors carry **zero** of the
300 targets — they are common-essential libraries and the VCC panel excludes common-essential
genes. Four genuinely usable ones exist outside `data/donors/` (`kolf2` 273/300, CD4 x3
~284/300), and coverage was never the constraint: the best transferable directional multiplier
from any added donor is **x1.026**, against x1.478 for the confidence-weight family already
built. The agreement key gets **worse** with more donors — it is sign unanimity, an
intersection, so the graded fraction falls 0.206 -> 0.104 (3) -> 0.029 (5) -> 0.015 (6).
Uniform donor weights have now survived three separate objectives on three surfaces.

A third harness was built in the process: `data/harness/orion-hek293t-thinned.h5ad`, depth-gate
PASS at 0.9999.

## The honest position

Calibration is spent; a15 is published and is the current live best entry.
**Winning needs a better predictor of which genes move and which way** — not another metric
lever. The donor-agreement key remains the only confidence signal that has ever transferred
across surfaces, and it captures about a fifth of a 1.97x per-perturbation oracle ceiling.
Experiments 080, 107, 109, 111, 114, 123, 125, 129, 132, 135, and 136 are
closed at their frozen or source-development gates. Experiments 137, 138,
and 139 have since closed at their development gates as terminal NO-GOs;
their untouched test rows remain sealed. The next experiment must use a
materially distinct source or mechanism rather than reopening those families.

## Experiments 101–103 — new directional sources

- **101, sparse wrong-sign correction: NO-GO.**
- **102, primary COMPASS donor-axis surrogate: NO-GO.** Its residual-only diagnostic was
  positive enough to justify one preregistered confirmation.
- **103, residual-only staged-a15 confirmation on KOLF2.1J: NO-GO, near-positive.** Exact
  staged-a15 baseline PDS was `0.620408`; the candidate reached `0.632124`, delta
  `+0.011716`. The target-cluster bootstrap was `[-0.003138, +0.027365]`, the conservative
  20-null maximum was `+0.011757`, and the add-one one-sided p-value was `0.0952`. The mean
  gate passed, but both uncertainty gates failed. KOLF access for experiment 103 is consumed.
  Do not tune projection strength, axes, target subsets, or null seeds against this outcome.

**Experiment 085, the source-only PerturbMatch inclusion gate: NO-GO before KOLF.** It was run
from the exact matching source at commit `2d5cd1d`. The fixed selector chose `always_D`; the
best nonidentity rule had macro delta `-0.000411` and was negative on HEK293T. The source-only
gate therefore stopped before the durable freeze or any KOLF authentication. KOLF remains
unused for this hypothesis. Do not reopen this displacement-percentile stump family without
a new external feature or mechanism.

## Historical campaign notes below this line

The sections below preserve earlier decisions and their reasoning. They are
not the current execution queue; use the sections above for live status.

## The finding that changed the project

For thirteen submissions we scored **exactly 0** on `expr_mse_unbiased_capped_norm`, as 376 of
385 teams still do, and `docs/STRATEGY.md` had written it off as "a solved constraint". It was
not a modelling wall. The metric *offers* a sampling-noise credit sized by how dispersed the
emitted cells are; ours were **0.67x the reference's**, so we forfeited most of it, and our
amplitude sat past the point where claiming it saturates. `VCC_ALPHA` 1.0796 -> 0.838 took raw
`mse` from **1.3524 to 0.9100** and the member from **0 to +0.0820**.

Two things made it findable, both in plain sight: the board publishes the **scaled** members
directly (their mean reproduces `score_avg` to 3e-16, where our own anchors reproduce it under
no clipping rule), and the board's `description` field, which **131 of 385 teams** populate —
including two teams' explicit metric-calibration runs that anchored the whole analysis.

## Two arms are built, gated and waiting

| arm | change | predicted | status |
| --- | --- | --- | --- |
| **a12** | dispersion to the `#247` correction cap + the cliff amplitude | raw `mse` 0.8882, `~+0.0038` | built, probed (`n_pred` 4440), pre-registered twice |
| **a13** | a12 + the two-part donor-agreement shrinkage, `--effect-scale 2.4004` | raw `mse` 0.8212, `~+0.0205` | built, probed (`n_pred` 4307), **DE-gated positive on both harnesses** |

`a13`'s configuration is now reproduced by passing `--shrink-key agreement --shrink-pert-key
donor-cos` alone; every other value is the default.

**The DE gate is the reason to trust a13 over the alternatives**, and it reversed both research
agents' own recommendation. Every arm more aggressive than a13 — including the one they were
recommending an hour earlier — goes **negative on orion**, and `nmae` does it every time.
`nmae` is the binding member on this family and was priced by nobody: it scores absolute
fold-change error over **the reference's** DE genes, which are not the genes our donors agree
about, so damping by donor agreement degrades it monotonically.

Also from that gate: **`reach` is unreadable on these arms.** Same-policy re-seeded twins differ
by `0.0207` / `0.0245` raw, an order of magnitude above every between-arm difference — which
retires the favourable reach readings as well as the unfavourable ones.

## The honest ceiling

Calibration tops out near raw `0.870`; the leader is at `0.8270` with an aggregate directional
cosine **at least 1.28x ours** (a late correction that reversed a more flattering reading).
**Calibration buys rank 3-5. First place needs the direction.**

**That lead is now closed, as a negative.** The whitening operator does cost `1.5-2.1x` on
directional accuracy, monotonically on both surfaces — and gated across all six members it
cannot be bought. `kappa` falls and `pds` rises monotonically in `gamma` on both surfaces, and
**no `gamma` is net positive on both**: every setting below 0.65 gains on h1-test and loses on
orion, crossing at ~0.45 (`-0.0001`, inside noise). At `gamma = 0` the arm reads `-0.0033`
(h1-test) and `-0.0080` (orion).

The exchange rate: **each `+0.01` of `kappa` costs `0.0181` raw `pds` on h1-test and `0.0517`
on orion.** `pds` alone more than eats the entire directional gain on both surfaces before any
DE member is counted. `nmae` disagrees in sign and reads as zero here; `fid` gains slightly;
`reach` was unreadable again. **The whole result is `pds`. Production's `WHITEN_GAMMA = 0.65`
is at or beside the two-surface optimum — do not change it.**

This was a deliberate re-opening of a lever closed under a different objective, per §9's
asymmetry note. It cost one measurement and returned a null, which is what that rule is for.

**So the direction cannot be bought from the operator, and has to come from somewhere else
entirely.** The donor-agreement family remains Experiment 097's main high-leverage result
(a13, `+0.0215 / +0.0040`) against a per-perturbation oracle ceiling of `1.97x`, most of which
is unbuildable. One small exception remains unresolved: target-specific
`normrel^-1.5` donor weighting was positive on two old analytic surfaces, but
has never been tested in the exact a15-shaped pipeline.

**Measured dead ends, each with the number that kills it:** shrinking toward our own grand mean
(captures 0.03% — our panel mean is uncorrelated with the reference's); the on-target knockdown
(the scorer drops each row's own target from **every** member, so the one signal transferring at
`rho = 0.647` is worth exactly zero); global donor re-weighting (tuned shares land on opposite
sides of 0.5 and each loses 8% on the other surface — which *validates* production's uniform weights);
and the Ahlmann-Eltze zero-shot transfer (**20 of 20 permuted maps beat the real fit**).

## The next two slots — decided, with the reasoning

**Slot 1: `a12` as built.** One coherent change (dispersion to the correction cap + the cliff
amplitude), already built, probed and pre-registered twice.

**Slot 2: `a13` = a12 + the PER-GENE shrinkage only** (`--shrink-key agreement`,
`--shrink-toward zero`), with its amplitude re-solved. **Not** the per-perturbation multiplier,
even though it is built and tested.

**REVISED, same evening, on measurement.** The attribution objection below is void and the
per-perturbation half goes in the same arm. `pds` is **exactly neutral** to `gamma_pert` —
structurally, because `pds_cosine` normalises each row before taking the cosine, and measured
identical to four decimals at `gamma_pert` 0 and 1 across every `gamma_gene`. So any `pds`
movement in a two-knob arm attributes **entirely** to `gamma_gene`, which is the thing that
would otherwise make it unreadable. The `mse` contributions are near-additive (`+0.00226`
alone, `+0.00323` stacked), so the knobs carry different information rather than duplicating.

**And the harness limit recorded below was the wrong diagnosis.** It is not that our surfaces
cannot choose between the unanimity and pairwise *keys*. It is that **`orion` cannot run a
three-donor key at all** — `orion_hct116` is itself a production donor and must be held out, so
that harness is structurally two-donor. Forcing h1-test down to two donors reproduces orion's
answer *and flips the sign of the pds result*. h1-test at three donors is leak-free and matches
the challenge panel closely (median non-zero genes per target 1,964 against 1,979). So the
disagreement is **donor count, not cell line**, orion is measuring a different estimator, and
rule 1's "take the smaller of the two harnesses" does not apply here.

**The original reasoning, kept for the record:** Combined, the two halves remove ~95% of the
amplitude (mean per-gene `alpha` 0.0474, mean `z^2` 0.2019) and need roughly 3x compensation
between them. That is the largest single transformation this pipeline has ever shipped, it has
never been on the board in any form, and the two halves interact through the amplitude solve
rather than stacking cleanly — so a disappointing result would be unattributable. The per-gene
half has a corrected, exact-identity value (`+0.0055` over a12); the per-perturbation half's
value is measured only on surfaces that **structurally cannot see it** (below). One at a time.

**A structural limit worth stating plainly.** The unanimity key grades 0.106 of genes at three
donors; the pairwise key grades 0.756. Both harnesses carry **two** donors, where the two keys
are provably identical and both grade 0.208. **Neither harness can choose between them, and
neither reproduces the estimator that would actually ship on the three-donor panel.** Any
choice between those keys made on harness evidence is unfounded — the open question is the
damping *level*, not the key, and it needs a realisation sweep over gamma rather than a
key comparison.

## Next, in order

1. **`a12-tier8000-cap` is built, checksum-verified and on disk.** Raises emitted dispersion to
   the metric's `#247` correction cap via a tiered spread policy
   (`--keep-top 8000 --keep-spread 0.25 --tail-spread 1.0 --overdispersion 4.0`) and trims
   `VCC_ALPHA` to 1.012. Predicted raw `mse` **0.878**, `+0.0192` gross against a *measured*
   DE cost of `-0.0013`. **Probe its `n_pred` before submitting** — the tail policy churns
   8-14% of the call set, and the churn was measured chance-in/chance-out on both harnesses.
2. **The donor-agreement shrinkage (Miller).** `donor_agreement` is the only confidence score
   that survives partialling out effect size on both harnesses (`+0.144` / `+0.229`); both
   self-referential scores fail or flip sign. Worth `+0.0044` on top of (1), and it is the
   only lever that touches the direction rather than the calibration.
2b. **The corrected shrinkage passes the leader on this member.** The donor-agreement
   shrinkage is now implemented in the builder (`--shrink-key agreement`,
   `--shrink-gamma`, `--shrink-toward {zero,mean}`), OFF by default and byte-identical when
   off, 30 tests. Two corrections were needed to get here and both matter:

   - **Shrink toward ZERO, not toward the panel mean.** I specified `mean`, reasoning that the
     metric's scored 0-anchor is the context's mean perturbation response. That assumed our
     predicted panel mean approximates the *reference's* mean effect; measured,
     `cos(ours, theirs)` is `+0.016 / -0.066` — no capture, inconsistent sign. `zero` wins
     12/12 across harnesses, gammas and keys. Default flipped against the specification.
   - **`MSE-MECHANISM.md`'s published `+0.0044` for this family was net NEGATIVE** (`-0.0030`).
     It was scored through a `sqrt(W/Wg)` shortcut that assumes `A ~ B`, on the one estimator
     family where `A - B` is exactly what moves (`(A-B)/A` goes 0.014 -> 0.077 under the mean
     target). Intervention 1 is unaffected — it was always scored through the exact identity.

   Corrected value: **raw `mse` 0.8203, `+0.0102` on top of a12** — **but treat that number as
   NOT ESTABLISHED.** Three independent routes agree that shrinking toward zero beats shrinking
   toward the mean, and that ordering is solid. All three, however, hold `claim` at its
   *unshrunk* value, and `claim` is a functional of the emitted cells that the optimum
   `s* = sqrt(claim/B)` sits directly on. **None of them re-emitted, so none of them could have
   disagreed about it** — rule 10. The ordering should survive; the levels, including "passes
   the leader's 0.8270", are unverified until `mech_challenge.py` is run on a built artifact.

   **The sharper risk: the harnesses understate this family's penalty, and understate it most
   on the arm we mean to combine it with.** `claim/B` is 0.183 (h1-test) and 0.369 (orion)
   against **0.879 on the challenge at the #247 cap** — so the board charges an eps-inflating
   estimator 4.8x / 2.4x what our surfaces report, and a12's entire design is to raise `claim`
   *to* that cap. Shrinkage is therefore most dangerous stacked on a12; the a11 base may be the
   safer host.

   Needs its amplitude re-solved before it can be built: mean surviving amplitude at gamma=1 is **0.108**, and `VCC_ALPHA` cannot
   supply the correction because PerturbMatch and CD4 are added at fixed magnitude after it —
   use the post-gate `--effect-scale`.

3. **Raising the directional cosine — still the only route to a durable first place.** Calibration tops out at
   raw `~0.870`; the leader is at `0.8270` with an aggregate cosine **at least 1.28x ours**
   (this reverses an earlier, more flattering reading — see the retraction in
   `MSE-MECHANISM.md` §(b)). Our directional signal is extremely concentrated: the top
   quartile of perturbations carries **80-85%** of it and **12-29%** point the wrong way.

## Closed — do not re-open without new evidence

The **emission channel is exhausted** (`experiments/094-reach-ordering/A10-BOARD-RESULT.md`):
`spread`, `phi`, `n_pred` and the call-set geometry now move the board by less than `0.001`.
a10 changed the call set by 24% and moved the score by `-0.00018`. Also closed: amplitude as a
*global pre-096 lever*, `min_cpm`, call-set truncation, `basal_beta`, basis-prune tuning,
spread-based ordering, per-gene `phi`, per-gene/per-perturbation shrinkage fitted on a slope.

**Note the asymmetry that cost us thirteen submissions:** "closed" means *it did not move the
members we were then measuring*. Amplitude was closed against the DE members while `mse` sat
in its clamp and was therefore invisible. A lever closed under one objective is not closed
under a new one.

---

## Leaderboard

- **Live artifact: `a8-nb-phi4`**, entry `lfnx5c9vp1gVhlNWcemT`, model name `Wildflower`,
  overall **`0.1543607`** — **our best ever** — submitted 2026-08-27 19:21 UTC.
  Config: `a2-prune` + `--overdispersion 4.0`.
- **Rank 7 of 363 teams**, up from 11. Board top: Aginglab.com `0.18991`,
  Ibrahim Mansour `0.17078`, Jurassic Park `0.16666`, Kate P `0.15617`,
  GISL `0.15602`, Vivai `0.15530`.
- **We hold the board maximum on `pds`** (`score_pds` 0.7366, raw cosine 0.8318).

### Today's two submissions, both measured exactly

| | overall | Δ | what changed |
| --- | ---: | ---: | --- |
| v11 (previous best) | 0.1501827 | — | — |
| `a2-prune` | 0.1503084 | +0.0001257 | exp091 pruned whitening basis |
| **`a8-nb-phi4`** | **0.1543607** | **+0.0040523** | gamma-Poisson over-dispersion φ=4 |

`a8` member by member: `pds` **+0.000000 exactly**, `nmae` −0.000038, `fid` +0.006393,
`reach` **+0.021398** (the largest single-member gain this project has recorded),
`jac` −0.003439.

### The gap to the leader, decomposed

| member | our raw | Aginglab raw | **they gain (overall pts)** |
| --- | ---: | ---: | ---: |
| `reach` | 0.1701 | 0.2674 | **+0.0183** |
| `nmae` | 0.9441 | 0.8923 | **+0.0140** |
| `fid` | 0.5021 | 0.5135 | **+0.0064** |
| `mse` | 1.5708 | 0.9594 | +0.0068 |
| `jac` | 0.0314 | 0.0287 | −0.0012 |
| `pds` | 0.8318 | 0.8196 | −0.0047 |
| | | | **+0.0396** |

(against the live `a2-prune`. `reach` is now the single largest gap by a clear margin,
which is what arm `a7` targets.)

## Submissions

- **UTC day 2026-08-27: slots in use — see the log below.** Slots reset 00:00 UTC.
- Rank uses a team's **latest** submission, not its best.
- Only the final entry after the Oct 22 test release is prize-eligible
  (Nov 5, 23:59 UTC). **Leaderboard rank is informational only.**
- 🚫 **NEVER spend a slot restoring rank, and never re-submit an already-scored package.**
  Ahmet's standing instruction (2026-09-01) after a19's regression was "fixed" by re-uploading
  a15 byte-identical, spending the day's second slot to buy zero information. If a bad arm
  drops our displayed row, leave it — the row is not scored. **A slot must buy information that
  moves us toward FIRST.** See the banner at `AGENTS.md` §6.

## The scoring map, solved against 780 board entries

`score_avg` is the plain mean of the six scaled members (max |diff| 3.3e-16).
`scaled = (raw − b)/(r − b)`, fitted at `R² > 0.9999` for four members:

| member | `b` | `r` | span | raw Δ for **+0.01 overall** |
| --- | ---: | ---: | ---: | ---: |
| `pds` | 0.50020 | 0.95228 | 0.452 | +0.0271 |
| `fid` | 0.51239 | 0.80964 | **0.297** | **+0.0178** |
| `reach` | 0.07902 | 0.96414 | 0.885 | +0.0531 |
| `jac` | 0.03038 | 0.39727 | 0.367 | +0.0220 |

`score_mse` and `score_nmae` are **not** functions of the reported aggregate raw
values — they are per-perturbation normalized then averaged, so do not do
arithmetic on them. `score_reach` is essentially exactly linear (rms residual
3.3e-4), so reach arithmetic is safe.

## Amplitude is RETIRED as a lever — see experiment 094

STRATEGY Addendum 5 §3 ("alpha's marginal direction is up") is **wrong**.

- Across 288 within-team consecutive submission pairs on the board,
  **`Δscore_avg` is uncorrelated with `Δlog mse`** (r = −0.07, p = 0.24).
- Between teams, `log mse` explains **9.7%** of reach variance. Our reach
  residual against the amplitude fit is −0.0056 (z = −0.15, 43rd percentile) —
  we are already ON that line; there is ~+0.001 there, not more.
- **We527** ships `STATE baseline` at raw mse `1.3577`, within 2.6% of our
  `1.3236`, and scores reach `0.2529` against our `0.1815`.
- **Aginglab** halved amplitude (`1.822 → 0.959`) and *gained* `+0.100` reach.

**v13 (α=2.0) is built and unsubmitted** at
`s3://<your-s3-bucket>/output/exp089-production-pm-cd4/20260825T193634Z/`.
On this evidence it should stay unsubmitted.

## The open research direction: the confidence ORDERING

`reach` walks the **reference's** significant genes in **our** confidence order
(our significant calls first, then ascending predicted `p_adj`) and measures the
deepest prefix that stays 90% sign-pure. The pool and each gene's correctness are
fixed by the model; **the only free variable is the order**, and ours is
inverted: `Spearman(our p_adj, sign correct) = +0.5009`.

- genes we **do** call significant: sign accuracy **0.3446**
- genes we **do not** call: **0.8259** (paired t = −54.0, inverted in 99.2% of perturbations)
- oracle ordering on the harness is worth **+0.110 overall**

Cause: a minimum-variance emitter puts 400 cells at the predicted **mean**, which
for right-skewed counts sits above the **median**, so the rank-sum reports
"significantly up" at zero effect. The reference's DE set is 90.8% **down**;
our most-confident decile is 83.5% **up**-predictions, and predicting up is right
10.8% of the time.

The order is a deterministic, **pre-computable** function of
`|F̄_g − ½| / sqrt(tc_g)` (Spearman +0.996 with the realized p-order), where `F̄_g`
is the mean of the control CDF over our emitted count law. `|predicted effect|` —
the key `--keep-top` ranks on — correlates with it at only **+0.071**.

**The tension to respect:** the same artefact HELPS `fid`. Reference
*non*-significant genes drift 54.8% **up**, and our call set is 76.5%
up-predictions. Removing the up-bias wholesale costs `fid` roughly −0.046 overall
at the extreme while buying `reach`. Arm a1 therefore changes the ORDER only and
leaves the call set alone.

## Today's arms (UTC day 2026-08-27)

| arm | config | package | status |
| --- | --- | --- | --- |
| `v11-control` | v11 exactly, via the fast path | sha `24a1282b…` — **byte-identical to `submissions/v11.vcc`** | control; validates the fast build. Held as the final-slot fallback. |
| `a1-order-k2000` | `--tail-spread 0.25 --keep-top 2000 --keep-spread 0.0` | — | **withdrawn before submission.** Its design assumed the HCT116 down-asymmetry, which h1-test refutes. |
| **`a2-prune`** | v11 + `--whitener-id union_pruned` (exp091) | **SCORED `0.1503084`, `+0.0001257` vs v11 — our best artefact** | **slot 1 spent**, entry `o8rCQvdgjb4I8rSLHKk1` |
| `a3-order-k2000` | v11 + `--tail-spread 0.25 --keep-top 2000 --keep-spread 0.0` | built, verified, held | superseded by a5 |
| `a5-prune-order` | a2-prune + the same ordering flags | built, verified | **slot 2 candidate** |
| `a6-graded-agree` | pruned + graded map on the agreement key | **build killed** — a fetch script stopped the worker while it was still running. Moot: its policy realises at −0.00193. | — |

### Slot 2 is no longer an easy call — the ordering channel is CLOSED

`experiments/095-clean-harness/realize.py` drives the **unmodified production
`emit_counts`** with each policy's spread vector and scores the exact DE path, with a
same-policy second seed as a noise control (**5e-5**). Every counterfactual in this
project's emission work evaporates under it:

| policy | counterfactual | **realised** |
| --- | ---: | ---: |
| block `\|effect\|` K=2000 (a5) | +0.00583 | **−0.00021** |
| graded `\|effect\|` | — | −0.00077 |
| graded agreement (a6) | +0.00678 | **−0.00193** |
| k500 / k1000 / k4000 / k6000 | — | −0.0001 … −0.0007 |
| the three positives (tail 0.40/0.60, graded wide) | — | +0.0001 … +0.0004, all by cutting `n_pred` 15–40% |

`spread` is a **weak lever on order and a strong lever on membership**: predicted p-values
span hundreds of orders of magnitude and `spread` moves them by a few, so a promoted gene
that was not already near the front does not travel. Anything that reorders enough to matter
shrinks the call set enough to lose `fid` — which is exactly how v14 and v15 lost.

So a5's expected board value is **≈0, possibly −0.001**, not the +0.005 it was built for.

### …but over-dispersion is NOT closed, and it is the arm

`spread` caps at Poisson. **Gamma-Poisson scatter does not**, and it is the only emission
change that moves the DE ordering. It is also what the real challenge data says we are
missing: our per-gene CV across 400 emitted cells is **0.38–0.40** against **0.84–0.92** in
the real control cells — **under-dispersed ~2.3×**.

Realised through `emit_counts(..., overdispersion=phi)` on the repaired harness, against a
`−0.00046` noise floor, with **`nmae` invariant to four decimals** (the exact-total repair
holding):

| policy | `n_pred` | reach | Δreach | Δfid | **ΔTOTAL** |
| --- | ---: | ---: | ---: | ---: | ---: |
| production | 3,934 | 0.0838 | — | — | — |
| **spread 0.25, φ=2** | **3,218 (−18%)** | **0.1378 (+64%)** | **+0.0102** | +0.0100 | **+0.0204** |
| spread 1.0, φ=0.30 | 692 (−82%) | 0.1961 | +0.0212 | +0.0167 | +0.0448 |

The big-φ rows' `fid` gains are the harness's coverage artifact (`n_real` 68) and would
invert on the board exactly as v14's did. **`spread 0.25, φ=2` is the operating point**
because its 18% call-set cut lands at `n_pred` ~4,880 on the challenge, where coverage is
still ~0.99. Board-anchored: ≈ **+0.016 overall**.

**Slot 2 is `a8-nb-phi4`** = `a2-prune` + `--tail-spread 0.25 --overdispersion 4.0`,
entry `lfnx5c9vp1gVhlNWcemT`, model name `Wildflower`, submitted 2026-08-27 19:21 UTC.

φ=4 beats φ=3 by `+0.0037` with two seeds each (mean reach `0.16052` against `0.15524`).
Its probe against the real challenge controls measured `n_pred` **3,836** against a
**predicted 3,815** — 0.5% out, through three measured transfers. Prediction pre-registered
at `experiments/094-reach-ordering/A7-PREDICTION.md`: **overall ≈ 0.1716**, `pds`
bit-identical, `nmae` within 0.0005.

`a7-nb-phi2` (φ=2, probed at `n_pred` 4,523) is built, verified and retired unused.

`a5` stacks the ordering flags on the pruned operator. That is **still** a clean
single-variable comparison, because slot 1 measured the operator exactly
(`+0.0001257`), so `a5 − a2-prune` isolates the ordering channel — and it makes the final
entry the best artefact available rather than a strictly worse one. `a3` (same flags on the
v11 base) is built, checksum-verified and held as the fallback.

Experiment 095 put K=2000 on the peak of the promotion curve (`+0.00583` on a repaired
harness), and probing the built artefact against the real challenge controls shows it grows
`n_pred` by only **3.9%** (5,961 → 6,193), which the coverage inversion prices at
**−0.0004**. Expected net ≈ **+0.005**.

## Real challenge `n_pred`, measured — every prior estimate was wrong

| artifact | harness estimate | **real challenge** | board raw `fid` | implied coverage |
| --- | ---: | ---: | ---: | ---: |
| v14 | ~2,317 | **1,563** | 0.49126 | 0.8950 |
| v15 | ~2,370 | **1,911** | 0.49654 | 0.9157 |
| v11 | ~5,340 | **5,961** | 0.50178 | **1.0000** |
| a3 / a5 | — | **6,193** | — | ~1.0 |

The entire `fid` rise across our three board points is **coverage recovery, and it is spent
at v11**. Past that, `fid` is pure precision and declining (−0.0004 at 6,193, −0.0019 at
7,000, −0.0037 at 8,000). Do not grow the call set further.

`submissions/v11.vcc` is on disk and byte-verified, so a restore needs no build and no
download.

## What experiment 095 settled, and what it could not

Rebuilt the production prediction against a depth-repaired harness:

- **The `|lfc|`-ordering gain is real.** `+0.00448` overall against `+0.00497` on the
  broken harness — agreement to 0.0005 across a repair that moved the reference's direction
  balance from 0.908 to 0.609. Not a confound artefact.
- **Block promotion beats full ordering.** Promoting the top-K ahead of the rest, keeping
  the emitted p-order inside the block, peaks at K=2000 (`+0.00583`) and exceeds a pure
  `|lfc|` sort (`+0.00448`).
- **No harness we hold can settle a call-set question.** `n_real` is 322 (broken), 68
  (repaired) and 222 (h1-test); all pin `fid`'s coverage term at 1, so every one of them
  says "truncate" — and the board refuted that twice. Our three board points are the only
  data in the challenge's regime, and they say **more calls is better** (raw `fid` 0.49126,
  0.49654, 0.50178 at `n_pred` 2317, 2370, 5340).

## Cloud

| instance | state | note |
| --- | --- | --- |
| `i-08300b395f88f6f55` | **running** | exp094 fast-build worker, r7i.8xlarge |
| `i-0bd65df670234f97e` | stopped | exp091 essential-prune worker |
| `i-06197bcff9375a2e1` | stopped | exp091 random-prune worker |
| `i-08d90ae88c7b210ed` | stopped | older clean-HCT scorer |

**A build is 18 minutes, not 3h20m.** `experiments/094-reach-ordering/aws/fastbuild.sh`
reuses the already-extracted tree on the worker; the exp089 ceremony runner spent
~3 hours re-downloading the donor tree, the controls bundle and the Python
environment around a build whose own `time -v` reads 16:55. Stop the worker when
idle — it burnt ~$95 idling 45h across 2026-08-25/27.

## Experiment gates

- **exp094** confidence ordering: **the live direction.** See its RESULTS.md.
- **exp091** basis prune: POSITIVE, not shipped. `+0.002409` (t=19.85) against a
  size-matched random control. Needs a whitener refit to ship.
- **exp093** call-set truncation: measured, shipped as v14/v15, both **negative**
  on the board (`0.147273`, `0.149436` against v11's `0.150183`). Its `fid` gain
  did not transfer because it moved `n_pred` the wrong way; coverage is already
  saturated at v11.
- **exp086** sparse-responder emission: re-openable.
- **Lineage-matched donors:** closed on data. `nadig_jurkat`, `rpe1` and
  `k562_essential` cover 0/300 of the panel.

## Corrections logged this session

- `experiments/093-spread-truncation/RESULTS.md` §2 claims the quantile key
  "saturates at 0/1 for thousands of genes at once". Measured: **0 of 10,258**.
  No gene saturates, `pred_p == 0` never occurs, and the `|lfc|` tiebreaker in the
  reach ordering is **dead code** on our submissions.
- The HCT116 harness is **reach-blind**: raw reach 0.0147 against the board's
  0.1815, `k*` = 1 on the median perturbation, 31.8% of perturbations score
  exactly 0. Trust its ordering *directions*, never its reach *magnitudes*.

---

## What to do next, in order

1. **Re-verify the v11 byte-reproduction.** `src/vcc2026/donor_transfer.py` and
   `src/vcc2026/multi_donor.py` both changed on 2026-08-27 (over-dispersion,
   `donor_agreement`), and both frozen implementation hashes were updated deliberately. The
   unit tests assert the default paths are unchanged and three arms built cleanly, but the
   end-to-end guarantee — an unflagged `fastbuild.sh` hashing to
   `24a1282b1ca69e17a5b1dfa3c4279817dedff725b0f885fdf785e4b28e5c22ef` — has **not** been
   re-checked since. One 18-minute run, and everything below assumes it.

2. **Read `a8`'s board result against `experiments/094-reach-ordering/A7-PREDICTION.md`.**
   The prediction (`overall ≈ 0.1716`) and its falsifiers were written before the slot was
   spent. Check `pds` first: it must be **bit-identical**. If it moved, the exact-total
   repair does not survive φ=4 and every over-dispersion number is void.

   The outcome that matters most is not the score — it is whether `realize.py` is
   calibrated. A validated local gate for the emission channel makes the next twenty
   experiments nearly free, and this project has spent three board slots on counterfactuals
   that the emitter would never have produced.

3. **Tomorrow's slot 1 is built and cross-validated: `a10-nb-s015-phi16`**
   (`--tail-spread 0.15 --overdispersion 16.0`, pruned operator). Package at
   `submissions/a10-nb-s015-phi16.vcc`, checksum verified against the worker.

   It **replaces** `a9-nb-s020-phi4`, which was built earlier the same day and then shown to
   be Pareto-dominated — a10 is better on *both* `n_pred` (3,553 vs 3,481) and `reach`
   (0.1565 vs 0.1473). `submissions/a9-nb-s020-phi4.vcc` is kept but **should not be
   submitted**; see `experiments/095-clean-harness/RESULTS.md` sections 21-22.

   The case for a10, in one line: **`spread` and `phi` are separable knobs** — `spread` sets
   `n_pred`, `phi` sets ordering quality at that `n_pred` — and the optimal `phi` rises as
   `spread` falls (0.25 -> 4, 0.15 -> 16, 0.10 -> 16-32). Tuning `phi` at the fixed
   `spread=0.25` we used through a8 could not find this.

   **Both harnesses were run, per a8's own rule, and they disagree about a8 in a readable
   way.** orion (`n_real` 92) ranks a8 first; h1-test (`n_real` 322) ranks it *last* and
   makes a10 the only positive arm (+0.00236 against a +0.00047 noise floor). The
   disagreement is the `n_real` artefact: orion rewards call-set cuts mechanically and a8
   cuts hardest. The board's `n_real` is ~2,000, so the trend runs against a8.
   Decisively, **a10 holds `jac` flat on h1-test (0.0682, unchanged) — the exact member a8
   lost on the board (-0.0034).**

   **Falsifier, pre-registered:** if a10 scores below a8's `0.1543607`, then `reach`
   magnitude beats the `n_pred` charge on the challenge, orion was right, and the next move
   is *higher* spread (`s0.30_phi4`), not lower.

4. **Do not re-open**, all closed by measurement on 2026-08-27: amplitude/alpha, `min_cpm`,
   global spread, call-set truncation, `basal_beta`, `mse` via low alpha or low gamma,
   basis-prune tuning, the `spread`-based ordering channel (nine policies, `5e-4` noise
   floor), and **per-gene φ matched to control dispersion** (loses to a global φ).

5. **Gate every harness with `scripts/check_harness_depth.py` before trusting a DE number**,
   and prefer `h1-test` for absolutes. `orion-hct116.h5ad` is depth-confounded;
   `orion-hct116-thinned.h5ad` is the repaired drop-in.

   **And the deeper rule:** `n_real` is far smaller on every harness we hold (92–322) than
   on the challenge, so **any metric with `n_real` in its denominator or its union is
   unreadable there** — that is `fid`'s coverage term, `jac`'s union, and `reach`'s
   normalisation. Every call-set cut looks free on all three. This artefact has now cost
   three board predictions: v14's `fid`, v15's `fid`, and a8's `jac`. Only the
   **ordering-driven** parts of those metrics transfer.

   **The repair is validated.** exp091's pruned basis is the one change measured on all
   three surfaces: the broken harness said `+0.000876` with `reach` going *up*, the repaired
   harness says `+0.0002` with `reach` going *down*, and the board delivered `+0.0001257`
   with `reach` down. A 7× over-statement with an inverted sign became a 1.6× over-statement
   with the signs right on four of five members.

   **Run every emission arm on BOTH harnesses.** For a8, orion-thinned said `+0.077` raw
   reach, h1-test said `−0.001`, and the board delivered `+0.019`. Neither bracketed it
   alone; their disagreement is the error bar.

6. **Force every emission claim through `realize.py` before it costs a slot.** A reordering
   that cannot be produced by an actual spread vector is not a result. And check the
   *mechanism*, not just the number — a per-gene φ arm read as a clean null on
   2026-08-27 because a units error had silently assigned φ=0 to 97.9% of genes.

7. **The remaining `reach` headroom is a model problem, not an emission one.**
   Over-dispersion captures 14% of the ordering ceiling; `+0.0875` overall sits between φ=4
   and a perfect confidence signal, and closing it means predicting per-gene *direction*
   better — which is also where `nmae` and `fid` live.

---

## 2026-08-30 16:30 UTC — current lead-research state

The live board has 474 entries. `a14` remains our latest submission at
`0.1704242443`, rank 11. The leader is `0.2112235154`. We still lead the
leader on scaled PDS (`0.738614` versus `0.708923`); the gap is in the
non-PDS members. `a15` remains built and unsubmitted.

### Experiment 117 — Pan direct residual: NO-GO

The repaired revision-2 evaluator passed 52 pretruth tests and then completed
one immutable held evaluation.

- cross-context raw-PDS delta: `+0.0011948386`;
- HEK293T: `+0.0000578302`;
- Jurkat macro: `+0.0023318469`;
- bootstrap lower 95%: `+0.0008382064`;
- largest target-label null: `+0.0004315281`;
- HEK293T directional-precision delta: exactly `0`.

The effect is real but below the frozen `+0.003` minimum and fails the strict
positive HEK293T directional gate. Do not tune its fraction against this
outcome.

### Experiment 118 — PDS-null panel residual: broad signal, frozen candidate NO-GO

CD4 source-only validation was strong, but the one-shot 64-target Pan
confirmation separated two claims:

- broad signed credit transferred: `+0.0237657838`, bootstrap lower 95%
  `+0.0114574525`, above both frozen credit nulls;
- confidence ordering did not: top-20% weighted sign precision `0.3547940640`,
  below the target-label null p95 `0.4037518492`;
- the direct unit-scale additive application produced a negative candidate
  mean, while every mass-conservation and PDS-invisibility invariant passed.

This exact candidate is NO-GO. Do not rescale or clip it using Pan. The result
does justify a new independently validated multi-context panel-response family:
the useful signal is broad and PDS-free, but its ranking and positivity rule
must be learned elsewhere.

### Experiments 119–121 — A549 sources acquired; direct donor family closed

- Altair A549: 1,000 targeting aggregates, 7,843 exact canonical response
  genes, 20 measured VCC-panel targets.
- Myllia A549: 53,511 raw-count cells, 175 QC-passing target effects, 16,477
  transferable non-panel genes, 14 QC-passing VCC-panel targets.
- Their preregistered 20-target cross-source test scored macro PDS `0.452632`
  against chance `0.5`, null p95 `0.575`, with negative directional credits in
  both directions.

Do not use either A549 source as an ordinary direct donor on the basis of
same-context reproducibility.

### Active work

- **Experiment 140:** terminal NO-GO. Target-total-UMI-weighted batch-matched
  controls selected exact identity on the full 613-target development set.
  Nested macro deltas were `-0.0000879634` PDS and `-0.0000017071` signed
  credit; only one of five folds selected nonzero eta. The frozen Orion-only
  diagnostic was also negative. All preparation, source-hash, target-axis,
  `union_pruned`, and eta-zero identity checks passed. The 410-target test
  remains sealed; do not integrate, retune, build, or submit this family.

- **Experiment 141:** terminal raw-gate NO-GO. The fixed cross-fitted
  Pan-coordinate primary was strongly positive on the reused 613-target
  surrogate (`+0.0119715` PDS, `+0.0076032` signed credit), but the scalar
  direct-blend control already achieved `+0.0118417` / `+0.0076233`.
  Coordinate calibration added only `+0.0001297` PDS over scalar versus the
  frozen `+0.001` distinctness gate. The expensive null stage and sealed
  410-target split remain unopened. Do not retune rho or reopen direct Pan
  interpolation; Experiment 117 already consumed that family.

- **Experiment 080:** cloud listwise adapter still running on
  `i-08300b395f88f6f55`; recover, hash-verify, and stop the instance immediately
  on completion.
- **Experiment 111:** DepMap context-program fitting is a NO-GO. The selected
  tune-only arm gained `+0.0103486`, but untouched nested-test gain was only
  `+0.00490196`; bootstrap lower 95% was `-0.000953`, geometry gain missed its
  gate, and the result did not beat the context-label null p95 `+0.006808`.
  Held-context access is not authorized.
- **Experiment 122:** source-only phase 1 is complete and selected
  `equal_mean_agreement`. Signed credit was positive on K562 (`+0.0220050`),
  Orion (`+0.0529608`), and CD4 (`+0.0149301`), with both frozen null
  p-values `0.004975` and all 300 targets supported. The source-only result is
  valid, but direct log-effect integration is not exactly PDS-null after
  production normalization.
- **Experiment 124:** the integration audit is a NO-GO. Naive
  log-effect addition changed non-panel abundance in all 900 A/B/C
  context/target profiles. Its proposed simplex transport preserved
  profile-level non-panel values, but used ordinary log-fold units instead of
  Experiment 122's `log1p(50000*p)` units and failed to preserve the summed
  emitted integer panel budget after per-gene flooring. The resulting PDS
  shift is tiny but nonzero. Its recovery-chain hashes are also not
  independently authenticatable.
- **Experiment 123:** exact SLIM/STRING reproduction is a NO-GO. Macro PDS gain
  was `+0.00119048`, below the `+0.003` gate; RPE1 was negative and both the
  PDS and precision bootstrap gates failed.
- **Experiment 116:** source integrity and all 33,372,300 author values verify,
  but the frozen broad-axis gate is a NO-GO: median split cosine `0.06963`,
  reliable-target fraction `0`, and only 13.91% broad-axis coverage. Do not
  open a held context under experiment 116. DLD1 remains the only newly
  identified source with at least 40 VCC targets, so any panel-only use would
  require a new preregistered sparse-response confirmation.
- **Independent-source scout:** HCASMC is the best orthogonal source but has
  only 27 VCC targets; primary macrophage has 12. Neither clears the proposed
  40-target confirmation floor.
- **Experiment 125:** terminal source-only NO-GO. The posterior improved
  top-20% precision by `+0.0121746`, but that gain did not exceed either null
  p95. More importantly, broad signed credit fell by `-0.0150499` and was
  negative in all three held biological families. All 1,800 transport/PDS
  invariants passed and the immutable 16-file result verified. Do not open
  TeloHAEC or HCASMC values for this candidate.
- **Experiment 126:** opaque HCASMC acquisition is complete. The exact
  1,531,381,676-byte SCEPTRE gzip matches published MD5
  `4d3067b3e934842c88b9ea52d74780f9` and local SHA-256
  `b87452c2db1f508f5013bc15dff8088828ab931b35492d77ed6aea569fb5d9b9`.
  Its 27-target overlap is metadata-only; the effect gzip has not been
  decompressed, tested, opened, or parsed. Value access remains locked pending
  a frozen candidate and a new explicit authorization.
- **Experiment 127:** panel-only oracle ceilings are small in aggregate.
  Mutable panel coordinates contain 2.283% of reference DE genes and 2.284%
  of absolute reference LFC mass across the three reused development
  surfaces. Perfect panel correction improves incumbent raw NMAE by at most
  0.024665 on average. The repaired-HCT116 metric-specific raw ceilings are
  +0.030948 fidelity, +0.025908 reach, and +0.000579 Jaccard; they are
  mutually incompatible optima and cannot be added. Panel work remains a side
  gain, not the route to first.
- **Experiment 128:** the clean replacement for 124 is a mechanical PASS.
  Across 14,400 A/B/C target/dose/library evaluations it used the exact
  `log1p(50000*p)` source units and preserved every non-panel integer total,
  target-self, mutable-panel budget, grand emitted total, and scorer-normalized
  non-panel linear/log/effect value bit-identically. Production is still
  unauthorized pending biological confirmation.
- **Experiment 114:** terminal NO-GO. On its already-authorized 35-target
  HEK293T PDS diagnostic, the frozen STRING+DepMap primary neighbor replacement
  changed raw PDS by `-0.120168`; STRING-only and DepMap-only were also strongly
  negative. Do not open the H1 surface.
- **Experiment 115 / 129:** the verified GSE295214 H1 compact donor has 2,354
  perturbations, 16,290 covered canonical genes, and 64 VCC targets, 57 new
  relative to the old H1 donor. It is not an untouched source: its direct
  residual was already tested in Experiment 117 and its 64 panel rows were the
  Experiment 118 Pan confirmation. Experiment 129 may proceed only as a
  materially different source-only learner on non-panel splits; it must not
  retune those closed residual, mask, dose, or ranking choices.
- **Fresh panel confirmation census:** no public source simultaneously clears
  40 exact VCC targets, broad response coverage, and complete independence
  from consumed outcomes. The best actionable compromise is IGVF TeloHAEC:
  65 VCC targets and 17,464 response symbols. Its broad direct-donor route
  failed Experiment 105, but its panel coordinates remain unopened and can
  support a separately frozen panel-response confirmation for a future
  source-only candidate. HCASMC's 27 targets are a secondary confirmation.

No challenge submission, count build, commit, or push was made.

---

## 2026-08-30 20:50 UTC — a15 promoted; Experiments 138 and 139 closed

### a15 — new best

`a15-cd4f010.vcc` published as entry `KpoZt5Rie0EjgMaZQECE` under
`The Fern Filed An Appeal`.

- overall: **`0.1740865087`**, rank 9 of 482;
- delta versus a14: **`+0.0036622644`**;
- scaled PDS `+0.007256`, MSE `+0.009441`, NMAE `+0.002266`,
  fidelity `+0.002190`, reach `+0.000826`;
- Jaccard changed by only `-0.0000066`.

All four preregistered a15 falsifiers passed. Full audit:
`experiments/096-mse-floor/A15-BOARD-RESULT.md`.

### Experiment 138 — DLD1 replicated sparse-sign gate: terminal NO-GO

The final arm was `k=200, beta=0.05`, but nested macro gains were only
`+0.0000995142` PDS and `+0.0000602495` signed credit. Both were statistically
above the frozen nulls, but 30-50 times below the `+0.003` practical gate.
RPE1 PDS was negative in two folds and PDS family concentration was `0.508929`.
The 410-target test remains sealed.

### Experiment 139 — high-MOI Jurkat donor-simplex routing: terminal NO-GO

Every development fold and the complete 613-target development set selected
exact identity. The 138-landmark high-MOI signal could not improve convex
weighting among the five allowed same-target source contexts. The 410-target
test remains sealed.

---

## 2026-08-31 18:45 UTC — Experiment 142 reliability weighting: terminal NO-GO

The corrected, hash-bound evaluator opened only the authorized 613-target
development split. All five nested folds and the complete development fit
selected exact identity (`eta=0`), so nested and fixed raw-PDS and
signed-credit deltas were all exactly zero. Only 2 of 59 promotion checks
passed.

The fixed nonidentity reliability arms showed a weak direction/PDS trade:

- `eta=0.25`: PDS `-0.000020`, signed credit `+0.001183`;
- `eta=0.50`: PDS `-0.000450`, signed credit `+0.001855`.

Non-promotable support controls at `eta=0.5` were positive but subthreshold:
cell count `+0.001961` PDS, total UMI `+0.001996`, and occupied-batch count
`+0.001917`. Their largest PDS gains were concentrated in Jurkat.
Pooled-effect norm selected identity.

The 410-target split remains sealed. Do not retune, rescue, integrate, build,
score, or submit Experiment 142. Full hash-bound result record:
`experiments/142-batch-reliability-weighting/RESULTS.md`.

---

## 2026-08-31 — Experiment 143 normrel replication: terminal NO-GO

The fixed leakage-safe `normrel613^-1.5` candidate was evaluated once on the
authorized 613-target development panel after a no-truth control-plane audit.
It used Experiment 142's exact PerturbMatch-aware equal-weight baseline and
the same 613 rows and 4,992 genes for both the relative-norm numerator and
denominator.

The candidate showed a small real retrieval signal but failed the frozen
promotion gate:

- macro raw-PDS delta: **`+0.00279128`**, bootstrap 95%
  `[+0.00112214, +0.00458825]`;
- macro signed-credit delta: **`-0.00093903`**, bootstrap 95%
  `[-0.00165080, -0.00021027]`;
- family PDS: K562 `-0.000529`, RPE1 `+0.002026`, Jurkat `+0.006877`;
- family signed credit was negative in all three families;
- PDS positive-gain family concentration was `0.77246`, above the `0.5` cap;
- 3 of 12 frozen checks passed;
- decision: **`TERMINAL_NO_GO`**.

The 410-target values were not loaded or interpreted; their containing donor
files were streamed bytewise only for preregistered SHA-256 verification
after the access lock. Reserved confirmation expression remains sealed. No
CD4 production-frame audit, count generation, `cell-eval2`, challenge build,
submission, commit, or push occurred. Do not retune, subset, blend, confirm,
or integrate this normrel family. Full record:
`experiments/143-normrel-production-replication/RESULTS.md`.

---

## 2026-08-31 — Experiment 145 spectral trace-bin attenuation: terminal NO-GO

The first pre-outcome implementation was withheld after an independent audit
found that its proposed development validator would hash complete donor
containers containing sealed 410-target rows. A corrected v2 package preserved
v1, reproduced every candidate and null array exactly, switched to metadata-
only container validation plus selected-row loading after the access lock,
enforced an 8 GiB process-tree limit, unified observed/bootstrap selection,
and passed 20 tests plus a second independent audit.

The one authorized 613-target evaluation showed a clean direction/retrieval
trade:

- `beta=0.25`: PDS `-0.00005109`, signed credit `+0.00986897`;
- `beta=0.50`: PDS `-0.00081566`, signed credit `+0.01920282`;
- `beta=0.75`: PDS `-0.00237279`, signed credit `+0.02740058`;
- `beta=1.00`: PDS `-0.00457143`, signed credit `+0.03406376`.

Nested selection chose identity in four of five folds. Nested macro PDS was
`-0.0000684160`, signed credit was `+0.0020382723`, the final selector chose
`beta=0`, and only 5 of 17 promotion checks passed. The low-variance tail
therefore carries retrieval signal even though attenuating it improves broad
signed direction.

Decision: **`TERMINAL_NO_GO`**. Do not retune, route, run the frozen nulls,
open the 410-target split, integrate, build, score, or submit this family.
Full record: `experiments/145-spectral-trace-bins/RESULTS.md`.

No challenge submission, count build, `cell-eval2`, commit, or push occurred.

---

## 2026-09-02 — a23 landed `-0.00073`; three research questions closed; the instrument moved

**Board: a23 `0.1733597079`, rank 25 of 578** (entry `IzQVfDpeDWufkPT8Lmd2`,
`Reviewer 2 Was Right About Everything`). a15 remains our best at `0.1740865087`. One Sep-2
slot remains and is deliberately unspent — see "why" below.

### a23: the profile channel works, `reach` charges a toll

`pds` **+0.001737** (above the optimistic prediction) and `mse` **+0.001525** both positive —
**the first time this project has moved `pds` and `mse` together with one profile change.**
`fid` and `jac` also correctly signed. `reach` **−0.007933** carried 182% of the net move.

Mechanism: the minimum-L2 mass correction is `eps ∝ −exp(b_ctrl)`, so it reduces the
highest-expressed genes — which carry the most counts, hence the smallest emitted p-values,
hence they populate exactly the head of the confidence order `reach` reads. **Cheapest in norm,
most expensive in ordering.** Three-arm law, one mechanism: a21 `−0.0185`, a22 `−0.0107`,
a23 `−0.0079`. **Price a ~0.008 `reach` toll into any change to the emitted profile.**

Also: the harness-cosine → board-`mse` transfer model **over-predicts by 1.7-7x**. A
+13.6%/+56.0% harness cosine gain bought a +1.2% relative board `mse` gain. Do not price a
profile arm by scaling a harness cosine ratio. Full record:
`experiments/156-delivered-geometry/A23-BOARD-RESULT.md`.

### The follow-up, and why the second slot is NOT being spent on it tonight

`weighted.py` measured seven correction directions and found the profile gain **insensitive** to
placement, and I wrote "there is no placement to tune." The board says placement is irrelevant
for the profile and decisive for the ordering. So: put `eps` on gated genes **outside the reach
head**, keep `pds` +0.0017 and `mse` +0.0015, drop the toll. Expected `+0.0005` to `+0.0018`.
**It must be priced on `reach` first, on `h1-test`, using the pre-computable emitted p-order
(`|Fbar_g − 1/2| / sqrt(tc_g)`, Spearman +0.996 with the realised order).** Shipping it unpriced
would repeat a19's and a23's exact mistake.

### 157 — `reach` ordering: NO-GO, and one real artefact

The exp-094 inversion is **gone** (it was a depth-confounded-harness artefact) and exp 154's
closure **survives** on a better surface with a better instrument. No buildable key beats the
shipped emitted order by enough to matter. **Do not spend a slot on reach ordering.**

But one finding is challenge-regime and not harness-limited: the shipped tier policy
(`--tail-spread 1.0 --keep-top 8000 --overdispersion 4.0`) puts `spread 1.0 / phi 4` on ~1,900
genes, where the gamma-Poisson sends most of a low-expression gene's 400 cells to zero. Measured
in closed form on the real challenge controls: that block is significant at **zero** predicted
effect for **73.9-78.9%** of its genes, **100% called DOWN**, median `|z| = 7.4` against the keep
block's 1.7 — so **~1,900 signal-free genes are placed at the very front of the reach ranking by
construction**, and 79.7% of the top-10 comes from them. The shipped call set is **28.8-31.4%**
zero-effect artefact, not the ~19% recorded for a8's flat law. `spread 1.0 / phi 0` on the tail
block would take the mixture to 14.6-20.7%. **Caveat: the board does not corroborate the size** —
the a11→a12 residual after the amplitude model is ≈ 0 against HCT116's −0.016, and HCT116
overstates it 4.4x against h1-test. Worth one properly-priced arm, not a rushed one.

### 158 — per-target amplitude: NO-GO, and the reason is durable

A build-time key **does** predict the mse-optimal per-target scale: production's own
`donor_cosine_key` restricted to scored-and-gated coordinates, Spearman **+0.395 / +0.465** to
`s_p*`, delivering **+16.2% / +29.6%** aggregate cosine strictly cross-surface at exactly zero
`pds` cost, cluster bootstrap positive, permutation null `p <= 0.0005` on both surfaces. It is
still a NO-GO: the deployable prize is **+0.0003 to +0.0029** overall (10-30x below the oracle
prior), and **`nmae` pays more than `mse` collects** — the recommended rule costs
`−0.0023 / −0.0080` on `nmae`, and the *oracle* scale costs `−0.0245`. **The two members
disagree about the per-target axis itself**, so no better key rescues it. This is rule 2 in its
original form: the metric the parameter serves is the one that cannot move. Per-target magnitude
joins the closed list. Note `ORACLE SIGN ONLY` is 17%/38% of the oracle gain, so the lever is
mostly relative confidence among correctly-signed rows, not the anti-correlated quarter.

### 159 — the paradigm is at ~30% of its ceiling, not 90%, but the headroom is expensive

Disattenuated, only ~30% of a target's target-specific response energy is shared across cell
lines (`kappa^2 ≈ 0.34`); ~70% is line-specific. **Production sits at 31.1% / 28.5% of the
infinite-donor context-independent ceiling** on the actual challenge panel. The ceiling is
approached as `1/k` in the number of *independent* cell lines, so reaching even half of it needs
about **ten** lines against the four we have. That reopens donor expansion as the main long-term
lever — `AGENTS.md` §12 closed it on *coverage* grounds (three prepared donors carry zero panel
targets), which is a different argument from *marginal value per independent line*.

**Verify before acting on the magnitudes.** 159's headline that the thinned Orion harnesses put
"86.5% / 86.9%" of their true effect energy in a target-independent vector against h1-test's
"5.4%" did **not** reproduce on my own re-measurement: I get **26.4% / 44.1% / 3.6%**. The
ordering and its consequence are confirmed and important; the magnitudes are 2-3x overstated.
The `1/k` extrapolation rests on the same code path and should be re-derived before ten cell
lines are acquired on the strength of it.

### The instrument finding, which outlasts all three

`h1-test` is the only harness validated against the board on `reach` (reproduces two board
points to within 0.014 raw; repaired HCT116 misses a15's by **2.6x**), and the Orion thinned
pair carries **26.4% / 44.1%** of its truth's effect energy in a single target-independent
vector against h1-test's **3.6%**. Two agents found this independently on different members.
A large shared component in the *reference* inflates decoy similarity, so absolute `pds` and
absolute transfer cosine on the Orion pair are depressed and not board-comparable.
**Prefer `h1-test`; use the Orion pair for relative direction only; never quote an absolute
from it.** Warning added at `AGENTS.md` §6.4.

### Corrections and traps from the three closures — read before the next screen

These came out of the closures' final reports rather than their result files, and three of them
invalidate arguments this repo has already made.

1. **`pds` is NOT exactly invariant to the shipped per-perturbation scale.** `docs/CURRENT_STATE`
   has recorded that "`pds` is exactly neutral to `gamma_pert`, structurally, because
   `pds_cosine` normalises each row". That is true only of a **post-render** scale. Production's
   `--shrink-pert-key` scales **pre-render**, and measured end-to-end it moves pds by
   `+4.05e-04 / -5.78e-04` — **sign-inconsistent across the two surfaces**. Verified numerically
   at `+0.00e+00` for a post-render scale, so the mechanism is the renderer, which is the same
   object experiment 156 is about. Any past attribution that leaned on "pds cannot move, so the
   movement attributes entirely to `gamma_gene`" needs re-reading.

   **And a23's mass fix does NOT repair it — I guessed it would and measured that it does not.**
   Applying a random per-row scale (sd 0.35-0.42) and then re-imposing the mass constraint moves
   pds by `-5.8e-05 / -5.2e-04`, the same order as the unfixed `+1.2e-04 / +2.9e-04`. The reason
   is the reason scaling cannot satisfy the constraint in the first place: the feasible set is
   **curved**, so moving along a ray and projecting back onto it changes direction. Since no
   post-render scale is feasible either, **an exactly pds-neutral per-target amplitude knob does
   not exist anywhere in this pipeline.** The structural claim should be stated as a property of
   the *unrendered* effect, which is not a thing we can emit.
2. **`debench.Surface(n_real_target=N)` must not be used to power-match `reach`.** Uniform
   `n_real` deletes the small-`n_real` tercile that carries **98%** of our reach — it reads the
   shipped arm at `0.008` instead of `0.183` for a 2.7x pool change. Scale each perturbation's
   pool proportionally. Exp 154's setup used the uniform form on the reach-blind surface; its
   conclusion survives the repair, but the instrument is wrong for anyone reusing it.
3. **`h1_2025` leaks 85/85 targets into h1-test.** See the AGENTS.md §6.4 warning. Audit
   donor↔surface overlap before every screen.
4. **A single-surface key screen will select a sign-flipping key and ship a regression.** Exp
   104's `score_agreement_mean` is the *strongest* per-target key on either Orion surface alone
   (|rho| 0.405 / 0.441) and **its sign flips between them**; used cross-surface it is the two
   worst rows in exp 158's table. Never select a per-target key on one surface.
5. **The reachable per-target ceiling is the whole oracle, not a fraction of it.** The
   positive-constrained oracle `max(s_p*, 0)` recovers **98% / 91%** of the full oracle gain, and
   sign-only recovers just 17% / 38%. So the lever was never "find the anti-correlated quarter";
   it is relative confidence among the correctly-signed rows. It is still a NO-GO on `nmae`.
6. **The tail-block artefact must be priced on `mse` before anyone touches it.** The
   `spread 1.0 / phi 4` block exists to feed `mse`'s jackknife credit (096), so removing its
   zero-effect call artefact (73.9-78.9% of ~1,900 genes, 100% DOWN) trades against the member
   it was built for. `spread 1.0 / phi 0` would take the shipped mixture's zero-effect call rate
   from 28.8-31.4% to 14.6-20.7%. Also note the oracle reach is exactly `sign accuracy / 0.9`
   and is invariant in pool size — the cleanest available statement that reach is model-limited.

### The next research object, named concretely

Exp 159's line-pair matrix says transfer is **strongly pair-dependent**, which no previous
experiment exploited: same-line control K562-gwps↔K562-essential `+0.2712`; K562↔RPE1 `+0.1164`;
Jurkat↔K562 `+0.1510`; but Orion HCT116↔HEK293T only `+0.0300` and H1↔HCT116 `+0.0332`. Our two
Orion donors are the *worst-transferring* pair in the matrix.

And the decisive asymmetry: knowing the **context** (same-line prediction, cos `0.2627`) is worth
**~3x** knowing the **target** (best cross-line single donor, `0.0897`). Production is worth 3.1
effective donors and 89% of the best 4-donor consensus buildable today; the leader is at 10.3.

> **All ~40 NO-GO experiments tried to improve the CONSENSUS. None attempted the RESIDUAL.**
> The object to build is a cross-line residual predictor: given the basal profile, the target,
> and the consensus effect, predict the **line-specific** part — which is ~70% of the energy and
> is exactly what a context-independent effect matrix cannot express by construction.

Caveat carried forward: 159's absolute magnitudes need independent re-derivation (its harness-blob
headline overstated by 2-3x against my own re-measurement), so treat the `1/k` extrapolation and
the "ten cell lines" figure as unverified. The pair-dependence and the context-vs-target
asymmetry are the parts worth acting on.

### 160 — the control cells cannot supply the context: NO-GO, and it is the fifth route into this wall

Experiment 159 said knowing the CONTEXT is worth ~3x knowing the TARGET, and our effect matrix
is context-independent, so the obvious move is to read the context's response geometry off its
own 18,400 control cells. **Containment is strongly positive and the payoff is negative.**

A line's own control basis captures **8.0x / 10.9x** (k=5) more of its target-specific response
residual than *another line's* control basis, and **153x / 339x** more than a random subspace
with the same per-gene sd envelope. Robust to regressing out library size and mitochondrial
fraction (Δ ≤ 0.0007) and reproducible across disjoint halves of the control cells.

But every way of using it loses `pds`: hard projection **−0.18 to −0.27** (a ≤200-dim projection
makes 132 predictions mutually similar, which is hubness by construction), and blend/amplify buy
+3 to +10% cosine for **−0.008 to −0.017** pds — roughly 5x worse than a23's measured board
exchange rate. Dropping PC1 and orthogonalising the basis against the panel mean do not rescue
it.

**The reason is one number:** on HEK293T, dropping PC1 alone takes residual capture from 0.1698
to **0.0525**. The context-specific capture is essentially a single direction, and that direction
is the line's *mean* response — the shared programme, which is `pds`'s zero and the hub whitening
exists to suppress. The own basis captures the shared mean at 0.71/0.89 against the residual's
0.11/0.17.

**The target-specific component of a perturbation response is not recoverable from unperturbed
variation.** That is a statement about biology, and it now agrees with 155 (per-gene regression,
split-half r 0.028), 067 (variance transport, −0.00025), 017 (basal modulation, +0.003 and
harmful by β=1) and Heisenberg's live board arm (added context-shared co-expression, pds
0.4509→0.3974). Five routes, one wall. Full record: `experiments/160-control-subspace/RESULTS.md`.

**So the line-specific 70% must be learned from donor line PAIRS, not from control geometry.**
That is experiment 161, and the unlock is that donor expansion was closed on *panel coverage at
inference* — `rpe1`, `nadig_jurkat` and `k562_essential` cover zero of the 300 panel targets but
carry thousands of shared non-panel targets, which is all a residual model needs to TRAIN on.

### Where the next slot should go — `fid` is precision, and a third of our call set is noise

`n_pred` (~4,400) far exceeds the challenge's `n_real` (~240, settled by a21), so the coverage
term is saturated and **`fid` reduces to precision on our own call set**. Board raw `fid`: ours
**0.5022**, board maximum 0.5175, leader cqawesome **0.4877** — *below* us. **Nobody in this
competition calls a DE gene's direction better than a coin flip.**

157 §5b measured in closed form on the real challenge controls (so NOT rule-9 limited) that the
shipped tier policy calls **28.8-31.4%** of gated genes significant at **zero** predicted effect,
100% DOWN, median `|z|` 7.4. If a third of the call set carries no model signal it dilutes
precision toward 0.5. Removing it should raise `fid` (precision), `reach` (those genes are 79.7%
of the top ten in the ordering) and `jac` (they are false positives) at once — three members
where we rank 179th-207th of 578. Precision 0.60 would be worth **+0.049 overall on `fid` alone**.

The binding constraint: the tail block exists to feed `mse`'s jackknife credit (096), so it must
be priced on `mse` too or it is a19 again. Experiment 164 is pricing it with a8→a15 as the
known-answer row.

---

## 2026-09-02 — THE COMPOSED SLOT GATE (use this before every future submission)

Experiment 162 built a panel-restricted leave-one-library-out validator (PLOLO) and measured
what it can and cannot predict, against **12 of our own board points** transcribed from primary
sources with the per-member decomposition asserted to sum to the overall delta. The answer is
narrow but usable, and it composes with two other instruments into a gate that **retrodicts the
a23 veto out of sample**.

### What each instrument can actually price

| channel | instrument | measured accuracy | verdict |
| --- | --- | --- | --- |
| `pds` | PLOLO, **challenge frame**, calibrated scale ≈ 1.011 | LOO σ **0.00149** raw pds = **±0.00025** overall; sign **6/6**; survives every fold drop | **TRUST** |
| `mse` | 096's `claim` machinery (`mech_claim_grid2.py`) | 3/3 known-answer `claim` gates PASS (0.11860 exact, +1.16%, +1.06%) | **TRUST for emission/dispersion arms** |
| `mse` | PLOLO | LOO **0.11239** | **USELESS** |
| `nmae` | PLOLO | LOO 0.00511, sign **3/6**, fails predict-zero null at every pool size | **UNPRICEABLE** |
| `reach` | PLOLO | LOO 0.00784, sign **2/6** (worse than chance), fails its null | **UNPRICEABLE** |
| `reach` | empirical emitted-profile toll | **−0.0146** from a21+a22; a23 realised **−0.00793** | **CHARGE IT AS A PRIOR** |
| `fid`, `jac` | — | carry `n_real` (rule 9); no reference DE pool exists on a donor library | **UNPRICEABLE** |

PLOLO cannot produce a `reach` number **at all**, by construction: a donor library yields
pseudobulk effects and never cells, so there is no reference DE pool. That is precisely the
member that killed a23.

### The numbers that set the bar

- Largest priced move PLOLO has ever produced: **+0.000997** on overall.
- Realised RMS of the **unpriced** channels (overall minus pds/6) across the six geometry arms:
  **0.00735**. **Unpriced/priced ratio 7.4×.**
- And the split that makes the gate usable:
  - **amplitude-changing arms** (a11, a12, a14, a19): unpriced RMS **0.00890**
  - **shape-only arms** (a15, a23): unpriced RMS **0.00188**

### The gate

```
predicted_overall  =  PLOLO Δpds / 6 · 1.011
                    + 096 Δclaim→Δmse / 6              (if amplitude / spread / phi / dispersion moves)
                    − 0.0146 / 6                       (if the EMITTED PROFILE or emission law changes)

SHIP only if  predicted_overall  >  unpriced RMS for the arm's class
              i.e.  > +0.0019  (shape-only)      or   > +0.0089  (amplitude-changing)
```

**Consequence, and it should change behaviour: stop shipping small arms.** An
amplitude-changing arm must predict better than **+0.0089** overall to be worth a slot. a19
predicted +0.002 and lost 0.01352 — it was never inside its own noise floor. Every arm since a15
has had a predicted gain below its class threshold.

### The gate is validated on the arm it would have saved

a23's veto, computed with a23 **held out** of every fit:

```
emitted-profile reach toll, from a21+a22 only        −0.014569
PLOLO Δpds for a23, scale fitted WITHOUT a23         +0.004060 raw → +0.000677 overall
composite prediction                                 −0.001752
realised                                             −0.000727      → VETOED, correctly
```

### Standing consequences

1. **Any arm touching the emitted profile is charged −0.0146 on `reach` unless it is priced on
   `h1-test` directly.** h1-test is the only surface that reproduces board `reach` (two points to
   within 0.014 raw; repaired HCT116 misses a15's by 2.6×).
2. **PLOLO must be read in the CHALLENGE frame, not the library frame.** Library frame
   over-predicts Δpds by 14× against the challenge frame's 5.5× (and the calibrated scale absorbs
   the rest). Challenge-frame render residual 2.88/2.90/2.87% independently reproduces the
   2.83–3.35% measured in `A23-PREDICTION.md`.
3. **`sgn_w` is unvalidated** — it fell 0.492 → 0.417 under a23 while board `nmae` was flat. Do
   not use it as an `nmae` proxy.
4. `mse` is priceable only through 096's `claim` route, and only for arms that move the emitted
   dispersion. A whitener/profile-shape change has **no** local `mse` instrument — which is
   exactly how a19 lost −0.0959 on that member.

Full record: `experiments/162-retrodicting-validator/RESULTS.md`.

## 2026-09-02 — 165 CLOSED NO-GO: no placement of the mass correction protects `reach`, and h1-test is blind to a23's toll

Head-protected variants of a23's mass fix (`eps = λ·u·w`, 14 placements × 2 seeds, 6 × 10-seed
confirmation, four instrument gates passed including exact reproduction of a8 0.20291 / a15 0.18283):
**no placement moves h1-test `reach` outside its error bar, and protection costs the profile gain
monotonically** (aggregate cosine +0.89% → −0.29% from `minnorm` to K = 9000).

**The instrument finding matters more than the arm.** `minnorm` — a23 itself — reads **+0.00041 ±
0.00046** raw h1-test `reach` where the board realised **−0.00793 scaled**. **h1-test reproduces board
`reach` for emission-law changes (a8/a15, spread policy only) but NOT for profile changes.** Rule 3
fails for `reach` on the mass-fix family. Conversely h1-test retro-predicts a23's *profile* channel to
1.5× against 12–48× for Orion, so it is the right surface for `mse`/cosine; it is only
`reach`-on-profile that is blind. Consequence for 171: its spread-schedule arms are inside the
validated regime **only while the pseudobulk stays bit-identical**.

## 2026-09-02 — the Orion pair shares one artefactual axis (cos 0.872) that K562 lacks; deflating it per donor helps

`175/deflate175.py`: HEK293T's top principal direction holds **19.7%** of its variance (K562 5.6%),
and HCT116's and HEK293T's top directions have **cosine 0.872** with each other against 0.23–0.31 with
K562's — two tissues from one atlas sharing their dominant axis is a **platform axis**, the
donor-table face of the 26–44% harness blob. Removing exactly one component per donor before the
consensus lifts H1 sign accuracy for every donor and for the consensus (**0.5697 → 0.5770**); two or
more removes signal. Whitening down-weights this direction but spends budget on it rather than
removing it. Handed to 178 as a fourth switch — it removes noise rather than trading `pds` for
accuracy, so unlike the other three riders its `pds` cost may be ~0. Distinct from 156's closure
(post-prediction deflation). **Identified:** the axis is essentially orthogonal to H1's truth (cos
0.05–0.13) and its top-30 loading genes — AUTS2, GPC6, NRXN3, PRKG1, WWOX, RAD51B, RORA, ERC2,
CACNA2D3, PARD3B, ADGRL3, FAT3, MAGI1, PTPRG … — are all very long intron-rich genes, all loading
negatively: **the nuclear-RNA-fraction / long-gene technical axis of single-cell data**, a property of
the X-Atlas libraries and not a response to any target. Production's whitener is spending budget
flattening a sequencing artefact. **Resolved:** Arc's controls carry the axis as within-population variation (cos 0.31–0.43 with
their PC1) but Arc's perturbation effects do not load on it (cos ≤ 0.13), while the **Orion harness
truth does — 36.3% of HEK293T's and 13.8% of HCT116's effect energy sit on it.** That is the mechanism
of the harness warning and of 168d's instrument rule: the Orion surfaces reward predicting a sequencing
artefact. **Production's whitener already removes it from the shipped effect (0.11% of energy vs 9.8%
in the raw consensus)** — so the deflation switch cannot move the profile directly; it acts only through
the whitener's refit and through any component that reads the raw consensus (sign repair: deflated
`cons` is +0.012 more accurate). Recorded permanently in AGENTS.md's harness warning.

## 2026-09-02 — 170 (the agent survived the kill and finished): `WHITEN_GAMMA = 0.65` is a LOCAL OPTIMUM. TERMINAL NO-GO on gamma.

Full production-staged sweep (`build_v9_pre_whiten_effects → apply_v9_whitening → … → emit_counts`
by calling the shipped functions), bit-identical known-answer gates (a15's `baseline_pre` /
`candidate_pre` sha256; a15's and a19's post-whitening audits to 1e-15). **Raw `mse` on the
challenge panel is V-shaped with its vertex at 0.65**: 0.926 (0.35) · 0.881 (0.60) · **0.874 (0.65)**
· 0.880 (0.70) · 0.893 (0.80). **Every other gamma in [0.35, 0.80] is net negative on all six
members, −0.0016 to −0.0071.**

**Mechanism — the ρ cliff.** Lowering gamma raises the effect's mean-square `A` by 14% (0.35 → 0.65)
and with it the across-perturbation budget `B`; a15 sits essentially on the cliff (ρ 0.99 / 0.94 /
0.88 per context); pushing `B` past `claim` saturates `min(claim, B)` and the extra amplitude is not
refunded. **097's h1-test `mse` estimate assumed the arm sits ON the cliff, and h1-test cannot
price `mse` at all (its raw `mse` is 0.96–1.02, above the clamp).** 097's "+0.0027 at gamma 0.30"
reads −0.005 to −0.007 once the `min()` is evaluated per context. My "sub-threshold live" reading
above inherited 097's error and is withdrawn. a19 re-attributed: gamma 21%, `lam` 79%.

**The design rule that falls out, and it is the most useful sentence of the day:** lowering gamma
does buy alignment (`d(log X)/d(gamma) = −0.576`, R² 0.974 — +16% X over [0.65, 0.35]) but pays in
amplitude the metric refuses to refund. **A transform that raised `X` without raising `A` — norm-
preserving, not norm-adding — would convert the full gain (≈ +0.0055 overall on `mse`, more at the
leader's cosine).** Sign repair (173) is exactly that shape: flips preserve `|H|`, cos +7.8%. Gamma
is not. 178's composite corrected accordingly: gamma dropped; sign repair (deflated source) +
K562-heavy weights + per-donor deflation, with `A`/`B`/ρ per context reported for every arm.

Side result: `pds` has a small seed-free interior maximum at gamma 0.50–0.55 on the *delivered*
object (+0.00046), so 097's "pds monotone in gamma" described the intended effect, not the delivered
one. Not actionable.

## 2026-09-02 21:05 — 178 CLOSED: the rider family is dead (gains compose sub-additively, costs do not)

Every rider and composite is **negative on all six, −0.005 to −0.030**, against the +0.0089 gate.
`pds` costs compose sub-additively (0.52–0.81) but `mse` costs compose additively-to-super-additively
(0.96–1.23) and DE gains sub-additively (0.29–0.84): **gains discount when composed, costs do not.**
Every rider except sign repair raises `A` 11–28% and pushes ρ to 1.000 in all three contexts (170's
cliff, three more measurements). Sign repair is the one norm-preserving rider (`A` unchanged, `mse`
**+0.0021** — 170's criterion confirmed) but on the board-validated PLOLO instrument its `pds` cost is
**−0.0245 raw, 6× what h1-test read** (173: −0.0042). Donor weighting: whitening consumes the sign
gain (+0.0043 → −0.0012 after the operator). Builder switches `VCC_DONOR_WEIGHTS` / `VCC_SIGN_REPAIR`
kept, default = shipped, a23's effects audit reproduced bit-for-bit (9/9 fields), 46 tests pass.
**Lesson for the instrument table: h1-test understates `pds` costs of un-whitening arms ~6×; PLOLO is
the `pds` instrument.** `experiments/178-composed-riders/`.

## 2026-09-02 20:46 — 177: THE ORDERING DEFECT IS SOLVED BY THE EMISSION LAW. 81–110% of the `|H|` counterfactual realised, `pds`/`mse` pinned. The arm's sign depends on the challenge `n_real`.

**The zero-effect null, on the three challenge contexts at production depth through the unmodified
production `emit_counts`:** a15's shipped law calls **4,302 / 4,312 / 3,883 genes (43.3% / 44.8% /
38.4%) significant at ZERO predicted effect**; the resampling law (`mode="resample_mix"`, real
control cells × fold) calls **1 / 18 / 19 (0.0–0.2%)**. Eleven known-answer gates pass, four of them
bit-identical fresh reproductions.

**With no schedule at all, resampling reverses the falsifier:** Spearman(−log p, ·) inside our call
set goes from a15's **+0.485 expression / −0.391 `|H|`** to **+0.001 / +0.601** (`hy_kt2000`).
h1-test reach, twin-averaged (resample twins are 5–60× tighter than the quasi pair):

| arm | board `n_pred` | Δ reach raw | % of `|H|` ceiling | precision | nmae |
|---|---:|---:|---:|---:|---:|
| a15 (`base_mf`) | 4,424 | — | 0 | 0.5079 | 0.9664 |
| `rs_all` (everything resampled) | 812 | **+0.0413** | 81% | 0.5184 | 0.9618 |
| `hy_kt4000` (top-4000 by `|H|` parametric, rest resampled) | 1,683 | **+0.0559** | **110%** | 0.5132 | 0.9704 |
| **`hy_kt8000` = a15 with its poisoned tail replaced** | **2,410** | **+0.0352** | 69% | 0.5131 | **0.9662** |

`pds` is an **identity**: one sha256 of the delivered per-gene panel totals is identical across all
twelve real-depth arms (85/85 targets). `mse` is exactly 0 through 096's claim route (every arm above
the knee). Five experiments (093, 095, 157, 171, the a15 tier) tried to deliver this ordering through
a spread schedule and 171 closed the family; **the emission LAW delivers it outright.** Also refuted
with a new mechanism: adding symmetric noise to a resampled sample *inflates* `|z|` — a rank-sum
has no variance channel in the needed direction on any emitter (`Var(F̄) > 1/12N` under
over-dispersion). The permutation control beats the schedule, as in 171.

**Where it pays: `n_pred` falls 4,424 → 812–2,410, and `fid`'s yield cap `min(1, n_pred/n_real)` plus
`jac`'s union now decide the arm.** Sum of the four DE members in overall units:

| arm | n_real 191 | **240** | 400 | 812 | 1,200 | **2,670** |
|---|---:|---:|---:|---:|---:|---:|
| `rs_all` | +0.0197 | **+0.0211** | +0.0239 | +0.0217 | −0.0807 | −0.2375 |
| `hy_kt4000` | +0.0075 | **+0.0076** | +0.0073 | +0.0030 | −0.0696 | −0.2479 |
| **`hy_kt8000`** | +0.0060 | **+0.0060** | +0.0056 | +0.0026 | −0.0158 | −0.1725 |

**What we know about the challenge `n_real`, four independent lines:** (1) board algebra bounds it
to [191, 4,716]; (2) 153's enrichment meeting-point and a21 read ~240; (3) h1-test — Arc's own 2025
data, same chemistry, the challenge's closest regime — has median `n_real` 222; (4) **our own board
`jac` (0.0301 at `n_pred` 4,424) implies recall 0.57–0.71 if `n_real` is 191–240 — matching h1-test's
0.67 — and recall 0.078 if `n_real` is 2,670.** 093's 2,670 is flagged circular by rule 9 and would
require our recall to be 9× worse on the challenge than on Arc's own data. **The weight of evidence is
~240, where every resample arm is positive.**

**Rule 13, applied literally, vetoes** (−0.00243: it charges a flat −0.0146 for touching the emission
law and treats reach/fid/jac as unpriceable). **The honest counter-argument, which the agent and the
coordinator both hold:** rule 13's "reach 2/6 signs" was fitted on arms whose reach deltas sat inside
the 0.0075–0.016 seed spread; these are +0.035–0.056, 3–7× that spread and 200× the arm's own twin
spread, on the one surface validated for reach on emission-law changes (157 §1, 165). This is a
different resolution regime from the one the gate was calibrated on. **And a22 — resampling the TOP
10% (the opposite assignment) — lost −0.0107 on reach, which is what this mechanism predicts for the
wrong-way arm.**

**CORRECTION (22:40) — the +0.006 was a single-`n_real` reconstruction; on h1-test's REAL `n_real`
distribution `debench` reads `hy_kt8000` at ≈ −0.005.** `fid` 0.5079 → 0.5061 (−0.0010 overall), **`jac`
0.0683 → 0.0541 (−0.0064)**, reach +0.027 (+0.0020), `nmae` ≈ 0. Recall falls 0.66 → 0.39. Decomposed by
stratum, **the entire `jac` loss is in the top tercile (`n_real` 456–5,484: −0.042); for the 61/85
perturbations with `n_real` ≤ 800 the two arms match to +0.0003.** Board calibration: our board `jac`
is 0.0301 against h1-test's 0.0683 for the same a15, between h1-test's "≤ 800 only" (0.0133) and its
full value — the challenge's heavy tail is roughly half as prevalent as h1-test's. Scaling the `jac`
loss by that puts **`hy_kt8000` at ≈ −0.0006 overall: a wash on expected score.** The manifest confirms
the regime (`ground_truth_cells` 138,400 = ~461 cells/perturbation vs h1-test's 400; controls 18,400 in
both). **What the slot buys is information, not score**: its board `fid`/`jac` read the call-set regime
directly, which gates whether the resample family (K between 4,000 and 8,000, or `rs_all`) can ever
be tuned for D/E/F. Honest range for `hy_kt8000`: −0.005 to +0.006. Not a score arm.

**179 BUILT AND VERIFIED IT (23:30). Package `submissions/hy-kt8000.h5ad`, sha256 `e2902589…`,
19.7 GB, 51.8% of the density cap, `vcc prep --dry-run` PASS, NOT submitted.** Rule 10: 78/78 effects
audit fields identical to a23 incl. `effects.final_sha256`; config diff exactly `emitter` +
`resample_keep_top`. **`pds`/`mse` pinned by identity on the two real artifacts: 900/900 delivered
pseudobulk blocks identical.** The resampled tail's zero-effect false-call rate collapses 33–49% →
0.06–1.9% on the built path. **The one surprise: the board `n_pred` cut is −21% (4,474 → 3,528), not
the −46% h1-test forecast** — the challenge's tail contributes less to our calls — so every DE delta
scales by 0.46 and the `fid` risk window narrows to `n_real ∈ (3,528, 4,474]`, which every estimate
sits below. **Six members: `pds` 0, `mse` 0, `nmae` ≈ 0, `reach` +0.0009, `fid` −0.0005 (or +0.0029
if `n_real` ≤ 3,528), `jac` −0.0030 → net −0.003 to +0.001. A wash; rule 13 vetoes.** The value of a
slot is (a) reading whether `n_real` exceeds 3,528 and (b) a board test of h1-test's reach transfer
on an emission-law arm — both modest. Ahmet's call; the package is ready either way. Builder now has
`--emitter` and `--resample-keep-top` (default = shipped, bit-identical), tests pass.

**Production change (built by 179):** `--emitter resample_mix --resample-keep-top 8000`, control
cells = the context's own `context_{A,B,C}.h5ad` controls passed as RAW counts with their OWN
libraries (never rescaled — the exact-total repair peels a rescaled block flat, §0.1), `control_profile`
= the pooled control profile. Exactly two meta differences from a15. **This is the first arm all day
with a transferable mechanism, a verified-pinned profile, and a positive expected value under the
best-supported `n_real`. It is also the one slot that measures `n_real` directly: `hy_kt8000`'s board
`fid` reads ≈ precision if `n_real` ≤ 2,410 and precision × 2,410/`n_real` otherwise.** The SUBMISSION
HOLD is Ahmet's; the decision is presented to him. `experiments/177-resample-emitter/`.

## 2026-09-02 evening — two agents in flight on the last two open questions (177, 178)

**177 — a resampling emitter as the delivery mechanism for the `|H|` ordering.** The anomaly that
motivates it: **Slurm Shadi holds our `pds` (0.842 vs 0.837) and our raw `mse` (0.883 vs 0.872) with
1.48× our reach.** So reach is not bounded by profile quality at our level; something in emission
delivers it. 171 closed spread-schedule delivery only through OUR emitter, whose even allocator and
zero-inflated gamma-Poisson both shift the median. `emit_counts(mode="resample_mix")` (a22's
primitive: real control cells × fold change) is exchangeable with the controls at zero effect, so
`z = 0` there exactly and symmetric noise reduces `|z|` monotonically — the regime in which INBOX-4's
expression-compensating schedule can work. Coordinator's first pass hit a depth confound (real 20k
cells rescaled to the harness's 100k convention → 10,776/10,779 genes called); agent relaunched with
real-depth emission, a quasi-at-real-depth control, and a zero-effect null gate.

**178 — the three riders composed as one production arm.** gamma 0.30–0.45 (+0.0027), sign repair
(+0.002), K562-heavy weights (+0.004 sign) were each priced alone; each costs `pds` and buys the
rest; none clears the +0.0089 profile gate alone. Whether their gains and `pds` costs add,
saturate or interfere has never been measured. Agent builds env-gated switches (default = shipped
behaviour, proven by reproducing a15's audit values), prices all six on h1-test + PLOLO + 096's
claim route, two matched seeds on the composites.

## 2026-09-02 — 176 CMap/LINCS shRNA: NO-GO — the perturbations are too weak to transfer (KO pending)

The agent never ran (Claude Code OAuth expired; **`claude login` needed before any further agents**).
Coordinator ran the gates directly on the 11.8 GB level-5 shRNA file. Coverage: 9 core lines ×
98–118 panel targets (union 127/300), 867 landmark + 7,101 inferred genes on our DE axis. **Sign
accuracy vs the H1 reference pool: CMap lines 0.48–0.55 each, 9-line mean 0.534 landmark / 0.522
inferred, against K562's 0.562 / 0.568.** Added to the trio: +0.005 on landmarks at best, nothing on
inferred, ≤ +0.001 overall on the third of targets covered. Platform control: CMap agrees with K562
at 0.50–0.53 on landmarks and **at chance on inferred genes** (our own two donors agree at 0.540).
**Mechanism: the target's own gene is knocked down only 0.5–0.9 SD in the consensus signatures and
25–34% of targets show no knockdown at all**, against CRISPRi's >80% in our donors. A perturbation
that barely moves its own target cannot carry a downstream response. **CRISPR-KO measured too: the 21-line KO mean adds +0.019 sign accuracy on the 867 landmark genes
(a real, transferable signal — bulk L1000 landmarks of a knockout do carry information about a
single-cell CRISPRi response in hESC) but subtracts −0.011 on the 7,101 inferred genes; net +0.0006
on covered pairs, and coverage (82/300 targets × ~13% of pool pairs) caps the challenge effect at
≈ +0.0007 pooled.** NO-GO for both modalities. **With 166 this closes public cross-line perturbation
data as a route to a better predictor at panel scale.**
`experiments/176-cmap-donors/`.

## 2026-09-02 — WHERE THE DAY LEAVES THE SCIENCE, and the one large source never examined (exp 176 launched)

**The unifying reading of every closure since 06:00.** Each lever found today trades `pds` against
every other member — gamma (097/170), whitening's sign cost (172), sign repair (173), donor
weighting (175) — because **`pds` rewards decorrelation and every other member rewards per-gene
accuracy, and production sits at the `pds`-optimal point of that curve.** The riders that buy
accuracy (gamma 0.30 +0.0027, sign repair +0.002, K562-heavy weights +0.004 sign) each cost `pds`
and none clears the profile gate; composed they might reach +0.005–0.008 and would still cost
`pds`, on a surface (h1-test) that 165 showed is blind to reach-on-profile. **The leader holds
better `pds` AND better `mse` AND better `reach` — a better direction, not a better point on our
curve.**

Denoising cannot supply it: 004 (the transferable signal is full-rank; truncation loses
monotonically) and 156 (deflation does not survive rendering) close both spectral directions. The
recipient's control cells cannot supply it: 160/161 (marginals, R² 1.7%) and 174 (pairwise
co-expression = the top-2 cell-state PCs) close both levels. **The predictor can only improve
through new information, and 166 exhausted genome-scale single-cell CRISPRi.**

**166 never looked at bulk perturbation signatures, and neither has anyone here** — `lincs`, `cmap`
and `l1000` appear nowhere in the repo. LINCS/CMap 2020 holds ~7,000 CRISPR-KO and ~4,000 shRNA
consensus signatures across ~9 cancer lines from tissues we lack, on 978 landmark genes (+ ~11k
inferred), bulk, public. Platform mismatch is a real risk (bulk L1000 at 96 h vs single-cell
CRISPRi at ~7 d) and the gene axis is narrow — but it is nine lines and the 1/k route on a modality
166 did not search. **Exp 176 (launched 14:34) inventories coverage, builds landmark-only donor
tables, and prices each line on h1-test with 175's sign instrument and 173's `pds` pipeline, with a
platform-agreement control against K562.** Stop-at-first-gate design; a clean NO-GO is a result.

## 2026-09-02 — 175 (coordinator): p53-status donor matching REFUTED; the Orion pair adds ~nothing to sign accuracy

Top-tier sign errors are target-level (SBDS 0.41, SOX4 0.46, CYCS 0.47). Tested whether donor TP53
status (K562 null, HEK293T large-T, HCT116 WT; H1 WT) explains them: **no** — on 31 canonical p53
targets the p53-WT donor reads 0.485, TP53-null K562 0.526, the consensus 0.444. What fell out
instead: **per-donor sign accuracy on the H1 pool is K562 0.567, HCT116 0.538, HEK293T 0.531; the
Orion pair alone 0.5425; all three 0.5697.** The Orion pair carries 2/3 of the weight and adds
+0.0027 over K562 alone; its effects are 1.7× larger and less accurate. (2,1,1) weighting buys
+0.004; the per-target oracle router ceiling is 0.603. Third instance of the pds-vs-sign trade:
019's equal weights are right for `pds`, wrong for sign. Rider at most.

## 2026-09-02 — 174 (coordinator): recipient co-expression is the shared programme in disguise — NULL

Tested the one recipient-side signal nobody had: pairwise `r(t, g)` across the 18,400 control
cells as a sign predictor. Raw: **0.5158** on the pool (beats `H` on the weak-donor half, 0.519 vs
0.508). **Remove two principal components and it is gone — 0.4885, below chance.** The signal is
the dominant cell-state axis (proliferation/stress) that target and gene both ride; not `t→g`
regulation. Adding it to `cons` makes `cons` worse (0.543 vs 0.552). With 161 this closes the
recipient's control cells as a source of the line-specific component at both the marginal and the
pairwise level. `experiments/174-coexpression/`.

## 2026-09-02 — 173 (coordinator): sign repair is real, targeted, and nets +0.002 — NO-GO on the profile gate

172's "whitening costs 0.011 sign accuracy" decomposed: **on the 15.7% of the reference's DE genes
where whitening flips the sign relative to the un-whitened consensus, `H` is right 45.6% — worse than
chance — and the raw consensus is right 54.4%.** Flips concentrate where `|cons|` is small (35% of the
bottom quartile, 1.4% of the top decile) and carry 1.6% of `H`'s energy. Setting those signs to
`cons`'s (magnitude kept) and running 169's full realisation pipeline on h1-test with matched seeds:

| member | Δ (overall units) | resolved? |
|---|---:|---|
| `fid` | **+0.00175** (precision +0.0031, call-set sign acc 0.5345 → 0.5454) | yes, twin spread 5e-5 |
| `nmae` | **+0.00116** (−0.0043 raw) | yes |
| `mse` | ≈ +0.0022 (cos +7.8% in delivered space) | estimate via 096 |
| `pds` | ≈ −0.0016 (−0.0042 raw, delivered space) | yes |
| `reach` | ≈ −0.0015, **inside seed noise** | no (matched-seed −0.020/−0.009 vs spread 0.016) |
| **net** | **≈ +0.002** — profile-class gate is +0.0089 | **NO-GO** |

The random-flip null (same count of entries) collapses `pds` by −0.034 and worsens `nmae` by +0.011,
so the targeting is genuine. **But `reach` does not respond to +0.014 pooled sign accuracy because
the gain lives in low-`|H|` entries at the *back* of the p-value order, and reach reads only the
front** (171's result from the other side). The `pds` cost is exactly gamma's trade: sign repair is
surgical partial un-whitening and lands on the same exchange curve (−0.054 vs 097's −0.068 raw `pds`
per +100% cos). It is a mechanism rather than a fitted exponent, and it is a composable +0.002 rider
for any future arm that buys `pds` back — not shippable alone. `experiments/173-sign-repair/`.

## 2026-09-02 — 169 / 171 / 172 ALL CLOSED NO-GO (07:45–07:54). The call-set channel is tapped out downstream of the predictor.

**171 (call set / reach) — NO-GO, three structural results.**
1. **"Our significant calls first" is a mathematical NO-OP.** BH-adjusted p is monotone in p, so
   sorting by `p_adj` already puts every called gene ahead of every uncalled one. Measured: reach
   with the call set replaced by *nothing* or by *everything* is 0.18283 either way. **The call set
   cannot move `reach` at all** — only the p-value order can. (My pole-position framing was wrong in
   mechanism; the ordering diagnosis stands and reproduced to the digit.)
2. **The spread-schedule family is closed by a sign error in its mechanism — mine.** A rank-sum
   test reads a **median shift**, not a variance. Our only noise channel is a zero-inflated
   gamma-Poisson (φ = 4) that moves the median as a side effect, so `|z|` is **U-shaped in `s` with
   its minimum at the shipped 0.25**: `n_pred` runs 3,581 → 6,034 → 9,472 as flat spread goes
   0.25 → tier → 1.0. "Suppress unearned significance by adding noise" is not available; past 0.25,
   noise *manufactures* significance. My expression-compensating schedule was built, was exactly
   `mse`-free as predicted, and **doubled** the expression dependence (+0.485 → +0.939); its
   permutation control scored higher than it did. **No setting of the per-gene spread vector delivers
   the `|H|` ordering.**
3. **a15 sits +0.01313 above the `mse` knee** (`s²B` = 0.23268 vs claim 0.24581), so `mse` is exactly
   flat for any arm above it. The one `mse`-free reach arm, `--keep-top 8300`, reads +0.0252 raw
   twin-averaged but 69% comes from five perturbations, it falls to +0.0094 on `n_real ≥ 20`, and
   the board's own `fid`-vs-`n_pred` series charges −0.0011: **net +0.0003 against a +0.0019 gate.**

**172 (learned sign-confidence key) — NO-GO.** Four model classes (22–2,900 parameters) trained
LODO on 783,000 labelled pairs from the three production donors read **0.2212–0.2429** raw reach on
h1-test, straddling `|H|`'s 0.23533 — every one inside `|H|`'s noise in both directions. The learned
keys beat the un-whitened consensus they were built from (+0.009 to +0.028) but **only recover what
spectral whitening already provides for free.** Pre-registered falsifier (top-5% accuracy > 0.75)
came out negative: 0.651 → 0.674. **The ordering gap is not learnable from AgingLab's feature
family.** Flagged for the record: **94.4% of h1-test's reach lives in 32 of 85 targets**; and
**whitening costs 0.011 pooled sign accuracy on the DE pool while buying `pds`** — the same
trade 097 measured as `kappa`.

**169 (sharpening, six members) — NO-GO, and `nmae` is why.** The `pds` lever is real and large —
**+0.019452 scaled on h1-test, positive on all four folds, 2.9× anything in this project's
history** — and the blanket reach toll rule 13 charges is **falsified** (five matched seeds:
+0.00061 overall, whole CI above the −0.00243 toll). `mse` is clean. It dies on `nmae`:
**−0.00023 on HCT116 vs −0.01289 on h1-test, a 56× disagreement**, and h1-test is the surface the
board validated on that member. Nearest precedent is a13 — the other per-gene magnitude
reweighting — gated positive on both harnesses and lost −0.0203 on the board with raw `nmae`
1.0185. Same family, same member, same shape. The optimum exponent is an exact identity (threshold
1/2), not a fit. 161's five-line `sign(C)·|C|²` result is the same lever measured without `nmae`.

**Where this leaves the project, stated plainly.** Eight experiments overnight, all NO-GO, no slot
spent, hold intact. The diagnosis is now sharp: our profile is world-class (`pds` 5th, `mse` 23rd of
581), our DE members are poor (176/169/192), and the reason is that **our design sign is right only
53.7% of the time on the reference's DE genes** — everything downstream of that (call set, ordering
key, spread schedule, learned confidence, magnitude sharpening) has now been measured and cannot
recover it, because the emitter cannot deliver an ordering and the ordering is not learnable
anyway. **The one remaining problem is the predictor's per-gene sign accuracy**, exactly as the
five closed direction families said. The open sub-threshold levers that would compose with a better
direction: gamma 0.30 (+0.0027, costs only `pds`), keep_top 8300 (+0.0003, `mse`-free).

## 2026-09-02 — THE ORDERING GAP IS WORTH 7 POINTS, AND IT SPLITS INTO TWO SOLVABLE HALVES (exp 172 launched)

Two more numbers on the same h1-test surface (a15 DE table, `H` from the 095 artifact, same
`Surface.reach(key=…)` instrument, `|H|` known-answer 0.23533 reproduced):

| ordering of the reference's DE pool | reach |
|---|---:|
| production (our emitted p-values) | 0.18283 |
| `|H|` — best simple key; **no composite** (× or ÷ expression, centred, z-scored, squared) beats it | 0.23533 |
| **ORACLE — our own signs, perfectly ordered** | **0.64553** |
| board leader `cqawesome` | 0.26063 raw |

**Our signs are already good enough to lead the field by 2.5×. The deficit is entirely ordering.**
The oracle-to-production gap is **0.46 raw reach ≈ +0.51 scaled ≈ +0.085 overall** — the largest
single quantity of locked score anywhere in this project. It splits cleanly:

1. **Delivery** (0.183 → 0.235): our emitted p-values do not track `|H|` because Wilcoxon power
   is expression-set (`σ² = s·m + φs²m²`, previous section). **Exp 171** owns this.
2. **Key quality** (0.235 → 0.646): `|H|`'s top-1% sign accuracy over the pool is 0.779, top-5%
   0.724 — never the 0.9 the metric needs, so `k*` stops after a few genes. **Exp 172** (launched
   07:55) learns a per-(target, gene) sign-confidence key from the donor lines by LODO across the
   three production donors, no gene or target identity, evaluated only as `reach(key=c)` and top-X%
   accuracy on h1-test. This is what `AgingLab` publicly describes ("578-parameter head … four
   measurable scalars … no gene identity and no target identity"), and they hold the field's best
   `nmae`.

**The learnable oracle is 0.621.** 94.5% of emitted signs equal the design sign `sign(H)`, so the
lucky-coin-flip component is small and the H-oracle (design sign correct, perfectly ordered) reads
**0.62119** — 96% of the headroom is a property of `H` and therefore learnable in principle. The
signal is thin: design-sign accuracy over the pool is **0.537**, rising monotonically by `|H|` decile
from 0.510 to only **0.624** at the top; the bottom half of the pool is at chance. Emission sign-flips
are confined to `|H|` < 0.013. And the leader's 0.26 is plausibly a modestly better direction
(0.56–0.60 accuracy) delivered properly, with no confidence model at all — so 171 is the surer half.

**A calibration that changes how 172 must be built:** `|H|` has pool `Spearman(|H|, sign-correct)`
of only **+0.069** yet reaches 0.235, while `sign-correct + N(0,σ)` needs Spearman **0.316** to
reach 0.224. Reach reads only the front, so **precision at the top is the objective, not rank
correlation** — a key mediocre on average but excellent in its top 5% wins. The two experiments
compose: 172's `c(t,g)` becomes the `rank_score` that 171's expression-compensating schedule
delivers through the emitter.

## 2026-09-02 — THE FALSIFIER FIRED: our p-value order is set by EXPRESSION, and the fix is not the one we kept trying

Measured on artifacts already on disk (a15's h1-test DE table from 164; the design effect `H`
recovered from the same 095 artifact by 169; axes aligned by name 10,779/10,779 and 85/85; the 095
artifact used `k562_gwps` + `orion_hek293t` only, so no `h1_2025` leak, and `H` comes from the
artifact's own counts so no differential leakage is possible).

**Our p-value order anti-tracks our own effect.** `Spearman(−log10 padj, ·)` median over perts:
vs `|predicted effect|` **−0.121** (universe) / **−0.391** (inside our call set); vs expression
**+0.444** / **+0.485**. Inside our significant block, the genes we rank most confident have the
*smallest* predicted effects. Sign accuracy at the front of our reach walk is **0.521** (chance);
ordering by `|effect|` would put **0.621** there with a clean monotone gradient to 0.499.

**`reach` on h1-test, same calls and lfc, only the order changed:**

| ordering | reach | Δ raw |
|---|---:|---:|
| production (sig first, then our padj) | 0.18283 | known answer ✓ |
| `|pred_lfc|` (157's key, the *emitted* lfc; Spearman 0.87 with `|H|`) | 0.22173 | +0.039 |
| **`|H|`** | **0.23533** | **+0.0525 (+29%)** |
| **sig-first, then `|H|` inside** = **`|H|` over called genes only** | **0.21914** | **+0.036** |

The last two are identical to five decimals: **the call set already contains the reachable genes;
the order inside it is what is broken.** +0.0525 raw would move board reach 0.1742 → ≈0.2267,
almost exactly the field line's 0.2242 — **this one defect is our whole deficit below the field
line.** ≈ **+0.0097 overall** at the full ceiling, ≈ **+0.0067** inside-block.

**Why every realisation attempt failed (093, 095, the a15 tier).** The metric reads our Wilcoxon
p-values against the fixed control cells: `z_g ≈ μ_g h_g / sqrt(σ²_ctrl + σ²_ours)`. At the shipped
constant spread, `z/h` runs ≈ 14 · 38 · 72 across `μ` = 1 · 10 · 100 — **expression moves `z`
five-fold**, which is the +0.485. `keep_top` tiers spread by `|effect|` rank, so the `μ` term runs
straight through it; that is the mechanism behind "counterfactual key gains do not survive the
emitter." **The fix is an expression-compensating spread schedule** — `s_g` increasing in `μ_g` so
that `μ_g / sqrt(σ²_ctrl + s_g²μ_g²)` is flat — which suppresses high-expression genes' unearned
significance and lets `h_g` set the order. We can only add noise on our side, so the low-`μ` Poisson
floor bounds how much of the +0.036 survives. **Measured: the exact emitted variance is
`σ² = s·m + φ·s²·m²` (`_quasi_block`), the crossover to the flat regime is `m* ≈ 1/(φs²) ≈ 200 CPM`
at the shipped setting, the matchable floor is ≈ 12–16 CPM, and ordering by `|H|` only above 12 CPM
recovers 73.4% of the inside-block ceiling (+0.0267 raw ≈ +0.0049 overall as a proxy lower bound)
while only 9.3% of the reference's DE pool sits below that floor.** The schedule is `s` INCREASING
in `m` (pull the plateau down; the floor cannot be raised), it requires `φ > 0` (pure Poisson has no
plateau), and it only raises `s`, so with a15 at 99.7% of the dispersion cap `mse` should be ≈
neutral — to be verified, not assumed. **It must be priced on
`mse` first** through 096's `claim` route: a15 sits at 99.7% of the dispersion cap and the tail block
exists precisely to feed the jackknife credit (157 §6). `pds`/`nmae` should be bit-identical and must
be verified, not cited.

This supersedes my "spread contrast at fixed call count" framing two sections down: contrast keyed
on `|effect|` is what the project has already tried three times. The key has to be `μ`.

## 2026-09-02 — THE CALL-SET MECHANISM, QUANTIFIED — and the decoupling it implies

Solving 169's measured `n_pred = 6034` / `sign = 0.53435` as a two-population mix:

| informative genes | implied accuracy |
|---:|---:|
| 600 | 0.8454 |
| **900** | **0.7303** |
| 1200 | 0.6727 |

**168b measured the top 5% of `|D|` at 69.8–73.1% of the alignment; a22's board probe read the top
decile at 0.7359. ~900 informative genes implies 0.730. Three independent measurements agree.**
900 is ~5% of the gated universe — exactly 168b's stratum.

**So our signs are not bad. ~85% of our call set carries no signal at all, and every call-set metric
averages over the whole thing.**

**Why the obvious fix fails.** `fid = precision × min(1, n_pred/n_real)`. At the inferred challenge
`n_real ≈ 2670` (093's number, flagged circular by rule 9):

| n_pred | precision | yield | fid |
|---:|---:|---:|---:|
| 900 | 0.7300 | 0.337 | 0.2461 |
| **3000** | **0.5690** | **1.000** | **0.5690** |
| 6034 (ships) | 0.5343 | 1.000 | 0.5343 |

Shrinking to the pure head collapses `fid` through the yield cap — which is why `keep_top` ships at
8000 although the builder's own docstring recommends **850**, and which re-reads 093's board
inversion (local +0.0202 → board −0.00291) as a yield-cap effect the depth-confounded harness could
not see. Note the `n_pred = 3000` row predicts `fid` 0.5690 against our 0.5343 at full yield: a mild
interior optimum may sit above the yield knee and below 6034.

**THE LEVER. `reach` has no yield term, and the order *inside* our significant block is set by our
own p-values — which the per-gene emission spread sets.** So `keep_spread`/`tail_spread` decide not
only who is called but **who leads the reach walk**. `fid` needs `n_pred ≳ n_real`; `reach` needs
only a pure *front*. **One parameter is serving two objectives that want opposite things.**

The arm is therefore **spread CONTRAST at fixed call count**, not a smaller call count: hold
`n_pred` near 6000 for the yield cap while driving the ~900 informative genes to far smaller
p-values than the other ~5,100. The builder's `--spread-map graded` with `--spread-lo/--spread-hi`
is the natural vehicle, since it makes p-value rank track `|effect|` rank continuously rather than
in two blocks. **This is a mechanism, not a hyperparameter, so it should carry to D/E/F.**

**The falsifier, which 171 must run first.** The chain rests on our p-value order tracking our
`|effect|` order. A Wilcoxon on counts gives highly expressed genes more power regardless of effect
size, so if expression dominates, the informative genes are not at the front and no contrast setting
will put them there. **171 must measure the rank correlation of our adjusted p-value against
`|predicted effect|` and against `b_ctrl` expression over the call set.** If expression wins, that
number is the headline and it redirects the repair.

## 2026-09-02 — 161 CLOSED NO-GO: the cross-line residual is not predictable, and the fifth wall

Experiment 161 attacked the ~70% of effect energy that is line-specific — the one thing 159 and 160
left open — by learning the held-line residual `y = R_L − C` from recipient control features and
`consensus × context` interactions, across five lines (K562/RPE1/Jurkat/HCT116/HEK293T; **1,883
shared targets, zero on the challenge panel**, so the trainable pool is disjoint from the deployment
set by construction). Eight known-answer rows passed on every fold, including the oracle-residual
ceiling at `pds` 0.99996.

**NO-GO, with a three-legged mechanism:**

1. **Context features carry nothing.** A context-only model is annihilated by the mandatory
   across-target centring: eval `R²` **−0.0024**, cosine to the true residual **−0.0033**, delivered
   Δpds statistically identical to its own permutation null on all five lines.
2. **Context slightly *loses* to target-side.** Cosine to the true line-specific residual: **T (no
   recipient information at all) 0.0606**, X (consensus×context) 0.0592, ALL 0.0550. The two
   scalars a donor table cannot express — per-gene control-cell detection rate and dispersion, built
   from 419k control cells — make it **worse**: dropping them moves worst-line Δpds from −0.00097 to
   +0.00382.
3. **The gain that survives is not prediction.** A zero-learning per-entry re-shaping
   `C ← sign(C)·|C|²` is worth **3.1×** the entire 27-feature learned model, and once it is in the
   baseline the learned model adds +0.0024 on `core4` but **−0.0027** on `core5` while costing
   cosine on every line of both surfaces.

**Do not reopen cross-line residual prediction without a mechanism that beats a per-entry magnitude
re-shaping of the consensus** — that is what a flexible learner finds when pointed at this
objective. This is the **fifth** independent route into the 159/160 wall and the first to attack the
residual rather than the consensus.

**One live thread it leaves:** leg 3 is global sharpening at q = 2, which 168d closed as a TERMINAL
NO-GO on h1-test. 161's five lines contain neither h1-test nor any low-blob surface, so the most
likely reading is that it **re-measured 168d's Orion overstatement on three more lines** — which
would strengthen the instrument rule from two surfaces to five. Handed to 169 (which owns the
family) as a one-row adjudication via the target-independent-energy discriminator; it must not delay
169's six-member table.

## 2026-09-02 — THE BOARD SAYS WE ARE A `pds` SPECIALIST WITH A `reach` HOLE (581 teams)

Pulled the live board at 06:30 UTC+3. **The field has grown from 129 to 581 scored teams and the
leader is now `cqawesome` at 0.2169** (was Vivai 0.1341 in §8's snapshot). We are **25th at
0.1734**. Snapshot: `data/leaderboard/2026-09-02T0630Z.json`.

Decomposed against the leader — gap **+0.0435** — with our rank in the full field:

| member | us | leader | their edge | share of gap | **our rank / 581** |
|---|---|---|---|---|---|
| mse | 0.1216 | 0.2375 | +0.01931 | 44.4% | 23 |
| **reach** | **0.1074** | **0.2049** | **+0.01625** | **37.3%** | **192** |
| pds | 0.7476 | 0.8165 | +0.01148 | 26.4% | **5** |
| nmae | 0.0941 | 0.1207 | +0.00444 | 10.2% | 68 |
| jac | −0.0005 | −0.0096 | −0.00151 | −3.5% | 169 |
| fid | −0.0301 | −0.0688 | −0.00646 | −14.8% | 176 |

**We are 5th in the world on `pds` and 192nd on `reach`.** 191 teams beat us there. Field best is
The_Way at 0.2379 scaled / 0.2898 raw against our 0.1074 / 0.1742. Matching the leader's reach
alone is **+0.0162 overall**; field best is **+0.0217**.

**The honest caveat:** the 14 teams above reach 0.20 average `pds` 0.6487 and overall 0.0929 — most
are buying reach by wrecking discrimination. **But `cqawesome` holds reach 0.2049 together with the
board's best `pds` 0.8165 and best `mse` 0.2375, so the trade is not forced, and the leader is not
winning by calibration — they beat us on `pds`, `mse`, `reach` and `nmae` at once, which is what a
genuinely better direction estimate looks like.** We beat them only on `fid` and `jac`.

**THE UNIFYING PATTERN, and it reframes the whole project.** Split the six members by what each
one *reads*:

| what it reads | member | our rank / 581 | us | field best | headroom |
|---|---|---:|---:|---:|---:|
| the pseudobulk **profile** | `pds` | **5** | +0.7476 | +0.8165 | +0.0115 |
| the pseudobulk **profile** | `mse` | **23** | +0.1216 | +0.2375 | +0.0193 |
| LFC **magnitudes** over the reference's DE genes | `nmae` | **68** | +0.0941 | +0.1838 | +0.0149 |
| **our own DE call set** | `fid` | **176** | −0.0301 | +0.0227 | +0.0088 |
| **our own DE call set** | `jac` | **169** | −0.0005 | +0.0196 | +0.0034 |
| **our own DE call set** | `reach` | **192** | +0.1074 | +0.2379 | +0.0218 |

**The three members that read our call set rank 176 / 169 / 192. The three that do not rank
5 / 23 / 68.** The two groups do not overlap, and the boundary is exactly whether the metric passes
through our emitted cells' Wilcoxon calls. **Our profile is world-class; our DE call set is not.**

Combined field-best headroom on the three call-set members is **+0.0339 overall**, which alone
moves 0.1734 → **0.2073** against the leader's 0.2169 — **with no improvement to the underlying
direction at all.** That is a larger, cheaper and better-evidenced target than anything in the five
closed direction families, and it is a *mechanism* rather than a hyperparameter, so it should carry
to the unseen D/E/F panel.

This also retro-explains the project's recent history: a22 lost `-0.01066` on `reach`, a23's `reach`
carried 182% of its loss, and 093's board result inverted — **every one of those was the call-set
path misbehaving**, and we read each as a separate accident rather than as one broken subsystem.

**AND THE DEFICIT IS A DEFECT, NOT A TRADE — this is the decisive evidence.** Across all 581
teams `raw_reach` and `raw_pds` are coupled at **corr +0.701** (181 corrected; +0.831 was wrong), with field regression
`raw_reach = −0.0751 + 0.3577 · raw_pds` (residual sd 0.0409). At our `raw pds` 0.8367 the field
predicts `raw reach` **0.2242**; we deliver **0.1742**, a residual of **−1.22σ**. Only 39 of 581
teams sit further below the line — **we are in the bottom 6.9% conditional on our own
discrimination**, and closing merely to the field line is **+0.0092 overall**.

Every team on the board with `raw pds` ≥ ours:

| team | raw pds | **raw reach** | raw mse | overall |
|---|---:|---:|---:|---:|
| pop | 0.8519 | 0.2781 | 2.661 | 0.1663 |
| cqawesome | 0.8682 | 0.2606 | **0.762** | 0.2169 |
| Slurm Shadi | 0.8418 | 0.2584 | **0.883** | 0.1918 |
| Garrrbage | 0.8454 | 0.2492 | 4.425 | 0.1612 |
| **us** | 0.8367 | **0.1742** | 0.872 | 0.1734 |

**All four reach ≥ 0.2492; we are the only one below 0.24.** The amplitude explanation dies here:
**Slurm Shadi holds `raw mse` 0.883 — statistically our own 0.872 — with 1.48× our reach**, and
**cqawesome is *better* shrunk at 0.762 with 1.50× our reach**. Our shrinkage is not what costs us
reach; four comparable pipelines do something we do not.

**Experiment 171 is now the top priority.** The reach spec makes the mechanism concrete: the
reference's DE genes are ordered by **our significant calls first**, then our adjusted p, then p,
then `|lfc|`, then **gene name** — and `P(k)` is cumulative, so poison at the FRONT is permanent.
157 measured that **74–79% of our calls carry zero predicted effect**; if those take the pole
positions they drag `P(k)` to chance immediately and `k*` collapses. That would explain rank 192
on `reach` alongside rank 5 on `pds`, and it locates the fault in our **call policy** rather than
in our ordering — a different conclusion from 157's "the ordering channel is closed."

## 2026-09-02 — gamma as an `mse` lever: RE-OPENED, CORRECTED, then TERMINALLY CLOSED by exp 170 (see the 170 entry below — it supersedes this section's 'sub-threshold live' reading)

I re-opened `WHITEN_GAMMA` on the argument that a19's board point makes gamma an 81×-larger `mse`
lever than a `pds` lever, that 097 measured lowering gamma raises `kappa`, and that **097 never
priced `mse`**. Experiment 170 was launched on it. **The third premise is false and I killed 170
seven minutes in.** 097 §9h contains a full six-member table with an `mse` column, measured
directly on both surfaces.

What 097 §9h actually measured on h1-test, against shipped gamma 0.65:

| gamma | mse | pds | net overall |
|---:|---:|---:|---:|
| 0.00 | +0.0108 | −0.0142 | −0.0033 |
| 0.15 | +0.0083 | −0.0069 | +0.0022 |
| **0.30** | +0.0055 | −0.0033 | **+0.0027** |
| 0.45 | +0.0029 | −0.0013 | +0.0023 |
| 0.85 | −0.0025 | +0.0006 | −0.0021 |

**Two things survive, and one of them is a genuine new measurement:**

1. **a19's implied slope overstates the truth by 6.4×, and that factor IS the `lam` confound.**
   a19 gives `d(raw mse)/d(gamma) ≈ +0.6075`; 097's direct h1-test measurement gives **+0.0947**.
   a19 moved gamma *and* `lam` 0 → 0.007, and this now quantifies how the −0.0959 `mse` loss split:
   **the ridge did roughly six-sevenths of the damage, not the exponent.** Any future operator arm
   should treat `lam` as the dangerous parameter and gamma as the mild one.
2. **The gamma lever is real, h1-test-positive at 0.15–0.45, and was closed ONLY by Orion** — the
   surface 168d has since discredited for sizing. Its peak is **+0.0027 at gamma 0.30**, which is
   above rule 13's shape-only gate (+0.0019) but **far below the +0.0089 gate that an operator- or
   amplitude-class arm must clear**, and gamma is squarely amplitude-class because it reshapes the
   profile the renderer reads. **So it stays unshipped** — but it is now a *sub-threshold live*
   lever rather than a closed one, and it is the obvious thing to compose with any future arm that
   buys back `pds`, since its cost is almost entirely `pds`.

**The methodological point, recorded because it nearly cost a slot:** I flagged the +0.016
extrapolation as "a hypothesis, NOT a result" and then checked the source document before letting
the agent build on it. The check took four minutes and refuted the premise. Rule 11's "record a
path and a sha, never a description" has a sibling: **never re-open a closure without reading the
closing document's own table.**

## 2026-09-02 — 161 (interim): the line-specific residual is NOT predictable from basal features

Experiment 161 fit the held-line residual `y = R_L − C` (what a line does that the equal-weight
donor consensus `C` does not) from recipient control-profile features and `cons × context`
interactions, on a 5-line surface (K562/RPE1/Jurkat/HCT116/HEK293T; 1,883 shared targets,
**0 on the panel**) with all three known-answer rows passing (oracle truth `pds` 1.000; oracle
residual `pds` 0.9998 / cos 0.781 — the ceiling of the whole family).

**`R²(y) = +0.015–0.017`.** The model explains under 2% of the residual. The orthogonal
component it predicts has cosine 0.13–0.14 with the true orthogonal residual. Delivered at
half weight it reads `pds` **+0.005** with **cosine −0.006**, and at full weight `pds` +0.005 /
cos −0.02; every heavier weight is catastrophic. A `pds` gain with a cosine *loss* is the
whitening pattern — hubness reduction, not prediction. Combined with 159 (production sits at
~30% of the context-independent ceiling) and 166 (donors exhausted), this says the 70% of
effect energy that is line-specific is, on present features, **unpredictable**. The agent died
before writing `RESULTS.md`; the watchdog is resuming it to finish the write-up from the
result files on disk (`experiments/161-crossline-residual/work/*.json`).

## 2026-09-02 — 164: the call-set artefact is real, at chance, and NOT worth removing

**NO-GO on the whole family, with the mechanism.** My own thesis (that ~30% of the call set is
signal-free artefact diluting `fid`'s precision, worth up to +0.049) is **refuted**, and the way
it failed is instructive.

The artefact is exactly as described. On h1-test, split by the shipped spread vector:

| block | calls | directional precision |
| --- | ---: | ---: |
| keep (spread 0.25) | 297,869 | 0.5122 |
| **tail (spread 1.0)** | 214,265 | **0.4999** |
| all | 512,134 | 0.5070 |

The tail block is 42% of the call set at **exactly chance** (per-perturbation sd of the mean
0.0011; paired keep − tail **+0.0122**, t = 4.33). Confirmed at scale on the **real challenge
controls**: `spread 1.0 / phi 4` is significant at zero predicted effect for **82–86%** of its
genes. **But 157's "100% DOWN" was the Wilcoxon z sign, NOT the scored direction** — the scored
`lfc` is up **42–44%**, near-balanced, which is exactly why the block reads precision 0.4999.
The shipped mixture's zero-effect call rate is **26.8 / 25.5 / 23.5%**, not 30.6/31.4/28.8%.

**But removing it is worth +0.00066 on `fid`, not +0.049.** The error was mine and it is
AGENTS.md rule 8: pooled precision rises 0.5070 → 0.5122, but the metric averages **per
perturbation**, where the gain collapses. `jac` then pays more than `fid` collects — overlap with
the reference falls 478 → 254, jac 0.0682 → 0.0538.

**And `mse` closes it decisively — harder than the interim numbers suggested.** Final
measurement: `claim` = **0.24581 against the cap 0.24653**, so **a15 is at 99.7% of the
dispersion ceiling** with `ρ = 0.9466`. Every arm drops `claim` below `s²B`, pushing raw `mse`
past 0.98786 where **`score_mse` clamps to zero** — so each arm forfeits the ENTIRE member,
not a slice of it:

| arm | Δpds | **Δmse** | Δfid | Δreach | Δjac | **NET overall** |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| armA tail φ=0 | 0.000000 | **−0.018520** | +0.001145 | +0.002394 | −0.000333 | **−0.015314** |
| armB tail spread 0.5 | 0.000000 | **−0.018520** | +0.000359 | +0.000301 | +0.000760 | **−0.017099** |
| armC a8 flat | 0.000000 | **−0.018520** | +0.000852 | +0.001488 | +0.000570 | **−0.015610** |

Known-answer row a8→a15: modelled a8 raw `mse` 1.3711 against board-observed **1.3627 (+0.61%)**,
accounting for 87% of the observed +0.019726. Re-tuning α per arm does not rescue any of them
(best case armA at −0.0008).

**And the board has already refuted this family five times.** `n_pred` 1,563 / 1,911 / 3,836 /
4,424 / 5,961 → raw `fid` 0.49126 / 0.49654 / 0.5021 / 0.5022 / 0.50178. **Shrinking the call set
has never once raised board `fid`.** Deleting the entire tail block perfectly and for free is
worth **+0.00089 overall**, and our gap to the board's best `fid` is **9.7x that**.

The tail block does two jobs priced on opposite sides of the ledger: **worth +0.0009 to remove as
a call set, worth +0.0185 to keep as dispersion.** a15 ships unchanged, zero flags different.

Two sub-hypotheses also dead: the 794/958/997-gene sub-gate band (DE-visible but below
production's CPM gate) is **never called significant** (0.0–0.2% at any policy), and
`--min-cpm 5` is `mse`-neutral but exp 088 already measured 10 as a sharp optimum.

**The unifying conclusion: `fid` is model-limited too.** Removing everything at chance takes
precision to 0.512, not 0.60. `fid` now joins `reach`, `nmae` and the per-target axis — five
independent routes (157, 158, 159, 160, 164) all ending at "we need a better predictor of
per-gene direction", not a better lever.

### Corrections to the composed gate, from 162's final run

The gate above was written off 162's interim log; its final run is stricter and adds two honest
failures that must travel with it.

- **A shape-only arm would need predicted Δ`pds` > +0.0126 raw to clear its own noise** — about
  2x the largest value in this project's history. **Under this gate none of the six historical
  geometry arms earns a slot.** So the gate's practical use is as a **veto**, not as
  authorisation. Read it that way.
- **`mse` fails PLOLO structurally, not by bad luck.** PLOLO's cross-line truth wants an
  amplitude 17-46x smaller than production (`gstar` 0.022-0.058), so its α-gradient always points
  down and every amplitude arm is mis-priced. A scale-sensitive `nmse` endpoint *did* recover
  a19's sign — the `sse` endpoint read a19's −0.0959 `mse` collapse as **+0.0000148** — but it
  still loses to the predict-zero null. Use 096's `claim` route for `mse`, never PLOLO.
- **Two failures in the instrument itself:** the `pds` result **depends on `k562_gwps`** being in
  the fold mean (that one three-fold subset fails the null), and 156's known-answer row
  (`intended` > `delivered`) reproduces on only **3 of 4 folds** — `k562_gwps` inverts.
- **`|effect|` is not a usable stand-in for the reach ordering.** The spec orders reach's pool by
  predicted adjusted p-value with `|lfc|` only fourth; substituting `|effect|` reads 2/6 signs at
  every pool size and prices a23's reach at **+0.000057** where the board read **−0.007933**.
  Any future reach claim must come from an emitted p-order, not a magnitude proxy.

### 162's amplitude anomaly — CLOSED the same night, it is a branch error (not a lever)

162's final message flagged that a11/a12/a14 "nominally differ only in α" yet their board `mse`
misses the quadratic-in-amplitude law by **8x the size of a real move**, and suggested either the
tier-8000 spread policy is not `mse`-neutral or per-context spans are wrong. **Neither. The law
was applied across a kink.**

Three things settle it, and the first two are already in this repo:

1. **a11→a12 is not a pure-α pair.** 162's own `board.py` labels it `channel="mixed"` (α
   0.838→1.0108 **and** the tier-8000 spread policy). a10→a11 is unreadable because a10's
   `score_mse` is **censored at the clamp floor** (0.0), so its raw value cannot be recovered. The
   only clean, uncensored pure-α pair in our whole history is **a12→a14**.
2. **`raw_mse(s)` is not a quadratic — it is a sum of three per-context piecewise functions with
   kinks at α = 0.982 / 1.008 / 1.042.** That is a14's own founding result: the scorer applies
   `min(claim, budget)` **inside each context**, `min` is concave, and the panel-mean evaluation
   is a different function. Below a context's cliff the metric refunds the differential part of
   the amplitude, so `raw_mse` *decreases* in `s`; above it, it increases. The minimum is the
   cliff itself, not an interior argmin.
3. **So a12 (α 1.0108) → a14 (α 0.982) spans two of the three kinks**, and a two-point slope
   across a kinked sum is not a quadratic slope.

I reproduced the failure to be sure it is the explanation: fitting the smooth `ρ = 1` branch to
that pair implies an aggregate cosine of **0.13389** against the independently board-fitted
**0.0722** (**1.85x**) and an interior argmin of **α\* = 0.21** against a shipped 0.982 (off by
79%). Both are nonsense, and they are nonsense in the same direction and for the same reason as
096's own retracted inversion — "inverting a formula whose omitted terms dominate produced a `c`
that was off by ~4.5x". The clean pair's observed direction is still informative and consistent:
lower α gave **lower** raw `mse` (0.901243 → 0.882761), i.e. both points sit on the *increasing*
side, which is exactly what "a14 is the true per-context cliff" predicts.

**Do not chase this, and do not re-derive an amplitude optimum from board slopes.** The only
valid instruments for `mse` are 096's `claim` route (per-context, cliff-aware) and the board
itself. Cost of the closure: one 20-line calculation.

### One inherited assertion that is now known false

`pred_lfc` invariance under emission changes **fails as stated**: 1.1% of predicted directions
flip between arms, and a re-seeded twin of the *same* policy flips 1.55%. So the honest statement
is "invariant to within re-seed noise", not "invariant". 095 saw this as `nmae` agreeing only to
4 decimals; it is the same fact.

---

## 2026-09-02 — WHY DONOR EXPANSION IS DATA-LIMITED, NOT EFFORT-LIMITED (structural, verified)

Experiment 159's `1/k` law says the consensus approaches its ceiling in the number of
**independent** cell lines, that we are worth **3.1** effective donors and the leader behaves
like **10.3**. That makes "acquire more lines" look like the top long-view lever. It is not
available, and the reason is structural rather than a matter of searching harder.

**Panel coverage of every donor library we hold or have just built:**

| library | perturbations | panel targets | verdict |
| --- | ---: | ---: | --- |
| `orion_hek293t` | 17,977 | **294 / 300** | STRONG — ships |
| `orion_hct116` | 17,767 | **290 / 300** | STRONG — ships |
| `k562_gwps` | 9,787 | **269 / 300** | STRONG — ships |
| `h1_2025` | 255 | 20 | training-only |
| `frangieh_melanoma` (built tonight, new tissue) | 235 | **5** | training-only |
| `nadig_hepg2` (built tonight) | 2,160 | **0** | ZERO |
| `nadig_jurkat` | 2,270 | **0** | ZERO |
| `rpe1` | 2,204 | **0** | ZERO |
| `k562_essential` | 2,021 | **0** | ZERO |

**The mechanism, verified rather than asserted.** The four zero-overlap libraries are 90-97% the
**same gene set** as each other (pairwise overlaps 1,828-2,123; union 2,389 genes) — they are all
the DepMap common-essential set. And the union of all four intersects the 300-gene challenge
panel in **exactly 0 genes**. The panel was constructed to exclude common-essential genes, so
**every essential-gene library is disjoint from it by construction**, no matter how many cell
lines it covers.

**Therefore only genome-scale screens can cover the panel** — note the perturbation counts: the
three panel-covering libraries carry 9,787-17,977 perturbations, the useless ones ~2,200. And
**exactly three genome-scale CRISPRi Perturb-seq datasets exist publicly, and we already ship all
three.** Replogle's RPE1 arm is essential-only; Nadig's HepG2 and Jurkat are essential-only;
Frangieh is a 235-gene immune-evasion library; Jiang is 218 targets (9/300); every other source
surveyed across experiments 105, 107, 116, 119-121, 126, 129, 131 is small-library or
access-limited.

### What this means, and it is worth being blunt about

1. **The effective-donor count cannot be raised for INFERENCE** unless a new genome-scale
   multi-line Perturb-seq atlas is published before Nov 5. Stop searching for donors on the
   strength of the `1/k` law; the search space is closed and the reason is the panel's own design.
2. **All of these libraries remain first-class TRAINING assets.** Experiment 161's residual model
   needs only targets shared *between pairs of lines*, and its `core5` surface already has 1,883
   shared targets with **zero** on the panel. Adding HepG2, melanoma and Jiang's six lines
   increases line diversity for a context-conditioned model, which is exactly what a universal
   residual head needs to generalise. **Coverage is the wrong axis for training and the binding
   axis for inference** — that distinction is what the old "donor expansion closed" note missed.
3. **The leader's implied 10.3 effective donors therefore cannot come from more lines.** Either
   159's effective-donor figure is overstated (its harness-blob headline was overstated 2-3x, so
   treat this as live), or the leader extracts far more from the same three sources, or they hold
   non-public data. The first two are actionable; the third is not. **Do not conclude they have
   more donors — conclude they have a better estimator on the same donors.**

---

## 2026-09-02 — 166: the panel-covering universe is EXHAUSTED. The 1/k route is closed on data.

The earlier structural argument (essential libraries are disjoint from the panel by
construction, so only genome-scale screens can cover it) is now backed by an exhaustive search:
Hugging Face, **GEO 258 series** across two queries, **Europe PMC 229 papers**,
BioStudies/ArrayExpress, the Arc Virtual Cell Atlas, the **IGVF portal (6,902 measurement sets,
374 construct libraries)**, and all 64 h5ad files already on disk (118 GB of Orion included).

**Every genome-scale human Perturb-seq screen that exists is already in `data/donors/` or already
logged as rejected.** Three closures worth recording because they retire open rows:

1. **X-Atlas/Pisces is an EMPTY repository** — 3 files, `used_storage: 0`, unchanged since
   2026-03-25. `AGENTS.md` §12 and `docs/CURRENT_STATE.md` both carried it as "unavailable";
   it is not merely unavailable, there is nothing in it. And even when it lands, its HepG2/iPSC
   arms are a 200-gene validation set.
2. **IGVF is closed by one number:** its largest gene-targeting guide library, across 374
   libraries, is **24 genes**.
3. **McFaline-Figueroa is closed** (494-gene kinome, 23/300) — the last open row of
   `AGENTS.md` §5.3.
4. **Zenodo** (166's one 403/WAF gap) is reachable with a browser user-agent; I ran three
   queries and its full-text search is keyword-OR'd noise (lightning K-leaders, dairy cattle,
   oil palm). Nothing relevant, and any Zenodo-hosted atlas would have surfaced in the Europe
   PMC sweep. **Treat the search as complete.**

**Two donor tables were built anyway, as TRAINING assets for 161:**

| built | tissue | targets | panel | shared targets vs held donors |
| --- | --- | ---: | ---: | --- |
| `nadig_hepg2` | **liver** | 2,160 | **0/300** | 2,143 k562 · ~1,950 Orion · 2,086 Jurkat · 2,039 RPE1 |
| `frangieh_melanoma` | **melanocyte** | 235 | 5/300 | 202 Orion · 156 k562 |

HepG2 is an 11th line, a **new tissue**, and an order of magnitude denser than Jiang's 171 shared
targets. **It does not move the effective-donor count by one unit** — 159's `1/k` law counts
lines carrying the panel *at inference*, and these carry 0 and 5. Production stays at 3.1.

**So: 161's residual model is the remaining path, and the coverage-vs-training distinction is the
operative one.** Do not spend further effort searching for donors.

### 166 also corrected 159, and the correction matters for how its numbers are quoted

159's line-pair estimator caps targets by a random draw **that also sets the gene axis**, so each
published row is one realisation. `orion_hct116 × orion_hek293t` reads **0.0285 ± 0.0011 over 12
draws** against 159's published 0.0306, and **+0.0231 on all 17,698 targets**. 159's conclusions
hold; **its values need a ±0.001–0.002 error bar.** That is the second independent correction to
159's magnitudes (the first was its 86.5%/86.9% harness-blob headline, which re-measured at
26.4%/44.1%). **Quote 159 for its orderings, never for a bare magnitude.**

---

## 2026-09-02 — 168: sharpening reopened, priced, and CLOSED. Plus the night's most important rule.

Experiment 030 measured magnitude sharpening (`sign(P)·|P|^q`) at **+0.072 macro PDS**, additive
with whitening, and parked it in August because the deciding axis (`nmae`) was only reachable
through a proxy. The instruments now exist. Reopened, and it is a **TERMINAL NO-GO**.

The full record is `experiments/168-sharpening/RESULTS.md` (168, 168b, 168c, 168d). The arc, and
it is worth reading as a process case study because two of my own claims were withdrawn:

1. **168** re-priced the global power on the Orion pair with exact `nmae`: q = 1.15 gave
   worst-surface Δpds **+0.01145** for Δnmae **+0.0036**, with a sharp `nmae` cliff past q ≈ 1.2.
   Marked INTERIM.
2. **168b** identified a mechanism — the top 5% of `|D|` carries **69.8% / 73.1%** of the
   aggregate alignment `Σ D·T` — and predicted a head-targeted arm, which measured
   **+0.02602** worst-surface on Orion, double the best global arm.
3. **168c** ran the control that could fail. `SHUFFLED` failed on both surfaces (good), but
   HEK293T was nearly indifferent to **REVERSING** the curve while HCT116 collapsed. The fitted
   curve *compresses* the magnitude range while the head arms *expand* it, and both appeared to
   help — so **168b's mechanism claim was withdrawn**. I proposed HEK293T's blob as the cause
   and **that was refuted too**: mean-removed scoring left REVERSED at +0.01897 there.
4. **168d** broke the tie on **h1-test** and closed the family.

### The closing numbers, and why h1-test is the only surface that should have been used

| arm | Δpds h1-test | Δpds Orion (worst) | |
| --- | ---: | ---: | --- |
| GLOBAL q=1.15 | +0.00112 | +0.01145 | 10.2x over |
| GLOBAL q=1.50 | +0.00784 | +0.02944 | 3.8x over |
| isotonic tau=0.5 | +0.00280 | +0.01845 | 6.6x over |
| **top 2% x2** | **−0.00042** | +0.02602 | **SIGN FLIP** |
| **top 2% x3** | **−0.00994** | +0.03452 | **SIGN FLIP** |

The tie-break favoured the method being *real* (h1-test: FITTED +0.00280, REVERSED −0.00378,
SHUFFLED −0.01737 — direction matters, agreeing with HCT116). **The magnitude is simply 4–10x
smaller than Orion said, and the head family inverts.** Every arm fails rule 13's +0.0019
shape-only threshold before the −0.0146 reach toll is charged, and the only arm with a real
h1-test `pds` gain (q = 1.50) carries an exact `nmae` cost of +0.2094 raw (≈ −0.058 overall).

## 🚨 STANDING INSTRUMENT RULE: never size an arm on the Orion pair. Use h1-test.

Measured four independent ways in one night:

| evidence | overstatement |
| --- | --- |
| 157: repaired HCT116 vs h1-test, tier-policy `reach` effect | **4.4x** |
| target-independent blob share: Orion 26.4% / 44.1% vs h1-test **3.6%** | — |
| 166: 159's line-pair values are one realisation of a random draw | ±0.001–0.002 |
| **168d: Orion vs h1-test on sharpening `pds`** | **3.8–10x, and the sign inverts** |

h1-test is also the only surface that can be scored through the **production** `union_pruned`
operator without leakage (it is not an Orion line), and the only one validated against a board
`reach` point. **Caveat that must travel with it:** `data/donors/h1_2025` covers **85/85** of its
targets and is pure oracle leakage — use `k562_gwps` + Orion only.

The Orion pair stays useful for arm-to-arm **ordering** on mean-removed statistics, and its
SHUFFLED/REVERSED controls did their job here. But a magnitude read off it can be wrong by an
order of magnitude and can invert. **This does not invalidate a23**: Orion *under*-stated that
small perturbation (predicted +0.00035/+0.00278, board +0.00075). The failure mode is specific to
**large reshapings**, which interact with the blob.

### Where this leaves the queue

Every lever I opened tonight is now closed with a mechanism: 160 (control geometry), 164
(call set), 166 (donor acquisition), 167 (target expression), 168 (sharpening), plus 157/158
earlier. **161's cross-line residual model is the last live thread**, and the coverage-vs-training
distinction 166 established is what makes it fundable: `nadig_hepg2` alone adds an 11th line, a
new tissue, and ~1,900–2,400 shared targets against four existing donors.
