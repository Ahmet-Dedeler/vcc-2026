# 657 — CLEAN4 octile 4 delivered UP, call-neutrally, priced on 618 (PREREGISTRATION)

Council session [7772dc], assigned by the coordinator after 625 Amendment 4 (ff3786b1) dropped 625's A/B/C readout. Written
and committed BEFORE any 657 run. Class: **transfer** (re-signing CLEAN4's weakest DOWN octile, from cross-lab evidence).
Code changes, if any ship, go to the post-10-08 apply list, never into the 10-06 chain.

## Why, and what is already known
- **Octile 4** = CLEAN4 P_down in [0.5058, 0.5251] (625's v8 edges), 793 covered genes. Its |P − ½| ≤ 0.025, below the hook's
  τ 0.05, so **today's hook never touches it**. The block prior signs these genes DOWN, with ~no confidence.
- **Agreement of CLEAN4's DOWN sign with each line's truth on octile-4 gate cells (octflip625b.py, and this file for 618):**
  - Broad-16 KO cancer lines: **0.399, wrong on 12/12 lines** (UP right at ~0.60);
  - H1 (h1-test): 0.502;
  - **618 (h1-trainval, the pricing surface): 0.524**, i.e. DOWN slightly right.
- **Consequence:** 618 is an H1-line surface where octile 4's UP sign is slightly WRONG. **It cannot show 657's value.** It
  can test call-neutrality and bound the downside on a line like H1. The value case rests on the 12 Broad cancer lines,
  and only an A/B/C read (Arc cancer lines) can test it.
- **Tries:** octile 4 was found by scanning 8 octiles (625 / octflip625b), so this is the **best of 8**. Its Broad 12/12 is
  partly selection.

## Routes (priced separately)
- **A, hook route:** the c4h hook's sign is set to +1 (UP) on every covered octile-4 gene of the gated axis. Elsewhere it
  is c4h's τ 0.05 rule, unchanged. The hook acts only on uncalled genes (rc_emit_k: `~head & ~blk`), as the MDE add's
  direction. Mode `tvc4h_o4` of `run657e.py` (a copy of 618's run618e with this one override; diff recorded).
- **B, 630 route:** the hook is unchanged. 630's UP pass also selects the uncalled octile-4 genes as UP (`sel = uncalled &
  (lfc > 0 | octile4)`). A variant suffix `_o4` in `rnup630_657.py` (a copy of d4cd5f's rnup630 with this one change;
  diff recorded); target t = max(l, K·MDE) as for every 630 gene.
- **AB:** both.

## Arms (618: h1-trainval, faithful key, `--zero-repair -5`, zs, seeds 0, 1, 2)
- **Emissions (`run657e.py`, `--dump 126`):**
  - `c4h` (base) — KA: its stats == 618's `hk_c4h` stats bit for bit (q, p, lfc, sums; pre-zs and zs);
  - `c4h_o4` (route A).
- **630 stacks (`rnup657.py` = rnup618's wrapper on the copy), on each emission's dumps:**
  - `b_k100_q20r3p10` (630);
  - `b_k100_q20r3p10_D200` (630 + DN, the coordinator's (b) stack);
  - each also with `_o4` (route B).
- **Read** (paired per seed, all vs the same base stack):
  - A = c4h_o4 + 630 vs c4h + 630;
  - B = c4h + 630_o4 vs c4h + 630;
  - AB = c4h_o4 + 630_o4 vs c4h + 630;
  - the same three with DN.
  - Route A without 630 (emission stats) is reported as a component.

## Measures
- Members at 618's regimes f 1 / 0.5 / 0.35 / cmp, with 618's own pricing functions (score661oct / stacks618b pattern).
- Δcalls per target; zero-total counts (`zerocount654.py`); octile-4 cells' delivered sign and mean |lfc|.
- Board price: nmae × 619's 1.0–1.2; fid / jac by 611's call-set rule if calls move (they must not).

## Gates and decision (fixed now)
**Hard gates, per route, every seed:**
- |Δcalls| ≤ 10 per target;
- |Δfid| ≤ 0.0005 and |Δjac| ≤ 0.0005 raw at f 0.35 and f 0.5.

**Bounded downside:** overall-eq Δ ≥ −0.0005 at f 0.35 AND f 0.5 (mean over seeds).

**Outcomes:**
- A route that passes both gates and the downside bound is **ELIGIBLE** for one A/B/C board arm with its own prereg. That
  arm is the only test of the value.
- If both A and B are eligible, the larger f 0.35 mean goes to the board.
- Otherwise the route is CLOSED.

## Predictions
- **A:**
  - |Δcalls| ≤ 10 (call-neutral: a small MDE add on uncalled genes);
  - on 618, Δnmae in +0.0000 … +0.0010 raw (slight harm: the UP sign is right only 0.476 there);
  - overall-eq in −0.0003 … 0, so eligible;
  - with 630, the harm grows a little, because 630 amplifies the newly-UP genes.
- **B:** a smaller effect than A, because 630 cannot raise a DOWN-delivered gene much without recipients above the control
  max. Δcalls 0.
- **Value on cancer lines (not measured here):** the Broad agreement favours UP. The size is bounded by the hook's delivery
  (612: the hook flips 2–6 % of signs). Expect ≤ +0.001 on A/B/C if it transfers.
