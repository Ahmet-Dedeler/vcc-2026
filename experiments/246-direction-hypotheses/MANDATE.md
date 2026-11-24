# 246 — The direction channel: audit every hypothesis that is NOT closed, then MEASURE the best open one

You are a research subagent. Work ONLY in `experiments/246-direction-hypotheses/`. Read-only
everywhere else. **No `vcc`, no build, no `cell-eval2`, no production change, no slot proposal.**
`uv run python`. Peak RSS ≤ 10 GB. **MEMORY GATE before any fold load:** `memory_pressure` ≥ 40%
free AND `vm_stat` free+inactive ≥ 20 GiB; 241 and 243 are running, so strictly serial, and record
peak RSS with `/usr/bin/time -l`.

## The position
Board 0.2118310, rank 3 of 656; leader `cqawesome` 0.22431. The whole gap is the direction
channel: implied transfer cosine **0.0995 vs 0.1336**, raw `pds` **0.836 vs 0.870**, and we are
AHEAD on `fid`/`reach`/`jac` (`docs/BOARD-LEDGER.md`, 222 `BOARD-RESULT.md`). `AGENTS.md` §8.1d:
**signs pay, magnitudes do not** — a direction arm must move |t|-weighted sign accuracy, and its
`pds` must not fall (every cosine gain that re-adds SHARED direction pays on `pds`: the
self-residual is +4.61% cosine and −0.0078 `pds`, 218 §3.1; `cd4f034` paid −0.00365 raw `pds`;
`h1union` is the only arm that raised it, +0.00344). `docs/SLOT-PROGRAM.md` revision 2: **the
scarce resource is now hypotheses.** The closed list is in `docs/CURRENT_STATE.md` ("Closed, do not
re-open without new evidence") and in the RESULTS headlines — treat a family as closed only when you
have read the closing experiment's headline and it closes THIS construction, not a neighbour.

## READ FIRST
`AGENTS.md` §8.1–8.4 and §11.1; `experiments/239-knowledge-audit/work/headlines.txt` (every
RESULTS headline); 218 §3, §7 (the residual construction and its controls); 222 RESULTS (the
staging: `stage222.py`, `work/effects-on028.npy`, `work/effects-off.npy`, `work/plolo/`); 216
`six216.py` (the exact six-member instrument, and note the `nmae` identity has mis-signed the board
twice — report it as a magnitude bound, not a sign); 099, 209, 212, 213, 224, 228, 232, 234
headlines; 162's `plolo.py` (import verbatim).

## Phase 1 — the audit (read-only)
Produce a ranked table of direction-channel hypotheses that are **OPEN**, each with: the
mechanism; **which information it adds that the three donors + CD4 + H1 do not already contain**
(the filter from `docs/STRATEGY.md` 2026-09-04: a hypothesis that cannot answer this is dead);
why it is sign-accuracy rather than magnitude (§8.1d); the expected sign on `pds`; the closing
experiment if CLOSED (id + headline number, cited); a PLOLO test design with its known-answer row;
and cost. Adjudicate at least these, explicitly, with citations:
1. **Specific-subspace relaxation.** The self-residual (099/218: +4.61% cosine, −0.0078 `pds`)
   with the component SHARED across the 300 perturbations projected out before it is added — per-
   gene centering across perturbations, or removal of the top-k right-singular directions of the
   un-whitened effect matrix — so the cosine gain is kept and the `pds` cost (re-adding shared
   direction) is not paid. Check 209 (Orion-pair excess subspace) and 212 exactly: was THIS
   construction measured, or a neighbour?
2. **Per-target dose** for the CD4 and H1 residuals from the target's cross-donor agreement, in
   place of one global `f`/`g` (the per-group dose edit spec'd on 2026-09-04, +0.00099 — cite it).
3. **Recipient-basal gene gating beyond the 10-CPM floor** (an effect on a gene the recipient does
   not express cannot be real): check 006, 213 and the `min_cpm` history.
4. **Sign-consistent consensus** (carry magnitude only where ≥ 2 donors agree in sign): state
   precisely what 232/234 closed and whether this variant is inside it.
5. **Per-source dose for the union** (Pan-2026 H1 and `h1_2025` separately) and a joint (`f`, `g`)
   grid: check 229 and 216.
6. Anything you generate yourself that passes the information filter.

## Phase 2 — measure the top OPEN candidate
On 162's PLOLO instrument, **on top of the shipped stack** (`cd4@0.34` + union `g` 0.28 — reuse
222's staging or 218's `compose218.py` `merged()` construction), three cancer folds with `kolf2`
reported separately. Endpoints: `cos_agg` (primary, ±0.20% relative resolution per 204), macro
cosine, `pds` (**must be ≥ 0**), `q_w`. Falsifiers, all mandatory: target-shuffled control at every
dose; the row-matched control where applicable; non-carrying rows must not move; leave-top-3-
targets-out retention on cosine and `pds` separately (220's lesson). Then price all six members
through `six216.py`'s route and 218 §7.4's rate card, never carrying a member at zero (§8.4).
**Known-answer row first:** reproduce 218's `cd40.34+h1PAIRCARRY0.28` on the three cancer folds —
macro **+7.6803%**, `pds` **+0.004887** — to the digit (222's `repro222.log` did), else stop.

## Deliverable
`RESULTS.md`: the Phase-1 table (CLOSED/OPEN with citations), the Phase-2 measurement with every
falsifier's number, the six-member price, and one line: **BUILD-WORTHY / NOT** against the
corrected gates (+0.0019 shape-only, +0.0089 amplitude-changing, `AGENTS.md` §11.1 rule 13 and
§8.4). A well-supported "nothing is open" is a complete result. Then stop.
