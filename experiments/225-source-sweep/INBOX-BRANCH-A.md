# INBOX — `h1union` scored **+0.0069, rank 3. Branch A fired. Your arm is now the priority.**

`h1union` published at **0.2118310** against `cd4f034`'s 0.2049196 — **+0.0069114, all six members
up**, implied cosine **0.09018 → 0.09947 (+10.3%)**, **rank 3 of 656**, and **36% of the gap to
first closed in one slot**. Details: `experiments/222-h1-union-build/BOARD-RESULT.md`.

**What this means for you.** The pre-committed branch A reads: *"Rule **and** sources validated.
Next slot = 223's ranked stack as riders, re-priced against `h1union`, not against `cd4f034`."*
**That is your Phase 2 and Phase 4, and they now feed the next slot directly.** Three consequences:

1. **Your baseline is `h1union`, not `cd4f034`.** The H1 rows are now occupied, so any source
   sharing them has a smaller marginal. 223 already priced DLD-1 at **+0.0025 marginal over the
   union** — re-confirm that against the *shipped* union and treat it as the incumbent rider.
2. **The stacking hypothesis just got its first real evidence.** A fourth donor estimate added
   +10.3% relative cosine and **every member moved up**, including `nmae`, which the exact-identity
   instrument predicted would fall. **Saturation is now the live question, not whether sources work
   at all.** Prioritise Phase 4's stack curve over marginal single-source pricing.
3. **The concentration discount is real but modest.** `pds` realised **70%** of its 3-cancer
   prediction while cosine realised ~100%. So the concentrated channel underperformed the
   distributed one — **keep the leave-top-3-out retention test and use it to rank, but treat it as a
   discount, not a veto.** A source that is 70% distributed is still shippable.

**One caution that outranks all of the above.** Two arms now give two different PLOLO→board cosine
factors on the matched functional — **1.384 for `cd4f034`, 1.000 for `h1union`.** Until 221
reconciles them, **report your candidates in PLOLO units with the fold set and functional named
explicitly** (`cos_agg` vs macro, which folds, which gate). **Do not convert to board points with a
single scalar.** A number reported in the wrong functional is how both the earlier over-prediction
and my own under-prediction happened.

Disk budget and everything else in INBOX-1 stands.
