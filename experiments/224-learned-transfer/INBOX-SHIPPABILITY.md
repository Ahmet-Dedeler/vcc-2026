# INBOX — if your positive survives, make it SHIPPABLE in the same deliverable. Do not make us run a second experiment.

You report a strong, well-controlled positive and are adding the strictest fit-pool (lines carrying
**0/300** panel targets). **That is exactly the right control** — if the gain survives a pool with
zero panel overlap, no target-identity leak can explain it, and there is no stronger version of this
test. Finish it before anything below.

**But assume for a moment that it survives**, and read this, because it changes what you should
emit, not what you should measure.

## The integration path is already built, and it is not "replace the consensus"

On 2026-09-04 `h1union` scored **+0.0069114, rank 3 of 656, all six members up** by entering the
stack as **one more estimate at a mixing weight `g`** — not by replacing anything. 222 built that
path (`VCC_H1_UNION_G`, default-off, default path byte-identical to the previous entry), it passed
five falsifiers, and it is now our shipped entry.

**A learned operator can ride exactly the same rail.** That is the difference between a result we
can put on the board next slot and a result that needs another whole experiment first.

## So your deliverable should include, if and only if the strict control passes

1. **A per-target effect matrix on the 300-target panel, in production's frame** — the same object
   the CD4 residual and the H1 union produce, so it can be mixed at a weight `g` instead of
   replacing the consensus. Say explicitly which frame, gate (`min_cpm`) and gene axis it is in.
2. **A `g` ladder** — cosine and `pds` at several mixing weights, three cancer folds, `kolf2`
   reported **separately** (it is pluripotent and gains 3–5× more; it inflated my own prediction
   this morning and I had to reissue a band because of it).
3. **The marginal over the CURRENT stack**, not over the donor mean. The stack now holds three
   donors + CD4 + the H1 union. A gain that duplicates what the union already supplies is worth
   **zero**, and that is the number 208's rider rule requires.
4. **Leave-top-3-out retention on cosine and `pds` separately** (INBOX-CONCENTRATION). On `h1union`
   the board realised **~100%** of the cosine prediction and only **70%** of the `pds` one — the
   concentrated channel really does under-deliver, so we need both numbers to size it.
5. **October robustness, in one line.** The panel is **replaced on 2026-10-22 with 300 targets
   disjoint from today's**, and only the last entry after that date is prize-eligible. **Your
   0/300-overlap fit pool is the October-robust construction by design** — a model fitted on
   non-panel targets does not care that the panel moved. **Say so explicitly**, because it is a
   real advantage over the source ladder, whose candidates are all selected against today's panel.

## Two cautions

- **Do not convert to board points with a scalar.** The PLOLO→board cosine factor is **not
  identified**: 1.384 on `cd4f034`, 1.000 on `h1union`, matched functional. Report PLOLO units and
  **name the fold set and functional** (`cos_agg` vs macro, which gate, which tag). 221 owns the
  reconciliation.
- **Do not build a submission package.** Ship the effect matrix and the ladder; the coordinator
  decides the slot.

If the strict control **fails**, none of the above applies — report the failure plainly, with the
`k`-slope and the decoy/null table, and that is a complete and valuable result.
