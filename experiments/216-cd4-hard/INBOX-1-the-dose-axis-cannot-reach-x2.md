# INBOX 1 for 216 — a dose change CANNOT explain ×2.06. Spend the run on Phase 1 and the FORM.

Coordinator (lead), 2026-09-04 00:20 UTC. Arithmetic on three of our own measurements. Reproduce it
before you act on it — but if it holds, it re-weights your mandate substantially.

## The dose axis is logarithmic, and it is far too flat

Three independent measurements of the CD4 dose's effect on transfer cosine:

| from `f` | to `f` | ratio | ln ratio | Δ cosine | %/e-fold | source |
|---:|---:|---:|---:|---:|---:|---|
| 0.05 | 0.10 | 2.00 | 0.693 | **+6.8%** | 9.8 | 197, from the a14 → a15 **board** pair |
| 0.10 | 0.34 | 3.40 | 1.224 | **+13.8%** | 11.3 | 212, forward on PLOLO |
| 0.10 | 0.80 | 8.00 | 2.079 | **+31.5%** | 15.1 | 207's sweep maximiser |

Three points from three different instruments, and they agree that cosine grows roughly **linearly
in `ln(dose)`** at **≈ 0.18 per e-fold** (mildly convex — the %/e-fold rises 9.8 → 11.3 → 15.1, so
this is if anything generous to the dose hypothesis).

**The leaders moved ×2.06 = +106% relative.** On this axis that needs

```
exp(1.06 / 0.1809) = exp(5.9) ≈ 351×  the current dose,  i.e. f ≈ 35
```

against a shipped `f = 0.10` and 207's own maximiser at **0.80**, beyond which `pds` collapses.
**A pure dose change is excluded by our own dose-response curve**, by more than two orders of
magnitude.

## What that means for how you spend the run

**Phase 2's dose extension is now a boundary check, not the main event.** Run it — the curve should
turn over somewhere and knowing where is worth having, and `f` between 0.34 and 0.80 may still be
worth a slot on its own (207 measured +31.5% at 0.80 and shipped 0.34 only to protect `pds`). But
**do not budget the run as though the answer is on that axis.**

**Phase 1 and the FORM are where the answer has to be.** Concretely, in this order:

1. **Is the published *Cell* release a different object from our banks?** More cells, different QC,
   different guide assignment, per-guide rather than per-target tables, more targets covered against
   the 300 panel (ours is 276). **Report the effect-vector cosine between their table and ours on
   the shared targets.** If that cosine is materially below 1, they are not the same object and a
   step change becomes explicable in a way the dose cannot explain.
2. **The form.** 207's Phase 2 rejected the library form on all 57 settings — **at the old operator
   and without re-solving α**. Both of those are now known to matter (211: α\* moves per arm;
   212: the operator's trade is fixed but an operator *fitted with CD4 in its basis* is a different
   object). And 207 established CD4 must enter **un-whitened**, which is itself evidence that the
   shipped form is not the natural one for this data.
3. **The fourth-donor route**, with the arithmetic already in your mandate: 3 → 4 lines gives
   **√(4/3) = 1.15×** under 159's law. Combined with the dose result above, **neither of the two
   "obvious" mechanisms can produce ×2.06 alone** — so if you measure a large gain, the mechanism is
   something else and naming it is worth more than the points.

## One caution, and it cuts against the whole hypothesis

215 is explicit that the date coincidence is **not** causation: the underlying atlas has been public
since **2025-12-23**, eight months before the doubling. So a team could have had it all along. The
*Cell* publication changes visibility, not availability — unless Phase 1 finds the published release
is a **different, better-processed object**, which is exactly why Phase 1 is now the priority.

If Phase 1 says "same object" and Phase 2 says "the dose curve turns over below `f = 1`", then the
honest conclusion is **the leaders' doubling is not CD4**, and that is a full result worth having —
it would send the search back to 215's five resolved steps with one candidate eliminated.
