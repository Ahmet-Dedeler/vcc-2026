# INBOX — add a concentration test to every candidate you price. 222 just showed why.

222 finished and found something that changes how a gain should be reported. On its `h1union`
arm, against three cancer folds:

- **Δ`pds` is concentrated: the top-3 targets per fold carry 99.8% of the mean, and
  leave-top-3-out retains only 14.6%.**
- **The cosine gain is NOT concentrated: leave-top-3-out retains 72.7%.**

Same arm, same folds, two channels, opposite fragility. A headline number alone cannot tell those
apart, and **a gain carried by three targets out of ~290 is a different object from a gain spread
over the panel** — the first is a handful of lucky retrievals and will not survive a new context;
the second is transfer.

**So: for every candidate that clears its gate, report leave-top-3-out retention on BOTH cosine and
`pds` separately.** It costs one extra evaluation per survivor and it is the difference between a
result that transfers to D/E/F and one that evaporates.

Two consequences for how you rank things:

1. **A smaller distributed gain beats a larger concentrated one.** If you must choose what rides a
   slot, rank by the retained fraction, not the headline.
2. **Concentration is diagnostic, not only a warning.** If a source's gain concentrates on targets
   that are in that source's own high-coverage set, it is measuring coverage; if it concentrates on
   targets unrelated to it, suspect the instrument.

For **225** specifically: 223's DLD-1 (+0.00228) has never had this test. **Run it on DLD-1 before
recommending it for a slot** — its two-library falsifier established that the gain is real, but not
that it is distributed, and those are different claims.

For **224** specifically: apply it to any rung that clears. You already caught one contaminated
hyperparameter selection through the `k=4` collapse — this is the same class of check on the
output side, and a learned operator whose gain sits on three targets is memorising, not learning.
