# 344 — Kill check for the empirical-null call-set key (construction A of the 2026-09-09 emission review)

Registered 2026-09-09 22:58Z before any output exists. Truth-free: reads only the challenge control files
(`data/challenge/context_{A,B,C}.h5ad`, 46 non-targeting constructs × 400 cells each). No perturbation truth,
no reserve, no H1, no emission, no package. Purpose: decide whether the H1 four-arm test (baseline / null key /
shuffled key / CPM-matched key) deserves a preregistration and a build for a 2026-09-10 slot.

**Construction under test.** Replace the 2,200-gene CPM-ranked block of the rc_callset with the top 2,200 by
per-gene mean |z| across the context's null constructs, where z is the scorer's own Wilcoxon (rc_emit.wilcoxon_z,
323's normalization order) of one construct's 400 cells against the pooled 18,400 controls — the same comparison
the scorer runs for a perturbation, so self-inclusion (400/18,400) is deliberate and matches the reference.

**Readings and stop rules, fixed now (emission review §2):**
1. Median BH-significant genes (q<.05, over the 5-CPM universe) per null construct. **Stop if < 10 in every context**:
   there is no construct-level null to model.
2. Spearman(null key, control CPM `ref_lin`) over the universe. **Stop if > 0.95 in every context**: CPM in disguise.
3. Leave-one-construct-out: key from the other 45 constructs; recall of the held-out construct's own significant
   set by its top-2,200 versus the top-2,200 by CPM (head excluded: no |H| exists for a null construct, so both
   keys rank the full universe). **Advance only if the null key beats CPM in ≥ 2 of 3 contexts with a 46-construct
   bootstrap 95% interval excluding zero.** Also recorded, no verdict: the same for BH-significant frequency as key.

Otherwise: scoped negative, recorded here, no H1 test.
