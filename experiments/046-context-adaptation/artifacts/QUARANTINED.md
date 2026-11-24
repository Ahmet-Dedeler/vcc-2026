# Do not use this Jurkat artifact

The artifact in this directory was created after the first outer-fold screen
passed, before the required label-dependence controls ran.

Those controls failed decisively:

- correctly paired map: **+0.04209** macro PDS;
- permuted target labels: **+0.05656**;
- best wrong-target map per fold: **+0.05786**.

The improvement is therefore not evidence of supervised target-context
adaptation. The files are retained only so the failed branch is reproducible.
They must not be used for context A or a challenge submission.
