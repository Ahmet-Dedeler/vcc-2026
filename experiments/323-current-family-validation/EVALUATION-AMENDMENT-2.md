# Implementation amendment 2 — September 9, 2026

The reference pass completed, but the known-answer pass stopped before emission
because the cached gene strings had NumPy object dtype. The cache deliberately
loads with `allow_pickle=False`.

Normalize all decoded labels to NumPy Unicode in the reader and extend the
encoding test to assert non-object dtype. Recompute the same deterministic
reference statistics; retain the first pass's JSON and log and compare its
DE counts with the corrected-cache pass. No reference cohort, candidate,
seed, metric, tolerance or gate changes. No candidate score has been computed.

Failure evidence: `work/known.log`; original implementation in commit `f22033be`.
