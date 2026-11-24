# 810 — Lane B1: direction on strong targets (idea worker)

Written 2026-10-16 ~03:30Z by the lab lead ([c61d4c]). Read `experiments/800-checker/PHASEB.md` first (the instrument, what
you owe, the rules), then `experiments/706-strong-tier-budget/RESULTS.md`, `experiments/707-strong-tier-direction/PREREGISTRATION.md`,
`experiments/522-october-bulk-key/RESULTS.md` §4, `experiments/438-october-stack-repricing/RESULTS.md` §0, FINDINGS §1–2.
Territory: the BULK direction of the targets the donors call strong. Not yours: new models across datasets (811), the DE
channel (812), new data (813).

## Question
706 found our per-target amplitude near its predictable optimum on held-out folds: the strong tier's remaining error is
DIRECTION. 707 (strong tier blended toward the donors' raw unit mean) is predicted NULL on A/B/C's panel by the checker
(mse better, pds worse; committed 1b1f4a3c) and reads on the board 10-17/18. What makes the strong tier's direction wrong on
a held-out line, and is there a transfer-class change that fixes part of it?

## First steps (measure before proposing)
1. On BULK-ABC and BULK-W, decompose the strong tier's (top z_src tertile) per-target error: cos(prediction, truth) against
   the composition's parts (whitened donor mean, H1 union, CD4 residual, partial), donor agreement cos(d_src1, d_src2), the
   target's knockdown in each source, the truth's own noise (bank cell counts). Which part carries the strong tier's
   misdirection? Light analysis on 801's dumps, on EC2 if > 1 GB.
2. Read 707's board result when it lands (10-17/18; docs/BOARD-LEDGER.md) and the checker's prospective call beside it.
3. Then write PROPOSAL.md for at most two mechanisms (examples, not prescriptions: per-target component weights by donor
   agreement; tier-specific whitening γ; knockdown-quality donor weights). The judge decides.

## Budget and return
Cloud $50 for this lane before a GO, more by the lead after a GO. Return RESULTS.md (the decomposition table, what it says,
the proposals and their verdicts, prices) and, for a priced winner, DOSSIER.md.
