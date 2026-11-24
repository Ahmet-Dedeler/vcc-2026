# Lab lead decision on 810 after the GO run (2026-10-16 08:55Z, commit time; first labelled 09:10Z)

decide810 returned **DE-H1 FIRST** (|dK_rel| > 5 % on 12 of 18 BULK-W folds; ≤ 2.2 % on the ABC folds). Past that gate the ladder
would CLOSE on C1 (sgu1's BULK-ABC point −0.00200 [−0.0053, +0.0011]) and C4 (sgu1 − 707 paired −0.0024 [−0.0048, −0.0002]).

- **DE-H1 is NOT funded.** DE side effects cannot rescue a negative bulk read on the board's panel (they would have to add
  ≥ +0.0024 over 707 with K moving ≤ 2.2 % on ABC), and the K trigger sits on BULK-W where the candidate already wins.
- **The family is CLOSED for the A/B/C-panel route**, with its mechanism: gating beats the closed global form on both surface sets
  (paired +0.0079 ABC, +0.0065 W), but on Arc's panel pds still eats 158 % of the mse gain.
- **One information read is funded (≤ $0.5): the same five arms on ABCp** (801's A/B/C folds with the DAY-recipe held-out operators)
  as soon as 801 registers them. Reason: the candidate is "undo the whitener", and BULK-ABC's folds whiten with 062's unpruned
  operator, which neither the board nor the day path uses (712). The W/ABC split (+0.0089 vs −0.0020) is the same divergence seen in
  522 / 706 / 707 / 811, and ABCp tells whether part of it is the operator. It is a read, not a reopening: decide810 is applied to it
  as frozen, and a DOSSIER would still need every G1 condition there.
- **The checker passed a new known answer here:** the closed global form (sgu1all) read −0.0099 on BULK-ABC with E 2.32, as 212 / 097
  found (STOP did not fire). Logged in LAB-LOG.
