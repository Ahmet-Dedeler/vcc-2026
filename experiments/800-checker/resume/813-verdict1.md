FROM THE LAB LEAD (2026-10-16 09:16Z): the judge's verdict on your PROPOSAL.md is **REVISE**. It is in
`experiments/813-public-data/VERDICT-1.md`. Read it in full. Your context is intact. The judge calls the idea "the right
next data move for the lab", so this is a fix list, not a rejection.

Fix all of R1–R7, meet C1–C4 ("Conditions that will carry into a GO"), and resubmit ONCE as a new frozen file,
`experiments/813-public-data/PROPOSAL-2.md`, that nobody edits after it is written. Keep it short. My calls on the open choices:

- R1 (KA-SG): use split-half reliability-corrected or Deming/TLS slopes, with the single-guide cells disjoint from A's carriers.
  Run the statistic on KA-SIM's simulated data and set [lo, hi] from that pass rate BEFORE any real row is read. That
  simulation is sub-second-to-minutes and small; if it needs > 100 MB or > 1 s locally, put it on EC2 via labcloud.sh, or
  state it as stage-0 work. Freeze the statistic by sha (C3).
- R2 (P2 must be able to fail): take the cheaper honest route unless PIE inference is clearly unseen and cheap. If any
  TeloHAEC data could be in PIE's training set, PIE is not an unseen source here: drop 702 on ABC/telohaec. Then name which
  board-pointed negatives in `experiments/800-checker/DESIGN.md` §9 the engine (`experiments/801-lab-folds/engine801.py`
  `named()` arms) can express on a TeloHAEC fold. If none can, say so and rename the bins: no "TRUSTED" on two near-certain
  positives.
- R3 (fold build): for L1–L6/telohaec, a fresh 3-donor fit per panel with 801's day recipe (`--basis-prune essential`, the 12k
  cap). Take its time and cost from 801's ABCp refit logs (`experiments/801-lab-folds/`, runs `lab801-abcp-*` in
  `experiments/800-checker/cloud/runs.jsonl`). Give the production operator sha for ABC/telohaec. You do NOT edit 801's files
  (lib801.py, registry.json, engine801.py). Write the registry/lib change as an UNAPPLIED patch in your own dir with the KA
  that every existing fold's stored rows stay bit-identical. 801 or I apply it after a GO.
- R4 (P3's donor-count confound): drop P3 as a blob sentinel unless the 2-donor cancer-subset pricing costs < $2 extra.
  Either is fine; say which.
- R5 (decision role and leakage): preregister this rule. TeloHAEC folds are **reported apart** (like kolf2), never pooled into
  BULK-W, and never a veto on their own. A candidate that gains on the cancer folds but has a TeloHAEC CI upper bound < 0
  carries a written flag in its dossier, with a mechanism, and the coordinator weighs it. Reason: all six challenge lines are
  cancer lines and TeloHAEC is not. Leakage: on a TeloHAEC-held fold every candidate's TeloHAEC rows (pilot and genome-wide)
  are dropped. Specify the declared-seen-surface refusal for check800 in your proposal; the lab lead implements it in
  `experiments/800-checker/check800.py`. You do not edit it.
- R6 (novelty): cite 105 and answer it for (a): P2/R2 is the test of whether a fold this far from every donor still orders
  arms like the board. **Delete P4 / use (b)**: the genome-wide rows as an 811 bank stay 811's own judged, pre-counted try.
- R7 (data terms): make it a STOP gate: fetch the terms text, record its sha256 and your reading in RESULTS.md, and fetch no
  data byte if the terms forbid this use.
- Lane cap: state a total for the lane (stage 0 ≤ $5; the full build as priced in R3). Per launch ≤ $50; Project=lab; on-demand.
- C1–C4 verbatim as conditions in the proposal.

Then stop. The lab lead sends PROPOSAL-2.md to the judge. Nothing is downloaded and nothing is priced before a GO.
