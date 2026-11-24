# Phase B: what every idea lane gets and owes (shared preamble for 810–813)

Written 2026-10-16 ~03:25Z by the lab lead ([c61d4c]). Read with the charter (`docs/LAB-2026-10-16.md` incl. Amendment 1),
`experiments/800-checker/DESIGN.md`, `JUDGE.md`, and the known-answer report in `docs/LAB-LOG.md`.

## The instrument you price on (what it can and cannot see)
- **BULK-ABC**: A/B/C's own panel, PLOLO folds held-out K562 / HCT116 / HEK293T (+ kolf2 reported apart), through
  `experiments/801-lab-folds/engine801.py` (production's fold arithmetic). The board's panel, so the closest read of an A/B/C
  board effect. Known answers: key change ρ 1.46, P2g sharpening ρ ~3.5, external rows (702 PIE) ρ 0.93. **Its mse side
  tracks the board (ρ_mse 0.9–1.65); its literal pds side does not** (over-read ~5× for sharpening, wrong sign for the key
  change). Report the board-emulated pds (pdsl, λ 0.5) beside the literal one.
- **BULK-W**: 18 fresh Arc-like folds (L1–L6 × the 3 cancer lines). The D/E/F-transfer read. Its panels carry heavier key
  tails than A/B/C (801), so key-class changes over-read 4.75×; external rows under-read (702: 0.19).
- **Not priceable on bulk folds:** per-target sharpening re-assigned by a statistic (479sp2's class; its board loss was in
  the emission). Such changes need **DE-H1** (`experiments/803-de-engine/de803.py`, the production emitter on h1-test and
  h1-trainval, six members from emitted cells).
- **Interface:** a bulk candidate is an `engine801.Arm(w, q, fin)` built from the fold's objects (`F.fin` composition rows,
  `F.kv` key, `F.carried`, `F.targets`, `F.sources` — never the held line), priced with `engine801.price(F, arm, base)`; write
  the STANDARD-FORMAT npz (`experiments/800-checker/stats800.py` docstring; `experiments/804-known-answers/ka804.py` is the
  worked example) and aggregate with `experiments/800-checker/check800.py aggregate <run> --cand <arm> --base inc
  --candidate <your candidate.json>` (this logs your try in the wear ledger). Heavy runs on EC2 via
  `experiments/800-checker/cloud/labcloud.sh` (ka804's job + sync files are the template: `experiments/804-known-answers/cloud/`).

## Which folds decide which class (lead, 2026-10-16 11:14Z; `docs/LAB-PANEL-REPORT-2026-10-16.md`)
The lab's fresh folds (BULK-W, L1–L6) were drawn by the arc800 rule, and they over-read against an Arc-drawn panel (A/B/C's
own 300 targets, on the same recipe = LABC) by: **key 3.08× [1.91, 6.76], direction ≈ 6.5× (lower bound 2.0), shape 1.19×
[0.81, 1.81]**. So:
- **key and direction candidates** are decided on the Arc-drawn panels: BULK-ABCp and BULK-LABC (and D/E/F-panel folds once
  the lab builds them after the release window). BULK-W is reported beside them, never as the decision number;
- **shape candidates** may still be priced on BULK-W (factor ≈ 1), with ABCp / LABC reported;
- any new class states which of these it belongs to before its first number, with the reason.

## What you owe before any compute (JUDGE.md)
`PROPOSAL.md` in your dir: mechanism, class (transfer / scorer; A/B/C-specific is never funded), a prediction with bins on
BULK-ABC and BULK-W (and DE-H1 if it touches emission), the falsifier, novelty against `docs/FINDINGS.md` §2 (cite the closure
it is NOT), cost and try budget, and the day-of path (what D/E/F's own inputs compute it from). The lab lead runs the judge;
only GO gets cloud. Every configuration you price is a try (best-of-N); plan > 10 only with nested folds or a held-out
subset of folds.

## What you deliver
A finished candidate becomes `DOSSIER.md` (mechanism and class, the checker's price with CI and tries, its known-answer
record, an UNAPPLIED implementation with its KA). The lab lead alone talks to production. Negative results are written with
their mechanism (they close a family for October).

## Rules (binding)
Write only in your own `experiments/8NN-*`; no git; no production code; nothing under `data/`; never the sealed reserve
(L7, L8) or production's reserves (R6, R7a T1–T4, R7b); local work only sub-second < 100 MB aggregations, anything else via
`experiments/602-council-gate/gate602.sh` or EC2, never while `bash experiments/800-checker/window800.sh` prints WINDOW;
cloud within your lane budget and the launcher's caps.

## Why the bulk folds price mse but not pds (read before leaning on pds)
pds is scale-invariant per target and a ranking among the panel's ~300 truths: a target's score moves only when its prediction
overtakes or falls behind ANOTHER target's truth in cosine, so it reads the target-specific part of the direction and how
crowded the competing truths are, never the size of the error. mse reads the cross term ⟨P, T⟩ against ⟨P, P⟩ and ⟨T, T⟩: it is
smooth and additive over targets, and truth noise enters only through ⟨T, T⟩, a per-fold constant that cancels in a relative
change (r = ΔX/X). Two consequences are measured (804): (1) a key change is a per-target scalar, to which pds is blind by
construction; what the folds show is second-order emission non-linearity (log1p render, ballast mixing), which the board
realises differently (re-levelled contexts, real cells), so 495's fold pds came out with the wrong sign (+0.0014 vs board
−0.0013) while its mse matched (ρ 0.89). (2) A direction change (sharpening, swaps, external rows) does move pds, but by how
many neighbours it crosses, and the fold truths (a public screen's bank rows, its own noise and gene axis) are spaced
differently from the board's 400-cell Flex truths: the folds read sharpening's pds about 5× the board's. So: price direction
on mse (r) with its calibration, report literal pds and board-emulated pdsl beside it, and treat a pds-only gain as unvalidated.

## Post-hoc rules (not validated; they may not decide anything until a later board read confirms them)
- **Call-set size is read at f 0.61** (DE-H1), not f 0.35: chosen after the frozen f 0.35 read mis-signed K20 / K16 (at f 0.61
  both land on the board, ρ 0.82 / 0.89). Label every price that uses it "post-hoc f 0.61"; a call-set candidate needs a board
  read before this rule supports a SHIP.
- **Board-emulated pds (λ 0.5)** got the key change's sign right where literal pds did not; same status.
