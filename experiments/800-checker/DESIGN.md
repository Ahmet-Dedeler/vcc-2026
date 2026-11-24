# 800 — the lab's checker: design (v1, 2026-10-15 ~22:40Z, lab lead [c61d4c]; committed 0e2b4e8e 22:47Z)

Charter: `docs/LAB-2026-10-16.md`. Log: `docs/LAB-LOG.md`. This file is the checker's contract. It changes only by a dated
amendment at the bottom, never silently.

## 0. What it is for, in one line
Give any candidate change a price on held-out data, in board units, member by member, with a target-bootstrap interval, a
count of how many tries its family has had, and a record of how often this same instrument was right on the board before.

## 1. The one command
```
uv run python experiments/800-checker/check800.py price experiments/8NN-x/candidate.json [--surfaces working] [--cloud]
uv run python experiments/800-checker/check800.py report <run-id>        # tables, bins, wear, calibration
uv run python experiments/800-checker/check800.py wear                   # configs seen per surface, N per family
uv run python experiments/800-checker/check800.py ka                     # the known-answer retrodiction table
```
`price` resolves the surfaces, writes one job per (surface, arm, seed) to the EC2 runner (`cloud/labcloud.sh`), collects the
per-target arrays, aggregates, appends the evaluation to the wear ledger and prints the report. Nothing is priced locally
above ~1 GB (charter: heavy compute on EC2; the Mac runs production builds).

## 2. A candidate is a file, not a description (AGENTS 11.1 rule 11)
`candidate.json`:
```
{ "name": "...", "family": "...",              # family = the switch it tries (tries are counted per family)
  "class": "transfer | scorer | abc",          # AGENTS 11.1c.1, declared before the numbers exist
  "prereg": "experiments/8NN-x/PREREGISTRATION.md", "prereg_sha256": "...",
  "base": "incumbent | <named arm>",           # what it is priced against (default: the October incumbent below)
  "bulk": { "module": "experiments/8NN-x/arm.py", "sha256": "...", "args": {...} },   # optional
  "de":   { "flags": [...], "module": "...", "sha256": "..." },                          # optional
  "predictions": { "uri": "s3://.../<fold>.npz" }                                         # optional: a whole model's rows
}
```
- **bulk.module** exposes `apply(fold, arm) -> arm'` on the fold's objects (composition rows `fin`, carriers, LODO key
  weights `w`, per-target sharpening `q`, frames). Everything else (render, ballast realisation, scoring) is the engine's.
- **de** is an emitter configuration on the real-emitter surfaces (run496/run618e flags, or a patched emitter module).
- **predictions** lets a new model hand in its own per-target rows for each fold; the engine blends or swaps them in as the
  module says, so a model trained without the held-out line is priced exactly like any other bulk change.
- The **incumbent** (the base every new idea is priced against) is the decided October recipe's bulk: LODO key, percentile
  cap P 96.64697522769474, tercilec sharpening (2.5, 2.0, 1.75), H1 ×2, CD4 ×1; DE: zs, CLEAN4 hook 0.05, UPRN, DNMAG
  (`docs/october-decisions.json`, read, never written).

## 3. Surfaces (fold registry: `surfaces/registry.json`, every entry pinned by sha256)
| set | what | folds | status |
|---|---|---|---|
| **BULK-W** (working) | 6 fresh Arc-like panels (§3.1) × held-out K562-GWPS / Orion HCT116 / Orion HEK293T, PLOLO (the held line is out of the donor mix, the operator, the key and the CD4/H1 composition's fit) | 18 | built by 801 |
| **BULK-ABC** | A/B/C's own panel, the same 3 cancer folds + kolf2 (438/522's) | 4 | worn; the board's panel, so the retrodiction surface |
| **BULK-Z** (sealed) | 2 more panels by the same rule, drawn and recorded by sha256 only; no operator, dump or score exists until the one access | 6 | 801 draws; lab lead opens once |
| **DE-H1** | h1-test (85 targets) and h1-trainval (618, 126 targets): the production emitter on H1, 3 seeds; the h1_2025 donor never enters | 2 | worn; board-validated (618 TRUSTGATE) |
| **DE-CD4** (A.2) | GWCD4i primary CD4 T cells, 10x Flex (the challenge's chemistry), ~500 cells per target: a fresh, unworn DE surface | 1–3 | 805, after Phase A's first pass |

### 3.1 The lab's panel rule "arc800" (AGENTS §2.1, as close to Arc's selection as public data allow)
- Pool: official axis ∩ carried by at least one cancer donor (efp's union pool), minus the K562 essential list, minus
  every target of every `pert_counts.csv` on disk (A/B/C, real, efp, efr, R6 S1–S4, R7a T1–T4, …) and earlier lab draws.
- **Coverage pattern matched to A/B/C** ("adequate coverage in public datasets"; A/B/C has K562 269/300): stratify by which
  of the three donors and the CD4 atlas carry the target, at A/B/C's shares.
- **Strength mix matched to A/B/C, measured on a line that is never a fold** (Arc binned by an H1 pilot screen; we bin by the
  GWCD4i Rest downstream-DE count, a genome-scale pilot in a fourth line; targets the atlas lacks are their own bin).
  Stratified draw: A/B/C's own share per bin.
- **Basal expression**: the target's control expression in the frames ≥ A/B/C's panel 5th percentile (Arc: "high basal
  target expression").
- 300 targets, 0 essential, seeds 800101…800108, drawn in order, each with its sha256, at draw time, by one script.
- Diagnostics recorded before any score: strength-bin shares, key/median distribution, floor share, K562 coverage, versus
  A/B/C and versus efr/efp (522 §4b's panel-mix test).
- **Amendment (2026-10-16 11:48Z, lead; judge's note on 814): the pool rule above BROKE a recorded rule.** FINDINGS §8.1g (and LOG
  l.200–212) state Arc's criterion: A/B/C's panel has 0 hits in the Replogle K562-essential, RPE1, Nadig HepG2 and Nadig
  Jurkat screens. A synthetic panel must satisfy it, not merely be disjoint. arc800 removed only the K562 list, so L1–L6 and
  the sealed LZ folds (L7, L8) carry 14–24 RPE1 ∪ Jurkat essential-screen targets each (814: 16.9 ± 3.7 expected per draw).
  Those targets carry most of the fresh panels' direction over-read (`docs/LAB-PANEL-REPORT-2026-10-16.md`; 814 RESULTS).
  **From now on every lab panel draw also removes `data/donors/rpe1/perts.txt` ∪ `data/donors/nadig_jurkat/perts.txt`
  (∪ a HepG2 list if it differs; 814 checks).** Existing L folds stay registered, but key and direction classes are not
  decided on them (PHASEB). The LZ access counts their essential-screen targets before its price is read.
- **Amendment 2 (2026-10-16 12:00Z, lead; 814 VERDICT-2):** the essential-screen cut is the screens' FULL perturbation library,
  `experiments/814-panel-rule/work/flags/esslib.txt` (2,394 genes, sha256 `24fffbfe…`). That is the set FINDINGS §8.1g's table
  counts, not the ≥ 20-cell donor banks (`rpe1/perts.txt` ∪ `nadig_jurkat/perts.txt`, 2,351). Inside the pool they differ by
  7 genes (CDK7, CEBPG, CMTR2, POLR3A, POT1, RBM22, ZNF830). 801's A1–A6 were drawn (run `lab801-draw-a`, 11:52Z) before this
  point reached the lead, with the 2,351-gene cut. They are NOT redrawn (a redraw would be a second try): their count of
  those 7 genes is reported as a recorded deviation. Every later lab draw uses esslib.txt.

## 4. Engines and what each member's number is
| member | bulk engine (BULK-*) | DE engine (DE-*) |
|---|---|---|
| pds | literal PDS on the fold (scorer's rule, panel genes out), and 414's board-emulated λ 0.5 | real (cell-eval2, or lib521 + KA vs cell-eval2) |
| mse | from the relative aggregate cosine r (518's identity: +1 % X ≈ −0.0041 raw; the levelled route), and, where the fold has cells, the scorer's numerator ratio | real |
| fid, jac, reach, nmae | **unmeasured** by the bulk engine. Mechanism that moves them for a bulk change: ballast K (re-roll class, sd 0.00025 overall). Bulk candidates that change K by > 5 % or reshape the low tercile are re-run through DE-H1 (479sp2's lesson: its board loss was nmae via ballast) | real-emitter members at f 1 and 618's board multipliers (619) |

Every report lists the members it measured and, for the rest, the mechanism that pins them (AGENTS 11.1 rule 12).

## 5. Units
Board units: Δscaled_m = Δraw_m × INV_m; Δoverall = Σ_m Δscaled_m / 6. INV is the board's own effective inverse span,
fitted on ledger pairs: pds 2.1614, mse 1.0466 (702 vs ALLSW; nmae/fid/jac/reach from 192: 1.64787 / 3.32482 / 2.72081 /
1.12927). A fold's raw Δ is read as a board raw Δ; the over-read ratio (§7) is what corrects for the fold not being the board.

## 6. Statistics
- **Paired target bootstrap per fold** (the same resampled targets for both arms; 4,000 draws; seed fixed per run), means
  re-computed per member (mse as a ratio of sums), averaged over seeds first on DE surfaces.
- **Pooled**: equal weight per fold; each draw resamples every fold independently. Report mean, 95 % percentile interval,
  P(Δ ≤ 0), folds positive n/N, and the leave-top-3-targets-out mean per fold.
- **Price from the lower end** (702 rule 1): the ship number is the 2.5 % bound after the best-of-N shift (§8).

## 7. Calibration: the over-read ratio, measured, not assumed
- Every known answer (§9) is run through the same frozen protocol. Over-read ratio ρ = checker Δoverall / board Δoverall,
  per class (bulk-shape, bulk-key, bulk-external, DE-magnitude, DE-sign). The board projection of a new candidate is
  Δ / ρ_class with ρ's leave-one-out spread added to its interval.
- A family with no board point gets ρ's full range across classes, so its interval includes zero (702 rule 4).
- 702 rule 3: every bulk report gives the share of the candidate's change energy outside the fold's gene axis and by
  control-CPM bin; above 25 % outside, the price is labelled PARTIAL.
- 702 rule 2: a surface the candidate has seen (its training lines, its own data) may bound a cost, never supply a gain.
  Candidates declare their training lines; the engine refuses a gain read on a declared-seen fold.

## 8. Wear and best-of-N
- `wear/ledger.jsonl`: one line per evaluation: UTC, candidate name, family, sha256 of candidate.json and every module, the
  surfaces (fold ids), engine version (git blob of the engine files), and the summary numbers.
- `check800.py wear`: per surface, distinct configurations and families seen; per family, N = distinct configurations
  priced on that surface set.
- **Best-of-N correction (AGENTS 11.1c.4):** a family's reported winner is shifted down by e(N) × its bootstrap sd, e(N) =
  E[max of N standard normals] (N 1: 0; 2: 0.56; 3: 0.85; 5: 1.16; 10: 1.54; 20: 1.87; 50: 2.25). The report shows both.

## 9. Validation (AGENTS 11.1 rule 3): the known answers come first
| row | board Δ overall (vs) | channel | must read |
|---|---|---|---|
| 495 cap 6 (LODO key cap 3 → 6) | +0.0016340 (446bj2) | bulk key | positive, ship |
| 498p2g (tercile sharpening on cap 6) | +0.0010083 (495) | bulk shape | positive, ship |
| 479sp2 (strength-tercile sharpening) | −0.0007141 (446bj2) | bulk shape (+ ballast nmae) | not ship |
| 702 (PIE blend w 0.25, seen targets) | −0.0025464 (ALLSW) | bulk external | not ship (HARM ideally) |
| 446bj2z zs | +0.0005493 (446bj2) | DE magnitude | positive, ship |
| 496a CLEAN4 hook 0.05 | +0.0005552 (446bj2) | DE sign | positive, ship |
| DN d_m200 | +0.0003736 (630arm) | DE magnitude | positive (PARTIAL) |
Extras where the engine can express them: 479p2g (+0.0010364), 630arm UPRN (+0.0009297), idall3 (−0.000501 vs the plain
hook), K20 / K16 call-set size (−0.0007 / −0.0021), 458d programme prior (−0.0085), 417a shape (+0.0070), 471c (~0).
- **Hit:** the checker's ship decision (lower bound > 0) and sign agree with the board. **Over-read ratio:** ρ per row
  with |board Δ| ≥ 0.0003. The checker is trusted only if it would not have shipped 702 or 479sp2 and ships ≥ 4 of the 5
  positives; otherwise the report says which class it cannot price.

## 10. Sealed reserve (BULK-Z)
Drawn and recorded (targets, seed, sha256) in `docs/CONFIRMATION-LEDGER.md` under "Lab reserve LZ", nothing built. One
access, for one preregistered final candidate, by the lab lead: claim file hashed into the ledger first, then operators,
dumps and the price, then the result appended whatever it is.

## 11. Cloud
`cloud/labcloud.sh`: code tarball + an S3 mirror of the data closure (`s3://…/input/lab/repo/<path>`), a job list run in
parallel on one instance, outputs to `s3://…/output/lab/<run>/`. Every instance and volume: tag Project=lab (charter
Amendment 1), on-demand for jobs over 30 min, per-launch cap $50 and per-UTC-day cap $150 enforced by the launcher, `instance-initiated-shutdown-behavior=terminate`, a user-data `shutdown -h +MAX` dead-man timer armed before
anything else, an EXIT trap that powers off; `scripts/aws_reap.sh` after each batch. Every launch and its cost go to
`docs/LAB-LOG.md`. Budget ≤ $800 through Nov 5.
