# 552 — Reserve R7 (running record)

Owner: session "Explore contribution areas". Declaration: `RESERVE.md` (+ amendments 1–2). Ledger: `docs/CONFIRMATION-LEDGER.md`,
"Reserve R7". Purpose: one sealed confirmation of the FINAL October recipe (bulk on R7a, DE on R7b), after R6's last access
(C10b, 551) and after the DE switches settle (not before the 10-07 reads), on the coordinator's message.

## Status (2026-09-30 22:05Z)
- **Declared** 21:43Z (b4dc26f9), before any draw. Amendment 1 (aae29b72): no mock / uninfected cells can be identified in the
  Sunshine file (no condition column, no viral features), so the infected state stands; asymmetric veto-only-on-HARM rule.
  Amendment 2 (b9fca152): control-only, prediction-side inputs are not reads.
- **R7a drawn** 21:44Z (85afb153): T1–T4, 300 targets each, 0 essential, disjoint from R6 and the 11 panels on disk;
  `work/draws.json` sha256 `418bb605…`. Controls are clones of 199's relabelled files. No operator, fold or score exists.
- **R7b surface:** `build_r7b.py` (construction only) is queued behind 551's efp KA and gate551 (`wait552.sh`); the machine's
  free+inactive sits near 20 GiB tonight because desktop apps hold ~12 GB in the compressor.
- **R7b instrument:** worker 553 (local Opus 5.5) generalises the real-emitter DE harness to a surface; KA on h1-test; the
  access script is written and tested on h1-test, never run on R7b. Deadline 10-06 12:00Z.

## Status update (2026-10-01 08:40Z)
- R7b surface (unfloored) built 00:27Z: 107 targets, 27,396 cells (ledger).
- **553 found R7b as built unfit for the production emitter** (bimodal control libraries, q75/q25 7.24; 107/107 placements
  fail; steering dose 0 under any class rule). **Amendment 3** (70affa19): floor 5,000 official-axis counts for controls
  and target cells, BST2 excluded, < 30 cells dropped. **Floored surface built 08:38Z** (f402a8c4): 106 targets, 18,759
  cells, 2,378 controls identical to 553's floored prediction side, 59–556 truth cells per target (median 127).
- 553's harness: all H1 known answers bit for bit; access tested end to end on h1-test; R7b-shaped rehearsal band sd
  0.0005–0.0015 overall-equivalent. Access: after the 10-07 reads, on the coordinator's message, with board-retrodiction rows.

## Candidates considered for R7b (metadata only; RESERVE.md §3 has the reasons)
| dataset | line | modality | depth (median UMI) | targets (≥100 cells) | verdict |
|---|---|---|---:|---:|---|
| SunshineHein2023 | Calu-3 (lung adenocarcinoma) | CRISPRi | 17,626 | 140 (128 carried by ≥ 2 donors; 22 essential) | **chosen** (SARS-CoV-2-infected) |
| XuCao2023 | HEK293 | CRISPRi | 2,155 | 169 | too shallow; HEK293 ≈ a donor line |
| FrangiehIzar2021 | melanoma | Cas9 KO, co-culture | 10,988 | 231 | knockout under immune co-culture; immune targets |
| Myllia / Altair A549, Jiang 5 lines | A549 etc. | CRISPR | — | — | used as banks; A549 may be an October line |
| RPE1 (Replogle) | RPE1 | CRISPRi | — | — | used by 430/516/517/521 |
| McFarland, Chang, Zhao, Gehring, Lotfollahi | various | drugs | — | — | not genetic |
| Wessels 2023 | THP-1 | Cas13 | 8,782 | 5 | too few targets |

## Requirement for R7b's claim preregistration (noted 2026-10-01 ~00:45Z, from worker 537's 10-01 read)
- **Board calibration point:** 446bj2zcp (the call-propensity block) read NULL on 10-01: −0.00079 overall vs 446bj2z; raw fid
  −0.00358, jac +0.00399, reach −0.00335. Against H1's real-emitter price the board gave jac 0.48× its predicted gain,
  fid 1.95× its predicted cost, reach −0.38×; 528's public-fold capture screen (K562 / HEK293T / HCT116) had predicted
  near-full transfer (0.87× H1's capture gain). Source: `experiments/537-arm-pipeline/RESULTS.md` §0 "10-01 READ", §4.
- **So R7b's claim must carry board-retrodiction rows** (AGENTS rule 3, "validate the validator"): inside the single access,
  before the claim's own rows are read, R7b prices the knobs the board has already read (the block: `lcap-kept` vs `cpm`; zs;
  and the hook if 496a has read by then) and reports each member's ratio to H1's price next to the board's (block: jac
  ~0.5×, fid cost ~2×). The preregistration fixes how far R7b may sit from the board's pattern before its claim rows are
  labelled "instrument not board-like" (reported, never used to veto). 553's access script already prices these MODEs.

## 553 final (2026-10-02 01:10Z): the R7b instrument is ready (`553/RESULTS.md` v6)
- Access: `access553.py r7bd5000 --owner-go <claim> …` runs the whole read in one command (~25 min, < 4 GB); its dry run
  verified the floored surface's 5 outputs by sha256, the target list and the controls against the prediction side, and
  the 10 control halves (10-02 00:54Z). The truth has not been opened.
- **The band to preregister is `band_total`** (emission seeds + split-half truths + control re-draws: the control sample
  moves an R7b-shaped reading as much as everything else together). On the h1r7 rehearsal `band_sd` 0.0004–0.0013 vs
  `band_total` 0.0006–0.0018 vs a between-draw spread 0.0004–0.0031.
- **Retrodiction must be against the SHAPE-matched rehearsal, not H1-full:** R7b's shape alone (H1 cut to R7b's cell and
  control counts) already gives the call-propensity block a board-like pattern (jac 0.45×, fid 1.4×, reach −0.05× of
  H1-full; board 0.48×, 1.95×, −0.38×). `--retro work/access553_rehearsal5` reports each member over the rehearsal mean.
- Shape-matched expectations: CP +0.0014, CPH +0.0015, K 2,000 −0.0016 overall-equivalent; the hook ~0.
- **Pricing other DE arms on a power-reduced H1:** `553/RESULTS.md`, last section (commands, cost, and why it is a regime
  instrument, not a board-shaped one: A/B/C's board has H1-full's 400 cells x 18,400 controls).
- **R7a's access machinery is ready (10-02, before C10b):** `551/r7a551.py` (plan / ka-cd4 / record / access) and
  `551/aggregate_r7a551.py` run 551's instrument unchanged on T1–T4 (its 28 imported files byte-identical to `e7be21ca`, the
  efp KA). `plan` passes (bundles == draws.json, no T output, not armed); KA-CD4 (efp) PASS; the aggregation's smoke test on
  efp's known-answer outputs runs (not a reading). The claim is drafted in `CLAIM-R7-DRAFT.md` (rules fixed; F and the DE
  arms filled after the DE reads); nothing is recorded in the ledger yet.
