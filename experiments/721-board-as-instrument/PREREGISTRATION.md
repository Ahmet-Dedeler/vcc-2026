# 721 — The public board as the instrument for the DE channel (preregistered 2026-10-15 02:40Z, before any fit)

**Question.** Our DE members (fid 0.0055, reach 0.233, jac 0.028 scaled) trail the 0.30 pack (fid 0.04–0.09) by ten times
on fidelity while our bulk is comparable. Our own instruments have twice been wrong in SIGN about the call-set size K
(K20, K16; 523). The board is the only instrument that knows the hidden quantity, the number of real DE genes per
perturbation (`n_conf`), and 18,634 scored entries from 1,841 teams across 252 snapshots are measurements of it. Is our DE
operating point dominated, what (call-set size, sign accuracy) explains the pack's fidelity, and does the implied truth
agree with the K-series reads?

**Raw material only:** `data/leaderboard/*.json` + the 10-14 capture (entries, raw and scaled members), the installed
cell-eval2 0.16.0 metric code (`metrics/direction.py`, `metrics/de.py`), Arc's spec constants (`docs/reference/
vcc2026-metrics-FULL.md`, b and r per context), and our three board entries whose call-set size is known by construction
(446bj2z K 2,200 `be90smZUDcIV7j5MgQ03`; K20 `kJfe2dP9LEJawDEf0pMm`; K16 `Imj3EXViNco8rT4eLLlE`). No FINDINGS numbers are
inputs; they are compared against at the end.

## Hypotheses, with the reading written first
- **H1 (dominance).** Among teams whose bulk sits within ±0.02 pds and ±0.03 mse (scaled) of our best, the median scaled fid is
  ≥ 0.03 and ≥ 70 % of them beat our fid. Reading: TRUE → the DE gap is an operating-point gap that a bulk-neutral change can
  close; FALSE (median < 0.015) → the pack's DE comes with its bulk and the gap is model quality.
- **H2 (regime).** The forward model transcribed from cell-eval2 says a submission's `fid_raw = mean_p[ acc_p · min(1,
  yield_p) ]` and `jac_raw = mean_p[ TP_p / (K_p + n_p − TP_p) ]`. For the 0.30 pack's (fid, jac) pairs, the feasible set
  under the inferred `n_p` distribution requires call sets of ≥ 2× ours (K ≥ 4,000) at sign accuracy ≤ 0.60, i.e. their
  fidelity is yield, not sign quality. Reading: TRUE → the lever is K (and the ordering that protects reach); FALSE (feasible
  only at K ≈ ours with acc ≥ 0.65) → the lever is sign accuracy and K is not the problem.
- **H3 (consistency).** The `n_p` distribution fitted on the field plus the K 2,200 anchor predicts the SIGN and roughly
  the size of our K16 and K20 reads (fid −0.0041 and −0.0017 raw; jac −0.0003 and +0.0004 raw) without using them. Reading:
  both signs right and sizes within 2× → the model is usable for a K decision; otherwise it is not.
- **H4 (our accuracy).** Our sign accuracy on called ∩ real genes is ≥ 0.60 (we call few, precisely); the pack's is
  0.53–0.58. If TRUE with H2, the final-recipe decision is a preregistered K dose arm upward (K ≈ 3,000–4,000) read on the
  board before Oct 21, with jac and reach costs priced from the same model.

## Analysis plan (fixed now)
1. Transcribe the three DE metrics from the code; verify the transcription on a synthetic truth through cell-eval2 itself
   (bit-level on fid/jac, 1e-6 on reach) before any inference.
2. Recover the panel-level scaling constants from (raw, scaled) pairs across all entries and check them against the spec
   (fid b 0.505–0.522 / r 0.795–0.832; jac b 0.021–0.037 / r 0.375–0.423; reach b 0.047–0.097 / r 0.958–0.978).
3. Infer `n_p`'s distribution (log-normal, two parameters) from: the baseline point (recall 1 by construction: `b_jac` =
   mean n_p/|gate|), the K 2,200 anchor, and the field's (fid, jac) envelope; report bootstrap intervals over entries.
4. For every latest entry, compute the feasible (K, acc) set; summarise the pack; locate ours.
5. Validate: H3 on K16/K20 (held out), the spec constants (step 2), replication by another session from the raw JSON.
6. Decision rule: adopt a K arm for the board only if H2 and H3 both read TRUE; otherwise record and stop.

## What would falsify the whole approach
If the field's (fid, jac, reach) triples cannot be reproduced by ANY (K, acc, ordering) under a single `n_p` distribution,
the panel-mean members hide per-perturbation structure the board cannot resolve, and the instrument is void for this
question. That too gets written down.
