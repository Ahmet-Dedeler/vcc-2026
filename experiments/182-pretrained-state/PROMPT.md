You are starting experiment 182 in ~/Code/virtual-cell-challenge-2026 (Arc Virtual Cell Challenge 2026). Work only in `experiments/182-pretrained-state/`; downloads under `data/external/state-pretrained/` (gitignored). **No builds, no `vcc`, no `cell-eval2`, no submission.** Stay under 40 GB of downloads and 20 GB RSS. `uv run python`, never bare `python`.

READ FIRST: `AGENTS.md` §4, §7 (State is allowed; confirm the non-commercial `MODEL_LICENSE.md` permits challenge use when you download it), §11.1 rules 1–13 incl. the 181 amendments (retrodict any board point of the family; size on h1-test not Orion; `h1_2025` is oracle leakage on h1-test); `docs/CURRENT_STATE.md` from the 2026-09-03 00:00 UTC `abc-knee` result downward — **our predicted sign is right only 53.7% of the time on the reference's DE genes, every downstream lever is closed, and the predictor's per-gene sign accuracy is the one remaining problem**; `experiments/106-state-baseline/RESULTS.md`, `110-*/RESULTS.md`, `113-*/RESULTS.md` (earlier State tests — ALL on a small `state_sm` trained locally from the Colab recipe, pds 0.50 = chance; 106 built a working State environment with an MPS patch — reuse it); `experiments/175-donor-pathway-match/RESULTS.md` (the per-donor sign-accuracy instrument, `perdonor175.py`); `experiments/161-crossline-residual/RESULTS.md` §0 (our LODO consensus on held-out K562: pds 0.885, agg cos 0.140 over 1,883 shared targets).

## Why this experiment exists

Arc published PRETRAINED State Transition checkpoints in Jan–Feb 2026 that nobody here has examined. HuggingFace: `arcinstitute/st-x-replogle-full` (and siblings `ST-Parse`, `ST-Tahoe`, `SE-600M` — inventory first). These are trained on four genome-scale CRISPRi lines (K562 gwps, RPE1, HepG2, Jurkat) with held-out lines. Our 106/113 tests used a tiny locally-trained `state_sm` and got chance pds — that is NOT a test of the published pretrained weights.

The remaining problem is per-gene sign. If a pretrained ST, conditioned only on a recipient's control cells + a target gene, produces signs that beat our 0.537 on a held line's DE pool, that is **new information** we do not currently have. If it does not, public pretrained State is closed as a donor.

## What to do

1. Inventory HuggingFace `arcinstitute/*` State checkpoints. Record id, size, training data, held-out lines, licence. Confirm challenge use is allowed before downloading weights.
2. Reuse 106's State env if it still works (`data/external/state-vcc/` or wherever 106 put it). Do not rebuild from scratch unless broken.
3. Download at most ONE ST checkpoint that was trained on CRISPRi / Replogle-Nadig and has a held-out line we can score. Prefer `st-x-replogle-full` if it exists and fits the 40 GB cap.
4. Run control-only inference on **h1-test** (or a held Replogle line that is NOT in the checkpoint's train set). Do not mix `h1_2025` into any donor comparison on h1-test.
5. Compare against our production consensus (K562+Orion, no H1) on the **same gene pool**: sign accuracy on the reference's DE genes, and pds in the challenge frame if you can do it without `cell-eval2`.
6. Known-answer row: on a line that IS in the checkpoint's train set, State should beat chance by a lot. If it does not, the inference path is broken — abort and write that.
7. Write `RESULTS.md` with a one-sentence verdict, the numbers, and the mechanism. Negative is a result.

## Verdict rule

- Sign accuracy on the held DE pool ≥ 0.56 (our K562-alone is 0.567) → interesting, price a blend.
- pds on a held line ≥ 0.55 → interesting.
- Chance (pds ≈ 0.50, sign ≈ 0.50) → NO-GO, close pretrained State as a source.
- If licence forbids challenge use → stop, write that, do not download weights.

Do not modify `scripts/`, `src/`, `AGENTS.md`, or production. No submit.
