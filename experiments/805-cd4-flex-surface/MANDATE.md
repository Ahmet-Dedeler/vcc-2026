# 805 — a fresh, unworn DE surface in Flex chemistry: GWCD4i primary CD4 T cells (worker mandate)

Written 2026-10-15 ~22:52Z by the lab lead ([c61d4c]). You own `experiments/805-cd4-flex-surface/` and write only there.
Read first: `docs/LAB-2026-10-16.md`, `experiments/800-checker/DESIGN.md` (§3 DE-CD4), then
`experiments/618-h1-trainval-surface/{PLAN,RESULTS,TRUSTGATE}.md` (the model for this whole job), `run618.py`, `run618e.py`,
`truth618.py`, `experiments/615-orion-hek-harness/RESULTS.md` (why public lines fail: power), `experiments/619*/RESULTS.md`,
`experiments/278-depth-inventory/RESULTS.md`, `experiments/273-flex-reference-instrument/RESULTS.md`,
`data/external/cd4_cell2026/data_sharing_readme.md`.

## Why
Every DE-channel price in this project comes from H1 (h1-test, 618's h1-trainval): one line, Arc's lab, and both surfaces
are worn by dozens of selections. Public cancer lines carry almost no DE power at Arc-like depth (615/616: median 1–2
calls per target). GWCD4i (Zhu/Dann/Marson, genome-scale CRISPRi in primary CD4 T cells) is **10x Flex, the challenge's
own chemistry, ~570 cells per target per condition**, cell-level on a public bucket in us-east-1
(`s3://genome-scale-tcell-perturb-seq/marson2025_data/D{1..4}_{Rest,Stim8hr,Stim48hr}.assigned_guide.h5ad`, ~110–160 GiB
each; read it from EC2 only, `--no-sign-request`, same region so the transfer is free). It is the one candidate for a
second, unworn, Flex DE surface in a third biology. Its CD4 *pseudobulk* rows are already a production input (the CD4 stack
residual), so the surface's bulk input must exclude CD4 (DESIGN §3: the held line is out of every donor mix).

## Stages (each gated; stop and report at a failed gate)
0. **Power gate (decides everything, run first, cheap):** on EC2, one condition (Rest), the four donors' files: take 100
   non-essential, knockdown-significant targets with ≥ 400 cells (random, seed 805), 400 cells each (pooled over donors),
   and 18,400 non-targeting controls (pooled, every donor in proportion), downsample every cell to a 20,000-UMI median on
   the official 18,533-gene axis (Arc's release protocol; record the axis overlap), and run the scorer's DE test (Wilcoxon,
   BH 0.05 per target, genes > 5 CPM in controls; cell-eval2's own DE code or 614's ref614.py route). **GATE: median
   n_real ≥ 20 and ≥ 60 % of targets with n_real ≥ 10** (H1: 89 and 68 %; K562: 5; Orion clean: 1). Also report the
   knockdown (target's own lfc) distribution and the blob test of 615 (split-half share ≤ 0.10). If the gate fails, try
   Stim48hr once (the same rule, written down before), then stop with RESULTS.md: a well-measured negative is a result.
1. **The surface** (only if stage 0 passes): one condition, a **300-target panel** drawn by a frozen rule written before the
   draw: non-essential, knockdown-significant, ≥ 400 cells, basal expression above the A/B/C panel's 5th percentile in CD4
   controls, strength mix matched to A/B/C's using a pilot that is NOT this line (K562-GWPS downstream DE counts, as Arc
   used an H1 pilot), and **excluding every target of every `pert_counts.csv` on disk** (A/B/C, real, efp, efr, R6, R7a,
   the lab's 801 draws). Truth = 400 cells per target + 18,400 controls, Arc-like depth, as stage 0. Build it as an h5ad in
   the challenge's format (`target_gene`, `non-targeting`), plus its cell-eval2 real bundle on EC2 (`cell-eval2 baseline
   --save-pred` → `prep-real-bundle --baseline <pred.h5ad>` → keep both on S3).
2. **The real-emitter harness on it**, 618's way: a reference dir (evaluate323.load_reference / score on this truth), a
   faithful key for this panel (run618.py's stage / operator / compose parameterised by panel, with **CD4 removed from the
   composition** and H1 kept), CLEAN4 banks (never CD4), CLASSES3 from this surface's own control libraries, and a
   `run805e.py` wrapper of run496 in 618e's style. KA first: the wrapper with no swaps on H1 reproduces 435's stats (618's
   ka618e row) on EC2.
3. **The trust gate, written down before any gate run (TRUSTGATE.md, bands like 618's):** the board's known DE pairs:
   zs (446bj2z − 446bj2, board nmae −0.00218), coarse MDE (board nmae −0.01679, reach +0.0024), CLEAN4 hook (496a, board
   nmae −0.00158), UPRN (630arm, nmae −0.00305), DN (nmae −0.001646), each 3 seeds. Report board / surface ratios per pair
   and regime (619's table format). PASS = every nmae pair the right sign and within 0.5×–2× of the board.

## Rules and budget
- Cloud budget: **$80** (stage 0 ≤ $10). Everything heavy on EC2 via `experiments/800-checker/cloud/labcloud.sh`;
  large intermediate files stay on S3 under `output/lab/` (record paths + sha256). Nothing > 2 GB on the Mac.
- Never use this surface's truth to choose anything you then price on it; every rule is written before its data is read.
- Write RESULTS.md as you go; stage 0's result goes to the lab lead the moment it exists (it decides whether this continues).

## Return
RESULTS.md: the power gate (numbers, sha256 of the subset), and if it passes, the surface's files, the KA rows, the trust
gate table, cost, and the driver interface (`run805e.py <mode> <name> <run496 args>` → stats files, members at f 1 / 0.5 /
0.35 / cmp, the same arrays 803 returns).
