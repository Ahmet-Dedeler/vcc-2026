# 803 — the lab's DE engine: the production emitter on the H1 surfaces, on EC2 (worker mandate)

Written 2026-10-15 ~22:47Z by the lab lead ([c61d4c]). You own `experiments/803-de-engine/` and write only there.
Read first: `docs/LAB-2026-10-16.md`, `experiments/800-checker/DESIGN.md` (your work is the DE half of §3–§4), then
`experiments/618-h1-trainval-surface/RESULTS.md` + `TRUSTGATE.md`, `experiments/619-board-multipliers-618/RESULTS.md`,
`experiments/496-uncalled-sign-hook/RESULTS.md`, `experiments/630-rank-neutral-up/RESULTS.md`, `experiments/631*/`,
`experiments/654-carrier-amplitude/RESULTS.md` (zero repair), `experiments/521-strength-nreal/lib521.py`.

## Objective
The checker must price DE-channel candidates (emitter / magnitude / sign / call-set changes) and the DE side effects of
bulk candidates (479sp2 lost on the board through nmae via ballast; no bulk-only fold can see that). Today that is done
by hand, on the Mac, chain by chain (run496 → run618e → score618 / score619). Turn it into an engine that runs on EC2:

1. **One driver** `de803.py run --surface {h1test,h1tv} --arm <arm.json> --seed S` that wraps run496 (h1-test) and run618e
   (h1-trainval) verbatim (import, never edit; redirect WORK into 803) and writes, per run,
   `work/runs/<surface>/<arm>/s<S>/members.npz`: PER-TARGET fid, jac, reach, nmae at f 1, f 0.5, f 0.35 and cmp (521/528's
   `rows_at` / lib521 member functions, the realized order, exactly as score619 computes them), the per-target ballast K,
   the per-target pseudobulk sums, plus a manifest (argv, git blob of every imported module, input sha256s, wall, peak RSS).
   If pds and mse can be computed from the sums against the surface's reference with the scorer's formulas (DESIGN §4;
   cell-eval2's own functions through evaluate323 if they expose them), add per-target pds and the mse numerator /
   denominator terms; if not, say why.
2. **Arms are files** (`arms/<name>.json`): the run496/run618e mode + flags + env (zero repair −5, production emitter,
   coarse MDE, zs, hook, UPRN, DNMAG, K / call-set CPM, bulk dump …). Write the arms for:
   - `inc`: the October DE incumbent on each surface (zs + coarse MDE + CLEAN4 hook 0.05 + UPRN + DNMAG + zero repair −5;
     where a component has no wrapper on a surface yet, say so and stop at the largest incumbent the wrappers support);
   - the known-answer pairs: zs (with / without, = 618's G1), coarse MDE (G2), CLEAN4 hook (c4h vs none), UPRN (630 vs
     none), DN (630+DN vs 630), idall3 vs c4h; and, if the wrappers can express them, call-set size K 2,000 / 1,600 vs
     2,200 (`VCC_RC_CALLSET_CPM`; board: both HURT via fid).
3. **Known answers first (abort on failure):** on EC2, (a) `run618e h1ka` == 435's H1 stats bit for bit (618's ka618e
   check), (b) 618's g1 and hk_c4h s0 stats files bit for bit. Linux vs macOS float or numba differences are possible:
   if a KA is not bit-equal, quantify it (max |Δ| on p, q, lfc; member Δ) and continue only if the member deltas between
   arms reproduce 618's to < 2e-5 raw. (c) score619's f 1 column equals 618's gate on every member of every pair.
4. **Run the known-answer pairs on both surfaces, 3 seeds each**, on EC2 (one r7i.8xlarge running ~8 emitter runs at a
   time; a run takes ~7.5 min on the Mac). Bring back only the members.npz + manifests (small).
5. **A real-scorer cross-check** (one run): if the emitted cells of one base run can be written as an h5ad without changing
   the emitter (a hook around the block writer), score it with the real `cell-eval2` (0.16.0, `--preset vcc2026`) on EC2
   against the surface's real h5ad and compare the four DE members with the engine's f 1 numbers. If writing cells needs
   more than a wrapper, skip it and say so.

## Rules and budget
- Cloud budget for this mandate: **$60**. No cell-eval2, emitter run or anything > 2 GB on the Mac.
- `data/donors/h1_2025` never enters an H1 run (oracle leakage).
- Do not aggregate into board units or decide bins; return per-target arrays. The lab lead's check800 does the rest.
- Aim to have the KA table and the known-answer pairs by ~10-16 20:00Z. Write RESULTS.md as you go.

## Return
RESULTS.md: the KA table, the arms (files + what each changes), per-pair mean member deltas per surface and regime with
seed SDs (as a sanity table, next to 618/619's), the cloud runs and cost, the driver's interface, and what is not
expressible yet.
