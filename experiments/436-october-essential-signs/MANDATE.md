# 436: the October-only DE sign lever. Mandate from the PI, 2026-09-23 10:15Z

**Why.**
- The essential-gene Perturb-seq screens we hold cover 0 of the 300 A/B/C targets but ~30 D/E/F targets each
  (`docs/DEF-RUNBOOK.md` §5.-2). The screens are K562-essential, RPE1, Nadig Jurkat and Nadig HepG2.
- On H1, their donor signs are much better than our delivered signs on the targets they cover.
  - 430 §0.3 found that H1's donor-sign advantage comes entirely from them.
  - 434g (`experiments/434-posterior-median-magnitude/RESULTS.md`, analytic): on H1's 20 covered targets, signing
    426's MDE add by the covering donors' mean buys nmae **−0.024** there, against −0.0075 with our own sign.
- On a 300-target panel with ~45 covered targets that is about **+0.0007 overall**.
- It is an exact no-op on A/B/C, so the board can never price it. H1 is the only instrument, and October is the only
  place it acts.

**Read first:**
- `AGENTS.md` (all of it);
- `docs/FINDINGS.md` §0–§2 (especially the "DE magnitude and delivery" block);
- `experiments/426-mde-magnitude/RESULTS.md` (the MDE and 426h `coarse`; note 426b's prior-sign MDE failure: class
  over-fill, "did not meet the per-cell totals");
- `experiments/434-posterior-median-magnitude/RESULTS.md` §434g and `ess434.py`;
- `experiments/430-cancer-fold-signs/RESULTS.md` §0;
- `docs/DEF-RUNBOOK.md` §5.

**Territory.** Only `experiments/436-october-essential-signs/`.
- No git, no `vcc`, no `cell-eval2`, no edits to `scripts/` or `src/`, no builder runs. The builder belongs to other
  runs today.
- H1 emission runs are allowed. Each is ~7 GB RSS, strictly serial, behind the memory gate (`memory_pressure` ≥ 40 %
  free, free+inactive ≥ 20 GiB).
- A ~20 GB production build may be running at the same time; the gate takes care of it.

**Tasks.**
- **(a) Harness.** Use `experiments/435-zs-on-coarse/run435.py`: 429's six-member H1 harness plus `--mde-coarse`,
  run with `RUN426_EMITTER=prod`, i.e. the production emitter `src/vcc2026/rc_emit_k.py` sha `edb5b8dc…`.
  - Copy it into your directory and add one option: an `mde_sign` spec vector (gated gene axis, ±1).
    - On targets that any essential screen (`k562_essential`, `rpe1`, `nadig_jurkat`; add `nadig_hepg2` if a table
      exists; check `experiments/371-consensus-block-prior-lolo/` and `data/external/nadig2025`) covers, the sign is
      that of the mean delta over ALL covering donor tables. Never use `h1_2025`.
    - Everywhere else, and wherever that mean is 0 or NaN, the sign stays as it is. The production emitter already
      has an `mde_sign` path (lines ~1015–1047); read it before you design anything.
  - Known answers:
    - with the option off, your copy must reproduce 435's `base_s{S}_bps_prodem_mde0.1cap0.1min0.3coarse` stats
      bit-for-bit (see `experiments/435-zs-on-coarse/ka435.py`);
    - with the option on for a target set that is empty, it must also reproduce them.
- **(b) The measurement.** `base S --bulk ps --mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse` with and
  without the donor-sign option, seeds 0/1/2.
  - Report all six raw members:
    - on all gated targets;
    - on the 20 essential-covered gated targets (per-target nmae from the per-target CSV);
    - on the 38 others, which must not move beyond seed noise.
  - Report `mde_lambda` per block and every refusal or capacity failure. 426b's prior-sign attempt crashed on class
    over-fill; `coarse` may or may not avoid that. If it does not, find the smallest fix: e.g. apply the donor sign
    only to non-block genes, or only where the donor mean agrees across ≥ 2 tables.
- **(c) October sizing.**
  - Use 199's synthetic D/E/F bundles (`experiments/199-def-runbook/work/def-bundle*`) and the pessimistic draws in
    DEF-RUNBOOK §5.-2 / 225.
  - Report how many D/E/F targets each essential screen covers, the union, and the covered gene counts.
  - Scale (b)'s per-covered-target gain to the panel.
- **(d) Production design, as text only.**
  - The exact builder change:
    - an env knob, e.g. `VCC_DE_MDE_SIGN=essmean`;
    - where the covering-donor signs are computed;
    - which donor tables and shas it pins;
    - how the config records it;
    - that unset is bit-identical.
  - The test cases that would prove it.
  - Do not apply it.

**Deliverable:** `experiments/436-october-essential-signs/RESULTS.md`:
- the known answers;
- the six-member table (all / covered / uncovered);
- the capacity story;
- the October sizing;
- the production design;
- a direct recommendation: ship it in the October stack or not, and at what expected overall gain on D/E/F.

## RESUME after a machine reboot (PI, 2026-09-23 10:45Z). Read this first.
- The Mac rebooted at ~10:39Z 09-23 and killed your previous run of this mandate mid-way, plus its background drivers
  and chains (`work/chain_q2.sh`, `work/chain_q3.sh`; `drive436-q2.log` stopped during `base 2 … --mde-coarse`).
- Everything you wrote is still here: `run436.py`, `rc_emit_k436.py`, `ka436.py`, `ana436.py`, `signs436.py`,
  `size436.py`, `drive436.sh`, the queues, the logs, `work/h1/` (20 files), `work/ka436.json`, `work/size436.json`,
  `work/tab436-dsall.json`.
- **Resume from those files.** Re-run only runs whose stats/metrics files are missing or truncated; a run is complete
  only if both files exist and its END line says rc=0. Re-verify your known answers from the stored files before you
  read any number.
- Production builder `2c90d9e1…` and emitter `edb5b8dc…` are unchanged.
- The machine may reboot again. Keep your state on disk as you go.
