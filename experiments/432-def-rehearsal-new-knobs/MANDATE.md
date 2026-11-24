# 432: D/E/F dry rehearsal of the September knobs on the synthetic bundle. Mandate from the PI, 2026-09-23 07:30Z

**Why now.**
- `docs/DEF-RUNBOOK.md` §5.-1 says the full rehearsal must be re-run once the stack is final, by Oct 15. The stack is
  not final.
- But every knob added since 410 has only ever run on A/B/C. On Oct 22, any integration bug in these knobs costs hours
  of a 14-day dark window:
  - ballast unwelding, the per-target key;
  - bulk sharpening, bulk dump/load (the stack), per-context re-level (`ntclevel417.py`);
  - the MDE (`VCC_DE_MDE`), the zero repair (`VCC_DE_ZERO_REPAIR`), and soon 429's zs (`VCC_DE_MAG`).
- Find those bugs now, on the synthetic D/E/F bundle.

**Read first:**
- `docs/DEF-RUNBOOK.md` (all of it, including the 410/417/426/428 additions at the end);
- `docs/FINDINGS.md` §0–§1;
- `experiments/199-def-runbook/RESULTS.md` (the rehearsal machinery and the synthetic bundles in
  `experiments/199-def-runbook/work/def-bundle*`);
- `scripts/def_pipeline.py` (its `SHIPPED_ENV` is stale: it still ships early September's h1union line).
- The current best configuration is in `submissions/428h-*.h5ad.meta.json` (and `428a`, `426a`, `417e`).

**Territory.** Only `experiments/432-def-rehearsal-new-knobs/` and scratch directories you create inside it.
- **Never edit `scripts/`, `src/`, `data/challenge/` or any other experiment.** Never write to `submissions/`. No
  `vcc submit`, no git.
- The driver's pin stage temporarily edits the production builder and restores it. That is allowed, but you must be
  the ONLY builder running. Check `pgrep -f make_submission_pm_cd4` before every build; if another builder runs,
  wait.
- Memory gate before every build (`memory_pressure` ≥ 40 % free, free+inactive ≥ 20 GiB). Builds are ~20 GB RSS and
  strictly serial.
- Leave the production builder byte-identical to what you found (sha before and after, logged).

**Tasks.**
- **(a) Base path.** Run `def_pipeline.py` on the realistic synthetic bundle (`def-bundle-r` or the one 199 names as
  realistic) with the CURRENT best config exported as env:
  - 428h's config minus the A/B/C-specific values;
  - the 410 key rebuilt for the new panel per DEF-RUNBOOK's 410 addition. If that needs a builder edit, STOP and
    document the exact edit;
  - no LOAD and no CTX scale in this first pass.

  Record every stage's timing, peak RSS and outcome. Keep the artifact only until (b)–(d) are done.
- **(b) Stack dump/load on D/E/F.** Reproduce 417e's stack (`VCC_BULK_DUMP` of the H1 ×2 and CD4 ×1.5 compositions →
  compose → `VCC_BULK_LOAD`) on the synthetic bundle, following `experiments/417-bulk-sharpen/run_dump417.sh` /
  `chain_loadarms417.sh`. Show that LOAD refuses an A/B/C dump (it must).
- **(c) Re-level.** Run `experiments/417-bulk-sharpen/ntclevel417.py <artifact> --bundle <synthetic> --contexts D,E,F
  --sd 0.01` and a rebuild with the printed `VCC_BULK_CTX_SCALE`. Check with `cliffcheck414.py` (or its D/E/F
  equivalent) that the energies move as intended.
- **(d) DE knobs.** Confirm `VCC_DE_MDE=0.1,0.1,0.3` and `VCC_DE_ZERO_REPAIR=-5` run on D/E/F:
  - `mde_lambda` distribution per context, via the meta's rc diagnostics;
  - `experiments/428-zero-gene-audit/audit428.py <artifact> --bundle <synthetic>` must report 0 zeroed gated genes;
  - if the zs patch (`VCC_DE_MAG`) is in the production builder by then, include it.
- **(e) Deliverable: `experiments/432-def-rehearsal-new-knobs/RESULTS.md`:**
  - a stage/knob table (ran / failed / time / RSS / what to look for);
  - every bug or trap found, with the exact fix, as text: you do not apply fixes to production files;
  - the exact day-of command sequence for this stack, ready to paste into DEF-RUNBOOK;
  - the new `SHIPPED_ENV` / flags line `def_pipeline.py` should carry.
  - Delete the large artifacts (h5ad, vcc) at the end; keep logs and metas.

## Update from the PI, 2026-09-23 09:50Z (read this; it overrides the text above where they differ)

**State of production.**
- Builder `scripts/make_submission_pm_cd4.py` is sha `2c90d9e1…` (commit `5568ee84`); emitter `src/vcc2026/rc_emit_k.py`
  is `edb5b8dc…`. Assert both at start and at end.
- DE knobs now in production:
  - `VCC_DE_MDE="c,cap,min_count[,coarse]"`;
  - `VCC_DE_ZERO_REPAIR=<floor>`;
  - `VCC_DE_MAG=dn_k250_zs` (429's zs pass);
  - `VCC_RC_AUDIT_DUMP=<path.jsonl>` (431: one JSON line per block, diagnostics only, no config key).
- **Rehearse the `coarse` form: `VCC_DE_MDE=0.1,0.1,0.3,coarse`.**
  - Without `coarse`, context A got zero MDE dose on the real data (431).
  - With it, A got 0.956 (433a's build). Whether a D/E/F context behaves like A is exactly what the dump tells you.
- The reference config is 433h, a 09-25 candidate: 417e's stack LOAD, coarse MDE, zero repair -5, and 417h's re-level.
  It is built by `experiments/433-coarse-mde-arm/chain433hk.sh`; its meta will be at
  `submissions/433h-*.h5ad.meta.json`. The re-level values are A/B/C-specific: derive D/E/F's with `ntclevel417.py`,
  as task (c) says.

**Builder schedule.**
- The PI's chain (`chain433hk.sh`: 433h then 433k) owns the builder until its log
  `experiments/433-coarse-mde-arm/work/chain433hk.log` contains a line ending in `END`, or until no process matches
  `pgrep -f chain433hk.sh`. Expect that around 12:15Z.
- Do not start any build before then. Use the time to read, plan, and write your command sequences.
- After that the builder is yours until you finish. Keep builds strictly serial, behind the memory gate.

**Overrides of the launcher's generic guardrails, for this mandate only.**
- You may run `scripts/def_pipeline.py` builds on the SYNTHETIC D/E/F bundles. A build peaks at ~20 GB RSS; keep them
  strictly serial and behind the memory gate. This overrides the launcher's generic ~10 GB note for these builds
  only; everything else you run stays under ~10 GB.
- You may run `vcc prep` on synthetic artifacts, with output inside your experiment dir, if a stage needs it.
- Everything else stands:
  - never write to `submissions/`;
  - never run `vcc submit` or `cell-eval2`;
  - no git;
  - leave the production builder byte-identical: log its sha before and after every build.

**Add to the deliverable.** Give the per-context `mde_lambda` distribution for D/E/F from the dump, and the
`zero_repaired` / `zs_n_moved` medians.

## RESUME after a machine reboot (PI, 2026-09-23 10:45Z). Read this first.
- The Mac rebooted at ~10:39Z 09-23 and killed your previous run of this mandate mid-way, plus its background chain
  (`work/chain432.sh`, which was waiting for chain433hk).
- Everything you wrote is still here: `build432.sh`, `chain432.sh`, `key432.py`, `cmpdump432.py`, `proposals/`, and
  `work/` (`state432.log`, the `ef_key` phase done at 10:12Z, `dryrun/`, `wraptest/`). **Resume from those files:**
  re-verify what they record instead of redoing finished phases, and do not trust a phase whose log has no "done" line.
- **The builder schedule changed.** 433h finished (e7495644). 433k is being rebuilt by
  `experiments/433-coarse-mde-arm/chain433k_resume.sh` (started 10:43Z, ~35 min).
  - Its completion signal is the same: `experiments/433-coarse-mde-arm/work/chain433hk.log` gets a line ending in
    `END` after 10:43Z.
  - There is no longer a `chain433hk.sh` process, so do not use `pgrep -f chain433hk.sh` as the signal. Use the log
    line, AND check that no process matches `pgrep -f chain433k_resume.sh` and no `make_submission_pm_cd4` python is
    running.
- Production builder `2c90d9e1…` and emitter `edb5b8dc…` are unchanged.
- The machine may reboot again. Keep writing your state to `work/state432.log` after every phase, as you were.

## RESUME 2 (PI, 2026-09-23 11:28Z): the builder is FREE now. Start building. Read this before anything else.
- The PI's 433k rebuild finished at 11:21:08Z. Its chain log ends with the END line, and no PI builder or chain is
  running. **Do not wait for any PI signal any more.** Your only gate before each build is the standard one: no other
  `make_submission_pm_cd4` python running, plus the memory gate.
- **Why your chain never started (a trap, now documented in AGENTS.md):**
  - `pgrep -f <name>` matches ANY process whose argument list contains that string.
  - Your own `claude -p` process has this whole mandate as its argument, and the RESUME note above names the PI's
    resume script. So `pgrep -f` on that name always found your own process, and your chain would have waited until
    its 20:00Z give-up.
  - PI killed that chain (pid 5501) and your previous run at 11:26Z; the log of that run is saved as
    `claude-run.run2-killed-1126Z.jsonl`.
  - Rule for your scripts: never `pgrep -f` a name that could appear in any agent's prompt. Match the interpreter
    and the script path instead, e.g. `ps -axo comm,args | awk '$1 ~ /^(bash|zsh|sh)$/ && /<script>/'`, and check
    `comm`, not only `args`.
- Your files and phase logs are intact (`work/state432.log`: `ef_key` done at 10:12Z; nothing else ran). Resume from
  them as before.
- After you finish, a PI chain waits for your claude process to exit and the builder to idle 10 min, then builds two
  arms. Do not wait for it or touch it.
