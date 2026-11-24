# 555 — adversarial review of the R7 machinery

Reviewed the blobs on `origin/master` at `21381a77` (the R7a driver, the aggregation, and the claim draft). Those three blobs are unchanged through `9725b988`. This checkout is `cloud/555-r7-review` at `54b8886d`, which does not contain the new files; the review ran against `repro/upstream/` copies whose `git hash-object` matches `origin/master`:

| file | blob |
|---|---|
| `experiments/551-r6-c10b/r7a551.py` | `19246a8dad622a24faa7c060aa3fae5e318d4b34` |
| `experiments/551-r6-c10b/aggregate_r7a551.py` | `2867b22b2602abc1ff5e1ac916bf25f10f53ccdb` |
| `experiments/552-r7-reserve/CLAIM-R7-DRAFT.md` | `757e502c15b2eea0b1b7370d808e9a42041329a8` |

`gate551.py`, `lib551.py`, `price551.py`, `designs551.py`, `full551.py`, `comp551.py`, `drive551.py`, `drive534.py`, `axes534.py`, `access553.py`, `band553.py`, `lib553.py` are the same blob on this HEAD and on `origin/master`. No `vcc`. No truth h5ad was opened. Patches under `patches/` were compiled from `/tmp` copies and were not applied.

**RAN** means a command in this review produced the log named below. **READ** means the conclusion is from the source (or from a tracked json) and was not executed as that program.

---

## Findings

### BLOCKER

**1. The instrument pin does not cover `ntclevel417`, and every T panel calls it.**

`experiments/551-r6-c10b/lib551.py:642` imports `ntclevel417` inside `ntc_claims`. `experiments/551-r6-c10b/full551.py:224-226` calls `ntc_claims` for every panel that does not start with `abc`. T1–T4 take that branch. The levelled claim is `cp46 * wsum` from that module (`lib551.py:648-649`).

`r7a551.py:89-106` finds the pin by importing `lib551`, `designs551`, `price551`, `gate551` and by exec'ing only the text of `full551` / `comp551` / `aggregate551` above the first `\ndef `. It never calls `ntc_claims`, so the import does not run. The committed plan (`experiments/551-r6-c10b/work/r7a/plan.json` on `origin/master`) lists 28 blobs, `differs_from_ka` empty, and does not name `ntclevel417`. `access()` (`r7a551.py:376-377`) only re-checks those blobs plus `external`.

Failure: after `record`, edit `experiments/417-bulk-sharpen/ntclevel417.py` (today's blob `3270a7cc1e3acdcf575a8a883f20325afdad7bb3` equals `e7be21ca`). `plan` / `access` still pass. T-panel level factors and every levelled price move, including PRIMARY. The efp known answer no longer describes the code that ran.

Smallest fix: in `instrument_files`, add that path to the hashed list so `blob_at(e7be21ca)` covers it. Do not edit `lib551.py` to import it at module level: that would change a file the pin requires to equal `e7be21ca`. Unapplied: `patches/01-r7a551-pin-floor-unlink-ledger.diff`.

**2. `access553.py --owner-go` does not look at the ledger. The draft's sha is not in the ledger, and the access would still arm.**

`access553.py:4-5` says the sha256 is in the ledger's R7 section. The code does not read `docs/CONFIRMATION-LEDGER.md`.

- `access553.py:532-535`: sealed and not `--dry-run` requires the path to be an existing file, then stores `sha256(file)` in `plan.json`.
- `access553.py:142-145`: each phase checks the file's sha256 still equals that plan entry, then `lib553.arm_access` (`lib553.py:243-244`) sets `_ACCESS_ARMED[name] = True`.
- `lib553.py:235-240`: `open_truth` only checks that flag.

`r7bd5000` copies the `r7b` spec and stays `sealed=True` (`lib553.py:174-189`). A control half `r7bd5000c<k><h>` copies that spec and does not clear `sealed` (`lib553.py:116-133`).

RAN `repro/test_owner_go.py` (`repro/logs/owner_go.txt`): `surface_spec("r7bd5000")["sealed"]` is true; `arm_access` is the one assignment above; the draft snapshot's sha256 `b20810cf7cfdd193c3a4dff66f44bc4da1c7220b3c4b42b9df3dbe3c26dfb9a9` is not in `docs/CONFIRMATION-LEDGER.md`; the flag was set anyway.

Failure: `access553.py r7bd5000 --owner-go experiments/552-r7-reserve/CLAIM-R7-DRAFT.md ...` (the draft's own first line says it is not the claim) opens the floored truth. Phases only prove the file did not change during the run.

Smallest fix: before `plan.json` is written and before phase `ctrl`, require that sha256 to occur in the ledger text under `## Reserve R7` and not merely somewhere in the file. Unapplied: `patches/03-access553-ledger-and-half-coverage.diff`.

**3. The published one-command access scores switch pairs that have no `band_total`, and `band553` hides that.**

`band553.py:50-52`: if an effect key from the main access is missing on any half, `continue`. No error, no row.

`access553.py:506` defaults `--half-arms` to `lcap-kept,cpm+c4h,lcap-kept+c4h,cpm-n1600,lcap-kept-n1600+c4h` and `--half-pairs` to the hook-on-CP pair and the K-on-CPH pair. Half readiness (`access553.py:591-596`) and the half subprocesses (`access553.py:619-620`) use those defaults, not `--arms` / `--pairs`.

`experiments/553-r7b-harness/RESULTS.md:510-518` runs exactly the claim's flags (`--ctrl-splits 5 --retro work/access553_rehearsal5 --halves 0,1,2,3,4`) with

```
ARMS=lcap-kept,cpm+c4h,lcap-kept+c4h,cpm-n1600,lcap-kept-n1600+c4h,cpm-n1200,lcap-kept+c4h-s
PAIRS=lcap-kept+c4h:lcap-kept,lcap-kept-n1600+c4h:lcap-kept+c4h,lcap-kept+c4h-s:lcap-kept+c4h,cpm-n1200:cpm
```

and does not pass `--half-arms` or `--half-pairs`. The claim draft (`CLAIM-R7-DRAFT.md:76-77`) is the same shape with `<FILL>` in `--arms` and `--pairs` only. The draft also says to emit every decided arm at seed 0 on the ten halves (`CLAIM-R7-DRAFT.md:56-58`). Emitting them does nothing unless those names are in `--half-arms` / `--half-pairs`: the half accesses and `band553` never load them.

Failure: main access writes `cpm-n1200 vs cpm` and `lcap-kept+c4h-s vs lcap-kept+c4h`. Readiness passes because the six core arms exist on the halves. `band553` drops the two extra effects. The process exits 0. The HARM rule ("≤ −2.5 × band_total of that pair") has no band, and a later reader can invent one.

RAN `repro/test_band_skip.py` (`repro/logs/band_skip.txt`): base keys `lcap-kept vs cpm` and `cpm-n1200 vs cpm`, halves contain only the first, the second is dropped, no exception.

`--pairs` itself is parsed only in `report()` (`access553.py:398-401`), which is the last phase, after `ctrl` and `ref` have opened the truth. A malformed pair spends the reserve and then raises. `RESULTS.md:527` also says a finished access can be re-reported with new `--pairs` / `--retro` without reopening the truth. That is a contrast chosen after the number exists.

Smallest fix: before phase `ctrl`, require every `--arms` arm and every `--pairs` side to be in `--half-arms`, and require every `--pairs` entry to be in `--half-pairs`. In `band553`, raise if a main-access effect is missing from a half. Unapplied: patches 02 and 03.

### SHOULD-FIX

**4. `arm_gate` does not patch `floor_gib`. A stale `GATE-OVERRIDE.json` replaces the 22 GiB floor.**

`gate551.py:166-177`: `floor_gib` returns the override's `min_free_gib` when `experiments/551-r6-c10b/GATE-OVERRIDE.json` exists, and only otherwise calls `thresholds()`. `gate()` (`gate551.py:185-186`) uses `floor_gib()` for the free+inactive floor and `thresholds()` for pressure and swap. `r7a551.py:279-281` replaces `window_ok`, `thresholds`, `GATE_LOG`, `HEAVY_LOG`. It does not replace `floor_gib` or `OVERRIDE`.

All of those names are loaded with `LOAD_GLOBAL` at call time (RAN `repro/test_gate_patch.py`, `repro/logs/gate_patch.txt`). Assigning `G.window_ok = ...` therefore takes effect for `gate` / `wait_gate` / `heavy`. `floor_gib` keeps calling the original `OVERRIDE` path.

RAN, same log: with the patch applied and an override file `min_free_gib: 10`, `floor_gib()` returned `(10.0, "owner override 10 GiB by stale")` while `thresholds()` returned 22 GiB. Without the override, `floor_gib()` returned 22.

Failure: a C10b override left on disk (legal range 10–27 GiB, `gate551.py:174`) is the floor R7a's heavy steps actually wait on.

Smallest fix: in `arm_gate`, replace `floor_gib` with a lambda that returns 22 GiB and does not read the override file. Same patch 01.

**5. A crash during operator deletion is not resumable. The fits stay "done" and are not rebuilt.**

`r7a551.py:361-366` (the same loop as `drive551.py:145-150` and `drive534.py:148-153`): for each npz, stat it, unlink it, and only then `done(ops_deleted)`. `heavy()` (`r7a551.py:284-287`) returns immediately when the step is already in `state.json`, before `G.heavy`.

RAN `repro/test_resume_and_seal.py` (`repro/logs/resume_and_seal.txt`): four fit steps marked done, sidecar json present, npz absent, comp jsons carrying `peak_rss_gib` and `op_axis_checked`. `panel()` skipped `full551` and `comp551`, then raised `FileNotFoundError` on `T1__held_k562_gwps.npz`. `ops_deleted` stayed absent. A whitener failure (nonzero rc) is different and is safe: `heavy()` does not call `done` (`r7a551.py:291-293`), and the whitener overwrites its output, so a re-run re-fits. READ of `scripts/build_union_whitener.py` save path, not executed.

Failure: kill the process after the first `unlink` and before `done`. The next `access` cannot rebuild the operator (the fit step is done) and cannot finish deletion. `aggregate` never starts.

Smallest fix: if the npz is already gone, record the sha256 from the sidecar json and continue. Same patch 01.

**6. `record` / `access` treat any occurrence of the hex string as the ledger entry, and `--config` is not read from the claim.**

`r7a551.py:246` and `r7a551.py:374`: `h in LEDGER.read_text()`. `r7a551.py:243` checks `a.config in designs551.CONFIGS` and never opens the claim markdown for the key. `record` (`r7a551.py:242`) refuses a second call only while `work/r7a/claim.json` exists. Deleting that file allows another record.

RAN, same resume log: a ledger file whose only mention was the sentence "Do not arm anything with `<hash>`…" was accepted; the claim text said M5; the recorded config was `M6+P2g`. A second `record` raised `AssertionError` while `claim.json` existed. An absent hash was refused.

Failure: a note that quotes the hash, anywhere in `CONFIRMATION-LEDGER.md`, arms R7a. The recorded F can disagree with the markdown the hash names. `aggregate_r7a551.py` then judges that F.

Smallest fix: take the hash only from the `## Reserve R7` section (stop at the next `\n## ` that is not `###`; `\n## ` does not match `\n### `). Require the claim to contain exactly one line `F = <key>` and require `--config` to equal it. Same patch 01. The frozen claim has to grow that line; the draft's `**F = FILL**` does not satisfy it, which is what you want before the fill.

**7. `aggregate_r7a551.py` has no run-once guard. A later invocation can change the verdict.**

`aggregate_r7a551.py:71-75`: `--panels` defaults to `T1,T2,T3,T4`, `--nb` defaults to 4000, `--test-config` is an alternate F allowed only with `--panels efp`. There is no check that `results_r7a.json` already exists. `verdict` (`aggregate_r7a551.py:66-67`) is TRANSFERS only when the point is positive and the 95% lower bound is positive. `access()` (`r7a551.py:382-388`) skips aggregation only when `state.json` contains `"aggregate"`, and it passes no flags. A direct re-run is outside that skip.

Failure: after the preregistered 4000-draw, seed-552, four-panel run, `aggregate_r7a551.py --nb 200` (or `--panels T1`) overwrites `results_r7a.json`. A bound that sat on zero moves between FRAGILE and TRANSFERS. `--test-config` is fenced to efp and writes `work/r7a/test/` only (`aggregate_r7a551.py:77-79`); that path is fine.

`retention` does not enter `verdict`. It can still change under you: `aggregate_r7a551.py:143-144` fills the A/B/C step from `Panel("abc")` when `work/abc/full551.json` and the three comp jsons exist. A later abc run changes the reported retention, not the verdict word.

Smallest fix: if `results_r7a.json` exists and `--test-config` is unset, refuse to write. Ignore `--nb` and `--panels` unless `--test-config` is set.

**8. Do not freeze the HARM line as written. −2.5 × `band_total` from zero is the wrong test for K-on-CPH, and the "1.3×" citation is a different ratio.**

RAN `repro/calibrate_harm.py` (`repro/logs/calibrate_harm.txt`) on `work/band553_h1r7.json`, `work/draws553.json`, and the five `access553_h1r7*` jsons. `band553_h1r7.json` has `between_draw_sd: null` because that run did not pass `--draws`. The between-draw sds below are from `draws553.json`. Normal tail `Φ(−2.5) = 0.0062`.

`band_total` in `band553.py:61-64` is `sqrt(band_sd² + sigma_c²)` with `sigma_c = sqrt(max(mean(d²) − 2 seed_sd², 0)) / 2`. That matches the claim (`CLAIM-R7-DRAFT.md:60-61`) and RESULTS "the band". It does not use `band_sd_ctrl`. `--ctrl-splits 5` writes `band_sd_ctrl` as well (`access553` report header). Substituting it after seeing which threshold fires is a post-hoc choice; the frozen rule has to name `band_total` and stop.

The claim's "between-draw spread … ~1.3× band_total (draws553: 0.83–3.0, 90%)" is CP's `ratio_ctrl` against the cells+controls band (`draws553.json`: `ratio_ctrl` 1.285, and RESULTS.md:267 prints 1.28 (0.83–3.05)). It is not between-draw sd / `band_total`. Those ratios on h1r7's `band_total` are CP 1.57, hook 0.42, CPH 1.44, K2000 0.90, CPHK20-vs-cpm 2.24, hook-on-CP 0.78, K-on-CPH 2.52.

Per effect, T = −2.5 × h1r7 `band_total`. Empirical vetoes are counts of the five draw values. `P0_b` is the null tail if the noise equals `band_total`. `P0_sd` / `Pmean` use the between-draw sd.

| effect | band_total | T | mean | between sd | emp | P0 if noise=band | P0 if noise=between | P(veto at the mean) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| CP | 0.00120 | −0.00299 | +0.00143 | 0.00188 | 0/5 | 0.006 | 0.055 | 0.009 |
| hook | 0.00089 | −0.00223 | +0.00002 | 0.00037 | 0/5 | 0.006 | ~0 | ~0 |
| CPH | 0.00128 | −0.00321 | +0.00152 | 0.00185 | 0/5 | 0.006 | 0.041 | 0.005 |
| K 2,000 | 0.00177 | −0.00442 | −0.00155 | 0.00158 | 0/5 | 0.006 | 0.003 | 0.035 |
| CPHK20 vs cpm | 0.00136 | −0.00340 | −0.00013 | 0.00305 | 0/5 | 0.006 | 0.133 | 0.142 |
| hook on CP | 0.00056 | −0.00139 | +0.00009 | 0.00043 | 0/5 | 0.006 | 0.001 | ~0 |
| K on CPH | 0.00092 | −0.00230 | −0.00165 | 0.00232 | 2/5 | 0.006 | 0.161 | 0.391 |

Family-wise, independent, noise = between-draw sd, true effect 0: CP+hook+CPH = 0.095 (empirical 0/5); hook + hook-on-CP = 0.001 (0/5); K2000 + K-on-CPH = 0.163 (empirical 2/5, both from the K-on-CPH draws).

Power against a K-size harm, μ = −0.0017, sd = 0.0023, T = −0.00230: **0.398**. At the rehearsal mean −0.00165 the same probability is 0.390. The h1r7 point itself is −0.00468, which is −5.1 × its own `band_total` and would already be a HARM. The sentence "a HARM there needs Calu-3 to show it more strongly than H1 did" (`CLAIM-R7-DRAFT.md:73-74`) is false for this rehearsal: two of five H1-shaped draws cross the line, and the expected value sits at −1.8 × `band_total`.

`band_total` is still the right *construct* for one frozen surface. RESULTS.md:307-310 is right that the between-draw sd mixes a reassignment of cells-per-target that R7b will not repeat, and that the all-controls ablation is still 1.25× for K-on-CPH. What is wrong is using that construct as a zero-centered 2.5-sd test for a switch whose preregistered prior is already −0.0017. Under noise = `band_total` and a true zero, the false-veto rate is the claimed ~0.6% (the draft's "1–3%" is that tail, not the between-draw tail). Under the between-draw sd it is 16% for K-on-CPH alone and about 10% for the three call-set switches together.

Proposed rule, to be written into the claim before the sha is recorded:

- CP, hook, CPH, hook-on-CP: veto iff effect ≤ −2.5 × that pair's own `band_total`. State the null tail as 0.6% only under noise = `band_total`. Also write the between-draw tail next to it (CP 5.5%, the three together ~9.5% if independent) so nobody quotes 0.6% as the family rate.
- K-on-CPH: freeze the prior −0.00165 from `draws553.json`. Veto iff effect ≤ prior − 2.5 × `band_total` = **−0.00395**. RAN: power against μ = −0.0017 with sd 0.0023 is 0.164, and the false-veto rate if the truth equals the prior and the noise is 0.0023 is 0.159. A zero-centered threshold cannot be reconciled with "a null or small effect never counts against a switch."
- If that power is not acceptable, K-on-CPH is not a veto candidate on one surface. Report the effect and `band_total`. A −2.5 × between-draw-sd threshold around the prior is −0.00740 and has power 0.007 against μ = −0.0017. One Calu-3 reading cannot resolve a K-sized harm; say so in the claim instead of implying the 2.5 line will.

**9. Strict `raw fid < 0` calls a rehearsal draw "not board-like" once in five.**

RAN, same calibration log. Board-like := raw jac > 0 and raw fid < 0 on `lcap-kept` vs `cpm`, which is the sign of `effect.mean`, not the retro ratio (the ratio is positive when both sides are negative; `access553.py:423-426`).

| draw | raw jac | raw fid | board-like |
|---|---:|---:|---|
| h1r7 | +0.004021 | −0.002852 | yes |
| h1r7s2 | +0.004160 | +0.000062 | no |
| h1r7s3 | +0.004353 | −0.001871 | yes |
| h1r7s4 | +0.003697 | −0.001553 | yes |
| h1r7s5 | +0.004350 | −0.002673 | yes |

4/5. A second command, same five files: jac mean +0.004116 (sd 0.000273, min +0.003697), fid mean −0.001777 (sd 0.001162, max +0.000062). The five-draw file `access553_rehearsal5` stores those means (jac +0.004116, fid −0.001777) and has no `retrodiction` key; the ratio is computed at R7b report time against this file.

Failure: R7b's CP fid lands at +6×10⁻⁵, the claim voids every veto (`CLAIM-R7-DRAFT.md:64-65`), and the surface was the same shape as the rehearsal. jac is not the fragile side.

Proposed: board-like iff raw jac > 0 and raw fid < +0.001 (one rehearsal fid sd, 0.00116). All five rehearsal draws pass. A fid of +0.002, the size of the rehearsal mean, still fails. Apply it to `here`, never to `ratio`.

**10. Retention vs A/B/C mixes the levelled instrument with hand-entered numbers. It does not move the verdict.**

`aggregate_r7a551.py:140-144`, READ, and the jsons were opened:

- R6 side: `results534B.json` `steps[s]["pooled"]["mean"]` for the six milestones. Same levelled overall as `price551`. `results551.json` does not exist yet. When it does, `config` is the whole c10b-config dict (`aggregate551.py:151`) and `claim["pooled"]["mean"]` is that step (`aggregate551.py:121`). `aggregate_r7a551.py:136-139` copies it only if `config.config == F`. If C10b's F differs, the per-target id is missing and `retention()` (`aggregate_r7a551.py:146-148`) returns None for any sum that includes it. Fail closed.
- A/B/C side: `stored[s]["abc"]["mean"]` when the key exists, else `text[s]["abc"]`. `stored` has M0→M1, M1→M2, M2→M3 only (levelled abc means, 0.01178 / 0.00238 / 0.00856). M3→M4, M4→M5, M5→M6 come from `text`. M5→M6's `src` says "522's unlevelled route, retention-scaled" and the value is 0.0021.

Failure: F = M6, PRIMARY is only M5→M6, `retention_vs_ABC` divides a levelled R7a mean by 0.0021 from a different route. The ratio is not a retention. The verdict does not read it (`aggregate_r7a551.py:122`). The draft lists it under "Reported, not judged." Leave it that way, and label the M3→M6 denominators as the record's text field, not as this instrument.

**11. "If October keeps plain tercile, price the carried-only form and say so" is not what the driver does.**

`CLAIM-R7-DRAFT.md:22-23`. `designs551.CONFIGS` (`designs551.py:98-102`) has six keys. P2g is the carried-only form. There is no key for plain tercile, for 459dr's H1×3, or for the 522 strength rider. The ledger's C10b amendment sends later bulk switches to R7; this aggregation cannot express them.

`r7a551.py:138-145` aborts if any of the 28 files differs from `e7be21ca`. A 598 edit that touches a pinned file makes `plan` / `access` abort. It does not "price the carried-only form and say so." 598 was not applied on the tree this review read; the conflict is with the sentence, and it becomes real the moment a pinned file moves.

Smallest fix: the frozen claim says the instrument stays at `e7be21ca` and F is one of the six keys. Anything else is a new pin and a new efp known answer before `record`, not a sentence after the score.

**12. F = M5 feeds the cumulative M0→M5 to the switch verdict.**

`aggregate_r7a551.py:84-86`: if `post` is empty, `primary_ids = MILESTONES[:5]` (M0→M1 … M4→M5). `verdict` is unchanged. C10a's cumulative M0→M6 was +0.0229 and M5→M6 was +0.00694 (ledger outcome C10a), so M0→M5 is about +0.016 and will be called TRANSFERS, whose consequence is "nothing changes" (`aggregate_r7a551.py:160-161`). The action happens to be right (there is no post-C10a switch) and the word is wrong: TRANSFERS is defined as the post-C10a switches transferring (`CLAIM-R7-DRAFT.md:33-36`).

Smallest fix: if F = M5, set the primary point to 0, verdict to a preregistered token meaning "no post-C10a switch", and keep M0→M5 as the secondary cumulative only.

### NOTE

**13. `mech_a12_audit` is a lazy import and is not on the R7a path.** `mse187.py:86-87`, called from `mse187.main` (`mse187.py:182`). `joint187.py:27` imports `CTX, BOARD_RAW, score, solve_X, mse_ctx` only. RAN the AST scan (`repro/logs/lazy_scan.txt`): uncovered repo modules are `ntclevel417` (BLOCKER 1) and this one.

**14. Non-repo lazy imports the pin cannot see:** `h5py` (`lib551.py:181`, `plolo.py:289`), `cell_eval2.metrics.delta._row_weights` (`lib551.py:643`, on the T-panel claim), `scipy.stats.norm` (`lib551.py:656`, every `level_rho`), and `ast` / `json` / `shutil` inside `gate551`. A venv upgrade of `cell-eval2` changes `wsum` with no blob failure. The pin is a git-blob pin; this is outside it. Worth one line in the frozen claim: the efp KA environment is the one that runs R7a.

**15. After `arm_gate`, C10b's window and the 2026-10-06 12:00Z limit do not apply.** RAN, gate log. Patched `window_ok` with `PHASE["name"] == "access"`: 2026-10-08 12:00Z open; 2026-10-04 08:00Z closed; 2026-10-07 closed; 2026-10-04 13:00Z closed; 05:59Z and 18:00Z closed (end exclusive). Unpatched `window_ok` still has the 10-06 limit and the C10b window (`gate551.py:67-77`). `PHASE`, `NO_START_AFTER`, `BLACKOUT`, `ACCESS` are read only inside the original `window_ok` / `thresholds` (`gate551.py:69, 73, 75, 146`). Child processes (whitener, `full551`, `comp551`, the axis check) do not inherit the in-process patch; they also do not call `gate551.gate`. The parent's `heavy()` does.

CD4 (`r7a551.py:323-327`), the `spec.json` rewrite (`r7a551.py:329-331`), and `aggregate` (`r7a551.py:382-388`) never call the gate. A re-run outside the window still executes an un-done CD4 and, once every heavy step is done, still runs the aggregation. CD4 and the aggregation read no perturbed PLOLO fold. The aggregation reads prices that were already computed. That is a light completion, not a second access, provided finding 7 is fixed so the command cannot change `--nb`.

**16. Axis check: what it proves, and what it does not.** `r7a551.py:297-316`, READ. The child `chdir`s to `experiments/534-generalization-audit/work/partB`, imports `axes534` (which imports `guard534`), recomputes `axis_sha(T bundle, DONOR_SETS[tag])`, and requires the sidecar json's `axis_sha256`, `n_cols`, and `donors` to match. `DONOR_SETS` (`axes534.py:26-27`) equals `r7a551.FITS` (`r7a551.py:56-57`). It then requires the digest not in `ABC_AXES` (two hardcoded hex strings, `axes534.py:28-29`), not in `tuned_on_axes()` (the nine paths in `TUNED_ON`, `axes534.py:30-33`, times four donor sets), and not in the four R6 panels recomputed from `partB/S*/bundle`. Not executed here: donor `genes.txt` and the bundles are not part of this run, and importing real `guard534` into the writer process would refuse writes under `experiments/555-cloud-r7-review`.

`axis_sha` / `tuned_on_axes` / `check_ops` (`axes534.py:43-84`) contain no write. `guard534` allows writes under `experiments/534-generalization-audit` and raises elsewhere in the repo. The child could write there if some called function wrote; these three do not. `check_ops` itself does not know about R6; the R6 exclusion is only in the child string. The check hashes gene names and the panel list, not the npz. Two T panels whose target lists differ only by genes outside the donor-gene union share an axis and both pass. RAN `repro/test_axis_predicate.py` (`repro/logs/axis_predicate.txt`): panels `{G4}` and `{G5}` shared a digest; `{G1}` did not. `TUNED_ON` is the 534 list, not every panel tuned since. R7's own bundles are not in it (if they were, the check would reject R7 for being R7).

**17. `plan`, `ka-cd4`, and `record` do not build a T-panel PLOLO fold.** READ of the call graph. `bundle_check` (`r7a551.py:114-128`) hashes control files, `pert_counts.csv`, `gene_names.csv`, and `manifest.json`, and compares the target column to `draws.json`. It does not call `build_fold`. `ka_cd4` (`r7a551.py:212-236`) runs the CD4 extractor on efp's bundle and deletes the copy. `record` calls `plan(write=False)`. `plolo.build_fold` (`plolo.py:514`) is not called. Importing `lib551` does import `plolo`, and `plolo.py:86-120` at import reads `data/challenge/pert_counts.csv` and `gene_names.csv` (the A/B/C panel, 300 and 18,533) and execs `experiments/062-cd4-post-whiten/run.py` and `experiments/064-pm-cd4-composition/run.py`, which define functions and do not build a fold at import. Those two files are in the 28 because that exec puts them in `sys.modules`; the AST closure missed them and the committed plan did not. The probe was not executed in this review (it imports `guard551`).

`access` is the step that runs `full551` / `comp551`, which do build folds. That is the reserved read. `access` cannot start with no `claim.json` (`r7a551.py:371`). It can start when the hash is only a substring of the ledger (finding 6). It cannot start when the hash is absent.

**18. Step differences from 534's R6 build, beyond the three intended ones.** Intended: no `--expect-sha256` (`r7a551.py:341-342` vs `drive551.py:124-126`; `drive534.py:128-129` also has no `--expect-sha256`); axis check in a child rather than `axes534.check_ops` in-process (`drive534.py:135-139`); no `--repro` (`drive551.py:140-144` passes `--repro`, and `r7a551.py:355-360` asserts `peak_rss_gib` and `op_axis_checked` only, not `ka.repro.worst`).

Same: donors and order of `FITS`; `--basis-prune essential`; `--scratch ops/.scratch`; `--max-rss-gib 14`; disk free > 10 GiB; spec fields `bundle`, `rules="def"`, `cd4_root`, `full_op`, `folds`, `synthetic=True`, `meta=None`, `fold_ops` (`r7a551.py:329-330`, `drive534.py:106-108`, `drive551.py:108-110`). `synthetic=True` means the controls are relabelled A/B/C frames (`lib551` bind), not that the biology is fake. CD4 flags match `cd4_534.py`: `--max-rss-gib 4.0`, `--overwrite`, gates on then `--min-targets 0 --min-union-targets 0` (`r7a551.py:168-172`). 534 writes `spec.json` once in `draw()` if absent (`drive534.py:110-112`); r7a rewrites it on every `panel()` call, including a resume. State also stores `axis_sha256` (`r7a551.py:346`), which 534 does not. Order: CD4, spec, four fits, axis check, full, comp, delete. 534's CD4 is a gated light step; r7a's is outside the gate (finding 15).

`comp551.py:302-306` skips an existing json when `peak_rss_gib` is present and the `ka` flag matches, and does not look at `config`. A re-run with a different F would skip. `aggregate_r7a551.py:91-94` then asserts `config == F` and stops. Fail closed at aggregation, not at comp. T panels do not start with `S`, so `comp551.py:296-299` does not require `--repro` or `c10b-config.json`.

**19. Plain and calibrated bootstraps see the same target draws.** `price551.pooled_bootstrap` (`price551.py:176-191`) builds `np.random.default_rng(seed)` and consumes it only in `rng.integers`. `cal` changes `k` / `k_pds` inside `step_prices` (`price551.py:155`) and does not touch the generator. Same seed, same `steps`, same panel order, same fold order ⇒ the same index matrix. RAN `repro/test_aggregate_r7a.py` (`repro/logs/aggregate.txt`): 240 index draws, plain and calibrated identical, calibrated prices exactly 10× the plain prices under a stub that scales by 10 when `cal` is set; seed 553 differed. `n_pos_folds` (`aggregate_r7a551.py:118`) sums `pooled[s]["per_fold"]` across steps. That object is a list (`price551.py:170-172`), so the count is the number of folds whose summed step effects are positive. A dict would hit `list(dict)` (the keys) and the sum would raise; the test provoked that TypeError. Per-panel `per_fold` is a dict and is not what this line reads.

PRIMARY ids from the same run, all six `CONFIGS` keys:

| F | primary | control |
|---|---|---|
| M5 | M0→M1 … M4→M5 | none |
| M6 | M5→M6 | none |
| M5+P2g | P2g@M5 | C-P@M5 |
| M6+P2g | M5→M6, P2g@M6 | C-P@M6 |
| M5+S_P2 | S_P2@M5 | C-S@M5 |
| M6+S_P2 | M5→M6, S_P2@M6 | C-S@M6 |

Verdict boundaries in that run match the claim: D>0 and lb>0 TRANSFERS; D>0 and lb≤0 FRAGILE; D≤0 FAILS.

**20. Retrodiction units.** `access553.py:412-427`: for `nmae`, `fid`, `reach`, `jac` only, `here` and `reference` are raw member effect means (arm − `cpm`), and `ratio = here/reference` unless `|reference| ≤ 1e-12`. `pds` and `mse` are absent. The reference path `work/access553_rehearsal5` is relative to `experiments/553-r7b-harness` (`access553.py:416-417`) and is the five-draw mean, not h1r7 alone. Board-likeness has to use the sign of `here`.

---

## A. `r7a551.py` against `drive551.panel`, `drive534.panel`, `cd4_534`, `axes534`

Checked, line by line, the panel body, the fit command, the spec fields, the post-comp asserts, the axis child, every read of the gate globals, resume, and the import pin.

Gate globals and whether the patch is visible: see findings 4 and 15 and `repro/logs/gate_patch.txt`. `window_ok`, `thresholds`, `GATE_LOG`, `HEAVY_LOG`, `floor_gib`, `OVERRIDE`, `PHASE`, `NO_START_AFTER`, `BLACKOUT`, `ACCESS` are module globals. `gate`, `wait_gate`, and `heavy` load `window_ok` / `thresholds` / `floor_gib` / the log paths at call time. Nothing in those functions binds them earlier. `arm_gate` covers `window_ok`, `thresholds`, and both logs. It does not cover `floor_gib` or `OVERRIDE`. With the patch applied, no path still applies C10b's window or the 10-06 limit. Without the patch, both still apply; `main` only installs the patch inside `access()` (`r7a551.py:378`), after `plan(write=False)`.

Resume and the second `record`: findings 5 and 6, `repro/logs/resume_and_seal.txt`. A fit that dies with rc ≠ 0 is not marked done. A crash after `G.heavy` returns and before `done()` re-runs the fit. A crash after `done()` and before the sha update (`r7a551.py:343-347`) skips the re-fit; the npz is already there and the axis check still runs if `axes` is not done. Operators deleted before every comp json exists: the comp step is not done, so resume re-enters `comp551`, which skips a json that already has `peak_rss_gib` (finding 18) and recomputes the missing fold. Then the unlink crash in finding 5 still applies.

Sealing of `plan` / `ka-cd4` / `record`: finding 17. `access` without the hash present in the ledger: refused (`r7a551.py:374`). `access` with the hash present only as a quoted warning: allowed (finding 6).

Lazy imports: RAN `repro/scan_lazy_imports.py`. Closure of module-level repo imports is 27 files. The committed pin has 28 because `plolo.py:119-120` loads `062` and `064` through `importlib` into `sys.modules` (the AST walker does not see that). Repo lazy imports inside functions:

| site | on the T-panel path? | in the 28? |
|---|---|---|
| `lib551.py:642` `ntclevel417` | yes (`full551.py:225`) | no — BLOCKER 1 |
| `mse187.py:87` `mech_a12_audit` | no (`mse187.main` only) | no — NOTE 13 |
| `full551.py:118` `lodo415.whiten_split` | only if that function runs; the module is imported at module level elsewhere, and the file is pinned | yes |
| `full551.py:221` `analyse418` | only `panel.startswith("abc")` | yes, via `price551` |

`aggregate551.py` is exec'd by the probe and then omitted from the 28, because the probe module is not inserted into `sys.modules` and the explicit add list is only `full551.py` and `comp551.py` (`r7a551.py:105`). R7a does not run `aggregate551.py`. Changing it does not change R7a. `aggregate_r7a551.py` is hashed in `external` only when the file exists (`r7a551.py:147-148`). The committed `plan.json` lists three external paths and not the aggregation, because the plan predates the file. `record` calls `plan(write=False)` and stores a fresh external map in `claim.json`. `access` compares against that map (`r7a551.py:377`). Not a hole if `record` is the thing that arms.

## B. `aggregate_r7a551.py`

Checked against `aggregate551.py` and `price551.py`. RAN `repro/test_aggregate_r7a.py` with the venv (system Python has no scipy; `price551` → `analyse418` → `joint187` imports `scipy.stats`). The test stubs `guard551`, loads `designs551` and `price551` from the tree, loads `aggregate_r7a551.py` from the snapshot by `importlib` so `HERE` stays `price551`'s, and replaces `Panel` with a fake whose `agg_fn` returns a pure function of the index vector.

(1) PRIMARY for each key: the table in NOTE 19. F = M5 falls back to cumulative M0→M5 (finding 12). F = M6 is only M5→M6, not the cumulative. The cumulative is `secondary` (`aggregate_r7a551.py:129`).

(2) Pairing: NOTE 19. The two `pooled_bootstrap` calls (`aggregate_r7a551.py:100-101`) pass the same `SEED` 552, the same `steps`, and the same `cum_ids`. `cal=True` does not consume the RNG.

(3) `n_pos_folds`: NOTE 19. It counts folds, and the denominator printed in the table is `3 * n_panels` (`aggregate_r7a551.py:168`).

(4) Retention: finding 10. Units on the R6 side are levelled overall means. Units on the A/B/C side are levelled abc means for M0→M1, M1→M2, M2→M3 and the record's text numbers after that. C10b's step is the same levelled overall, included only when the recorded config equals F.

(5) Verdict inputs that can be chosen after the number: finding 7 (`--nb`, `--panels`, a second write). The verdict function itself (`D>0 and lb>0`) has no other switch. Calibrated rows, targeting, and retention are not inputs to `verdict`. `--test-config` cannot write the real result.

## C. The claim draft against RESERVE.md and the ledger

RESERVE.md §4: one claim, once, sha256 in the ledger before the first R7 number, no fitting on an R7 score. Amendments 1–3: infected state stands; veto only on a preregistered HARM beyond the surface's own re-roll band, null or small never counts against a switch; sealed versus control-only inputs; floor 5000, BST2 excluded, <30 cells dropped. The ledger's R7 section (from line 1024) records the declaration and the three amendments and does not record a claim sha. The draft's first line matches that.

(1) PRIMARY is the right single question for `lodo_hi_pct` and `sharpen_pt`, which are the only post-C10a bulk switches the six `CONFIGS` keys can express: F versus M5, the sum of the steps those switches add. It is not the right question for a bulk switch outside those keys (finding 11), and it is the wrong estimand when F = M5 (finding 12).

(2) Board-likeness and the HARM veto: findings 8 and 9. The code does not apply either rule. A person reads `access553.json` and `band553_*.json`. That is a post-hoc choice unless the frozen sentence is the one that is applied and the ledger hash is checked before the truth opens (finding 2).

(3) Choices the draft still leaves open after a number exists, or that the code will accept after a number exists:

- F and the arm list are marked FILL. That is legitimate only if the file is frozen and hashed before any R7 score. The driver does not check that the hashed file contains the F it was given (finding 6).
- The tercile sentence (finding 11).
- `--nb` / `--panels` / a second aggregate (finding 7).
- New `--pairs` or `--retro` on `--phase report` after the access (`RESULTS.md:527`; `access553.py:398-427`).
- `band_total` versus `band_sd_ctrl` (finding 8). Both will be on disk.
- Board-like applied to `ratio` instead of `here` (NOTE 20).
- Whether a missing `band_total` (finding 3) is "no veto" or "go get another band."
- `retention` returning None when C10b's F differs, and someone substituting a different R6 number by hand (finding 10). The code's answer is None.
- The window timestamps are arguments to `record`, not lines of the claim. That is before the score if `record` precedes `access`. It is not fixed by the claim hash.

## D. `access553.py r7bd5000 --owner-go … --pairs … --ctrl-splits 5 --retro … --halves 0,1,2,3,4`

Static review plus the stubs above. The floored truth is not on the path this review executed.

(1) `--owner-go`: finding 2. It verifies that the path exists, stores sha256 in the plan, and re-hashes the file at each phase. It does not verify a ledger entry. Yes, it runs against a claim whose sha256 is not in the ledger. `--dry-run` does not require `--owner-go` (`access553.py:532`) and does not open the truth (`access553.py:600-602`). A sealed dry-run still calls `check_r7b_outputs(check_truth_sha=True)` (`access553.py:536-537`), which hashes the unfloored `truth_r7b.h5ad` and does not parse it (`lib553.py:220-222`). That is a hash, not a score.

(2) An arm named in a switch pair and not emitted on the halves: finding 3. If the arm is in `--half-arms` and the stats file is missing, readiness raises `SystemExit` before the truth (`access553.py:595-596`). If the arm is only in `--arms` / `--pairs`, readiness passes, the main access scores it, and `band553` skips it.

(3) `band_total`: `band553.py:58-64`, as in finding 8. The claim's formula matches. `band_sd` is the main access's seeds + split-half of target cells (`effect_scaled.overall.band_sd`). `sigma_c` is the control-half difference, seed-corrected. `band_sd_ctrl` is not an input. h1r7's json has no between-draw sd because `--draws` was not passed; the claim's between-draw numbers live in `draws553.json`.

(4) `--retro`: NOTE 20. Units are raw member deltas. The reference `work/access553_rehearsal5` is the five-draw mean (CP raw jac +0.004116, raw fid −0.001777). It is not a single draw and it is not H1-full (`--retro work/access553_h1` is the other file RESULTS names).

---

## Commands

```text
/workspace/.venv/bin/python experiments/555-cloud-r7-review/repro/calibrate_harm.py
/workspace/.venv/bin/python experiments/555-cloud-r7-review/repro/test_owner_go.py
/workspace/.venv/bin/python experiments/555-cloud-r7-review/repro/test_band_skip.py
/workspace/.venv/bin/python experiments/555-cloud-r7-review/repro/test_axis_predicate.py
/workspace/.venv/bin/python experiments/555-cloud-r7-review/repro/test_aggregate_r7a.py
/workspace/.venv/bin/python experiments/555-cloud-r7-review/repro/test_resume_and_seal.py
python3 experiments/555-cloud-r7-review/repro/test_gate_patch.py
python3 experiments/555-cloud-r7-review/repro/scan_lazy_imports.py
```

Logs: `repro/logs/`. The aggregate test needs the venv (scipy). The gate test and the AST scan do not. `test_resume_and_seal.py` writes under `repro/resume_out/` and that directory was removed after the log was saved.
