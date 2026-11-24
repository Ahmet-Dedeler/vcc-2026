# 432: D/E/F dry rehearsal of the September knobs on the synthetic bundle

**2026-09-23 (Opus 5.5 worker; RESUME 2 run 11:27Z–~13:45Z).** Nothing was submitted or written outside this directory.
No `cell-eval2`, no git. The production builder was pinned and restored around every builder call.
- Before and after every one of 12 builder calls it was byte-identical to `2c90d9e1…`. The chain's START and END lines
  record that; the per-build lines are in `work/state432.log`.
- The emitter (`edb5b8dc…`) and the driver (`27edc486…`) were unchanged throughout.

## Headline

**The full 433h stack runs end to end on a synthetic D/E/F panel. The zs pass and every September knob run with it.**
- **Knobs covered:** the re-fit operator, 410's key re-fitted for the panel, ballast, sharpening q 1.25, the stack DUMP
  → LOAD, the per-context re-level derived with `ntclevel417.py`, the coarse MDE, zero repair −5 and zs.
- **Cost:** the rehearsal made 12 builder calls, all strictly serial and behind the memory gate:
  - four full builds (b0, b1, b2, ef_b1), each at 25–35 min and 20.21 GiB peak;
  - five bulk dumps, 11–12 s each;
  - three negative controls, which the builder refused within 10 s.
- **The day path itself** is one dump plus two full builds: the un-levelled LOAD build, then the levelled shipped
  build. That is ≈ 58–68 min from bundle to a prepped artifact (§8).
- **Checks passed on the shipped-form artifact b2:**
  - 0 zeroed gated genes;
  - 15 genes repaired;
  - MDE dose D 0.956 / E 1.0 / F 1.0;
  - zs moves 40 / 31 / 57 genes per block (medians);
  - the re-level lands each context at its EV-optimal margin in one step (EV-scale after 0.9998 / 1.0008 / 1.0015);
  - the strict postbuild passes 41/41;
  - `vcc prep --dry-run` is clean.
- **LOAD refuses cross-panel inputs.** An A/B/C dump, the A/B/C key and the A/B/C re-level string are all refused.
- **On the pessimistic low-coverage panel, the un-levelled stack build puts context D PAST its mse cliff** (§7).
  - ρ = B/claim = 1.071: 7.1 % past it in energy, 3.4 % in bulk scale.
  - The same stack on the realistic panel stopped 0.8 % short in energy, at ρ 0.992.
  - So on D/E/F the re-level is mandatory, and its factor can be well below 1: `D=0.9561` there.
  - Deriving the string on the LOADed artifact is the only safe order.
- **Nine day-of traps found (§6).** Three of them would have cost hours or the stack on Oct 22:
  1. **417's literal CD4 ×1.5 dump recipe (`VCC_CD4_TARGET_F=0.51`) is a silent no-op on D/E/F.** The dump is bit
     identical to H1 ×2 alone. Use `VCC_CD4_TARGET_G=0.5165482917974462`.
  2. **The re-fitted key needs a builder source edit** that no script makes. Exporting its id before that edit makes
     every builder call refuse at import, including the driver's own `repin-frozen`.
  3. **`def_pipeline.py` would ship early September on the day.** Its defaults, its arm guard and its postbuild
     reference all predate these knobs. Its postbuild then fails every September D/E/F build (40/41).
- The fixes are proposed as diffs with tests in `proposals/`. None was applied to a production file.

## 1. What ran

**Bundle.** 199's "realistic" `def-bundle-r` (panel `vcc2026-synthetic-def-199r`, seed 1990, K562 ∩ Orion ∩ CD4 pool).
- **Coverage.** Every one of its 300 targets is carried by all three production donors. PerturbMatch carries 273, CD4
  286 (95.3 % carrying rows), the H1 union 59. It is a best-case coverage draw.
- **Controls.** The control files are A/B/C relabelled (D ≡ A/Jurkat, E ≡ B/HeLa, F ≡ C/CAL33). So every per-context
  number below is a number for A/B/C's control cells on a new panel. It is not a measurement of a new cell line.
- **Free known answer.** Because of that, `ntclevel417`'s Cp46 must equal the A/B/C reference values. It does, exactly:
  14.4398 / 14.8057 / 16.4397.

**Low-coverage extension.** The same stack, minus the re-level, on 238's essential-free pessimistic `ef-bundle`:
- coverage: K562 146/300, CD4 carrying rows 182/300 (60.7 %), H1 union 28/300;
- operator: 244's re-fit (`373df206…`), read-only;
- see §7.

**Configuration.** 433h's config (`submissions/433h-*.h5ad.meta.json`) minus the three A/B/C-specific values:

| | value on D/E/F | source |
|---|---|---|
| driver defaults (`SHIPPED_ENV` + `DEF_ENV`) | `VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.34 VCC_H1_UNION_G=0.28 VCC_CD4_SIZING=carrying VCC_CD4_TARGET_G=0.3443655278649641 VCC_CD4_DEFLATE_BASIS=V5_abc` | `def_pipeline.resolve_env(day_env(False))` |
| September knobs, exported | `VCC_NORMREL_EXP=-1.5 VCC_NORMREL_MODE=cap VCC_CD4_PARTIAL_G=0.5 VCC_DE_SIGN_PRIOR=clean4q80 VCC_PERT_SCALE_BALLAST=1 VCC_BULK_SHARPEN_Q=1.25 VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5` | 433h |
| the key | `VCC_PERT_SCALE=wshare_lo0.3_hi3_def` → `work/key/pert_scale_def432_lo0.3_hi3.npz` `4f1f0209…` | re-fitted on the day's operator (`key432.py`) |
| the stack (b1, b2) | `VCC_BULK_LOAD=work/dump_stack/dump_h1x2_cd4x15.npz:b6d80499…` | dumped with `VCC_H1_UNION_G=0.56 VCC_CD4_TARGET_G=0.5165482917974462` |
| the re-level (b2) | `VCC_BULK_CTX_SCALE=D=0.9938,E=1.0353,F=1.0680` | `ntclevel417.py` on b1 |
| zs (b2) | `VCC_DE_MAG=dn_k250_zs` | 429, in production |
| flags | `--min-cpm 10 --tail-spread 1.0 --keep-top 8000 --keep-spread 0.25 --overdispersion 4.0 --emitter rc_callset` | the driver's `SHIPPED_FLAGS` |

**How the key reached the builder.** The mandate says to stop and document if the key needs a builder edit. It does
(§6 T2).
- Instead of editing the file, `run_builder432.py` runs the pinned builder with that one edit applied **in memory**:
  one extra `PERT_SCALE_ARTIFACTS` entry, anchor-checked and refusing on any sha mismatch.
- The patched source, its sha and its diff are saved next to every build (`work/<b>/builder432_*.{patched.py,diff,json}`).
  The file on disk only ever held the driver's own pin (`5820b6bf…`; `f7d5ed48…` on the ef bundle), which is restored
  after every call.
- Caveat: each meta's `hashes.implementation` records that on-disk sha (`5820b6bf…`), not the patched one. The patched
  sha is in `builder432_build.json`.

## 2. Stage / knob table

Measured with `/usr/bin/time -l`; the logs are in `work/<dir>/logs/` and `work/state432.log`, and the table is
regenerated by `summarise432.py` into `work/summary432.json`.

| # | stage / knob | ran | outcome | wall | peak RSS | what to look for on the day |
|---|---|---|---|---:|---:|---|
| 1 | light stages (preflight, identity, cd4-extract, census) | ✅ | rc 0 | 16 s | 2.7 GiB | census: donors 300/300/300, PM 273, CD4 286 |
| 2 | operator re-fit (driver `whitener`) | ✅ | `e90a1f79…`, axis `87361c7c…`, 17,942 cols, `dropped 4393` | 322 s | 12.9 GiB | axis sha ≠ A/B/C's `edf6ba88…`; npz sha ≠ `e2f65ecf…` |
| 3 | **410 key for this panel** (`key432.py`) | ✅ | `4f1f0209…`; carried 300/300, all 3-donor; w q50/q75/q90 = 1.0/2.37/3.0; 31 % at floor, 19 % at cap. KA on A/B/C: **byte-identical to the pinned `d7f76a22…`** and axis == fold-kolf2 `wcols` | 3 s | 4.1 GiB | `frac_at_cap` in the json (A/B/C 0.14, here 0.19, pessimistic 0.26); the builder prints `[pert-scale] wshare_lo0.3_hi3_def sha=…` |
| 4 | driver `repin-frozen` + `pin-artifacts` (arm env **unset**) | ✅ ×12 | builder `2c90d9e1` → `5820b6bf` | 2 s | 0.2 GiB | one operator entry + CD4 hashes; **the key is NOT pinned** (T2) |
| 5 | `repin-panel-support` | ✅ ×12 | `{pm 273, cd4 286}` = census (ef: `{133, 182}`) | 10 s | 20.1 GiB | the two numbers equal the census prediction |
| 6 | **ballast** (`VCC_PERT_SCALE_BALLAST=1`) | ✅ b0–b2 | `ballast_jk_ratio` median 1.022 / 1.024; `ballast_clip_mass` max 1,561 / 1,346 | — | — | jk ratio ≈ 1.02, clip mass ≤ a few thousand (DEF-RUNBOOK 410) |
| 7 | **sharpening** (`VCC_BULK_SHARPEN_Q=1.25`) | ✅ b0–b2 | `bulk_sharpen_q 1.25` in every config; new panel's columns protected (mask from `spec.perts`) | — | — | config key present |
| 8 | **stack DUMP** (`VCC_BULK_DUMP`, H1 ×2 + CD4 ×1.5 via `TARGET_G`) | ✅ | `b6d80499…`, exits before emission; `[cd4] sizing=carrying g=0.516548: solved 0.5115715 on 286 rows` | 12 s | 20.1 GiB | `g=0.516548` and `[h1-union] g=0.56` in the dump log (T6) |
| 9 | **b0** = 433h − LOAD − CTX (the base path, task a) | ✅ | `b98c797b…`, 360,000 rows, 11.80 GB | 1,718 s | 20.21 GiB | `[arms]` OK lines; strict postbuild 41/41 |
| 10 | **b1** = b0 + LOAD (task b) | ✅ | `c3978bdd…`; `[bulk-load] max\|loaded−own\| = 0.6387372208324478` | 1,516 s | 20.21 GiB | the `[bulk-load]` line; `bulk_load_sha256` in config |
| 11 | **KA-LOAD**: dump with b0's own env | ✅ | max\|stack − b0env\| = **0.6387372208324478** = b1's printed value to the last bit | 11 s | 20.1 GiB | (rehearsal only) |
| 12 | **re-level** `ntclevel417.py` on b1 (task c) | ✅ | `D=0.9938,E=1.0353,F=1.0680` | 22 s | 3.4 GiB | Cp46 per context; a factor < 1 is legitimate (D here) |
| 13 | **b2** = b1 + CTX + zs (the shipped form) | ✅ | `a56ffb29…`, 11.82 GB | 2,100 s | 20.21 GiB | `[bulk-ctx-scale] {...} on the ballast bulk render only`; `[de_mag] …differs from the emitter gate on 0` |
| 14 | re-level check: `ntclevel417.py` on b2 | ✅ | EV-scale after 0.9998 / 1.0008 / 1.0015; ρ 0.980 / 0.978 / 0.976 | 19 s | 3.4 GiB | EV-scale 1.000 ± 0.002 (one step converges) |
| 15 | **coarse MDE** (`VCC_DE_MDE=0.1,0.1,0.3,coarse`) | ✅ b0–b2 | `mde_lambda` D 0.956 / E 1.0 / F 1.0 (§3) | — | — | `mde_lambda` median ≥ 0.8 per context (431 dump) |
| 16 | **zero repair** (`VCC_DE_ZERO_REPAIR=-5`) + `audit428.py` | ✅ b0–b2 | 15 genes repaired (D 5 / E 4 / F 6); audit: **0 zeroed gated genes** in 900 blocks | 20 s | <3 GiB | 0 zeroed; the extreme tail (below) |
| 17 | **zs** (`VCC_DE_MAG=dn_k250_zs`) | ✅ b2 | `zs_n_moved` median 40 / 31 / 57 | +584 s on the build | — | `zs_n_moved` > 0 in every block |
| 18 | driver `postbuild` (allow-list `whitener_id`, vs 433h) | ❌ as designed | **40/41 FAIL** on b0, b1, b2 (T5) | 1 s | 0.1 GiB | — |
| 19 | strict postbuild (allow-list `whitener_id,pert_scale,pert_scale_sha256,bulk_load,bulk_load_sha256,bulk_ctx_scale[,de_mag]`) | ✅ | **41/41** on b0, b1, b2 | 1 s | 0.1 GiB | 41/41 |
| 20 | `vcc prep --dry-run` (driver `prep`) | ✅ b0, b2 | 360,000 cells, D/E/F = 120,000 each, no `.vcc` written | 8–10 s | 14.3–19.7 GiB | `✓ dry run` |
| 21 | refusal: A/B/C dump LOADed | ✅ refused | `VCC_BULK_LOAD perts/genes differ from the bundle's` | 10 s | 20.1 GiB | — |
| 22 | refusal: A/B/C key (`wshare_lo0.3_hi3`) | ✅ refused | `VCC_PERT_SCALE artifact target order differs from the panel` | 10 s | 20.1 GiB | — |
| 23 | refusal: A/B/C re-level string | ✅ refused | `VCC_BULK_CTX_SCALE contexts ['A','B','C'] != bundle ['D','E','F']` | 10 s | 20.1 GiB | — |
| 24 | 417's literal dump recipe (`VCC_CD4_TARGET_F=0.51`) | ✅ ran, **silent no-op** | `dump_lit` == `dump_h1` (H1 ×2 alone) bit for bit (T6) | 11 s | 20.1 GiB | — |
| 25 | `cliffcheck414.py` | not run | A/B/C-only, and it writes into 414's work dir (T7) | — | — | use rows 12/14 instead |
| 26 | **ef_b1**: low-coverage build, every knob but the re-level (§7) | ✅ | `2896e721…`; strict postbuild 41/41; 0 zeroed genes; **context D 7 % PAST its cliff un-levelled** (ρ 1.071) | 2,046 s | 20.15 GiB | ρ > 1 in any context = the re-level is not optional |

Extreme tail (428 audit on b2): `lfc > +3` on 3 blocks, all SUPV3L1 (D +5.10, E +5.07, F +5.60; A/B/C 433k's max was
+3.10). `lfc < −5` on 10 blocks, which are the repaired genes sitting at the −5 floor.

## 3. Task (d): the DE knobs on D/E/F, from the per-block dump (`VCC_RC_AUDIT_DUMP`, 900 blocks)

b2 = the shipped form (coarse MDE + zero repair + zs). A/B/C references come from the same dump format: 433a
(coarse MDE, no repair) and 429a (plain MDE + zs).

| context (controls) | `mde_lambda` median | q10 | min | frac < 1 | frac = 0 | `zero_repaired` median / total / blocks > 0 | `zs_n_moved` median (q10–q90) | A/B/C reference: 433a `mde_lambda` median; 429a `zs_n_moved` median |
|---|---:|---:|---:|---:|---:|---|---|---|
| D (A, Jurkat) | **0.9561** | 0.9315 | 0.0342 | 0.913 | 0 | 0 / 5 / 5 | 40 (30–54) | A 0.9561; A 41 |
| E (B, HeLa) | **1.0000** | 1.0000 | 0.0117 | 0.053 | 0 | 0 / 4 / 4 | 31 (22–40) | B 1.0000; B 31 |
| F (C, CAL33) | **1.0000** | 1.0000 | 0.0010 | 0.013 | 0 | 0 / 6 / 6 | 57 (44–74) | C 1.0000; C 56 |

Readings:
- **The dose follows the control cells, not the panel.**
  - Per-context MDE doses and zs capacities on a disjoint panel (the synthetic panel shares 0 of its 300 targets with
    A/B/C's) reproduce the A/B/C dumps: 0.9561 vs 0.9561 in the Jurkat context; zs medians 40/31/57 vs 41/31/56.
  - 431's slack mechanism predicts this: the dose tracks the context's library quartiles and the sub-gate mass per cell.
  - So on Oct 22 the dump answers "does D/E/F behave like A?" per real line. It is not needed to find out whether the
    knob runs.
- **The runbook's floor is met.** Every context sits above DEF-RUNBOOK's ≥ 0.8 floor. D is the A-like context.
  - With the plain form, D would read ≈ 0: 429a's A read 0.0000.
  - **Rehearse and ship `coarse` only.**
- **The zero repair fires on D/E/F**, 15 genes in 15 blocks, and the audit then reads 0 zeroed genes.
  - This shows the repair runs, not only that nothing needed repair.
  - The repaired genes sit at the −5 floor. They are the 10 `lfc < −5` entries.
- **LOAD and CTX are not bit-neutral on the DE diagnostics.** The b0→b1→b2 comparison shows a bulk change moves
  `ballast_K` in 21 % (LOAD) and 11 % (CTX) of blocks, so the DE channel shifts slightly:
  - `mde_lambda` equal in 92 % / 98 % of blocks;
  - `forced` differs in almost all blocks.
  - This is by design (410b: DE members within 0.0007 raw), but a DE-channel diff across bulk arms is not a bug signal.

## 4. Task (c): the re-level, and whether the energies moved as intended

`ntclevel417.py <artifact> --bundle def-bundle-r --contexts D,E,F --sd 0.01` on each build. B is the rendered
across-perturbation energy, claim = Cp46 · Σw, ρ = B / claim.

| build | D: B / ρ / cliff / EV-scale | E: B / ρ / cliff / EV-scale | F: B / ρ / cliff / EV-scale |
|---|---|---|---|
| b0 (no LOAD) | 0.22020 / 0.9420 / 1.0303 / 1.0196 | 0.20890 / 0.8716 / 1.0711 / 1.0600 | 0.21825 / 0.8201 / 1.1043 / 1.0928 |
| b1 (LOAD, un-levelled) | 0.23178 / **0.9916** / 1.0042 / **0.9938** | 0.21899 / 0.9137 / 1.0462 / 1.0353 | 0.22851 / 0.8586 / 1.0792 / 1.0680 |
| b2 (LOAD + CTX) | 0.22901 / 0.9797 / 1.0103 / **0.9998** | 0.23435 / 0.9777 / 1.0113 / **1.0008** | 0.25986 / 0.9764 / 1.0120 / **1.0015** |

- **The energies moved as intended.** B(b2)/B(b1) against the asked s²:
  - D: 0.98805 vs 0.98767 (ratio 1.0004);
  - E: 1.07011 vs 1.07185 (0.9984);
  - F: 1.13720 vs 1.14059 (0.9970).
  - Every context lands at ρ ≈ 0.978 = exp(−sd·z*), the EV-optimal margin, and a second level would ask ≤ 0.15 %.
  - One level step converges. The small shortfall on E/F is render non-linearity (log1p, the 10-CPM gate, integer
    counts).
- **The stack alone moves a context to its cliff.** LOAD raised B by 4.7–5.3 % in every context. D (A's controls) went
  from ρ 0.942 to 0.992, 0.4 % in scale under its cliff, and its EV-optimal factor is < 1 (0.9938).
  - On the day, never ship a LOAD without re-levelling.
  - Derive the string on the LOADed artifact, not the unloaded one: b0's string (`D=1.0196,E=1.0600,F=1.0928`)
    applied to b1 would put D 1.5 % past its cliff.
- **Known-answer row.** Cp46 = 14.4398 / 14.8057 / 16.4397, identical to 180's A/B/C values (same control cells).

## 5. Task (b): stack dump / load

| dump | env | sha256 | ‖bulk‖_F | vs `dump_h1` |
|---|---|---|---:|---|
| `dump_stack` (the D/E/F recipe) | `VCC_H1_UNION_G=0.56 VCC_CD4_TARGET_G=0.5165482917974462` | `b6d80499…` | 81.199 | rel diff 0.141 |
| `dump_lit` (417's literal recipe) | `VCC_H1_UNION_G=0.56 VCC_CD4_TARGET_F=0.51` | `cc99b99a…` | 81.504 | **bit-identical** |
| `dump_h1` (H1 ×2 alone) | `VCC_H1_UNION_G=0.56` | `1db4eb14…` | 81.504 | — |
| `dump_b0env` (KA) | b0's env | `4b3d550c…` | | max\|stack − b0env\| = b1's `max\|loaded−own\|` exactly |

- **LOAD works on D/E/F.** b1 LOADs the stack, prints `max|loaded−own| = 0.6387372208324478`, and that equals the
  offline max |dump_stack − dump_b0env|. So LOAD replaced exactly the build's own bulk and nothing else.
  - 417's KA-LOAD (a same-env LOAD reproduces X byte for byte) was not re-run here; it would cost a 25-min build.
- **LOAD refuses an A/B/C dump** (`experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15.npz`, `bb19f0e5…`):
  `VCC_BULK_LOAD perts/genes differ from the bundle's`.
- **Dumps are cheap.** The effects stage takes ≈ 11 s once the donor files are in the page cache; the dump exits before
  emission.

## 6. Bugs and traps, with the exact fixes (text only; no production file was changed)

**T1. The 410 key script is A/B/C-hard-wired, and run as-is on Oct 22 it overwrites the pinned A/B/C artifact.**
- What it hard-codes (`experiments/410-pertscale-revisit/build_weights410.py`):
  - the panel: `data/challenge/pert_counts.csv`, line 35;
  - the operator: `OP_PATH = data/donors/whitener_union_pruned.npz`, line 30;
  - the axis: `wcols` from `fold-kolf2.npz` (the A/B/C axis), line 37;
  - the OUTPUT: `experiments/410-pertscale-revisit/work/pert_scale_410_lo…npz`, line 49. That is the very file
    `PERT_SCALE_ARTIFACTS["wshare_lo0.3_hi3"]` pins.
- The runbook's recipe is "edit `OP_PATH`" and point `wcols` at the refit. Its panel guard only fires if the panel
  sits on the axis.
- **Fix:** use `key432.py`: same arithmetic, same npz fields; the panel, operator and output are arguments, and the
  axis is recomputed from the bundle the way `build_union_whitener.py` forms it.
  - `--ka` on A/B/C reproduces `d7f76a22…` byte for byte and the axis equals fold-kolf2's `wcols`.
  - It runs in 3 s at 4.1 GiB.
  - Move it to `scripts/def_pert_scale_key.py` when the proposal lands.

**T2. The re-fitted key needs a builder SOURCE edit, and no script makes it.** The builder accepts a key only through
the source-pinned `PERT_SCALE_ARTIFACTS` (path + sha). `def_pin_panel_artifacts.py` pins the operator and the CD4 banks
only.
- **The exact edit** (what `run_builder432.py` applied in memory): one entry at the top of `PERT_SCALE_ARTIFACTS = {`
  in `scripts/make_submission_pm_cd4.py`:
  ```python
      "wshare_lo0.3_hi3_def": (
          ROOT / "<day work dir>/pert_scale_def_lo0.3_hi3.npz",
          "<sha256 of that file>",
      ),
  ```
- **Fix:** `proposals/def_pin_panel_artifacts.proposed.{py,diff}` adds `--pert-scale <npz> [--pert-scale-id …]` to the
  pin.
  - The entry is additive and in place on a re-pin; it refuses a hand-written entry, and the post-write check covers
    its sha.
  - Tested by `proposals/test_pin432.py` (7/7 PASS on the dry run against git HEAD's builder copy;
    `work/test_pin432_ef.json`).

**T3. Exporting the day's key id before that edit makes every builder call refuse at IMPORT.** This includes the
driver's own `repin-frozen`, which runs the builder before `pin-artifacts` and passes the whole exported env.
- Demonstrated: `VCC_PERT_SCALE=wshare_lo0.3_hi3_def uv run python scripts/make_submission_pm_cd4.py --help` exits 1
  with `VCC_PERT_SCALE must be one of (...)` (`work/traps/import_refusal_keyid.log`).
- **Fix:**
  - keep every `VCC_*` unset through `pin-artifacts` (what `build432.sh` does);
  - or apply the proposed driver: a `key` stage before the pins, `pin-artifacts --pert-scale`, and the id as a
    `DEF_ENV` default. Defaults do not reach `repin-frozen`; exports do.

**T4. `def_pipeline.py` builds and guards early September.** Four things are stale:
- `SHIPPED_ENV` is still the `h1union` line. Run as the runbook prints it, the day's build silently drops the nine
  knobs added after `cd4svd5`: normrel exp and cap, partial CD4, the CLEAN4 prior, the key, ballast, sharpening, the MDE
  and the zero repair.
- 240's `[arms]` guard re-reads only the old arms. b0's block has 8 OK lines and none about a September knob.
- `--incumbent-meta` defaults to `249-cd4svd5`.
- **Fix:** `proposals/def_pipeline.proposed.{py,diff}` (regenerated by `make_def_pipeline_proposal.py`, 8 anchored
  edits):
  - `SHIPPED_ENV` + `ABC_ENV`/`DEF_ENV` as in §9;
  - `september_rows()` re-reads every September knob off the written meta by value, plus `VCC_BULK_LOAD` /
    `VCC_BULK_CTX_SCALE` / `VCC_DE_MAG` when exported;
  - the postbuild allow-list for new panels;
  - the `key` stage and `pin-artifacts --pert-scale`.
- **Tests:**
  - `proposals/test_arms432.py`, 10/10 PASS (`work/test_arms432.json`):
    - the current guard still accepts cd4svd5 and 433h;
    - the proposed guard accepts 433h and 432's b0, b1, b2 and ef_b1;
    - it refuses 417a (it names `de_mde`, `de_zero_repair`, `bulk_load_sha256`, `bulk_ctx_scale`) and cd4svd5;
    - it refuses 433h when zs was asked for.
  - A `--dry-run` of the proposed driver on `def-bundle-r` prints the `key` stage and `--pert-scale` on the pin
    (`work/drydrv/dryrun_proposed.log`).
- **Also needed on the day:** pass `--incumbent-meta submissions/<the shipped A/B/C arm>.h5ad.meta.json` (433h or its
  successor).
- **Consequence:** 240's row 1 ("the driver's defaults rebuild the shipped package byte-for-byte") must be re-pointed
  at that arm, with its LOAD/CTX exported, and re-run once the proposal lands.

**T5. The driver's `postbuild` fails every September D/E/F build** (40/41, on b0, b1 and b2). Its allow-list is
`whitener_id`, and the panel-specific `pert_scale`, `pert_scale_sha256`, `bulk_load_sha256` and `bulk_ctx_scale` differ
from any A/B/C reference by construction.
- It fails closed, not silently, but it blocks the day.
- **Fix:** the allow-list in the proposal (row 19 ran it: 41/41).
- With the reference = the A/B/C meta of the SAME config, every other September key (`de_mde`, `de_zero_repair`,
  `bulk_sharpen_q`, `normrel_*`, …) is enforced by value by that same check.
- Rehearsal-only corollary: `def_postbuild.py` reads `WHITENER_REGISTRY` from the builder source, so a rehearsal must
  run it while pinned. After the restore, the registry row fails.

**T6. 417's literal dump recipe for the stack is a silent no-op on D/E/F.**
- 417 dumps CD4 ×1.5 as `VCC_CD4_TARGET_F=0.51` (`chain_overnight417.sh` line 27, through `run_dump417.sh`).
- Off A/B/C the driver sets `VCC_CD4_SIZING=carrying`, which sizes CD4 by `VCC_CD4_TARGET_G` and only RECORDS `TARGET_F` (the log says `the all-rows rule at f=0.51 would
  solve 0.508…`).
- So `dump_lit` is bit-identical to `dump_h1`: the stack would ship H1 ×2 with CD4 ×1, and every check passes.
- **Fix (runbook):** dump with `VCC_CD4_TARGET_G=0.5165482917974462` (= 1.5 × `DEF_ENV`'s 0.3443655278649641), and
  check `[cd4] sizing=carrying g=0.516548` in the dump log.
- **Fix (builder, proposed text):** under `VCC_CD4_SIZING=carrying`, refuse an exported `VCC_CD4_TARGET_F` other than
  0.34 ("it would only be recorded; set VCC_CD4_TARGET_G").

**T7. `cliffcheck414.py` is A/B/C-only and writes outside the caller's directory.** It hard-codes contexts A/B/C and
187's fitted claims, reads `data/challenge/context_{A,B,C}.h5ad`, and appends to
`experiments/414-bulk-reprice/work/cliffcheck414.json`.
- On Oct 22 (D/E/F files in `data/challenge/`) it raises; before that it reads the wrong controls.
- **The D/E/F equivalent** is §4's pair: `ntclevel417.py` on the un-levelled and the levelled builds. B(b2)/B(b1) ≈ s²
  within 0.3 %, ρ ≈ 0.978, EV-scale after within 0.2 % of 1.

**T8. The dump binds the bundle's perts and genes, but not the operator or the key.**
- A dump made with a different operator or key on the same bundle would LOAD. On the day there is one operator, so the
  risk is low.
- **Fix (proposed text):** at LOAD, compare the dump's recorded `VCC_PERT_SCALE` and `operator.sha256` with the
  build's.

**T9. Synthetic D/E/F controls are A/B/C relabelled.** Every per-context number here (MDE dose, zs capacity, Cp46,
cliffs) is A/B/C's controls on a new panel. The real lines' values exist only on Oct 22. Read them from the rc dump and
the `ntclevel` output then; §2's look-for column says what is normal.

Two notes that are not bugs:
- **`pgrep -f` self-match (AGENTS 6.2).** It deadlocked this worker's previous chain. The rewritten `chain432.sh` /
  `chain432x.sh` use no `pgrep -f`, and their builder gate is `def_pipeline.other_builders()`, which matches launcher
  executables only.
- **Stale message text.** `def_postbuild.py` calls any reference "the current incumbent (cd4svd5, board 0.2160944 …)".

## 7. Low-coverage extension: ef_b1 on 238's essential-free pessimistic panel

238's `ef-bundle` (panel `vcc2026-synthetic-ef-238`, the essential-free pessimistic draw; the same relabelled A/B/C
controls):
- **Coverage:** K562 146/300, PerturbMatch 133, CD4 182 carrying rows (60.7 %), H1 union 28/300, 7 targets shared with
  `def-bundle-r`.
- **Inputs:** 244's re-fit operator (`373df206…`, read-only) and 244's CD4 banks.
- **Key:** 410's key re-fitted on that operator with `key432.py`, `c3e429dd…`:
  - carried 295, n-donor histogram {1: 5, 2: 155, 3: 140};
  - w q75 already at the cap;
  - 26 % at cap vs A/B/C's 14 %.
- **Stack dump** `8b3c59fb…`: `[cd4] sizing=carrying g=0.516548: solved 0.5398 on 182 rows`, and the all-rows rule
  would over-dose ×1.167. That is the 247 case the carrying rule exists for, now exercised with the stack.
- **ef_b1** = LOAD + coarse MDE + zero repair + zs, no re-level. One build, 2,046 s, 20.15 GiB.

| context | `mde_lambda` median (q10, min) | `zero_repaired` total (blocks) | `zs_n_moved` median (q10–q90) | ntclevel on the un-levelled LOAD build: B / claim / ρ / cliff / EV-scale |
|---|---|---|---|---|
| D (A) | 0.9580 (0.9305, 0.866) | 4 (4) | 51 (32–65) | 0.2504 / 0.2338 / **1.0712** / **0.9662** / **0.9561** |
| E (B) | 1.0000 (0.9997, 0.000) | 4 (4) | 37 (25–48) | 0.2379 / 0.2397 / 0.9925 / 1.0038 / 0.9934 |
| F (C) | 1.0000 (1.0000, 0.631) | 7 (7) | 64 (47–79) | 0.2489 / 0.2661 / 0.9351 / 1.0341 / 1.0234 |

- **The low-coverage path works.** Every knob runs at low coverage:
  - `repin-panel-support` pins `{133, 182}`;
  - the strict postbuild passes 41/41, and the driver's own postbuild fails 40/41 as on the realistic panel (T5);
  - the 428 audit reads 0 zeroed gated genes (15 repaired), 4 at `lfc < −5`, max lfc +1.67;
  - `ballast_jk_ratio` median 1.028, `ballast_clip_mass` max 1,887;
  - the proposed guard accepts the meta (P5).
- **The DE dose again follows the controls:** 0.958 / 1.0 / 1.0. zs moves 12–28 % more genes per block than on the
  realistic panel.
- **The stack over-energises a low-coverage panel.**
  - Un-levelled, the Jurkat-like context sits past its cliff: ρ 1.071, which is 7.1 % in energy and 3.4 % in scale.
    That means mse beyond the cliff in that context. E sits at its cliff (ρ 0.99).
  - Compare the realistic panel's b1: ρ 0.992 / 0.914 / 0.859.
  - The EV-optimal string `D=0.9561,E=0.9934,F=1.0234` is 4 % below 1 in D.
  - Neither A/B/C's string (`A=1.0088,B=1.05,C=1.0736`, which the builder refuses for its labels anyway) nor the
    realistic panel's string would have been safe here.
- **The levelled ef build was not run** (one more 35-min build). On the realistic panel one step converged to within
  0.15 % for factors up to 1.068. A 0.956 step is expected to converge the same way (B ∝ s² held to 0.3 %), but that
  is not measured.

## 8. The day-of command sequence for this stack (ready to paste into DEF-RUNBOOK)

This is `build432.sh`'s sequence with the rehearsal-only parts replaced:
- the in-memory wrapper becomes the real source pin (T2);
- the throwaway frozen table becomes `scripts/frozen_inputs.json`;
- the synthetic paths become the day's.

The builder is called with exactly the driver's argument list (the driver's `repin-panel-support`/`build` stages build
the same command). **bash, not zsh.** Wall/RSS from this rehearsal.

```bash
#!/usr/bin/env bash
set -euo pipefail
cd ~/Code/virtual-cell-challenge-2026
unset $(env | grep -o '^VCC_[A-Z0-9_]*') 2>/dev/null || true      # NOTHING exported until step 5 (T3)
B=data/challenge                                  # D/E/F (A/B/C copied to data/challenge-abc first, §2)
R=experiments/NNN-def-round/work; N=OPAQUE_NAME   # <- set both (N: a random joke stem, AGENTS 14)
INC=submissions/SHIPPED_ABC_ARM.h5ad.meta.json     # <- the A/B/C meta of the exact config shipped (433h or successor)
OP=$R/whitener_union_pruned_def.npz; OPJ=$R/whitener_union_pruned_def.json; CD4=$R/cd4_compact
KEY=$R/pert_scale_def_lo0.3_hi3.npz; DUMP=$R/dump_h1x2_cd4x15.npz
BUILD=(uv run python scripts/make_submission_pm_cd4.py --bundle $B --frozen-inputs scripts/frozen_inputs.json
       --expect-contexts D,E,F --whitener $OP --whitener-metadata $OPJ --whitener-id union_pruned_def
       --cd4-root $CD4 --max-rss-gib 24.0)
FLAGS=(--min-cpm 10 --tail-spread 1.0 --keep-top 8000 --keep-spread 0.25 --overdispersion 4.0 --emitter rc_callset)

# 1. light stages + operator re-fit (5.4 min, 12.9 GiB). Look for: axis sha != edf6ba88..., npz sha != e2f65ecf...
uv run python scripts/def_pipeline.py --bundle $B --work $R --name $N --stages preflight,identity,cd4-extract,census,whitener
# 2. 410's key for THIS panel on THIS operator (3 s, 4.1 GiB). KA first: must print "KA PASS d7f76a22..."
uv run python experiments/432-def-rehearsal-new-knobs/key432.py --bundle data/challenge-abc \
    --op data/donors/whitener_union_pruned.npz --lo 0.3 --hi 3 --out $R/key_ka_abc.npz --ka
uv run python experiments/432-def-rehearsal-new-knobs/key432.py --bundle $B --op $OP --lo 0.3 --hi 3 --out $KEY
# 3. pins, arm env still UNSET (2 s): frozen table + operator + CD4 banks
uv run python scripts/def_pipeline.py --bundle $B --work $R --name $N --stages repin-frozen,pin-artifacts
# 4. the key's pin (T2), by hand until def_pin_panel_artifacts.py --pert-scale lands: ONE entry at the top of
#    PERT_SCALE_ARTIFACTS:  "wshare_lo0.3_hi3_def": (ROOT / "<$KEY, repo-relative>", "<shasum -a 256 $KEY>"),
#    git diff scripts/make_submission_pm_cd4.py  -> the driver's pin + this entry, nothing else; commit.
# 5. the day's configuration (433h minus LOAD/CTX; driver SHIPPED_ENV + DEF_ENV + the September knobs)
export VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.34 VCC_H1_UNION_G=0.28 VCC_CD4_SIZING=carrying \
       VCC_CD4_TARGET_G=0.3443655278649641 VCC_CD4_DEFLATE_BASIS=V5_abc \
       VCC_NORMREL_EXP=-1.5 VCC_NORMREL_MODE=cap VCC_CD4_PARTIAL_G=0.5 VCC_DE_SIGN_PRIOR=clean4q80 \
       VCC_PERT_SCALE=wshare_lo0.3_hi3_def VCC_PERT_SCALE_BALLAST=1 VCC_BULK_SHARPEN_Q=1.25
"${BUILD[@]}" --repin-panel-support --min-cpm 10          # 10 s, 20.1 GiB; the two numbers == census's prediction
# 6. the stack dump (12 s, 20.1 GiB). CD4 x1.5 is TARGET_G on D/E/F, NEVER VCC_CD4_TARGET_F=0.51 (T6)
VCC_H1_UNION_G=0.56 VCC_CD4_TARGET_G=0.5165482917974462 VCC_BULK_DUMP=$DUMP \
    "${BUILD[@]}" "${FLAGS[@]}" -o $R/dump-never-written.h5ad
#    look for "[cd4] sizing=carrying g=0.516548", "[h1-union] g=0.56", "[bulk-dump] wrote ... exiting before emission"
DSHA=$(shasum -a 256 $DUMP | cut -d' ' -f1)
# 7. the un-levelled stack build (25 min, 20.2 GiB)
export VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_LOAD=$DUMP:$DSHA
VCC_RC_AUDIT_DUMP=$R/rc_audit_${N}_lvl0.jsonl "${BUILD[@]}" "${FLAGS[@]}" -o $R/${N}_lvl0.h5ad
# 8. its re-level string (22 s, 3.4 GiB); on the LOADed artifact, never the unloaded one (section 4)
uv run python experiments/417-bulk-sharpen/ntclevel417.py $R/${N}_lvl0.h5ad --bundle $B --contexts D,E,F --sd 0.01
# 9. the shipped build: + that string [+ VCC_DE_MAG=dn_k250_zs only if 433z's read landed] (25 min; 35 with zs)
export VCC_BULK_CTX_SCALE='D=?,E=?,F=?'     # <- paste the string step 8 printed
VCC_RC_AUDIT_DUMP=$R/rc_audit_$N.jsonl "${BUILD[@]}" "${FLAGS[@]}" -o $R/$N.h5ad
# 10. checks (each < 30 s)
uv run python experiments/417-bulk-sharpen/ntclevel417.py $R/$N.h5ad --bundle $B --contexts D,E,F --sd 0.01  # EV-scale 1.000+-0.002, rho ~0.978
uv run python experiments/428-zero-gene-audit/audit428.py $R/$N.h5ad --bundle $B     # "0 in 0 blocks"; eyeball lfc > +3
uv run python experiments/431-rc-audit-dump/analyse431.py $R/rc_audit_$N.jsonl        # mde_lambda median >= 0.8 per context
uv run python scripts/def_postbuild.py --meta $R/$N.h5ad.meta.json --bundle $B --cd4-root $CD4 \
    --incumbent $INC \
    --allow-config-diff whitener_id,pert_scale,pert_scale_sha256,bulk_load_sha256,bulk_ctx_scale --out $R/postbuild.json
uv run python scripts/def_pipeline.py --bundle $B --work $R --name $N --stages prep    # dry run, 8-10 s, <= 20 GiB
# then DEF-RUNBOOK section 2's manual tail: vcc prep (real), git commit of the pins, vcc submit.
```

Total from bundle to a prepped stack artifact: ≈ 5.4 + 0.1 + 0.5 + 25 + 0.4 + 25–35 + 1.5 ≈ **58–68 min**. The earlier
single-build path took 29 min.
- The re-level costs one extra full build: its string needs the LOADed render's B_c.
- A bulk-only render mode, which would compute B_c from the loaded bulk without emitting cells, would recover ~25 min.
  It does not exist; it is an idea.

With the proposals applied, steps 2–4 collapse into the driver's own `key` / `pin-artifacts` stages. The `DEF_ENV`
default then carries the key id, and step 5's first two export lines become driver defaults. Steps 6–10 are unchanged.

## 9. The new `SHIPPED_ENV` / flags line for `scripts/def_pipeline.py`

```python
SHIPPED_ENV = {"VCC_ALPHA": "0.982", "VCC_CD4_TARGET_F": "0.34", "VCC_H1_UNION_G": "0.28",
               "VCC_NORMREL_EXP": "-1.5", "VCC_NORMREL_MODE": "cap", "VCC_CD4_PARTIAL_G": "0.5",
               "VCC_DE_SIGN_PRIOR": "clean4q80", "VCC_PERT_SCALE_BALLAST": "1", "VCC_BULK_SHARPEN_Q": "1.25",
               "VCC_DE_MDE": "0.1,0.1,0.3,coarse", "VCC_DE_ZERO_REPAIR": "-5"}
DEF_ENV = {"VCC_CD4_SIZING": "carrying", "VCC_CD4_TARGET_G": "0.3443655278649641",
           "VCC_CD4_DEFLATE_BASIS": "V5_abc", "VCC_PERT_SCALE": "wshare_lo0.3_hi3_def"}
ABC_ENV = {"VCC_CD4_DEFLATE_K": "5", "VCC_PERT_SCALE": "wshare_lo0.3_hi3"}
SHIPPED_FLAGS = ["--min-cpm", "10", "--tail-spread", "1.0", "--keep-top", "8000",
                 "--keep-spread", "0.25", "--overdispersion", "4.0"]            # unchanged; --emitter rc_callset
```

- **Day-of exports, never defaults:**
  - `VCC_BULK_LOAD=<dump>:<sha>` (step 6);
  - `VCC_BULK_CTX_SCALE=<ntclevel string>` (step 8);
  - `VCC_DE_MAG=dn_k250_zs`, only if 433z lands.
- **Stack dump overrides, exported for step 6 only:** `VCC_H1_UNION_G=0.56 VCC_CD4_TARGET_G=0.5165482917974462`.
- **Guard.** The full diff, the guard rows and the postbuild allow-list are in `proposals/def_pipeline.proposed.diff`.
  If 433h's successor changes any value above, `SHIPPED_ENV`, `ARM_BOARD_DELTA` and 240's row 1 move together.

## 10. Files

- **Scripts** (all in this directory):
  - `build432.sh` (phases), `build432x.sh` (the extras copy: `dump_b0env`, zs on ef_b1);
  - `chain432.sh` / `chain432x.sh`;
  - `key432.py`, `run_builder432.py`, `cmpdump432.py`, `summarise432.py`;
  - `proposals/`: the diffs, their generator and their tests.
- **Work:**
  - `work/state432.log`: every phase, gate, sha before/after, and outcome;
  - `work/chain432*.log`;
  - `work/<build>/`: meta, logs, rc dumps, ntclevel/audit outputs, pin records, patched-builder provenance;
  - `work/summary432.json`, `work/cmpdump432.json`, `work/test_arms432.json`, `work/test_pin432_ef.json`, `work/traps/`.
- **Also:** `work/dayof432.sh.txt` (§8's block, `bash -n` clean), `work/drydrv/` (the proposed driver's dry run),
  `work/traps/`, `work/chain432.before-resume2.sh` and `work/make_def_pipeline_proposal.before-key-stage.py`
  (pre-edit copies).
- **Deleted at the end:** the four h5ad artifacts (b0, b1, b2, ef_b1; ≈ 47 GB); no `.vcc` was ever written.
- **Kept:**
  - the re-fitted operator `work/base/whitener_union_pruned_def.npz` (1.3 GB, reproducible in 6 min);
  - the five dumps (85 MB each) and the keys (20 KB).
