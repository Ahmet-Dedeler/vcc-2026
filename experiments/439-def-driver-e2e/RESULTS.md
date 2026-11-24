# 439 — the installed October driver runs end to end on an October-like panel with nothing exported

Status: 2026-09-24 03:32Z, PI. `run439.sh`; the driver log is `work/efr/driver.log`; the artifact was deleted after
the audit (the meta is kept).

**What ran.**
- Components: `scripts/def_pipeline.py` `5463fe28`, pin script `10af47e7`, key script `0e375b77`, builder `63dfd7d0`.
- Panel: 238's essential-free, realistic-coverage synthetic D/E/F bundle (`ef-bundle-r`).
- **Zero `VCC_*` exported.** Every knob came from the driver's `SHIPPED_ENV` + `DEF_ENV`.
- Stages: `preflight, identity, cd4-extract, census, whitener, key, repin-frozen, pin-artifacts, repin-panel-support,
  build, postbuild`, with `--restore-builder-after`.

| stage | wall | peak RSS | rc |
|---|---:|---:|---|
| preflight / identity / cd4-extract / census | 14 s | 2.7 GiB | 0 |
| whitener (the day's operator) | 316 s | 12.2 GiB | 0 |
| **key** (`scripts/def_pert_scale_key.py` on the day's operator) | 3 s | 4.1 GiB | 0 |
| repin-frozen (bare env) | 1 s | 0.2 GiB | 0 |
| **pin-artifacts** (operator + CD4 banks + the key via `--pert-scale`) | 0.5 s | 0.06 GiB | 0 |
| repin-panel-support | 10 s | 20.1 GiB | 0 |
| **build** | 1,428 s | 20.2 GiB | 0 |
| postbuild | 1 s | 0.1 GiB | 1 (see below) |
| **restore-builder** (in `finally`) | 0 s | — | 0: `56d608b3 → 63dfd7d0`, identical |

**What it proves.**
- **Every September knob lands from the defaults.** The driver's own arms guard passes 66/66. The written config
  carries:
  - `pert_scale wshare_lo0.3_hi3_def` (the day's key, sha `5d98344f…`, registered by the pin);
  - `de_mde` 0.1 / 0.1 / 0.3 with `coarse: true`;
  - `de_zero_repair -5`, `bulk_sharpen_q 1.25`;
  - `cd4_sizing carrying` with `cd4_target_g 0.3443…` (CD4 ×1, 438's October default), `cd4_deflate_basis V5_abc`;
  - `normrel_exp −1.5`, `de_sign_prior clean4q80`, ballast (`kmax 16`, `lb 995000`), `h1_union_g 0.28`;
  - `expect_contexts D,E,F`.
- **MDE dose:** `mde_lambda` median 1.0 (min 0.14). Zero-gene audit: **0 zeroed gated genes** in 900 blocks.
- **The restore works on failure.** The postbuild stage exited 1, and the builder was still restored byte-identical
  in `finally` (437's T10 fix).
- **The one postbuild FAIL is the expected config-drift line** (40/41). With no `--incumbent-meta`, postbuild compares
  against the September-5 incumbent (cd4svd5), which predates every September knob.
  - **On the day, pass `--incumbent-meta submissions/<the shipped A/B/C arm>.h5ad.meta.json`**, as 432 T4 says. Then
    only the panel-specific keys differ, and those are on the allow-list.
  - Until then, read the 40/41 as the pass condition for a run with the default incumbent.

**Not covered here** (432 covered them with its own wrappers):
- the stack LOAD on D/E/F (use 438's dump recipe, H1 ×2 / CD4 ×1 by default, then `dayof438.py`);
- the re-level (`ntclevel417.py` on the LOADed artifact, then a rebuild with `VCC_BULK_CTX_SCALE`);
- zs;
- `vcc prep`.
