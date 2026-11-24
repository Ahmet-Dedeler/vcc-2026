## 3. Ranked tables (all vs 417a's bulk = PP + 410 key clip(0.3, 3) + sharpening q 1.25)

**Conventions**
- Every Δoverall is (2.2138 Δpds_bd − 1.051 Δraw mse)/6 on the primary route (r417a: ballast-realised, KMAX 16, anchored
  at 417a's measured energy). Values are per fold in the order k562_gwps / orion_hct116 / orion_hek293t, with the
  cancer-fold mean and n+ = folds positive.
- "418 route" is 418's design-energy route. On 418's own rows it reproduces 418's table exactly (KA-A).
- "A/B/C energy vs 417a" is the realised render energy per context. Flags: **HI** > 1.03 and **LO** < 0.90 in any
  context.
- The builder column lists only what changes vs 417a's env: `VCC_BULK_SHARPEN_Q` 1.25, `VCC_PERT_SCALE=wshare_lo0.3_hi3`,
  PP composition.
  - `LOAD[…; PS …]` means `VCC_BULK_LOAD=<dump.npz>:<sha256>`, with a dump made by `VCC_BULK_DUMP` under 417a's env plus
    the listed composition knobs (`VCC_H1_UNION_G`, `VCC_CD4_TARGET_F`) and that `VCC_PERT_SCALE`. The key lives inside
    the dump, and sharpening applies at render time, so dumps are q-agnostic.
  - `CTX=` means `VCC_BULK_CTX_SCALE=`.
  - `*` means an artifact that is not wired into the builder yet. `lodo415_lo0.5_hi3` = `work/pert_scale_415pp_lo0.5_hi3.npz`,
    sha256 `8c82bef37d0a…` (full sha in `work/pert_scale_415pp_lo0.5_hi3.json` and `analyse420.json`). It is 415's key
    vector re-clipped at (0.5, 3), i.e. w = max(w_lo0.3, 0.5).

The full 150-design × 4-fold grid (6 compositions × 4 keys × 6 q, plus own-key and control rows) is in
`work/analyse420.log` ("ALL DESIGNS") and `work/analyse420.json`.

### 3.1 The one-variable steps from 417a (the ranking that decides the next arm)

{{T-one}}

Reading T-one:
- **Class-calibrated.** Composition steps use 418's history constants. The stack uses H1 then CD4, each with its own
  family constants. Key steps use k 1.407. Sharpening uses the kolf2-like transfer, which is pessimistic.
- **The composition steps carry a real board-emulated pds cost.** H1 ×2: P(Δpds_bd > 0) = 0.06 / 0.02 on the Orion
  folds. It is paid for by mse. kolf2 reads the H1 arms 5–10× higher; that fold is biased (iPSC truth vs the hESC source,
  418 §3).
- **Keys.** The floor-0.5 variants gain only through the emulated pds (literal +0.0001 / +0.0002).

{{T-a}}

{{T-b}}

{{T-c}}

### 3.2 (d) Best combinations

Top of the whole grid (no re-levelling). The q ≥ 1.75 rows are all **LO**, so read them with §4 and §6.

{{T-d}}

The same ranking restricted to q ≤ 1.5 with the energy inside [0.90, 1.03] of 417a:

{{T-d2}}

With re-levelling on top: each design at its own EV-optimal `VCC_BULK_CTX_SCALE`, with the pds and cos of the scaled design
measured on PLOLO (realised, K 16):

{{T-lev}}

- **The level-only gain is nearly identical across folds by construction.** It is an A/B/C-energy computation; only X
  differs by fold. So its 3/3 is not independent evidence. Its risk sits in 187's cliff positions (claim_c), not in the
  folds.
- **The scaled designs cost a little direction.** Realised cos falls 0.3–0.9 % and pds_bd −0.0002 to −0.0007 at
  factors ≥ 1.09. The render is nonlinear and more blocks cap. Both costs are already inside the numbers.
- **q 2 cannot be fully levelled.** Its factors hit the builder's 1.2 cap.

**417b re-levelled from its MEASURED artifact energy** (`cliffcheck414.json`, "my-thesis-advisor-is-a-cat" =
`submissions/417b-…-sharpen150`; ρ = 0.916 / 0.847 / 0.813):
- the EV-optimal factors are `VCC_BULK_CTX_SCALE=A=1.045,B=1.092,C=1.115`, at k = 1 and at k = 1.407;
- worth −0.0264 raw mse, **+0.0046 overall** (the mse part; PLOLO puts the pds change at +0.0002);
- these, not the design-route factors, are the ones to build with.

## 4. (b) Sharpening dose: ballast realisability, KMAX 24 / 32, and the DE channel

{{T-kmax}}

**Bulk side.**
- Raising KMAX recovers part of each dose's realisation loss: +0.0005 / +0.0008 at q 1.25 and +0.0009 / +0.0014 at q 1.5
  (K 24 / K 32).
- It cannot recover more. Among capped blocks, most have r_min = 0: the sharpened (or keyed) design sends a gene to zero
  that the carriers express at > 1e-9. No finite K realises that; the builder clips the gene at the carriers' total.
- The capped share falls only 20.1 → 17.2 → 15.8 % (q 1.25) and 28.8 → 26.1 → 24.3 % (q 1.5).

**DE side (the carrier count drops to 376 / 368 on capped blocks only).** Two channels:
1. **Arithmetic mean CPM** (nmae's log2FC, fid's direction).
   - Here the ballast leaks L_K (p_realised − p_de) into the mean. The leak for 417a's design is 0.0318 of 410's full
     weld at K 16, 0.0346 at K 24 and 0.0371 at K 32. For q 1.5 it is 0.0336 / 0.0371 / 0.0403.
   - 410b sat at 0.0333, and the board moved its DE members **+0.0002** overall vs 395.
   - At the weld's measured −0.0171 overall, a linear bound puts the extra leak of KMAX 32 at ≤ −0.0001 (q ≤ 1.5).
   - Re-levelling barely moves the leak: 417a re-levelled 0.0328, q 1.5 re-levelled 0.0346, q 1.75 re-levelled 0.0413
     (`abc420_leak.json`).
2. **Wilcoxon ranks** (fid/jac/reach calls).
   - The ballast share of the rank sample rises from mean K/400 = 1.4 % (q 1.25, K 16) to 2.0 % (K 32). On capped blocks
     it goes from 4 % to 8 % of cells.
   - The emitter sizes its call set on 400 − K carriers, so the carriers themselves adapt. What 8 % of deep,
     bulk-composition cells do to the U statistic of marginal genes is **not bounded here**.

**Verdict.** KMAX 32 is worth ≤ +0.001 on the bulk at the doses that matter, against a mean-channel cost ≤ 0.0001. It
needs a builder constant change and an H1 check of the rank channel (416/419 territory). It is a rider for a q ≥ 1.5 arm,
not an arm.

## 5. Stacking arithmetic: what adds and what doesn't

{{T-int}}

- **The key and the composition are orthogonal to the dose.**
  - The V1 key's gain grows with q: +0.0019 → +0.0024.
  - The stack's gain is flat to q 1.5 (+0.0037 → +0.0038) and falls only at q ≥ 1.75.
  - The CD4 half fades with dose (+0.0017 → +0.0013 at q 2). The H1 half does not.
- **The stack is 87 % additive at q 1.25.** Its CD4 increment on top of H1 ×2 is +0.0012 (3/3).
- **The key composes with the stack at ≈ 95 %.** Stack + V1_5_3 is +0.0064, against +0.0037 + +0.0031.
- **Along any path, cumulative step prices match the direct end-vs-417a price to ±0.0002** (the PATHS block of
  `analyse420.log`). Stacks here price as sums, with the stated redundancy.

## 6. Recommendation

### 6.1 The single next bulk arm, one variable from 417a

**The stack via `VCC_BULK_LOAD`:**

    417a's env + VCC_BULK_LOAD=experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15.npz:<its sha256>

- **Dump.** 417's overnight chain makes it with `VCC_H1_UNION_G=0.56 VCC_CD4_TARGET_F=0.51 VCC_PERT_SCALE=wshare_lo0.3_hi3`.
- **Gate.** 417's KA-LOAD (same-env dump reproduces 417a's X) must PASS first.
- **Config diff vs 417a:** exactly `bulk_load`.
- **Predicted vs 417a:**
  - Δoverall **+0.0037 (k = 1) / +0.0056 (k = 1.407) / +0.0044 (history-calibrated)**, 3/3 cancer folds;
  - Δraw mse −0.026 … −0.037, Δpds_bd −0.0025;
  - render energy 0.995 / 0.991 / 0.992 of 417a;
  - realisation: 18.6 % of blocks at K 16 (417a: 20.1 %), cos q1 0.991;
  - DE members unchanged by construction; leak 0.0353 of the weld against 417a's 0.0318.
- **Why this arm and not the others:**
  - it is the largest step whose value does not depend on the sharpening transfer (unmeasured until 417a reads);
  - it needs no builder edit;
  - its class has two board-calibrated steps (418 §5).
- **Where it could disappoint:** hek293t (+0.0009) and a significant pds cost on the Orion folds.

**Exception, decided by 417a's board row:** if the bulk members imply k_sh ≥ 0.75, send 417b (q 1.5, built) first. That
means raw mse ≤ 0.703 or raw pds ≥ 0.863 (T-read below), and then q 1.5 is worth +0.0037 … +0.0091 vs 417a. Between
k_sh 0.55 and 0.75 the two arms are tied (q 1.5 +0.0021 … +0.0037); below 0.55 the stack wins.

### 6.2 The next three bulk arms (each one variable from the previous)

{{T-seq}}

- **Default (WEAK-alt): stack → + q 1.5 → + re-levelling.** One env key per step: `bulk_load`, then `bulk_sharpen_q`,
  then `bulk_ctx_scale`.
  - Predicted steps +0.0037 / +0.0063 / +0.0049 at k = 1.
  - **Cumulative +0.0148 (k = 1), +0.0182 (k = 1.407), +0.0091 at kolf2-like sharpening.**
  - At zero sharpening transfer, the path still ends ≈ +0.006: the re-levelling refunds q's amplitude cost.
  - Re-derive step 3's factors from step 2's measured artifact energy with the same EV rule. The design route gives
    `A=1.047,B=1.095,C=1.119` for this bulk; for 417b's own bulk the measured route gives `A=1.045,B=1.092,C=1.115`.
- **If 417a reads k_sh ≥ 0.75 (STRONG): 417b → + re-levelling → + stack.** It reaches the same bulk and the same totals,
  collecting the transfer-dependent step first.
- **Fourth arm, either branch: the key.**
  - LODO V1 at floor 0.5 adds +0.0033 … +0.0042 on the stacked q 1.5 bulk. It needs the new artifact wired in and a new
    dump.
  - `lodo415_lo0.3_hi3` (already wired) adds +0.0021 … +0.0029, and does not lean on the β emulation.
- **Not recommended now.**
  - q ≥ 1.75: LO energy, factors ≥ 1.10 to relevel, and the least-known transfer. Take it only after 417b's read confirms
    k_sh.
  - KMAX 24/32 as a standalone arm.
  - The own-composition LODO key: it adds ≤ +0.0004 over the fixed PP key and would need a new artifact per composition.

### 6.3 Reading 417a's board row (for the PI at ~00:30Z)

{{T-read}}

- **The pds read is the cleaner one.** k_pds = (raw pds − 0.851156) / 0.01713 does not go through the cliff model, and
  k_sh ≈ k_pds / 0.895.
- **The mse read goes through 187's claims.** k_sh = (solve_X(417a measured A/B, raw mse) / (X_410b × 0.9798) − 1) / 0.1013.
- **If the two reads disagree by more than ~0.3 in k_sh,** the cliff model (and therefore every re-levelling number here)
  is what failed. Re-price the levelled rows before building a CTX arm.

## 7. Limits

- **Sharpening transfer is unmeasured.** Every q row above k = 1 / 1.407 is an assumption. 414's kolf2 reads it at 0.37
  (rel) and 0.33 (pds) of the cancer folds. 417a is the first board point.
- **Composition and key transfer rest on thin calibration.** 418's history constants come from one board step per family.
  k_key 1.407 comes from 410's single key step. FINDINGS H.5 says to treat a once-measured constant as a lower bound on
  magnitude.
- **Re-levelling is a model prediction.** It rests on 187's claim_c (retrodicted 410/410b to +0.0002 raw) and on the
  realised-energy map. It is A/B/C-specific; nothing here calibrates D/E/F cliffs.
- **The energy route has a small residual.** The realised route over-predicts sharpening's energy drop by 0.2–0.9 %
  (417a/410b, 417b/417a, 417b/410b). The error is in the conservative direction for the unlevelled q rows (it
  overstates their amplitude cost). For the levelled rows the factors are ~0.3 % high: use the measured-artifact factors
  (§3.2).
- **The β emulation (board-emulated pds) is fitted to one board point** (410b − 395). The floor-0.5 key's value is
  entirely inside it.
- **The ballast emulation is 414's.**
  - It emulates the design profile, not ratio mode on the carriers' delivered composition. It skips the jackknife boost
    and the 999k trim, and uses a constant carrier total.
  - The energy validation above is its end-to-end check.
- **kolf2 cannot vote**, and it over-reads H1 arms. The cancer folds are two-donor folds, blind to the normrel family
  (418 §6); normrel was not varied here.
- **Parallelism.** Up to three fold processes (each one fold, ≤ 6.6 GiB) and one A/B/C process ran concurrently. The
  memory gate (≥ 40 % free, ≥ 20 GiB free+inactive) was checked before each and never went below 67 % / 22.7 GiB.

## 8. Files and reproduction

```
cd experiments/420-bulk-stack            # everything writes only into ./work (PYTHONDONTWRITEBYTECODE=1 everywhere)
uv run python abc420.py                  # A/B/C: KAs, energy-route validation rows, 146 designs x {design, realised K16/24/32}
bash chain420.sh k562_gwps orion_hct116 orion_hek293t kolf2   # comp420.py <fold>: KAs + the 150-design grid (+own-key, control)
uv run python analyse420.py              # prices, KA-A, validation, tables, re-levelling factors -> work/level_specs.txt
uv run python abc420.py --level <designs...>   # realised energy maps at f = 0.94..1.20 -> work/abc420_level.json
bash chain420L.sh k562_gwps orion_hct116 orion_hek293t kolf2  # comp420.py <fold> --level <specs>: scaled designs' pds/cos
uv run python abc420.py --leak <designs...>    # DE leak in units of 410's weld -> work/abc420_leak.json
uv run python analyse420.py && uv run python table420.py > work/table420.md && uv run python assemble420.py
```

- Libraries: `lib420.py` has the KMAX-parameterised `realise_k`, the realised rows, the key artifacts and the new LODO
  (0.5, 3) artifact. 418's `lib418.py` / `analyse418.py` are imported read-only.
- Outputs:
  - `work/comp420-<fold>.json` / `arrays420-<fold>.npz` (per-target arrays for the bootstrap);
  - `work/comp420L-<fold>.json`;
  - `work/abc420.json`, `abc420_level.json`, `abc420_leak.json`;
  - `work/analyse420.{json,log}`, `work/table420.md`;
  - `work/pert_scale_415pp_lo0.5_hi3.{npz,json}` (candidate artifact).
- Peak RSS (`/usr/bin/time -l`, max resident):
  - comp420: 6.65 / 7.00 / 7.04 / 5.26 GB (k562 / hct116 / hek293t / kolf2);
  - comp420 --level: ≤ 6.56 GiB (internal peak);
  - abc420: 6.09 GB; --level / --leak ≤ 5.43 GB.
- **Writes.**
  - Nothing outside this directory. Fold caches are APFS clones of 418's in `work/plolo/`.
  - `find -newer` over every imported directory lists only other agents' files: the builder and its `.pyc` from 417's
    builds, and 414's `cliffcheck414.json`, which 417 updated with 417b.
