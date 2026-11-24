# 811 — Lane B2: a cross-dataset model, always leave-one-line-out. Inventory, baseline, and M811 (GO run)

## 00000000000. X2 banks fix (coordinator SIXA 344c4c31 items 1–2; written 2026-10-22T08:09Z by the clock)

**Why.** Readiness's preflight of dumphook811c0.sh on r0's stack dump refused: 136 of D/E/F's targets were not in the pinned banks
(37e34335: ABC + L1–L6). The DEFo/DEF composites were priced on EC2 with kit811ext's banks_ext, which was never brought to the Mac.

**What changed.**
- One code change: `hook811.py` BANKS_MANIFEST_SHA 37e34335… → **4cc6d82a9e09f5cf6a0edad28a3f704e33db3af3e4ea06dca0ecf8ed12e7bba2**. The file is now
  c7cff26d… and the diff is exactly one line (`work/defbanks/hook811.diff`).
- New data: `work/banks` = kit811ext.py (unchanged, KIT811_PANELS=DEF), D/E/F's rows appended to the pinned rows. Old rows are byte-identical
  (asserted) and every `measured` mask is unchanged.
  - The 18 npz files are copied byte for byte by the new adapter `defbanks811.py` into ext811's manifest schema, with the DEF panel added.
  - The old set is kept at `work/banks.pre-def-37e34335`. Manifest copy: `work/defbanks/banks.manifest.4cc6d82a.json`.
- Runs: `20261022T074226Z-defbanks811` ($0.18) and `20261022T075850Z-kaabc811b` ($0.03).

**Known answers.**
1. KA-same-C0 (kasame811c0, unchanged) PASS on 10 folds (DEFo ×3, DEF ×3, ABCp ×3, L1): the re-pinned hook's transform == composite811's C0.
2. KA-priced PASS: composite811 re-run with the new banks reproduces the priced DEFo/DEF arrays exactly (402 arrays, 0 differ).
3. A/B/C hook, same platform (EC2), PASS: old banks + old hook vs new banks + re-pinned hook on dump 2dfe327d at λ 0.25 are byte-identical
   (both ba7a6477…). Against the Mac's 10-20 file (96204332…) both differ by ≤ 1.3e-15 in `bulk`: platform arithmetic. λ 0 == input content.
4. **Mac file-sha KA PASS** (2026-10-22T08:11Z; run beside X1 on the coordinator's OK, 8.9 s, peak RSS 0.75 GB): the re-pinned hook (c7cff26d) with the
   new banks (4cc6d82a) at λ 0.25 on 2dfe327d writes 9620433220459420…, byte for byte the 10-20 file. λ 0 content == input.
   `kamacdef811.sh`, `work/kamac/def1022/summary.txt`.

## 0000000000. THE LZ SIGN CHECK (sealed reserve, one access; written 2026-10-22T07:48Z by the clock): **PASS**

- Claim `experiments/800-checker/LZ-CLAIM-811.md` (3a8f7252…, ledger e9b4bfd6). The deciding increment is C0 − B0, fixed from DEFo before the access
  (Amendment 4b; decide6a).
- Override: `experiments/800-checker/lzaccess811.sh` (e990d562…). No guard or engine was edited; tags LZ7 = L7, LZ8 = L8; A/B/C frames.
- Kits 20261022T070218Z / 070608Z (18 jobs rc 0). Composites 20261022T072630Z / 072928Z: every KA passed before any score was read.
- Aggregation: agglz811.py (327b7d02, committed before any LZ score), 4,000 draws, seed 800, N = 1, 6 folds pooled.
- **C0 − B0 (deciding): +0.00999 [+0.00659, +0.01273], 6/6 folds > 0, LT3 +0.00571** (LZ7 K562 / HCT116 / HEK293T, LZ8 K562 / HCT116 / HEK293T: +0.00929 / +0.00912 / +0.01651 / +0.00611 / +0.00811 / +0.01080; pdsl +0.01044). C1 − B1 (reported): +0.00645 [+0.00373, +0.00882], 6/6.
- Essential-screen counts: L7 9, L8 20. As DESIGN predicts for arc800 draws, LZ reads 811 ~2.2× DEFo.
- Record: `experiments/800-checker/ACCESS-LZ-811.md`. With runbook §5 items 1–5 all holding, the 6a dossier goes to the coordinator.

## 000000000. THE D/E/F READ (runbook §5, swap-ON branch per 708; written 2026-10-22T07:01Z by the clock)

**Inputs.**
- Kits: DEFo `20261022T061527Z-lab801-kit-defo` (D/E/F's own control frames, shas in LAB-LOG) and DEF `20261022T061958Z-lab801-kit-def`
  (A/B/C frames). All 18 jobs rc 0. The operators are identical across the two tags.
- Composite runs: `20261022T064006Z-lab811c-kit-defo`, `20261022T064311Z-lab811c-kit-def`. composite811 6ec25788, m811 e00192e6.
- Aggregation: `20261022T065734Z-aggtag811-defo`, `20261022T065756Z-aggtag811-def`. aggtag811 af68906a, agg811 d34d1e26, 4,000 draws, seed 800, N = 1.
  These are §000000's exact commands with C1:C0 added.
- Arrays + jsons: `work/run-comp811kit-DEFo/`, `work/run-comp811kit-DEF/`.

**Known answers: all PASS on both tags.**
- KA-705 and KA-707 hold on every tag fold (worst |d| 0.0), on ABCp/k562 (≤ 1.1e-13) and on L1 (0.0).
- KA-811 holds on L1 (0.0), and KA-L0 is 10/10 per tag. L1's CD4 compact equals 801's.

**The read** (literal pds unless named; per fold K562 / HCT116 / HEK293T):

| pair | tag | overall [95 %] | folds > 0 | per fold | LT3 | pdsl |
|---|---|---|---:|---|---:|---:|
| C1 − B1 | DEFo | **+0.00299** [+0.00054, +0.00566] | 3/3 | +0.00351 / +0.00195 / +0.00351 | +0.00142 | +0.00297 |
| C0 − B0 | DEFo | **+0.00446** [+0.00175, +0.00733] | 3/3 | +0.00446 / +0.00371 / +0.00521 | +0.00291 | +0.00442 |
| C1 − C0 | DEFo | **+0.00033** [-0.00071, +0.00143] | 3/3 | +0.00048 / +0.00035 / +0.00015 | -0.00036 | +0.00050 |
| C1 − B0 | DEFo | **+0.00488** [+0.00158, +0.00840] | 3/3 | +0.00496 / +0.00418 / +0.00550 | +0.00299 | +0.00501 |
| C1 − B1 | DEF | **+0.00222** [-0.00024, +0.00488] | 3/3 | +0.00328 / +0.00171 / +0.00167 | +0.00069 | +0.00256 |
| C0 − B0 | DEF | **+0.00363** [+0.00099, +0.00647] | 3/3 | +0.00429 / +0.00333 / +0.00328 | +0.00212 | +0.00397 |
| C1 − C0 | DEF | **+0.00027** [-0.00079, +0.00135] | 2/3 | +0.00055 / +0.00039 / -0.00014 | -0.00047 | +0.00043 |
| C1 − B0 | DEF | **+0.00398** [+0.00075, +0.00743] | 3/3 | +0.00486 / +0.00384 / +0.00324 | +0.00216 | +0.00448 |

**§5 items, swap-ON branch (708: swap ON, rider ON, weights h3c1.5):**
1. DEFo C1 − B1 has its 2.5 % bound at +0.00054 > 0, with 3/3 folds and LT3 +0.00142. **PASS.**
2. DEF C1 − B1 point is +0.00222 > 0. **PASS.**
3. No Arc-drawn panel reads HARM: ABCp and LABC points are > 0, and the A/B/C board reads were C1 ABOVE and C0 ABOVE. **PASS.**
4. K: |dK_rel| vs inc on DEFo stays within 5 % for C0 (max +4.46 %, HCT116) but not for C1 (+6.09 % on HCT116). The DE-H1 read stands (803's hook:
   no DE cost), and both A/B/C board reads, which carry all six members, came out ABOVE. **Clean: YES.**
- **Amendment 4b:** C1 − C0 on DEFo is +0.00033 [−0.00071, +0.00143], with the bound ≤ 0. The swap does not clear its bound on top of 811.
  **X2 = C0** (R0 + 811, no swap): C0 − B0 passes its own bars, DEFo +0.00446 [+0.00175, +0.00733] 3/3 with LT3 +0.00291, and DEF +0.00363.
- decide6a.py (read-only, `--lz NONE`): swap ON, X2 construction **C0**, hook `dumphook811c0.sh`, LZ increment **C0 − B0**. X2 is not eligible until LZ.
- **Projection (reported, not decided on):** C0 − B0 DEFo +0.00446 ÷ ρ. The preregistered ρ ∈ [1.44, 1.72] gives board +0.0026 … +0.0031.
  ka800's current bulk-direction|ABCp ρ 0.660 (n 3) gives +0.0068.
- **Note:** on D/E/F's panel 811 reads larger than on A/B/C's (C0 − B0 ABCp +0.00269 → DEFo +0.00446). As on ABCp, the swap's residual on top of
  811 is small and not significant.
- Next: the sealed LZ access (LZ-CLAIM-811.md 3a8f7252…, ledger e9b4bfd6), deciding increment C0 − B0.

## 00000000. C1:C0 for SIXA Amendment 4: format test (written 2026-10-20T14:21Z by the clock; no D/E/F data exists)

**Why.** The coordinator's decide6a.py (`670-coord-1015`, b2bd8ba9) reads `"C1|C0|DEFo|pds"` and its `per_fold[].overall` values on a swap-ON day
(LAB-DEF-RUNBOOK Amendment 4). §000000's exact 10-22 loop priced only C1:B1, C1:B0 and C0:B0. C1:C0 is now in the loop (same arms, no new build).

**Format test (PASS).**
- Inputs: `work/run-comp811kit-LABCk/` (the LABCk dry run, SHA256SUMS `b4b14322…`), byte-cloned into the scratchpad with the LABCk folds relabelled
  DEFo and DEF.
- Run: agg811.py (`d34d1e26…`, unchanged) through aggtag811's three lines, with agg811's `HERE` pointed at the scratchpad so the test's wear rows stay out
  of `work/wear811.jsonl`. Pairs C1:B1, C1:B0, C0:B0, C1:C0; `--boot 400` (the day uses 4,000); gate602; peak RSS ≈ 12 MB.
- The json carries all eight keys per tag (`<cand>|<base>|<TAG>|pds` / `pdsl`).
- `decide6a.load_agg` + `c1_over_c0` + `bars` parse them: "C1:C0 DEFo point +0.00023, folds+ 2/3 → swap stays". C1:B1 and C0:B0 reproduce §000000's LABC
  points (+0.00142, +0.00235).
- The full `decide6a.py` run on the 10-16 rehearsal's decision708 (swap OFF) completes and prints X2 "C0:B0 … FAIL", as LABC's bound says.

**By-product, an information read (decides nothing; LABC is worn, n_tries unchanged at 12 since C1 and C0 are not new arms).**
- C1 − C0 on LABC: **+0.00023**, 2/3 folds > 0 (K562 / HCT116 / HEK293T +0.00027 / +0.00117 / -0.00074). Bootstrap interval at 400 draws: [-0.00052, +0.00108].
- The board's exact bulk members gave C1 − C0 = −0.0003 on A/B/C. So LABC, like ABCp (+0.0005 to +0.0006 by difference), reads the swap's
  residual as weakly positive, and the board read it as ≈ 0 or negative. Amendment 4's 2/3-folds rule on DEFo would have kept the swap on LABC.
  That is the over-read pattern the prospective grades showed for composition/direction changes. Amendment 4's bar is the coordinator's; the lab
  only reports.
- **2026-10-20T15:20Z:** the coordinator adopted the tighter bar as SIXA Amendment 4b (1e9a9695): C1:C0 DEFo ci95[0] > 0, else C0. Re-run on this
  test's jsons, decide6a's updated `c1_over_c0` returns False on LABC (ci95[0] −0.00052), i.e. C0, as the board would pick.

## 0000000. C1 board-read hook (built 2026-10-18 08:32–08:57Z by the clock): **READY, with one caveat (below)**. UNAPPLIED.

**What it is.** For production's offered board read of C1 on A/B/C: R-ABC2H's recipe (R0 + 707's swap via a dump hook) plus
c0_l25. It uses 690g's `DEF_DAY_DUMP_HOOK` interface, the same as `experiments/707-strong-tier-direction/dumphook707.sh`.
Nothing in production or in 707 is edited.

| file | sha256 | role |
|---|---|---|
| `dumphook811.sh` | `6cccd69a0651faa100d4aaa12bd68054ff7bd8e3a67d6f5e638221ee84ef4378` | the hook: `<bundle> <work> <dump.npz> <identity.json> <hookdir>` |
| `hook811.py` | `c5a5614ad74c3dab8d975351ed5a4a808fd56de8831c2e0506885a3dc3772f16` | step 2: c0_l25's blend on the swapped dump |
| `kasame811.py` | `44960229508a27a732af1fa700df125488e01036ac33e575a856090b31bccc6b` | KA-same-construction / KA-scale on engine folds |
| `cloud/ka-hook811.sh` | `44c1214b99bbc7fa1a79a03ab0100c806347df3fdd5db7c827f0aa3c4ea6e236` | the KA driver (EC2) |

**What the hook does, in order.**
1. **707's swap first.** It calls `swap707.py --beta 0.25` unchanged, exactly as dumphook707.sh does, and writes
   `<hookdir>/swap707/<name>`.
2. **Then 811's blend.** `hook811.py` blends the swapped dump's bulk rows and writes `<hookdir>/<name>`. Carriers, perts,
   genes and env are copied unchanged, in swap707's key order and dtypes.
   - The construction is the one composite811's C1 applies on engine folds. base811's `bank_rows`, `unit`, `consensus` and
     `blend` and its GROUPS / COMPOSITION literals are exec'd by ast from the pinned base811.py (`c29f9f97…`), the way 803's
     c0811h does: importing base811 would pull in engine801's write guard.
   - m811's day rule is restated (DAY = GROUPS minus KOLF2) and asserted against m811.py's text (`e00192e6…`).
   - The blend is λ 0.25, centring off, norm-preserving on the non-panel genes. Targets without U, and the panel columns,
     are kept byte for byte.
   - **Safeguards:**
     - 811's banks manifest (`37e34335…`) and every bank file's sha256 are checked;
     - the hook refuses any target that is not in the banks' extraction panels (A/B/C + L1–L6), so it cannot run silently
       on an un-extracted D/E/F panel;
     - it asserts that the dump's genes are the official 18,533-gene order.
   - `HOOK811_LAMBDA=0` returns the swapped rows unchanged (the KA-off path).

**How the engine fold maps to the dump (stated explicitly).**
- Same gene axis (official, asserted) and the same effect space (the stack's composition rows).
- The engine blends the unkeyed composition rows; the dump's bulk rows carry the key, and sharpening is a render knob
  applied to the loaded dump. The blend is norm-preserving and depends only on each row's direction, so a positive per-row
  scale commutes with it, and the key's Frobenius constant does not change. **KA-scale measures this at ≤ 3.8e-16 relative.**
- **Held line.** On engine folds the held line's group is out of U. On the board there is no held line, so hold = none and
  U uses all 14 day groups, including the H1 group. The D/E/F identity rule is a no-op here: the one matching bank line,
  Jurkat (context A), carries 0 of A/B/C's targets, and a single dump cannot be per-context anyway.

**Known answers: all PASS.** Run `20261018T084656Z-kahook811b`, r7i.4xlarge, 2 BLAS threads; outputs in
`work/hook811/ka2/` (SHA256SUMS `8f069763…`).

| KA | result |
|---|---|
| **KA-off.** dumphook707.sh run unchanged (its hard-coded repo path supplied by a symlink) vs dumphook811.sh with `HOOK811_LAMBDA=0`, on a `cp -c` clone of A/B/C's dump from 707's KA log (`495-lodo-cap6/work/dump_h1x2_cd4x15_lodo05hi6.npz`, `d58750b8…`) | **byte-identical on the same machine:** dumphook707 `24a88608…` == dumphook811 λ 0 `24a88608…` == its own 707 step `24a88608…` |
| **KA-same-construction.** `hook811.transform` on composite811's B1 rows with the fold's held group out vs composite811's C1 rows, on ABCp/k562_gwps, ABCp/orion_hct116, ABCp/orion_hek293t and L1/k562_gwps | **max \|d\| 0.0 on all 4**; the same targets change |
| **KA-scale.** transform(s·B1 rows) vs s·C1 rows, s a random positive per-row scale | max relative \|d\| 1.9e-16 … 3.8e-16 |
| **690g's validator** (`scripts/def_dumphook_check.py`, read-only) on the λ 0.25 output | PASS: perts, genes, env identical; bulk and carriers finite with the same shape and dtype; mean bulk cosine to the stack dump 0.9752 (min 0.9341) |

- **Caveat: the same-machine byte check ran on EC2.** On EC2, swap707's own output is `24a88608…`; on the Mac, 707's log
  records `77e9de2d…`. Both hooks share that difference: it is swap707's floating-point result on each platform. The Mac
  check was not run, because LAB Amendment 1 forbids lab local jobs while `watch691.sh` / `r0ready692.sh` run (both are
  armed).
  - **Recommended rehearsal for readiness, on the Mac** (seconds), before the build:
    ```
    bash experiments/707-strong-tier-direction/dumphook707.sh data/challenge <w> <dump> <w>/identity.json <w>/h707
    HOOK811_LAMBDA=0 bash experiments/811-lolo-model/dumphook811.sh data/challenge <w> <dump> <w>/identity.json <w>/off
    shasum -a 256 <w>/h707/<name> <w>/off/<name>          # must be identical (and == 77e9de2d on 707's KA dump)
    ```
  - If you want it, the lab lead can run that with the watchers' leave.
- A first KA run (`20261018T083702Z-kahook811`) is **VOID for its KA-off line.** It compared two empty strings, because
  `shasum` does not exist on Amazon Linux. Its other KAs equal the v2 run's.

**What the hook does to A/B/C's dump** (λ 0.25, from the sidecar `work/hook811/ka2/dump_h1x2_cd4x15_lodo05hi6.npz.hook811.json`):
- **Rows moved:** all 300 targets carry U and all 300 rows move, against 100 for the swap alone. Moved-row cosine to its
  input is 0.965 / 0.981 / 0.984 / 0.986 / 0.994 (quantiles).
- **Change energy:** 2.7 % of the bulk.
- **Norm preservation:** rows keep their norm to 9.4e-15 relative; panel columns and unmoved rows are unchanged.
- **Bank coverage of A/B/C's 300 targets** (no held exclusion): lines per target min 2, median 5, max 11.
  - Histogram: 2 → 2, 3 → 10, 4 → 122, 5 → 112, 6 → 39, 7 → 5, 8 → 1, 9 → 2, 10 → 3, 11 → 4 targets.
  - Lines: K562 299 (GWPS + VIPerturb + Jiang), HCT116 290, HEK293T 294, CD4 291, **H1 80** (the H1 group is IN, as on
    the day), TeloHAEC pilot 65, DLD1 67, A549 29, BXPC3 / HAP1 / HT29 / MCF7 9 each; RPE1 and Jurkat 0.
- **Board rule vs fold rule** on the engine rows: the change directions agree at cosine 0.52–0.97 (median 0.79–0.84). The
  board includes the held line's own group, which the folds exclude, so the board's U is not the folds' U.

**Readiness.**
- **Env:** `DEF_DAY_DUMP_HOOK=experiments/811-lolo-model/dumphook811.sh`, together with R-ABC2H's own settings. The hook
  replaces 707's dumphook, because it runs the swap itself; do not chain both.
- **Commit first:** def_day refuses an uncommitted hook (`git ls-files`, a clean `git status`), so the lab lead must
  commit dumphook811.sh and hook811.py.
- **Local inputs:** the hook reads `experiments/811-lolo-model/work/banks/*.npz`, which is git-ignored and local on the Mac
  (430 MB; every file sha-checked).
- **Runtime:** 17 s wall on EC2 (the swap ~5 s, hook811 12.1 s), peak ≈ 1.2 GB.
- **D/E/F:** the hook refuses D/E/F's panel until its banks are extended. That is the kit's job (kit811ext), not this hook's.

**The prediction input** (for the lead's preregistration; board overall units, deflated as the lab does):

| | checker C1 − B1 [95 %] | folds > 0 | ÷ ρ 1.72 | ÷ ρ 1.44 | interval with ρ ∈ [1.44, 1.72] folded |
|---|---|---:|---|---|---|
| ABCp | +0.00153 [−0.00096, +0.00428] | 2/3 (leave-top-3-out −0.00022) | **+0.00089** [−0.00056, +0.00249] | +0.00106 | [−0.00067, +0.00297] |
| LABC | +0.00142 [−0.00104, +0.00412] | 2/3 (leave-top-3-out −0.00013) | +0.00083 [−0.00060, +0.00240] | +0.00099 | [−0.00072, +0.00286] |

Differences between the board read and these folds:
- the board's U includes every line, with no held group;
- the board's base is R-ABC2H's own built bulk.

**Cost.** EC2 $0.28 (kahook811 $0.133, kahook811b $0.143), within the $3. Lane total ≈ $2.54. No 811 instance is running
(`labcloud.sh status` 08:56Z).

## 000000. Composite arms for LAB-DEF-RUNBOOK Amendment 1 (built 2026-10-17 22:35–23:13Z by the clock; ready for 10-22)

**Why.** 6a allows one replacement attempt, and the likely one is R0 plus 707's swap. So c0_l25 is priced as an INCREMENT on
top of what would ship, using the four arms Amendment 1 defines.

**The script.** `composite811.py`, sha256 **`6ec25788276372e824d067b943268ff0e5cbea6b477494d3d53a2b9f7fa28bda`**. `m811.py` is
byte-identical (`e00192e6…`), and so are `base811.py` (`c29f9f97…`), `ka804.py` (`d5d304db…`) and `engine801.py`
(`7442433b…`). On any engine801 fold the script builds:

| arm | definition |
|---|---|
| B0 | inc (pct+P2g) + 705's θ 0.25 rider (P2g re-cut on w705): R0's bulk |
| B1 | B0 + 707's swap (β 0.25, top z_src tertile) |
| C0 | B0 + c0_l25's blend (λ 0.25, day bank list, centring off) |
| C1 | B1 + c0_l25's blend, applied to B1's rows (swap first, then 811) |

All four use B0's keys.

**Where each piece comes from.**
- From ka804 (imported): `zsrc_fold`, `donor_unit_mean`, `p2g_for`, and the constants FLOOR, THETA and BETA.
- The 705 expression lives inline in `ka804.main`, so it cannot be imported. composite811 carries ka804's two lines and
  asserts at run time that they appear verbatim in ka804.py's source.
- The 707 loop and 811's blend are `base811.blend`. Earlier KAs showed it bit-identical to ka804's loop.
- U comes from `m811.fold_U(F, centre=False)`.
- `--mode day` (used for every read) sets inc = pct+P2g and cuts P2g with p2g_oct on every fold, ABCp included.
  `--mode ka804` follows ka804's own per-panel choice and exists only to reproduce its stored arrays.

**Known answers: all PASS** (r7i.8xlarge, 2 BLAS threads). Run `20261017T224004Z-comp811`, 7 folds, every job rc 0.

| KA | folds | worst \|d\| |
|---|---|---|
| KA-705: ka804-mode B0 == ka804's stored 705 | L1/k562_gwps (804's run), ABCp × 3 (801 `ka-abcp`), LABC × 3 (801 `ka-labc`) | L1: 0.0. ABCp / LABC: ≤ 1.7e-13 (references from 801's runs) |
| KA-707: ka804-mode B1 with the rider OFF == ka804's stored 707 | same 7 folds | L1: 0.0. ABCp / LABC: ≤ 2.8e-13 |
| KA-L0: C0 at λ 0 == B0 and C1 at λ 0 == B1, bit-exact | 7 folds | rows identical, priced arrays 0.0 |
| KA-811: C0 with the rider OFF == m811 stage 1's c0_l25 | L1/k562_gwps | 0.0 |

- The KA references are byte copies in `work/refcomp/<fold>/{705,707}.npz` (SHA256SUMS `3d281b54…`).
- **How KA-L0 is exact.** The composite path returns the base rows when λ = 0. The blend arithmetic itself, at λ 0, differs
  from the base rows by at most 8.9e-16 (recorded, not required).
- **One reading of the brief, stated.** "B1 reproduces ka804's 707" can only hold with the rider OFF, because ka804's 707 arm
  carries inc's keys. With the rider ON, B1 = B0's keys plus 707's rows, which is the arm Amendment 1 defines. Its pieces
  are each KA-checked: rows by KA-707, keys by KA-705.

**Information read** (counted in wear, decides nothing): day mode, vs the arm named. agg811 with N = 1, so no shift. Arrays
`work/run-comp811/` (SHA256SUMS `15dd2c49…`); json `work/run-comp811/agg811-811-composite-*.json`.

| pair | surface | overall [95 %] | folds > 0 | per fold K562 / HCT116 / HEK293T | leave-top-3-out | pdsl | r | Δpds |
|---|---|---|---:|---|---:|---|---:|---:|
| **C1 − B1** (811 on top of R0 + swap) | ABCp | **+0.00153 [−0.00096, +0.00428]** | 2/3 | +0.00113 / −0.00001 / +0.00348 | −0.00022 | +0.00193 | +4.49 % | −0.00467 |
| C1 − B1 | LABC | +0.00142 [−0.00104, +0.00412] | 2/3 | +0.00100 / −0.00036 / +0.00362 | −0.00013 | +0.00186 | +5.15 % | −0.00629 |
| C1 − B0 (swap + 811 together) | ABCp | +0.00333 [+0.00036, +0.00666] | 3/3 | +0.00278 / +0.00294 / +0.00426 | +0.00143 | +0.00414 | +7.71 % | −0.00607 |
| C1 − B0 | LABC | +0.00265 [−0.00036, +0.00589] | 3/3 | +0.00207 / +0.00211 / +0.00378 | +0.00092 | +0.00335 | +7.89 % | −0.00831 |
| **C0 − B0** (811 without the swap) | ABCp | **+0.00269 [+0.00005, +0.00564]** | 3/3 | +0.00215 / +0.00124 / +0.00467 | +0.00098 | +0.00335 | +5.92 % | −0.00430 |
| C0 − B0 | LABC | +0.00235 [−0.00019, +0.00522] | 3/3 | +0.00178 / +0.00083 / +0.00445 | +0.00078 | +0.00311 | +6.46 % | −0.00630 |
| B1 − B0 (the swap alone, for scale) | ABCp | +0.00169 [+0.00041, +0.00307] | 3/3 | +0.00161 / +0.00283 / +0.00063 | +0.00075 | +0.00210 | +3.07 % | −0.00140 |
| B1 − B0 | LABC | +0.00113 [−0.00016, +0.00253] | 3/3 | +0.00103 / +0.00232 / +0.00002 | +0.00022 | +0.00139 | +2.60 % | −0.00203 |

**Reading. The expected size for the day.**
- **811's increment over R0 + swap (C1 − B1) is about +0.0015** on both Arc-drawn panels, about 55 % of its size without the
  swap (C0 − B0 ≈ +0.0025–0.0027). The swap and the consensus overlap partly: 707 already moves the strong tier toward the
  sources' mean.
- **Applied to these two panels, Amendment 1's rule would NOT pass on C1 − B1.** The bound is below 0 (ABCp −0.00096, LABC
  −0.00104), 2/3 folds are > 0, and leave-top-3-out is < 0. On HCT116 the increment is ≈ 0.
- **On C0 − B0 (the no-swap branch)** ABCp's bound is +0.00005, with 3/3 folds.
- The day's own read, on D/E/F's panel with N = 1, decides. These numbers decide nothing.
- R0's 705 rider is ≈ 0 on these panels, as recorded (B0 vs inc, r −0.5 … +0.6 %).

**The kit (step 4): copies in 811's dir; 801's files untouched.**
- `kit/run-kit811c.sh` (sha256 `442e867a…`) is a copy of 801's `kit/run-kit811.sh` (`9ccd12f1…`). It does five things:
  1. calls 801's `kit811ext.py` (`2a545547…`) unchanged, to extend the banks to <TAG>'s targets;
  2. copies the kit's ka804 references from `801/work/ka-<tag, lower case>/bulk/` into `work/refcomp/`;
  3. runs `composite811.py --mode day --ka` on L1/k562_gwps and ABCp/k562_gwps (KA-705/707/L0/811);
  4. runs it on <TAG>'s three folds (`--ka` when the kit's references exist, so KA-705/707/L0 gate them too);
  5. writes `work/run-comp811kit-<TAG>/`. Any rc ≠ 0 fails the run.
- `kit/make-jobs811c.sh` (`c5b40246…`) writes `kit/jobs811c-<TAG>.txt` and `kit/sync811c-<TAG>.txt`. It takes 811's
  A-only sync list, drops the A bundles, and adds <TAG>'s bundle, ops, registry and dumps, plus 811's run-m811s1, refcomp
  and 801's `ka-<tag>`.
- **One bug, found and fixed before the dry run:** kit801 names its ka dirs by the lower-case tag (`ka-labck`).
- `aggtag811.py` (`af68906a…`) is agg811 (unchanged) restricted to one tag's folds, so the KA folds sitting in the same run
  dir are not pooled.
- **Dry run on LABCk** (run `20261017T225259Z-lab811c-kit-labck`, $0.199; aggregation `20261017T230951Z-aggtag811-labck`, $0.015):
  - kit811ext's KA PASS: L1's CD4 compact is byte-identical to 801's;
  - the kit's ka804 references were found;
  - all 5 folds rc 0;
  - KA-705/707: 0.0 on the LABCk folds, ≤ 1.1e-13 on ABCp, 0.0 on L1;
  - KA-L0: 10/10 PASS; KA-811: 0.0;
  - **the kit path on LABCk reproduces the direct LABC read exactly:** 15 arm files at worst |d| 0, and aggtag811's pooled
    C1 − B1 / C1 − B0 / C0 − B0 on LABCk equal agg811's LABC numbers to every printed digit.
  - Arrays: `work/run-comp811kit-LABCk/` (SHA256SUMS `b4b14322…`).

**The exact 10-22 commands** (after LAB-DEF-RUNBOOK §2's kit801 + kit801-collect for the tag; window800 CLEAR; repeat with
DEF in place of DEFo):
```
bash experiments/811-lolo-model/kit/make-jobs811c.sh DEFo
bash experiments/800-checker/cloud/labcloud.sh launch --name lab811c-kit-defo --jobs experiments/811-lolo-model/kit/jobs811c-DEFo.txt \
  --type r7i.4xlarge --minutes 45 --parallel 1 $(sed 's/^/--sync /' experiments/811-lolo-model/kit/sync811c-DEFo.txt) \
  --collect experiments/811-lolo-model/work/run-comp811kit-DEFo
bash experiments/800-checker/cloud/labcloud.sh wait <run-id>          # then copy out/repo/.../run-comp811kit-DEFo into 811's work/
R=experiments/811-lolo-model/work/run-comp811kit-DEFo
for p in C1:B1 C1:B0 C0:B0 C1:C0; do echo "OPENBLAS_NUM_THREADS=1 uv run python experiments/811-lolo-model/aggtag811.py DEFo $R --family 811-composite --tries 1 --pairs $p"; done > experiments/811-lolo-model/cloud/jobs-aggtag-DEFo.txt
bash experiments/800-checker/cloud/labcloud.sh stage $R
bash experiments/800-checker/cloud/labcloud.sh launch --name aggtag811-defo --jobs experiments/811-lolo-model/cloud/jobs-aggtag-DEFo.txt \
  --type c7i.2xlarge --minutes 20 --parallel 4 --root-gb 30 --sync $R --sync experiments/811-lolo-model/work/ka800.copy.json --collect $R
```
- **10-20 edit (LAB-DEF-RUNBOOK Amendment 4):** C1:C0 added to the loop and `--parallel` 3 → 4. Format test in §00000000.

**How to read the result** (Amendment 1's rule, from the json):
- On DEFo: the C1:B1 pair's literal-pds `ci95[0] > 0`, AND (`folds_pos == 3` OR `leave_top3_out > 0`).
- On DEF: C1:B1's `point > 0`.
- If the swap is OFF on R0's day path, read the same bars on C0:B0.
- Swap ON (SIXA Amendment 4b, coordinator 1e9a9695): X2 = C1 only if C1:C0's DEFo literal-pds `ci95[0] > 0`; else C0 if C0:B0 passes the same bars, else C1 (decide6a applies this). Amendment 4's earlier bar (point > 0 on ≥ 2/3 folds) is superseded. Pass all eight jsons (four pairs × DEFo, DEF) to `decide6a.py --agg`.
- Gates before reading: the run's KAs pass and rc is 0 everywhere. If D/E/F's panel has no ka804 references, the kit read
  is gated by the L1 and ABCp KA folds only. The kit produces those references, so this is not expected.

**Cost and wear.**
- Cost: EC2 $0.60 for this build (comp811 $0.340, agg811comp $0.045, the LABCk dry run $0.199, aggtag $0.015), within the
  $2 budget. Lane total ≈ $2.26.
- Instances: none of 811's is running; `labcloud.sh status` at 23:13Z lists none.
- Wear: 8 rows appended to `work/wear811.jsonl` (now 78 rows, sha256 `2af7b9e6…`), with `n_tries` 11 on ABCp and 12 on
  LABC, `surface_set` "… composite information read, Amendment 1", and `decides` "nothing".
- The LABCk dry-run aggregation is a known answer (byte-identical inputs), not a new read, so it adds no wear row.

## 00000. A-only instrument read (2026-10-16 12:45–13:16Z; an instrument read, no decision)

**The folds.** 801's 18 `A<k>/<held>` folds: A1–A6 × {k562_gwps, orion_hct116, orion_hek293t}, with registry status
`a-only`. They are arc800 draws with the essential screens removed, built on the L recipe, with inc = pct+P2g.

**Banks had to be extended first (a data input, not a code change).**
- 811's compact banks covered only the A/B/C + L1–L6 target union. The A panels' 1,390 distinct targets were mostly
  outside it, and no CD4 bank for the A panels existed on disk or S3. So every extra line would have been silently
  missing on A folds.
- `ext811a.py` (sha256 `d8b4c205…`) ran on EC2. It used 047's `extract_cd4.py`, with the same arguments as 801's
  `cd4_801.py`, on A1–A6. It then took the other extra banks' rows for the A targets with `ext811.read_bank`, unchanged,
  and merged everything with the old banks into `work/banks_ext/` (manifest sha256 `09255f5f…`).
- Rows added per bank: viperturb +542; CD4 rest / stim8hr / stim48hr +509 / +509 / +506; dld1 +146; telohaec +93;
  h1_pan +56; jiang +9 each; altair +15; h1_vcc +4; kolf2, rpe1 and jurkat +0. Old rows were kept byte-for-byte.
- **Known answers for the extension:**
  - re-extracting L1's CD4 compact gives delta, genes and perts byte-identical to 801's L1 compact, in all 3 conditions;
  - every extra bank's source-file sha256 equals the Mac-built manifest's (`work/banks/manifest.json`);
  - every old bank file's sha256 matches (0 mismatches).
- The 1.7 GB of CD4 compacts stayed in the launcher's run directory
  (`experiments/800-checker/cloud/runs/20261016T125901Z-aonly811/out/…/cd4a`). Only their logs were copied to `work/cd4a/`.
- Cosmetic bug: ext811a's printed "peak 0.00 GiB" is a units slip (Linux reports ru_maxrss in KiB). The job's
  `/usr/bin/time` shows 1.34 GB.

**Code changes, with every sha recorded.**
- `m811.py` is unchanged (`e00192e6…`) and `labc811.py` is unchanged (`3072b668…`).
- `aonly811.py` (`05c056d5…`) is labc811.py with **one detection line changed**. The line
  `assert fold.startswith("LABC/") or fold == KA_FOLD, …` became
  `assert fold[:1] == "A" and fold[1:2].isdigit() or fold == KA_FOLD, …; B.BANKS = HERE / "work" / "banks_ext"`.
  The env-var name also changed (LABC811_OUT → AONLY811_OUT), and so did one docstring path.
- `agga811.py` (`479ab5bd…`) is agg811.py (unchanged, `d34d1e26…`) with one extra surface, "Aonly". agg811's surfaces only
  match fold ids starting "L"/"ABC", and its first aggregation run (`…130915Z-agg811aonly`, $0.033) printed nothing for
  that reason.
- The EC2 driver is `cloud/run-aonly.sh` (`f135d5c0…`).

**Known answers on the run** (r7i.8xlarge, 2 BLAS threads):
- Engine KAs (PT0 cap3 / pct) are 0.0 on all 18 A folds.
- On the L-fold known answer, L1/k562_gwps priced with the EXTENDED banks: inc and src_l25 equal base811's arrays, and
  c0_l25 and c1_l25 equal m811 stage 1's, all at |d| = 0.0 (4/4). So the extension leaves every L-target row intact.
- All 19 fold jobs ended rc 0.
- Arrays: `work/run-m811aonly/bulk/A<k>__<held>/{inc,c0_l25,c1_l25,src_l25}.npz`, listed in
  `work/run-m811aonly/SHA256SUMS.txt` (`12757fb8…`). The KA fold sits in `work/run-m811aonly/ka/`.

**Prices** (18 A folds, vs inc = pct+P2g, literal pds with pdsl beside; N = 1 on this surface set, so there is no shift)

| arm vs base | overall [95 %] | folds > 0 | pdsl overall [95 %] | r | Δpds | Δmse raw |
|---|---|---:|---|---:|---:|---:|
| c0_l25 vs inc | **+0.00638 [+0.00454, +0.00829]** | 18/18 | +0.00769 [+0.00590, +0.00957] | +11.95 % | −0.00603 | −0.04901 |
| c1_l25 vs inc | +0.00659 [+0.00480, +0.00852] | 18/18 | +0.00808 [+0.00627, +0.00995] | +11.92 % | −0.00538 | −0.04887 |
| src_l25 vs inc | +0.00452 [+0.00229, +0.00691] | 11/18 | +0.00580 [+0.00362, +0.00823] | +11.61 % | −0.01050 | −0.04760 |
| **c0_l25 vs src_l25 (many-line increment)** | **+0.00193 [+0.00091, +0.00304]** | **14/18** | +0.00197 [+0.00094, +0.00305] | +0.45 % | +0.00447 | −0.00186 |
| c1_l25 vs src_l25 | +0.00215 [+0.00105, +0.00330] | 14/18 | +0.00235 [+0.00129, +0.00349] | +0.42 % | +0.00512 | −0.00172 |
| c1_l25 vs c0_l25 | +0.00021 [−0.00014, +0.00058] | 9/18 | +0.00038 [+0.00004, +0.00072] (15/18) | −0.04 % | +0.00065 | +0.00015 |

By held line (mean over A1–A6, from `per_fold`):

| held line | c0_l25 vs inc | r | Δpds | src_l25 vs inc | many-line increment | its Δpds | folds > 0 (increment) |
|---|---:|---:|---:|---:|---:|---:|---:|
| K562 | +0.0027 | +4.8 % | −0.0021 | −0.0019 | +0.0046 | +0.0076 | 6/6 |
| HCT116 | +0.0049 | +13.2 % | −0.0125 | +0.0034 | +0.0016 | +0.0041 | 5/6 |
| HEK293T | +0.0115 | +17.9 % | −0.0036 | +0.0121 | −0.0003 | +0.0017 | 3/6 |

**Reading (descriptive, for 800's over-read work):**
- **On A-only panels the many-line increment works through pds, not mse.** The increment is +0.0019, 14/18 folds, with
  Δpds +0.0045 and r only +0.45 %. On L panels it ran mostly through mse (W-tune r +5.6 %, pds +0.006). On these panels
  the sources-only blend already captures nearly all the aggregate-cosine gain (r +11.6 % vs c0's +12.0 %). What the
  extra lines add here is target discrimination.
- **It again comes from the K562 folds**, where the held line's group removal leaves the most non-source lines:
  +0.0046, 6/6. HEK293T gets nothing (−0.0003, 3/6). This matches the per-target R_tgt read in §00.
- **Against inc, the full arm sits between the fresh L panels and A/B/C's panel:** +0.0064 here, against W-tune +0.0085,
  W-conf +0.0068, LABC +0.0023 and ABCp +0.0026. Removing the essential screens does not bring a fresh panel down to
  A/B/C's level for this arm.
- **Centring (c1 − c0):** +0.0002 literal (NULL) and +0.0004 pdsl, with 15/18 folds positive. That is the same small
  size as everywhere else.
- No decision changes: M811's §5 outcome stays WEAK (§000).

**Wear and cost.**
- Wear: 6 rows appended to `work/wear811.jsonl` (now 70 rows, sha256 `82663b4b…`), with `n_tries` 1 and `surface_set`
  "A-only (A1-A6 x 3; instrument read)". The family's count goes up by 1.
- Cost: EC2 $0.355 for this read (aonly811 $0.285, the empty aggregation $0.033, agg811aonly2 $0.037). The lane total is
  **$1.66**.
- Instances: no 811 instance is running (`labcloud.sh status` at 13:16Z lists none).

## 0000. LABC instrument read (2026-10-16 10:57–11:08Z; an instrument read, not a reopening, and no decision)

**The folds.** 801's `LABC/<held>` is A/B/C's 300 targets built with the fresh-fold (L) recipe, under D/E/F rules, with
inc = pct+P2g. Registry entries are status `labc`; the dumps' sha256 are `eba095c1…`, `e97f1fbd…` and `34ecd80e…`.

**How it ran: a wrapper, not an edit.**
- `m811.py` is byte-identical at `e00192e6…`. As frozen, it does recognise `LABC/` (the panel starts with "L"), but it
  then treats LABC as an L fold and requires 804 and base811 reference arrays that do not exist for LABC, so every fold
  would be voided.
- A single detection line cannot fix that. At least two conditions would change, and so would the frozen sha. So I added
  a wrapper instead: `labc811.py` (sha256 `3072b668…`). It adds no arithmetic. It calls `m811.fold_U` (c0, c1),
  `base811.consensus` over the fold's sources (src, exactly as base811 built `src_lX`), `base811.blend`,
  `engine801.named("pct+P2g")` and `engine801.price`.
- **Wrapper KA, on the L fold L1/k562_gwps:** inc and src_l25 equal base811's arrays, and c0_l25 and c1_l25 equal m811
  stage 1's arrays, all at |d| = 0.0 (4/4).
- **Engine KAs:** PT0 cap3 and PT0 pct are 0.0 on all three LABC folds.

**Run and files.**
- EC2 `20261016T105755Z-labc811`: r7i.8xlarge, 2 BLAS threads, 4 jobs, rc 0, actual $0.149 (its launch cap was $0.89 on
  15 + 10 min). Aggregation `20261016T110424Z-agg811labc`: $0.014.
- Arrays: `work/run-m811labc/bulk/LABC__<held>/{inc,c0_l25,c1_l25,src_l25}.npz`; list in
  `work/run-m811labc/SHA256SUMS.txt` (`fae8725c…`). The KA fold sits apart in `work/run-m811labc/ka/`.
- agg811 labels these folds "W" because they start with "L". Its board projections use W's ρ 0.19 and are not meaningful
  here.
- Tries: N = 11 (the BULK-ABC count of 10 plus 1), e(11).

**Prices** (3 LABC folds, vs inc = pct+P2g, literal pds; pdsl beside)

| arm vs base | overall [95 %] | folds > 0 | after e(11) | pdsl overall [95 %] | r | Δpds | Δmse raw |
|---|---|---:|---|---|---:|---:|---:|
| c0_l25 vs inc | +0.00230 [−0.00028, +0.00521] | 3/3 | +0.00006 [−0.00252, +0.00297] | +0.00293 [+0.00046, +0.00562] | +6.47 % | −0.00646 | −0.02653 |
| c1_l25 vs inc | +0.00248 [−0.00007, +0.00522] | 3/3 | +0.00033 [−0.00223, +0.00306] | +0.00325 [+0.00074, +0.00593] | +6.38 % | −0.00577 | −0.02616 |
| src_l25 vs inc | −0.00135 [−0.00439, +0.00204] | 1/3 | −0.00391 | −0.00108 [−0.00402, +0.00203] | +3.75 % | −0.01119 | −0.01535 |
| **c0_l25 vs src_l25 (many-line increment)** | **+0.00360 [+0.00175, +0.00546]** | **3/3** | +0.00207 [+0.00021, +0.00392] | +0.00395 [+0.00220, +0.00570] | +2.65 % | +0.00473 | −0.01086 |
| c1_l25 vs src_l25 | +0.00379 [+0.00181, +0.00572] | 3/3 | +0.00219 [+0.00022, +0.00412] | +0.00428 [+0.00244, +0.00613] | +2.56 % | +0.00542 | −0.01052 |
| c1_l25 vs c0_l25 | +0.00019 [−0.00043, +0.00081] | 2/3 | −0.00030 | +0.00032 (3/3) | −0.08 % | +0.00069 | +0.00034 |

Per fold vs inc (from the job logs):

| fold | c0_l25 r | c0_l25 Δpds | src_l25 r | src_l25 Δpds | c0 dK |
|---|---:|---:|---:|---:|---:|
| K562 | +3.94 % | −0.0029 | +0.89 % | −0.0063 | +3.3 % |
| HCT116 | +6.95 % | −0.0118 | +3.11 % | −0.0224 | +4.0 % |
| HEK293T | +8.52 % | −0.0047 | +7.24 % | −0.0049 | +4.6 % |

**The same arms on the other surfaces** (λ 0.25, vs each surface's inc):
- c0_l25 vs inc:

  | surface | c0_l25 vs inc |
  |---|---|
  | W-tune L1–L3 | +0.00853 |
  | W-conf L4–L6 | +0.00683 |
  | LABC | +0.00230 |
  | ABCp | +0.00264 |
  | 062-ABC | +0.00092 |

- many-line increment over sources-only:

  | surface | increment |
  |---|---|
  | W-tune | +0.00621 |
  | LABC | +0.00360 |
  | 062-ABC | +0.00305 |

- **Reading (descriptive, for 800's over-read work):**
  - Holding the recipe and changing only the panel, from fresh L panels to A/B/C's own Arc-drawn targets, cuts the gain of
    the full arm against inc about 3.7× (W-tune → LABC) and the many-line increment about 1.7×.
  - A/B/C's panel under the L recipe (LABC +0.0023) reads like the day-recipe A/B/C folds (ABCp +0.0026). So the panel
    effect dominates, and the operator or recipe difference between LABC and ABCp is small for this arm.
  - The increment over sources-only stays positive on all three A/B/C instruments (3/3 folds on each).
  - The sources-only blend itself turns negative on A/B/C's panel (LABC −0.0014; 062 −0.0022 in the baseline), while it
    was positive on L panels (+0.0024).
  - No decision changes: M811's §5 outcome is still WEAK (§000).

**Wear and cost.**
- Wear: 6 rows appended to `work/wear811.jsonl` (now 64 rows, sha256 `b05b9c31…`), with `n_tries` 11 and `surface_set`
  "LABC (A/B/C truths, L recipe; instrument read)".
- Cost: $0.16 for this read; the lane total is **$1.30**.
- Instances: no 811 instance is running. The live lab instance is 814's `feat814a`.

## 000. ABCp decision read (2026-10-16 09:18–09:29Z): **WEAK. No dossier. ABCp alone would NOT pass §5.**

**Inputs and checks.**
- `m811.py` still hashes `e00192e6…` (checked 09:18Z). The frozen job file `cloud/jobs-s2b-abcp.txt` (`b32f4cf1…`) was
  run unchanged. Launch `20261016T091854Z-m811s2b` (r7i.8xlarge, 2 BLAS threads; actual $0.204, cap $1.43 on 30 + 10 min)
  and aggregation `20261016T092610Z-agg811s2b` ($0.017).
- **ABCp folds**: 801's registry (`work/registry.json` sha256 `ce2ced4c…`) holds three ABCp folds: `ABCp/k562_gwps`,
  `ABCp/orion_hct116` and `ABCp/orion_hek293t`. Their dumps' sha256 equal the registry's: `16927940…`, `10f5860c…` and
  `dc7f4df6…`.
- **No ABCp/kolf2 fold exists.** The frozen job for it exited rc 1 (`FileNotFoundError: …/abcp/dump518-kolf2.npz`), so
  the run's STATUS is `FAILED 3`. The three registered folds all completed rc 0. This is not a failed known answer on a
  fold that exists.
- **G1 detection worked**: every meta file records `variant: ABCp`, `inc_arm: pct+P2g`, and the ABCp note.
- **Engine known answers** pass on every fold: PT0 cap3 0.0 and PT0 pct 0.0, plus the dump-stage KA marker. No 804
  reference exists on ABCp, so none was required.
- **Centring is mechanically as on L**: shared share drops from 0.021–0.022 to 0.0003–0.0009.
- Arrays: `work/run-m811s2b/SHA256SUMS.txt` (`e609b42d…`).

**BULK-ABCp** (3 folds, vs the October incumbent pct+P2g; N = 10 branch, e(10) 1.547; literal pds, pdsl beside)

| arm | overall [95 %] | folds > 0 | after e(10) [95 %] | bin | pdsl overall [95 %] | pdsl after e(10) | r | Δpds | Δmse raw |
|---|---|---:|---|---|---|---|---:|---:|---:|
| **c0_l25** (winner) | **+0.00264 [−0.00003, +0.00561]** | 3/3 | +0.00042 [**−0.00224**, +0.00339] | WEAK | +0.00331 [+0.00089, +0.00593] | +0.00128 [−0.00114, +0.00391] WEAK | +5.91 % | −0.00441 | −0.02424 |
| c1_l25 (C-partner) | +0.00273 [+0.00017, +0.00556] | 3/3 | +0.00058 [−0.00198, +0.00341] | WEAK | +0.00365 [+0.00120, +0.00632] | +0.00159 [−0.00087, +0.00425] WEAK | +5.82 % | −0.00398 | −0.02386 |
| c1_l25 vs c0_l25 | +0.00009 [−0.00052, +0.00070] | 2/3 | −0.00039 | NULL | +0.00035 (3/3) | −0.00008 NULL | −0.09 % | +0.00043 | +0.00036 |

**The 062 read beside it** (vs the board's ALLSW bulk, cap6+P2g(abc); §00):
- c0_l25: +0.00092 [−0.00164, +0.00364], 2/3 folds;
- c1_l25: +0.00095 [−0.00145, +0.00354].

ABCp reads about +0.0017 higher. Two things differ between the reads: the operators (062 vs the day recipe) and the
incumbent (cap6 vs pct). This run does not separate them.

**Many-line increment over sources-only on ABCp: not produced.** The frozen jobs price only c0_l25 and c1_l25, and
base811 never ran `src_l25` on ABCp. On 062 it was +0.00305 [+0.00121, +0.00492], 3/3 folds (§00).

**§5 on ABCp, item by item:**
1. **FAILS.** The 2.5 % bound after e(10) is −0.00224 (literal pds); pdsl's is −0.00114. The point is positive and all
   3 folds are > 0, but the interval after the best-of-10 shift spans 0.
2. Holds: c0 vs all_l25 on 062 is +0.00186 by point (KOLF2 removal does not hurt).
3. Holds: the W-conf bound after e(7) is +0.00234.
4. Holds: F1 does not fire (89 % of the increment kept).
5. Off-axis share is 0.26 on ABCp/k562 and ~0.01 elsewhere. The fold mean is ≈ 0.09, under 25 %, so no PARTIAL label on
   the mean.
6. **FAILS.** The winner's |dK_rel| is +5.86 % on L6/k562 (W-conf); on ABCp it is at most 4.13 %. DE-H1 is required and
   803 has no bulk_transform hook yet.

- **Outcome under §5:** the point is > 0 and the bound is ≤ 0, so the result is **WEAK: no dossier.**
- **The family is not closed:** F2 does not fire, because the point is positive.
- **Reopening condition (frozen):** only with the 813 lines (KOLF2 full rows, TeloHAEC genome-wide), as one counted try.
- **ABCp alone would not pass §5, even if rule 6 were cleared:** rule 1 fails by itself. Rule 6 blocks independently.

**On λ 0.10, which this read did not include: I do not recommend reading it now.**
- In W-tune it ranked below both λ 0.25 arms (bounds +0.00243 and +0.00232 against +0.00422 and +0.00397), so the frozen
  selection did not choose it.
- Choosing it after seeing A/B/C reads would pick a dose on A/B/C. That is the A/B/C-specific tuning AGENTS 11.1c.3 forbids
  for October.
- Its only argument is the baseline's smaller pds cost on A/B/C. But c0_l25 already reads +0.0026 on ABCp, and λ 0.10
  scales the mse gain down by about 0.4×.
- If the lab lead wants it anyway, it is a new configuration on BULK-ABCp at N = 11 (e 1.59), not part of this decision.
- DE-eligibility is not an argument for it now. λ 0.10 moved K by ≤ 1.6 % on W-tune, but rule 6 is about to be answered by
  803's hook, not by changing the dose.

**Wear:** 3 rows appended to `work/wear811.jsonl` (now 58 rows, sha256 `0a7d6791…`), with `n_tries` 10 and
`surface_set` "BULK-ABCp (N=10 branch, after 062)".

**Cloud:** this read cost $0.22. The lane total is **$1.14** of the $3 granted. No 811 instance is running; the one live
lab instance is 803's `de803-b9`.

## 00. M811 outcome (PROPOSAL-2 + ADDENDUM, GO under VERDICT-2; run 2026-10-16 08:33–08:55Z)
Frozen inputs: `PROPOSAL-2.md` sha256 `8336f260…`, `PROPOSAL-2-ADDENDUM.md` `2a3c8e95…`, `m811.py` `e00192e6…` (all in
`work/SHA256SUMS.txt`, recorded 08:33:22Z, before the first launch at 08:33:36Z).

**Decision: NO DOSSIER. The decision rule cannot yet be applied in full, because the decision set (ABCp) is not
available.** 801's `lab801-abcp-b` run ended `FAILED 3` and no ABCp dump is registered (checked 08:55Z). Even so, the
outcome on dossier is already fixed by rule 6, independently of ABCp:
- **Rule 6 fails for the winner.** The winner `c0_l25` moves ballast K by +5.86 % on L6/k562_gwps (a W-conf fold). Rule 6
  sets the bar at 5 %, and 803 has no external-rows hook to clear it on DE-H1. So whatever ABCp reads, this run produces
  no dossier.
- **What ABCp still decides:** closure or reopen. If its point is ≤ 0 (F2), the family is closed for October for c0/c1
  at λ 0.25. If its point is > 0, it is WEAK: no dossier, and the family reopens only with the 813 lines.
  - The run uses the N = 10 branch (e 1.547), because 062 was read first.
  - Job file ready: `cloud/jobs-s2b-abcp.txt` (sha256 `b32f4cf1…`). It runs `c0_l25,c1_l25` on `ABCp/*` and syncs
    801's `work/dumps/abcp`.
- **Closure scope** (VERDICT-2 G3): λ 0.25, C off and C on, on the day bank list. The read was 062-BULK-ABC; ABCp is
  pending. No other configuration is closed by this run.

### Known answers (all on r7i.8xlarge, 2 BLAS threads)
**71 KA rows, every one exact at |d| = 0.0.**
- Stage 1 (36 rows): `inc` and `707` equal 804's arrays on L1–L3 (18 rows); `c0_l10` and `c0_l25` equal the baseline's
  `all_l10` and `all_l25` (18 rows).
- Stage 2 (35 rows): `inc` and `707` equal 804's arrays on L4–L6 and the four 062 A/B/C folds (26 rows); `c0_l25`
  equals `all_l25` on L4–L6 (9 rows).
- The engine fold KAs passed on every fold.
- The `src_l*`, `all_l25` arrays copied from the baseline for the pairings were each sha-checked against
  `work/run-base811/SHA256SUMS.txt` (all match).

### Stage 1, W-tune (L1–L3, 9 folds; vs inc; N = 10, e 1.547)
Source: `work/run-m811s1/agg811-*.json`.

| arm | overall lit. pds [95 %] | best-of-10 2.5 % bound | pdsl overall | r | Δpds | Δmse raw | max \|dK_rel\| |
|---|---|---:|---|---:|---:|---:|---:|
| c0_l10 | +0.00434 [+0.00330, +0.00547], 9/9 | +0.00243 | +0.00454 | +5.89 % | +0.00035 | −0.02415 | 1.6 % |
| **c0_l25** | **+0.00853 [+0.00617, +0.01114], 9/9** | **+0.00422** | +0.00946 | +14.03 % | −0.00417 | −0.05751 | 4.9 % |
| c1_l10 | +0.00413 [+0.00314, +0.00521], 9/9 | +0.00232 | +0.00456 | +5.65 % | +0.00026 | −0.02315 | 1.1 % |
| c1_l25 | +0.00804 [+0.00581, +0.01052], 9/9 | +0.00397 | +0.00936 | +13.38 % | −0.00426 | −0.05486 | 4.2 % |

- **Eligibility (G3):** every arm's max |dK_rel| is ≤ 5 %, so all four are eligible.
- **Winner (§3.3):** `c0_l25`. Its bound +0.00422 beats c1_l25's +0.00397 by 0.00025, which is more than the 0.0002 tie
  window.
- **The winner adds nothing new on L folds:** c0_l25 is bit-identical to the baseline's all_l25 on every L fold (KA), so
  its W-conf read repeats an already-priced configuration.
- **Centring does what it should mechanically.** It cuts U's across-target shared share from 0.022–0.024 to 0.0002–0.0008
  per fold (meta `U.c0/c1.shared_share`).
- **Restated F1 (G2) does NOT fire.**
  - c1's r gain over src_lX is 1.93 % against c0's 2.16 % (λ 0.10) and 5.01 % against 5.60 % (λ 0.25).
  - So centring keeps **89 % / 89 %** of the many-line increment, where 30 % was the line.
  - On fresh folds the increment is therefore not a shared component, and the class stays transfer.
  - Caveat (addendum G2): centring also touches the source lines.
- **P1 (revised): all four checks hold.**
  - The inc-based ratio is 96 % (λ 0.10) and 95 % (λ 0.25), against "≥ 70 %".
  - |Δoverall| for c1 − c0 is 0.00020 at λ 0.10 (bar 0.0010) and 0.00044 at λ 0.25 (bar 0.0020).
  - Literal Δpds is −0.0001 at both λ, inside the predicted [−0.001, +0.003].
  - c1 − c0 at each λ: −0.00020 [−0.00042, +0.00002] (λ 0.10) and −0.00044 [−0.00092, +0.00007] (λ 0.25), literal;
    with pdsl, +0.00003 and −0.00004.
- **c1 predictions: all inside their bins.** c1_l10 overall +0.0041 (bin +0.003 … +0.005), r +5.65 %. c1_l25 overall
  +0.0080 (bin +0.006 … +0.009), r +13.4 %. Both pdsl overalls also fall inside their bins.

### Stage 2, W-conf (L4–L6, 9 folds; winner c0_l25 vs inc; N = 7, e 1.364)
- **Price:** overall **+0.00683 [+0.00421, +0.00960], 9/9 folds**. After e(7): **+0.00497 [+0.00234, +0.00774]**, so the
  2.5 % bound is > 0.
  - pdsl: +0.00791.
  - Members: r +14.00 %, Δpds −0.00883, Δmse raw −0.05741.
- **Selection over-read:** W-conf ÷ W-tune = 0.80, inside [0.7, 1.3].
- **Against the predictions:** overall, pdsl, r and Δmse fall inside their bins. **Δpds misses**: it came out −0.0088
  against a predicted −0.008 … −0.003.
- **dK_rel:** max +5.86 % (L6/k562_gwps), which triggers rule 6.
- **Off-axis change share:** 0.25–0.29 on the K562 folds and ~0.01 elsewhere.

### Stage 2, 062-BULK-ABC (3 cancer folds, kolf2 apart; N = 9 for this read)
Source: `work/run-m811s2/agg811-*.json`. In the N = 10 branch ABCp is pending.

| pair | overall lit. pds [95 %] | folds + | best-of-9 | pdsl overall | r | Δpds | Δmse raw |
|---|---|---:|---|---|---:|---:|---:|
| **c0_l25 vs inc** | **+0.00092 [−0.00164, +0.00364]** | 2/3 | −0.00107 [−0.00364, +0.00164] | +0.00109 [−0.00125, +0.00360] | +4.95 % | −0.00727 | −0.02032 |
| c1_l25 vs inc | +0.00095 [−0.00145, +0.00354] | 2/3 | −0.00097 | +0.00133 (3/3) | +4.88 % | −0.00705 | −0.02002 |
| c1_l25 vs c0_l25 (C mechanism) | +0.00003 [−0.00057, +0.00062] | 2/3 | −0.00042 | +0.00024 | −0.07 % | +0.00022 | +0.00028 |
| c0_l25 vs baseline all_l25 (KOLF2 dropped) | +0.00186 [−0.00010, +0.00394] | 2/3 | +0.00031 | +0.00046 | −0.13 % | +0.00540 | +0.00051 |
| c0_l25 vs src_l25 (many-line increment, day list) | **+0.00305 [+0.00121, +0.00492]** | 3/3 | +0.00163 [−0.00021, +0.00350] | +0.00353 [+0.00178, +0.00534] | +2.33 % | +0.00383 | −0.00957 |
| ABC/kolf2: c0_l25 vs inc / vs src_l25 | +0.00123 / +0.00864 [+0.00467, +0.01286] | 1/1 | −0.00237 / +0.00553 | +0.00185 / +0.00935 | +4.22 / +9.86 % | | |

- **Decision rule on 062** (reported only; 062 alone never decides):
  - rule 1: the 2.5 % bound after e(9) is −0.00364, so it fails;
  - rule 2 (c0 ≥ all_l25 by point, +0.00186): holds;
  - rule 3: holds;
  - rule 4 (F1): holds;
  - rule 5: PARTIAL by fold (0.25 on the K562 fold), not on the mean;
  - rule 6: fails.
- **Bins:** predicted c0_l25 −0.002 … +0.001 (literal) and −0.001 … +0.002 (pdsl). Read +0.00092 and +0.00109: inside both.
- **Dropping KOLF2 helps A/B/C.** c0_l25 vs the baseline's all_l25 reads +0.0019 literal (Δpds +0.0054). This matches
  align811's raw-delta read that the KOLF2 z-score bank hurt every A/B/C cancer truth.
- **The many-line increment also holds on A/B/C** (day list, λ 0.25): +0.00305 [+0.00121, +0.00492], 3/3 folds; best-of-9
  WEAK (literal) and SHIP (pdsl). Against the incumbent the whole blend is still only +0.0009: the sources-only blend it
  sits on costs pds on A/B/C (baseline src_l25 vs inc −0.0022).

### Mechanism read R_tgt (unweighted per-target rendered cosine, c0_l25 − src_l25)
Script `rtgt811.py` (sha256 `93b150c4…`); output `work/rtgt811.txt` (`4ce0eefa…`).

| held line | R_tgt by fold | share of targets that improve |
|---|---|---:|
| K562 | +0.0016 … +0.0034 (L1–L6), +0.0015 (A/B/C) | 0.55–0.62 |
| HCT116 | −0.0009 … +0.0011 | 0.39–0.51 |
| HEK293T | −0.0005 … +0.0002 | 0.44–0.50 |
| KOLF2 (A/B/C) | +0.0019 | 0.63 |

- **The prediction holds:** R_tgt > 0 on K562 folds and |R_tgt| < 0.002 on HCT116/HEK293T folds.
- The many-line gain on held HCT116/HEK293T (r +5 % on W) is therefore energy-weighted: strong targets gain, and the
  median target does not.

### What this run settles
- **(1) More lines and the shared component.** On fresh folds the many-line increment over sources-only is real (W-tune
  9/9 folds, both λ). It is not a shared component: F1 does not fire, and centring removes 97 % of the shared energy while
  keeping 89 % of the gain. On A/B/C the increment is also positive, +0.0031, 3/3 folds.
- **(2) Against the October incumbent: not an A/B/C gain, and not shippable at λ 0.25.**
  - On 062 it reads +0.0009 with an interval spanning 0.
  - Ballast K exceeds the 5 % DE-H1 trigger on one W-conf fold.
  - The pds cost of blending raw consensus into the whitened composition on A/B/C is what eats the mse gain.
- **(3) KOLF2 z-scores should never enter a raw consensus without per-line handling.** Dropping them gains +0.0019 on A/B/C.
- **What is still open:** λ 0.10 was never read on A/B/C in this run. In the baseline, λ 0.10 was the dose that broke even
  on A/B/C (+0.0005, KOLF2 in). A day-list, λ 0.10 read on ABCp is the obvious next try, as one counted try (N = 11 on
  BULK-ABC). So is the same model with the 813 lines (KOLF2 full rows with per-line centring, TeloHAEC genome-wide).

### Cost and tries
- **Cost:** EC2 **$0.47** for this run (m811s1 $0.204, agg811s1 $0.044, m811s2 $0.190, agg811s2 $0.035). Lane total
  **$0.92** of the $3 granted. No 811 instance is running (`labcloud.sh status` 08:55Z: none).
- **Wear:** 25 rows appended to `work/wear811.jsonl` (now 55 rows, sha256 `332844590e44…`), each with `n_tries` and
  `surface_set`: W-tune 10, W-conf 7, BULK-ABC (062) 9.
- **Files:**
  - run arrays `work/run-m811s1/SHA256SUMS.txt` (`27192af3…`) and `work/run-m811s2/SHA256SUMS.txt` (`a8a036ac…`);
  - job files `cloud/jobs-s1.txt` (`fe2d9479…`) and `cloud/jobs-s2.txt` (`fdc2b524…`).

---

Worker: Claude Opus 5.5 under `MANDATE.md` (lab lead [c61d4c]), 2026-10-16 07:06–08:25Z. Writes only in
`experiments/811-lolo-model/`. The lab's own launcher appended its usual rows to `experiments/800-checker/cloud/{runs.jsonl,
staged.tsv}` and `docs/LAB-LOG.md`. No git, no `vcc`, no production code, nothing under `data/`, no sealed panel (L7/L8)
read. `check800.py aggregate` was NOT run, because it appends to 800's wear ledger and calls git. `agg811.py` runs
stats800's own statistics instead, and the would-be wear rows are in `work/wear811.jsonl` (30 rows) for the lab lead to
append.

## 0. Answer
- **The baseline is a W-only win against the incumbent** (preregistered: `PREREGISTRATION.md`, sha256 `884405075a57…`).
  - The equal-weight consensus of all non-held public lines, blended into the composition, gains on BULK-W at every λ:
    **+0.0039 [+0.0032, +0.0048]** (λ 0.10), +0.0077 (0.25), +0.0093 (0.50); 18/18 folds.
  - On BULK-ABC it is **NULL at λ 0.10** (+0.0005 [−0.0007, +0.0017]; best-of-6 −0.0003) and HARM at 0.50 (−0.0074).
  - mse improves everywhere and literal pds pays. This is the same picture 706 drew for the sources-only mean.
- **The family's real finding: more lines carry more transferable direction** (this falsifies preregistered prediction 2).
  - At the same λ, the all-lines consensus beats the sources-only consensus on BULK-W: **+0.0018 / +0.0050 / +0.0105**
    (λ 0.10 / 0.25 / 0.50), 17/18 folds, literal pds positive too.
  - On A/B/C: positive with pdsl (+0.0011 / +0.0031 / +0.0058, 3/3 folds at 0.25 and 0.50). The literal-pds intervals span 0.
  - On ABC/kolf2: strongly positive.
- **Mechanism of the pds cost** (align811, raw delta space, nothing priced): the all-lines consensus carries 2–4× the truth's
  across-target shared energy (0.021–0.046 against 0.005–0.017). The KOLF2 z-score bank (shared share 0.245) lowers U's raw
  cosine with every A/B/C cancer truth.
- **PROPOSAL.md (M811)**: per-line centring plus reliability weights from inner leave-one-line-out among the non-held lines,
  with λ selected on L1–L3 and confirmed once on L4–L6 and BULK-ABC. It is waiting for the judge; nothing beyond the
  baseline was priced.
- **Known answers: all pass.**
  - incumbent and 707 arrays equal 804's on 22/22 folds at |d| = 0.0 (44 rows);
  - 707's pooled price equals 804's grade to every printed digit;
  - every fold's engine KA passed.
- **Cloud: $0.45** of the $50 lane budget. Five launches: two failed in align811 (an OOM, then a dtype bug), both fixed and
  logged. No 811 instance is running (§6).

## 1. Inventory (step 1; `inv811.py` → `work/inventory811.json` sha256 `efffd2d8…`, target lists only)
Coverage of each panel's 300 targets (A/B/C and the WORKING L1–L6; the sealed L7/L8 never read):

| bank | line | targets | ABC | L1 | L2 | L3 | L4 | L5 | L6 | on disk as rows? |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| k562_gwps | K562 | 9,787 | 269 | 269 | 269 | 269 | 269 | 269 | 269 | yes (composition donor) |
| viperturb_k562 | K562 (Flex, ~48 cells) | 17,489 | 284 | 271 | 274 | 276 | 273 | 267 | 272 | yes |
| perturbmatch_k562 | K562 (GWPS re-analysis, same cells) | 9,471 | 251 | 250 | 250 | 255 | 249 | 245 | 254 | yes; excluded (not independent) |
| orion_hct116 / orion_hek293t | HCT116 / HEK293T | 17,767 / 17,977 | 290 / 294 | 291 / 293 | 290 / 295 | 290 / 294 | 290 / 294 | 290 / 294 | 291 / 293 | yes (composition donors) |
| cd4_gwcd4i (3 conditions) | primary CD4 T | 11,526 | 291 | 283 | 284 | 283 | 288 | 284 | 283 | yes (panel compacts) |
| h1_pan2025 (115) / h1_vcc2025 | H1 hESC | 2,354 / 300 | 64 / 25 | 42 / 3 | 40 / 2 | 35 / 1 | 33 / 1 | 33 / 0 | 39 / 4 | yes |
| kolf2 (local) / kolf2 official list | KOLF2.1J iPSC | 273 / 10,167 | 273 / 273 | 0 / 252 | 0 / 244 | 0 / 244 | 0 / 250 | 0 / 247 | 0 / 248 | A/B/C rows only; **full rows not on disk** |
| telohaec_pilot / telohaec genome-wide (813) | TeloHAEC | 2,167 / 19,276 | 65 / 300 | 55 / 298 | 56 / 299 | 51 / 299 | 47 / 300 | 52 / 298 | 51 / 299 | pilot yes; **genome-wide not on disk** |
| dld1_lowmoi (223) | DLD-1 | 4,715 | 67 | 65 | 69 | 82 | 81 | 80 | 73 | yes |
| altair_a549 (223) | A549 | 997 | 20 | 22 | 12 | 15 | 17 | 16 | 19 | yes |
| rpe1 / nadig_jurkat / k562_essential | essential screens | 2,204 / 2,270 / 2,021 | 0 / 0 / 0 | 22 / 24 / 0 | 14 / 14 / 0 | 17 / 19 / 0 | 14 / 12 / 0 | 15 / 14 / 0 | 15 / 16 / 0 | yes |
| jiang (A549, BXPC3, HAP1, HT29, K562, MCF7) | 6 lines | 218 each | 9 | 4 | 5 | 7 | 7 | 8 | 3 | yes |
| cmap_xpr / cmap_sh_cgs | L1000 (inferred genes) | 5,157 / 4,342 | 82 / 122 | 89 / 89 | 79 / 91 | 81 / 96 | 87 / 84 | 92 / 103 | 83 / 95 | excluded (other assay; 176 closed) |
| igvf_cardiomyocyte | iPSC-CM | 1,888 | 87 | 28 | 40 | 39 | 38 | 31 | 27 | excluded (significant hits only, median 4 genes/row) |

Reading:
- On a fresh Arc-like panel, the lines on disk with near-full coverage are four: K562 (two banks), HCT116, HEK293T and CD4.
  Everything else covers fewer than ~80 targets.
- The two missing near-full lines are TeloHAEC genome-wide (813's #1) and the KOLF2 full rows (813's #3; the official list
  has 244–252 per panel).
- **Extra held-out folds:** only lines with full panel rows can be folds. Today that means the three cancer lines on every
  panel, plus kolf2 on A/B/C. TeloHAEC genome-wide would add a fourth fold line on every panel once 813 wires it.
  KOLF2 would add one on L folds with its full rows. CD4 cannot be a clean fold, because it is inside the composition.

Compact banks for the pricing jobs (`ext811.py` via gate602; peak RSS 1.20 GiB, 22 s): `work/banks/<bank>.npz`, one per
extra bank, holding the union of the A/B/C and L1–L6 targets on the official axis. 18 banks, 430 MB, git-ignored.
`work/banks/manifest.json` (sha256 `37e34335…`) records every source file's sha256 and every output's sha256. The CD4
rows of a target present in several panel compacts were asserted identical.

## 2. Baseline (step 2), preregistered before any number
Construction (frozen in PREREGISTRATION.md):
- **U_all(t)**: the mean over non-held LINES of the line's unit row (a line = the normalised mean of its banks' unit rows).
  The held line's whole group is out: on K562 folds that includes viperturb and jiang_K562; on ABC/kolf2, KOLF2.
- **U_src(t)**: the same over the fold's composition sources only.
- **Blend**: norm-preserving into the composition rows at λ (707's formula, applied to all targets), then the incumbent's
  key, P2g, render and ballast unchanged.
- Engine: `engine801.price` vs `inc`. Script: `base811.py` (sha256 `c29f9f97…`). Run on EC2 `20261016T072449Z-base811`
  (r7i.8xlarge, OPENBLAS 2 threads, 22 jobs, 72–130 s each, peak 1.55 GB).
- Arrays: `work/run-base811/bulk/<fold>/<arm>.npz`, list in `work/run-base811/SHA256SUMS.txt` (sha256 `05643d07…`).
  Per-fold records in `work/run-base811/meta/`.
- Pooled with stats800 (paired target bootstrap, 4,000 draws, seed 800; equal weight per fold) by `agg811.py` (sha256
  `d34d1e26…`) on EC2 `20261016T074655Z-agg811`. Tables in `work/table-base811.md`, json in
  `work/run-base811/agg811-811-consensus-*.json`.
- Calibration: `bulk-external` from a byte copy of `ka800.json` (`work/ka800.copy.json`, sha256 `a7085382…`). Best-of-6 for
  the family's 6 configurations.

### 2.1 Known answers (read before any number)
| row | required | measured |
|---|---|---|
| engine801 fold KAs (every fold, `E.load(check=True)`) | PT0, 540c / 522 / 546o / 551 rows ≤ 1e-12 | PASS on 22/22 folds |
| `inc` arrays vs 804's `inc.npz` (byte copies in `work/ref804/`, SHA256SUMS sha256 `31a4ad5b…`) | every array ≤ 1e-12 | **0.0 on 22/22 folds** |
| `707` rebuilt from ka804's own functions vs 804's `707.npz` | every array ≤ 1e-12 | **0.0 on 22/22 folds** |
| 707's pooled price vs 804's `grade804.txt` | same numbers | W +0.00485 [+0.00360, +0.00639], ABC +0.00044 [−0.00092, +0.00178]: **identical** |
| U_src direction vs 707's U (positive rescaling of the same rows) | min cos = 1 | recorded per fold in meta (`U_src_vs_707U_min_cos`) |

**Instrument note (the Mac smoke, `work/smoke811-abc-k562.log`).** On the Mac, 707 failed its KA on ABC/k562_gwps:
3 targets moved by |dx| up to 0.34. The cause is ballast-K flips: the key's Frobenius constant differs at the 1e-16
level between Mac and EC2 BLAS, and that flips K at boundary targets. The incumbent matched (8.5e-14). On EC2 (804's
instance type and thread count) every row is exact. **Engine KAs must be judged on the same instance type and thread
count as their reference.** The Mac smoke is not used for any number.

### 2.2 Prices (board units; literal-pds overall is the checker's standard; pdsl = board-emulated pds)
| arm vs base | surface | overall [95 %] | folds + | best-of-6 → bin | r | Δpds | Δmse raw | overall with pdsl |
|---|---|---|---|---|---:|---:|---:|---|
| **all_l10 vs inc** | W | **+0.00393 [+0.00316, +0.00476]** | 18/18 | +0.00341 SHIP | +5.86 % | −0.00072 | −0.02404 | +0.00424 SHIP |
| all_l10 vs inc | ABC | +0.00050 [−0.00069, +0.00170] | 2/3 | −0.00029 NULL | +2.33 % | −0.00322 | −0.00954 | +0.00145 [+0.00025, +0.00271] WEAK |
| all_l10 vs inc | ABCk | +0.00126 [−0.00082, +0.00343] | 1/1 | NULL | +1.93 % | −0.00034 | −0.00790 | +0.00138 WEAK |
| all_l25 vs inc | W | +0.00768 [+0.00596, +0.00958] | 18/18 | +0.00650 SHIP | +14.02 % | −0.00650 | −0.05746 | +0.00869 SHIP |
| all_l25 vs inc | ABC | −0.00092 [−0.00359, +0.00178] | 1/3 | HARM | +5.09 % | −0.01267 | −0.02088 | +0.00064 HARM (bo6) |
| all_l50 vs inc | W | +0.00932 [+0.00621, +0.01265] | 18/18 | +0.00719 SHIP | +24.10 % | −0.02198 | −0.09882 | +0.01139 SHIP |
| all_l50 vs inc | ABC | −0.00736 [−0.01218, −0.00256] | 0/3 | HARM | +6.82 % | −0.03396 | −0.02796 | −0.00550 HARM |
| src_l10 vs inc | W | +0.00211 [+0.00126, +0.00311] | 14/18 | +0.00150 SHIP | +3.89 % | −0.00188 | −0.01596 | +0.00258 SHIP |
| src_l10 vs inc | ABC | +0.00042 [−0.00097, +0.00199] | 2/3 | HARM (bo6) | +1.44 % | −0.00169 | −0.00589 | +0.00035 HARM (bo6) |
| src_l25 vs inc | W / ABC | +0.00242 / −0.00216 | 9/18 / 1/3 | WEAK / HARM | +8.67 / +2.57 % | −0.0105 / −0.0111 | | |
| src_l50 vs inc | W / ABC | −0.00172 / −0.01109 | 8/18 / 0/3 | HARM / HARM | +12.5 / +0.9 % | −0.0296 / −0.0326 | | |
| **all vs src, λ 0.10** | W | **+0.00179 [+0.00127, +0.00225]** | 17/18 | +0.00147 SHIP | +1.92 % | +0.00116 | −0.00786 | +0.00162 SHIP |
| all vs src, λ 0.10 | ABC | +0.00008 [−0.00115, +0.00125] | 2/3 | HARM (bo6) | +0.88 % | −0.00153 | −0.00360 | +0.00110 [+0.00004, +0.00215] WEAK |
| **all vs src, λ 0.25** | W | **+0.00505 [+0.00388, +0.00607]** | 17/18 | +0.00433 SHIP | +5.04 % | +0.00400 | −0.02068 | +0.00488 SHIP |
| all vs src, λ 0.25 | ABC | +0.00120 [−0.00155, +0.00384] | 2/3 | HARM (bo6) | +2.46 % | −0.00157 | −0.01009 | +0.00307 [+0.00053, +0.00565] WEAK |
| all vs src, λ 0.25 | ABCk | +0.00864 [+0.00467, +0.01286] | 1/1 | SHIP | +9.86 % | +0.00441 | −0.04043 | +0.00935 SHIP |
| **all vs src, λ 0.50** | W | **+0.01046 [+0.00834, +0.01242]** | 17/18 | +0.00911 SHIP | +10.79 % | +0.00761 | −0.04425 | +0.01049 SHIP |
| all vs src, λ 0.50 | ABC | +0.00371 [−0.00104, +0.00823] | 2/3 | WEAK | +5.88 % | −0.00137 | −0.02410 | +0.00575 [+0.00116, +0.01044] WEAK |
| oracle_vip (K562 folds; not shippable) | W (6) + ABC | r −34 … −64 % per fold, Δpds −0.14 … −0.20 | 0/7 | — | | | | |

Full table with every surface: `work/table-base811.md` (sha256 `ffeee455…`). "HARM (bo6)" means the best-of-6-shifted
point fell below −0.0005, so the decide rule labels it HARM. The raw point itself was near zero.

### 2.3 Against the preregistered bins
1. **Baseline vs incumbent: as predicted.** W is positive at every λ (+0.0039 … +0.0093). ABC is ≤ +0.0005 at λ 0.10 and
   negative at 0.50 (−0.0074, inside the predicted −0.001 … −0.004 bin's direction but larger). ABC/kolf2 was predicted
   negative at every λ: **wrong**. It is +0.0013 / +0.0012 / −0.0003.
2. **More lines: prediction FALSIFIED.** I predicted |Δ| < 0.0005 (NULL), leaning negative. It reads +0.0018 … +0.0105 on W
   (17/18 folds, literal pds positive), positive-leaning on A/B/C, and large on kolf2. The lane's premise survives its
   falsifier.
3. **Same-line ceiling: outside the bin, and not a ceiling.** VIPerturb's own K562 rows lose −34 … −64 % r against the
   composition. A ~48-cell bank is far noisier than the stack, so this probe cannot bound line-specific headroom. A
   same-line bank with depth comparable to the truth's is what it would need. None exists on disk.
- dK_rel: λ 0.10 within −0.3 … +2.3 %; λ 0.25 up to +5.9 %; λ 0.50 up to +15.4 %. Under DESIGN §4, λ ≥ 0.25 needs DE-H1
  before any dossier.
- Change energy off the held line's axis: 0.00–0.30 per fold, at or near 702's 25 % PARTIAL line on some folds. Read the
  bulk price as partial.

## 3. Diagnostic for the proposal (align811; raw delta space; prices nothing)
`align811.py` (sha256 `1e8b1264…`), EC2 `20261016T080512Z-align811c` (43 s). Output `work/align/align811.json` (sha256
`28b88bd2…`). It reports the per-target cosine of a consensus with the held line's own rows, on the held line's genes,
mean over the panel's targets, plus the across-target shared share.

| panel | held | cos src | cos all | cos non-src only | shared U_src | shared U_all | shared truth |
|---|---|---:|---:|---:|---:|---:|---:|
| ABC | K562 | 0.0278 | 0.0278 | 0.0129 | 0.010 | 0.043 | 0.011 |
| ABC | HCT116 | 0.0381 | 0.0331 | 0.0109 | 0.012 | 0.043 | 0.005 |
| ABC | HEK293T | 0.0376 | 0.0238 | −0.0015 | 0.010 | 0.046 | 0.010 |
| ABC | KOLF2 | −0.0248 | −0.0109 | 0.0142 | 0.013 | 0.020 | 0.245 |
| L1–L6 (range) | K562 | 0.016–0.018 | 0.023–0.026 | 0.016–0.020 | 0.008–0.014 | 0.021–0.024 | 0.010–0.011 |
| L1–L6 (range) | HCT116 | 0.024–0.032 | 0.025–0.032 | 0.010–0.015 | 0.011–0.015 | 0.023–0.025 | 0.005–0.006 |
| L1–L6 (range) | HEK293T | 0.024–0.029 | 0.023–0.031 | 0.009–0.017 | 0.009–0.010 | 0.021–0.023 | 0.010–0.017 |

- Ablations (the mean cosine lost when a line is dropped from U_all):
  - the cancer partners carry most of it (HCT116↔HEK293T +0.011 … +0.015: the same lab, Orion);
  - K562, CD4 and H1 help every held line (+0.001 … +0.006);
  - **KOLF2 hurts every A/B/C cancer line** (−0.009 / −0.009 / −0.015);
  - TeloHAEC pilot and A549 are slightly negative (≤ −0.0007).
- The engine read (all > src on A/B/C, pdsl) and the raw read (all ≤ src for A/B/C HCT116/HEK293T) disagree. The engine's
  aggregate cosine is energy-weighted and rendered. The raw mean cosine is per target. The proposal's centring factor is
  the experiment that separates a shared-component gain (mse mechanics) from a target-specific direction gain (transfer).

## 4. Tries and wear (AGENTS 11.1c.4)
- Family `811-consensus`: 6 priced configurations (all/src × λ 0.10/0.25/0.50) on 22 folds (18 W, 3 ABC, kolf2).
- Plus 2 known-answer arms (`inc`, `707`), already priced by 804, so not new configurations. Plus 1 diagnostic
  (`oracle_vip`, K562 folds; not a candidate).
- Wear rows for 800's ledger: `work/wear811.jsonl` (sha256 `18244b3a…`, 30 rows: 10 pairs × 3 surfaces).
- align811 read raw delta rows of all panels' fold lines. It chose nothing and fed no priced configuration. M811's factors
  were motivated by it, so count it as the reason for the proposal's grid, not as a try.

## 5. Files (sha256)
| file | sha256 |
|---|---|
| PREREGISTRATION.md | `884405075a57ad1bf7d8b57a9dfc0072ca980ddea5f0ef7495265d6bdce0f567` |
| PROPOSAL.md | `cc5b8a8fdc2ee07bf790d5e25f4f6b8c3886e0f33155efe1bb303e8e70e2b7f9` |
| inv811.py / ext811.py | `1e21a06f…229c` / `ccd606f2…e843` |
| base811.py / agg811.py / align811.py | `c29f9f97…2652` / `d34d1e26…c9f0` / `1e8b1264…9d9f` |
| work/inventory811.json | `efffd2d8af17f715cb6e6ee48f5e5c36e330fd6d73565689eeda4007b4632037` |
| work/banks/manifest.json | `37e343351000b9b65ceb4adc6a18004563219908cf03b805ed8b27215e57d83b` |
| work/run-base811/SHA256SUMS.txt | `05643d07c5566b3aff9b3e258e32233985525b36b65d7f978362b1e5052ce02b` |
| work/ref804/SHA256SUMS.txt | `31a4ad5bd9d51ad65a163a5bc4555d7327376ab3b606fbdfe4854d9d3c6dd5ed` |
| work/align/align811.json | `28b88bd22b7c81428be0acf18dd95191beb8f43ad6d04d0ccf167b7db1fc0786` |
| cloud/jobs-base.txt / cloud/sync.txt | `9fcce647…` / `58c530c1…` |

## 6. Cloud (all through `labcloud.sh`, on-demand, Project=lab)
| run id | type | outcome | actual |
|---|---|---|---|
| 20261016T072449Z-base811 | r7i.8xlarge, 22 jobs | DONE 0 | $0.319 |
| 20261016T074655Z-agg811 | c7i.4xlarge, 10 jobs | DONE 0 | $0.036 |
| 20261016T075142Z-align811 | r7i.xlarge | FAILED: OOM (an npz array re-read per target) | $0.039 |
| 20261016T080038Z-align811b | r7i.xlarge | FAILED: a bool & float mask (my bug) | $0.041 |
| 20261016T080512Z-align811c | r7i.xlarge | DONE 0 | $0.015 |

Total **$0.45**. After the last batch, `labcloud.sh status` shows no 811 instance. The one running Project=lab instance
belongs to 801 (`lab801-abcp-a`).

## 7. What is left
- **The judge reads PROPOSAL.md.** On GO: 8 configurations on L1–L3, then one confirmation on L4–L6 and BULK-ABC, about
  $0.5. DE-H1 for the winner if its dK exceeds 5 %.
- **Data (813)**: TeloHAEC genome-wide and the KOLF2 full rows are the two near-full lines missing on fresh panels. Both
  would also add fold lines.
- **For the lab lead**:
  1. Append `work/wear811.jsonl` to 800's ledger.
  2. 707 (810's lane) and M811 both blend toward donor means. Do not stack them without a joint price.
  3. The Mac-vs-EC2 K-flip caveat for engine KAs (§2.1).
- **One honest sentence:** the consensus of all public lines points the bulk better than the two cancer sources alone on
  fresh folds (W +0.005 at λ 0.25, 17/18 folds). Against the incumbent it pays in pds and is not yet an A/B/C gain, so it is
  a lead for M811, not a candidate.
