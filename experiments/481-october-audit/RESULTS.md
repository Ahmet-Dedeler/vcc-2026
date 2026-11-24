# 481 — October transfer audit: which A/B/C-tuned settings change meaning on a new panel?

Status: **complete** (2026-09-26 21:04Z – 21:55Z, resumed 21:14Z after the memory-rule stop). Analysis only: no builds, no
`def_day.sh`, no `vcc`, no `cell-eval2`, no git, no edits outside this directory. Every process ≤ 2.6 GB peak RSS
(largest: `whiten481.py`, 2.5 GB, float32 operator). Scripts here, outputs in `work/`, proposed patches in `patches/`
(see `patches/README.md`).

## Summary: ranked by expected October impact

| # | item | kind | A/B/C vs new panel (measured) | expected October impact | fix and its known answer |
|---|---|---|---|---|---|
| 1 | **F1 CD4 extraction gate (275 / 290) fires; `def_day.sh` cannot relax it** | confirmed, operational | A/B/C 283/287/284, union 291 (passes by ONE target); Arc's 2025 H1 panel 262/265/262, union 267 (FAILS); efr passes only because its pool required CD4 | **P(fires) 0.29 if D/E/F is as CD4-covered as A/B/C, ~1.0 if like the 2025 H1 panel.** The day stops in stage 1 with "setup did not complete" and no documented way through the script | `DEF_DAY_CD4_THRESHOLDS=relaxed` pass-through + explicit message (`F1_F5_def_day.diff`); the census floors still gate and PASS on the H1-2025 stand-in (union 267, complete 258) |
| 2 | **Prior gate q33 cuts at 1 on an Arc-style panel, not only on efr** (the known shift, now sized for D/E/F) | confirmed mechanism; calibrated prediction | strength ≤ 1: A/B/C 0.277 (cut 2), **essential-free 2025 H1 panel predicted 0.38–0.40 (95 % 0.32–0.45) → cut 1 in 100 % of draws**, efr 0.540 (cut 1) | if the programme + gate switches go on with q33: the strength-1 class (H1's whole gate benefit, +0.0043 H1 overall-equiv, target-bootstrap CI ≈ [−0.003, +0.013]) silently stays on the programme | 472's absolute `le1` (in flight); this makes it necessary, not cosmetic |
| 3 | **F2 re-level cliff assumes 400-cell NTC constructs** | confirmed mechanism; risk on trigger | Cp46 ∝ 1/n: 800-cell constructs ×0.500, 200-cell ×2.003, uneven 200–700 ×1.105 | if D/E/F controls are not exactly 46 × 400: **−29 % in bulk scale (≈ −0.015 overall, silent)**, **+5 % past the cliff for uneven guides (≈ −0.01, silent)**, or a misleading stop at 200 cells. Probability low (A/B/C exactly 46 × 400) | ×n/cells_per_pert in `ntclevel417.cp46` (`F2_ntclevel417_cp46.diff`): A/B/C Cp46 14.439751080335085 **bit for bit**; corrected to 1.0003 / 1.0015 / 0.9998 of the truth |
| 4 | **F5 the pre-upload gate (`def_postbuild.py`) cannot pass on a `def_day.sh` artifact** | confirmed, procedural | 38/41 on 465's final: the day's operator is gone from the restored builder, and `bulk_load` / `sameline` are not allowed diffs | the rule-10 diff that caught 239's −0.0199 will read FAIL on every correct day-of package and be skipped | postbuild inside `def_day.sh` before the restore (`F1_F5_def_day.diff`): **41/41** on 465's artifact with the builder as pinned during that run |
| 5 | F4 programme dose falls with the weak share | measured | median `prior_lambda` 0.441 (A/B/C) → 0.278 (efr); the low-λ blocks are the weak targets (469: 0.953 once they are on CLEAN4) | the 09-30 read prices the programme at A/B/C's dose, an upper bound for D/E/F | `le1` removes most low-λ blocks |
| 6 | F6 "absolute" τ 0.2 is read on carriers scaled to the panel median | measured, small | median pre-whitening norm 4.1637 (A/B/C) vs 3.8232 (efr); ¼ of efr's weak-share excess (54.0 → 47.7 % at matched scale) | small: the gate's share and DE magnitudes move ~8 % with the panel median | recorded only |
| 7 | F7 programme source degenerates when the panel overlaps it | measured, fail-open | rows left: A/B/C 275, efr 294; a panel inside the 2025 H1 panel: **0 rows → all-down, `programme_ok` true** | low probability (A/B/C shares 25 targets with that panel); bounded (all-down ≈ CLEAN4) | refuse < 150 rows (`F7_def_sameline_programme_rows.diff`): rebuilds 465's `sameline.npz` byte for byte (`b4eb2a6d…`), refuses the overlap case |
| 8 | labels hard-coded (`D,E,F` in `def_day.sh`; `is_abc` by label in the driver) | reasoned from code | — | fail-closed but ~1 h late; low probability | contexts from the manifest (in `F1_F5_def_day.diff`) |
| 9 | `dayof438.py` q grid stale (1.25 / 1.5; shipped 1.75) | confirmed from code | — | only if the optional day-of re-pricing is run; `--board-q 1.75` is refused | not patched |
| 10 | H1 union carrying-row share | measured, small | carried targets stronger than the panel median on A/B/C (×1.4 proxy), typical on efr | effective H1 share per carried row +10–30 % on efr; dose curve flat | none |

Checked with no meaning shift (measured): ballast K/KMAX/LB, the 410 and LODO key clips, normrel cap weights, the
context-relative DE constants on a new line (H1), CLEAN4 block coverage, the pinned CD4 deflation basis on an Arc-style
panel, and every delivery diagnostic in the rc meta. Details below.

**New instrument: `transfer_check481.py <def_day work dir> --bundle <bundle>`** prints every quantity above next to its
A/B/C value with a flag. Known answers: on 453's A/B/C-relabelled run every row equals the A/B/C reference (strength ≤ 1
0.277 / q33 2, NTC 46 × 400, CD4 283/287/284 u291, key floor 0.279, ballast-at-cap 0.402, carrier anchor 4.1637); on
465's efr run it flags exactly the strength / q33, programme-dose and anchor rows (`work/transfer_check_{453,465,469}.json`).
Run it on the day after `def_day.sh`, before the manual `vcc prep`.

---

## Materials (measured)
- The three synthetic bundles (`180/def-bundle`, `238/ef-bundle-r`, `238/ef-bundle`) carry byte-identical context files
  (sha256 prefixes D `a3b66213…`, E `c6f58df4…`, F `81af8449…`): A/B/C's controls relabelled. They move the PANEL only;
  context-relative knobs need a new-line stand-in (H1 Flex controls, `data/harness/h1-test.h5ad`).
- Control statistics (`ctrl_stats481.py`, ≤ 1.2 GB): D/E/F = A/B/C 18,400 cells, 46 NTC × 400, median lib 20,109 /
  19,946 / 20,034; genes > 5 CPM 9,923 / 9,610 / 10,124. H1 Flex controls: 18,400 cells, one ntc label, median lib 20,190.
- Arc's other 300-target panel: the 2025 VCC H1 panel (`data/donors/h1_2025_full/perts.txt`; 78 K562-essential
  targets, 222 not). It shares 25 targets with A/B/C. It is the only stand-in we have for how Arc builds a panel.

## Known answer (PASS): the prior gate's q33 shift is rediscovered
`strength481.py` recomputes the builder's design strength (count of render-gated |carrier| > 0.2; gate = pooled control
CPM > 10 and b_ctrl > 0) from each run's stack dump `carriers`, then the q33 threshold with `np.quantile`:

| run (dump) | panel | q33 threshold D/E/F | weak (< thr) D/E/F | strength ≤ 1 | strength = 0 | quartiles (D) |
|---|---|---|---|---|---|---|
| 453dr, 469dr | A/B/C relabelled | **2 / 2 / 2** | 83 / 78 / 76 | 26–28 % | 21 / 21 / 22 | [0, 1, 3, 8, 430] |
| 465, 452 efr, 448 efr | efr | **1 / 1 / 1** | 68 / 73 / 73 | 54–58 % | 68 / 73 / 73 | [0, 1, 1, 3, 70] |

The builder's own A/B/C counts (469dr meta) are 86 / 83 / 83: the dump holds the H1 ×2 stack carriers rather than the
level builds', so the recomputation is ~5 % low, far smaller than the efr shift. Same numbers as the 09-26 review's §2.

---

## F1 (CONFIRMED, operational, HIGH): the CD4 extraction gate fires on an Arc-style panel, and `def_day.sh` cannot relax it
`extract_cd4.py` refuses a panel with < 275 targets in any CD4 condition or < 290 in their union (047's defaults), and
`def_day.sh` calls the driver without `--cd4-thresholds relaxed` (no option or env passes it). Measured with
`extract_cd4.py --dry-run --min-targets 0 --min-union-targets 0` (0.2 GB) on every panel we hold:

| panel | Rest / Stim8hr / Stim48hr | union | gate (275 / 290) |
|---|---|---|---|
| A/B/C (Arc 2026 validation) | 283 / 287 / 284 | 291 | pass by 1 target on the union |
| **2025 VCC H1 panel (Arc's other 300-target panel)** | **262 / 265 / 262** | **267** | **FAIL** |
| efr (238, realistic) | 295 / 298 / 296 | 300 | pass (its pool REQUIRED CD4 coverage) |
| efp (238, pessimistic) | 183 / 185 / 184 | 186 | FAIL |
| 199 realistic / pessimistic | 291–296 / 168–173 | 300 / 174 | pass / FAIL |

Every `def_day.sh` rehearsal (452, 453, 457, 459, 460, 461, 465, 469) ran on a panel that passes by construction, so none could
see this. The only documented fix (DEF-RUNBOOK §0.0: re-launch with `--cd4-thresholds relaxed` after checking the mapping)
is not reachable through `def_day.sh`.

How likely on D/E/F? (`cd4cover481.py`, extract_cd4's own target resolver, split by K562-essential status.) The H1-2025
shortfall is NOT an essential-gene artefact: 21 of its 33 uncovered targets are non-essential (non-essential union rate
0.905 vs A/B/C's 0.970). Bootstrap of a fresh 300-target essential-free draw (20,000 draws):

| targets drawn like | P(gate fires) | union q05 / q50 / q95 |
|---|---|---|
| A/B/C's own targets (the most optimistic Arc-style reference) | **0.29** | 286 / 291 / 296 |
| the 2025 H1 panel's non-essential targets | **1.00** | 263 / 272 / 280 |
| both pooled | 0.96 | 276 / 283 / 289 |

Downstream of a relaxed extraction the real fail-closed gate still stands: `def_census.py` on a stand-in bundle carrying
the 2025 H1 panel (CD4 extracted with thresholds 0 into the scratchpad; `work/census_h1arc.json`) reads **GATE PASS**:
K562 278, HCT116 289, HEK293T 292 (floor 100), CD4 union 267 (floor 120), complete 258 (floor 100), 1 target with no
production donor (≤ 5). The pinned CD4 deflation basis also behaves there (captured fraction 0.132, see negatives).

Fix (`patches/F1_F5_def_day.diff`, `bash -n` clean, not rehearsed): `DEF_DAY_CD4_THRESHOLDS=relaxed` passes
`--cd4-thresholds relaxed` to the setup stages; when the gate fires, the script now dies with the coverage line and the
exact next step (mapping check with `244/mapping244.py`, then rerun with the variable). Default unchanged, so every
rehearsal's command line is unchanged. Whether the default itself should become `relaxed` (with the census as the gate)
is the coordinator's call; this audit's numbers put the strict default's firing probability between one in three (targets like A/B/C's) and near certainty (targets like the 2025 H1 panel's).

## F2 (CONFIRMED mechanism, RISK on trigger; silent in one direction): the re-level's cliff assumes 400-cell NTC constructs
`ntclevel417.py` (every re-level step of `def_day.sh`) sets claim_c = Cp46_c × Σw_c, with Cp46 = the mean over the
context's NTC constructs of cell_eval2's delete-1 jackknife ON EACH CONSTRUCT'S OWN CELLS. That proxies the truth's
per-perturbation noise only if a construct has as many cells as a reference group (cells_per_pert, 400). Nothing checks
it: the manifest carries `control_cells` and `n_ntc_ids` PER CONTEXT, and the preflight compares the file with the
manifest (self-consistent), never cells per NTC id with cells_per_pert. `cpsize481.py` (1.8 GB): as shipped **14.4398
= 180's cp_all exactly**; first 200 cells of each construct 29.82 (×2.07); constructs merged in pairs (800 cells) 7.222
(×0.500). The cliff scale goes as sqrt(claim).

**Fix, tested (`patches/F2_ntclevel417_cp46.diff`, `patches/cp46_fix481.py`, `patches/test_cp46_481.py`, 1.8 GB):**
multiply each construct's jackknife by n / cells_per_pert (from the manifest); refuse constructs < 100 cells; print the
sizes when they are not all equal. At n = 400 the factor is exactly 1.0.

| context D controls relabelled as | current rule: Cp / true | cliff-scale error | fixed: Cp / true | cliff-scale error |
|---|---|---|---|---|
| **as shipped, 46 × 400 (KA)** | 1 | 0 | **14.439751080335085 = CP46_ABC["A"], bit for bit** | 0 |
| 23 constructs × 800 | 0.500 | **−29.3 %** (silent) | 1.0003 | +0.02 % |
| 92 constructs × 200 | 2.003 | **+41.5 %** (the loop leaves the [0.7, 1.4] bound: "lower q and rerun", the wrong diagnosis) | 1.0015 | +0.08 % |
| 46 constructs of uneven size 200–700 (mean 400) | 1.105 | **+5.1 %** (silent: past the cliff) | 0.9998 | −0.01 % |

Impact if triggered: 800-cell constructs put the bulk ~30 % below its cliff (0.31 raw mse per unit scale below: ≈ 0.09
raw mse per context, ≈ −0.015 overall); uneven guides put it ~4 % past (ρ* is ~1 % below in scale; 1.78 raw mse per unit
scale above: ≈ −0.01 overall). Every check passes in both cases because the check reads the same claim.

## F3 (MEASURED): an essential-free Arc-style panel is predicted to cross the q33 boundary
Is efr (a random draw) a fair stand-in for how Arc builds a panel?
1. Raw donor proxy (`proxy481.py`: agreeing-donor count |mean raw delta| > τ on context D's gate): Spearman with the
   builder's strength only 0.67 (A/B/C) / 0.45 (efr); it under-predicts efr's weak share (0.27 vs 0.54). Too crude; it
   ranks Arc's two panels alike and above every random draw (median 7 / 7 vs 2–4 at τ 0.2).
2. **Whitened base-only proxy (`whiten481.py`, `isomap481.py`)**: 438's pre-whitening composition (3 donors, normrel cap,
   PM50), the panel's columns zeroed, whitened with the production operator at γ 0.65 with the on-axis median preserve,
   ×0.982; strength = count of gated |value| > 0.2. Spearman with the builder 0.84 (A/B/C) / 0.76 (efr); the operator's
   panel does not matter (efr with its own re-fit vs the A/B/C operator: weak share 0.723 vs 0.730). An isotonic map
   proxy → builder strength **cross-validates both ways**: fit on A/B/C predicts efr's weak share 0.533 (actual 0.540)
   and q33 1 (actual 1); fit on efr predicts A/B/C's 0.250 (actual 0.277) and q33 2 (actual 2).

| panel | strength ≤ 1 (target-bootstrap 95 %) | q33 threshold | q33 weak set |
|---|---|---|---|
| A/B/C (builder, measured) | 0.277 | 2 | {0, 1} |
| 2025 H1 panel, all 300 targets (predicted) | 0.30 (0.25–0.35) | 2 | {0, 1} |
| **2025 H1 panel, 222 non-essential = the essential-free Arc-style stand-in (predicted)** | **0.378 (0.315–0.446)** | **1** | **{0} only** |
| efr (builder, measured) | 0.540 | 1 | {0} only |

The isotonic map predicts conditional means, which blur the 0 / 1 classes. A distributional version (`knnmap481.py`:
each target's strength sampled from its 15 nearest reference targets in proxy space, 2,000 draws) validates on the
quantities that decide the gate: A/B/C references predict efr's ≤ 1 share 0.49 [0.45, 0.54] (actual 0.54) and a q33 cut
at ≤ 1 in 100 % of draws (actual 1); efr references predict A/B/C's 0.31 [0.28, 0.35] (actual 0.277) and a cut at ≤ 1 in
only 11 % of draws (actual 2). It under-predicts efr's zeros (0.13 vs 0.23), so read the 0 / 1 split as approximate.
For the **essential-free 2025 H1 stand-in: strength 0 = 0.135 [0.10, 0.18], strength 1 = 0.261 [0.21, 0.31], ≤ 1 = 0.396
[0.36, 0.43], q33 cut at ≤ 1 in 100 % of draws** (all 300 of its targets: ≤ 1 0.317, cut ≤ 1 in 17 %). So on such a
panel q33 would send only the ~13 % strength-0 targets to CLEAN4 and leave the ~26 % strength-1 class (the class behind
all of H1's reach loss) on the programme; `le1` would send ~40 % (~119 per context, vs 86 on A/B/C) to CLEAN4.

q33 flips from threshold 2 to 1 once more than 33 % of targets have strength ≤ 1. A/B/C sits just below (27.7 %); an
essential-free panel drawn the way Arc drew its H1 panel is predicted above it. So the review's efr finding is **likely
on the real D/E/F panel**. `le1` fixes it either way. efr is the weak extreme (54 %); an Arc-style panel is predicted at
~38 %, so shifts measured on efr (F4, the key's floor share) overstate D/E/F by roughly half.

## F4 (MEASURED): the programme's delivered dose falls with the weak share
With `DEF_DAY_PROGRAMME_SOURCE`, the emitter's feasible prior dose `prior_lambda` (class capacity, 0..1) has median
**0.441 on A/B/C relabelled (461dr) and 0.278 on efr (465)**: same controls, programme and builder. The gated A/B/C run
(469dr: the weakest 86/83/83 targets on CLEAN4, blocks otherwise byte-identical to 461dr's) has median 0.953, so the
low-λ blocks ARE the weak targets. The 09-30 read (458d/458g) prices the programme at A/B/C's delivered dose, an upper
bound for a panel with more weak targets; `le1` removes most low-λ blocks. No separate fix.

## F5 (MEASURED, procedural, MEDIUM): the documented pre-upload gate cannot pass on a `def_day.sh` artifact
DEF-RUNBOOK §6 item 6 requires `def_postbuild.py` PASS before each upload; `def_day.sh` never runs it. By hand on 465's
final (efr, same-line on D, programme source) against 446b's meta with the driver's own allowed-diff list: **38/41**. All
three FAILs are structural: the operator and its metadata against `WHITENER_REGISTRY['union_pruned_def']`, which after
`def_day.sh`'s EXIT restore no longer carries the day's operator (today it read a running chain's pin, e2f65ecf); and the
config diff, because `bulk_load` (the dump's file name) and `sameline` (and `de_prior_gate` when on) are not allowed. A
gate that fails on every correct day-of artifact gets skipped on the day. Fix (in `F1_F5_def_day.diff`): a postbuild step
inside `def_day.sh` BEFORE the restore, `--incumbent` = the shipped A/B/C arm (`DEF_DAY_INCUMBENT_META`, default 446b),
allowed list + `bulk_load` (+ `sameline` / `de_prior_gate` only when this run set them). Known answer on 465's artifact
with a scratch builder copy carrying 465's own operator pin (what the builder held during that run): **41/41 PASS**;
40/41 with the driver's current list.

## F6 (MEASURED, small, hidden relativity): τ 0.2 is read on carriers whose scale is anchored to the panel median
The whitening's "candidate-specific on-axis preserve" rescales the whitened rows so their MEDIAN norm equals the median
pre-whitening row norm OF THIS PANEL (`whitening.before_operator_median_norm`: 4.1637 on A/B/C, 3.8232 on efr), and the
later preserves are median-relative too. So every carrier on efr is ~8 % smaller than the same target would be on A/B/C.
Recomputing efr's strength at τ × 3.8232 / 4.1637 = 0.184: strength ≤ 1 **54.0 % → 47.7 %** (A/B/C 27.7 % either way):
about a quarter of the efr–A/B/C gap is this anchoring. `le1` at τ 0.2 still moves mildly with the panel median; the
same anchoring shrinks the DE carriers' magnitudes ~8 % on an efr-like panel (the MDE floor partly covers it). Recorded,
not proposed (it would need its own H1 known answer, and 472 is in flight).

## F7 (MEASURED, fail-open, LOW probability): the programme source degenerates when the panel overlaps it
`def_sameline.py --programme-source h1_2025_full` computes P_down over the source's perturbations MINUS the panel's
targets, with no minimum. With the script's own `programme()`: A/B/C leaves 275 of 300 rows (gene coverage 97.5 %),
efr 294; random overlaps of 100 / 200 / 250 / 280 leave 200 / 100 / 50 / 20 rows with sign agreement to the full
programme 0.946 / 0.872 / 0.820 / 0.728; a panel inside the 2025 H1 panel leaves **0 rows, coverage 0, every gene down,
`programme_ok` true**. (All-down ties CLEAN4 on H1 per 470, so the damage is bounded, but the switch would silently not be
what it claims.) Guard (`F7_def_sameline_programme_rows.diff`): refuse < 150 rows or < 50 % coverage. The patched copy
rebuilds 465's `sameline.npz` byte for byte (sha256 `b4eb2a6d…`) and refuses the overlap case.

## Checked, no meaning shift (measured unless marked)
- **Ballast `BALLAST_KMAX` 16** (`ballastK481.json`; ballast cells counted per block from `indptr`, nnz > 8,000): blocks at
  the cap 120 / 118 / 124 of 300 on A/B/C relabelled vs 108 / 101 / 101 on efr; median K 5/4/6 vs 3/3/2; clip mass max
  14,305 vs 18,971 (medians 0). The cap binds in ~40 % of A/B/C blocks already (board-validated) and less on efr.
- **410 and LODO key clips** (median-relative by design): the key median drops 2.25× / 2.4× on efr, but after the re-level
  only ratios matter. Floor share 27.9 → 33.3 % (410), 33.0 → 39.0 % (LODO), driven by targets whose shared energy is ≤ 0
  (15.6 → 24.3 %, 14.6 → 27.3 %), which any clip floors. The day-of 410 fit on the relabelled A/B/C bundle reproduces the
  pinned `d7f76a22` byte for byte (453's key json).
- **normrel cap** (`normrel481.py`): the between-donor median ratios that set two donors' relative weights are stable,
  HCT116/K562 1.045 (A/B/C), 1.105 (H1 2025), 1.032 (efr), 1.066 (efp); HEK293T/K562 0.906 / 0.943 / 0.918 / 0.982; rows
  down-weighted by more than half 6–14 % on every panel.
- **Context-relative DE constants on a new line** (`context_consts481.json`, H1 Flex controls): gated genes 9,135 / 8,668 /
  9,127 vs 9,598; control CPM at gated rank 2,200 125–131 vs 125; library q25/q50/q75 ~14.0k/20.0k/26.7k vs
  16.1k/20.2k/25.3k; MDE add median 0.088–0.093 vs 0.093; genes below min_count 0.3 6.5–8.3 % vs 8.7 %.
- **CLEAN4 block coverage** (uncovered genes default to down): 87 / 79 / 77 % of the top-2,600 gated genes on A/B/C,
  83 % on H1. The default is the weak-response lean; no fix.
- **Pinned CD4 deflation basis V5_abc** (`v5capture481.py`, lib438's verbatim `deflate_pinned`): captured fraction of the
  centred residual 0.19824 on A/B/C (= 453's meta), 0.12494 on efr (= 465's meta), **0.132 on the 2025 H1 panel
  stand-in**, inside the 0.10–0.14 range the runbook calls measured territory.
- **Delivery diagnostics** (rc meta, 453/461/469 vs 452/465): MDE dose median 1.0 everywhere; call-set size, sign-flip
  fraction, clip fraction and ballast jackknife ratio within a few % between panels.
- **Call-set size** (400 head + 2,200 CPM block per target): a count, context-relative through the CPM ranking (block =
  26–27 % of the > 5 CPM universe on A/B/C, 24 % on H1). Its n_real-relative meaning (over-calling on weak targets, a jac
  cost) grows with the weak share (F3: ~38 % vs 28 %), but per-target K is closed (385) and no count rule was re-priced
  here.
- **By construction (code read):** sharpening (sign|h|^q with the row norm restored is scale-invariant per row); q80
  protection (a per-block quantile by design); MDE law (the context's own control counts); zero-repair floor, same-line
  bank classes, identity unanimity, bank registry (absolute); re-level ρ*, tolerance and bound (absolute; the A/B/C slopes
  0.31 / 1.78 move ρ* by < 1 % in ρ if the below-cliff slope halves); census floors (loose: 100 / 120 / 100 / ≤ 5).
- **Stack weights H1 ×2 / CD4 ×1**: chosen by 438 on efp/efr, which are weaker than Arc's panels (F3). The 09-29 board read
  of 457dr (CD4 ×1.5) vs 453dr (CD4 ×1), both on the October path over A/B/C, is at least as relevant to D/E/F.
- **Zero-gene audit** (`audit428.py`): prints, never exits non-zero; `def_day.sh` logs its first line only. A report, not a
  gate (the repair itself runs in every build).
- **Labels** (code): `def_day.sh` hard-codes `D,E,F`; the driver picks the A/B/C env by label. Fail-closed, but ~1 h late
  (after the dump and lvl0); fixed in `F1_F5_def_day.diff` by reading the manifest.
- **`dayof438.py`** (code): `--board-q` accepts 1.25 / 1.5 only and `recommend438` prices the H1 × CD4 surface at q 1.25;
  the shipped q is 1.75. `def_day.sh` does not call it, so it matters only for the optional day-of re-pricing.
- **H1 union carrying-row dose** (`h1share481.py`): g is sized against the base's median over ALL rows; on A/B/C the 77
  H1-carried targets are stronger than the panel median (proxy ×1.43, strength 6 vs 3), on efr the 35 are typical
  (proxy ×1.0, strength 2 vs 1), so efr's effective H1 share per carried row is ~10–30 % larger. 438 found H1 ×2 holding
  on efr/efp; no fix.

## Files
- Scripts: `ctrl_stats481.py`, `strength481.py` (KA), `census481.py`, `cd4cover481.py` (F1), `cpsize481.py` (F2),
  `proxy481.py`, `whiten481.py`, `isomap481.py`, `knnmap481.py` (F3), `h1share481.py`, `normrel481.py`, `v5capture481.py`,
  `transfer_check481.py` (the instrument).
- Outputs (`work/`): `ctrl_*.npz`, `strength_*.json`, `census481.json`, `cd4cover481.json`, `cpsize481_D.json`,
  `test_cp46_481.json`, `proxy481.json`, `whiten481.json`, `isomap481.json`, `isomap481_classes.json`, `knnmap481.json`, `h1share481.json`, `normrel481.json`,
  `ballastK481.json`, `context_consts481.json`, `clean4_block_cov481.json`, `v5capture481.json`, `census_h1arc.json`,
  `postbuild465_vs446b.json`, `transfer_check_{453,465,469}.json`, `panels/*.csv`.
- Patches: `patches/README.md`.
- Scratchpad only (temporary, outside the repo): the 2025 H1 panel stand-in bundle and its CD4 extraction (222 MB), a
  builder copy carrying 465's operator pin for the F5 known answer, the patched-`def_sameline` test outputs. One stray
  `/tmp/abc_*.txt` from a mistyped redirect was deleted at once.

## Not done / open
- None of the patches is applied or rehearsed. F1_F5 needs one rehearsal on efp (CD4 union 186) with the variable set and
  one on ef-bundle-r that must retrace 452/465.
- F3's prediction rests on one stand-in panel (Arc's 2025 H1 panel, 222 non-essential targets) and a base-only proxy; it
  is a calibrated estimate, not a build. The day itself should read `transfer_check481.py`'s strength row.
- The coordinator's 470b / 471 / 472 chains were not touched; 480's nmae question and 482's literature hunt not duplicated.
