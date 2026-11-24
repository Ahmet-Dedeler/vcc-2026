# 553 — R7b's instrument: the real-emitter DE harness on a sealed Calu-3 surface

Worker 553 (Claude Opus 5.5, headless), mandate `MANDATE.md` (2026-09-30 ~22:00Z). Writes only in this directory; no git,
no `vcc`, no `cell-eval2` CLI, nothing built for submission. **The R7b truth file has not been opened, loaded or parsed by
anything in this directory** (only its sha256 may be checked; `lib553.open_truth` refuses it — and 552's floored truth — unless `access553.py
r7b|r7bd5000 --owner-go <claim>` armed the access in-process).

**Status: v6 (2026-10-02 01:10Z), final.** Items 1, 2 and 4 are DONE and pass every known answer: the surface-parametric
harness reproduces 537's / 528's / 517's recorded H1 runs bit for bit (5 arms, pre-zs and zs, plus R7b's own code path
on H1; KA-FINAL re-ran the final code: bit-equal to 528's `cpm` and 537's CPH), 323's donor-effect construction and
control reference; the access reproduces 323's H1 reference and anchors exactly. Item 3 is DONE: R7b's prediction side
from its controls and the three public donors (operator refit under the mandate's own gate, sha `d846b222…`; no
forbidden source). **Item 5: R7b as built is unfit for the production emitter** (107 / 107 placements fail; under any
class rule the CLEAN4 prior / MDE / hook dose is 0); **a 5,000-count human-library floor restores it** — the floor sweep
shows no floor restores everything (higher floors give the full MDE / hook dose but lose the call set and squeeze the
K-size switch), 5,000 keeps the call set and the K-size switch as on H1. **552 adopted it (Amendment 3) and built the
floored surface; the dry run verifies it against the prediction side.** The 8 arms x 3 seeds are emitted on it.
**The band, measured:** the mandate's band (emission seeds + split-half re-draws of the truth) understates how far an
R7b-shaped reading moves — five independent rehearsal draws spread 1.4–2.0x wider for the call-set switches — because
the CONTROL sample (prediction side and truth) moves the reading as much as everything else together. A control re-draw
(the surface's controls split into halves, each half's prediction side rebuilt and emitted again) measures that part on
the surface itself; with it the band matches the rehearsal's between-draw spread (validated on h1r7). R7b's 10 control
halves (5 splits) are emitted and the dry run confirms them; the owner's access runs everything in one command (tested
end to end on the rehearsal and reproducible to 1e-15; ~25 min, < 4 GB on R7b). Shape-matched
expectations (five draws): CP +0.0014, CPH +0.0015, K 2,000 −0.0016 overall-equivalent (H1: +0.0051, +0.0065, −0.0036),
hook 0. The R7b truth has not been opened, loaded or parsed by anything in this directory (sha256 checks only).

## For the owner: what the claim needs

1. **The surface is `r7bd5000`** — 552's floored build (`552-r7-reserve/work/r7b_floor5000`, 08:38Z): 106 targets
   (BST2 excluded), 16,381 target cells (59–556 per target, median 127, 22 under 100), 2,378 controls. `access553.py
   r7bd5000 --dry-run ... --halves 0,1,2,3,4` verifies its 5 outputs by sha256 (the truth hashed only), its target
   list and its controls (column sums and count) against the prediction side, and the 10 control halves' readiness
   (passed 10-02 00:54Z).
2. **The band to preregister is `band_total`** (`work/band553_r7bd5000-xBST2.json`, computed inside the access from
   R7b's own control halves), not the per-arm `band_sd` alone: on the rehearsal `band_sd` was 0.0004–0.0013 and
   `band_total` 0.0006–0.0018, against a between-draw spread of 0.0004–0.0031 (the largest for K-size effects, whose
   between-draw spread also carries the draws' different cells-to-targets assignment, which R7b does not have).
3. **Board-retrodiction rows** (552's requirement): `--retro work/access553_rehearsal5` gives each raw member's effect
   over the SHAPE-MATCHED rehearsal mean (five h1r7 draws). Against H1-full, R7b's shape alone produces a board-like
   pattern for the call-propensity block — rehearsal mean jac 0.45x H1's, fid 1.4x, reach −0.05x (board: 0.48x, 1.95x,
   −0.38x) — so R7b-over-H1 would mix the shape into the cell-line comparison; R7b-over-rehearsal does not.
4. **Expect the call-set switches at about a third of their H1 size** (CP +0.0014, CPH +0.0015 overall-equivalent on the
   rehearsal mean), the hook alone and the hook on CP at ~0, and jac's LEVEL below its baseline anchor for every arm (the
   shape, not the switch).

## Files

| file | what |
|---|---|
| `lib553.py` | the surface library: `surface_spec` (h1, h1cb, h1r7[sN], h1t7[sN], r7b, r7bd<floor>, <base>c<k><A|B> control halves, fx), `controlref` (evaluate323.reference's control half), `truth_reference` (its truth half, any cell counts), `score_stats` (evaluate323.score with the reference passed in), `anchors` (anchors323's baseline + replicate, group by group), `scaled` (cell-eval2's own `score_one` at the replicate anchor), the sealing guard |
| `run553.py` | the harness: run537's MODE grammar on a SURFACE; drives run496.main unchanged (production emitter 0f56209d asserted) |
| `effects553.py` | 323's stage -> operator -> compose for any panel (cap15 arm) |
| `prep553.py` | the prediction side's control-only steps: control reference, call-propensity key (528's callblock528.py), CD4 banks (047's extract_cd4.py), audit; `subset r7bd<floor>` (the floored control file), `h1r7 build` (the R7b-shaped rehearsal surface from h1-test), `h1t7 build` / `h1t7 <h1r7 draw>` (its ablation: a draw's target cells + all h1-test controls), `chalf <base> <k>` (split k of a surface's controls into halves A / B, each a surface of its own) |
| `access553.py` | the access: truth reference, anchors, arm x seed scoring, split-half reference re-draws, paired effects vs `cpm`; `--ctrl-splits` (cells + controls re-draw, truth side), `--pairs` (switch-on-recipe effects), `--retro` (board-retrodiction ratios), `--halves` (the control-half accesses + band553), 552's floored-build checks, the halves' readiness in the dry run |
| `ka553.py` | KA-H1 comparisons (stats bit-equality + |d| per member) |
| `test553.py` | fixture tests of the group-by-group functions against the full-AnnData ones |
| `emitdiag553.py` | item 5: emission-side diagnostics of a surface's arms (K, call sets, identity checks vs `cpm`); reads only predicted blocks' own statistics |
| `kascore553.py`, `katable553.py`, `fx553.py` | KA-S0 (scorer on 323's stored reference), the KA table, the synthetic `fx` surface for access tests |
| `diag553.py` | emission-side diagnostic: run496.main's per-target setup with the emitter's `place_k` instrumented (class layout, residuals, prior / MDE dose λ) under any class rule (`--classes q25,q50,q75`, `spread<r>`, `h1ratio`, `r7bratio`, …); failures caught; reads no truth |
| `libq553.py` | control-library quantiles per context (A / B / C / H1 / R7b, and R7b under a floor) |
| `cmp553.py` | side-by-side of access reports (h1-test vs the rehearsal: effects, seed / split-half / band sd, bootstrap, anchors, cost) |
| `draws553.py` | the rehearsal's between-draw spread: per arm and pair, the sd of the effect over N independent R7b-shaped draws against the access's mean within-surface band (cells; cells + controls) |
| `band553.py` | the control re-draw band of a surface from its control halves (`prep553.py chalf`): sigma_c = rms(effect difference between halves) / 2, band_total = sqrt(band_sd^2 + sigma_c^2) |
| `prune553.py` | deletes regenerable `.h5ad` intermediates under `work/` only (records each in `work/PRUNED.md`) |
| `BLOCKED.md` | RESOLVED: the R7b operator fit ran under the mandate's own 27 GiB gate (no override written) |
| `drive553.sh` | serial gated driver (mandate gate + blackout 10-04 05:00Z – 10-05 18:00Z + 8 GB RSS watchdog); queue `work/queue553.txt`, log `work/drive553.log`; 402 jobs ok, 1 expected failure (`r7b_cpm_s0`: R7b as built, the emitter's AssertionError); stopped 10-02 01:04Z with an empty queue (`work/STOP553`; delete it to restart). Pruned regenerable files: `work/PRUNED.md` |

## Known answers

| KA | result |
|---|---|
| T1–T4 fixtures (`test553.py`): `read_rows` == AnnData slicing; group-by-group `anchors` == anchors323's `stats` / `build_baseline_prediction` + `fast_score` / `replicate` on 323's 3-target fixture and on a 5-target interleaved fixture with odd group sizes and a categorical label; `controlref` == evaluate323's control half; the sealing guard refuses R7b's truth | **PASS**, max \|d\| 0.0 on every member of baseline, replicate mean and each split |
| KA-CB: 528's `callblock528.py --ka` (the October key formula from h1-test's controls == keys528's `lcap`) | **PASS**: max \|diff\| 8.9e-16, blocks identical on 85/85 targets |
| KA-CD4: 047's `extract_cd4.py` re-run for 273's H1 panel == 273's `cd4_compact_h1` | **PASS**: 36/36 array/label files equal |
| KA-S0: `lib553.score_stats` on 323's stored reference == recorded metrics json (528 `cpm` s0 zs, 528 `lcap-kept` s1 zs, 537 `lcap-kept+c4h` s2 pre-zs); `scaled` (cell-eval2's `score_one`) == report323's formula | **PASS**: \|d\| 0, 1.1e-16 (nmae), 0; scaled vs formula ≤ 6e-17; 0.3–1.4 s per scoring |
| KA-C: `controlref` on h1-test's 18,400 controls == 323's reference.npz control arrays (keep, cs, tie, ref_lin, emit_cs, emit_tie, emit_ref_lin, libraries, b_ctrl, the control row's sums / jackknife / count) | **PASS**, 13/13 bit-equal (17 s, peak RSS 7.87 GB; a first try with a view-retaining reader was killed by the 8 GB watchdog at 8.46 GB) |
| KA-E: `effects553.py stage` + `compose` on H1's 85 targets (current production builder, 323's model_env, 323's operator a46fd323…, 273's H1 CD4 banks) == 323's cap15 | **PASS**: stage — baseline / candidate pre-whitening shas, the staged cap15 matrix, the 18,157 operator columns and the CD4 arrays equal 323's (6.6 s, 7.32 GB); compose — final sha `565d57ab…` equal, effects.npz cap15 bit-equal (2.5 s, 4.35 GB) |
| KA-H1 (`ka553.py`, table below) | **PASS** on all 12 compared stats files |
| KA-FINAL (10-01 04:30Z): the FINAL `run553.py` (after the 10-01 `--exclude-targets` / `--class-rule` / `--opt=value` / NaN-safe-json edits; cached bulk key; `--tag=-kafinal`) re-run on h1 and compared with the RECORDED runs: `cpm` s0 vs 528's `h1_cpm`, `lcap-kept+c4h` s0 vs 537's `h1_lcap-kept+c4h` | **PASS**: p, q, lfc, sums, jk, targets bit-equal, pre-zs and zs, both arms (296–299 s, peak 5.2–5.5 GB: the cached key spares the operator load) |
| T5 / T2b fixtures: `truth_reference` == evaluate323.reference's arithmetic; the access's lean anchor path (control block and control moments reused from the truth reference, control mean streamed) == the standalone path == 323's functions | **PASS** |
| `access553.py fx` (a synthetic 5-target surface, 3 fake arms x 3 seeds): the whole access end to end, report format | runs (1.6 s, 0.33 GB) |
| access on h1-test, KA-R: `lib553.truth_reference` (control rows + target rows, any counts) == 323's reference.npz (sums, counts, jk, p, q, lfc, keep) and reference_de.parquet | **PASS**, all bit-equal |
| access on h1-test, KA-S: `score_stats` on 18 recorded stats files (6 arms x 3 seeds; 435 / 528 / 537 / this directory) == their metrics json | **PASS**, max \|d\| 1.1e-16 |
| T7 (`test553.py`): a floored prediction side (`r7bd<floor>`) without 552's floored truth is refused by the access's orchestrator before any phase, and by `truth_rows` | **PASS** (stand-in spec; no R7b file touched) |
| KA-SWEEP: `diag553.py r7bd5000 / r7bd8000 cpm 0 --exclude-targets BST2` == the full `run553.py` runs' per-target emitter record (call set prescribed / delivered / missed / unsuppressible, ballast K, prior and MDE λ) | **PASS**, all 106 targets equal at both floors |
| KA-PAIRS: the report phase re-run with `--pairs` on the h1, h1t7 and h1r7 accesses leaves every vs-`cpm` row (raw, scaled, sds, bootstrap) unchanged | **PASS**, `arms` sections identical |
| KA-ACCESS-FINAL (10-01 22:37Z): the FINAL `access553.py` (with the csplit phase, `--pairs`, `--retro`, `--halves`, the floored-build checks) re-run on h1-test | **PASS**: KA-R all arrays bit-equal, KA-S max \|d\| 2.2e-16 on 18 files, KA-A max \|d\| 0.0; its arm rows equal the first h1-test access's within 7.4e-16 (425 s, peak 7.51 GB) |
| KA-E2E: the one-command access with `--halves 0,1,2` on h1r7 (main access + cells-and-controls re-draw + retrodiction + 6 control-half accesses + band553), from inputs rebuilt after pruning (all 7 rebuilt files sha-identical to the first builds) | **PASS**: every effect, band, anchor and the band553 output equal the earlier step-by-step runs within 1.1e-15 (cell-eval2's nmae jitter); 864 s, peak 3.75 GB |
| access on h1-test, KA-A: the group-by-group anchors (cell-eval2's generic profile + dispersed seed-0 baseline + fast_score; the 5-split replicate, one half per process) == 323's h1-anchors.json, baseline, replicate mean and all 5 splits x 6 members | **PASS**, max \|d\| 2.2e-16 |

### KA-H1 detail (`work/ka553.json`)

| surface | mode | seed | form | recorded run | arrays bit-equal (p, q, lfc, sums, jk, targets) | max \|d\| member (which) | pass |
|---|---|---|---|---|---|---|---|
| h1 | `cpm` | 0 | pre-zs | `528-callprop-block/work/h1_cpm/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.npz` | True | 1.1e-16 (nmae) | **PASS** |
| h1 | `cpm` | 0 | zs | `528-callprop-block/work/h1_cpm/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.npz` | True | 0.0e+00 (pds) | **PASS** |
| h1 | `lcap-kept` | 0 | pre-zs | `528-callprop-block/work/h1_lcap-kept/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.npz` | True | 0.0e+00 (pds) | **PASS** |
| h1 | `lcap-kept` | 0 | zs | `528-callprop-block/work/h1_lcap-kept/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.npz` | True | 1.1e-16 (nmae) | **PASS** |
| h1 | `cpm+c4h` | 0 | pre-zs | `517-flex-programme/work/h1_c4-c4t0.05/stats-base_s0_bps_hookem_mde0.1cap0.1min0.3coarse.npz` | True | 2.2e-16 (nmae) | **PASS** |
| h1 | `cpm+c4h` | 0 | zs | `537-arm-pipeline/work/h1_cpm+c4h/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.npz` | True | 2.2e-16 (nmae) | **PASS** |
| h1 | `cpm-n1600` | 0 | pre-zs | `537-arm-pipeline/work/h1_cpm-n1600/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.npz` | True | 1.1e-16 (nmae) | **PASS** |
| h1 | `cpm-n1600` | 0 | zs | `537-arm-pipeline/work/h1_cpm-n1600/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.npz` | True | 1.1e-16 (nmae) | **PASS** |
| h1 | `lcap-kept+c4h` | 0 | pre-zs | `537-arm-pipeline/work/h1_lcap-kept+c4h/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.npz` | True | 0.0e+00 (pds) | **PASS** |
| h1 | `lcap-kept+c4h` | 0 | zs | `537-arm-pipeline/work/h1_lcap-kept+c4h/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.npz` | True | 1.1e-16 (nmae) | **PASS** |
| h1cb | `lcap-kept` | 0 | pre-zs | `528-callprop-block/work/h1_lcap-kept/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.npz` | True | n/a (sealed stub scores nothing) | **PASS** |
| h1cb | `lcap-kept` | 0 | zs | `528-callprop-block/work/h1_lcap-kept/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.npz` | True | n/a (sealed stub scores nothing) | **PASS** |

`h1cb` = h1's inputs through R7b's code path (the sealed scoring stub and the October call-propensity key from h1-test's own
controls, 528's callblock formula): bit-equal to 528's `lcap-kept`.

Peak RSS of the H1 emission runs 7.02–7.47 GB (`/usr/bin/time -l`; 537's were 6.99–7.61), 299–324 s each.

## The access, tested end to end on h1-test (`work/access553_h1/access553.md`)

h1-test as the stand-in truth (85 targets x 400 cells, 18,400 controls; n_real median 89, q10 2.4, q90 2,691), 6 arms x 3
seeds of zs-form stats (seed 0 from this directory's KA runs, seeds 1–2 the recorded runs), 5 split-half reference
re-draws (10 half-truths), anchors computed inside the access (baseline pds 0.510 / mse 0.971 / nmae 0.960 / fid 0.035 /
reach 0.120 / jac 0.036; replicate 0.954 / 0.0087 / 0.379 / 0.852 / 0.976 / 0.461 — 323's values exactly).
Overall-equivalent = mean of the six members scaled at these H1-local anchors (cell-eval2's own `score_one`; it equals
report323's formula to 6e-17 here). These are H1-LOCAL units: 192 / 537's board-fitted units weight fid ~3x more (H1's own
fid baseline is 0.035 against the board's 0.51), so K-size arms read larger there.

| arm (vs `cpm`, zs) | overall-equiv (3 seeds) | target 95 % CI | seed sd | split-half sd (half-truth) | band sd (suggested) | raw members (fid / jac / reach / nmae) | record |
|---|---:|---|---:|---:|---:|---|---|
| `lcap-kept` (CP) | **+0.00505** | [+0.0036, +0.0067] | 0.00027 | 0.00114 | 0.00082 | −0.0013 / +0.0091 / +0.0096 / +0.0004 | 528: +0.0051 [+0.0033, +0.0070]; jac +0.0091, reach +0.0096, fid −0.0013 |
| `cpm+c4h` (hook alone) | +0.00019 | [−0.0017, +0.0020] | 0.00084 | 0.00043 | 0.00057 | +0.0000 / −0.0001 / −0.0011 / −0.0015 | 537 P1: +0.0002 [−0.0017, +0.0019] |
| `lcap-kept+c4h` (CPH) | **+0.00647** | [+0.0045, +0.0084] | 0.00027 | 0.00106 | 0.00076 | −0.0006 / +0.0093 / +0.0131 / −0.0014 | 537 P1: +0.0067 [+0.0047, +0.0088] |
| `cpm-n1600` (K 2,000) | −0.00364 | [−0.0062, −0.0013] | 0.00030 | 0.00129 | 0.00093 | −0.0100 / −0.0029 / −0.0020 / +0.0003 | 537 P2 (board units): −0.0073 at f 1 |
| `lcap-kept-n1600+c4h` (CPHK20) | +0.00020 | [−0.0028, +0.0033] | 0.00036 | 0.00189 | 0.00135 | −0.0130 / +0.0062 / +0.0006 / −0.0011 | 537 P5 (board units): −0.0040 at f 1 |

What the H1 test says about the reading rule's band (the claim's to fix, not mine): on 85 targets x 400 cells the
reference re-draw dominates the emission re-roll for the block arms (split-half sd of the seed-averaged paired effect
0.0011–0.0019 per half-truth, i.e. ~0.0008–0.0013 at full truth size, vs a seed sd of 0.0003); for the hook alone the
seed term dominates (reach). The suggested band sd = sqrt(seed sd^2 / n_seeds + (split-half sd / sqrt 2)^2) is
0.0006–0.0014 overall-equivalent on H1. R7b has fewer truth cells per target and fewer controls; on the R7b-shaped
rehearsal below the band does not widen (0.0005–0.0015), the call-set arms' effects shrink, and the band itself turns
out too narrow at that shape: it cannot see the control sample (five draws and the control re-draw, below).

Phase cost on h1-test: ctrl 9 s / 4.3 GB, ref 30 s / 6.5 GB, score 8 s / 1.7 GB, split 152 s / 7.5 GB, anchors' baseline
36 s / 6.6 GB, ten replicate halves 14.5 s / 5.0–5.7 GB each, combine 4 s, report 2 s: 390 s in all, measured peak
7.53 GB (the 8 GB watchdog never fired). An earlier monolithic version was killed at 8.01 GB and a single-process
replicate at 9.31 GB (macOS RSS, see findings); both were restructured before any number was read.

## The R7b-shaped dress rehearsal (`h1r7`): what an R7b access can resolve

h1-test is not a reserve, so it can be cut to R7b's shape and accessed freely: `prep553.py h1r7 build --controls 2095
--cells-scale 0.59` keeps 2,095 of its controls (R7b ≥ 8,000's count) and gives each of the 85 targets R7b's
cells-per-target distribution x 0.59 (build_r7b.json's COUNTS, a metadata read; 0.59 = the controls' kept fraction at
8,000, an assumption for target cells): 61–400 truth cells per target, median 110, 29 targets under 100; seed 55307.
Its prediction side is rebuilt from those controls only (control reference, call-propensity key; classes 16,128 /
20,275 / 25,325, prior and MDE dose ~1 as on H1), 537's 8 arms x 3 seeds emitted (285 s, ~2 GB each), and the access
run end to end exactly as on R7b (`work/access553_h1r7/access553.md`; side by side: `cmp553.py`,
`work/cmp553_h1_h1r7.md`). Effects are in each surface's own overall-equivalent units (anchors computed inside the access).

| arm | surface | overall-equiv effect | seed sd | split-half sd (full-est) | band sd | effect / band sd | target-bootstrap 95 % | raw nmae / fid / reach / jac |
|---|---|---:|---:|---:|---:|---:|---|---|
| `lcap-kept` (CP) | h1 | +0.00505 | 0.00027 | 0.00081 | 0.00082 | +6.1 | [+0.0036, +0.0067] | +0.0004 / −0.0013 / +0.0096 / +0.0091 |
| | **h1r7** | **+0.00270** | 0.00062 | 0.00063 | 0.00073 | **+3.7** | [−0.0023, +0.0075] | −0.0005 / −0.0029 / +0.0070 / +0.0040 |
| `cpm+c4h` (hook) | h1 | +0.00019 | 0.00084 | 0.00030 | 0.00057 | +0.3 | [−0.0017, +0.0020] | −0.0015 / +0.0000 / −0.0011 / −0.0001 |
| | **h1r7** | +0.00045 | 0.00041 | 0.00038 | 0.00045 | +1.0 | [−0.0003, +0.0016] | +0.0001 / +0.0000 / +0.0022 / +0.0000 |
| `lcap-kept+c4h` (CPH) | h1 | +0.00647 | 0.00027 | 0.00075 | 0.00076 | +8.5 | [+0.0045, +0.0084] | −0.0014 / −0.0006 / +0.0131 / +0.0093 |
| | **h1r7** | **+0.00219** | 0.00056 | 0.00069 | 0.00076 | **+2.9** | [−0.0024, +0.0064] | −0.0005 / −0.0031 / +0.0048 / +0.0041 |
| `cpm-n1600` (K 2,000) | h1 | −0.00364 | 0.00030 | 0.00091 | 0.00093 | −3.9 | [−0.0062, −0.0013] | +0.0003 / −0.0100 / −0.0020 / −0.0029 |
| | **h1r7** | −0.00265 | 0.00047 | 0.00125 | 0.00128 | −2.1 | [−0.0067, +0.0007] | −0.0011 / −0.0007 / −0.0150 / +0.0008 |
| `lcap-kept-n1600+c4h` (CPHK20) | h1 | +0.00020 | 0.00036 | 0.00133 | 0.00135 | +0.1 | [−0.0028, +0.0033] | −0.0011 / −0.0130 / +0.0006 / +0.0062 |
| | **h1r7** | −0.00249 | 0.00051 | 0.00114 | 0.00118 | −2.1 | [−0.0098, +0.0036] | −0.0022 / −0.0043 / −0.0203 / +0.0049 |
| `cpm-n1200` (K 1,600) | h1r7 | −0.00296 | 0.00036 | 0.00151 | 0.00152 | −1.9 | [−0.0082, +0.0016] | −0.0019 / −0.0028 / −0.0153 / +0.0005 |
| `lcap-kept+c4h-s` | h1r7 | +0.00271 | 0.00059 | 0.00069 | 0.00077 | +3.5 | [−0.0026, +0.0086] | +0.0009 / −0.0030 / +0.0095 / +0.0041 |

| surface | targets | controls | truth cells / target (min / median / max) | anchors baseline (pds / mse / nmae / fid / reach / jac) | anchors replicate | `cpm` overall (scaled) | access s | peak RSS GB |
|---|---:|---:|---|---|---|---:|---:|---:|
| h1 | 85 | 18,400 | 400 / 400 / 400 | 0.510 / 0.971 / 0.960 / 0.035 / 0.120 / 0.036 | 0.954 / 0.009 / 0.379 / 0.852 / 0.976 / 0.461 | 0.2990 | 384 | 7.53 |
| h1r7 | 85 | 2,095 | 61 / 110 / 400 | 0.483 / 1.002 / 0.963 / 0.032 / 0.180 / 0.071 | 0.873 / 0.034 / 0.491 / 0.810 / 0.963 / 0.477 | 0.3224 | 157 | 3.69 |

What it says about a floored R7b access (the same targets and arms, only the shape changed; the shape is the 8,000
floor's, chosen before the floor sweep — the recommended 5,000 keeps ~13 % more controls and, est., cells, so the
rehearsal is mildly conservative):

- **The noise band does not widen; the effects shrink.** Band sd 0.00045–0.00152 (h1: 0.00057–0.00135). With 2,095
  controls and ~110 cells per target the seed term grows (CP 0.00062 vs 0.00027) but the split-half term does not
  (0.00063 vs 0.00081): the paired differences are smaller and so is their re-draw noise in raw units (CP's jac 0.00011
  vs 0.00015, reach 0.0044 vs 0.0060), and the scaling at h1r7's own anchors is close to h1's (jac 1/0.406 vs 1/0.425).
  But the call-set arms' effects shrink — on this first draw CP +0.0027 (h1 +0.0051), CPH +0.0022 (h1 +0.0065) —
  because jac's gain halves (+0.0040 raw vs +0.0091) and fid turns more negative: a smaller truth call set leaves less
  overlap to win. (This first draw read CP / CPH at ~3–4 band sd; the five draws below put the shape's mean at
  +0.0014 / +0.0015 and show the reading moving ~0.0019 between draws, so one draw's "3–4 band sd" was luck.)
- **jac sits below its baseline anchor for every arm** (`cpm` raw 0.032 vs baseline 0.071; scaled −0.096; h1: 0.075 vs
  0.036). A 2,600-call predicted set against a truth call set this small loses jac to the generic baseline. The member
  is still informative as a paired difference (its band sd is 0.0002), but its LEVEL on R7b will look poor; that is the
  shape, not the switch.
- **The target bootstrap is 3–6x wider than the band** ([−0.0023, +0.0075] for CP): with ~110 cells per target the
  per-target effects are noisy and heterogeneous. The band (emission seeds + truth re-draws) answers "is this panel's
  effect distinguishable from noise"; the bootstrap answers "would another panel of targets see it". A HARM rule on the
  band can detect a CP-sized harm; one on the bootstrap CI cannot.
- **Cost.** The whole access (8 arms x 3 seeds, 10 half-truths, anchors inside) took 157 s with a 3.69 GB peak (the
  split phase); R7b ≥ 5,000 has 106 targets and ~16k cells (est.), so expect ~3 min and ~4 GB.

**Ablation: which half of the shape costs what (`h1t7`).** `prep553.py h1t7 build` takes h1r7's target cells (the same
h1-test rows) with ALL 18,400 h1-test controls, and the access scores h1's own emitted stats against it (`--stats-map
h1test`, the 18 files of the h1-test access, sha-identical; 219 s, 6.19 GB). h1 → h1t7 changes only the truth's cells
per target; h1t7 → h1r7 changes the control count (the prediction side's control reference and key, and the truth's
DE controls). `work/cmp553_h1_h1t7_h1r7.md`:

| arm | h1 (400 cells, 18,400 controls) | h1t7 (61–400 cells, 18,400 controls) | h1r7 (61–400 cells, 2,095 controls) |
|---|---:|---:|---:|
| `lcap-kept` (CP) | +0.00505 (6.1 band sd) | +0.00397 (7.4) | +0.00270 (3.7) |
| `cpm+c4h` (hook) | +0.00019 (0.3) | +0.00008 (0.1) | +0.00045 (1.0) |
| `lcap-kept+c4h` (CPH) | +0.00647 (8.5) | +0.00447 (4.4) | +0.00219 (2.9) |
| CPH − CP (the hook on top of CP) | +0.0014 | +0.0005 | −0.0005 |
| `cpm-n1600` (K 2,000) | −0.00364 (−3.9) | −0.00114 (−0.8) | −0.00265 (−2.1) |
| `lcap-kept-n1600+c4h` (CPHK20) | +0.00020 (0.1) | +0.00024 (0.2) | −0.00249 (−2.1) |

The CP / CPH shrinkage splits roughly evenly between the two halves of the shape (CP: −0.0011 from the cell counts,
−0.0013 from the controls; CPH: −0.0020 and −0.0023; each surface in its own anchor units, but the raw members agree:
CP's jac gain +0.0091 → +0.0049 → +0.0040, reach +0.0096 → +0.0100 → +0.0070). The hook's increment on top of CP does
not survive the shape (+0.0014 → +0.0005 → −0.0005, inside a band of ~0.0008). So **a floored R7b access cannot separate
CPH from CP**: a null hook increment there is the expected reading, not evidence against the hook. (Whether it can see
CP's own gain is settled by the five draws below: at R7b's shape CP's mean gain is ~1 band — R7b can veto a large
harm, not confirm the gain.)

**Switch-on-recipe pairs (`--pairs`, the same machinery with the reference arm swapped; the report phase re-run on the
three accesses, their vs-`cpm` rows unchanged bit for bit).** The reading rule judges switches on top of a recipe, so
the access also reports A vs B for any two emitted arms:

| pair | h1 | h1t7 | h1r7 | raw nmae / fid / reach / jac on h1r7 |
|---|---:|---:|---:|---|
| `lcap-kept+c4h` vs `lcap-kept` (the hook on top of CP) | +0.00142 (+4.3 band sd; band 0.00033) | +0.00050 (+0.6) | −0.00051 (−1.3; band 0.00038) | +0.0000 / −0.0003 / −0.0022 / +0.0001 |
| `lcap-kept-n1600+c4h` vs `lcap-kept+c4h` (K 2,000 on top of CPH) | −0.00626 (−6.2) | −0.00422 (−3.3) | −0.00468 (−6.4; band 0.00073) | −0.0017 / −0.0012 / −0.0251 / +0.0008 |
| `lcap-kept+c4h-s` vs `lcap-kept+c4h` | | | +0.00052 (+0.6) | +0.0014 / +0.0001 / +0.0047 / +0.0000 |

A paired difference between two block arms that share their seeds has a small band (0.0003–0.0004 for the hook on CP).
On this first draw the hook's increment read −1 band sd and K-on-CPH −6.4 band sd; the five draws below put K-on-CPH's
mean at −0.0017 with a between-draw sd of 0.0023, so a K-size harm of 537's kind reads at about −1 to −2 of R7b's full
band — R7b can catch it only if Calu-3 shows it at least as strongly as H1, and only if the floor delivers the K-size
switch (5,000 does, 556 of 600 calls; 8,000 squeezes it to 248: the floor sweep below).

**A second draw breaks the rehearsal's reading — the within-surface band misses the control sample.** `h1r7s2` is the
same shape drawn again (seed 55308: another 2,095 controls, other cells, another assignment of R7b's counts to targets;
3,845 of its 13,034 cells shared with h1r7 by chance), with its own prediction side, 8 arms x 3 seeds and the same access;
`h1t7s2` is its ablation (draw 2's target cells with all 18,400 controls, h1's stats). `work/cmp553_draws12.md`:

| arm | h1t7 (draw 1 cells, all controls) | h1t7s2 (draw 2 cells, all controls) | h1r7 (draw 1, 2,095 controls) | h1r7s2 (draw 2, 2,095 controls) |
|---|---:|---:|---:|---:|
| `lcap-kept` (CP) | +0.00397 (band 0.00053) | +0.00318 (0.00094) | +0.00270 (0.00073) | **−0.00068** (0.00134) |
| `cpm+c4h` (hook) | +0.00008 (0.00079) | +0.00003 (0.00104) | +0.00045 (0.00045) | −0.00058 (0.00054) |
| `lcap-kept+c4h` (CPH) | +0.00447 (0.00102) | +0.00382 (0.00122) | +0.00219 (0.00076) | **−0.00077** (0.00136) |
| `cpm-n1600` (K 2,000) | −0.00114 (0.00138) | −0.00265 (0.00090) | −0.00265 (0.00128) | −0.00366 (0.00069) |
| `lcap-kept-n1600+c4h` (CPHK20) | +0.00024 (0.00139) | +0.00113 (0.00170) | −0.00249 (0.00118) | −0.00203 (0.00146) |
| pair: CPH vs CP (the hook on CP) | +0.00050 (0.00087) | +0.00064 (0.00098) | −0.00051 (0.00038) | −0.00009 (0.00085) |
| pair: CPHK20 vs CPH (K on CPH) | −0.00422 (0.00129) | −0.00269 (0.00202) | −0.00468 (0.00073) | −0.00126 (0.00119) |

With all controls the two cell draws agree within their bands (CP +0.0040 / +0.0032, CPH +0.0045 / +0.0038): the
split-half band does capture the truth cells' sampling. With 2,095 controls the same two draws disagree by 0.0034 (CP)
and 0.0030 (CPH) — 2.5–4.5 of their own band sds — and K-on-CPH moves from −6.4 to −1.1 band sd. So the dominant
uncertainty at R7b's shape is **the control sample**, which the access's band cannot see: the split-half re-draw halves
target cells against the full, fixed control block, while a different control sample changes both the prediction side
(control reference, call-propensity key, the emitter's control block) and the truth's DE. On R7b the control sample is a
single draw of Calu-3 controls; a reading's real uncertainty is larger than its band by a factor this rehearsal measures:

**Five draws (`draws553.py`, `work/draws553.json`).** Three more independent draws of the same shape (h1r7s3–s5, seeds
55309–55311; the 6 core arms x 3 seeds each), every access also running the cells + controls re-draw (`--ctrl-splits
5`, truth side only), each with its all-controls ablation (h1t7s3–s5). Per effect: the spread of the reading over the
five draws against the access's mean within-surface band:

| effect | mean over 5 draws | between-draw sd | band, cells (the mandate's) | ratio | band, cells + controls | ratio (90 % interval) | all-controls ablation: ratio to its band |
|---|---:|---:|---:|---:|---:|---|---:|
| CP (`lcap-kept` vs `cpm`) | +0.00143 | 0.00188 | 0.00125 | 1.50 | 0.00146 | 1.28 (0.83–3.05) | 1.08 |
| hook (`cpm+c4h` vs `cpm`) | +0.00002 | 0.00037 | 0.00054 | 0.70 | 0.00092 | 0.41 (0.26–0.97) | 0.75 |
| CPH (`lcap-kept+c4h` vs `cpm`) | +0.00152 | 0.00185 | 0.00136 | 1.36 | 0.00172 | 1.08 (0.70–2.56) | 1.08 |
| K 2,000 (`cpm-n1600` vs `cpm`) | −0.00155 | 0.00158 | 0.00110 | 1.44 | 0.00116 | 1.37 (0.89–3.25) | 1.48 |
| CPHK20 vs `cpm` | −0.00013 | 0.00305 | 0.00149 | 2.04 | 0.00205 | 1.48 (0.96–3.52) | 1.50 |
| hook on CP (CPH vs CP) | +0.00009 | 0.00043 | 0.00051 | 0.84 | 0.00091 | 0.47 (0.31–1.13) | 0.58 |
| K on CPH (CPHK20 vs CPH) | −0.00165 | 0.00232 | 0.00119 | 1.95 | 0.00135 | 1.72 (1.12–4.08) | 1.25 |

Reading: (1) for the call-set switches the mandate's band understates how far an R7b-shaped reading moves by 1.4–2.0x
(geometric mean 1.6); the truth-side controls re-draw recovers part (1.1–1.7x, geometric mean 1.4); (2) with all
18,400 controls the CP / CPH bands are right (1.08x) — the excess there is the control sample's — while the K-size
effects keep ~1.25–1.5x even with all controls (the cells' assignment to targets moves them beyond the split-half's
reach); (3) the hook's band is already conservative. (4) The five draws' mean is the R7b-shape expectation: CP +0.0014,
CPH +0.0015, K2000 −0.0016 (H1: +0.0051, +0.0065, −0.0036) — one R7b reading can sit 2 sd from it either way (draw 2:
CP −0.0007). The remaining control component (the prediction side's own control sample) is measured by control
halves:

**The control re-draw band, validated on the rehearsal (`prep553.py chalf`, `band553.py`).** Split a surface's controls
into two random halves A / B (3 splits); each half is a surface of its own — its prediction side (control reference,
call-propensity key, the emitter's control block and classes) rebuilt from that half, 537's 6 core arms emitted again
(seed 0), and the access run against the same target cells with that half of the controls. A and B are independent
draws of half the controls on BOTH sides, so an effect's difference between them, d, has E[d²] = 4 sigma_c² + 2 seed
variance (variance ∝ 1 / controls): sigma_c = sqrt(mean d² − 2 seed_sd²) / 2 is the control sample's share of the full
surface's uncertainty, and band_total = sqrt(band_sd² + sigma_c²). On h1r7 (`work/band553_h1r7.json`; each half access
~100 s, 2.9 GB):

| effect | h1r7 effect | band_sd (seeds + cells) | sigma_c (controls, 3 splits) | band_total | between-draw sd (5 draws) | ratio |
|---|---:|---:|---:|---:|---:|---:|
| CP | +0.00270 | 0.00073 | 0.00095 | 0.00120 | 0.00188 | 1.57 |
| hook | +0.00045 | 0.00045 | 0.00077 | 0.00089 | 0.00037 | 0.42 |
| CPH | +0.00219 | 0.00076 | 0.00103 | 0.00128 | 0.00185 | 1.44 |
| K 2,000 | −0.00265 | 0.00128 | 0.00122 | 0.00177 | 0.00158 | 0.90 |
| CPHK20 vs `cpm` | −0.00249 | 0.00118 | 0.00067 | 0.00136 | 0.00305 | 2.24 |
| hook on CP | −0.00051 | 0.00038 | 0.00041 | 0.00056 | 0.00043 | 0.78 |
| K on CPH | −0.00468 | 0.00073 | 0.00056 | 0.00092 | 0.00232 | 2.52 |

The control component is as large as the whole of the mandate's band (sigma_c 0.0006–0.0012 against band_sd
0.0004–0.0013), and adding it brings the band to the between-draw spread for CP, CPH, K 2,000 and the hook (ratios
0.4–1.6; with 3 and 4 degrees of freedom a ratio's 90 % interval is roughly x0.45–x2.5, so these are consistent with 1).
The K-size effects on top of a recipe keep 2.2–2.5x: their between-draw spread also carries the draws' different
assignment of cell counts to targets (the all-controls ablation already shows 1.25–1.5x there), a source a real surface
does not have — R7b's counts per target are what they are. **So the right band for R7b is its own band_total**: the
access's band plus the control re-draw from R7b's own control halves. R7b's 10 halves (5 splits of its 2,378 controls,
1,189 each — 5 splits rather than the rehearsal's 3 so that sigma_c has 5 degrees of freedom) are built and their
prediction side emitted before the access (controls only): all 10 run in the production regime (106 targets; prior dose
≥ 0.93; MDE / hook dose median 0.85–0.96; missed calls median 116–202 per target against 84 on the full surface, as
halving the controls does on the rehearsal: 47 → 84). The owner's access runs them with `--halves 0,1,2,3,4` (one
command; section "How to run").

## R7b: the production emitter cannot emit this surface as built — two depth constants break; a depth floor restores most of it

**Finding (emission side only; nothing read from the truth).** Every R7b prediction-side input was built (audit below).
The first R7b emission (`cpm` s0) stopped inside the production emitter: `rc_emit_k.py:1273 AssertionError: rc_hybrid_k
did not meet the per-cell totals`. Two depth-dependent mechanics of the K-class emitter fail on R7b as built:

1. **The class layout.** Production's classes are the context's control-library q25 / q50 / q75. R7b's are 2,816 /
   11,418 / 20,379 (q75/q25 **7.24**; A / B / C / H1: 1.95 / 1.57 / 1.86 / 1.58). The design then over-fills the low
   class (residual −550 … −1,413 counts per block in class 0, the matching surplus in class 2), the sub-gate repair stalls
   and the within-class rebalance has nowhere to go: **107 / 107 targets fail.**
2. **The steering margin.** The block-prior / MDE steering is feasible only if every class keeps a FIXED 6,000-count
   margin (`margin_ = 6000` in `rc_hybrid_block_k`), which at dose 0 already needs each class's share of the sub-gate
   (< 5 CPM) mass to reach 6,000 counts. H1: ~111 sub-gate counts per cell → ~14,700 per class. R7b: ~43 per cell →
   median 5,712 per class (min 4,557). So even when the placement succeeds, the CLEAN4 block prior's dose is λ = 0 —
   and with it the coarse MDE and the CLEAN4 hook (which rides on the MDE add) are not delivered.

Cross-tests (`diag553.py`: the same per-target calls, failures caught, the placement and the steering dose recorded):

| surface (controls) | emitter classes | q75/q25 | placement failures | prior dose λ (placed targets) | MDE dose λ, which carries the hook (q10 / median; targets < 0.9) |
|---|---|---:|---:|---|---|
| R7b as built (3,536) | production quartiles 2,816 / 11,418 / 20,379 | 7.24 | **107 / 107** | — | — |
| R7b as built (3,536) | `flex`: R7b median x H1's relative levels, 9,093 / 11,418 / 14,325 | 1.58 | 0 / 107 | **0 on 106 / 107** | not applied (no MDE, no hook) |
| R7b, controls ≥ 8,000 counts (2,095) | production quartiles 13,544 / 18,777 / 24,970 | 1.84 | 1 / 107 (BST2) | **1.0 on all 106** | **1.0 / 1.0; 0 / 106** |
| R7b, controls ≥ 5,000 counts (2,378) | production quartiles 11,265 / 17,310 / 23,929 | 2.12 | 1 / 107 (BST2) | **0.963–1.0 on all 106** | 0.77 / 0.90; 54 / 106 (the emission run without BST2, s0: 0.78 / 0.93; 44 / 106) |
| h1 (18,400) | production quartiles 16,079 / 20,190 / 25,330 | 1.58 | 0 / 85 | 0.93–1.0 (every KA run) | 1.0 / 1.0; 2 / 85 |
| h1 | median x √r spread, r = 2.5 | 2.5 | 0 / 85 | (not recorded) | (not recorded) |
| h1 | r = 3.5 | 3.5 | 0 / 85 | 0.93–1.0 | 0.96 / 1.0; 1 / 85 |
| h1 | r = 5 | 5.0 | 0 / 85 | **0 on 85 / 85 (silent)** | not applied |
| h1 | median x R7b's relative levels, 4,979 / 20,190 / 36,035 | 7.24 | **49 / 85** | — | — |

Why R7b's depth looks like this (control cells only, `libq553.py`; per-cell totals on the official axis): its controls
are two populations — a shallow mode (~1/3 of controls, mostly 0.5k–3k human-gene counts) and a broad Flex-like one;
above ~4k the density is flat. The cells are SARS-CoV-2-infected and viral reads are not among the 33,538 human features,
so heavily infected cells plausibly keep small human libraries (not tested). With the shallow mode removed the controls
look like A and C (floor 8,000: quartiles 13.5k / 18.8k / 25.0k).

| context | controls | min | q1 | q10 | q25 | q50 | q75 | q90 | max | q75/q25 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| A | 18400 | 3,275 | 4,470 | 9,027 | 13,986 | 20,109 | 27,234 | 34,328 | 52,420 | 1.95 |
| B | 18400 | 710 | 1,007 | 10,553 | 15,722 | 19,946 | 24,709 | 29,600 | 42,666 | 1.57 |
| C | 18400 | 3,447 | 4,685 | 9,972 | 14,414 | 20,034 | 26,742 | 33,875 | 50,965 | 1.86 |
| R7b | 3536 | 230 | 450 | 1,225 | 2,816 | 11,418 | 20,378 | 28,368 | 54,906 | 7.24 |
| R7b ≥ 5,000 | 2378 | 5,000 | | | 11,265 | 17,310 | 23,929 | | 54,906 | 2.12 |
| R7b ≥ 8,000 | 2095 | 8,000 | | 10,254 | 13,544 | 18,777 | 24,970 | 32,544 | 54,906 | 1.84 |
| H1 | 18400 | 7,327 | 8,377 | 12,774 | 16,079 | 20,190 | 25,330 | 31,039 | 124,278 | 1.58 |

BST2 fails the placement under every floor tested (5,000 / 6,000 / 8,000: residual −4,874 / −5,187 / −5,099 counts in
class 0, one cell 425–1,067 over) although its class sub-gate mass is ample (12.5k). It is a near-null target here
(effect-row norm rank 83 / 107, max |effect| 0.23, bulk weight clipped at 0.3; its own control CPM 82, unremarkable): an
emitter edge case on this surface, not a strong profile. It is excluded from emission and scoring with
`--exclude-targets BST2`; the panel (operator, carriers, bulk key) stays 107 targets, and the access treats its truth
cells as absent (fixture T6: identical to a file without them).

**The floor sweep: no floor restores everything; 5,000 keeps what the switches need.** `diag553.py r7bd<floor> <MODE> 0
--exclude-targets BST2` (run553's own stream: at 5,000 and 8,000 it reproduces the full runs' per-target call-set
records, K and doses exactly — the sweep's known answer), `cpm` (K 2,600) and `cpm-n1600` (K 2,000):

| floor | controls kept | q75/q25 | prior λ | MDE λ (carries the hook) q10 / median | K 2,600 calls missed per target q10 / median / q90 | K2600 − K2000 calls delivered, median (q10–q90) of 600 |
|---|---:|---:|---|---|---|---|
| H1 (no floor) | 18,400 | 1.58 | 0.93–1.0 | 1.0 / 1.0 | 0 / 31 / 63 | 566 (534–585) |
| h1r7 (H1, 2,095 controls) | 2,095 | 1.57 | 0.92–1.0 | 1.0 / 1.0 | 8 / 47 / 85 | 556 (527–581) |
| **R7b ≥ 5,000** | 2,378 (67 %) | 2.12 | 0.963–1.0 | 0.76 / 0.91 | 48 / 84 / 138 | **556 (522–570)** |
| R7b ≥ 6,000 | 2,281 (65 %) | 2.00 | 1.0 | 1.0 / 1.0 | 78 / 169 / 415 | 497 (272–558) |
| R7b ≥ 7,000 | 2,184 (62 %) | 1.90 | 1.0 | 1.0 / 1.0 | 107 / 323 / 612 | 362 (143–544) |
| R7b ≥ 8,000 | 2,095 (59 %) | 1.84 | 1.0 | 1.0 / 1.0 | 167 / 451 / 695 | 248 (124–493) |

Raising the floor buys the MDE / hook dose (full from 6,000) and costs the call set: the emitter misses a growing share
of its 2,600 prescribed calls (3 % at 5,000, 6.5 % at 6,000, 17 % at 8,000; H1 1 %), and because a smaller prescribed
set (K 2,000) is still delivered almost fully, the K-size switch is squeezed — to 41 % of its H1 size at 8,000. The
loss is not the control count (h1r7, H1 cut to 2,095 controls, misses 47); its mechanism in the emitter is not
diagnosed here. The October switches act through the call set (CP chooses which genes are called, K how many), and the
rehearsal below shows the hook's increment is unresolvable at R7b's shape anyway; so the floor that keeps the call set is
the one to use.

**Verdict.** R7b as built **cannot be priced faithfully with the production emitter**: under production's rule it does not
emit at all, and under any class rule its shallow third makes the steering infeasible, so the incumbent's DE construction
(CLEAN4 block prior, coarse MDE) and two of the October switches' mechanisms (the hook; the block's `kept` steering) would
be absent — a regime D/E/F (Flex, ~20k counts per cell) will not be in. **A human-library floor restores the production
emitter unmodified, with one trade-off and one exclusion**: at 5,000 (67 % of controls kept) the placement runs on 106 /
107 targets, the prior's dose is ~1, the call set and the K-size switch are delivered as on H1, and the MDE add (with
the hook riding on it) is delivered at a median 0.91 (q10 0.76; sub-gate mass 81.5 counts per cell against H1's 111, so
the fixed 6,000-count margin still bites); higher floors give the full MDE dose and lose call-set fidelity (sweep
above). The floor is an amendment of RESERVE §3's construction ("no cell downsampling; depth as measured"), so it is the
owner's (552's) decision: it would be applied by `build_r7b.py` to controls AND target cells before the access, at a
value fixed now. (RESERVE §6 clears the CONTROL cells' library sizes; selecting TARGET cells by their per-cell total
reads a scalar of perturbed cells' expression — no per-gene or per-target statistic, but the owner's to rule on; I did
not compute it.) The value I suggested is **5,000** (6,000 if the owner weighs the hook's full dose above the call set);
**552 adopted 5,000 (RESERVE §7, Amendment 3)**, ruled the per-cell library read admissible as QC, and built the floored
surface: 16,381 of the target cells kept (59–556 per target, median 127), 2,378 controls identical to the prediction
side's. Its cost, as 552 states it: the floor selects cells by a quantity the infection screen itself moves, so
per-target truths are those of less-infected cells (closer to D/E/F's uninfected state, but a selection that can differ
by target). The harness reads the floored truth as an ordinary sealed surface (`r7bd5000` picks up
`552-r7-reserve/work/r7b_floor5000/`; the access refuses to floor a truth itself).

**A production robustness fact for October (H1 sweep and floors above).** The K-class emitter has two depth/shape cliffs
and one slope: the steering dose collapses to 0 WITHOUT any error somewhere between class spread 3.5 and 5 (h1: λ ≥ 0.93
at 3.5, λ = 0 on 85 / 85 at 5); the placement crashes by 7.24 (49 / 85); and the MDE / hook dose erodes quietly when the
sub-gate mass per cell is low (R7b ≥ 5,000, 81.5 counts per cell: median 0.90–0.95, q10 ~0.78; H1, 111: 1.0) — the
capacity margin that zeroed context A's plain MDE (429a's rc dump, `make_submission_pm_cd4.py:711`), which the `coarse`
form used here shrinks but does not remove. A / B / C sit at 1.57–1.95 spread, and D/E/F are Flex contexts from the same
Arc run, so the risk is low; but a silent λ = 0 build
would ship without its CLEAN4 prior, MDE and hook. On R7b's floors a fourth quantity moves too — the call set's
delivery (median missed calls 84 → 451 from 5,000 to 8,000; H1 31), by a mechanism not diagnosed here. Cheap guards for
the day: q75/q25 of each context's control libraries (flag above ~3), and per context the median of the emitter's
per-target `prior_lambda` and `mde_lambda` (flag below ~0.9) and `call_set_missed` (flag above ~100) diagnostics
(`rc_hybrid_block_k`'s record, the fields this harness reads).

## Item 5: the floored candidates' emissions (emission side only; no truth read)

`emitdiag553.py` reads only the predicted blocks' own statistics against the CONTROL cells (RESERVE §6): per arm the
predicted call set per target (q < 0.05, own gene excluded), the emitter's call-set record, the ballast K, the steering
doses, and the two identity checks a DE-only knob must pass against `cpm` at the same seed (per-gene block sums = the
pds / mse inputs, and pass-1 ballast K). All 537 arms x seeds 0–2 on `r7bd5000` (BST2 excluded): 24 runs, rc 0, 332 s
median (342 max), peak RSS ≤ 2.14 GB; `access553.py r7bd5000 --dry-run` hashes the 24 stats files, plans 264 scorings
and verifies 552's floored build against them.

R7b ≥ 5,000 (`work/emitdiag553_r7bd5000-xBST2.json`, means over seeds 0–2):

| arm | called / target (median) | up frac | prescribed / delivered | missed | unsuppressible | hook signs | calls new vs `cpm` | sums == `cpm` | pass-1 K == `cpm` | prior λ min / median | MDE λ min / median |
|---|---:|---:|---|---:|---:|---:|---:|---|---|---|---|
| `cpm` | 2510 | 0.369 | 2600 / 2519 | 90.1 | 1.1 | | | | | 0.963 / 1.000 | 0.025 / 0.907 |
| `lcap-kept` | 2518 | 0.398 | 2600 / 2531 | 77.1 | 1.1 | | 422 | True | True | 0.969 / 1.000 | 0.025 / 0.919 |
| `cpm+c4h` | 2510 | 0.368 | 2600 / 2518 | 90.6 | 1.0 | 2548 | 5 | True | True | 0.963 / 1.000 | 0.580 / 0.926 |
| `lcap-kept+c4h` | 2518 | 0.393 | 2600 / 2530 | 78.1 | 1.0 | 2548 | 418 | True | True | 0.969 / 1.000 | 0.653 / 0.946 |
| `cpm-n1600` | 1959 | 0.376 | 2000 / 1962 | 38.3 | 0.0 | | 37 | True | True | 0.975 / 1.000 | 0.581 / 0.962 |
| `lcap-kept-n1600+c4h` | 1967 | 0.407 | 2000 / 1971 | 29.5 | 0.0 | 2548 | 272 | True | True | 0.981 / 1.000 | 0.748 / 0.978 |
| `cpm-n1200` | 1583 | 0.381 | 1600 / 1584 | 16.5 | 0.0 | | 40 | True | True | 0.993 / 1.000 | 0.719 / 0.994 |
| `lcap-kept+c4h-s` | 2518 | 0.396 | 2600 / 2532 | 77.0 | 1.1 | 2548 | 425 | True | True | 0.969 / 1.000 | 0.535 / 0.920 |

H1 for comparison (this directory's KA runs, seed 0, `work/emitdiag553_h1.json`):

| arm | called / target (median) | up frac | prescribed / delivered | missed | unsuppressible | hook signs | calls new vs `cpm` | sums == `cpm` | pass-1 K == `cpm` | prior λ min / median | MDE λ min / median |
|---|---:|---:|---|---:|---:|---:|---:|---|---|---|---|
| `cpm` | 2566 | 0.352 | 2600 / 2569 | 32.4 | 5.2 | | | | | 0.926 / 1.000 | 0.012 / 1.000 |
| `lcap-kept` | 2567 | 0.389 | 2600 / 2570 | 31.2 | 23.7 | | 611 | True | True | 0.945 / 1.000 | 0.011 / 1.000 |
| `cpm+c4h` | 2567 | 0.352 | 2600 / 2569 | 31.8 | 3.1 | 2550 | 5 | True | True | 0.926 / 1.000 | 0.012 / 1.000 |
| `lcap-kept+c4h` | 2573 | 0.383 | 2600 / 2574 | 26.0 | 6.7 | 2550 | 616 | True | True | 0.945 / 1.000 | 0.159 / 1.000 |
| `cpm-n1600` | 1998 | 0.341 | 2000 / 1999 | 2.4 | 24.1 | | 22 | True | True | 0.942 / 1.000 | 0.011 / 1.000 |

Reading: every arm's mechanism is present on the floored surface and acts as on H1 — the DE-only knobs leave the block
sums and pass-1 K bit-identical to `cpm` (so pds / mse can only move through the zs pass), the call-propensity block
swaps ~420 calls per target (H1 ~610), the hook signs 2,548 genes (H1 2,550) and moves ~5 calls (H1 5), K2000 / K1600
shrink the call set as prescribed. The ballast K histogram (bulk side) is H1-like: 14 / 106 targets at the cap K = 16
(H1 10 / 85). Two differences: the call set misses more prescribed calls (90 vs 32 per target for `cpm`) and the MDE dose
(which carries the hook) is lower — median 0.91–0.95 on the 2,600-call arms against 1.0 on H1 — the 5,000 floor's
partial restoration (cross-test table above).

R7b ≥ 8,000, all arms at seed 0 (`work/emitdiag553_r7bd8000-xBST2.json`; seeds 1–2 held once the sweep showed the call-set
loss; `cpm` and `lcap-kept` s1 done): every mechanism is present and both doses are 1.0, but `cpm` delivers 2,150 of its 2,600 calls
(missed 447 per target) while `cpm-n1600` still delivers 1,910 of 2,000 — the K-size switch squeezed to ~240 calls.

| arm | called / target (median) | up frac | prescribed / delivered | missed | unsuppressible | hook signs | calls new vs `cpm` | sums == `cpm` | pass-1 K == `cpm` | prior λ min / median | MDE λ min / median |
|---|---:|---:|---|---:|---:|---:|---:|---|---|---|---|
| `cpm` | 2114 | 0.380 | 2600 / 2150 | 447.2 | 0.1 |  |  |  |  | 1.000 / 1.000 | 1.000 / 1.000 |
| `lcap-kept` | 2122 | 0.404 | 2600 / 2186 | 422.5 | 0.0 |  | 378 | True | True | 1.000 / 1.000 | 1.000 / 1.000 |
| `cpm+c4h` | 2108 | 0.380 | 2600 / 2152 | 447.8 | 0.0 | 2546 | 16 | True | True | 1.000 / 1.000 | 1.000 / 1.000 |
| `lcap-kept+c4h` | 2114 | 0.397 | 2600 / 2182 | 423.6 | 0.0 | 2546 | 376 | True | True | 1.000 / 1.000 | 1.000 / 1.000 |
| `cpm-n1600` | 1892 | 0.381 | 2000 / 1910 | 116.9 | 0.0 |  | 264 | True | True | 1.000 / 1.000 | 1.000 / 1.000 |
| `lcap-kept-n1600+c4h` | 1923 | 0.407 | 2000 / 1936 | 87.7 | 0.0 | 2546 | 490 | True | True | 1.000 / 1.000 | 1.000 / 1.000 |
| `cpm-n1200` | 1550 | 0.386 | 1600 / 1553 | 45.5 | 0.0 |  | 236 | True | True | 1.000 / 1.000 | 1.000 / 1.000 |
| `lcap-kept+c4h-s` | 2120 | 0.402 | 2600 / 2180 | 423.5 | 0.0 | 2546 | 386 | True | True | 1.000 / 1.000 | 1.000 / 1.000 |

The `r7bd6000` alternative's emissions were held once 552 adopted 5,000 (its prediction side exists: control file,
control reference, call-propensity key).

## R7b prediction-side audit (`work/r7b/audit553.json`; controls + public donors only)

| item | value |
|---|---|
| surface (552's `build_r7b.py`, 2026-10-01 00:27Z) | 107 targets, 23,860 target cells (103–826 per target, median 187) + 3,536 controls; 18,333 of 18,533 official genes mapped (200 all-zero); every bundle file's sha256 == `build_r7b.json` (the truth by sha256 only: `bd35f450…`) |
| control reference (`work/r7b/controlref553.npz` sha `a0fb10d3…`) | 3,536 controls; 5-CPM universe 9,923 genes (H1 10,779); control libraries min 230 / q25 2,816 / median 11,418 / q75 20,379 / max 54,906; 3,868 all-zero genes; control jackknife 3.165 |
| donors (`k562_gwps`, `orion_hct116`, `orion_hek293t`) ↔ targets | 92 / 105 / 107 of 107; 90 targets carried by all three, 17 by two; 0 K562-essential targets |
| PerturbMatch K562 replacement | 91 targets modified (89 with 3 votes), 14 supported by neither, 1 + 1 one-sided |
| CD4 banks (047's `extract_cd4.py`, relaxed thresholds, `work/r7b/cd4_compact_r7b`) | 98 of 107 targets in any condition, 96 in all three; 10,113 official genes |
| whitening operator (`build_union_whitener.py` as 323 ran it; 330 s, 13.66 GB, under the mandate's 27 GiB gate) | `work/r7b/whitener_r7b_pruned.npz` sha256 **`d846b222…`** (full: d846b22273fb46d91ea1c28d4811f36127c36f613c1645596dac86930e1c84f5); axis sha `a6ba78e5…` over 18,135 columns (the panel's 107 genes off the axis); basis 29,394 rows (7,784 + 10,807 + 10,803), 0 panel targets among the pruned rows |
| carriers (cap15 under 323's model_env) | final sha256 `cec09402…`; effects.npz `89e4d340…`; CD4 k-5 deflation on 98 carried rows (26.1 % of the centred energy) |
| 410b bulk key (the panel's, live = cache bit-equal: KA-B) | 107 carried; median key 0.00275; weights clipped to [0.3, 3]; energy constant 0.656 |
| CLEAN4 q80 block prior on R7b's axis | covers 61.8 % of the 9,923 gated genes (53.5 % down among covered); H1 57.2 % of 10,779; R7b ≥ 8,000 64.2 % of 9,488 |
| CLEAN4 hook (τ 0.05) | 2,558 signed gated genes (1,497 down); H1: 2,550 |
| call-propensity key (528's `callblock528.py`, the October form) | excess-dropout median 0.00974 → scale 0.901 to H1's; the key's block swaps ~518 genes vs the CPM block (no-head approximation) |
| forbidden sources | none: the 12 input paths and this directory's code never name the Sunshine study or a Calu-3 bank; the CLEAN4 banks are K562-GWPS, K562-essential, Nadig Jurkat, Nadig HepG2 |

## How to run the access (the owner, after the claim is preregistered; never before the 10-07 reads)

What it needs, all in place for `r7bd5000` (552's floored surface, Amendment 3): (a) the claim file (its sha256 goes into
the plan and every phase re-verifies it); (b) the arms' stats — 537's 8 arms x seeds 0–2 on `r7bd5000`, BST2 excluded
(emitting reads no truth); (c) 552's floored truth `experiments/552-r7-reserve/work/r7b_floor5000/` — the dry run checks
its 5 outputs by sha256, its 106 targets and its 2,378 controls against the prediction side; (d) for the control
re-draw band, the 10 control halves `r7bd5000c{0..4}{A,B}` (1,189 controls each; prediction side + the 6 core arms at
seed 0, emitted before the access); (e) nothing else running heavy.

```bash
cd ~/Code/virtual-cell-challenge-2026
ARMS=lcap-kept,cpm+c4h,lcap-kept+c4h,cpm-n1600,lcap-kept-n1600+c4h,cpm-n1200,lcap-kept+c4h-s
PAIRS=lcap-kept+c4h:lcap-kept,lcap-kept-n1600+c4h:lcap-kept+c4h,lcap-kept+c4h-s:lcap-kept+c4h,cpm-n1200:cpm
# 0. readiness, never opens the truth (stats hashed, 552's build by sha256, controls equal, control halves ready)
uv run python experiments/553-r7b-harness/access553.py r7bd5000 --dry-run --exclude-targets BST2 --arms $ARMS \
    --seeds 0,1,2 --halves 0,1,2,3,4
# 1. the access + the cells-and-controls re-draw + board retrodiction + the control-half band, one command
/usr/bin/time -l uv run python experiments/553-r7b-harness/access553.py r7bd5000 --owner-go <claim file> \
    --exclude-targets BST2 --arms $ARMS --seeds 0,1,2 --pairs $PAIRS --ctrl-splits 5 \
    --retro work/access553_rehearsal5 --halves 0,1,2,3,4
```

Outputs: `work/access553_r7bd5000-xBST2/access553.{json,md}` — per arm (vs `cpm`) and pair: raw and scaled member effects,
seed sd, split-half sd, `band_sd` (the mandate's band: seeds + target-cell re-draws), `band_sd_ctrl` (cells + controls,
truth side), the target-bootstrap CI, the retrodiction table (each raw member's effect over the shape-matched rehearsal's
mean, `work/access553_rehearsal5`, five h1r7 draws; use `--retro work/access553_h1` for H1-full); the 10 half accesses
`work/access553_r7bd5000c<k><h>-xBST2/`; and **`work/band553_r7bd5000-xBST2.json`: per effect `band_total` = sqrt(band_sd²
+ sigma_c²), the band the rehearsal validates** (section "the band"). Switch-on-recipe pairs need both arms in `--arms`;
a finished access can be re-reported (`--phase report` with the same arguments, new `--pairs` / `--retro`; the report
reads only saved phase outputs, never the truth; on a sealed surface it still re-verifies `--owner-go`). Option values
that begin with a dash need `=` (`--tag=-flex`).

Cost (measured on the rehearsal, which is ~10 % smaller): the same one command with 3 splits took 864 s at a 3.75 GB
peak (KA-E2E: main access with csplit ~6 min, each half access ~100 s at 2.9 GB); with 5 splits on R7b ≥ 5,000 expect
~25 min and < 4 GB. Emitting
another arm first (~5.5 min, ~2 GB with the cached bulk key, reads no truth): `uv run python experiments/553-r7b-harness/
run553.py r7bd5000 <MODE> base <seed> --bulk ps --mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse --zs
dn_k250_zs --bulk-cache use --exclude-targets BST2` (and on each half `r7bd5000c<k><h>` at seed 0 for the band).

Other surfaces: `r7b` as built cannot be emitted by the production emitter and `r7b --class-rule flex` emits without the
steering — neither is a faithful pricing; `r7bd8000` (seed 0 pre-emitted) squeezes the K-size switch (the floor sweep).

## What was generalised (file:line) — and what was not touched

The per-target loop is `experiments/496-uncalled-sign-hook/run496.py:main` itself (imported, unmodified): profile /
spread / rank axis / 410b ballast K / carrier block (the production emitter `src/vcc2026/rc_emit_k.py`, sha asserted
`0f56209d…` at `run553.py:106-107`) / 426's MDE add / 423-429's zs pass / `evaluate323.block_stats`. What run553 replaces,
one input each (the same points run537 patched, plus the four frozen-H1 inputs); line numbers are of the final files:

| H1-frozen input (where the H1 instrument reads it) | surface input (553) | on h1 |
|---|---|---|
| `evaluate323.load_reference()` → 323's reference.npz (18,400 controls) | `run553.py:164-174` (`surface_reference`, from `CR = lib553.load_controlref`, `:142`; h1: 323's arrays; r7b / r7bd / h1r7: `controlref553.npz` built by `lib553.controlref`, `lib553.py:354`, the same `rc_emit.control_block_from_counts` / `exact_de323.control_block` / jackknife calls as evaluate323.reference, any control count) | identical arrays (KA-C) |
| run400's `CLASSES3` constant (16078.75 / 20190.5 / 25330.5) | `run553.py:156-161` → `lib553.library_classes` (`lib553.py:386`): the context's control-library q25/q50/q75, production's `de_sign_classes` rule; `--class-rule flex` → `lib553.flex_classes` (`:396`), a deviation for R7b as built | == the constant (asserted, `run553.py:158`); flex == production on h1 (KA-flex) |
| `run400.native_cap15` → 323's effects.npz | `run553.py:177-189` (`surface_cap15`, same compose-sha check; with `--exclude-targets` the kept rows) → `effects553.py` (323's stage / operator / compose for any panel) | 323's file (KA-E reproduces it) |
| `run496.ps_bulk` → `H1Op` (323's operator, 18,157 cols asserted) | `run553.py:192-262` (`SurfaceOp`, `BULK_FILE` `:206`, `ps_bulk` `:210`): the same code on the surface's operator + columns, the key and energy constant on the WHOLE panel when targets are excluded, cached with provenance (targets, carriers sha, operator sha) | 323's operator (KA-B: cached == live) |
| run537's `keys528.npz` `lcap` | `run553.py:275-282` → 528's `callblock528.py` artifact from the surface's controls (the October form) | keys528 (and KA-CB2 via `h1cb`) |
| run537's `_clean4_P_full` on h1-test's axis | `run553.py:330`, same code on the surface's axis | same |
| `evaluate323.score` inline | h1: unchanged; sealed: a stub returning null (`run553.py:266-269`); NaN-safe metrics json (`write_json_finite`, `:133`) | unchanged |
| `evaluate323.reference` / `score` / anchors323 (access side) | `lib553.truth_rows` (`:463`, the only truth reader, via `open_truth` `:235`; a control half keeps the truth's controls by name), `control_moments` (`:500`, checks the truth's controls == the prediction side's), `truth_reference` (`:513`), `score_stats` (`:547`, prediction counts 400 per target, `:561`), `anchors_base` / `replicate_halves` / `anchors_rep` / `anchors` (`:642` / `:694` / `:744` / `:769`, group by group), `scaled` (`:778`, cell-eval2's `score_one`) | KA-R / KA-S / KA-A |
| (new) the band beyond the mandate's | `access553.py:220` (`csplit`: cells + controls re-draw, truth side, control blocks on the full gene universe via `lib553.control_block_keep` `:332`), `:344` (`paired`: any arm vs any arm, `--pairs`), `:582-613` (`--halves`: the control-half accesses + `band553.py`); the half surfaces themselves are `lib553.surface_spec` `:110-133` (`CTRL_HALF_RE`) and `prep553.py chalf` | T8 (`control_block_keep` == `control_block`); KA-PAIRS |

Not generalised: `+mzg` (527's Arc-sourced hook) — its strong / weak split reads 467's H1 strength file and its source
excludes H1's targets; `run553.py` refuses `+mzg` off h1 (537 P12: October keeps CLEAN4 unless CPMZ beats its r 0.37
prediction by ~2×; if CPMZ enters the recipe, its production form — `le1` weak gate on count(|key| > 0.2) and 527's
artifact — has to be added first).

## Findings on the way (instrument facts)

- **cell-eval2's nmae dispatch is not bit-deterministic.** Six identical `score_stats` calls on one stats file (528 `lcap-kept`
  s1 zs) returned three distinct nmae values (…816 / …817 / …818, ±1.1e-16); the other five members were identical. So
  "bit for bit" is a property of the emitted STATISTICS (p, q, lfc, sums, jk), not of every member: this directory's
  known answers require bit-equal arrays and member |d| ≤ 1e-12, and report the exact |d|.
- **A process gate must read argv, not `ps comm`.** `ps -o comm` truncates to 16 characters, so a venv python shows as
  `~/Cod` and a `comm ~ python` test never fires; my first driver draft had that bug (fixed before any heavy job;
  537's own `alive` helper is args-based and correct). Chains are now matched as argv[0] = a shell AND its first non-flag
  argument = the chain script, which also keeps headless workers whose prompts name the scripts from matching.
- **macOS keeps freed large allocations in a process's RSS.** Freeing 1.5 GB of numpy arrays left `ps` RSS unchanged and
  `malloc_zone_pressure_relief` released 0 bytes; KA-C read 7.87 GB against a ~4.7 GB live peak. Under an RSS cap the live
  high-water mark is what counts, so the access was rebuilt to never hold the whole truth matrix twice: the truth is held
  once as its control rows + target rows; the control jackknife (whose kernel copies its input twice) runs first and its
  moments are reused by the anchors; scipy's sparse `.mean(axis=0)` (which multiplies every element by 1/n into a copy
  before summing) is reproduced bit for bit by a row-streamed replica (`lib553.csr_mean0`, equal to scipy on 3,000
  h1-test control rows and a 400-cell group).
- **The K-class emitter has several depth-dependent mechanics, not one.** (1) The class layout: production's quartile
  classes fail the per-cell placement when the control libraries' q75/q25 is large (R7b 7.24: 107 / 107 fail; h1 at
  R7b's relative levels: 49 / 85). (2) The block prior's dose: `rc_hybrid_block_k`'s steering keeps a FIXED 6,000-count
  margin per class, so the dose drops to 0 — silently — when a class's sub-gate (< 5 CPM) mass is too small (R7b under the
  flex rule; h1 at class spread 5). (3) The MDE add's dose (which carries the CLEAN4 hook) shares that budget: at R7b
  ≥ 5,000 (81.5 sub-gate counts per cell vs H1's 111) it is delivered at a median 0.90–0.95 and below 0.9 on 44–54 of
  106 targets, while the prior's dose stays ~1. (4) On R7b's floors the call set's delivery falls as the floor rises
  (missed calls per target, median: 84 / 169 / 323 / 451 at 5,000 / 6,000 / 7,000 / 8,000; H1 31; h1 cut to 2,095
  controls 47), which squeezes the K-size switch; mechanism not diagnosed. Placement errors are loud; the two doses and
  the call-set loss are not, and the emitter's own per-target record carries all three.
- **One R7b target crashes the production emitter at any floor tested (BST2).** Forced placements are normal (~5,000
  counts per block on R7b ≥ 5,000 / 8,000, ~9,000 on h1 at spread 3.5; the repair leaves |residual| ≤ ~42 per cell);
  after both repairs BST2's block still has one cell 425–1,021 counts over its total and its low class ~4,900–5,100 over
  in sum (the sub-gate repair stalls on over-full cells that hold too few sub-gate counts; the fallback rebalance moves
  counts only within a class). The failure is an AssertionError (loud), not a silent mis-emission.
- **argparse needs `--opt=value` for values that begin with a dash** (`emitdiag553.py r7bd5000 … --tag=-xBST2`;
  `--tag -xBST2` is read as two options). `run553.py` parses its own options and takes either form (tested on five
  argv shapes before the edit was written; KA-FINAL ran with `--tag=-kafinal`).
- **The cached 410b bulk key spares the operator load on H1**: peak RSS 5.2–5.5 GB instead of 7.0–7.5 GB, same bits
  (KA-FINAL; KA-B / KA-flex earlier).
- **A within-surface noise band can miss the dominant noise.** The mandate's band (emission seeds + split-half re-draws of
  the truth's target cells against the fixed control block) is exact on H1's shape for the call-propensity arms (all
  18,400 controls: between-draw / band 1.08) but understates an R7b-shaped reading's spread by 1.4–2.0x, because with ~2k
  controls the control SAMPLE — through the prediction side (control reference, call-propensity key, the emitter's
  control block) and the truth's DE — moves the paired effects as much as the cells do. Re-drawing the controls on the
  truth side only recovers part; re-drawing them on both sides (control halves, each re-emitted) closes the gap on the
  rehearsal. Two draws were not enough to see this reliably (one ratio's 90 % interval spanned 1x–30x); five were.
- **A rehearsal's between-draw spread is not the same thing as a real surface's uncertainty.** The rehearsal draws also
  re-assign R7b's cell counts to targets, which moves the K-size effects (the all-controls ablation keeps 1.25–1.5x for
  them); a real surface's counts per target are fixed, so its band (seeds + cells + controls) is the comparable quantity.

## Pricing any DE arm on a power-reduced H1 (`h1r7`-style) surface: how, what it costs, and what it is not (owner, 10-02)

Added by the owner at the coordinator's request (10-02). **What it is not: a board-shaped instrument.** A/B/C's board
references have exactly 400 cells per target and 18,400 controls, h1-test's own shape, so H1-full is already
shape-matched to the board. The h1r7 cut (R7b's approximate shape, built before the floor was chosen: 2,095 controls and 61–400 cells per target, median 110, against the floored R7b's 2,378 controls and 59–556 cells, median 127) LOWERS power.
That it gives the call-propensity block board-like ratios (five-draw mean jac 0.45×, fid 1.4×, reach −0.05× of H1-full,
against the board's 10-01 read of 0.48×, 1.95×, −0.38×) points at a power or regime effect, like 537's f-regimes (f 0.5 /
0.35 / cmp), not at shape. Use it as a REGIME instrument. It is validated only as far as it retrodicts known board reads
(the 10-01 block read; 446bj2z's zs, nmae −0.0022 raw; 496a's hook on 10-05; the K pair on 10-06). Here it matches jac,
under-reads fid's cost (1.4× vs 1.95×) and misses reach's sign. **Prices already on disk** (`work/draws553.json`; five
independent h1r7 draws, 3 seeds each, overall-equivalent at each surface's own anchors, between-draw sd in brackets):
CP +0.00143 (0.0019), CPH +0.00152 (0.0019), K 2,000 −0.00155 (0.0016), hook +0.00002 (0.0004), CPHK20 −0.00013 (0.0031),
hook on CP +0.00009 (0.0004), K on CPH −0.00165 (0.0023); one draw: K 1,600 −0.00296. Read them as "one regime's mean",
and use the between-draw sd as their noise, not the band.

**How, for an arbitrary arm (537's MODE grammar), all from h1-test (not a reserve; no R7b data beyond its count metadata):**
```bash
cd ~/Code/virtual-cell-challenge-2026; E=experiments/553-r7b-harness
# a new draw of the reduced-power surface (or reuse h1r7 / h1r7s3..s5): N controls, cells per target = R7b's counts x f
uv run python $E/prep553.py h1r7 h1r7s6 --seed 55312 --controls 2095 --cells-scale 0.59      # ~1 min, < 3 GB
uv run python $E/prep553.py controlref h1r7s6 && uv run python $E/prep553.py callblock h1r7s6  # prediction side, controls only
# emit each arm x seed (~5 min, ~2 GB each; the cached bulk key spares the operator load)
uv run python $E/run553.py h1r7s6 <MODE> base <seed> --bulk ps --mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 \
    --mde-coarse --zs dn_k250_zs --bulk-cache use
# score (unsealed surface: no --owner-go); --ctrl-splits re-draws cells + controls on the truth side
uv run python $E/access553.py h1r7s6 --arms cpm,<MODE>[,...] --seeds 0,1,2 --pairs <a>:<b>[,...] --ctrl-splits 5
```
Cost: per arm 3 seeds × ~5 min at ~2 GB; the access ~3–6 min at ≤ 3.8 GB; a new draw adds ~3 min. Run one draw per
(controls, cells-scale) setting and several draws per setting (`draws553.py` summarises them) before reading a sign: at
this power one draw sits up to ~2 sd from the mean (draw 2 read CP −0.0007 against a mean of +0.0014). The surface name
must start with `h1r7` (`lib553.surface_spec`). To vary only the cells, `prep553.py h1t7 build` keeps all 18,400
controls (the ablation).

**Reusing the five existing draws (after the 10-02 prune, `work/PRUNED.md`).** Each draw keeps its prediction side
(`controlref553.npz`, callblock) and the base arm `h1r7[sN]_cpm` at seeds 0–2. Its truth h5ad and the other arms' stats
were pruned. Rebuild the truth first, with the draw's own arguments (the defaults build a different surface): `prep553.py h1r7
build --controls 2095 --cells-scale 0.59` for h1r7 (seed 55307), `prep553.py h1r7 h1r7s<N> --seed <55306+N> --controls
2095 --cells-scale 0.59` for s2–s5 (seeds 55308–55311). The rebuild overwrites `work/<draw>/<draw>.json`, which is
committed: `git diff experiments/553-r7b-harness/work/<draw>/<draw>.json` must be empty, or the file is not the priced
draw (KA-E2E rebuilt seven files this way sha-identically). Then emit only the new arm and run the access
as above, so a five-draw price of one new arm costs 5 × 3 emissions plus 5 accesses, about 1.5–2 h serially.
