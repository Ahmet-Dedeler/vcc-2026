# 658b — 524dr2, the October dress rehearsal on the applied post-10-06 tree (PREREGISTRATION)

Council session [7772dc], assigned by the coordinator 10-04 06:18Z. Written 10-04 from 06:40Z and committed before 524dr2 is
built, before 524dr (10-06 s1) is scored, before 496a (10-05 s1) is read, and before tercilec is decided (10-06 step 9).
It carries over 524's PREREGISTRATION (bins, falsifiers, Σ rule) and amends it only where marked **[A]**.
**Reader:** `read658.py dr` (tests: `test658.py`, 23/23). With `--prev` omitted, the same reader is 524's own reader for 524dr
(653 R6 needed one, and none existed).

## 1. The arm
- **524dr2 = 524dr's procedure, run on the applied tree.** That is `recipe_dr524.sh`'s steps with the post-10-06 builder,
  emitter and def_day (APPLY DAY v2 steps 1–8 committed, every step-6 KA PASS):
  - resolve `october_recipe.py command --allow-pending --json` at build time and save it as `<w>/recipe524.json` in 524dr2's
    own work dir `<w>` (the name `gate524.py --work <w>` reads), with a copy of `docs/october-decisions.json`;
  - 180's relabelled A/B/C bundle; 453's day operator (66610be7);
  - a FRESH work dir (def_day treats an existing h5ad as built, 485 F1);
  - `def_day.sh` with exactly the resolved exports and args; relabel453; `vcc prep` dry run then real; novelty;
  - `gate524.py --work <w> --stem <s>` vs 453dr. Its SWITCH_KEYS already cover `DEF_DAY_SHARPEN_PT` (tercile or tercilec)
    and `DEF_DAY_MDE_SIGN`.
- **[A] Switches that can be new vs 524dr:**

  | switch | how it gets in | channel |
  |---|---|---|
  | `mde_sign` | 496a TRANSFERS on 10-05 and is DECIDED after 524dr's 10-05 02:00Z resolve. 524dr runs it OFF if 496a was not yet read | DE |
  | `sharpen_pt` → `tercilec:2.5,2.0,1.75` | step 9 decides it after C10b | bulk, 4 A/B/C targets |
  | `h1_dose` (653 R10) | only if 459dr reads ≥ +0.0015 before the resolve. The coordinator wires dayof438's H1 choice; then a `gate524.py` SWITCH_KEYS row and `read658.py --alias <ENV>=h1_dose --h1 <459dr>` are added before the build | bulk |

  - Nothing else may be new.
  - 630's UPRN, DNMAG, idall3, 642, 604 and K16/K20 wait for their own reads. If any is DECIDED before 524dr2 resolves, this
    prereg is amended first; the reader refuses any unwired switch (exit 2).
- **[A] Branch N0, no switch new:** do NOT submit.
  - Same recipe and an applied tree whose new knobs are all OFF should give the same content, and identical content must
    never take a slot.
  - Run `bulkid658.py <524dr2>.h5ad <524dr restored>.h5ad --exact` (restore with `scripts/restore_prediction_h5ad.py`).
  - IDENTICAL → the applied tree carries the October path bit for bit. Recorded; the queue moves on.
  - DIFFERENT → a defect found before Oct 22: diagnose locally.
- **Board:** `-m` from the MODEL-NAMES queue (`--claim 524dr2`), `-d` empty.

## 2. Prediction (fixed rule; numbers from captured entries)
**Σ rule (524, unchanged):** Σ = the sum, over the ON switches in the resolved recipe, of each switch's own A/B/C board read
(the deciding entry minus its reference). Per member, raw:

| switch | read | overall | pds | mse | nmae | fid | reach | jac |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| `lodo` | 446bj2 − 446b | +0.0022338 | +0.000208 | −0.014473 | +0.000060 | −0.000022 | −0.001963 | +0.000016 |
| `lodo_hi_pct` | 495 − 446bj2 | +0.0016340 | −0.001289 | −0.011991 | +0.000615 | +0.000058 | +0.000699 | +0.000021 |
| `sharpen_pt` | 479p2g − 446bj2 | +0.0010364 | +0.001587 | −0.003739 | +0.000120 | +0.000004 | −0.000836 | −0.000011 |
| `de_mag` | 446bj2z − 446bj2 | +0.0005493 | 0 | 0 | −0.002179 | 0 | 0 | 0 |
| **Σ (524dr's recipe)** | | **+0.0054536** | +0.000505 | −0.030202 | −0.001384 | +0.000039 | −0.002101 | +0.000026 |
| `mde_sign` [A] | 496a − 446bj2 | read on 10-05 (prereg: nmae −0.0015 … −0.0035 raw) | 0 | 0 | | pinned | unpriced | pinned |
| `tercilec` [A] | none on the board | **0** (4 A/B/C targets re-tier) | | | | | | |
| `h1_dose` [A] | 459dr − 453dr | only if R10 fired | | | | | | |

- **Secondary prediction (vs 453dr):** 524dr2 = 0.2663613 + Σ. With 524dr's switch set that is **0.2718149**; add the
  [A] rows that are ON.
- **[A] Primary prediction (vs 524dr):** 524dr2 − 524dr = the sum of the NEW switches' rows only. The rest of the path is the
  same recipe on a tree whose KAs proved its new knobs inert.
- The bulk-switch reads' DE parts are re-rolls (sd 0.00025 each; 537 §2.5), so Σ's DE part carries √3 × 0.00025 ≈ 0.00043
  of noise. 524's DE tolerances (below) are far wider.

## 3. Bins
**Secondary (vs 453dr): 524's bins, verbatim.**
- |Δ − Σ| ≤ 0.0020: **PASS**.
- Δ < Σ − 0.0020: **LAG**.
- Δ > Σ + 0.0020: **SUPER-ADDITIVE** (recorded, no action).
- **Falsifiers (diagnose before reading Δ):**
  - a DE member off its Σ by more than nmae 0.002, fid 0.003, reach 0.012 or jac 0.001 raw;
  - pds or mse moving while every ON switch is DE-only.

**[A] Primary (vs 524dr), by branch.** Only when 524dr itself read PASS; otherwise see §5.

| branch | bulk channel (pds, mse) | DE channel (nmae, fid, reach, jac) |
|---|---|---|
| **N1:** only `mde_sign` new | **pinned:** \|Δpds\|, \|Δmse\| ≤ 1e-6 raw, else FALSIFIER (tree / determinism defect) | **CARRIED** if \|ΔDE − hook's ΔDE\| ≤ 0.00025 (draws shared). **LAG** if ΔDE ≤ 0.25 × hook's ΔDE (the day does not deliver the hook). **PARTIAL** in between. Hook leak: \|Δfid − 496a's\|, \|Δjac − 496a's\| ≤ 0.001, \|Δreach − 496a's\| ≤ 0.01 raw, else FALSIFIER |
| **N2:** only bulk switches new | tercilec alone: **CARRIED** if \|Δbulk\| ≤ 0.0003 overall, else FALSIFIER (more than 4 targets moved, or the tree changed the re-level). H1 dose: **CONFIRMED** if Δbulk ≥ ½ × (459dr − 453dr)'s bulk, **NOT CONFIRMED** (October keeps ×2) if Δbulk < 0, **WEAK** (keep ×2) in between | a bulk step only re-rolls the DE: \|ΔDE\| ≤ 0.0005 (2 sd), else FALSIFIER |
| **N3:** both | as N2 | as N1, with tolerance 0.0005 (the bulk change re-rolls the draws) |
| **N0:** none new | not submitted (§1) | |

## 4. Who moves what: the culprit map for a LAG or a falsifier
| member off its prediction | owners among the ON switches | what to check first |
|---|---|---|
| **mse** short | `lodo`, `lodo_hi_pct` (and `h1_dose`) | `bulk_load` / `bulk_load_sha256` in the meta vs the hi6 dump; def_day's "LODO key ON" line; the hi-pct export |
| **pds** short | `sharpen_pt` / `tercilec` (and `lodo_hi_pct`'s −0.0013) | meta `bulk_sharpen_pt` (mode, edges, tier counts: 4 targets differ between tercile and tercilec on A/B/C); def_day's "per-target sharpening" line |
| pds and mse both off | the re-level | def_day's CONVERGED line, ρ per context vs 0.9793, `bulk_ctx_scale` vs 524dr's |
| **nmae** short | `de_mag` (zs), `mde_sign` (hook) | meta `de_mag` / `de_mde_sign`; the DEF_DAY_DE_MAG / DEF_DAY_MDE_SIGN pass-through in def_day.log |
| fid / jac / reach beyond tolerance | **no ON switch moves calls** | `rc_callset_*` keys in the meta, the emitter sha, the zero-repair record: a call-set change is a defect |

`read658.py` prints this map, member by member, whenever a bin is not PASS / CARRIED or a falsifier fires.

## 5. If 524dr (10-06) was not PASS
- On LAG or a falsifier, R6 makes 10-07 s1 a diagnosis arm or holds it.
- 524dr2 then also confirms the fix. Its primary reference is **453dr**, with the secondary bins and the full Σ, and the defect
  named in 524dr's diagnosis must be gone: its member back inside its Σ tolerance.
- The primary read vs 524dr is skipped (`read658.py dr` without `--prev`).

## 6. Class, tries
- **Class:** scorer mechanics / integration. It tests that the day script carries every decided switch on the applied tree,
  with no model change. Its [A] rows add the first October-path reads of the hook and tercilec, and R10's H1-dose
  confirmation if it fires.
- **Tries:** 1 arm. No knob is tuned on its read.
- A switch whose October-path contribution reads ≤ 0 comes out of the recipe, unless the defect is found and fixed with a
  rehearsal (524's rule).

## 7. Read commands
```bash
# 524dr (10-06), 524's own bins (653 R6):
uv run python experiments/658-prereg-1007-builds/read658.py dr <524dr id> \
  --recipe experiments/524-october-recipe-dr/work/recipe524.json [--hook <496a id>]
# 524dr2:
uv run python experiments/658-prereg-1007-builds/read658.py dr <524dr2 id> --recipe <w>/recipe524.json \
  --prev <524dr id> --prev-recipe experiments/524-october-recipe-dr/work/recipe524.json [--hook <496a id>] \
  [--alias <ENV>=h1_dose --h1 <459dr id>] --out experiments/658-prereg-1007-builds/work/read_524dr2.json
```
Exit codes: 0 = read, 2 = members hidden or a refusal, 3 = a falsifier.
