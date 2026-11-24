# 707: the strong-tier direction swap, an A/B/C board read (preregistration)

Written 2026-10-15 ~21:35Z by the science session [d80a10] (was [bf003b]), on the coordinator's YES (10-15 ~05:40Z, [400b6e];
order confirmed by [2d6f46] 21:3xZ). At this time the 707 dump exists (below), but no build, package or entry does.

## The arm (one variable vs ALLSW, entry `4WNmsNfdQhYx6bfib0dR`, 0.2760851)
- **707sw** = ALLSW's exact env, with `VCC_BULK_LOAD` = `work/dump_h1x2_cd4x15_lodo05hi6_sw025.npz`
  (sha256 77e9de2de8d04947de39e9eb4c43fd2bccd2f4029666ac60dd9bc4cae43d4f45). Re-levelled r1 at ALLSW's factors →
  ntclevel417 → r2 (rho within 0.5 % of 0.9793; at most one r3), as 702 / 705.
- **Construction** (`swap707.py`, the same code the day-of hook calls):
  - The 100 targets in the top tertile of z_src (522's production-donor strength rank; on A/B/C it equals 518's
    abc518.zprod exactly, KA-z) get their bulk row's non-panel part blended toward the production donors' equal-weight
    unit mean U = Σ_j r_j/|r_j|. The donors are K562-GWPS, Orion HCT116 and Orion HEK293T; the rows are raw delta rows.
    The blend is norm-preserving, at β 0.25.
  - Panel columns, the other 200 rows, the carriers (the DE design), perts, genes and env are unchanged.
  - Moved rows keep cos 0.972 … 0.991 to their input (median 0.982); the change is 2.4 % of the bulk's energy.
- **Known answers (PASS):** KA-β0 reproduces ALLSW's dump file byte for byte (sha d58750b8); KA-z PASS.
- **Expected config diff vs ALLSW:** r1 `bulk_load`, `bulk_load_sha256`; r2 adds `bulk_ctx_scale` (+
  `de_sign_emitter_sha256` if the emitter pin moved).
  - No `bulk_sharpen_pt`: every row keeps its norm to 4.5e-15. The P2g tercile edges are bitwise equal to ALLSW's (checked
    with the builder's arithmetic), and 0/300 targets change tier.
- **Where the change sits:** by A/B/C control CPM it follows the base's energy ([0,5) 0.153 vs 0.170; [5,20) 0.145 vs 0.141;
  [20,100) 0.448 vs 0.407; ≥ 100 0.254 vs 0.282). This is unlike 702, which piled onto 5–20 CPM. None of it falls outside
  the donor-gene union, and 30.6 % falls outside K562's gene axis.

## Class and the held-out evidence (706; AGENTS 11.1c)
**Transfer** (a better bulk direction for the strong targets from the same public donors, re-weighted; nothing fitted to
A/B/C). Held-out reads from experiments/706-strong-tier-budget (522's exact fold surfaces; 522 board units, k 1):

| surface | point | CI-low | folds > 0 |
|---|---|---|---|
| efp (3 October-like folds) | +0.0032 | −0.0001 | 3/3 |
| efr (3) | +0.0026 | +0.0011 | 3/3 |
| real (3) | +0.0043 | +0.0013 | 3/3 |
| **A/B/C cancer folds (3)** | **+0.0010** | −0.0002 | 3/3 |
| **A/B/C kolf2 fold** | **−0.0019** | — | — |
| h1-test (Arc lab; ew-cos) | top +8.5 % [−0.3, +22.6] | | mini-pds −0.009 |

- **Tries (11.1c.4): about 50 configurations** (3 selections × 3 β × centred / essential-free × 13 folds, plus h1-test).
  This one (top tertile, β 0.25) was picked as the best minimax form.
- **Best-of-N correction:** the configurations are strongly correlated (same folds, nested selections), so I take
  effective N ≈ 8–10, i.e. about 1.5 sd of each panel's price.
  - Corrected: real +0.0021, efr +0.0014, efp +0.0007, A/B/C cancer +0.0001.
  - With kolf2, A/B/C reads about −0.0004.
- **So the corrected A/B/C prediction is ≈ 0** (−0.0004 … +0.0001). The board read is the decisive, unbiased test. It is
  one fresh draw, so its bins stand uncorrected.
- 438 §4.4 (`abcdef` ≈ `abc` within 0.00014) shows that the October-like gain comes from the panel, not the day-of
  rules. D/E/F's panel is Arc-selected like A/B/C's, so the A/B/C board read is the best available proxy for D/E/F's value.

## Primary statistic and bins (board, vs ALLSW's own entry)
- **Primary: Δ(pds + mse) in overall units** = (scaled Δpds + scaled Δmse)/6. These members read exactly (rule 13).
- **TRANSFERS, ≥ +0.0010:** adopt for October, but only if the day-of hook (below) is applied with its known answers and
  one rehearsal by 10-20. Without both, it ships off.
- **WEAK, +0.0005 … +0.0010:** not adopted (below the bar the coordinator set for a late, code-bearing switch); recorded.
- **NULL, −0.0005 … +0.0005:** not adopted. The October-like gain is then a panel artefact, as 706's A/B/C folds suggest.
- **HARM, ≤ −0.0005:** not adopted. Record what the October-like folds and h1-test over-read (702's rules).
- **Falsifiers (STOP and diagnose before reading):**
  - DE members move by more than 0.0015 overall in total (a bulk-only arm moves them only by re-roll, σ ≈ 0.00025);
  - a config diff beyond the expected keys;
  - a final rho outside 0.5 % of 0.9793.

## Day-of construction (for October, only after TRANSFERS)
`dumphook707.sh` is a committed DEF_DAY_DUMP_HOOK script (690g's interface:
`bash <hook> <bundle> <work> <dump.npz> <identity.json> <work>/dumphook`). It calls `swap707.py` on the day's stack dump
with the day's bundle (z_src and U are recomputed for D/E/F's panel from the donor banks on disk). It is given to [17adab]
as an UNAPPLIED switch: `DEF_DAY_DUMP_HOOK` stays unset unless 707 reads TRANSFERS. The day-of z uses 522's
def_strength_key522 zsrc, the same function the rider (705) uses on the day.

## What it does not test
The October-like panels' larger gain (+0.003 … +0.004) cannot be read on A/B/C. This entry tests the Arc-panel proxy
only.

## AMENDMENT 1 (2026-10-15T21:35Z, [d80a10]; written before any 707 build, package or entry): the D/E/F-panel route
The coordinator [2d6f46] approved 708 (post-release re-pricing of panel-dependent switches on the REAL D/E/F panel, with
held-out public lines as truth). It adds a second October route for this switch:
- **TRANSFERS (board ≥ +0.0010):** unchanged. Adopt, with dumphook707 applied and one rehearsal.
- **WEAK or NULL (board −0.0005 … +0.0010):** the switch is decided by 708's read on D/E/F's own panel after release.
  It ships only if top-tier β 0.25 passes 438's rule verbatim on the D/E/F panel's PLOLO folds: a mean gain ≥ +0.0010
  (522 board units, k 1, ret 1), positive on 3/3 cancer folds, and a positive mean after leave-top-3-out. It then ships
  as a later final entry with its own build and gates (never the R0 entry).
- **HARM (board ≤ −0.0005):** off for October; 708 does not reopen it.
- Why this route is legitimate (AGENTS 11.1c): D/E/F's panel is Arc-selected, and 438 §4.4 showed these prices follow the
  panel, not the day-of rules. 708 measures that panel directly, on public lines the model never saw for those targets,
  and fits nothing to D/E/F's truth (there is none). The A/B/C read still vetoes a harm on Arc's own lab and truth.
- Tries for the D/E/F-panel read: one configuration (top tertile, β 0.25, the form frozen here). No re-selection on D/E/F.

## AMENDMENT 2 (2026-10-15T21:58Z, before any 707 build or entry): October adoption follows 708 AMENDMENT 2
For October this switch is adopted iff 708 AMENDMENT 2's relaxed rule holds on D/E/F's panel: the A/B/C read is not HARM,
and the D/E/F folds give a mean ≥ +0.0005, ≥ 2/3 folds positive and a positive leave-top-3-out. That replaces this file's
TRANSFERS route and AMENDMENT 1's 438-rule route. The A/B/C board bins above stay as the A/B/C read's labels.
