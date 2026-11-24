# 615 — Trust gate for the Orion raw HEK293T real-emitter DE harness (fixed BEFORE any harness run)

Written 2026-10-02 ~18:20Z by [21dd4c], after the blob PASS (395c186f). The truth, reference and keys were still being
built, and no 615 emitter run existed. The coordinator's instruction: the SAME gate as 614, bands committed before any
run. AGENTS 11.1 rule 3: the harness prices a lever class only after it retrodicts that class's board pairs. If this
gate FAILS, every 615 number is descriptive only.

## The harness (stages in RESULTS.md)
- **Truth:** Orion raw HEK293T on the official 18,533-gene axis.
  - 100 seeded panel targets (`work/targets615.txt`), ≤ 250 cells each (median 220). These are the exact cells the
    blob test read (`harness615_rows.npz`).
  - 18,400 controls drawn with the preregistered P′ batch mixture (`truth615.py`).
  - Native depth ~15k. The full axis has its own sub-gate genes, so no reservoir is needed (614's A1 does not arise).
- **Reference:** `614/ref614.py`, evaluate323's `reference()` logic → `work/ref_hek`.
- **Scorer:** evaluate323's `score()`, i.e. cell-eval2's own DE member code (vcc2026 preset).
- **Emitter:** run496's driver, production emitter, via `run615.py`, a path-only copy of run614. It swaps:
  - the reference dir;
  - the key (`work/key615.npz`);
  - CLASSES3 = the reference's control-library quartiles.
  - All four CLEAN4 banks are kept (none is Orion-derived; asserted).
- **KA for the wrapper:** `run615.py h1ka` must reproduce 435's H1 stats (`ka615.py`, the 614 KA). If it fails, nothing
  below is read.
- **Design key = PROXY** (coordinator 18:2xZ: proxy first, the faithful key is the fallback).
  - k562_gwps delta rows × one scalar s fitted on H1 only (`key615.py`).
  - Both Orion lines are excluded. The Nadig lines and k562_essential cover 0 / 300 panel targets, so the proxy is a
    single-line key, weaker than production's 3-line whitened consensus.

## The gate: 614's pairs and bands, copied UNCHANGED (3 seeds each, f 1 = the HEK truth's own BH set)
| pair | harness construction | board raw Δ (611) | required on HEK |
|---|---|---|---|
| **G1 zs** (446bj2z − 446bj2) | the same run, pre-zs vs `--zs dn_k250_zs` (prodem, coarse MDE, zero repair −5) | nmae −0.00218 | Δnmae in **[−0.0065, −0.00073]** |
| **G2 coarse MDE** (433a − 417e) | prodem without MDE vs prodem + MDE 0.1 / 0.1 / 0.3 coarse (no zs, no zero repair, as both board arms) | nmae −0.01679 | Δnmae in **[−0.0504, −0.0056]** |
| **G3 call set** (446bj2zcp − 446bj2z) | 528's call-propensity block ("kept") vs the CPM block, on the zs base (zero repair −5) | fid −0.00358, jac +0.00399 | **Δfid < 0 AND Δjac > 0** |

- **PASS iff G1, G2 and G3 all hold** (mean of seeds 0–2). `gate615.py` is a copy of gate614.py.
- Reported with no rule attached:
  - each pair at f 0.5 / 0.35;
  - G2's and G3's reach;
  - per-seed values;
  - the count of gated cells with lfc < −30 per arm (7772dc's `654/zerocount654.py`).

## The low-power diagnostic, read BEFORE the gate runs (coordinator's step 2; `diag615.py` → `work/diag615.json`)
- **Reported:** median n_real, its q90, the targets with n_real ≥ 10 (the ones nmae scores), and the fraction of DE calls
  that are negative.
- **Comparators:** K562 (614: median 5, 34 / 100 targets ≥ 10) and H1 (h1-test: median 89, 58 / 85).
- **"K562-like" means** median n_real ≤ 10 OR fewer than 45 of 100 targets with n_real ≥ 10. It is said in RESULTS.md
  and to the coordinator before any gate run starts. The gate still runs, with unchanged bands.

## If the gate FAILS (fixed now, before any number)
- **Descriptive depth control first:** H1@15k. That is h1-test binomially thinned to the HEK truth's median control
  library (614's `thin614.py` with p = HEK / H1), with the same pairs. It runs before any fallback decision.
- **Faithful-key re-run (the coordinator's "key-attributable failure"):** when the power diagnostic is NOT K562-like AND
  H1@15k's G2 Δnmae is ≤ −0.0056 (depth does not shrink the lever).
  - The faithful key (run323's stage / operator / compose with K562-only donors, CD4 kept, the whitener refit without
    both Orion lines and pruned of the 100 targets) is then built.
  - This gate is re-run once, unchanged.
- **Otherwise the failure is attributed** to power or depth, with that diagnostic, and no re-run is made.

## After a PASS (not before): register, then price
- Each lever's prediction is committed before its HEK run (native depth, plus 633's d20 rescale beside it where the
  lever is add-based).
- The coordinator's queue: CPH, DN, K16/K20-class call sets, 630 v3 `b_k100_q20r3p10`, 631 / v7.
- Also the hook forms, plain CLEAN4 hook (c4h) vs idall3, on the zs base.
- Every lever row carries `--zero-repair -5`.

## STOPPED before any gate row (coordinator, 2026-10-02 19:06Z)
The coordinator's written reason, verbatim:
> "the preregistered pre-gate power diagnostic (diag615.json) reads worse than K562, whose gate failed on power (median 1
> DE gene per target, 18/100 nmae-scored). The gate's 2–3 h of a shared slot would measure power, not the harness.
> Stopped by the coordinator before any gate row, 10-02."

- chain615b was killed at 19:05Z, after its first step had been queued.
- The wrapper KA (`run615 h1ka` == 435, `ka615.py`) is kept and finishes on its own. It validates run615 for any later
  use of this harness.
- No G1 / G2 / G3 run exists. 615 is not a trusted instrument, and no 615 lever is priced. PRED-615-levers.md stays
  registered and unread.
