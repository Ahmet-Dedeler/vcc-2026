# 614 — A real-emitter DE harness on K562: built and KA'd end to end, but its pre-registered trust gate FAILS, so every 614 number is descriptive

2026-10-02, council session [21dd4c], the coordinator's assignment (spec: 631/SPEC-614-K562-harness.md, [d4cd5f]).
Gate fixed first: `TRUSTGATE.md` (cf5f202f), amendment A1 (e0a3ed54) and reading rule A1-R, all before any valid run.
Verdict script: `gate614.py` → `work/gate614.json`.

## What was built (each with a known answer)
| piece | KA |
|---|---|
| truth: K562 GWPS, 100 panel targets (seed 614), 100–250 cells, 18,400 NTC, 8,248 genes (`k562harness480` → work/k562-gwps-c250.h5) | integer counts; 223 eligible targets |
| A1 sub-gate reservoir (`sub614.py`): 7,301 synthetic sub-5-CPM genes, H1's share 0.58 % and profile (the emitter needs them; GWPS has none) | 20 of them land just over 5 CPM (null, recorded) |
| reference (`ref614.py` = evaluate323.reference minus H1's asserts) | **its H1 build is IDENTICAL to 323's reference.npz** (arrays, DE parquet sha) |
| proxy key (`key614.py`): Orion-only donor deltas × s = 1.124 fitted on H1 only | on H1 the donor mean tracks cap15 at r 0.73 |
| wrapper (`run614.py`): swaps the reference loader, key, CLEAN4 banks (K562 banks dropped) and library classes (H1's are its quartiles) | **h1ka mode == 435's stats, max \|Δ\| 0** |
| G3 call-propensity key (`key614cp.py`, 528's lcap) | **reproduces 528's H1 lcap bit for bit** |
| scorer: evaluate323.score() (cell-eval2's own DE member code, vcc2026 preset) | — |

## The gate (3 seeds, f 1; raw Δ)
| pair | K562 Δ (per seed) | board Δ | band / rule | result |
|---|---|---|---|---|
| G1 zs | nmae **−0.00409** (−0.00429 / −0.00395 / −0.00403) | −0.00218 | [−0.0065, −0.00073] | **PASS** (1.9× the board) |
| G2 coarse MDE | nmae **−0.00075** (−0.00197 / −0.00069 / +0.00043); reach −0.0034; fid +0.0016 | nmae −0.01679, reach +0.0024 | [−0.0504, −0.0056] | **FAIL** (≈ 1/22 of the board; not even resolved across seeds) |
| G3 call set (lcap-kept) | fid **+0.00009** (+0.00013 / +0.00040 / −0.00026); jac +0.00030 (all three seeds) | fid −0.00358, jac +0.00399 | fid < 0 AND jac > 0 | **FAIL** (jac's sign is right but ~1/13 the size; fid is noise) |
- **Verdict: FAIL**, so per TRUSTGATE every 614 number is descriptive only.
- The first gate614 "FAIL" (crashed runs, before A1) was VOID and is kept as `work/gate614.VOID-crashed-runs.json`.

## Why (diagnosis, not a rescue)
- **The truth has little power.** Median n_real is 5 DE genes per target (q90 64, max 386). nmae scores only the 34 of 100
  targets with ≥ 10 gated genes. jac's level is ~0.006 against the board's ~0.04.
- The coarse MDE add raises |lfc| of weakly-detected genes toward the minimum detectable effect, and those are exactly
  the genes a low-power truth does not gate. So the lever has almost nothing to act on. A call-set composition change has
  equally little room in a tiny DE set.
- zs works on the delivered magnitudes of strong calls, which K562's gate does see, so G1 passes.
- So the harness can price levers that act on STRONG genes, but not MDE-type or call-set levers. And the gate requires
  all three classes.
- K562 GWPS has at most ~250 cells per target and ~11k UMIs, so more power is not available from this source.

## Still running (descriptive rows only)
- **H1@11k control** (chain614d: h1-test thinned to K562's depth; G1 / G2 on it).
  - If H1@11k reproduces the board's MDE size, K562's failure is TRUTH POWER or LINE.
  - If H1@11k also shrinks it, DEPTH is a factor.
- **Hook forms** (chain614h: c4h, idall3 on G1's base), against `PRED-614-levers.md`. Descriptive only now.

## The pre-registered faithful-key fallback
- TRUSTGATE says: if the gate fails on the PROXY key, build the faithful key (run323.stage(): Orion-only donors, CD4,
  the whitener pruned of K562's targets) and re-run the same gate once.
- This diagnosis points at the truth, not the key: a better key cannot add gated genes. So the expected value of that
  re-run is low and its cost is hours.
- The decision is the coordinator's. It is recorded here either way.

## Forward question: is a higher-power 614 possible before Oct 22? (`inventory614.py`, measured on disk, 18:05Z)
Bar (coordinator): ≥ 100 cells per target AND ≥ 20k median UMI per cell, in a cancer line.
| dataset (on disk) | line | median cells / target | targets with ≥ 100 cells | median UMI | meets the bar |
|---|---|---|---|---|---|
| **X-Atlas / Orion, raw** (332 batch parquets) | **HEK293T** | **200** | **14,719** (of 18,310) | **20,699** | **yes** (transformed line, not tumour-derived) |
| X-Atlas / Orion, raw | HCT116 (colorectal carcinoma) | 150 | 12,711 | 17,292 (q75 23.8k) | close |
| Replogle K562 GWPS (614's truth) | K562 | 178 | 7,990 | 11,304 | no (depth) |
| Replogle K562 essential | K562 | 121 | 1,299 | 13,207 | no |
| Nadig 2025 | HepG2 | 45 | 256 | 15,641 | no (cells) |
| Nadig 2025 | Jurkat | 83 | 906 | 10,160 | no |
| Replogle RPE1 | RPE1 (not cancer) | 72 | 783 | 11,940 | no |
| Frangieh 2021 (Cas9 KO, Perturb-CITE) | melanoma | 659 | 230 | 10,988 | no (depth, modality) |
| Myllia A549 | A549 | (label column not found) | — | 13,578 | no (depth) |
- **Yes, via Orion's RAW data.** The repo's Orion harnesses were thinned to 12–14k; the raw cells are 17–21k with 150–200
  cells per target, about 1.6–1.8× K562's depth.
- Caveats that a higher-power 614 would have to design for:
  - (1) Orion's thinned truths carry a target-independent "technical blob" (26–44 % of effect energy; FINDINGS H.4).
    With 109 / 223 batches, batch-matched controls (controls drawn from each target's own batches) are the first thing
    to test against it.
  - (2) The design key must exclude both Orion lines, the truth's own and its same-platform sister, so donors are K562
    (+ the Nadig lines).
  - (3) CRISPRi is the same modality as Arc's. HEK293T is transformed, not tumour-derived, and HCT116 is the
    tumour-derived but slightly shallower option.
- The 614 code is generic, so the new work is an Orion extractor (parquet → the CSR-h5 truth, with batch-matched
  controls) plus the same pre-registered gate.
- **CORRECTION (18:10Z, measured in 615's extraction).** The Orion median-UMI figures above are `total_counts` over all
  38.6k Orion features.
  - On the scorer's 18,533-gene axis, raw HEK293T is **15.0k** for panel cells and 15.7k for controls (615
    `work/hek-raw.json`). HCT116 is **14.4k** (the reservoir build `data/harness/orion-hct116.h5ad.meta.json`).
  - So **neither meets the ≥ 20k bar.** HEK is ~0.75× A/B/C's 20.1k and ~1.33× K562's 11.3k (K562's truth library
    median, thin614).
  - The other rows' UMI come from each file's own count column, so compare them with Orion's all-feature figure, not
    with 15.0k.
  - HEK raw: 235 of the 300 panel targets have ≥ 100 cells (median 200 cells).

## Hook forms on K562 (descriptive: the gate FAILED; `hook614.py`, `work/hook614.json`; chain614h, crash fixed and relaunched 18:05Z)
Base = G1's zs construction (prodem, coarse MDE, zero repair −5). CLEAN4 here is Nadig Jurkat + HepG2 only (both K562
banks dropped): the weakened-CLEAN4 stress case.

| contrast (f 1, raw) | s0 / s1 / s2 nmae | mean nmae | mean reach | registered (PRED-614-levers.md) |
|---|---|---|---|---|
| c4h − no hook | +0.00021 / +0.00019 / +0.00040 | **+0.00027** | −0.0001 | > 0: **MATCH** (613d predicted +0.0014; the real size is ~1/5) |
| idall3 − c4h | −0.00018 / −0.00023 / −0.00038 | **−0.00026** | −0.0059 (+0.0009 / −0.0147 / −0.0037) | > 0: **MISS** (613d predicted +0.0005) |
| idall3 − no hook | | +0.00001 | −0.0060 | |

- The plain hook costs a little on K562 (3 / 3 seeds), and idall3 removes that cost (3 / 3 seeds), so idall3 ≈ no hook in
  nmae.
- 613d's K562 "stress tail" (idall3 worse than c4h) does NOT appear on the real emitter. 613d over-states the K562 tail,
  as the registration's miss branch says.
- The reach drop is one noisy seed on a truth with a median of 5 DE genes per target. Not read.
- Descriptive only (gate FAIL, low-power truth). For idall3's case: the K562 row that "opposes" idall3 came from 613d.
  On the real emitter it is neutral to slightly favourable in nmae.

## Faithful-key re-run (TRUSTGATE 3d851ae7 + amendment (8); coordinator: run, no waiver) — key built 10-03 17:03Z
- `fkey614.py`: 618's KA'd run323 path with the 8 scoped differences. Inputs were re-frozen after the donor-weights fix
  (e0eed20d).
  - stage 0 s, operator 5.4 min (quiet admission), compose 4.4 GB.
- **The key, described before any gate read:**
  - 100 / 100 targets have a row; 98 % of entries are nonzero (whitened + CD4 residual);
  - per-target Pearson with the proxy key (key614.npz: the two Orion donors' raw consensus × s) has median 0.83
    (q10 0.74).
  - So the re-run tests a whitened, CD4-adjusted version of mostly the same donor information. That fits the scope note's
    expectation that the key is not the binding constraint, but the 12 gate runs decide.
- The 12 gate runs (`run614 k562f / k562fcp`, the same flags as the first gate) started at 17:03Z. Verdict:
  `gate614f.py` → `work/gate614f.json`.

## Faithful-key re-run: **FAIL**, so 614 is CLOSED (17:39Z)
| pair (3 seeds, f 1, raw Δ) | proxy key (first gate) | **faithful key** | band / rule | board |
|---|---|---|---|---|
| G1 zs Δnmae | −0.00409 (PASS) | **−0.00299 (PASS)** | [−0.0065, −0.00073] | −0.00218 |
| G2 coarse MDE Δnmae | −0.00075 (FAIL) | **−0.00186 (FAIL)** (reach −0.0103, fid +0.0021) | [−0.0504, −0.0056] | −0.01679 |
| G3 cp Δfid / Δjac | +0.00009 / +0.0003 (FAIL) | **+0.00147 / +0.00034 (FAIL)** | fid < 0 AND jac > 0 | −0.00358 / +0.00399 |

- The faithful key (production's whitening and CD4, Orion-only donors) lifts G2 2.5× but stays far outside the band, and
  G3's fid moves the wrong way.
- So the K562 failure is not the key. It is truth power and / or line biology (median 5 DE genes per target), which this
  re-run cannot separate (scope note).
- **614 is closed as an instrument.** Its hook rows stay descriptive.
- Side note, unexplained and not pursued: the first (proxy-key) gate's G1 pair showed mse −0.047 between the pre-zs and zs
  arms, while the faithful-key runs show mse 0. zs is a DE-only pass, so the proxy-key value is a bulk artefact of the
  first runs' construction (sub-gate reservoir era). No reading depended on it.
