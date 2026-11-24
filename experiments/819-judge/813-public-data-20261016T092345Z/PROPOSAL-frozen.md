# 813 PROPOSAL-2 — TeloHAEC genome-wide CRISPRi as a fourth held-out BULK fold line (resubmission; frozen once written)

Written 2026-10-16 09:21–09:35Z by the 813 worker. It answers `VERDICT-1.md` (R1–R7, C1–C4) and the lab lead's calls of
09:16Z. It replaces PROPOSAL.md. Nobody edits it after it is written. Nothing has been downloaded or priced; spend $0.
Artifacts it freezes:

| artifact | sha256 |
|---|---|
| `kasg813.py` (KA-SG statistic; selftest PASS, 0.06 s) | `b94bf49ed92f89904a465418cf01370d81ece6edf835bc4aff322401435dddbe` |
| `patches/lib801-telohaec.patch` (UNAPPLIED; `patch --dry-run -p0` clean against lib801.py `41035aeb…`) | `67621126fc3910975c46c3c5ce9a5b492e53d594e73f5c8b7b0e24e612ea0d29` |
| `patches/ka_patch813.py` | `b53d34a35d55cdba901c457fc36c5d560a3f4d5d558dea375017cbaf515057c5` |

## 1. What and why
**What:** IGVF TeloHAEC genome-wide CRISPRi Perturb-seq, released 2026-10-12 (inventory: `RESULTS.md` §2).
- 19,276 genes, 114,885 guides.
- The A/B/C panel is 300/300 in the library, and L1–L6 are 298–300/300.
- In resting rep 1 alone, 295–299 targets per panel have ≥ 100 cell-incidences.

**Use:** a held-out truth line only. That gives L1–L6/telohaec_gw (6 folds) and abc/telohaec_gw (1 fold). It adds no donor
to any stack.

**Use (b) is deleted (R6).** The genome-wide rows as an 811 bank are 811's own judged, pre-counted try
(811 PROPOSAL-2 l.24–26, l.157). This lane only writes the bank that the folds need as truth.

**Mechanism (falsifiable):**
- Today every held-out read rests on three cancer lines, two of them from one lab (Orion). That pair is 811's strongest
  mutual partner (+0.011…+0.015, 811 RESULTS l.291).
- A line far from every donor tests whether an arm's gain is about the panel's targets or about lines like the donors.
- **105 is the record's warning:** the TeloHAEC pilot rows, used as a source for HCT116/HEK293T, read pds −0.141 / −0.127
  below K562 (`experiments/105-telohaec-residual/RESULTS.md` l.14–17).
  - Part of that deficit is pilot noise (67 cells per target).
  - The rest is distance. A fold that far away may compress arm differences toward zero, or reorder them.
  - **P2 (§5) is the test of whether such a fold still orders arms like the board.** If it does not, the line is
    rejected as an instrument.

**Class:** transfer (an instrument for pricing transfer). Nothing in it is fitted to A/B/C.

## 2. Truth rows from multi-guide cells (R1)
Rep 1: mean 4.70 guides per cell, median 4, 10 % with 0.

**Estimator A (primary), stratum-matched background.**
- Strata s = (n_guides_total, pool, rep). Carriers C_t are the cells that carry any guide of target t.
- P_bg(t) = Σ_s (|C_t∩s| / |s∖C_t|) · S(s∖C_t), with S the summed counts.
- δ̂_t = log1p(5e4·P_t/ΣP_t) − log1p(5e4·P_bg/ΣP_bg): the donor-table and pds space.
- The other guides and guide-load technicals cancel to first order. What is left is δ_t − μ̄, the mean effect of one
  random guide, the same for every target.

**Estimator B (check), sparse deconvolution.** Per-stratum intercepts plus one coefficient per target, solved as
XᵀX β = Xᵀy by conjugate gradient. XᵀX is 19 k × 19 k with ~3.5e7 nonzeros.

**KA-SIM, stage 0a: before any real expression row is read; it reads only obs/obsm metadata.**
- Inputs:
  - rep 1's real guide assignment and real `total_counts`: about 60 MB of range reads, already shown to work in
    `work/obsq813.py` / `cov813.py`;
  - base profile: the PILOT's control profile (`data/external/igvf-telohaec/b_ctrl.npy`, sha `cbc0f2cc…`, 70 KB). It is a
    TeloHAEC profile that is not the new data;
  - truth: K562-GWPS donor rows for the ~9.8 k targets it carries (S3 mirror `input/lab/repo/data/donors/k562_gwps/`),
    zero for the rest.
- Group sums are drawn Poisson from the summed per-cell rates L_c · softmax(log p + Σ_{t∈c} δ_t), on a 3,000-gene axis.
- The design is simulated three ways:
  1. the real multi-guide design;
  2. its single-guide twin (one random guide kept per cell);
  3. the real design with an injected guide-load offset: a fixed random direction of norm 0.5 × median‖δ_t‖ per guide.
- R = 200 replicates. The A/B choice and the KA-SG bounds are fixed here. Pass rules are unchanged from PROPOSAL:
  - A's median cos(δ̂, δ) ≥ 0.9 × the twin's;
  - ‖mean_t(δ̂_t − δ_t)‖ ≤ 0.10 × median‖δ_t‖ on designs 1 and 3;
  - noise-free recovery cos > 0.9999;
  - B replaces A only if B's median cosine is ≥ +0.02 higher.

**KA-SG, stage 0b, on real data (`kasg813.py`, frozen above).**
- Cell sets: A is computed on carriers and background with n_guides ≥ 2. S comes from cells with n_guides == 1. The two
  are **disjoint**, so their noises are independent.
- Two noise-unbiased statistics over the 100 stage-0 targets:
  - **β** = Σ⟨A,S⟩ / Σ⟨S^h1,S^h2⟩: A's signal slope on S's signal, not attenuated by noise. E ≈ 1 if A is unbiased.
  - **dsh** = ⟨m^h1,m^h2⟩ / mean⟨A^h1,A^h2⟩ with m^h = mean_t(A^h − S^h): the energy of a shared A−S offset, noise-free in
    expectation. E ≈ 0.
- **The bounds are set in stage 0a, before any real row:** [lo, hi] = the central 95 % of each statistic over the R
  unbiased replicates.
- **Power row:** the naive zero-guide estimator on design 3 must fall outside [lo, hi] in ≥ 80 % of replicates. If it
  does not, KA-SG is recorded as **uninformative at this depth** and is not a gate. KA-SIM stays the estimator KA.
- **Blind spot, stated:** dsh cannot see μ̄, because A (k-matched background) and S (whose k = 1 background carries one
  other guide) both carry −μ̄. It tests only k-dependent technical offsets. μ̄ is G0-BLOB's to catch.
- **Freeze (C3):**
  - Stage 0a writes `freeze813.json` (estimator, β and dsh bounds, power verdict, R, seeds) to S3. Its sha256 is
    recorded in RESULTS.md and committed by the lab lead before stage 0b is launched.
  - Stage 0b refuses to start if the json's sha differs.
  - No re-tune after a failed gate.

## 3. Depth, the 5-CPM gate, the blob test
- **Depth** (median 5,452 UMI): because each cell counts toward ~4.7 targets, per-target pseudobulk depth is ~3.0 M UMIs
  over both resting reps (EST). K562-GWPS is ~2.0 M and Orion HCT116 ~2.3 M. The other guides add variance, not bias (§2).
  G0-POWER measures what is left.
- **5-CPM gate:**
  - The bank's gene axis is the official-axis genes whose **mean per-cell CPM over TeloHAEC's background cells is > 5**:
    the scorer's `filter_gene_min_cpm_cell` (`docs/reference/vcc2026-metrics-BRIEF.md` l.83) on the line's own
    control-equivalent.
  - The fold's scored set is that axis minus panel columns. lib438's `scored` takes the held bank's on-axis genes, so the
    gate lives in the bank and no engine change is needed.
  - Stage 0b reports the gene count and the split-half reliability per CPM bin.
- **G0-BLOB (615's statistic, unchanged):** split-half share on the 100 targets with A's control, bar ≤ 0.10. The
  zero-guide and NTC-only controls are reported only.

## 4. Stages, gates and costs (C1: each stage runs alone, and the next starts only after the lab lead reads its RESULTS)

| stage | what | instance, cap | cost |
|---|---|---|---|
| **T (terms, R7)** | **DONE 09:18Z, PASS.** IGVF External Resource Sharing Policy: "External users may freely download, analyze and publish results using any … IGVF resources without restrictions"; citation requested. Text sha `fc17237b…` and the reading are in `RESULTS.md` §6. Re-checked as the first step of stage 0b: if the sentence has gone, STOP with no byte fetched. | — | $0 |
| **0a** | KA-SIM, bounds, A/B choice; metadata range reads only | r7i.xlarge, 45 min | ≤ $0.5 |
| **0b** | **terms re-check**; fetch resting rep 1 (IGVFFI2522NCFR, 66.3 GB) EC2-side; one CSR pass; G0-KD, G0-POWER, G0-BLOB, KA-SG on the 100 seeded A/B/C panel targets (seed 813, sorted list, ≥ 100 incidences, ≤ 250 carrier cells each) | r7i.2xlarge, 200 GB gp3, 120 min | ≤ $1.5 (≤ $3 with one retry) |
| 1a | bank: estimator A over all 19,276 targets on resting reps 1+2 (+44.5 GB); `build_donor_table` format to `s3://…/output/lab/813/bank/telohaec_gw/` | r7i.2xlarge, 200 min | ≤ $2 |
| 1b | **R3 fold build**, below | r7i.8xlarge, 75 min | ≤ $3 |

**Lane cap: stage 0 ≤ $5, stage 1 ≤ $6, total ≤ $11.** Per launch ≤ $50; `labcloud.sh` only; Project=lab; on-demand (C4).
Nothing > 200 MB to the Mac, and no Mac transfer at all in a WINDOW. The stimulated conditions are dropped with use (b).

**Stage-0 gates (all must pass; a failure is a recorded negative with its mechanism):**
1. **G0-KD:** own-gene residual median ≤ 0.60, with ≥ 60 % of targets < 0.70 (targets with own CPM ≥ 10).
2. **G0-POWER:** median split-half residual cosine (271's definition) ≥ +0.028, with share ≤ 0 at most 34 %.
3. **G0-BLOB:** split-half share ≤ 0.10.
4. **KA-SIM passes,** and **KA-SG** is inside its frozen bounds (or recorded as uninformative per §2).

**R3, the fold build:**
- **L1–L6/telohaec_gw:** one new full 3-donor operator per panel. No full-donor operator exists for any L panel; 801
  built held-out fits only.
  - Day recipe: `scripts/build_union_whitener.py --bundle <L bundle> --donors k562_gwps,orion_hct116,orion_hek293t
    --basis-prune essential --max-rss-gib 14`, with the 12k cap as in 801's runs.
  - Cost basis: 801's day-recipe ABCp fits (2 donors) took **322 s** each at 8 threads, peak **12.0 GiB**
    (`s3://…/output/lab/20261016T080732Z-lab801-abcp-a/jobs/2.log`; run cost $0.70 for 3 fits plus dumps on r7i.8xlarge).
    A 3-donor fit has a ~1.5× basis on the same 17.8 k² covariance. Estimate ≤ 10 min and ≤ 18 GiB each, so six run in
    parallel at 5 threads within the 75-min cap.
  - **sha recording:** each operator is stored at `s3://…/output/lab/ops801/L<k>/full_telohaec_gw.npz` with its sha256
    and axis sha in the new registry entry, under the same fields as existing entries.
  - **Reproducibility rule:** fits are not bit-reproducible across machines (801 RESULTS l.30). A fold is only ever
    rebuilt from the stored operator, pulled by sha, and never refit. If a refit is ever forced, it is a new fold id.
- **abc/telohaec_gw:** production's operator `data/donors/whitener_union_pruned.npz`, sha256
  `e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123`, as ABC/kolf2 uses (registry). It is the day recipe
  (LAB-LOG 08:05Z).
  - Built through lib801.bind("abc") and dump801 as a registry DumpFold.
  - Route KA: the same route re-dumps abc/kolf2 and must match ABC/kolf2's stored rows to ≤ 1e-12. 801 measured this
    route against 518 at 8.7e-15.
- **The lib801 change is an UNAPPLIED patch** (`patches/lib801-telohaec.patch`). It adds only new keys: `FOLD_SOURCES
  ["telohaec_gw"] = MAIN`, a runtime `plolo.DONOR_PATHS.setdefault` (plolo.py untouched), and
  `PANELS["abc"]["fold_ops"].setdefault("telohaec_gw", PROD_OP)`. `PANELS["abc"]["folds"]` is unchanged, so no loop over
  existing folds changes.
  - **KA (`patches/ka_patch813.py`):**
    - `snapshot` before and after applying it; `check` PASSes only if every pre-existing registry entry, dump sha,
      FOLD_SOURCES key and fold_ops key is byte-identical;
    - `rows` re-dumps L1/k562_gwps and abc/kolf2 through the patched lib and requires bit-identical rows (same instance type
      and threads as the stored dump), else ≤ 1e-12 with max |d| reported.
  - 801 or the lab lead applies it after a GO. 813 never edits 801's files.
- Then the engine's own fold KAs, and the P2 rows (§5).

## 5. Predictions and bins (R2, R4)
- **P1, stage 0** (expected values): G0-KD median 0.35–0.60; G0-POWER +0.03…+0.08; G0-BLOB 0.03–0.10; KA-SIM A ≥ 0.9× twin.
  P(all pass) ≈ 0.6. Likeliest failure: G0-BLOB (guide-load technicals that stratum matching does not absorb).
- **P2, retrodiction:**
  - **702 is dropped on TeloHAEC folds.** Its fold arm blends the HELD line's own PIE rows, and no TeloHAEC PIE rows
    exist.
    - The line is unseen: PIE's training and evidence sets are tahoe, jiang, arc_vcc_25 and orion, with Replogle at
      weight 0 (`experiments/700-pie-transfer/RESULTS.md` l.49).
    - But building a PIE query for a new context (a preprocessed control set plus context text) is an untested
      construction with choices of its own. And 702 reads only NULL (1/6) even on fresh K562 folds (LAB-LOG known-answer
      table), so it would discriminate little here.
  - **Board-pointed negatives in DESIGN §9 that the engine's named arms can express on a TeloHAEC fold: only 479sp2**
    (`cap3+S_P2` vs `cap3`). Every lab surface misses its class: bulk 18/18 and 3/3, and DE-H1 (LAB-LOG 03:20Z, 07:55Z). It
    is **reported, not graded**. The others (458d, K20/K16, idall3, DN) are DE-channel rows that bulk folds cannot express.
  - **So the graded rows are positives only:** 495 (`cap3` → `cap6`), 479p2g (`cap3` → `cap3+P2g`), 498p2g (`cap6` →
    `cap6+P2g`).
  - **Bins, renamed (no "TRUSTED"):**
    - **SIGN-CONSISTENT** = all three positive, with a lower bound > 0 on abc/telohaec_gw, and each positive on ≥ 4/6 L folds.
      This is a necessary condition only.
    - **INCONSISTENT** = any of the three with an upper bound < 0 on abc/telohaec_gw, or ≤ 3/6 L folds positive on ≥ 2 rows.
      → **REJECTED as an instrument.**
  - **The folds carry the label "uncalibrated for negatives" until a discriminating board row is graded on them.** Rule,
    fixed now: the first bulk-expressible rows that get a board read after the build (707's strong-tier swap; 710's
    compositions) are graded on the TeloHAEC folds exactly as on the cancer folds. ρ per class is reported beside the
    cancer folds'.
- **P3 is dropped as a blob sentinel (R4).** A TeloHAEC fold has three donors and a cancer fold has two, so the level
  comparison is confounded. G0-BLOB is the blob gate. Levels are reported descriptively only.

## 6. Decision role and leakage (R5, preregistered)
- **Role:**
  - TeloHAEC folds are **reported apart** (like kolf2). They are never pooled into BULK-W and never a veto on their own.
  - A candidate that gains on the cancer folds but whose TeloHAEC CI upper bound is < 0 carries a written flag in its
    dossier, with a mechanism. The coordinator weighs it.
  - Reason: all six challenge lines are cancer lines (AGENTS.md l.23), and TeloHAEC is not.
- **Leakage:**
  - On a TeloHAEC-held fold, every candidate's TeloHAEC rows (the pilot `telohaec_pilot` and `telohaec_gw`) are dropped,
    as PLOLO drops the held line.
  - **Refusal to implement in check800 (the lab lead implements; 813 does not edit it):**
    - candidate.json gains a required `banks` list: every external bank or line the candidate reads, per DESIGN §7's
      "candidates declare their training lines".
    - On any fold whose held line is `telohaec_gw`:
      - if `banks` is missing, the price is **refused** (fail-closed);
      - if it names `telohaec_pilot`, `telohaec_gw` or any bank whose line is TeloHAEC, the fold is marked
        **declared-seen**: its read may bound a cost and never supplies a gain (702 rule 2), and the report prints
        "seen: cost-bound only".
    - The wear entry records `surface=telohaec`.
- **C2:** no TeloHAEC fold, dump or score on an LZ panel (L7/L8). Bank rows for LZ targets may exist on S3 but are not
  read.

## 7. Novelty (R6)
Not these closures:
- **The pilot as a stack donor** (223/225 `telohaec@0.28` −0.00120; FINDINGS l.885 "sources saturate").
- **122/135's atlas** (135 terminal NO-GO).
- **105's direct source** (pds −0.13…−0.14 vs K562). Here TeloHAEC is never a source; 105's distance is what P2 must
  survive.

New: no RESULTS on disk uses TeloHAEC as a fold truth (VERDICT-1 check 1), and the genome-wide release is four days old.

## 8. Tries and day-of
- **Tries:**
  - Estimator: 2, chosen on simulation only.
  - KA-SG statistic and bounds: frozen before real rows.
  - Control, gate, targets: 1 each.
  - Priced arms on BULK-W from this lane: 0. P2's rows are the engine's named known-answer arms; they choose nothing.
- **Day-of (Oct 22 → Nov 5):** the bank is public and fixed and covers ~all genes. A telohaec_gw-held fold on D/E/F's real
  panel is built the same day (a full 3-donor fit on that panel, as 1b), next to the cancer folds in 708's re-pricing. It is
  reported apart (§6), with nothing fitted to A/B/C.

## 9. Conditions carried into a GO (verbatim from VERDICT-1)
- **C1:** stage 0 runs alone. Stage 1 is launched only after every stage-0 gate passes and the lab lead has read the
  stage-0 RESULTS.
- **C2:** no TeloHAEC fold, dump or score on an LZ panel. Bank rows for LZ targets may exist on S3, but they are not read.
- **C3:** the estimator (A or B) and the KA-SG statistic are frozen by sha before the first real row is read. A failed gate
  is recorded as a negative with its mechanism, with no re-tune. The DLD-1 fallback (frozen §9) is a new proposal, not an
  automatic branch.
- **C4:** cloud via `labcloud.sh` only, on-demand, Project=lab, ≤ $50 per launch, the lane cap as resubmitted. Nothing
  > 200 MB to the Mac, and no Mac transfer at all in a WINDOW.
