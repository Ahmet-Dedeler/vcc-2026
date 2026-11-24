# 590 — The day-of identity call is right whenever the line is in CCLE; a line ABSENT from CCLE gets a false unanimous call. A similarity floor (min top-1 correlation ≥ 0.65) removes every false switch-on

2026-10-01, session "Explore contribution areas [57577c]" (numbers 590–599). Analysis only: no build, no `vcc`, no
`cell-eval2`; every process ≤ 4.1 GB (`/usr/bin/time -l`). Preregistered before any query was scored
(`PREREGISTRATION.md`). Class (AGENTS §11.1c): **day-of procedure** (it fits nothing to A/B/C's board; A/B/C enter only
as the Flex reference values of the measure).

## Answer
- `scripts/def_identity.py` (041's matcher) identifies a line correctly **whenever that line is in the CCLE matrix**:
  98.5 % top-1 over 390 labelled profiles, 99.5 % when the 8 methods agree. Every wrong call on an in-CCLE line is a
  Cellosaurus-documented duplicate or same-patient sister: PC-14 = PC-9, KCI-MOH1 = HPAC, KPL-1 = MCF-7 ("Problematic cell
  line: Contaminated. Shown to be a … derivative"), GP2d/GP5d (same individual), A549 → A549H3.
- **146 "sister" profiles** (in CCLE, same lineage as a registry line, not themselves registered: colon, lung, breast,
  pancreas, liver, blood lines) produced **0 false switch-ons**.
- **The failure is a line absent from CCLE.** Both HEK293 profiles (Orion HEK293T, Xu HEK293) are called HAP1_ENGINEERED
  by all 8 methods. HAP1 is a registry line (Jiang's bank), so today's rule would switch the wrong rows on.
  Preregistered rule: FP_reg > 0, so a gate is proposed.
- **What separates them is absolute match quality.** The measure is the lowest top-1 correlation over the 8 methods
  (`min_top1_score`):
  - the false calls: 0.554 and 0.606 (Flex-simulated: 0.506 and 0.561);
  - correct unanimous calls: q01 0.681, q05 0.719, median 0.769;
  - A/B/C, measured exactly as D/E/F will be (Flex, 18,400 cells): **0.772 / 0.719 / 0.733**.
- **Proposed gate (UNAPPLIED, for the coordinator): same-line rows switch on only if the call is registered, unanimous
  AND `min_top1_score` ≥ 0.65.**
  - It removes both false switch-ons in the plain and the Flex-simulated read.
  - It costs 2 of 19 correct registry switch-ons, none of them D/E/F-shaped. Plain: Nadig HepG2's b_ctrl, which carries
    only 8,828 genes, 0.540; sci-RNA-seq3 A549, 0.636. Flex-simulated: the same HepG2 profile, 0.396; ECCITE THP-1, 0.636.
  - A/B/C clear it by ≥ 0.069.
  - Patch + 8 tests ready: `patch_identity590.py`, `test_identity590.py` (ALL PASS).

## Instrument and known answer
- `profiles590.py`: control pseudobulks per labelled line on the official 18,533-gene axis (NaN where a dataset lacks a
  gene), read in h5py chunks.
- `match590.py`: imports `def_identity`'s own functions and runs Oct 22's arithmetic on every profile. That means the same
  CCLE cache, the official ∩ CCLE axis (17,900 genes), the same 8 methods, the same three anchors (so the same method
  weights), the same Borda consensus and the same unanimity flag.
- **Known answer: PASS on all 4 runs.** The A/B/C rows reproduce a fresh `def_identity.match()` on `data/challenge`
  exactly: method_quality, top-25 consensus, call and unanimity (JURKAT / HELA / CAL33, 8/8 each).
- **Flex-like read (the peer's caveat: MIX-seq is 3', D/E/F are Flex).**
  - `flexbias590.py` estimates a per-gene "Arc Flex − 3'" log2 bias. Each profile's deviation from its own CCLE profile
    cancels the line identity; the bias is A/B/C's mean deviation minus the MIX-seq lines' mean deviation.
  - It is cross-fitted in 2 folds over all 200 MIX-seq lines: a fold's bias never uses its own lines.
  - The bias is large: median |b| 0.95 log2 per gene, sd 1.5. It is shared across A/B/C (their deviations correlate
    0.67–0.78), so it is a platform/lab signal, not line noise.
  - It is added to every non-Flex profile before re-scoring.

| profiles (plain read) | n | true line in CCLE | top-1 acc | P(unanimous) | acc \| unanimous | switch-ons | FALSE switch-ons | median cells |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| A/B/C (Flex) | 3 | 3 | 1.000 | 1.00 | 1.000 | 1 | 0 | 18,400 |
| MIX-seq, control cells (10x 3' v2) | 166 | 165 | 0.994 | 0.96 | 1.000 | 3 | 0 | 152 |
| MIX-seq, all conditions | 200 | 199 | 0.990 | 0.96 | 0.995 | 4 | 0 | 521 |
| other single-cell sets (3', 5', sci-RNA-seq3, TAP-seq, CROP-seq, scifi) | 23 | 17 | 0.824 | 0.48 | 0.889 | 9 | **1** (Xu HEK293 → HAP1) | 9,627 |
| donor b_ctrl profiles (anchors, K562, Orion, RPE1, H1, HepG2) | 10 | 6 | 1.000 | 0.40 | 1.000 | 4 | **1** (Orion HEK293T → HAP1) | — |

| read | top-1 acc | acc \| unanimous | switch-ons | false | correct registry switch-ons (of 30 registry rows) | sister false |
|---|---:|---:|---:|---:|---:|---:|
| plain, today's rule | 0.985 | 0.995 | 21 | **2** | 19 | 0 / 146 |
| plain, + gate 0.65 | 0.985 | 0.995 | 17 | **0** | 17 | 0 / 146 |
| Flex-simulated, today's rule | 0.974 | 0.995 | 20 | **2** | 18 | 0 / 146 |
| Flex-simulated, + gate 0.65 | 0.974 | 0.995 | 16 | **0** | 16 | 0 / 146 |

(402 readable profiles in both reads, 390 with the true line in CCLE; `work/summary590.json`.) The registry-line rows that
stay off under today's rule are non-unanimous K562/Jurkat sets with partial gene axes (8–9k genes) or odd chemistries,
and two duplicate-line calls (MCF-7 → KPL-1). A D/E/F context carries the full axis, like A/B/C, all three unanimous.

**Not in CCLE (the abstention cases), plain read.**
- HEK293T (Orion) → HAP1, 8/8, 0.606; HEK293 (Xu) → HAP1, 8/8, 0.554.
- Primary T cells (Shifrut) → CCLF_HEME_0004_T, 8/8, **0.719**.
- Arc H1 (Flex) → TERA2, 7/8, 0.665; iPSC (Tian) → SUSA, 4/8, 0.671; iPSC neurons → CHP-134, 7/8, 0.642.
- Patient melanoma → GAMG, 5/8, 0.661; patient GBM → SNU-201, 6/8, 0.565; RPE1 → SW982, 7/8, 0.466.

## Reading for D/E/F
1. **Which rows resemble D/E/F.** Only A/B/C (Flex, Arc lab, 18,400 cells, full axis) are D/E/F-shaped. The MIX-seq rows
   share neither chemistry nor depth (median 152 control cells). The Flex-simulated rows fix chemistry, not depth, and
   fewer cells only make a correct call harder. So the 0.65 floor is set against A/B/C's 0.719–0.772 (Flex), with the 3'
   q01 (0.681) as a second check.
2. **The residual risk the gate does not cover.** A line absent from CCLE whose nearest CCLE line is a close biological
   relative can score as high as a correct call: primary T cells reach 0.719 against a T-cell lymphoma line. A non-CCLE
   T-ALL line would likely be called JURKAT at a similar value and pass. The gate catches distant non-CCLE lines (HEK293:
   0.55–0.61), not close relatives. CCLE's 1,600 lines cover most commercial cancer lines; HEK293 and RPE1 fail because
   they are not cancer lines. **On Oct 22, write each context's `min_top1_score` into the round's RESULTS. For a value in
   0.60–0.70, read the consensus top-5 and its lineages before accepting a same-line switch.**
3. **The size of what is guarded.** A false call to a thin bank (HAP1: 9/300 panel targets in Jiang, no programme) is
   small. A false call to HCT116 (290/300) or K562 (272/300) would put 270+ wrong rows in at w 0.5. Neither occurred
   anywhere here: no colon or blood line was called HCT116, DLD-1 or K562 unless it was that line.

## The patch (UNAPPLIED) — `patch_identity590.py`
- `def_identity.py`: `match()` also returns `min_top1_score` per context (min over the 8 methods of the top-1 score it
  already computes) and prints it. No other key or table changes; `--verify` output is identical apart from that key.
- `def_sameline.py`: `--min-top1` (default 0.65; 0 disables).
  - A registered, unanimous call below the floor leaves the context unchanged; the value is recorded in its provenance
    `mode`.
  - `--force` rows are untouched. An identity file without the key (pre-590) skips the gate with a WARNING.
  - A passing call writes nothing new: same npz arrays and provenance as today.
- `def_day.sh` needs no edit (the default lives in `def_sameline.py`).
- **Order:** after 538's applier and 554 on 10-05, or any later day. 530 also edits `def_sameline.py`; 590's anchors are
  disjoint, and `--check` proves they hold on both the live and the post-530 file.
- **Tests (`test_identity590.py`, on copies; ALL PASS):**
  - T1: patched `--verify` PASS; output equals the live script's apart from the new key; values A/B/C 0.7716 / 0.7191 /
    0.7329.
  - T2a–f: blocked at 0.554 with the 590 mode; on at 0.70; on with a WARNING when the key is missing; on with
    `--min-top1 0`; `--force` untouched; the non-unanimous rule unchanged. Every "on" case's npz equals the live script's.
  - T3: `--check` on the live and the post-530 file.
- **Inert on every known-answer path.** A/B/C and their relabelled bundles reach the gate only for A (JURKAT, a registry
  line), at 0.772, so it passes and writes nothing new. 451's KA rows use `--force`.
- Commands:
  `uv run python experiments/590-identity-validation/patch_identity590.py --check`, then `--apply`, then
  `uv run python experiments/590-identity-validation/test_identity590.py`.
- The test run included one real catch: the first draft crashed on a pre-590 identity file (a log line formatted
  `None`). It was fixed before anything was handed over.

## Not read, and why
`sciplex3` (gene symbols stored as `nan:<k>`), Cui / Liang / Gehring (≤ 15 official genes mapped), Xie (no var names),
KOLF2.1J (b_ctrl gives zero scores). All are listed in `analyse590.py`'s `UNREADABLE`; none is a registry or a D/E/F-shaped
case.

## Files
| file | what |
|---|---|
| `PREREGISTRATION.md` | question, instrument, decision rule (written before scoring) |
| `profiles590.py` | control pseudobulks per labelled line (h5py, chunked) → `work/profiles/*.npz` |
| `match590.py` | Oct 22's matcher on every profile + known answer; `--flexbias` for the Flex-like read |
| `flexbias590.py` | cross-fitted per-gene Arc Flex − 3' bias → `work/flexbias590_fold{0,1}.npz` |
| `analyse590.py` | the preregistered read → `work/summary590.json` |
| `patch_identity590.py`, `test_identity590.py` | the unapplied gate and its 8 tests |
| `work/match590*.csv` | one row per profile: call, unanimity, `min_s1` (= `min_top1_score`), gaps, true line |

## Addendum (2026-10-01 ~15:40Z, found while reading 596's Flex inventory)
- `donor:viperturb_k562` is Satija's genome-wide VIPerturb-seq K562, measured on **10x GEM-X Flex v2**. It is a real
  Flex profile, but its chemistry label ("VIPerturb-seq") lacks "Flex". So `match590.py --flexbias` added the
  simulated Flex-minus-3' bias on top of real Flex, and its Flex-simulated row (min_top1 0.684) is an artefact.
- Its plain read, **0.765, unanimous K562**, is the true Flex value. That makes a fourth genuine Flex data point for a
  correct in-CCLE call, next to A/B/C's 0.772 / 0.719 / 0.733.
- No conclusion changes. The row was a correct switch in both reads, the floor's two false calls are untouched, and
  removing the artefact row moves the Flex-simulated correct-unanimous quantiles by less than 0.001.
- The real-Flex picture for the floor (0.65):
  - correct calls 0.719–0.772 (four profiles, two labs, Flex v1 and GEM-X Flex v2);
  - the only Flex line absent from CCLE (H1) gets a split call at 0.653–0.665, which never switches anyway.

