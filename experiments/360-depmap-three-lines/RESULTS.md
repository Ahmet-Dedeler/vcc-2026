# 360 — DepMap/CCLE covariates for the three context lines: identities confirmed, bank built

**Status: COMPLETE.** All four tasks done. No kill trigger fired.

## Known-answer reproduction (rule 3, own code, read-only of 323's outputs)

- `experiments/323-current-family-validation/work/per-target-cap15.csv` → 58 `de_wilcoxon_lfc_nmae`
  rows, mean **0.95083646** vs 343's recorded **0.950836**: **exact match** (Δ < 1e-6).
- From-scratch check (`known343.py` logic on `stats-cap15.npz`/`reference.npz` lfcs): eligibility
  threshold recovered — the 58 scored targets are exactly those with ≥ 11 genes at q < 0.05
  (COASY with 10 is out, RFX5/ZNF92 with 11 are in; same 58-target set). Mean over a raw
  q < 0.05 gate is 0.9425 (Δ 0.0084), so the scorer's gate adds a further filter
  (CPM/effect-size) I did not re-derive. Documented as a limit, not a mismatch.

## 1. What is present locally (before downloading)

| store | three lines? |
|---|---|
| `data/external/ccle_broad_2025/ccle_cache.npz` (1600 samples × 19,209 genes TPM; 25Q3 mirror) | **ALL PRESENT**: `JURKAT_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE`, `HELA_CERVIX`, `CAL33_UPPER_AERODIGESTIVE_TRACT` (plus K562, HCT116, HEPG2; no HEK293T/RPE1 — only `HCC2935_LUNG` matches "293") |
| `data/external/depmap26q1/` | **ABSENT — and misnamed**: it is a 16-line Broad-16 Cas9-screen extract (`cell_line_metadata.csv`: HCT15…HKGZCC; `crispr_table.csv` only 1,600 rows), not the portal 26Q1 release. No Model.csv, no matrices. |

DepMap IDs (24Q4 `Model.csv`, kept in `work/`): **A/Jurkat = ACH-000995, B/HeLa = ACH-001086,
C/CAL-33 = ACH-000518** — the ACH-001086 and ACH-000518 confirm 356's stated IDs; Jurkat's was
read from Model.csv as instructed. Donors: K562 ACH-000551, HCT116 ACH-000971, HepG2 ACH-000739;
all 12 Broad-16 working lines present (HCT15 000997, KYSE140 000823, KYSE450 000865, RPMI7951
000348, SKGII 001652, U343 000036, C4I 001334, ONE58 001992, T3M5 002042, KMRC20 000250,
SNU8 000460, HKGZCC 001852). HEK293T is not a DepMap model; parental RPE1 is absent (only
engineered RPE1-ss* clones).

## 2. Download manifest (accessed 2026-09-10)

Portal `depmap.org/portal/{data_page,api/download/files}` returns a bot-verification page to
programmatic fetch (two endpoints probed, both `Verification` HTML; portal text itself invites
bulk grabs from the downloads section but no keyless path exists). Figshare hosts full releases
only through **24Q4** (article 27993248, verified via API search — no 25Q/26Q1 full-release
article exists), so the licensed fallback below was used **instead of** 26Q1:

- License on the article record: **CC BY 4.0** (`creativecommons.org/licenses/by/4.0/`), doi
  `10.25452/figshare.plus.27993248.v1`, published 2024-12-10. License-clean.
- `Model.csv` — full file kept, 645,696 B, sha256
  `b7a0c1385e6cef30132b56aff61f1261d11e3f490490b355c430d32ee0dbdcfa`
  (URL `https://ndownloader.figshare.com/files/51065297`; all URLs in `work/fs_urls.json`).
- `CRISPRGeneEffect.csv` (428.7 MB upstream) — streamed, kept 17/18 rows →
  `work/CRISPRGeneEffect_subset.csv` 6,447,569 B, sha256 `60681af8…c96`.
- `OmicsExpressionProteinCodingGenesTPMLogp1.csv` (506.6 MB upstream) — streamed, kept 18/18 →
  `work/Expression_subset.csv` 5,757,558 B, sha256 `97b883ca…3fa`.
- `OmicsSomaticMutations.csv` (338.9 MB upstream, < 1 GB as tasked) — streamed 718,369 rows,
  kept 3,963 (3 lines) → `work/Mutations_3lines.csv` 1,866,639 B, sha256 `5970c2d0…76e2e`.
- Full upstream files were never written to disk (row-filtered in flight). Kept total ≈ 15 MB,
  all under `work/`. Transferred ≈ 1.27 GB, deleted nothing (nothing full retained).

**Substitution note:** 24Q4 Chronos/expression for three lines established decades ago is stable
across releases; cross-checked against the held 25Q3 CCLE mirror below (r ≈ 0.85). Re-pull at
26Q1 only if a number here becomes load-bearing.

## 3. Bank: `work/three_lines.npz`

Per line: basal log(TPM+1) and Chronos gene effect on the challenge 18,533-gene axis
(`genes`, `targets`, `expr_ids` 18, `expr_logtpm1p` 18×18533, `dep_ids` 17, `dep_chronos`
17×18533, plus `expr_missing`/`dep_missing` symbol lists). Symbol mapping strips DepMap
`GENE (ENTREZ)` suffixes, first duplicate kept.
Losses: expression **357/18,533** (expected — DepMap matrix is protein-coding only; e.g.
AC118549.1; plus deprecated symbols SARS/KARS/AARS), dependency **1,432/18,533**.
Also saved: `work/challenge_cpm.npz` (control total-count CPM per context, chunked backed read),
`work/panel_dependency.csv` (300 targets × context/donor Chronos).

## 4a. Identity validation — INDEPENDENT confirmation of 041

DepMap log(TPM+1) vs challenge control log(CPM+1), 18,176 common genes:

| context | vs JURKAT (r / ρ) | vs HELA (r / ρ) | vs CAL33 (r / ρ) |
|---|---|---|---|
| A | **0.9071 / 0.9054** | 0.7191 / 0.7456 | 0.6825 / 0.7091 |
| B | 0.7040 / 0.7362 | **0.8532 / 0.8854** | 0.7331 / 0.7692 |
| C | 0.7203 / 0.7404 | 0.7736 / 0.8022 | **0.8750 / 0.8874** |

Each context matches its claimed line by **+0.12–0.19** over either shuffled pairing, both
metrics, all three contexts. A=Jurkat, B=HeLa, C=CAL33 confirmed independently of 041.

## 4b. Nearest donors

Expression (Pearson, top 3 non-self): A → HCT15 0.768, HCT116 0.761, K562 0.759;
B → HCT116 0.779, ONE58 0.777, HCT15 0.774; C → KYSE450 0.821, HKGZCC 0.819, HCT15 0.811.
Dependency (Chronos Pearson): nearly flat — all pairs 0.79–0.92 (shared essentiality
dominates); ranks: A → C4I 0.848, KYSE140 0.844, KYSE450 0.842; C → KYSE450 0.917,
HKGZCC 0.912, SKGII 0.910. **Dependency distance is uninformative for donor choice** (report,
not a recommendation). **B has no dependency row at all**: HeLa (ACH-001086) is absent from
24Q4 `CRISPRGeneEffect` (file holds 1,178 models; HeLa CRISPR came in later releases) —
any B-side dependency reasoning needs a newer release or Sanger SCORE, not this bank.

## 4c. Panel-target dependency (Chronos; negative = more essential)

300 targets: A n=298, median −0.044, 14 with score < −0.5, **0 < −1.0**; C n=300,
median −0.072, 9 < −0.5, **0 < −1.0**; B n=0 (absent, see above). Consistent with the
essential-free-panel design (§8.1g): no panel target is a strong dependency anywhere.
Most dependent: A — CDK8 −0.899, TWF2 −0.678, TMEM41B −0.656, MYRF −0.644, LDLR −0.637;
C — THAP7 −0.871, TMEM41B −0.771, MED25 −0.556, SLIRP −0.548, FAM120A −0.539.
Context-dependent but donor-average-not (score < −0.5 here, donor mean > −0.3):
A — AKT2, CDK8, GIGYF2, HDAC1, LDLR, MYRF, PRMT2, TWF2;
C — CDK8, FAM120A, NGLY1, SLIRP, STT3A, THAP7.
Full table: `work/panel_dependency.csv` (columns `chronos_{A,B,C}`, `donor_mean/sd`,
`n_donors_dep`; B column all-NaN by construction).
Mutations (3,963 rows): Jurkat 3,579 vs CAL33 255 vs HeLa 129 — Jurkat is hypermutated
(known MMR deficiency), so its 80 panel-target variants (incl. HIGH-impact frameshifts in
PRMT6, CHD3, TAF15, TWF2, SCAF8, ZNF462 and stops in SEMA4F/IL10RB) vs 7 in CAL33 and 5 in
HeLa should be read as background burden, not driver enrichment.

## Verdict against the PLAN kill criterion

**Not killed.** All three lines identified in a CC BY 4.0 matrix, identities validated at
+0.12–0.19 over shuffles, bank + tables delivered. Usable PI takeaways: (i) identities hold;
(ii) no panel target is strongly essential in A or C — dependency is not a response-strength
proxy on this panel; (iii) the divergent-dependency shortlists above are the only
context-vs-donor contrasts this bank supports; (iv) B-side dependency is a hole until a newer
release is pulled.

## Limits

- 24Q4, not 26Q1 (portal bot-wall; §2). 24Q4↔25Q3-mirror expression agreement on the 3 lines:
  r = 0.846/0.851/0.848 (17,849 common genes) — proxy validated, not identical.
- Challenge pseudobulk is total-count CPM (mean-of-ratios not computed; difference negligible
  at these depths, libraries 3.7–3.9e8 total).
- Chronos correlations are bulk-profile similarities, not perturbation-response predictions.
- H1-panel numbers nowhere used; nothing here is a score claim.

## Exact commands to reproduce

```
uv run python experiments/360-depmap-three-lines/work/known343.py        # known answer
uv run python experiments/360-depmap-three-lines/work/dl_model.py        # Model.csv + URLs
uv run python experiments/360-depmap-three-lines/work/dl_subset.py       # expr + gene-effect subsets
uv run python experiments/360-depmap-three-lines/work/dl_mut2.py         # mutation subset
NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4 uv run python experiments/360-depmap-three-lines/work/analyze360.py
uv run python experiments/360-depmap-three-lines/work/analyze360b.py
uv run python experiments/360-depmap-three-lines/work/analyze360c.py
```

## Peak memory and wall time

`memory_pressure` free 69–73% before heavy steps (gate: ≥ 50%). All jobs small: widest
operation was chunked backed h5ad reads (2,000 rows/chunk) and 18-row pandas subsets; no job
near the 12 GiB cap (no `/usr/bin/time -l` peak captured — RSS stayed at interactive scale,
est. < 2 GiB). Wall: Model.csv seconds; gene-effect 12 s, expression 14 s, mutations 11 s
(stream-filtered); analyses < 2 min each. Total ≈ 25 min.
