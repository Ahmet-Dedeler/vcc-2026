# 645 — HPA's bulk cell-line atlas as non-CCLE identity references: it independently confirms 642; a broad veto panel is rejected

Date 2026-10-02, session [537757]. Analysis only: no build, no `vcc`, no `cell-eval2`. Peak RSS 3.35 GB, 294 s.

**Download** (Ahmet's standing OK, 10-02):
- `rna_celline.tsv.zip`, https://www.proteinatlas.org/download/tsv/rna_celline.tsv.zip, 205,861,141 bytes, sha256
  `e5b11779…`. HPA v23 (2023-05-24), 1,206 cell lines × 20,151 genes, TPM, CC BY-SA 4.0.
- Kept in `data/external/hpa/` (gitignored). The cache is `hpa_cache645.npz` (log2(1 + TPM), symbol-summed), built by
  `load_hpa645.py` in 3 s.

**Class:** day-of instrument (642's second stage). Nothing is fitted to A/B/C.

## Why

- 642's HEK293 reference was built from single-cell, non-Flex HEK293 controls (XuCao 3′ + Orion HEK293T FiCS).
- A bulk reference from a different source tests whether the real-Flex HEK293 call depends on how the reference was
  made.
- HPA's non-cancer lines (hTERT-RPE1, HaCaT, hTERT-HME1, RPTEC/TERT1, …) also allow a test of a broader "veto panel"
  for D/E/F lines absent from CCLE.

## Results (`match645.py`; 590/642's matcher with stage-1 weights frozen; `work/match645.json`)

**1. HPA bulk is a CCLE-compatible platform.**
- 973 HPA profiles of lines also in CCLE, scored CCLE-only: top-1 accuracy **0.997**, unanimous 0.998.
- min_top1 of the correct calls: q01 0.761, q05 0.795, median 0.825.
- The 3 misses are sister or duplicate lines: EB2 → EB1, NCI-H1882 → NCI-H1876, SuSa → TERA2.

**2. HEK293 from HPA: the same answer as 642's row.**

| reference row | Orion HEK293T (FiCS) | XuCao HEK293 (3′) | **real-Flex HEK293 (Satija)** | false accepted (of 423) |
|---|---|---|---|---:|
| HPA bulk "HEK293" (all three held out) | 8/8, 0.785 | 8/8, 0.734 (below the 0.75 stage-2 floor) | **8/8, 0.822** | 0 |
| 642's shipped non-Flex row (XuCao + Orion) | 8/8, 0.970 (in-sample) | 8/8, 0.972 (in-sample) | **8/8, 0.820** | 0 |

- On the D/E/F platform (Flex) the two independent references give the same call and the same score (0.822 vs 0.820).
- The real-Flex HEK293 call does not depend on the reference's source.
- HPA's row scores single-cell 3′ lower (0.73): a bulk-vs-3′ platform gap that does not matter for Flex contexts.
- **No change to 642.** Its row stays: equal on Flex, wider margins on single-cell profiles, already tested and
  sha-pinned.

**3. A broad veto panel of HPA non-CCLE lines: REJECTED.**
- 17 rows were added at once: HEK293, HEK TE, hTERT-RPE1, HaCaT, hTERT-HME1, RPTEC/TERT1, HBEC3-KT, TIME, HUVEC/TERT2,
  three BJ variants, fHDF/TERT166, two PODO, hTEC/SVTERT24-B, HBF/TERT88.
- Accepted calls (8/8, ≥ 0.75): the two held-out HEK293 profiles → HEK293 (correct).
- **But there are 2 FALSE vetoes:** MIX-seq BICR31 → HaCaT (0.755) and SNU1076 → HaCaT (0.772). Both are head-and-neck
  squamous carcinomas, CAL-33's (context C's) tissue.
- An immortalized keratinocyte row attracts squamous carcinomas. Close-relative confusion is exactly what such a panel
  would add.
- No registry line was vetoed, but the false-veto rate on true CCLE lines must be 0 for a panel that can override a CCLE
  call. HEK293 alone stays the only non-CCLE row (642).

## One sentence

An independent bulk HEK293 reference gives the same unanimous 0.82 call on real-Flex HEK293 as 642's shipped row,
confirming 642. Adding many non-cancer lines at once is unsafe: HaCaT pulls squamous carcinomas, CAL-33's tissue.

## Reproduce

```bash
uv run python experiments/645-hpa-references/load_hpa645.py   # needs data/external/hpa/rna_celline.tsv (unzip the download)
uv run python experiments/645-hpa-references/match645.py
```
