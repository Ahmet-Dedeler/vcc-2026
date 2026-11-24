"""533: append the Round 14 reference rows to REFERENCES.md (from the untouched v14 copy). Nothing was fetched from the web
for this revision: each row is transcribed from the experiment that read the source, as its row says."""
import pathlib
HERE = pathlib.Path(__file__).resolve().parent
src = (HERE / "REFERENCES.v14.orig.md").read_text()
ADD = """

## Added for Draft v15 (533, 2026-09-30)

Nothing was fetched from the web for this revision. Each row is transcribed from the experiment that opened the source; the
route, access date and the sentences relied on are in that experiment's `RESULTS.md`.

| key | reference | read by |
|---|---|---|
| **VCC26-CELL** | *Virtual Cell Challenge 2026: Benchmarking zero-shot generalization across cellular contexts.* Cell (2026), doi:10.1016/j.cell.2026.08.004, PMID 42648290, published online 2026-08-26, CC BY-NC-ND 4.0. Six CRISPRi Perturb-seq lines from different tissues on 10x Flex; three validation, three test; no line names. | 520 (full text via the Cell HTML and the Europe PMC abstract, 2026-09-28) |
| **JIANG25** | Jiang *et al.* (Satija laboratory, NYGC; the author list was not transcribed by 520). *Systematic reconstruction of molecular pathway signatures using scalable single-cell perturbation screens* (the GEO series title). bioRxiv 10.1101/2024.01.29.576933 (Nature Cell Biology version the same study); GEO GSE281048, Zenodo 10.5281/zenodo.14518762. Six cancer lines (A549, MCF7, HT29, HAP1, BxPC3, K562), CRISPRi (dCas9-KRAB-MeCP2), 44–61 regulators per pathway. The "Jiang 2025" banks of §3(j) and F-10; the closest published analogue of the challenge's design. | 520 (full text and supplement, 2026-09-28) |
| **ATLASSHIFT** | *AtlasShift* emitter by Kaipeng Zheng (board author, as 525 records), github.com/kaipengm2/Virtual-Cell-Challenge-2026 at commit `d24ce4fdae0cd7cba3ba29546cd8737094eae98a` (2026-09-23); `predict.py`, `model.py` (dual-moment count fitting). The most reproduced public emission design on the 2026 board (§3(h)). | 525 (code read at the pinned commit, blob sizes matched, 2026-09-29) |
| **DEPMAP-CHRONOS** | DepMap 24Q4 CRISPR gene effect (Chronos), `CRISPRGeneEffect` subset for K562 (ACH-000551) and HCT116 (ACH-000971), CC BY 4.0. §3(h)'s dependency test. | on disk: `experiments/360-depmap-three-lines/work/CRISPRGeneEffect_subset.csv`; used by 514 |
"""
out = src.rstrip("\n") + ADD
(HERE.parent / "paper" / "REFERENCES.md").write_text(out)
print(f"REFERENCES v14 {len(src):,} -> v15 {len(out):,} chars")
