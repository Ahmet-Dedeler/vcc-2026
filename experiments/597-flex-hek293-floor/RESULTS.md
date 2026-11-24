# 597 — On real Flex the 0.65 floor is CONFIRMED for HEK293, but by 0.003; HEK293 gets a split call (5/8), so unanimity also blocks it; the simulated Flex bias in 590 was not conservative for absent lines

2026-10-01, session "Explore contribution areas [57577c]", approved by the coordinator. **Ahmet OK'd the download
(17:2xZ).** Preregistered before the download (`PREREGISTRATION.md`, e2b4e80e), with the columns recorded before any
pseudobulk (782d962c). Light: max RSS 3.4 GB (the first attempt was stopped at 4.9 GB by my watchdog: it hashed the 1.58
GB file in one read, now streamed). No refit; outside the coordinator's windows.

## Data
- `multimodal_cell_line_mixing_pilot.rds`, Zenodo 18460279 (Satija lab, VIPerturb-seq), CC BY 4.0.
  - 1,579,169,440 bytes, md5 `3f910a17baa164df2bd66bfac3197360` (equal to Zenodo's record).
  - Kept in this session's scratchpad, not committed.
- 10x Flex v1. Three KRAB-dCas9-MeCP2 lines pooled and separated by hashing.
- Parsed with `rds_stream597.py`, which is 424's streaming reader plus R byte code (two byte-code items skipped; clean
  end of stream).
- Non-targeting cells (`Gene_assignment == "NO.TARGET"`) per `hash.ID`:
  - K562 498 of 3,437 cells, 2.62 M UMI;
  - HAP1 478 of 3,379, 3.10 M;
  - HEK293 486 of 4,297, 5.56 M.
  Integer counts, RNA `counts` layer.

## Known answer: PASS
- The 8 methods' weights computed here equal 590's recorded `method_quality` exactly.
- A/B/C re-scored in the same run reproduce 590's calls, unanimity and min_top1 (0.772 / 0.719 / 0.733).

## Result (590's matcher unchanged)
| line | call | methods agreeing | **min_top1** | consensus top 5 |
|---|---|---:|---:|---|
| K562 (in CCLE, registry) | K562 ✓ | 8/8 | **0.693** | K562, P2URK562, KYO-1, JURL-MK1, NCO2 (all CML) |
| HAP1 (in CCLE, registry) | HAP1_ENGINEERED ✓ | 8/8 | **0.747** | HAP1, SNU-398, C-33A, EN, TOV-112D |
| **HEK293** (absent from CCLE) | SNU398_LIVER | **5/8 (split)** | **0.6472** | SNU-398, **HAP1**, TOV-112D, HOS, C-33A |

## Reading (the preregistered rule)
- **HEK293 min_top1 0.6472 < 0.65 → the floor is CONFIRMED on real Flex.** The margin is 0.003. The call is also split
  (5/8, SNU-398 over HAP1), so the unanimity requirement blocks it independently. On real Flex the exact 3′ failure (a
  unanimous HAP1 call) does not occur.
- HAP1 (0.747) and K562 (0.693) are called correctly and clear the floor, so on this object the floor costs nothing.

## What it changes in what we know (beyond the binary reading)
1. **Real Flex raises an absent line's match quality.**
   - HEK293 reads 0.647 here, against 0.554 / 0.606 in 3′, 0.616 in Orion's FiCS 5′ file (595), and **0.506 / 0.561
     under 590's simulated Flex bias**.
   - The simulated bias moved absent lines the wrong way, by about 0.09–0.14. 590's Flex-simulated read was therefore
     NOT conservative for the false-positive side. Its conclusion survives only because the real value still sits under
     0.65.
2. **The real-Flex gap between correct calls and absent lines is narrow.**
   - Correct real-Flex calls: A 0.772, B 0.719, C 0.733, VIPerturb GEM-X v2 K562 0.765, and here HAP1 0.747, K562 0.693.
     The lowest is 0.693.
   - Absent lines on real Flex: HEK293 0.647 and H1 0.653–0.665. Both are split calls.
   - Above 0.65 there is room up to ~0.69 without blocking any correct real-Flex call.
   - The preregistered rule does not call for raising the floor. The coordinator may judge whether to raise it to
     ~0.67–0.68 for margin. The cost would be degraded non-Flex profiles such as 595's THP-1 at 0.665, which a Flex D/E/F
     context should not resemble.
3. **The day-of reading rule matters.** "Log min_top1 and read the top 5 when it is 0.60–0.70" would flag both HEK293
   (0.647, with HAP1 second) and K562 v1 (0.693) for a human look. Together with the unanimity requirement, it is what
   actually protects the switch in this chemistry.

## Caveats
- Flex v1 with custom probes and a non-standard fixation (1 % formaldehyde, 10 min).
- CRISPRi machinery is present in every cell.
- About 480–500 non-targeting cells per line, against 18,400 on the day: fewer cells make a profile noisier, which
  lowers rather than raises the scores.
- HEK293 here is HEK293, not HEK293T.

## Files
| file | what |
|---|---|
| `PREREGISTRATION.md` | the rule, written before the download; the columns, recorded before any pseudobulk |
| `rds_stream597.py` | 424's streaming reader plus R byte code |
| `score597.py` | `inspect` (structure) and `score` (pseudobulks + 590's matcher + known answer) |
| `work/inspect597.json`, `work/score597.json`, `work/*.log` | structure, results (incl. the file's sha256), logs of both attempts |

## Addendum (2026-10-02): the pseudobulks themselves are saved
- 597's score run kept only the scores. `pseudobulk597.py` rebuilds the per-line NTC pseudobulks with the same RDS, parser
  and selection, for exp 620 (council session, identity across chemistries).
- Known answer PASS: the RDS sha 3218665c…, the cells 498 / 478 / 486 and the totals 2,623,903 / 3,103,451 / 5,555,320
  equal the 597 run. Peak 0.92 GB.
- `work/pseudobulk597.npz` (516 KB) holds:
  - `lines`, `n_cells`;
  - `genes` (16,336 RNA features) and `counts` (3 × 16,336 integer sums);
  - `log1p_cp50k` over the object's own features;
  - `official_genes` and `counts_official` (3 × 18,533, NaN where the probe set lacks a gene; 16,333 present). This is
    exactly what 597 passed to 590's matcher.
  - The sidecar `work/pseudobulk597.json` carries its sha256.
- Derived from CC BY 4.0 data (Satija lab, Zenodo 18460279). Caveats as above: Flex v1 with custom probes; HEK293, not
  HEK293T.
- **Used by 620d/e** (council session [7772dc], `experiments/620-bank-identity/RESULTS.md`, c0afc398).
  - Renormalised over shared genes and lab-centred, Satija K562 → Weissman k562_gwps matches: +0.28–0.29, margin
    0.16–0.23.
  - Satija HEK293 → Orion HEK293T does **not** clear the negatives: +0.17–0.20, rank 1–2, inside the negatives'
    0.13–0.21.
  - So a bank's control profile cannot identify a HEK293(T) context on Flex. The 0.65 CCLE floor plus unanimity stays
    the guard: rows off.
  - Caveat: D/E/F are cancer lines (FINDINGS H.1, from the 2026 Cell paper), and every non-CCLE bank line (HEK293T,
    PSC, RPE1, TeloHAEC, CD4 T) is non-cancer. HEK293 is therefore a stress test of the absent-line branch, not a likely
    D/E/F identity.

