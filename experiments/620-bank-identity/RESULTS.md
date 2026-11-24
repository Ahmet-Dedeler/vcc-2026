# 620 — Bank-line identity across chemistries: the day-of matcher can find a line's bank without CCLE

2026-10-02, council session [7772dc]. Analysis only: no build, no `vcc`, no `cell-eval2`, no download. Peak RSS < 1 GB,
runtime seconds. Scripts: `match620.py` (raw), `match620b.py` (with chemistry correction). Outputs: `work/match620.json`,
`work/match620b.json`.

## Question
`def_identity.py` matches a context's controls against CCLE bulk RNA-seq, and 590/597's gate (registered AND unanimous
AND min_top1 ≥ 0.65) refuses lines absent from CCLE. Several non-CCLE cell types have large perturbation banks on disk
(Orion HEK293T genome-wide, Arc's H1, Pan's H1, RPE1, TeloHAEC, CD4 T). Can a context be matched to a bank line's own
control pseudobulk (`b_ctrl.npy`, log1p CP50k) across study, lab and chemistry?

## Method
- References: every bank whose `b_ctrl` is a real control profile. Placeholders are dropped: jiang_*, feng_ipsc,
  altair_a549, dld1_lowmoi and igvf_cardiomyocyte are all-zero; kolf2 is all-one.
- Common gene axis: genes in ≥ 85 % of profiles. Genes are standardised by the reference set's mean and SD (mean > 1,
  SD > 0.1). Score: Spearman.
- Flex probes are corrected by a per-gene Flex-minus-other bias, fitted leave-one-pair-out on the known same-line pairs:
  K562 (viperturb_k562 vs k562_gwps), PSC (h1_2025 vs pan_h1), Jurkat (challenge A vs nadig_jurkat).

## Results
| probe | truth | top match | score | 2nd | same-line rank |
|---|---|---|---:|---|---:|
| challenge A (Flex, Arc) | Jurkat | nadig_jurkat | +0.424 | k562_gwps +0.150 | 1 |
| h1_2025 (Flex, Arc) | H1 | pan_h1 | +0.383 | k562_gwps +0.191 | 1 |
| viperturb_k562 (Flex) | K562 | gasperini_k562 / k562_gwps | +0.406 / +0.392 | nadig_jurkat +0.168 | 1 |
| challenge B (Flex) | HeLa, no bank | rpe1 | +0.196 | telohaec +0.102 | none |
| challenge C (Flex) | CAL33, no bank | telohaec | +0.196 | sunshine_calu3 +0.173 | none |
| xucao_hek293 (3′; 620 raw) | HEK293 | orion_hek293t | +0.517 | kolf2 +0.190 | 1 |

- Without the chemistry correction (620 raw), every Flex profile matches the other Flex profiles first (A/B/C →
  h1_2025, viperturb_k562, cd4). The correction is necessary.
- With it, all four known-answer probes match at rank 1, with margins 0.19–0.27. The two no-bank negative controls stay
  at ≤ 0.196. A floor near 0.30 separates matches from non-matches on this set (n = 4 positives, 2 negatives: small).

## Reading
A non-CCLE D/E/F line that has a bank (HEK293(T) → orion_hek293t, PSC → h1_2025/pan_h1, RPE1, TeloHAEC, primary T →
cd4) could be identified on the day and switched to same-line rows. Today the pipeline gives such a line nothing. That
is the 450 mechanism (H1: +0.053 … +0.079 overall-equivalent at 35 % coverage). Class: day-of instrument plus transfer;
nothing is fitted to A/B/C. Not yet decisive for HEK293 in Flex. Next: a Flex HEK293 probe (597's Satija data), more
negative controls, a relabelled-bundle rehearsal, and a real kolf2 control profile. Council proposal H-620 in
`experiments/600-breakthrough-council/7772dc.md`.

## 620c — negative controls with the truth removed (`neg620c.py`, `work/neg620c.json`)
Each Flex probe whose line has a bank is re-scored with every reference of its own line removed, and with its own pair
removed from the bias fit. This is the regime a D/E/F line with no bank would face.

| probe (truth removed) | best wrong line | score | 2nd | margin |
|---|---|---:|---|---:|
| challenge A (Jurkat) | k562_gwps | +0.192 | pan_h1 +0.160 | 0.03 |
| h1_2025 (PSC) | k562_gwps | +0.224 | nadig_jurkat +0.138 | 0.09 |
| viperturb_k562 (K562) | nadig_jurkat | +0.257 | pan_h1 +0.226 | 0.03 |
| challenge B (HeLa) | rpe1 | +0.196 | telohaec +0.102 | 0.09 |
| challenge C (CAL33) | telohaec | +0.196 | sunshine_calu3 +0.173 | 0.02 |
| cd4_rest (CD4 T) | wessels_thp1 | +0.306 | nadig_jurkat +0.268 | 0.04 |

- True matches score +0.383 … +0.424 with margins 0.19–0.27 to the next line. Wrong lines score ≤ +0.306 (the
  hematopoietic neighbours of CD4 are the closest call) with margins ≤ 0.09.
- **Candidate day-of rule:** top ≥ 0.33 AND margin to the next different line ≥ 0.15. It separates every row here, but
  n is small: 3 cross-chemistry positives and 6 negatives. Widen it (Flex HEK293, more Flex lines) before trusting it.

## 620d — renormalised over shared genes, plus Satija's Flex probes (`match620d.py`, `work/match620d.json`)
[36a81a] rebuilt 597's Satija Flex NTC pseudobulks (`597/work/pseudobulk597.npz`: K562, HAP1, HEK293; Flex v1). Every
profile is renormalised over the shared gene axis, and the bias, scoring and rule are as in 620b/620c.

| probe | kind | top (score) | margin | rule (≥ 0.33 and margin ≥ 0.15) |
|---|---|---|---:|---|
| challenge A | positive (Jurkat) | nadig_jurkat +0.430 | 0.285 | fires, correct |
| h1_2025 | positive (PSC) | pan_h1 +0.393 | 0.188 | fires, correct |
| viperturb_k562 | positive (K562) | gasperini_k562 +0.455 | 0.296 | fires, correct |
| satija_K562 | positive (K562) | k562_gwps +0.280 | 0.157 | quiet (below floor): miss |
| **satija_HEK293** | near-positive (HEK293 vs 293T) | **k562_gwps +0.201; orion_hek293t +0.199** (rank 2) | 0.002 | quiet: miss |
| satija_HAP1 | negative (no real bank profile) | k562_gwps +0.236 | 0.065 | quiet |
| B, C, cd4, four truth-removed probes | negatives | ≤ +0.296 | ≤ 0.089 | quiet |

- No false positives in 9 negatives. 3 of 5 positives fire. **The HEK293 case, the one with a genome-wide bank, does
  not match.** All three Satija probes lean toward k562_gwps, which points to a lab/chemistry offset the cross-lab bias
  does not remove.

## 620e — within-lab centring (`lab620e.py`, `work/lab620e.json`)
Each lab group with ≥ 2 lines is centred on its own mean: Arc {A, B, C}, Satija {K562, HAP1, HEK293}, Weissman
{k562_gwps, rpe1, nadig_jurkat, nadig_hepg2}, Orion {hct116, hek293t}. On Oct 22 the Arc group is A–F (six lines).

| probe | top | score | 2nd | margin |
|---|---|---:|---|---:|
| challenge A | nadig_jurkat | +0.454 | k562_gwps +0.062 | 0.39 |
| satija_K562 | k562_gwps | +0.286 | orion_hct116 +0.060 | 0.23 |
| satija_HEK293 | orion_hek293t | +0.173 | rpe1 +0.120 | 0.05 |
| challenge B (neg) | rpe1 | +0.205 | orion_hek293t +0.083 | 0.12 |
| challenge C (neg) | nadig_hepg2 | +0.185 | rpe1 +0.116 | 0.07 |
| satija_HAP1 (neg) | orion_hct116 | +0.132 | nadig_hepg2 +0.065 | 0.07 |

- Within-lab centring removes the lab offset: Satija K562 now matches Weissman K562 cleanly. HEK293 → Orion HEK293T
  reaches rank 1, but at +0.173, inside the negatives' range (0.13–0.21). **Not decisive.** HEK293 is not 293T, and
  Orion's group has only two lines, so its contrast is just (HEK293T − HCT116)/2.

## Reading (updated)
- Bank-line identity works across labs and chemistries for Jurkat, K562 and H1 (PSC). With ≥ 3-line lab groups or a
  fitted chemistry bias, true matches score 0.29–0.45 with margins ≥ 0.16; negatives stay ≤ 0.31 with margins ≤ 0.12.
- **For the HEK293(T) case (the genome-wide bank) the control-profile route is not reliable.** A better route is to put
  a bulk HEK293 / HEK293T profile into def_identity's own reference, beside CCLE. The Human Protein Atlas cell-line RNA
  atlas covers many non-cancer lines; it is a small public download, which needs Ahmet's OK. The existing, validated
  matcher would then call it, and the registry maps HEK293 → orion_hek293t.
- The PSC branch (Arc H1 → Pan H1) is decisive enough to wire: same-line rows from h1_2025_full + pan_h1 for a
  pluripotent context (coverage 25 + 64 of A/B/C's 300).
