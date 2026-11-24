# Board intelligence, 2026-09-26 (peer-review session)

Source: all 176 snapshots in `data/leaderboard/` (2026-08-30 → 2026-09-26 19:44Z), deduplicated by entry id: 11,731
entries from 1,252 teams, 3,934 with a non-empty `description`. The website's note icon beside a model name shows
exactly the API's `description` field (checked on CellWM-57: tooltip "v57"), so the snapshots hold everything the page
shows. Every description is self-reported and unverified; the member columns are measured.

## 1. Who discloses

- The top six are mostly silent: Illumina AI, cqawesome and nailong king post no descriptions; CellWM posts only
  version strings ("v56").
- Two top teams disclosed real detail earlier and stopped: **Testh** (#3; silent since 09-15, every entry named
  "test") and **chill** (#6; silent since 09-13).
- Everyone else who describes a method scores ≤ 0.24. Arc's Stack model, fine-tuned (Yifan SONG, Advaita Labs):
  0.206. The public repo `github.com/kaipengm2/Virtual-Cell-Challenge-2026` (Kaipeng Zheng; reproduced by Yao Li, Vu
  Mai, Jia): 0.180. Flow matching / CellOT: ≤ 0.19.
- Our 27 names are jokes and our descriptions are empty (AGENTS §14 holds).

## 2. Where we stand among teams with a comparable bulk

Among the 95 teams whose latest entry has pds ≥ 0.84, our bulk ranks 3rd (pds) and 4th (mse), and 86 % have a worse nmae
(median 0.923 vs our 0.896). 98 % have a worse jac, 68 % a worse fid, 72 % a worse reach. The DE channel is our relative
weakness against the top six and a few specialists, not against the field.

## 3. Lead 1: the shelved reliability arm (410d2), which Testh's board history supports

- Testh's disclosed steps, board raw members:

| date | change (their description) | pds | mse |
|---|---|---:|---:|
| 09-11 | all-source combined delta, hybrid pinned emitter | 0.8054 | 0.8418 |
| 09-12 | × fixed per-(target, gene) reliability weight `0.15 + 0.85·z²/(z²+8)·n/(n+4)` | 0.8203 | 0.8097 |
| 09-14 | × learned per-(target, gene) gate (sigmoid of 8 truth-free source features: log z², log n_src, log control expression, log\|μ\|, logit source-sign agreement, log cross-target variance, log n_cells, log source MAD; fitted on 6 leave-one-context-out folds for energy-weighted cosine) × per-target weight (ridge on 37 per-target features) | 0.8356 | 0.6985 |
| 09-15 | − 0.5 × projection on the top-14 singular directions of LOCO truth deltas; library-size tilt off | 0.8450 | 0.7055 |

- Their per-target weight is the analogue of our 410 key, which we already ship. The per-(target, gene) reliability part
  is not in our pipeline.
- Ours was priced and shelved: 414 found 228's per-gene reliability gain real in direction, with an amplitude excursion
  (SYVN1-C). The bounded variant (m clipped to [1/2, 2], amplitude-matched, GO rule cliffcheck ratio 0.97–1.03) was
  priced at **+0.0025 to +0.0046 overall** and shelved as "fourth in the queue" on 09-22 (NEXT-SLOT-DECISION l.199).
  It was never built. The objection that closed 224/228 (nmae pays) no longer applies under ballast.
- Suggested: re-price 410d2 on the current stack (q 1.75 sharpening may already take part of it) and, if it holds,
  consider Testh's feature set (cross-source sign agreement, n_src, control expression) for the gain.

## 4. Lead 2: a separately learned DE-membership model (chill, Testh)

- chill posts raw nmae 0.833–0.838 in every entry we have, from its first (09-11), when its bulk was far weaker than
  ours (mse 0.845, pds 0.824). Its nmae therefore does not come from bulk information.
- chill's descriptions (09-11/09-12): "signed conditional response estimates, separate DE membership/direction
  ranking"; "DE membership trained using equal total weight for each source study family, separate directional
  ordering, and bounded joint integer population refinement… Each cell and gene total is preserved, cell-mean LFC drift
  is bounded by 0.01". Testh lists a "DE-membership policy" as a component too.
- Reading: chill predicts which genes will be DE as its own model, sets arithmetic-mean LFCs to the response given DE,
  then refines integer counts for the Wilcoxon ranks with the means pinned. Our call set is the design's own |z| order,
  and our delivered means are ~0.17× the truth plus the MDE add (434, 426).
- Related closed items to check it against: order keys (406), per-target K (385), 434's learned L1 magnitude. A learned
  membership model is a different object from those, but the coordinator should judge the overlap.

## 5. Lower-value observations

- **The pack** (Ava Brown, Anna Xie, Vivian Xu, SQ W, HY H; near-identical members). **CORRECTED 2026-09-26 21:20Z
  (found by 480):** the first version of this bullet attached SQ W's 09-07 description ("public K562 supervised response
  means, scale 1… Gamma-Poisson decoder… mean-cellCPM inverse") to the pack's later numbers. That described entry
  (`SourceResponse-GammaPoisson-v2`, 09-07 13:20Z) scored nmae **0.9818**, pds 0.736, fid 0.470, worse than our 0.9380
  that day; H1 agrees (480 §2–3). The pack's later entries (nmae 0.878–0.89 from 09-18, fid 0.528, reach 0.29, jac 0.033;
  ours 0.896 / 0.515 / 0.279 / 0.041) carry **no description**, so their method is unknown. There is no board evidence
  that a plain full-scale donor decoder beats our DE channel.
- **Testh's 09-23 jump** was DE-only (nmae 0.925 → 0.877, reach 0.226 → 0.270, pds/mse flat), undisclosed. Its shape
  resembles our 433h (zero repair + MDE).
- Already in our pipeline: rank-invisible magnitude moves (black_handsome_egg's "q-space amplitude transport" = our
  `VCC_DE_MAG`, the 09-28 zs riders); cis-neighbour boosting ("1 Undergrad", `cis:3:2.5`; 425); source weighting,
  whitening, shrinkage, decoupled rendering (Zemin Li's "per-cell total cap 900k").
- Le0 submits under other teams' model names (Coral Sea, nailong, sunny_cpro_v19f, hh-v32, …): misdirection, not
  information.

## 6. Reproduce

The aggregation script is in the reviewer session's scratchpad; the logic is: load every `data/leaderboard/*.json`,
take `entries`, dedupe by `id`, keep `model_name`, `description`, `organization` and the six raw/scaled members.
