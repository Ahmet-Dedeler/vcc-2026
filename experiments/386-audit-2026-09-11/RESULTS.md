# 386 — Independent re-derivation of the PI's 2026-09-11 numbers: all rows MATCH except two stale gap members (mse/pds vs 319p, not 319p2) and a mislabeled −216

Status: COMPLETE

## Known-answer reproduction (numbers, match)
- 323 cap15 legs ratio from on-disk `metrics-cap15.json`, my code: ΣN/ΣD = **0.9749462**;
  recorded `raw.mse` = **0.9749462**. Exact match (7 dp).
- 319p2 `score_avg` = mean of its six `score_*` from `data/board-entries/lqWNXX5YaYArdglu6Stn.json`:
  diff **−2.8e−17**. 319p likewise **0.0**. `score_avg = mean of six` holds exactly.
- PI `mse_alloc.py` copied verbatim (md5-identical) and run: reproduces every PI intermediate
  (alignment 18,077; top-10 share 0.58; D-weighted cos 0.103; global α\* = 0.0914;
  `(S_or−base)/ΣD = −216.0358`). My independent implementation agrees on all of them.

## What was measured (inputs, estimator, controls)
- Own scripts only: `audit_board.py` (board JSONs), `mse_alloc_mine.py` (independent CSR load,
  own pseudobulk/alignment/oracle algebra), `callset.py` (closed-form fid/jac algebra),
  plus an H1 `reference_de.parquet` n_real distribution check. No PI code reused except the
  verbatim `mse_alloc_pi.py` copy used as the control.
- Control for row 1: latest-by-mtime snapshot verified as `2026-09-11T1910Z.json`
  (mtime 1789153845, newest of all `data/leaderboard/*.json`); snapshot semantics checked
  (893 entries = 893 distinct `team_id`s: one row per team).
- Control for row 3: exact per-target-G conversion computed alongside the PI-style single-G
  conversion, so the units question is answered by comparison, not by assertion.

## Results

Main table (row, PI value, my value, verdict):

| # | row | PI | mine | verdict |
|---|---|---|---|---|
| 1a | latest snapshot / count | 2026-09-11T1910Z, 893 | same file by mtime, 893 entries, 893 team_ids | MATCH |
| 1b | Illumina AI iCell avg / reach raw | 0.2872 / 0.367 | 0.2872087 / 0.3670132 | MATCH |
| 1c | cqawesome avg | 0.2506 | 0.2505626 (`m0911a`) | MATCH |
| 1d | our row 0.2201159, rank 15/893 | live board | JSON exact 0.2201159; exactly 14 snapshot rows above ⇒ implied rank 15 | MATCH (consistency; rank itself is live-only, see Limits) |
| 1e | gap/6: nmae / reach / mse / fid / pds / jac | +0.0224/+0.0170/+0.0156/+0.0079/+0.0067/−0.0023 | vs 319p2: +0.0224/+0.0169/**+0.0152**/+0.0079/**+0.0069**/−0.0023; vs 319p: all six exact | mse, pds MISMATCH vs 319p2 (stale: computed vs 319p); rest MATCH |
| 1f | Illumina 09-09→09-11 step | mse +0.151, rest within ±0.011 | mse +0.150836; pds +0.010680, nmae −0.001819, fid −0.003406, reach −0.008261, jac −0.003031 | MATCH |
| 1g | builders ≡ Illumina AI, one team_id | 0b01… | `0b01OX4HQx2BB1tRiny3` on both; names builders→Illumina AI, models PromoterAI_precomputed→iCell | MATCH |
| 2a | 319p2−319p Δ (bin 2) | +0.000297 | +0.00029691, inside [−0.0005, +0.0005] | MATCH |
| 2b | member deltas pds / mse | −0.00099 / +0.00204 | −0.000994 / +0.002037 | MATCH |
| 3a | oracle signal Δ | −216 (gene-sum) | gene-sum dS = **−33.1**; −216.04 = dS/ΣD | LABEL MISMATCH (number −216 reproduced; units mislabeled) |
| 3b | ΣD·G ≈ 18,077·ΣD | 18,077 | aligned genes exactly 18,077/18,533; ΣD = 0.153073; h1 universe 18,080 | MATCH |
| 3c | ≈ −0.012 raw mse | −0.012 | **−0.011950** exact per-target-G; −0.011951 PI-style | MATCH |
| 3d | top-10 share 58% / Spearman 0.19 / D-cos 0.103 | as stated | 0.580 / 0.1873 / 0.1035 | MATCH |
| 4a | I for n_real ∈ {100, 288, 765} (J=0.0409, n_pred=2,599) | — | 106.05 (**infeasible**: exceeds n_real; min feasible n_real = 106.3) / 113.44, recall 0.394 / 132.18, recall 0.173 | computed; n_real=100 ruled out by the definitions alone |
| 4b | implied sign accuracy a on intersection (fid 0.508) | — | **not identified**: uniform-accuracy ⇒ 0.508 (= fid) for every n_real ≤ n_pred; coin-flip-outside ⇒ 0.68/0.66. k = 1,320 vs I ≈ 113: ≤10% of fid credit can lie in the intersection | computed; no single a implied |
| 4c | fid at exactly 400 calls, n_real 288 vs 1,000 | — | ceiling min(1, 400/n_real): **1.0** (no yield penalty) vs **0.40** (hard under-call cap) | computed |
| 4d | "over-calls 30×, fid/jac unclaimed 2/6" | — | 30× is TRUE vs the median (H1: median n_real = 90 ⇒ 28.9×; mean 766 ⇒ 3.4×). "Unclaimed 2/6" NOT SUPPORTED by §4/§6: fid has no over-call penalty term (yield cap binds only when K < n_real; marginal extra call at a≈F is neutral), and jac break-even at K=400/n_real=288 needs 6.76% top-400 precision vs current 4.36% (1.55× enrichment) just to tie 0.0409 | slogan half-true; prize claim rejected |
| 5 | "flat ⇒ no third dose" per prereg | bin 2 | Δ=+0.000297 ∈ bin 2 ⇒ "saturated; the family closes at 0.34–0.5 with the incumbent staying 319p (ship the lower dose)"; top bin grants a 0.75 point only "if mse still carries it" (conditioned, not automatic) | MATCH with nuance |

Quote for row 5 (`experiments/319-cap-variant/PREREGISTRATION-partial-dose2.md`, bins table):
"| **−0.0005 … +0.0005** | saturated; the family closes at 0.34–0.5 with the incumbent staying 319p
(ship the lower dose). |" and "| **Δ ≥ +0.0005** | under-dosed at 0.34; 319p2 becomes the incumbent;
a third point (0.75) earns a slot only if mse still carries it. |"

## Verdict against the PLAN.md kill criterion
No kill triggered: the audit is the deliverable and every row is adjudicated. Nothing overturns
a PI decision — the dose family stays closed (our Δ lands bin 2 cleanly) and the −0.012 ceiling
stands (my exact number agrees to 6 dp). Two findings need PI attention but change no verdict:
(1) the gap decomposition is computed against 319p, stale for mse/pds now that 319p2 is best
(+0.0152/+0.0069, not +0.0156/+0.0067); (2) the "−216 (sum over genes)" label is wrong —
the gene-sum is −33.1 and −216 already contains a /ΣD; the −0.012 survives only via the extra /G.

## Limits
- H1-panel numbers are development evidence, not A/B/C forecasts; the H1 n_real distribution
  (mean 766, median 90, deciles 1/3/8/90/1020/2692/5485) grounds only the heavy-tail premise.
- Neither 319p nor 319p2 appears in ANY on-disk leaderboard snapshot (snapshots are one-row-per-team
  and lag the live board; our on-disk team row is the older fern entry at 0.21854). Rank 15 is
  implied (14 rows above 0.2201159), not directly re-derived.
- "Best public-data team" for cqawesome is not verifiable from the JSON schema (no public/private flag).
- n_pred = 2,599 is task-given; its provenance (which entry/panel mean) was not re-derived.
- The "lead calls it customized PromoterAI_precomputed" quote is external; on disk only the model-name
  string `PromoterAI_precomputed` (09-09) → `iCell` (09-11) is confirmed.

## Anything the PI did not ask about
The verbatim PI script's own printed lines read "delta raw mse −216.0358" and "overall /6 ->
+37.7026": the first omits the /G the task text applies afterwards (×18,077 off as a raw-mse
statement), and the second additionally flips the sign, printing a +37.7 overall "gain" for what
is a −0.012 raw improvement. Anyone re-running `mse_alloc.py` naively will quote −216/+37.7, not
−0.012/+0.0021. Related: the script's proxy deltas (−76.9 norm-match, −112.5 sqrt(D)/|d|) are in
the same /ΣD units, i.e. −0.0043/−0.0062 raw after /G — consistent with the "proxies ≤ −0.006"
gloss only for the best proxy, barely. Separately: k = 1,320 correct-direction calls vs I ≈ 113
means ≥90% of 319p2's fid credit comes from sign agreement on genes the reference did NOT call
significant — fid is mostly a background-direction meter, which is why shrinking K cannot move it
much and why the "2/6 unclaimed" framing misreads §4.

## Exact commands to reproduce
- `uv run python experiments/386-audit-2026-09-11/audit_board.py` (rows 1–2; needs repo root CWD)
- `uv run python experiments/386-audit-2026-09-11/mse_alloc_pi.py` (verbatim PI control; repo root CWD for its relative paths)
- `.venv/bin/python experiments/386-audit-2026-09-11/mse_alloc_mine.py` (row 3, independent)
- `uv run python experiments/386-audit-2026-09-11/callset.py` (row 4)
- H1 n_real distribution: parquet read of `experiments/323-current-family-validation/work/reference_de.parquet`
  (p_adj < 0.05 per-target nunique; see command in session log)

## Peak memory and wall time
- `mse_alloc_mine.py`: peak RSS **3,992,977,408 B (≈3.7 GiB**, `/usr/bin/time -l` maximum resident set
  size), wall **6 s**. Under the 8 GiB cap; memory gate at start 86% free, no other heavy job run.
- Board/call-set/H1-distribution checks: seconds each, negligible memory. Total session well inside
  the 3 h budget. No scratch files retained (only the five audit scripts + PLAN + RESULTS in
  `experiments/386-audit-2026-09-11/`); nothing written outside it except the md5-verified PI copy
  inside it.
