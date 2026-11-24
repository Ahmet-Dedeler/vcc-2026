# 330 — Panel-independent CD4 gene-scale estimation

September 9, 2026. Freeze code, protocol and input hashes before extracting
the broader source-response coordinates or constructing candidates.

## Hypothesis and distinctness

The current CD4 preparation estimates a separate mean-absolute scale for each
gene from only the source rows matching the requested panel: 78–79 rows for
H1, versus hundreds for A/B/C. A scale estimated on a fixed larger non-panel
source population may transfer better and be less query-panel-dependent.
This is an assumption to test, not a claim that broader normalization is
automatically more biological or better for retrieval.

This experiment does **not** use the 328 posterior or tune the 329 candidate.
All raw panel source coefficients remain unchanged. It changes only the
per-gene denominator's estimation population; it adds source variation
information absent from the compact panel-only preparation.

History checked: 062's original normalization; 125/203's posterior and
confidence models; 188's delivered magnitude maps; 194's donor re-estimation;
213's transport laws (CD4 not re-expressed); 246/249/250/251's centering and
projectors; 271/286's noise/cross-product estimators; 273's explicit panel-
dependent CD4 preparation; 322's KOLF-source breadth. None is this construction.
Do not call it a new statistical method or repeat the KOLF breadth claim.

329's mixed primary point row has been seen while its controls are running.
The present hypothesis comes from the source-normalization implementation,
not a selection among posteriors or a rescue change to 329. H1 remains reused
development throughout.

## Data and exact construction

Use the first **1,000 fit-target IDs already fixed in 328's metadata selection**,
all complete across Rest/Stim8hr/Stim48hr and excluding H1/A/B/C target symbols
and mapped IDs. Do not use the 250 source-development targets.
The 256 response genes of 328 have been inspected; the broader coordinates
are newly extracted from the already-used local DESeq2 source, not a reserve
or independent study.

Response axis: the existing compact gene map, restricted to official genes
and excluding **both** H1 and A/B/C target coordinates. Hash the exact mapping
and selected target rows. Read only source `log_fc`, convert exactly as 047:
float32 source coefficient times Python-float ln(2), nonfinite coefficients
zeroed. Require at least 99% finite entries per condition and no all-zero
response column; otherwise stop, without replacing targets or columns.

For a condition and background subset X:

1. subtract its own per-gene mean;
2. compute mean absolute centered effect per gene;
3. floor at the 25th percentile of positive gene means, as production does.

The resulting positive scale vector is used **only on the eligible response
coordinates** of the current compact panel banks. Other gene scales retain
their baseline values. Keep the panel banks' original centering, gene support,
row counts and raw coefficients. Recompute their median normalized row norm
with the new scale and the existing `EFFECT_SCALE` rule. Thus the normalization
rule retains its scalar calibration; only its gene-scale estimate changes.

Afterwards run identical mean3, complete-row freeze, k=5 day-of deflation,
f=.34 dose solving and global preservation. Main-donor cap15/PerturbMatch,
alpha .982, gamma .65 and the corrected 323 H1 whitener stay fixed.
H1 union OFF; mass fix unset; unchanged actual-count emitter and scorer.
The DESeq2 baseMean proxy is not used as recipient expression.

## Frozen arms

1. `cap15`: exact original preparation, known-answer row.
2. `background1000`: the only candidate, scale vector from all 1,000 targets.
3. `background80`: same rule on the first 80 targets of that fixed ordering.
   This tests sample-size sensitivity on the same non-panel population; it is
   not exactly matched to each bank's 78/79 rows.
4. `shuffled_scale`: background1000's scale values independently permuted
   across eligible genes in each condition (seed 20260909), then the same
   panel median-norm calibration. Preserves scale marginals, not biological
   gene identity.
5. `shuffled_targets`: complete target-row derangement of background1000's
   final prediction (seed 20260909); harness negative control, not candidate.

No scale blend, new quantile floor, sample-size grid, selected gene subset,
posterior estimates, new projector or dose adjustment.

## Checks and gates

Synthetic checks: centered scale's translation invariance, exact production
formula, floor handling and rejected degenerate data; unchanged ineligible
scales; median normalization under a fixed scale; identity override; shuffle
marginal preservation/determinism. Broader extraction must reproduce the
256 previously stored source gene coordinates wherever eligible, for all
1,000 targets/three conditions, exactly at the float32 coefficient level.

Before new scoring, cap15 must reproduce the prepared CD4 arrays, U5 hash,
final effect hash, count/panel digests, and all six metrics at unchanged 323
tolerances. Support counts stay 80 carried, 77 complete; incomplete rows and
H1 target response coordinates stay zero through deflation. Hash the operator,
raw banks, source selection, scale vectors, production imports and instrument.

Actual counts and all six members on 85 H1 targets. Primary contrasts:
background1000 minus cap15 and shuffled_scale; also report background80.
10,000 paired target bootstrap draws, seed 20260909, fixed gallery/models/
anchors; recompute panel-ratio MSE and global claim cap in every draw.
Report raw MSE even when local scaled MSE is clamped. Report source-scale,
prepared-residual, composed and emitted pseudobulk changes.

Earn further multi-context development only if both primary overall lower
95% bounds exceed zero, background80 does not beat background1000's overall
point estimate, and no primary profile-member interval is clearly harmful
(PDS upper >= 0; raw MSE lower <= 0; raw NMAE lower <= 0). Deranging targets
must reduce candidate raw PDS by >.10. All construction checks must pass.
Otherwise close this construction; do not retune it to pass.

These are conditional development intervals, not multiplicity-adjusted
discoveries or A/B/C submission clearance.

## Bound

One guarded heavy process, launch only after 329 exits. Extraction and all
scores within 90 minutes, ≤22 GiB tree RSS, ≤5 GiB scratch, repository memory
launch gates. No new source download, reserve, paid service, production file
change, package, submission, public disclosure or push.
