# September 9 research sprint — final handoff

Research cutoff: **20:58 UTC / 23:58 Istanbul, September 9, 2026**.
Final read-only live refresh finished **21:10:44 UTC, September 9 /
00:10:44 Istanbul, September 10**; closure operations followed the cutoff.
No research or new worker was launched after it. Ahmet alone decides submissions.

## Current decision

Through experiment 342, **no new candidate has cleared its frozen gates**.
The best recorded A/B/C result remains 319 at **0.21865973944346428**.
The authenticated/public check at **21:10 UTC, September 9** shows **rank 13 of 839**, unchanged
score. Rank was 6 at 08:19 UTC; no model or submission changed.

Two Claude Code Fable 5.1 reviews completed: general strategy and emission
construction. Both were launched with xhigh effort on the existing Max
login, with only Read/Glob/Grep, no MCP servers and `apiKeySource=none`.
Their observed rate records say `isUsingOverage=false`. A third, scoped
Fable worker implemented 342 without command execution. Codex reviewed
and corrected that code, passed twelve pure tests and froze it at
**5cb7dcfc**; numerical execution/verification remains Codex's responsibility.
Those three workers are terminal. A fourth restricted read-only critique
also completed; Codex adjudicated it in
`docs/sprint-2026-09-09/FINAL-REVIEW-DECISIONS.md`. Only the corrected 341 label-regime pilot
and 342 rounding remeasurement were selected from their recommendations.

No submission, reset-time queue, production-default promotion, push,
paid API/cloud use or new reserve access has been performed by this sprint.
This handoff is not permission for any of those actions.

## Experiment ledger

Each experiment's directory contains its authoritative registration,
input manifest, result and verification record. The table distinguishes
source/engineering evidence from actual-count H1 results.

| experiment | hypothesis and key controls | measured outcome | terminal status |
|---|---|---|---|
| 325 source provenance | Accession repeats might provide independent residual labels; official GEO/ENCODE sample provenance | 19 records inspected; no pair verified as fully independent same-state repeats; KHSRP reuses ENCODE material | Metadata complete; no label pilot |
| 326 cap30 diagnostic | Existing cap30 versus cap15, with plain-dose interaction | H1-local overall +.0001924 [-.0046193,+.0031757]; raw MSE worsens | Descriptive, not candidate clearance |
| 327 coordinate median | Robust donor pooling; magnitude-only, power and shuffled controls | H1-local overall +.0008450 [-.0020704,+.0032102]; controls not beaten, raw MSE cost | No advancement |
| 328 joint CD4 estimation | Joint source coefficient posterior; diagonal and target-shuffled covariance | After a disclosed post-result numerical amendment, source log-density +.01337746 over diagonal, CI [+.01034952,+.01685980] | Source gate passed only; tested in 329 |
| 329 joint count stage | 328's fixed posterior on current cap15; all six members and fitted controls | H1-local overall -.0056597 [-.0169085,+.0028336]; resolved NMAE/fidelity harm | No advancement |
| 330 background gene scales | Non-panel source scales; shuffled/scalar controls | All six raw point estimates favorable; overall +.00083497 [-.00030217,+.00216558] | Unresolved; no advancement |
| 331 TabICL runtime | Local synthetic regression; mean and shuffled-label controls | Repeatable CPU inference near 1 GiB peak RSS | Engineering feasible only |
| 332 conditional median | TabICL response median; mean, recipient-free and wrong-recipient controls | Absolute error improves; retrieval .675662 versus source mean .732673; recipient advantage unresolved | No advancement |
| 333 GenePT features | Released vector coverage and exact identity checks | Coverage passes, but 34 challenge targets participate in exact-vector collisions across distinct descriptions | Vectors quarantined |
| 334 direct text | TF-IDF/SVD text plus source responses; diagonal, source-only, STRING and shuffles | Retrieval -.0342726 versus source mean; -.0131545 versus source-only; shuffled text not beaten | No advancement |
| 335 biomedical encoder | Local fixed-revision encoding; identity, token and repeatability checks | 1,528 available rows; duplicate vectors all have identical consumed tokens | Engineering feasible only |
| 336 semantic transfer | Same 334 design with biomedical features; matched prior and shuffle controls | Retrieval -.0221851 versus source-only, -.0059327 versus shuffled semantic features, both resolved | No advancement |
| 337 availability consensus | Normalize zero-filled donor support; count-only and shuffled controls | Overall +.00322188 [-.00555539,+.00643962]; raw MSE improves, NMAE/fidelity/reach worsen | No advancement |
| 338 rank/moment feasibility | Fixed-U/zero/column-total histogram LP; same-depth and no-rank controls | 128/128 constrained pairs solve, but useful symmetric room fails in 3/4 targets; 768 certificates verified | No integer construction earned |
| 339 marginal calibration | Gene-wise empirical marginal mapping; diagonal, Gaussian, clipping and shuffled curves | Source retrieval +.00293255 [crossed +.00013215,+.00580927], below floor; diagonal/Gaussian not beaten, cosine cost | No advancement |
| 340 neighborhood transfer | Sparse response-gene coupling; diagonal, own-only, random/correlation graphs and label shuffle | Retrieval -.00520658 versus diagonal [crossed -.00748176,-.00282211], negative on all six lines despite cosine gains | No advancement |
| 341 known-line calibration | Own-line training labels; matched cross-line full/diagonal models and whole-target label shuffle | Weighted sign +.03254973, but retrieval -.04829784 versus diagonal; every sign gate passes, both retrieval guards fail; independently audited | No advancement |
| 342 Hamilton rounding | Existing rider remeasured on current cap15; exact floor and mass-matched random-up control | Overall -.00238441 [-.00472537,+.00056315], no advantage over random, resolved fidelity harm; exact baseline and both audits pass | No advancement |

Important scope: 326/327/329/330/337/342 use actual emitted counts and the
current-family, H1-union-OFF 323 evaluator. Their “overall” is on H1-local
anchors, **not an A/B/C score forecast**. The raw MSE cost is retained even
when a local point score is clamped to zero.

The six-line response pilots use selected coordinates and unit-RMS source
responses; their squared/absolute errors are **not competition MSE/NMAE**.
Their paired target and crossed line/target intervals are conditional on
the fitted models. Reused development data, correlated targets, study/line
confounding and the selected gene axis limit transfer. No new independent
biological confirmation is claimed.

### Source-stage uncertainty, without treating it as a leaderboard forecast

332 uses **12 Broad Cas9 lines and 82 targets**, with 512 source-selected
coordinates per fold. The shared frame in 334, 336, 339, 340 and 341 is
**six CRISPRi lines, 1,024 targets and 512 coordinates**. Their retrieval
numbers are not directly comparable across those two frames.

| fixed primary contrast | point | target 95% interval | crossed line/target 95% interval |
|---|---:|---|---|
| 332 absolute-error improvement vs no-recipient | +.00000416 | [-.00026110,+.00028280] | [-.00064268,+.00062722] |
| 332 retrieval vs source mean | -.05701094 | [-.07620697,-.03768851] | [-.08854073,-.02834143] |
| 334 retrieval vs source-only | -.01315448 | [-.01525956,-.01110246] | [-.01893270,-.00805888] |
| 334 retrieval vs shuffled text | -.00057420 | [-.00271954,+.00154441] | [-.00441454,+.00346421] |
| 336 retrieval vs source-only | -.02218506 | [-.02477252,-.01970046] | [-.03072895,-.01439978] |
| 336 retrieval vs shuffled semantic | -.00593272 | [-.00810338,-.00382780] | [-.01018319,-.00157457] |
| 339 retrieval vs source mean | +.00293255 | [+.00193574,+.00390466] | [+.00013215,+.00580927] |
| 340 retrieval vs diagonal | -.00520658 | [-.00649882,-.00395619] | [-.00748176,-.00282211] |
| 341 weighted sign vs cross-line diagonal | +.03254973 | [+.02995022,+.03519118] | [+.01553181,+.05421251] |
| 341 retrieval vs cross-line diagonal | -.04829784 | [-.05445425,-.04226641] | [-.07344306,-.02750191] |

These are the saved, previously audited paired intervals, not new
resampling or selection. 339's positive source-mean contrast misses its
.005 point gate and does not beat the stronger calibration controls.
341 deliberately allows the recipient's **training-target** labels;
it is known-line development, not an unseen-context validation. Its
sign improvement does not override the failed retrieval safeguards.

## Technical readiness versus scientific readiness

**319:** retained best-scored package. Its bytes were rehashed in this
continuation and match the saved pin:

`submissions/319-normrel15-cap.vcc`

`57585d558f9777f066c02f5b299cacccca0dd4954b4201b8fa653bd0f6816dd2`

**319d3:** pre-existing, unsubmitted diagnostic package, not a sprint-built
improvement. Its bytes were also rehashed and match:

`submissions/319d3-cap30.vcc`

`6ddbf76d2c7b914fae09f606f5da44262f25e415870552f211d606fa08c33d74`

The intended distinction is cap-mode exponent -3.0 versus -1.5, not a
different biological source or an isolated emitted-amplitude contrast.
Fresh package/operator/sidecar comparison passes with explicit provenance
qualifications in `docs/sprint-2026-09-09/package-verification-2005.json`:

- Configs differ in exactly `normrel_exp`.
- Three CD4 metadata hashes differ solely because `peak_rss_gib` differs;
  consumed response/basal/axis files match and were rehashed.
- The driver temporarily embedded those pins into each builder. Both
  historical builder hashes are reconstructed exactly in memory; all
  code outside the two pin registries is identical. Current source
  matches both preserved pre-pin copies.
- `frozen_inputs.json` differs only in its `repinned_at` timestamp.

No package, source, sidecar or original build receipt was changed.
Neither package was rebuilt in this continuation. A final live operational
refresh is recorded in `docs/sprint-2026-09-09/live-verification-closure.json`.

Earlier September 9 storage verification recovered and checked the
prediction-content digests. Original deleted `.h5ad` byte identities are
not promised by restoration. See `docs/STORAGE-AUDIT-2026-09-09.md`.

No new artifact is scientifically cleared merely because its code tests,
hash receipts or local scorer pass. Experiment 323 reproduced known
board point signs but did not resolve the small contrasts; it is not a
standalone selector for small improvements.

The complete retained-record index is
`docs/sprint-2026-09-09/ledger-verification-2005.json`: **639 declared input
hash checks, 346 unique files, no mismatches**. One reserve-labelled path
was deliberately not opened. This verifies retained artifacts and flat
hash tables, not every nested manifest field or a new model rerun.
328's original failure and 334's original text-coverage failure are
explicitly indexed alongside their recoveries.

## Final operational state and interpretation

`docs/sprint-2026-09-09/live-verification-closure.json` records the authenticated
entry and public board. `can_submit=true`, `limit_reached=false`; the
endpoint does not return a numeric remaining allowance. The known entry
is terminal, but that is not an exhaustive search for other in-flight
entries. No upload is authorized by those availability flags.

At that checkpoint, leader score .2630424866 versus our .2186597394
leaves a .0443827472 gap. Reading the public scaled columns directly:

| member | leader-minus-us contribution to overall |
|---|---:|
| PDS | +.00490029 |
| MSE | -.00885597 |
| NMAE | +.02282486 |
| fidelity | +.00860005 |
| reach | +.01868322 |
| Jaccard | -.00176970 |

These contributions sum to the observed gap; they do not estimate anyone's
cosine or a causal route to improvement. They reinforce a limitation of
source-only proxy gains: neither own cosine nor retrieval nor a weighted
sign proxy alone prices the six-member outcome.

Four Claude workers completed through the existing subscription: three
read-only reviews and one restricted experiment implementation. All
numerical execution and integration were coordinator-owned. The observed
streams show no API-key source, MCP servers, recursive workers or overage.
No new Claude login was needed.

## Submission recommendation

1. **Keep 319 as the best known result.** Do not resubmit it to restore rank.
2. Do not upload any failed/unresolved sprint construction as a claimed
   improvement.
3. 319d3 remains an optional **information** experiment if Ahmet wants the
   frozen cap-dose contrast, **not an improvement recommendation**.
   Its local overall is unresolved, and raw MSE worsens **+.0151952
   [+.0029296,+.0290382]**. The H1 point clamp hides that loss; 319's A/B/C
   MSE member is positive. That is a concrete downside risk, not a
   calibrated prediction of a board loss. The local cap-versus-plain
   interaction is already measured; an upload would test its transfer.
   The leaderboard shows the latest entry, so a lower result replaces
   the visible best row. Read the September 9 coordinator corrections
   alongside the original dose registration: no automatic first slot,
   no claim that magnitude was excluded, no pre-authorized higher dose.
4. 341 does not earn cross-study/count-stage advancement; 342 does not
   earn a new package or dose. The proposed
   empirical-null upweighting and block-sign repair were not selected.
5. No new candidate is ready for an improvement submission. The result is
   better-verified evidence and preserved artifacts, not a score gain.

The numerical jobs and four Claude workers are terminal. The closure receipt
`docs/sprint-2026-09-09/handoff-final-verification.json` records the final
package/sidecar checks, unchanged production/driver/submission ledger,
preserved snapshot PID 82906, and the paused research heartbeat with its
other fields preserved. No automatic submission or reset-time queue exists
from this sprint.
