# 342 — current-family Hamilton remeasurement does not advance

September 9, 2026. **NO_LOCAL_ADVANCEMENT.**
Frozen at **5cb7dcfc** before emission. No retuning, additional arm,
package, production change or submission.

Hamilton is an existing rider from 196/208/230, not a newly discovered
method. Its older positive price does not transfer unchanged to this
current cap15 composition and exact 323 evaluator.

## All six actual-count results

Same 85 H1 targets, 323's H1-excluding operator, H1 union OFF, cap15
effects and rc_callset. Each changed arm restores the same integer mass:
7,761–7,933 added counts per target, median 7,868, one count per selected
gene. Hamilton selects largest fractional remainders; the control selects
uniformly from positive fractions using its own generator.

| arm | raw PDS | raw MSE | raw NMAE | fidelity | reach | Jaccard | H1-local overall |
|---|---:|---:|---:|---:|---:|---:|---:|
| exact floor | .88963585 | .97494624 | .95083646 | .49755452 | .25434726 | .07475983 | .28095315 |
| Hamilton | .88333333 | .97064645 | .95025181 | .49724789 | .25367803 | .07476164 | .27856874 |
| random up | .88375350 | .97338276 | .95001986 | .49729918 | .25879935 | .07476265 | .27979416 |

Lower raw MSE/NMAE is better. Overall uses H1-local anchors, not the
leaderboard's A/B/C anchors.

Paired 10,000-target bootstrap, seed 20260909:

| contrast | local overall change | 95% interval |
|---|---:|---|
| Hamilton minus floor | **-.00238441** | [-.00472537,+.00056315] |
| Hamilton minus random up | **-.00122542** | [-.00405121,+.00127225] |
| random up minus floor | -.00115898 | [-.00399794,+.00244423] |

Hamilton versus floor:

- Raw MSE improves **-.00429979**, interval [-.00572617,-.00336141].
- Raw NMAE improves **-.00058465**, [-.00116248,-.00006173].
- Raw PDS falls **-.00630252**, [-.01302521,+.00028011], unresolved.
- Fidelity falls **-.00030663**, [-.00058781,-.00002971], resolved.
- Reach -.00066924 remains unresolved; Jaccard +.00000180 has a
  lower bound exactly zero.

Both overall lower-bound gates and the +.001 point floor fail; the
fidelity-harm guard also fails. Better agreement with the intended
fractional profile is not sufficient for a competition gain. The
historical forecast of a roughly -.006 raw NMAE improvement is not
reproduced on this composition/evaluator.

No metric was assumed invariant. Both changed arms alter all 85 panel
total vectors. The report recomputes the MSE ratio and across-target
jackknife budget in every resample; although point rho is 1, bootstrap
rho reaches below 1 in all arms.

## Construction and independent checks

Claude wrote the initial experiment-local implementation without
command execution. Codex reviewed every module, passed its eleven
initial tests, and corrected one pre-freeze integrity check: floor
must retain its intended truncation deficit, not restore rounded mass.
A twelfth regression test covers that distinction. The scientific
gates and rounding rules did not change.

The adapter changes exactly one assignment in production `emit_counts`
inside the quasi branch. It compiles a separate function object; no
NumPy or production-module monkeypatch. The other floor-allocation
sites remain unchanged. New totals pass through pass-1 statistics,
rank design, placement and per-cell welding, never a post-emission
count addition.

- **Exact baseline:** all 85 floor blocks reproduce 323's cap15 whole
  count digest, panel digest, gene sums, jackknife/DE arrays, per-target
  rows, MSE legs and all six metrics. Differences are **zero**; even
  the saved statistics archive SHA matches.
- All 255 emitted and scored total vectors satisfy their rule.
  Both changed arms move each coordinate by 0 or 1, conserve the same
  rounded mass and select only positive fractional remainders.
- The reviewed scorer replay reproduces all six metrics and target
  rows. An alternate dense gene-wise MSE budget reproduces all three
  contrasts and member intervals within **4.39e-16**.
- A separate coordinator audit, importing none of 342's rule/seam/
  report/audit modules, reconstructs intended counts from frozen
  profiles and library picks. An explicit lexicographic Hamilton
  ordering and separate random selection reproduce **all 255 total
  vectors exactly**, their panel digests and the frozen verdict.
- Both changed-arm full count matrices were streamed, not retained;
  their complete emission was not independently rerun. Stored
  sufficient statistics, totals and digests are retained.

These are implementation/arithmetic checks, not independent biological
replication. This is a reused H1 development panel. Bootstrap intervals
condition on one fixed emission/control seed and fixed reference
gallery/anchors; they do not establish A/B/C transfer or cross-seed
stability. No submission is cleared.

## Execution and receipts

Run **19:01:37.130–19:08:54.016 UTC**, 431.22 seconds after setup.
Sampled process-tree peak **4.80 GiB**; driver high-water RSS
5,217,615,872 bytes. Scorer/dense-bootstrap audit peak 3.72 GiB.
Coordinator allocation audit finished **19:15:18.201 UTC**.
All processes exited successfully; production files remain unchanged.

- Inputs SHA:
  `c40fc9be90839e593b14e3de3d905c45ba7ae1c68ddc859c9c8e2d20fbe2e5f7`
- Result SHA:
  `6c4f1cc9d280b6b5230490c488416f8fe0c1aed680cab83e8c31d7e36feea37c`
- Scorer/dense-bootstrap audit SHA:
  `e2a3eb7ee88540883bc325685749c0b25c04092b13b60b2292c4c64ad586fedb`
- Coordinator allocation audit SHA:
  `f36ba9ac0cbafe897cfc9fa2eb70793fe0a12ed9eb25f4624a29ac06e0187e74`
