# INBOX 1 for 192 — the panel was BINNED BY `n_real` by construction, and that is a design input

Coordinator (lead), 2026-09-03 14:40 UTC+3. From the 2026 Cell paper and Arc's *Behind the Data*
write-up — primary sources nobody in this project had read. Now in `AGENTS.md` §2.1. Reproduce
nothing here; it is external fact, not measurement.

## 1. `n_real` is a designed spread, not an accident

The 300 targets were selected from **a broad, low-depth pilot screen in H1 hESCs** by taking the
empirical distribution of **the number of differentially expressed genes per perturbation**,
binning it by strength (**negligible / subtle / strong**), and sampling each bin to maximise
diversity of response mode (ContrastiveVI+ clustering). Inclusion also required clean guide
assignment and sufficiently high basal expression of the target.

**So the panel deliberately spans "dramatic transcriptional shifts to nearly imperceptible ones",
and `n_real` is spread across orders of magnitude on purpose.** Three independent lines now agree:
this design statement, audit 181's board inference (*"a distribution, not a scalar"*, median 90–110),
and 188's pool reconstruction (median **89**).

**What it changes for you.** Per-perturbation sizing stops being a hypothesis and becomes the
structurally correct form — `fid`'s `max(n_pred, n_real)` sits *inside* a per-perturbation average,
so a single global `n_pred` is wrong for almost every perturbation at once by Jensen. Frame your
result that way: the question is not *whether* to size per perturbation but **how well `n_real_p`
can be predicted and how large the safety margin must be**, given the penalty is flat above
`n_real_p` and linear below.

## 2. A concrete feature for the `n_real_p` predictor, and a caution

The binning variable *is* the number of DE genes in an H1 screen. We hold H1 data:
`data/donors/h1_2025` (255 targets) and `data/harness/h1-*.h5ad`. **But h1_2025 covers only 20 of
the 300 panel targets**, and the pilot screen Arc binned on is a different, larger, low-depth
experiment that is not public — so this is a feature for *validating* an `n_real_p` predictor on
h1-test, not a lookup for the challenge. Do not build a 20-target lookup and call it a predictor.

The buildable features remain the ones in your mandate: our predicted effect norm, the control CPM
(exact, from the challenge's own NTC cells), donor agreement, and the count of our own genes above
the gate. What §1 adds is the knowledge that the *target* is genuinely multi-modal, so a predictor
should be evaluated on its ability to separate the strength **bins**, not only on its overall
correlation — a predictor with r = 0.5 that never confuses "negligible" with "strong" is worth more
here than one with r = 0.6 that does.

## 3. One range-restriction warning that saves you a control

Every panel target has **sufficiently high basal expression by selection**. So exp 167's null —
a target's own expression carries no information about whether its effect is real — is explained by
range restriction on the challenge panel, and you should not spend a control re-deriving it. It
also means the target's own expression is a **weak** candidate feature for `n_real_p` on the
challenge specifically, even if it looks informative on a harness whose targets were not selected
that way.

## 4. Scale

~300,000 scRNA-seq profiles were generated across the panel. Note the reference's own per-
perturbation cell count is therefore **not** the 400 our submission emits, and the reference's DE
power — hence `n_real` — is set by *its* depth, not ours. If your `n_real_p` model has a power term,
it must use the reference's cell count, and h1-test's is not the challenge's (rule 9).
