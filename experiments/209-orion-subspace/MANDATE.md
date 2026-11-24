# 209 — The Orion pair shares an excess SUBSPACE, not a rank-1 axis. Equal weighting double-counts a platform.

You are a research subagent. Work ONLY in `experiments/209-orion-subspace/`.
No `vcc`, no `cell-eval2`, no submission build, no writes outside your directory.
Peak RSS ≤ 12 GB. `uv run python`.

## READ FIRST
`AGENTS.md` §4.2 **in full**, **§8.2** (direction is the entire remaining gap), the standing
HARNESS WARNING on the long-gene technical axis, §11.1 rules 1–13. Then:
1. `experiments/204-repression-weighted-mixing/RESULTS.md` — **the measurement this mandate
   executes, and its explicit hand-off.** Verify it; do not inherit it.
2. `docs/LOG.md` 2026-09-02, the entry *"Orion pair shares a protocol axis (cos 0.872, HEK293T
   19.7% variance); per-donor rank-1 deflation lifts consensus sign 0.5697→0.5770; handed to 178"*,
   and `experiments/178-composed-riders/RESULTS.md`, which closed the rider family.
3. `experiments/159-transfer-ceiling/RESULTS.md` — the `1/k` law and κ².
4. `experiments/197-cosine-inversion/RESULTS.md` — the exchange rate, **+0.00365 overall per +0.01
   of cosine**, and the field's implied `k_eff`.
5. `experiments/206-cluster-method/RESULTS.md` §5 — why nobody is on the `1/k` curve.

## The argument

Production averages three donors **with equal weights**, and `AGENTS.md` §4.2 lists that as one of
its four load-bearing findings. 204 has now closed re-weighting for good: equal weights are the
**maximum of the transfer cosine over all 29 alternatives**, the permutation falsifier **ties** the
best real arm, and the negative has power (the minimum-variance model's own inputs predicted +8.5%
to +16.8% cosine; measured **−1.50%, CI [−2.38%, −0.67%]**).

**But every one of those 29 arms was a SCALAR per donor.** 204's own mechanism section says why that
cannot work, and hands over what it found instead:

> a one-factor test on 4–6 lines **rejects the shared-component model wherever it has power**
> (RMS residual 19.7–22.4%, |z| up to 11), its same-cell-line positive control fires at **z = +11.2**,
> and **the second-largest excess is the Orion↔Orion pair (z = +5.1 to +7.6)** — the product the
> three-donor triple has zero degrees of freedom to avoid.

and, as its closing hand-off: **"the Orion↔Orion excess is a *subspace*, not the blob."**

**So the three donors are not three independent observations.** Two of them are the same platform,
same chemistry, same depth, run by the same lab, and they share structure *beyond* the shared
biological signal and *beyond* the long-gene technical axis. Equal weighting therefore counts the
Orion platform twice and K562 once. **A scalar weight cannot fix that** — the excess lives in a
subspace, so removing it requires a projection, not a coefficient. That is precisely why 204's
29-arm scalar sweep found equal weighting optimal and still left an effect on the table.

**Why this is not the closed rank-1 result.** On 2026-09-02 a **rank-1** deflation of the shared
Orion axis (cos 0.872, 19.7% of HEK293T's between-perturbation variance) was measured to lift
consensus sign accuracy **0.5697 → 0.5770** and was handed to 178, which closed the rider family on
composition grounds — *"gains compose sub-additively, costs do not"* — not on this arm's own merits.
204's finding is that the excess is **not rank-1**. A rank-1 deflation therefore removed part of it
and a subspace deflation has never been tried at all.

**And it is the right shape for the problem.** 206 §5 measured that the leaders' implied `k_eff` of
10–14 **exceeds the ~11 CRISPR-perturbed human cell lines that exist in total**, so nobody is on
159's `1/k` curve. Being "off the curve" means extracting more independent information per line than
naive averaging does — and de-correlating two donors that are secretly one platform is exactly that
operation. Our `k_eff` is 3.09 with three lines; if the Orion pair is partly one line, our true `k`
is nearer 2, and removing the shared excess should move it toward 3.

## Work

### Phase 0 — known-answer rows (mandatory, abort on failure)
- Reproduce 204's Orion↔Orion excess: the **z = +5.1 to +7.6** figure and the same-cell-line positive
  control at **z = +11.2**. If your one-factor test does not fire on the positive control, it has no
  power and nothing below is readable.
- Reproduce the 2026-09-02 rank-1 numbers: shared-axis **cos 0.872**, **19.7%** of HEK293T's
  between-perturbation variance, and consensus sign **0.5697 → 0.5770** under rank-1 deflation.
- Reproduce PLOLO's a15 `pds` **0.795784** and a23 **0.789121**, and the a15→a23 mass-fix price
  **+0.004014** scaled (204 reports all three exact).
- **Rule 11:** path and sha256 for every operator. Shipped: `29,394 × 17,942`, sha **`e2f65ecf`**
  (193's correction to §4.2's stated dimensions, which match no file on disk).

### Phase 1 — characterise the excess subspace
Estimate the subspace shared by `orion_hct116` and `orion_hek293t` **beyond** what each shares with
`k562_gwps`. Report its dimension (by a criterion you state in advance, not chosen on the endpoint),
its energy share in each donor, and — mandatory — **its overlap with the long-gene / nuclear-fraction
technical axis**, which the standing harness warning documents as the Orion pair's dominant shared
direction (top-PC cosine 0.872, and 36.3% of HEK293T's *harness* effect energy). 204 reports blob
deflation moves its own numbers by <1%, so the excess is claimed to be distinct — **verify that
independently, because if the subspace is just the blob at higher rank, production's whitener
already removes it (0.11% of shipped energy) and the arm is void.**

### Phase 2 — deflate it, and price the direction
Remove the excess subspace from the two Orion donors (per donor, symmetric), re-form the equal-weight
consensus, and measure:
- **Δcos on PLOLO — the primary endpoint.** 204 measured this instrument as resolvable to **±0.20%
  relative ≈ ±0.00006 overall**, which is 25–200× below any effect worth shipping.
- PLOLO `pds` (162's four folds, challenge frame). **`debench.Surface.pds()` must never size a `pds`
  arm** (188 Phase 3).
- **|t|-weighted `q_w`** (187 §5) and pooled consensus sign accuracy, against the rank-1 arm's
  0.5697 → 0.5770 so the two are directly comparable.
- The implied **`k_eff`** by 159's route, before and after. That number is the mechanism claim.

Sweep the deflation strength from 0 (identity, your known-answer row) to full, and report the
maximiser of cosine and of `pds` **separately** — 062's dose was chosen on `pds` and 197 showed that
is the member least sensitive to direction, which is the mistake this project has already made once.

### Phase 3 — six members and the gate
Rule 12. `mse` through 096's `claim` route — a donor-side change moves the effect's shape, so `mse`
is **not** exempt (a19 lost −0.01352 exactly here) and we ship `claim` at **100.0% of the cap** with
only 5.6% margin, so report what the arm does to saturation. `nmae`/`reach`/`fid`/`jac` on h1-test
with a **significance-selected pool** (188). Donors on h1-test: `k562_gwps` + `orion_hek293t` only
— **note that this surface holds out one of the two lines whose correlation you are removing, so it
cannot see the arm.** Say so rather than reporting a silent zero; 200 caught exactly that failure in
the CD4 residual today.

Rule 13's `reach` toll does not apply to a profile-class arm. Gate **+0.0019** shape-only,
**+0.0089** if the arm changes amplitude — re-solve `VCC_ALPHA` and say what that cost.

### Phase 4 — falsifiers
1. **Permutation control**: deflate a random subspace of the same dimension and re-score. 204's own
   repression falsifier *tied* its best real arm, so this is the outcome to watch for.
2. **Blob control**: the arm's gain must survive deflating the long-gene axis first. If it does not,
   the arm is the blob and is void.
3. **Leak-free fold**: report all four PLOLO folds including `kolf2`, which 162 flags as the only
   leak-free one, and never quote an absolute from the Orion pair.
4. **The K562 sanity check**: deflating an Orion↔Orion subspace must not change `k562_gwps`'s own
   rows at all. Assert it.

## Deliverable
`RESULTS.md`: known-answer rows, the subspace's dimension/energy/blob-overlap, the deflation sweep
with **cosine primary**, the `k_eff` before and after, the six-member price, all four falsifiers,
and a verdict in overall units. **"The excess is the blob at higher rank" or "deflation ties its
permutation control" are full results** and would close the last structural idea about our donor set.

## Territory
207 owns the CD4 residual. 208 owns emission riders and rounding. 199 owns D/E/F readiness.
**You own the donor-donor correlation structure.** Do not re-weight donors by a scalar (204 closed
it), do not tune the operator (193), and do not re-estimate the donor tables (194).
