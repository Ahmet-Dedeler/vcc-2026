# 212 — Whitening costs 1.5–2.1× on cosine. Our whole remaining gap IS 1.66–1.76× on cosine. Re-sweep gamma with cosine primary.

You are a research subagent. Work ONLY in `experiments/212-gamma-cosine/`.
No `vcc`, no `cell-eval2`, no submission build, no writes outside your directory.
Peak RSS ≤ 12 GB. `uv run python`.

## READ FIRST
`AGENTS.md` §4.2 **in full**, **§8.2**, §11.1 rules 1–13 with the rule-13 table. Then:
1. `experiments/170-gamma-mse/RESULTS.md` **in full** — the closure this reopens, and **its
   instrument and known-answer gates are yours.**
2. `experiments/197-cosine-inversion/RESULTS.md` — the exchange rate and the field's cosine.
3. `experiments/195-retrieval-transform/RESULTS.md` — γ = 0.65 is the *retrieval* optimum by the
   direct objective. **Read why that is not the same claim as this mandate's.**
4. `experiments/208-free-riders/RESULTS.md` Phase 4 — the Ĉ/C_p by-product.
5. `experiments/178-composed-riders/RESULTS.md` and `docs/CURRENT_STATE.md`'s 178 entry — the
   ρ-cliff mechanism that closed the whole rider family.
6. `experiments/202-callset-rider/BOARD-RESULT.md`, `185-…/BOARD-RESULT-rc_hybrid.md` — the two
   board points that replaced the emission underneath 170's algebra.

## The argument, and it is arithmetic

`AGENTS.md` §4.2 records, measured: **"the whitening operator does cost 1.5–2.1× on directional
accuracy, monotonically on both surfaces."**

197 measured our transfer cosine at **0.079249** and the leaders' at **≥0.1312** (`cqawesome`
**0.1392**) — **1.66× to 1.76× ours.** 206 established they are on **public data with a method**
(three mutually independent pipelines above cosine 0.11, two of them lone individuals; implied
`k_eff` 10–14 against the ~11 CRISPR-perturbed human lines that exist *in total*, so nobody is on
159's `1/k` curve).

**The cost of our own whitening and the size of our entire remaining gap are the same number.**
§8.2: `mse` + `pds` + `nmae` are 135% of the gap and we already *lead* on `fid` and `jac`. At 197's
rate of **+0.00365 overall per +0.01 of cosine**, recovering even a third of the whitening cost is
worth more than every arm shipped today.

**Why 170 does not close this, and why 195 does not either.**
- **170 swept gamma on all six members and found every setting in [0.35, 0.80] net negative** — but
  its binding mechanism was the **ρ cliff**: *"lowering gamma raises the predicted effect's
  mean-square `A` and with it `B`; raising `B` past `claim` pushes all three contexts above it,
  where `min(claim, B)` saturates and the extra amplitude stops being refunded."* 178 closed the
  entire rider family on the same mechanism: *"every rider except sign repair raises `A` 11–28% and
  pushes ρ to 1.000 in all three contexts."* **That algebra was measured on a15's emission.** Two
  board submissions have since replaced it, and 208 measured `rc_callset`'s **Ĉ/C_p at 2.940 against
  a15's 1.125**. Whether that moves the ρ = 1 cliff, the `claim`-saturation cliff, or neither is
  **your Phase 1 and it is the crux** — do not assume it, and do not assume it does not.
- **195 showed the direct retrieval objective picks γ = 0.65 on its own.** That is a statement about
  **`pds`**, a rank statistic. 197 measured that we extract **26% more retrieval per unit of cosine**
  than the median shrunk team — which is exactly what whitening buys and exactly why the `pds`
  optimum is not the cosine optimum. **170 and 195 both optimised the member least sensitive to the
  quantity that is now the entire gap.** This is the same error 062 made with the CD4 dose, which
  207 corrected today for **+0.0044**.
- 181's audit adds: **`pds` has a seed-free interior maximum at gamma 0.50–0.55 on the *delivered*
  object** (+0.00046). So a step down from 0.65 may not cost `pds` at all.

## Work

### Phase 0 — known-answer rows (mandatory, abort on failure)
Re-run 170's instrument unchanged. Its eight gates are yours:
- `baseline_pre` / `candidate_pre` sha256 bit-identical to a15's shipped audit
  (`7678bab9…`, `0601336f…`) — **note 170's own trap: this only reproduces at `--chunk-rows 128`**,
  the builder's default, because PerturbMatch's centring mean accumulates per chunk;
- a15's post-whitening audit at γ 0.65 / λ 0 — 8 scalars, worst rel 1.1e−15;
- a19's at γ 0.80 / λ 0.007 — 6 scalars, rel 0.00e+00;
- `msekit` on the frozen 095 h1-test artifact = **2.055101829973639** against cell-eval2's own
  **2.0551018299736388**;
- `claim(0.25, φ0)` **0.11860** and the cap **0.24653**; challenge-panel `claim` at a15's policy
  **0.24610** against 164's 0.24581; challenge-panel ρ **0.938** against 164's 0.9466;
  per-context `Kc` **0.005084 / 0.005228 / 0.005919**.
- **Rule 11:** path and sha256 for every operator. Production is **29,394 × 17,942**, sha
  **`e2f65ecf`** (193's correction — §4.2's stated 12,000 × 7,760 matches **no file on disk**), and
  **nine** near-identically-named operators exist with max eigenvalues that differ in the 5th
  decimal. 170's gate 1 failed once for exactly this reason.
- **Then reproduce 208's Ĉ/C_p independently**: a15 **1.125**, `rc_callset` **2.940**. If you cannot,
  say so — the premise is then unverified.

### Phase 1 — re-locate BOTH cliffs on the current emission
For a15, `rc_hybrid` and `rc_callset`, per context: `A`, `B`, `claim`, `Ĉ/C_p`, `ρ = s²B/claim`, and
the distance to (a) the `ρ = 1` cliff and (b) `claim` leaving saturation. **Say explicitly which
cliff 208's "2.6× further away" refers to and whether it changes what 170 and 178 measured.** If it
does not, the honest answer is that gamma stays closed on the old mechanism and you should report
that — the sweep in Phase 2 is still worth running because of the endpoint change, but its verdict
will be different.

### Phase 2 — the sweep, with COSINE primary and `pds` secondary
Sweep γ over at least [0.30, 0.80] on the current emission (`rc_callset`), holding λ = 0 (193: there
is no ridge; production ships `WHITEN_LAM = 0.0`). At each γ report:
- **Δcos on PLOLO — primary.** 204 measured this instrument as resolvable to **±0.20% relative
  ≈ ±0.00006 overall**, which is 25–200× below anything shippable. Report its interval.
- **PLOLO `pds`** — 162's four folds, challenge frame. **`debench.Surface.pds()` must never size a
  `pds` arm** (188 Phase 3).
- **|t|-weighted `q_w`** (187 §5).
- `mse` through 096's `claim` route, per context, evaluating `min(claim, s²B)` **at every γ** —
  170's own correction was that linearising the `min` around the shipped point inverted its result
  (+0.0021 became −0.0039 at γ 0.50).
- **The maximiser of cosine and of `pds`, separately.** The whole hypothesis is that they differ, as
  they did for the CD4 dose.

Then **re-solve `VCC_ALPHA` at each γ** and report the price at the re-solved α, not at 0.982 —
170's V-shape in `mse` is partly an α-mismatch artefact if α is held fixed while `A` moves 14%.

### Phase 3 — six members, the gate, and the composition question
Rule 12, all six. `nmae` on h1-test with a **significance-selected pool** (188: `|lfc|·√cpm` at the
reference's BH depth, median ~89 genes/row; a top-N-by-magnitude pool **inverts** `nmae`'s sign).
Donors on h1-test: `k562_gwps` + `orion_hek293t` only — `data/donors/h1_2025` covers 85/85 and is
pure leakage. Rule 13: profile-class, so **no `reach` toll** (see the rule-13 table).

Gate **+0.0089** if the arm changes amplitude, **+0.0019** if you can show it does not. Lowering γ
raises `A` by ~14% over [0.65, 0.35], so expect the higher bar and say so.

**Coordinate, do not collide:** experiment **211** is re-solving the global α optimum against the
same moved cliff, and **210**'s `cd4f034` (built, submitting at the 00:00 UTC reset) already spends
some amplitude headroom in the favourable direction (delivered `A` **0.97970×**). Read all three
directories. Report your γ optimum **on `rc_callset`** and, if cheap, **on `cd4f034`'s delivered
profile** — if they differ, say which depends on what.

### Phase 4 — falsifiers
1. **170's own two-surface criterion**: no γ was net positive on both h1-test and Orion. Re-check on
   the current emission, and never quote an absolute from the Orion pair (its truth carries the
   long-gene artefact at 26.4% / 44.1% of effect energy).
2. **The leak-free fold**: report all four PLOLO folds including `kolf2`.
3. **A permuted-operator control**: apply a γ-swept operator built from permuted donor rows. 209's
   random control was the maximiser over all 36 of its arms — that is the outcome to watch for.
4. **The identity row**: γ = 0.65 must reproduce `rc_callset` exactly. 170 warns the identity check
   at production settings reproduces the baseline by orthonormality and *cannot* detect a mismatched
   basis — so assert the operator **sha**, not just the numbers.

## Deliverable
`RESULTS.md`: the known-answer rows, the two-cliff table across the three emissions, the γ sweep with
**cosine primary** and the two maximisers side by side, the six-member price at the re-solved α, all
four falsifiers, and a verdict in overall units. **"The cliff did not move and gamma stays closed" is
a full result** — it would close the operator for the third and final time, on the current emission
and on the right endpoint.

## Territory
210 owns the CD4 dose. 211 owns the global α scalar. 208 owns emission riders. **You own
`WHITEN_GAMMA` and the operator's exponent.** Do not re-fit the operator's basis (193 closed the
eigenvalue map; 195 closed learned transforms), do not touch the donor tables (194, 204, 209).
