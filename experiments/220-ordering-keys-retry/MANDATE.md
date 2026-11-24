# 220 — The ordering-key closures were noise-limited. The instrument just got 77× tighter.

You are a research subagent. Work ONLY in `experiments/220-ordering-keys-retry/`.
No `vcc`, no `cell-eval2`, no submission build, no writes outside your directory.
Peak RSS ≤ 10 GB. `uv run python`.

## READ FIRST
1. `experiments/214-cheap-positives/RESULTS.md` — **its by-product 1 is why this experiment exists.**
2. `experiments/189-per-gene-reliability/RESULTS.md` and `190-regulon-ordering/RESULTS.md` — the two
   closures being re-opened. **Know exactly what each measured and on what instrument.**
3. `experiments/185-rank-controlled-emission/BOARD-RESULT-rc_hybrid.md` and
   `experiments/202-callset-rider/BOARD-RESULT.md` — the two board points that built the new
   instrument and proved the channel converts.
4. `experiments/215-leader-decomposition/RESULTS.md` — its by-product on how little of this channel
   the field uses.
5. `AGENTS.md` §8.1, §8.4, §11.1 rules 1–13 with the rule-13 table.

## The argument

**Both ordering-key closures were noise-limited, and the noise just fell 77×.**

214 measured that `reach` under the shipped `rc_callset` emitter has a **three-seed sem of 0.00009**,
against 200's **0.00695** on the older emitter — a **77× tightening**. The ordering channel finally
has a low-variance instrument, and it did not exist when 189 and 190 ran.

Re-read their verdicts against that:

| experiment | measured | its own uncertainty | verdict then |
|---|---:|---|---|
| **189** per-gene cross-line reliability key | **+0.00408** raw reach | **± 0.00424** | NO-GO — inside noise |
| **190** regulatory-prior key | **+0.00137** raw reach | CI [−0.0014, +0.0052] | NO-GO — and its shuffled control scored *higher* |

**189's effect was 0.96σ on its instrument. On an instrument 77× tighter it would be ~80σ** — if the
effect is real and not an artefact of the old emitter. That is the question.

**And the channel is nearly unexploited by anyone.** 215 measured **9 arms in 707 board arm-pairs**
that moved the strict ordering channel — **two of them ours**. We are 1st in the world on `jac`
and 29th on `reach`, and our `reach` residual against the field's own conditional went
**−1.92σ → +0.70σ** in one day. This is the cheapest live territory we have.

**A caution that must shape the work.** 189's own closure was **structural, not statistical**: it
argued `reach` needs **90% purity**, the median pure prefix is **3 genes**, and lifting sign accuracy
inside `|H|`'s top decile from 0.587 to 0.669 **does not cross 0.90**, capping the whole family at
**≤ +0.0005 raw reach**. 190 derived the same cap independently. **If that structural argument
holds, a tighter instrument changes nothing** — it just measures a small number precisely.
**Deciding between "noise-limited" and "structurally capped" is the entire experiment.**

## Work

### Phase 0 — known-answer rows (mandatory, abort on failure)
- Reproduce the instrument itself: `reach` under `rc_callset` at three seeds, and its **sem
  0.00009**. If you cannot reproduce the tightening, stop — everything here rests on it.
- Reproduce the incumbent rows: production `reach` **0.18283**, `|H|` key **0.23533**, sig-first
  **0.21914**, oracle **0.64553**, H-oracle **0.62119** (171/172).
- Reproduce 189's key at its published **+0.00408** and 190's at **+0.00137** on the *old* emitter,
  so the old and new instruments are directly comparable on the same arms.
- Donors on h1-test: `k562_gwps` + `orion_hek293t` **only** — `data/donors/h1_2025` covers 85/85 of
  h1-test's targets and is pure oracle leakage. Scale pools **proportionally**, never
  `Surface(n_real_target=N)`.

### Phase 1 — settle "noise-limited" vs "structurally capped"
Re-measure 189's and 190's keys on the **new** instrument and report the effect with its new
uncertainty. Then, separately and decisively, test the structural claim:
- What **is** the median pure prefix under `rc_callset` (189 measured 3 under the old emitter)?
  `rc_callset` changed the call set and the delivered order, so this number may have moved.
- What sign accuracy would be needed at the front to cross `P₀ = 0.90` at depth 5, 10, 20?
- Does any available key reach it **for a non-trivial share of perturbations** — not on average?
  189's cap argument is about the *median* perturbation; `reach` averages `k*/n_real` over
  perturbations, and 172 measured that **94.4% of h1-test's reach lives in 32 of 85 targets**.
  **A key that helps only the perturbations where the front is already nearly pure could pay even
  if the median is untouched.** That is the loophole in the structural argument and nobody has
  tested it.

### Phase 2 — if anything survives, price it
Six members (rule 12). The ordering channel's great virtue is that `pds`, `mse` and `nmae` are
**pinned by construction** when the per-gene panel totals are held — board-verified to **sixteen
significant figures** across `rc_hybrid` and `rc_callset`. **Prove that with a digest for your arm,
never assert it.** Rule 13: an ordering arm that holds the totals identical carries **no** `reach`
toll (the rule-13 table). Gate **+0.0019**.

Convert with the factor the board fixed for this class: **h1-test → board `reach` ×1.290** (185's
board point), **not** 164's stale 0.4277.

### Phase 3 — falsifiers
1. **Permutation control at the new precision.** 190's shuffled regulon scored *higher* than the
   real network; 189's shuffled control tied its call-set arm. At 77× tighter noise these controls
   become sharp — **run them, and if a control still ties or beats the arm, the family is dead for
   good and that is the result.**
2. **The old-emitter comparison** must reproduce the published numbers, or the two instruments are
   not measuring the same thing.
3. **Do not select the key on the surface you score it on.** 189's own λ was selected on its
   evaluation surface, which it flagged; do better.

## Deliverable
`RESULTS.md`: the reproduced instrument and its sem, 189's and 190's keys re-measured old-vs-new,
the structural test with the per-perturbation stratification, and a verdict. **"Structurally capped,
and the tighter instrument just measures the cap precisely" is a full result** and would close the
ordering-key family permanently — which is worth having, because it is currently the only channel
where we lead the field and we should know whether there is anything left in it.

## Territory
216 owns the CD4 dose. 217 is building. 218 owns acquisition. 219 owns the closure audit —
**it will hand you 189/190 rather than re-price them itself.** 213 owns the transport law.
**You own ordering keys.**
