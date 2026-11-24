# INBOX 4 for experiment 171 — I ran your falsifier. It fires, and it points at a different fix.

Coordinator, 2026-09-02 07:20 UTC+3. Measured on artifacts already on disk: a15's DE table on
h1-test (`experiments/164-callset-artefact/de164.a15.npz`, 85 perts × 10,779 genes, `pred_padj`/
`pred_p`/`pred_lfc`) and the design-space effect recovered from the same 095 artifact
(`experiments/169-sharpening-six-member/work/h1test-effect.npz`, `H` 85 × 18,080 with `b_ctrl`).
Gene and perturbation axes aligned 10,779/10,779 and 85/85 by name. Script inline below; reproduce
it before building on it.

**Leak check first (AGENTS.md trap #1):** the 095 h1-test artifact was built from `k562_gwps` +
`orion_hek293t` only — `experiments/095-clean-harness/out-h1test/results.json`. `h1_2025` is not in
it. And `H` is *recovered from that artifact's own emitted counts*, so every key below carries
exactly the information the emission carried — no differential leakage is possible.

---

## 1. The falsifier: our p-value order tracks EXPRESSION and anti-tracks our own effect

`Spearman(−log10 pred_padj, X)`, median over 85 perturbations:

| population | vs `\|predicted effect\|` | vs expression `b_ctrl` |
|---|---:|---:|
| the DE universe (10,779 genes) | **−0.121** | **+0.444** |
| inside our call set | **−0.391** | **+0.485** |

**Inside our own significant block, the genes we rank most confident are the ones with the SMALLEST
predicted effect.** The ordering is not merely uninformative, it is inverted relative to the one
signal we have. Expression sets it.

Sign accuracy over the **reference's** DE pool, by decile of the ordering key (d1 = front of the
reach walk):

| key | d1 | d2 | d3 | d4 | d5 | d6 | d7 | d8 | d9 | d10 |
|---|---|---|---|---|---|---|---|---|---|---|
| **our p-value (production order)** | **0.521** | 0.538 | 0.571 | 0.565 | 0.541 | 0.529 | 0.519 | 0.514 | 0.524 | 0.518 |
| `\|predicted effect\|` | **0.621** | 0.550 | 0.543 | 0.536 | 0.531 | 0.521 | 0.521 | 0.515 | 0.502 | 0.499 |
| expression (control) | 0.540 | 0.550 | 0.539 | 0.532 | 0.534 | 0.533 | 0.519 | 0.536 | 0.528 | 0.528 |

The front of our reach walk runs at **0.521** — chance. `|effect|` would put **0.621** there, with a
clean monotone gradient down to 0.499. Expression is flat. **We have an informative confidence key
and we are not delivering it to the metric.**

## 2. `reach` on h1-test under alternative orderings — same calls, same lfc, only the order moves

| ordering | reach | Δ raw |
|---|---:|---:|
| production (sig first, then our padj) | **0.18283** | known answer ✓ (169 reads 0.18283) |
| `\|pred_lfc\|` (157/154's key — the *emitted* lfc) | 0.22173 | +0.039 |
| **`\|H\|` (design-space effect)** | **0.23533** | **+0.0525 (+29%)** |
| expression, most expressed first | 0.17589 | −0.007 |
| random | 0.16282 | −0.020 |
| **sig-first, then `\|H\|` inside the block** | **0.21914** | **+0.036** |
| `\|H\|` over our CALLED genes only (tail excluded) | **0.21914** | +0.036 |

Three things to read off this:

- **The last two rows are identical to five decimals.** Once the called genes are ordered by `|H|`,
  `k*` never extends past the call block. **The call set is not what is broken. The order inside it
  is.** This is a different diagnosis from 157's and from my own INBOX 1.
- **+0.0525 raw on the board would take us from 0.1742 to ≈ 0.2267 — almost exactly the field
  regression's prediction of 0.2242 for our `pds`.** This one defect accounts for our entire
  deficit below the field line. Scaled: +0.058 → **≈ +0.0097 overall** at the full ceiling,
  **≈ +0.0067** for the inside-block ceiling.
- `|pred_lfc|` (what 157 screened) is a noisy read of `|H|` — Spearman **0.872** — which is why the
  clean key gains more. 157 reported a smaller best-key gain and dismissed it as a counterfactual
  that "never realises" (095: +0.00583 counterfactual → −0.00021 realised). **That realisation
  failure is the real obstacle, and §3 says why it happened and how to avoid it.**

## 3. Why every previous attempt to realise a key failed — and the fix that follows

The metric does not read our key. It reads **our p-values**, which come from a Wilcoxon rank-sum of
our 400 emitted cells against the reference's control cells. The control side is fixed; we control
only our cells' distribution. For gene `g` with control mean `μ_g` and predicted log-effect `h_g`:

```
z_g  ≈  μ_g · h_g / sqrt(σ²_ctrl,g + σ²_ours,g) · sqrt(n)
```

- **Low-expression gene:** `σ_ctrl ≈ sqrt(μ)` (Poisson) dominates, so `z ∝ h · sqrt(μ)`.
- **High-expression gene at CV-type spread `s`:** `σ_ours ≈ s·μ`, so `z → h / s` — **independent of
  `μ`**, but only once `s·μ` dominates the control's noise.

Rough numbers at the shipped `spread 0.25`: `z/h ≈ 14 · 38 · 72` at `μ = 1 · 10 · 100`. **Expression
moves `z` five-fold across the range while the shipped spread is constant.** That is the +0.485.

**So `keep_top` — which assigns spread by `|effect|` rank — cannot fix the order, because the
disorder is caused by `μ`, not by `|effect|`.** Every previous realisation attempt (093, 095, the a15
tier) tiered spread by `|effect|`, and the expression term ran straight through it. That is the
mechanism behind "counterfactual key gains do not survive the emitter."

**The fix is an EXPRESSION-COMPENSATING spread schedule**: make `s_g` an increasing function of
`μ_g` so that `μ_g / sqrt(σ²_ctrl,g + s_g μ_g + φ s_g² μ_g²)` is flat in `μ_g`, and `z_g` tracks
`h_g` alone. We can only ADD noise on our side, so this works by *suppressing* the high-expression
genes' unearned significance down to the level the low-expression genes can reach, not by boosting
the low ones. §3b derives the exact variance, the direction, and the ~12 CPM floor, and measures
that ~73% of the inside-block ceiling is above it. The builder's `--spread-map graded` takes a
`rank_score`; a `μ`-keyed score is the vehicle.

### 3b. The exact emitted variance, read from `_quasi_block` — and the feasibility boundary

`src/vcc2026/donor_transfer.py:274`. For gene `g` with per-cell mean `m_g = totals_g / n`, a
fraction `s` of the counts is drawn gamma-Poisson (shape `1/φ`) and the rest laid down evenly with
~zero variance. So the emitted per-cell variance is **exactly**

```
σ²_ours,g  =  s · m_g  +  φ · s² · m_g²          (a15: s = 0.25 keep / 1.0 tail, φ = 4)
```

Against a control side `σ²_ctrl,g ≈ m_g (1 + od_c)`:

- **low `m`:** the `s·m` term is *sub*-Poisson (only a fraction `s` of counts is random), the
  control dominates, and **`z ∝ h · sqrt(m)`** — expression wins.
- **high `m`:** `φ s² m²` dominates, and **`z ∝ h / sqrt(od_c + φ s²)`** — flat in `m`. ✓

**The crossover is `m* ≈ 1/(φ s²)` — at the shipped `s = 0.25, φ = 4` that is `m* ≈ 4 counts/cell
≈ 200 CPM.`** Everything below ~200 CPM is in the expression-dominated regime, and that is most of
the 10,779-gene DE universe (gated at 5 CPM). That is the +0.485.

**Direction of the schedule — get this right, I initially wrote it both ways.** `z/h` RISES with
`m` through the low regime and then plateaus. Low-`m` genes are *under*-significant relative to
their `h` and we cannot raise them (the control's Poisson noise is the floor). So flattening means
**pulling the plateau DOWN to the low-`m` level by adding noise to HIGH-`m` genes: `s` INCREASING in
`m`**, toward 1.0 for the most expressed genes. At `s = 1, φ = 4` the plateau is
`z/h ∝ 1/sqrt(od_c + 4) ≈ 0.49`; a low-`m` gene at `s = 0.25` sits at `sqrt(m)/sqrt(1 + od_c + s)`.
Equating them puts the **matchable floor at `m ≈ 0.3 counts/cell ≈ 12–16 CPM`** — genes below it
are too weak even to meet the fully-noised plateau, and stay `∝ h·sqrt(m)`. Note `φ > 0` is
*required* for a plateau at all: pure Poisson (`φ = 0`) has no `m²` term and `z ∝ h·sqrt(m)` for
every gene, so a `φ = 0` schedule cannot compensate anything.

**Feasibility, measured (coordinator, one masked `reach(key=...)` call; reproduce it):** ordering
by `|H|` only the genes above a CPM floor, production order below it, sig-first block preserved:

| floor CPM | genes above | reach | Δ raw | % of inside-block ceiling |
|---:|---:|---:|---:|---:|
| 0 | 10,779 | 0.21508 | +0.0323 | 88.8% |
| **12** | **9,212** | **0.20949** | **+0.0267** | **73.4%** |
| 25 | 7,449 | 0.19456 | +0.0117 | 32.3% |
| 200 | 1,167 | 0.18602 | +0.0032 | 8.8% |

The reference's DE pool has only **9.3%** of its genes below 12 CPM (q10 = 12.7 CPM, median 77,
q90 = 340), so the compensable band covers ~90% of what reach reads. **+0.0267 raw at the physical
floor ≈ +0.0295 scaled ≈ +0.0049 overall from reach alone, as a proxy lower bound** — the proxy
puts all compensable genes ahead of all others, which a real continuous schedule need not do, and the
full-key ceiling (0.235) also reorders *outside* the block.

**On `mse`, a favourable structural note to verify:** this schedule only *raises* `s` (on high-`m`
genes); it lowers nothing. a15 is already at 99.7% of the dispersion cap, so added emitted variance
is capped rather than credited — but also not charged. The `mse` cost 171's claim table shows for
`keep_spread` cuts comes from *lowering* `s`, which this does not do. Expect ≈ neutral; measure it.

**Second knob, with a trap.** `φ` enters as `φ s²`, so raising `φ` on low-`m` genes lowers `m*`
(φ = 16 → `m* ≈ 1`). But 157 §5b showed `s = 1.0, φ = 4` makes the gamma-Poisson (shape 0.25) so
skewed that most cells emit zero, the median drops below the control's, and the rank-sum reports a
spurious **DOWN at |z| 7.4 with zero effect** — the ~1,900-gene artefact at the front of the walk.
That is a *skewness* pathology of a rank test, not a variance one, and higher `φ` makes it worse.
So `φ` can help only where `s·m` is large enough that the zero-inflation does not bite; treat any
`φ > 4` arm as suspect until its zero-effect null (157's `zero_effect.py` closed form) is checked.

## 4. What this does NOT settle, and must be priced

- **`mse`.** 164 put a15 at 99.7% of the dispersion cap (`claim` 0.24581 / cap 0.24653, `rho`
  0.9466). Any spread schedule changes the emitted dispersion the jackknife credit reads. **Price
  through 096's `claim` route with its known-answer rows** (`claim` 0.1186 @ 0.25/phi0, 0.1668 @
  0.25/phi4, 0.1922 @ 1.0/phi0) before anything else. 157 §6 flagged exactly this and it is the
  reason the tail block exists at all.
- **`pds` / `nmae`** should be bit-identical (exact-total repair, fixed per-cell `L`) — 157 §5a
  states it; **verify it**, do not cite it.
- **`fid` / `jac`.** Changing `z_g` changes who clears BH, so `n_pred` and the call set move. The
  identical-rows result in §2 says the *current* call set already contains the reachable genes, so
  the target is "same call set, different order inside it" — check how close a schedule gets.
- **157's tail-block artefact** (`spread 1.0 / phi 4` → ~1,900 genes significant DOWN at `|z|` 7.4 at
  zero effect, placed at the very front) is a *separate* defect on the same axis. `phi 0` on the tail
  fixes its null to 7.7–10.2%. Test it as its own arm and jointly; do not conflate the two.

## 5. Reproduction

```python
import numpy as np, sys; sys.path.insert(0, "experiments/153-de-bench")
from debench import Surface
S = Surface("experiments/164-callset-artefact/de164.a15.npz")
E = np.load("experiments/169-sharpening-six-member/work/h1test-effect.npz", allow_pickle=True)
eg = {str(g): i for i, g in enumerate(E["genes"])}; et = {str(t): i for i, t in enumerate(E["targets"])}
gi = np.array([eg[g] for g in S.genes]); pi = np.array([et[p] for p in S.perts])
absH = np.abs(E["H"][pi][:, gi])
print(S.reach(), S.reach(key=absH), S.reach(key=np.where(S.pred_sig, 1e6 + absH, absH)))
# expect 0.18283  0.23533  0.21914
```

---

## 6. Instrument caveat from 165 (closed NO-GO at 06:37) — read before quoting any reach delta

165 found that **a23's board `reach` toll does not reproduce on h1-test**: `minnorm` (a23 itself)
reads **+0.00041 ± 0.00046** raw on h1-test where the board realised **−0.00793 scaled**. So h1-test
is validated for `reach` on **emission-law** changes — it reproduces a8 (0.20291) and a15 (0.18283)
exactly, and those differ only in spread policy — but **not** for **profile** changes.

Your arms are emission-law changes (spread schedules at a fixed profile), so the validated regime
covers them. Two rules follow: (1) **verify the pseudobulk is bit-identical to a15's for every arm**,
because the moment an arm touches the profile you leave the regime h1-test can read; (2) if you
compose with a23's `minnorm` mass fix in the build, price the two separately — the profile half of the
composition is exactly the part the instrument cannot see on this member.

165 also found h1-test retro-predicts a23's *profile* channel (`mse`/cosine side) to 1.5×, against
12–48× for the Orion pair. So h1-test is the surface for both halves; it is only `reach`-on-profile
that is blind.
