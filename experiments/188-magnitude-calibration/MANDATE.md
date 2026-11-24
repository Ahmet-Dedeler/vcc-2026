# 188 — L1-optimal magnitude calibration: keep 169's `pds` gain, delete its `nmae` cost

You are a research subagent. Work ONLY in `experiments/188-magnitude-calibration/`.
No `vcc`, no `cell-eval2`, no submission build, no writes to `submissions/` or outside your
directory. **Peak RSS ≤ 8 GB** (four other agents are running on this machine). `uv run python`.

## READ FIRST, in this order
1. `experiments/169-sharpening-six-member/RESULTS.md` **in full**. This experiment exists to
   resurrect its arm. §1, §2, §5, §6.1–§6.4, §7, §10, §11 all matter.
2. `docs/CURRENT_STATE.md`, the `2026-09-03 02:15` entry — the `nmae` landscape identity.
3. `AGENTS.md` §4, §11.1 rules 1–13 with the audit-181 amendments (especially rule 12,
   six-member pricing, and rule 13's arm-class split).
4. `experiments/168-sharpening/`, `experiments/098-nmae/` §4 (the gate-depth artefact),
   `experiments/162-retrodicting-validator/` (PLOLO, the gate).

## The situation, precisely

169 measured the largest single signal this project has ever produced. Sharpening the per-gene
effect magnitude, `h_g → sign(h_g)·|h_g|^q` with the row norm restored, reads on PLOLO — the only
`pds` instrument validated against the board (6/6 signs, LOO RMS 0.001494 scaled) — a **Δ`pds` of
+0.019452 scaled at q = 1.15**, positive on all four folds including the leak-free `kolf2` fold, at
blob shares of 0.0–0.6% so it is **not** an Orion artefact. That is **2.9× a15** and **1.54×** the
Gate-4 bar that no other geometry arm in this project has ever cleared. `mse` is neutral
(Δ`claim` −0.00001, `rA` 0.99903). The `reach` toll rule 13 charges was **falsified** for this arm
class (+0.000605 overall measured, whole CI above the toll). `fid` and `jac` are noise.

**It died on one member, `nmae`, and §5 gave the exact mechanism.** The 02:15 identity states it in
closed form: for a reference DE gene with true magnitude `|t|`, ours `|p| = k|t|`, sign accuracy `q`,

```
nmae(k, q) = q·|k − 1| + (1 − q)·(k + 1)
```

minimised at **k = 1**, with **asymmetric slopes**: `−(2q − 1) ≈ −0.074` below and **`+1` above**.
A power law with exponent > 1 inflates the largest magnitudes, and on h1-test **9% of gated genes
already have `|h| > |r|` and carry 22.6% of the L1 error**. Sharpening charges those at slope +1
while earning slope −0.074 on the other 91%: a 13:1 asymmetry. That, and only that, is the NO-GO.

169 also found the transform space matters: the builder's `h` space is `log1p(TS·p)`, linear below
~20 CPM and logarithmic above, so a power of `h` concentrates on **high-expression** genes, not on
**large-fold-change** genes. Rebuilt in **fold-change space** the arm is better on every member —
`nmae` 5.4× cheaper — and its six-member sum with `nmae` from h1-test is **+0.000930**, against a
gate of +0.0019. **It misses by 0.001.**

## The hypothesis

**A power law is the wrong magnitude map.** `nmae` is an **L1** loss on magnitudes, so the map that
minimises it is the **conditional median of the true magnitude given the predicted one**, not a
monomial. Fit that map and you get, by construction, the magnitude reweighting the data supports —
expansive where our estimate is over-shrunk, saturating where it is not — instead of one that is
expansive everywhere including on the 9% of genes that are already over-predicted.

If the conditional median is convex over the bulk of the range, it delivers 169's `pds` contrast
for free; if it is concave at the top, it delivers less `pds` but stops paying the +1 slope. Either
way the arm becomes *correctly specified* rather than *approximately right and catastrophically
wrong in one tail*. Our predictor is heavily shrunk — 169 §6.2 measures median `|h|/|r|` at
**0.1266** (HCT116) and **0.2825** (h1-test) — so there is a lot of room for an expansive map.

## What to do

### Phase 0 — known-answer rows (mandatory, abort on failure)
- The identity map must reproduce 169's q = 1.00 rows: 16/16 bit-identical PLOLO folds, and the
  a15 → a23 mass-fix PLOLO price **+0.004014** scaled.
- Reproduce 169's q = 1.15 numbers in **both** spaces: PLOLO Δ`pds` **+0.019452** (`h` space) and
  **+0.017631** (`lfc` space), and h1-test emitted Δ`nmae` **+0.04640** / **+0.00860**.
- `claim(a15 as shipped)` = **0.24581**; a15 raw `mse` / `score_mse` = **0.88228** / **0.11112**.

If you cannot reproduce 169's numbers you cannot improve on them — stop and report.

### Phase 1 — fit the map, leakage-safe
On the donor lines where truth exists, learn `m(·)` mapping our predicted per-gene magnitude to the
conditional **median** of the true magnitude, **leave-one-line-out**, in **fold-change space**
(169 §6.4 shows why). Use a monotone, low-capacity estimator — isotonic regression or a binned
median with monotonicity enforced — not a parametric family you have to choose. Condition on
predicted magnitude alone first; only then consider one extra covariate (control expression decile,
or donor agreement), and report whether it earns its complexity on held-out lines.

Report the fitted map's shape explicitly: is it convex, concave, S-shaped? Where does it cross the
identity? What fraction of the gated mass does it move up versus down? **The shape is the finding**,
whatever the arm's fate.

### Phase 2 — price all six (rule 12), on the instruments each member is validated on
| member | instrument | incumbent |
|---|---|---|
| `pds` | **PLOLO**, challenge frame, 162's pre-committed 4-fold unweighted mean | a15/a23 rows |
| `mse` | 096's `claim` route on `data/challenge/context_{A,B,C}` | `claim` 0.24581 |
| `nmae` | h1-test **emitted**, ≥5 matched seeds, **and** HCT116 at BH 0.35 — report both, take the worse (169 §6.2: the board-validated surface is also the worse one) | a15 board raw 0.944251 |
| `reach` / `fid` / `jac` | h1-test emitted, matched seeds | 169 §4 |

Predict `nmae` **from the identity as well as from the surfaces**, and compare. The identity needs
the fraction of DE-pool L1 mass your map pushes above `k = 1`; you can compute that directly on the
donor folds. **If the identity and the surfaces disagree, that is a finding and you must report the
disagreement rather than pick one.**

### Phase 3 — resolve the open `pds` contradiction 169 flagged
169 §6.4 records that the h1-test DE-table `pds` reads **−0.0228** for the lfc-space arm where
PLOLO reads **+0.017631**. It trusted PLOLO (correctly — it is the board-validated instrument in
the scorer's own space with 269–294 decoys, against the DE table's 85), but flagged the sign
disagreement as something to resolve before the family is revisited. **Resolve it.** The likely
cause is that the two statistics are computed in different spaces on different gene sets with
different decoy pools; demonstrate that rather than assert it, e.g. by computing an `h`-space
`pds` on the h1-test panel with a matched decoy count.

### Phase 4 — the gate
Composed slot gate, with rule 13's arm-class split (this is a **profile-class** arm, so the
`−0.0146` `reach` toll does **not** apply — 169 falsified it for this class with matched seeds):
```
predicted = PLOLO d_pds/6 × 1.011  +  096 d_mse/6  +  measured nmae/reach/fid/jac terms
SHIP only if predicted > +0.0019 (shape-only)
```
Note the standing warning: **a13** is this family's nearest board point and it lost **−0.0203
overall** with raw `nmae` 1.0185, and **both harnesses were off by 0.075 raw and in the wrong
direction on it**. Audit 181's rule 1 applies — if you cannot retrodict a13, say so explicitly and
treat your `nmae` number as unvalidated. 169 notes a13 is "not reconstructible"; if that is still
true, state it as a limitation rather than working around it.

## Deliverable
`RESULTS.md`: known-answer table, the fitted map's shape, the six-member price on both surfaces and
from the identity, the Phase-3 resolution, the gate verdict, and — if it clears — the exact builder
invocation plus a preregistered per-member prediction written **before** anything is built. Do not
build. Do not submit.

## Territory
184 owns Tahoe / cross-line direction. 185 owns the emission law and the DE call set. 186 owns
competitive intel. 187 owns the global amplitude scalar α and the `mse`/`nmae` joint solve — **you
own the per-gene magnitude SHAPE**; coordinate by reading their `RESULTS.md` if they land before
you, and do not re-solve their scalar.
