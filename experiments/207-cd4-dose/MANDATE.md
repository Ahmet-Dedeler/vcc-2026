# 207 — The CD4 residual dose was chosen by maximising `pds`, the member least sensitive to direction. Sweep it against cosine.

You are a research subagent. Work ONLY in `experiments/207-cd4-dose/`.
No `vcc`, no `cell-eval2`, no submission build, no production changes, no writes outside your
directory. Peak RSS ≤ 12 GB. `uv run python`.

## READ FIRST
`AGENTS.md` §4.2, **§8.2** (direction is the entire remaining gap), §11.1 rules 1–13. Then:
1. `experiments/206-cluster-method/RESULTS.md` **§6 and its H1** — the argument this mandate
   executes. **Verify it; do not inherit it.**
2. `experiments/197-cosine-inversion/RESULTS.md`, especially §1.2 — the CD4 residual is measured
   there as **the only resolvable cosine gain in eight board arms, +6.8%**.
3. `experiments/062-cd4-post-whiten/`, `047-cd4-direct/`, `050-cd4-transfer-confidence/`,
   `051-cd4-breadth-blend/`, `054-cd4-union-hct/`, `060-cd4-union-k562/`, `085-kolf-pm-gate/` —
   the CD4 family. **062 is the one that set the shipped dose.**
4. `experiments/041-ccle-match/RESULTS.md` — context A = **JURKAT**, under all eight feature/rank
   variants, which is why a T-lineage donor is matched at all.
5. `experiments/096-mse-floor/A15-BOARD-RESULT.md` — a15's board point, the arm that shipped CD4.

## The argument

**Everything else has closed and this one thing was never tuned against the channel it pays in.**

- 197 §1.2: across eight board arms, **the CD4 residual is the only profile change that measurably
  raised our transfer cosine — by +5 to +7%.** Nothing else did.
- 206 §6: the residual ships at effective fraction **0.0528**, and **062 selected that dose by
  maximising LODO `pds`** — the member 197 proves is *least* sensitive to direction, because `pds`
  is a rank statistic on which whitening buys us 26% more retrieval per unit of cosine than the
  field. **The dose was optimised on the wrong objective.**
- 206 §5: three *mutually independent* pipelines sit above cosine 0.11 (pairwise distance
  0.43–1.02), two of them lone individuals with no organisation, and the implied `k_eff` of 10–14
  **exceeds the ~11 CRISPR-perturbed human cell lines that exist in total.** So nobody is on 159's
  `1/k` curve; the edge is not a private atlas, and **the documented exit from that curve is a
  lineage-matched donor.**
- Nine teams name the Marson/Zhu genome-scale CD4⁺ T-cell atlas publicly, and **GISL published a
  board-measured +0.01325 overall for shipping it as a LIBRARY** rather than as a small residual.
  It covers **284–291/300** of the panel. We already hold it.

## Work

### Phase 0 — known-answer rows (mandatory, abort on failure)
- **The built-in row: at the shipped fraction 0.0528 your pipeline must reproduce a15's cosine
  0.07925** (197's value for the shipped arm) and PLOLO's a15 `pds`. 206 names this as the
  known-answer row and it is the one that proves your dose axis is the same axis 062 used.
- Reproduce PLOLO's a15 → a23 mass-fix price **+0.004014** scaled (162).
- Reproduce 062's own headline: HCT116 **+0.00482**, HEK293T **+0.00242**, frozen K562 **+0.00441**
  at fraction 0.05, and state which of those numbers your axis reproduces and which it cannot.
- **Rule 11:** record a path and a sha256 for every operator. The shipped operator is
  **29,394 × 17,942**, sha **`e2f65ecf`** (193's correction — §4.2's stated 12,000 × 7,760 matches
  no file on disk), and nine near-identically-named operators exist.

### Phase 1 — sweep the dose against the RIGHT objective
Sweep the CD4 effective fraction over a wide grid (at least 0 → 0.40, log-spaced, including the
shipped 0.0528 as the identity row). At each dose report:
- **implied transfer cosine** by 197's route — **this is the primary endpoint**;
- **PLOLO `pds`** (162's four folds, challenge frame; `debench.Surface.pds()` must never size a
  `pds` arm — 188 Phase 3);
- **|t|-weighted `q_w`** (187 §5, the quantity `nmae` moves through);
- and the **maximiser of each**, separately. The whole hypothesis is that these maximisers differ.
  If the cosine maximiser is also 0.0528, the hypothesis is dead and that is a clean result.

### Phase 2 — the FORM, not just the dose
GISL ships it as a **library** (a fourth donor entering the consensus) and measured +0.01325 on the
board; we ship it as a **post-whitening residual** at 5%. Those are different estimators, not two
doses of one. Price both:
- residual at its cosine-optimal dose;
- CD4 as a **fourth donor line** in the equal-weight consensus, whitened with the rest;
- and, if it is cheap, the two composed.

Note what 206 §5 implies about *why* a lineage-matched donor escapes 159's `1/k` law: the law counts
*independent* lines contributing to a context-independent mean, and a donor matched to the recipient's
lineage contributes to the *line-specific* component the law treats as noise. **Say whether your
measurements support that mechanism**, because it decides whether this generalises to D/E/F (where
the lineages are unknown until Oct 22) or is specific to A/Jurkat.

### Phase 3 — the six-member price and the gate
Rule 12. `mse` through 096's `claim` route (a profile change; `mse` is **not** exempt — a19 lost
−0.01352 exactly here, and we currently ship `claim` at **100.0% of the cap**, so say what the dose
does to that). `nmae`/`reach`/`fid`/`jac` on h1-test with a **significance-selected pool** (188:
`|lfc|·√cpm` at the reference's BH depth, median ~89 genes/row — a top-N-by-magnitude pool inverts
`nmae`'s sign). Donors on h1-test: `k562_gwps` + `orion_hek293t` only; `data/donors/h1_2025` covers
85/85 of h1-test's targets and is pure leakage.

Rule 13's `reach` toll does **not** apply to a profile-class arm (169 falsified it there; the
`rc_hybrid` board point confirmed the toll belongs to profile-*changing emission* arms). Gate:
**+0.0019** shape-only, **+0.0089** amplitude-changing — a dose change is amplitude-changing unless
you re-solve `VCC_ALPHA`, so **re-solve it and say what that cost**, exactly as 203 had to.

### Phase 4 — falsifiers
1. **Target-shuffled CD4 control** at every dose — 206 names this and it is the decisive one. If a
   shuffled donor buys the same cosine, the gain is amplitude or geometry, not biology.
2. **Lineage specificity**: if CD4 helps because it is T-lineage-matched to A/Jurkat, the gain should
   be **larger on context A than on B (HeLa) or C (CAL33)**. Measure per context. If it is uniform,
   the lineage story is wrong and the mechanism is something else — report that.
3. **The dose must not be selected on the surface it is scored on.** Select on donor folds, report
   on the held frame, and say which is which.
4. 062's LODO number and yours must be reconcilable; if the cosine-optimal dose is far from 0.0528,
   show that 062's `pds` curve is genuinely flat there rather than assuming it.

## Deliverable
`RESULTS.md`: the known-answer rows, the dose sweep with **cosine as the primary endpoint** and the
three maximisers side by side, the residual-vs-library comparison, the six-member price at the
chosen operating point, all four falsifiers, and a verdict in overall units. If it clears the gate,
write the exact builder invocation and a preregistered per-member prediction, **and stop** — the
coordinator builds and submits.

206's H1 estimate is **+0.004 to +0.0133 overall**. The gap to the cluster is +0.01895 and to first
place +0.03740. **This is currently the highest-value open arm in the project.**

## Territory
You own the CD4 residual: its dose, its form, and its interaction with the operator's amplitude.
Do not tune the whitening operator (193 closed it), the donor tables (194), or the emission law
(185/191/202).
