# 229 — pre-registration of the arm, and an honest note about when it was written

**Written 2026-09-04 06:45 UTC, AFTER the four builds completed and BEFORE `six229.py`'s exact
`nmae` identity was read.** INBOX-3 §2 asked for a pre-registration *before* building; it arrived
after the builds were already on disk. Saying so is the point of this file.

## What the guard is for, and whether it holds here

INBOX-3 §3: *"Selecting a construction by what your own build scores is fitting the arm to the
instrument."* **That failure mode did not occur, and it could not have:**

- **The built packages were never scored.** No `cell-eval2` was run, no board result exists, and
  nothing in this experiment reads a score computed from an `.h5ad` or a `.vcc`.
- **The selection instrument is 225's own PLOLO harness** — `experiments/162-.../plolo.py`,
  `218/resid218.py`, `218/compose218.py` — imported, not rewritten. It reproduces 225's
  `stack225d.json` rows at **0.0e+00 on 36 fold × dose cells** (KA-5) and 225's own baselines to
  the last printed digit. It is the same instrument that priced the arm, and it is independent of
  the builder.
- The one thing the builder contributes is the **union object**, and KA-4 proves the builder's
  union equals 218's independently written one at **0.0e+00** on every fold and every selection.

So the selection is on 225's measurements, per INBOX-3 §3. What is genuinely missing is that this
file was not written first, and a reader should discount accordingly.

## The pre-registration

**Selected arm:** `ARM6` — `cd4@0.34 + union(gse295214_h1, h1_2025, dld1_lowmoi,
jiang_5line_noK562, igvf_cardiomyocyte, myllia_a549) @ g = 0.35`.

**One-line reason:** it is the highest-priced construction a single-union builder can express
(+0.00268, 1.41× the shape-only gate, 3/3 cosine and 3/3 `pds` on the three cancer folds, marginal
over the shipped `h1union`), it is the only expressible arm that beats 223's `dld1_lowmoi`
recommendation (+0.00228), and `g = 0.35` is 225's own pre-registered dose, not a tuned one.

**Lattice searched** (every single-union construction, 6 doses, 4 folds, priced before selection):
`ARM1` = dld1 only · `ARM3` = H1 pair + dld1 · `ARM4` = 225's four names · `ARM5` = H1 pair +
dld1 + jiang + igvf · `ARM6` = all six.

## What would falsify it — declared here, measured next

1. The exact `nmae` identity (216's `six216.py` / 188's split), which the rate card cannot see.
   **If `nmae` falls materially, the arm is NOT recommended regardless of cosine and `pds`**
   (INBOX-4 §3). Break-even efficiency on 188's split is `e = ΔM/ΔU = 0.5`; 224's arm died at
   0.072.
2. The five package falsifiers (config diff, digest, non-vacuity, carrying separation, prep).

## Outcome

🔴 **Falsifier 1 FIRED.** `six229.py` measured `nmae wald89` **0.98569 → 1.00291** (raw
**+0.01722**, overall **−0.004783**) with efficiency **`e` = −0.033**. The arm was withdrawn.
See `RESULTS.md` §0 and §9. **The pre-registered falsifier did its job on the first arm it was
applied to.**
