# Spec for the next source-ladder build — read before building anything for the 2026-09-05 slots

Written 2026-09-04 04:40 UTC by the coordinator, from a read of
`scripts/make_submission_pm_cd4.py` §`VCC_H1_UNION_G`. **No code was changed.** This is a
readiness note so the next builder does not discover the gap at build time.

## 1. The good news: the machinery is already general

`VCC_H1_UNION_G`'s implementation is **not** H1-specific. `prepare_h1_union` takes a list of
`(name, path)` sources in a standard compact-donor layout — `delta.npy`, `b_ctrl.npy`, `genes.txt`,
`perts.txt` — unions them, sizes the result on its **carrying rows**
(`median(‖added_i‖ : ‖added_i‖ > 0)`), merges it with the CD4 residual into **one** add, and runs
production's single global preserve rescale. That construction is 218 §7.2 exactly, it passed five
falsifiers, and it produced our shipped entry (`h1union`, board **0.2118310, rank 3**).

**Every part of that is source-agnostic.**

## 2. The gap: exactly one hardcoded tuple

```python
H1_UNION_SOURCES: tuple[tuple[str, Path], ...] = (
    ("gse295214_h1", ROOT / "experiments/115-gse295214-h1/work/compact-donor"),
    ("h1_2025",      ROOT / "data/donors/h1_2025"),
)
```

To ship any source 225 finds — `sunshine_calu3` and `dld1_lowmoi` are the current front-runners —
**this list must become configurable.** Nothing else needs to change.

## 3. How to do it, and the known-answer row that governs it

- Add a **registry-driven** selector (a name→path table plus an env var naming which entries are
  active), or extend the env var to carry `name:path` pairs. **Prefer a registry file** — rule 11
  says a path and a sha, never a description, and a registry keeps the shas in one auditable place.
- **`VCC_H1_UNION_G` unset must remain bit-for-bit inert**, and with the H1 pair selected the build
  must reproduce **`f13ce40bca908af7cb85b8e84887af473aae4588410c238529139b7d1f245b1e`**
  (the shipped `.vcc`, 2,824,785,920 bytes) — or, on the `.h5ad`,
  `cd4f034`'s `7036e89f8a03c82d69e79a5c6690b275adae52dbb3ff72bb6ab23e7ca622bf7a` with the flag
  unset. **Those two are the known-answer rows. A build that cannot reproduce them does not ship.**
- **Do not touch `solve_cd4_fraction`.** CD4 sits at 93–98% coverage, its own solve is deliberately
  untouched by this flag (216/217's territory), and 222 proved the fraction did not move by one bit.

## 4. Sizing the new arm

- **One variable per slot.** If the new source rides *alongside* the H1 pair, the H1 dose `g` stays
  at **0.28** and only the new source's presence changes.
- **`g` for the new source comes from 225's ladder on the three cancer folds**, with `kolf2`
  reported separately — it is pluripotent, gains 3–5×, and inflated a prediction this morning.
- **Do not convert PLOLO to board points with a scalar.** The factor is **not identified** (1.384 on
  `cd4f034`, 1.000 on `h1union`, matched functional). Quote PLOLO units, name the fold set and the
  functional. 221 owns the reconciliation.
- Carry the **leave-top-3-out retention** on cosine and `pds` separately; `h1union` realised ~100%
  of its cosine prediction and **70%** of its `pds` one.

## 5. Operational

`vcc prep` peaks at **22.33 GiB RSS** — bigger than the build that feeds it (222 measured it; it is
recorded nowhere else). **Never run it alongside a build**, and check what else is running first.

---

# 6. THE NEXT PRODUCTION EDIT AFTER 229 — per-group doses. Worth +0.00099 PLOLO, and it compounds.

Added 2026-09-04 05:27 from 229's §11. **Do not start this until 229 has finished and committed**
— it owns `scripts/make_submission_pm_cd4.py` and a second writer would corrupt both.

## The gap 229 measured

The builder carries **one** extra residual at **one** dose. 225 priced its stack as a **second**
residual riding on top of the shipped one: `cd4@0.34 + h1pair@0.28`, **then** `+ stack4@0.35`.
Naming all six sources in one union at one dose is a *different construction*, and it is worth less:

| construction | buildable today | PLOLO vs shipped `h1union`, 3 cancer folds |
|---|---|---:|
| **`CUM4`** — h1pair@0.28 **then** stack4@0.35 (two residuals, two doses) | **no** | **+0.00367** (1.93× gate) |
| **`ARM6`** — all six in one union at one dose | yes | **+0.00268** (1.41× gate) |
| `ARM4` — 225's four names, H1 pair dropped | yes | **−0.00102** |

**The buildable form leaves +0.00099 PLOLO on the table**, ~37% of the arm.

## The change, and why it stays a ONE-VARIABLE arm

Extend `UNION_SOURCE_REGISTRY` from a flat list to **named groups, each with its own dose**, applied
as successive post-whitening residuals before the single global preserve rescale.

**Fix the H1 group at its shipped `g = 0.28`** — it is on the board at 0.2118310 and is not the
thing under test — **and let only the new group's dose vary.** That is **one** new variable, not two,
and it keeps `AGENTS.md`'s one-variable-per-slot discipline intact. **Do not re-tune the H1 dose in
the same arm**; if it needs re-tuning that is a separate slot.

## Known-answer rows that govern it

1. **No groups declared** → byte-identical to `cd4f034`
   (`7036e89f8a03c82d69e79a5c6690b275adae52dbb3ff72bb6ab23e7ca622bf7a`).
2. **H1 group alone at 0.28** → byte-identical to the shipped `h1union`
   (`f13ce40bca908af7cb85b8e84887af473aae4588410c238529139b7d1f245b1e`, 2,824,785,920 bytes).
3. **A single group containing all six at one dose** → reproduces 229's `ARM6` package exactly.
   That is the row that proves the generalisation is a strict superset of what 229 shipped.
4. **`solve_cd4_fraction` untouched**, fraction `0.32207201067578145` unmoved to the last bit.

## Two cautions from 229's own run

- **Order is part of the arm.** `union_residual` accumulates a float64 sum per coordinate, so a
  permutation of the same members can differ in the last bit. **Record group order and
  within-group order in `config`.**
- **Keep `verify_union_source_hashes()`.** Three of the six banks live under `experiments/*/work/`,
  a tree other experiments write to. 229 added a pre-build re-hash of every selected bank's four
  files that **raises** on mismatch, before a byte of a 25-minute build is spent. **That guard must
  survive the refactor** — it is the only thing standing between a silently-edited bank and a
  corrupted submission.
