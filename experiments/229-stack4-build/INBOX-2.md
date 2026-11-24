# INBOX 2 — while you are building the registry, REGISTER four sources you will not ship. October needs them.

**This does not change your arm.** Build exactly what 225 recommended:
`dld1_lowmoi ∪ jiang_5line_noK562 ∪ igvf_cardiomyocyte ∪ myllia_a549` at `g ≈ 0.35`. **Nothing
below alters the package.** It is a registry addition only.

## Why

225 §8.3 projected every source onto `DEF-RUNBOOK` §5's two synthetic October draws (three
known-answer rows against §5's own numbers — `k562_gwps` 269→**150**, Orion **291**/**293** — exact
on all three). **Four sources we already hold carry ZERO of today's panel and ~30 targets each of
October's:**

| source | today /300 | pessimistic October /300 |
|---|---:|---:|
| `nadig_hepg2` | **0** | **33** |
| `nadig_jurkat` | **0** | **33** |
| `rpe1` | **0** | **32** |
| `k562_essential` | **0** | **30** |

**PLOLO cannot price them at all** — it scores today's panel and they carry none of it. Their first
measurable price is on the new panel. If they are not in the registry, we discover them on Oct 22,
under time pressure, with the operator being re-fit and the banks re-extracted at the same time.

For contrast, the two sources that look best by today's overlap — **`kolf2` (273/300)** and
**`feng_ipsc` (182/300)** — project to **EXACTLY ZERO**. They are panel-*fitted* libraries and
October's panel is disjoint by construction.

## What to do

**Add all four to `UNION_SOURCE_REGISTRY` with their paths and shas, marked as not-active.** Same
treatment as any other entry — rule 11, a path and a sha, never a description. They must be
selectable by name via the env var and must **not** be selected in your build.

**Then add one known-answer row proving they are inert**: with your arm's source list selected, the
package hash is unchanged whether or not the four are present in the registry. **Registry membership
must not be able to affect a build** — that is the property that makes it safe to pre-register
sources we cannot yet price, and it is worth one assertion.

If any of the four is not in the standard compact-donor layout (`delta.npy`, `b_ctrl.npy`,
`genes.txt`, `perts.txt`), **say so in `RESULTS.md` with what is missing and do not fabricate it** —
knowing which of the four need a build step is itself the useful output, and it is far better to
learn it now than on the day.
