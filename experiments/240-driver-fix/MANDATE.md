# 240 — FIX the Oct-22 driver. It discards both board-validated arms and the post-build check passes. −0.0199671 on the only prize-eligible entry.

You are a build subagent. Work in `experiments/240-driver-fix/`. **You MAY modify
`scripts/def_pipeline.py` and `scripts/def_postbuild.py`.** You may **NOT** run `vcc submit`, choose
a model name, or write a description. `uv run python`.

**MEMORY GATE:** proceed when `memory_pressure ≥ 40% free` AND `vm_stat` free+inactive **≥ 20 GiB**
AND swap free **≥ 400 MiB**. Heavy jobs strictly serially; `/usr/bin/time -l` for peak RSS. `prep`
peaks at **22.33 GiB** and must never run beside a build. Two other agents are live (221, 238).

## The defect (239 §2.1 — read it first, it is exact)

`docs/DEF-RUNBOOK.md:90` gives **one** command for Oct 22, `scripts/def_pipeline.py`, under
*"Follow it top to bottom. Do not improvise on the day."* That driver contains:

```python
# scripts/def_pipeline.py:58
INCUMBENT_ENV = {"VCC_ALPHA": "0.982", "VCC_CD4_TARGET_F": "0.1"}   # the ONLY env the incumbent sets
```

pinned, by its own comment, to **`build191.sh` — the 2026-09-03 incumbent, two shipped arms ago**.
And it is applied as an **override**, not a default:

```python
env={**os.environ, **(env or {})}      # INCUMBENT_ENV goes LAST and wins
```

**Consequences on the one entry that can win the prize:**

| arm | the one variable | what the driver does | Δ overall lost |
|---|---|---|---:|
| `cd4f034` | `VCC_CD4_TARGET_F` 0.10 → **0.34** | **silently resets to 0.1**, even if you export 0.34 | **−0.0130557** |
| `h1union` | `VCC_H1_UNION_G` = **0.28** | **never set**; would pass through only on top of `f = 0.10`, a combination never built or scored | **−0.0069114** |
| | | | **−0.0199671** |

That is **board 0.2118310 → 0.1918638**, and **larger than our entire gap to first place (+0.01248)**.
**And `def_postbuild.py` passes**, because it diffs against **`rc_callset`'s** meta — the incumbent
from before both arms.

## The fix

1. **`INCUMBENT_ENV` must become a DEFAULT, not an override.** Use `{**INCUMBENT_ENV, **os.environ,
   **(env or {})}` — or an explicit `setdefault` — so an operator-supplied value **wins**. Then
   re-point the defaults at the **currently shipped** configuration: `VCC_ALPHA` **0.982**,
   `VCC_CD4_TARGET_F` **0.34**, `VCC_H1_UNION_G` **0.28**, and the registry's default source
   selection (the H1 pair) as landed by 229 (`0204168`) and hardened by 235 (`b540148`).
2. **`def_postbuild.py` must diff against the CURRENT incumbent**, not `rc_callset`. Make the
   reference explicit and **fail loudly** if the built config does not carry the arms the runbook
   says we ship.
3. **Add a guard that makes this class of failure impossible to repeat**: the driver must **assert**
   that the config it produced contains the arms it intended, and **raise** otherwise — the same
   discipline 235 applied to the operator axis. **An assertion never observed to fire is not a
   guard: demonstrate it firing.**

## Known-answer rows — all must pass
1. **Driver defaults, on the A/B/C bundle, reproduce the shipped `h1union`**: `.vcc` sha256
   `f13ce40bca908af7cb85b8e84887af473aae4588410c238529139b7d1f245b1e`, 2,824,785,920 bytes. **This
   is the row that proves the fix.** (If a full build is too expensive, reproduce the **config
   block** exactly and say so explicitly — `effects.cd4.fraction` **0.32207201067578145**,
   `target_f_requested` **0.34**, `f_per_unit_fraction` **1.0556645369046553**,
   `fully_supported_targets` **276**, `operator.sha256` `e2f65ecf…`, γ 0.65, λ 0.0.)
2. **An operator-supplied override now WINS**: `VCC_CD4_TARGET_F=0.10 … def_pipeline.py` must
   produce `f = 0.10`, not be silently overridden. Show the before/after.
3. **The old behaviour is reproduced once, to document what was at risk**: with the pre-240 driver,
   show `f` forced to 0.1 and `VCC_H1_UNION_G` absent, and **`def_postbuild.py` PASSING** on it.
4. **The new assertion FIRES** on a config missing either arm.
5. **`solve_cd4_fraction` untouched.**

## Then update the runbook
`docs/DEF-RUNBOOK.md` §2's step table and §3's decision table both quote the old flags. **Bring them
to the shipped configuration**, and add a line at the top of §2 naming the two arms the day must
carry and their board deltas, so a future reader can check the driver against the record in one
glance.

## Discipline
- **Do not change any shipped number**, and do not re-tune anything. This is a wiring fix.
- **Commit both scripts and the runbook** when the rows pass, staging those paths explicitly —
  never `git add -A`, two other agents are writing.
- Deliverable: `RESULTS.md` with the five rows (row 3 quoted verbatim — the old driver passing while
  wrong), the diff summary, and a one-line **FIXED / NOT FIXED** verdict.
