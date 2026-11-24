# 187 — The joint amplitude optimum: `mse` and `nmae` want different α, and nobody has solved them together

You are a research subagent. Work ONLY in `experiments/187-joint-amplitude/`.
No `vcc`, no `cell-eval2`, no submission build, no writes outside your directory. RSS ≤ 10 GB.
`uv run python`. This is mostly algebra plus board arithmetic — it should be fast.

## READ FIRST
`AGENTS.md` §4.1 / §4.1b (the magnitude lever and its correction), §11.1 rules 2, 10, 12, 13 with
the audit-181 amendments. Then:
- `docs/CURRENT_STATE.md` — the `2026-09-03 02:15` nmae-landscape entry, the `2026-09-03 00:00`
  `abc-knee` board result, and the composed-slot-gate entry.
- `experiments/170-gamma-mse/RESULTS.md` (the ρ knee), `096-mse-floor/` (the `claim` route),
  `180-def-readiness/BOARD-RESULT.md` (abc-knee), `162-retrodicting-validator/`.

## The situation

Two of the six members are pure functions of the predicted-effect amplitude `α`, and they are
optimised at **different** α. Nobody in this project has solved them jointly; every arm has moved α
to serve one member and paid on the other.

**`nmae`** (exact identity, `docs/CURRENT_STATE.md` 02:15): for a reference DE gene with true
magnitude `|t|`, our magnitude `|p| = k·|t|`, and per-gene sign accuracy `q`,
```
nmae(k, q) = q·|k − 1| + (1 − q)·(k + 1)
```
minimised at **k = 1 exactly**, value `2(1 − q)`, with **asymmetric slopes**: `−(2q − 1) = −0.074`
below the optimum and `+1` above it. Our board `nmae` 0.9432 at `q = 0.537` reads as **k ≈ 0.75** —
we sit on the cheap side of the kink, 25% short of the optimum.

**`mse`** is minimised at the ρ knee `s* = √(claim/B)` (170). `claim(a15 as shipped)` = **0.24581**,
the knee is `s²B` = **0.23268**, so a15 sits **+0.01313 of claim headroom above the knee** and `mse`
is flat above it. The board-validated ρ-cliff algebra holds to **±1%** on a15's global α (audit 181
lists it as SOLID).

**And the local `mse` proxy just failed a board test.** `abc-knee` (entry `N3ui4hCkKx9bxExG2pxE`,
2026-09-03 00:00 UTC) raised B and C toward a locally computed per-context knee — local alphas
0.98586 / 1.01275 / 1.04568 against a15's global 0.982 — and the board returned, as **Δ scaled per
member** (overall = sum/6):

| pds | mse | nmae | fid | jac | reach | overall |
|---:|---:|---:|---:|---:|---:|---:|
| +0.00012 | **−0.00226** | +0.00208 | +0.0006 | +0.0001 | +0.00347 | **+0.00068** |

Raw `mse` 0.87235 → 0.87449. The prereg predicted `mse` **+0.0050**; the board delivered **−0.00226**.
**The `C_p` NTC-jackknife proxy does not transfer.** The per-context knee is therefore NO-GO as a
D/E/F default until a proxy retrodicts this point.

## Where the board says we actually stand

Regressing each raw member on raw `pds` across the 206 board teams with `pds > 0.70` (and, for
`mse`, the 49 shrunk teams with raw `mse < 1.2`), our residuals conditional on our own
discrimination are:

| member | ours | field fit at our pds | residual | in σ | overall points to reach the line |
|---|---:|---:|---:|---:|---:|
| `mse` | 0.8745 | 0.8392 | +0.0353 (worse) | +0.59σ | **+0.0062** |
| `nmae` | 0.9432 | 0.9291 | +0.0141 (worse) | +0.43σ | **+0.0039** |
| `fid` | 0.5038 | 0.5013 | +0.0025 | +0.06σ | — |
| `jac` | 0.0301 | 0.0298 | +0.0003 | +0.08σ | — |
| `reach` | 0.1773 | 0.2206 | −0.0433 | −1.04σ | +0.0082 (185 owns this) |

So **+0.010 overall sits in `mse` + `nmae` alone**, and the field proves it is attainable at our
discrimination. That is a quarter of the 0.0438 gap to first place, in the two members that a single
scalar moves.

## What to do

1. **Build a two-member model `overall(α) = [mse_scaled(α) + nmae_scaled(α)] / 6`** from first
   principles: the ρ-cliff algebra for `mse`, the `nmae(k, q)` identity with `k ∝ α`, and our
   measured `q = 0.537`. State every assumption.
2. **RETRODICT THE BOARD (mandatory, audit-181 rule 1).** Your model must reproduce, as its
   known-answer rows, the measured Δ of **abc-knee vs a23** given abc-knee's actual per-context
   alphas — `mse` −0.00226 and `nmae` +0.00208 scaled. Also retrodict a19's `mse` **−0.01352** and
   a15's raw `mse` **0.88228** on h1-test / 0.87235 on the board. **If your model cannot retrodict
   abc-knee it is the same broken proxy and you must say so and stop**; a model that only fits the
   points it was built from is worth nothing here.
3. **Solve for the joint optimum** — the global α that maximises the sum, and separately whether a
   per-context α helps once `nmae` is in the objective (abc-knee only priced `mse`). Report the
   predicted six-member delta and put it through the composed slot gate:
   `predicted = PLOLO d_pds/6 × 1.011 + 096 d_mse/6 (amplitude arm) − 0.0146/6 (emission-law arm
   only)`; SHIP only if predicted > **+0.0089** for an amplitude-changing arm.
4. **Explain the −0.0141 `nmae` residual mechanistically.** Is our `nmae` deficit entirely `k ≈ 0.75`,
   or is part of it a `q` deficit? The identity separates them: `nmae` at `k = 1` is exactly
   `2(1 − q)`, so from a measured `nmae` and a measured `k` you can solve for the implied `q` and
   compare it to 0.537. Do that for us AND for the leader (0.9269) and `Aginglab.com` (0.8892).
   Which of the two is the field beating us on? That single number decides whether `nmae` is an
   amplitude problem (cheap, ours today) or a direction problem (expensive, 184's territory).
5. Note the constraint from the same 02:15 entry: any magnitude increase must not push DE-pool
   magnitudes **above** their true values, because the wrong-sign 46% are charged `|p| + |t|` there.
   A uniform `k` increase toward 1 is safe on the identity; check whether the distribution of
   per-gene `k` makes it unsafe in practice.

## Deliverable
`RESULTS.md`: the model, the retrodiction table (pass/fail per board point), the joint optimum with
its predicted six-member delta, the gate verdict, and the `k` vs `q` decomposition for us, the
leader and Aginglab. If the answer is "the joint optimum is where a15 already sits", say so — that
closes amplitude for good and is worth knowing.

**Do not build or submit anything.** If an arm clears the gate, write the exact builder invocation
and the preregistered per-member prediction, and stop.

## Territory
184 owns Tahoe / cross-line direction. 185 owns emission and the DE call set. 186 owns intel.
You own amplitude, `mse` and `nmae`.
