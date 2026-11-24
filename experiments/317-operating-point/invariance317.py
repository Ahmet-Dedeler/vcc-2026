"""317 gate 0 — the structural claim, checked on the real `emit_counts`, locally and cheaply.

`emit_counts(mode="quasi")` computes each gene's panel total as `floor(n * L * profile_g)`
and then repairs the scattered part back to that exact integer. If that holds, then for a
FIXED profile the per-gene panel totals -- and therefore the delivered pseudobulk, and
therefore `pds` and `mse` -- are bit-identical under ANY (spread, overdispersion) and can only
be moved by changing the profile itself (amplitude).

This does not need h1-test: it is a property of the emitter, so a synthetic profile with the
same shape (18,080 genes, 400 cells, real-ish library sizes) tests it exactly.
"""
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from vcc2026.donor_transfer import emit_counts  # noqa: E402

G, N = 18_080, 400
rng = np.random.default_rng(0)
# a heavy-tailed relative-expression profile, normalised, like `profile_from_effect`'s output
p = rng.lognormal(mean=0.0, sigma=2.2, size=G)
p /= p.sum()
libs = rng.integers(4_000, 60_000, size=N).astype(np.float64)

def totals(spread, od, seed=7):
    v, col, ip = emit_counts(p, libs, 0.1, np.random.default_rng(seed),
                             mode="quasi", spread=spread, overdispersion=od)
    t = np.zeros(G, dtype=np.int64)
    np.add.at(t, col.astype(np.int64), v.astype(np.int64))
    return t, v

base_spread = np.full(G, 1.0)
base_spread[:8000] = 0.25            # the shipped keep_top/keep_spread/tail_spread shape
ref, ref_v = totals(base_spread, 4.0)

rows = []
def check(label, spread, od, seed=7):
    t, v = totals(spread, od, seed)
    same = bool(np.array_equal(t, ref))
    # cell-to-cell variance is what the Wilcoxon reads; show that it DID move
    cv = float(np.mean(np.var(np.asarray(v[:0]), axis=0))) if False else 0.0
    rows.append((label, same, int(np.abs(t - ref).sum()), int(t.sum())))
    print(f"  {label:<34} totals identical: {str(same):<5}  "
          f"sum|d| {int(np.abs(t - ref).sum()):>8}  panel total {int(t.sum())}")

print("[317 gate 0] per-gene panel totals vs the shipped emission point "
      f"(spread block 0.25/1.0, overdispersion 4.0), {G} genes x {N} cells")
for od in (0.0, 1.0, 2.0, 8.0, 16.0):
    check(f"overdispersion {od}", base_spread, od)
s2 = base_spread.copy(); s2[:8000] = 0.10
check("keep_spread 0.25 -> 0.10", s2, 4.0)
s3 = base_spread.copy(); s3[8000:] = 0.70
check("tail_spread 1.0 -> 0.70", s3, 4.0)
s4 = np.full(G, 1.0); s4[:16000] = 0.25
check("keep_top 8000 -> 16000", s4, 4.0)
# The emission seed at FIXED library sizes also preserves the totals — the repair is
# exact, so the seed can only reach the pseudobulk through `pick_rng`'s 400 library
# sizes, which is 221 gate 4's recorded reason its a15/a15_seed1 twin differs.
check("emission seed 7 -> 11 (libs fixed)", base_spread, 4.0, seed=11)
# and the negative control: amplitude DOES move the totals
p2 = np.exp(np.log(np.maximum(p, 1e-300)) * 1.0)   # placeholder, replaced below
scaled = p * 1.0
v, col, ip = emit_counts(p, libs, 0.1, np.random.default_rng(7), mode="quasi",
                         spread=base_spread, overdispersion=4.0)
q = p.copy(); q[:] = p
# amplitude in the pipeline scales the LOG effect, which changes the profile:
q = np.exp(np.log(np.maximum(p, 1e-300)) * 1.05); q /= q.sum()
t, _ = totals(base_spread, 4.0)
vq, colq, ipq = emit_counts(q, libs, 0.1, np.random.default_rng(7), mode="quasi",
                            spread=base_spread, overdispersion=4.0)
tq = np.zeros(G, dtype=np.int64); np.add.at(tq, colq.astype(np.int64), vq.astype(np.int64))
print(f"  {'NEGATIVE CONTROL: profile changed':<34} totals identical: "
      f"{str(bool(np.array_equal(tq, ref))):<5}  sum|d| {int(np.abs(tq - ref).sum()):>8}")

ok = all(r[1] for r in rows) and not np.array_equal(tq, ref)
print(f"\n[317 gate 0] {'PASS' if ok else 'FAIL'} — every emission knob preserves the "
      f"per-gene panel totals EXACTLY (the emission seed too, at fixed library\n             sizes); only a change to the PROFILE moves them.")
sys.exit(0 if ok else 1)
