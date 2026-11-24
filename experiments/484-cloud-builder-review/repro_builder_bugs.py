"""Reproduce production-builder correctness bugs without challenge data.

    uv run python experiments/484-cloud-builder-review/repro_builder_bugs.py

Prints one section per finding. Exits 0 if the q80 shrink reproduces, 1 otherwise.
The q80 section execs the production lines in src/vcc2026/rc_emit_k.py (the prior
block), then optionally runs rc_hybrid_block_k on a tiny control. It does not
edit that file.
"""
from __future__ import annotations

import json
import textwrap
import traceback
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src/vcc2026/rc_emit_k.py"
# 1-based line range of the prior target / shrink block (A0 through `move =`).
PRIOR_LO, PRIOR_HI = 941, 982


def section(title: str) -> None:
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def production_prior_lines() -> str:
    lines = SRC.read_text().splitlines()
    chunk = lines[PRIOR_LO - 1 : PRIOR_HI]
    text = textwrap.dedent("\n".join(chunk))
    if "side = blk & (tgt < 0)" not in text:
        raise SystemExit(f"{SRC} lines {PRIOR_LO}-{PRIOR_HI} no longer contain the shrink mask")
    if "the strongest keep their delivered LFC" not in text:
        raise SystemExit("protection comment moved; re-pin PRIOR_LO/PRIOR_HI")
    return text


def run_prior_math() -> dict:
    """Exec the production prior block on a hand-built gated axis."""
    n_keep = 24
    # Two head genes (largest |key|, low control CPM so they are not the block).
    # Then high-CPM genes: two strong down, two strong up, and weak ups the prior flips.
    # Magnitudes are chosen so the mass residual changes sign inside (0, 1), which is
    # the production regime (recorded kf median ~0.60, never 1).
    ref_lin = np.array(
        [30.0, 40.0] + [2_000.0] * 2 + [800.0] * 2 + [400.0] * 16 + [20.0, 15.0],
        dtype=np.float64,
    )
    assert ref_lin.size == n_keep
    e0 = np.array(
        [0.2, -0.2]
        + [-1.20, -1.05, 0.90, 0.80]
        + [0.12] * 16
        + [0.05, -0.05],
        dtype=np.float64,
    )
    grand = 48 * 12_000
    a0 = ref_lin * np.power(2.0, e0)
    ms = np.maximum(1, np.rint(a0 * grand / 1e6)).astype(np.int64)
    # Recompute E0 from the integer totals the emitter actually sees.
    key_gated = np.zeros(n_keep)
    key_gated[0] = 5.0
    key_gated[1] = 4.0
    key_gated[2:] = 0.3

    class CS:
        n_head = 2
        n_cpm = 18

    ns = {
        "np": np,
        "n_keep": n_keep,
        "Ms": ms,
        "grand": grand,
        "ref_lin": ref_lin,
        "key_gated": key_gated,
        "call_set": CS(),
        "lfc_target": {"prior": -np.ones(n_keep), "shrink_floor": 0.002, "protect_q": 0.8},
    }
    exec(production_prior_lines(), ns)  # noqa: S102 — the production source, not user input
    e0v = ns["E0v"]
    blk = ns["blk"]
    thr = float(np.quantile(np.abs(e0v[blk]), 0.8))
    protected = blk & (np.abs(e0v) > thr)
    side = ns["side"]
    kf = float(ns["kf"])
    t2 = ns["_t2"](kf)
    print(f"source: {SRC}:{PRIOR_LO}-{PRIOR_HI}")
    print(f"R(0)={ns['_R'](0.0):.4f}  R(1)={ns['_R'](1.0):.4f}")
    print(f"kf={kf:.6f}  block={int(blk.sum())}  protected={int(protected.sum())}  thr={thr:.4f}")
    print(f"protected & side (down genes the shrink still multiplies): {int((protected & side).sum())}")
    rows = []
    for i in np.flatnonzero(protected):
        kind = "DOWN" if e0v[i] < 0 else "UP"
        ratio = t2[i] / e0v[i] if e0v[i] != 0 else float("nan")
        in_side = bool(side[i])
        print(
            f"  gene {i:2d} {kind}  E0={e0v[i]:+.4f}  delivered_target={t2[i]:+.4f}  "
            f"ratio={ratio:.4f}  in_side={in_side}"
        )
        rows.append((kind, ratio, in_side, float(e0v[i]), float(t2[i])))
    down = [r for r in rows if r[0] == "DOWN"]
    up = [r for r in rows if r[0] == "UP"]
    if not down or not up:
        raise SystemExit("synthetic did not put both a down and an up gene in the protected set")
    # Production bug: protected DOWN ratios track kf; protected UP ratios stay 1.
    # kf in (0, 1): `_t2` returns tgt*kf, and tgt was restored to E0, so the ratio is kf.
    # The floor only binds when |E0|*kf < shrink_floor.
    down_bad = all(abs(r[1] - kf) < 1e-6 and r[2] and abs(r[3] - r[4]) > 0.05 for r in down)
    up_ok = all(abs(r[1] - 1.0) < 1e-9 and not r[2] for r in up)
    print(f"protected DOWN scaled by kf (not left at E0): {down_bad}")
    print(f"protected UP left at E0:                     {up_ok}")
    if not (down_bad and up_ok and 0.05 < kf < 0.99):
        raise SystemExit("q80 shrink did not reproduce")
    return {"kf": kf, "n_protected_down_in_side": int((protected & side).sum())}


def audit_shrink() -> None:
    section("F1b. prior_shrink_k on recorded clean4q80 rc audits (if present)")
    paths = [
        ROOT / "experiments/432-def-rehearsal-new-knobs/work/b2/rc_audit_b2.jsonl",
        ROOT / "experiments/432-def-rehearsal-new-knobs/work/b0/rc_audit_b0.jsonl",
        ROOT / "experiments/429-zs-production/work/rc_audit_429a.jsonl",
    ]
    for path in paths:
        if not path.exists():
            print(f"missing {path}")
            continue
        ks = []
        ctx = {}
        for line in path.open():
            d = json.loads(line)
            if "prior_shrink_k" not in d:
                print(f"{path.name}: no prior_shrink_k")
                break
            ks.append(d["prior_shrink_k"])
            ctx[d.get("context", "?")] = ctx.get(d.get("context", "?"), 0) + 1
        else:
            a = np.asarray(ks, dtype=np.float64)
            print(
                f"{path.relative_to(ROOT)}  n={a.size} contexts={ctx}\n"
                f"  min/p10/p50/p90/max {a.min():.4f} {np.quantile(a, .1):.4f} "
                f"{np.median(a):.4f} {np.quantile(a, .9):.4f} {a.max():.4f}\n"
                f"  fraction kf<0.99 {float((a < 0.99).mean()):.3f}  kf<0.5 {float((a < 0.5).mean()):.3f}"
            )


def try_full_emitter() -> None:
    section("F1c. end-to-end rc_hybrid_block_k on a 48 x 28 toy (optional)")
    try:
        from vcc2026.donor_transfer import emit_counts
        from vcc2026.rc_emit import CallSet, control_block_from_counts, de_from_block
        from vcc2026.rc_emit_k import LAST_PRIOR_REC, rc_hybrid_block_k
    except Exception:
        traceback.print_exc()
        print("SKIP full emitter (import failed)")
        return
    rng = np.random.default_rng(0)
    n_ctrl, n, g = 36, 48, 28
    # Control: genes 0..21 expressed, libraries vary so mean-CPM is well above 5.
    ctrl = np.zeros((n_ctrl, g), dtype=np.int64)
    for i in range(n_ctrl):
        lib_target = 8_000 + (i % 5) * 1_500
        weights = np.array([80, 70, 40, 30] + [12] * 18 + [0] * 6, dtype=np.float64)
        weights = weights / weights.sum()
        ctrl[i] = rng.multinomial(lib_target, weights)
    control = control_block_from_counts(ctrl, min_cpm=5.0)
    n_keep = int(control.keep.sum())
    print(f"DE universe {n_keep}/{g}")
    if n_keep < 16:
        print("SKIP full emitter (gate kept too few genes)")
        return
    # Profile: suppress two high-CPM genes, boost the weak ones, so E0 spans signs.
    ref = np.zeros(g)
    ref[control.keep] = control.ref_lin
    order = np.argsort(-ref)
    profile = np.maximum(ref, 0.0)
    profile[order[2]] *= 0.15   # strong down, high CPM, not the |key| head
    profile[order[3]] *= 0.20
    profile[order[4]] *= 3.0    # strong up
    profile[order[5]] *= 2.5
    profile[order[8:20]] *= 1.4
    profile = profile / profile.sum()
    key = np.zeros(g)
    key[order[0]] = 4.0
    key[order[1]] = 3.5
    key[order[2:]] = 0.2
    libs = np.full(n, 12_000.0)
    n_head, n_cpm = 2, min(14, n_keep - 2)
    LAST_PRIOR_REC.clear()
    try:
        block, diag = rc_hybrid_block_k(
            profile, libs, 0.0, np.random.default_rng(1),
            spread=0.0, overdispersion=0.0, control=control, key=key, jitter_seed=7,
            call_set=CallSet(n_head, n_cpm, 1.90),
            classes=[(1 / 3, 8_000.0), (1 / 3, 14_000.0), (1 / 3, 20_000.0)],
            lfc_target={"prior": -np.ones(n_keep), "shrink_floor": 0.002, "protect_q": 0.8},
            emit_counts=emit_counts,
        )
    except Exception as exc:
        print(f"full emitter raised {type(exc).__name__}: {exc}")
        if LAST_PRIOR_REC:
            print("prior record before the failure:", LAST_PRIOR_REC[-1])
        return
    rec = LAST_PRIOR_REC[-1] if LAST_PRIOR_REC else {}
    _z, _p, _q, lfc = de_from_block(block, control)
    print(f"prior record: {rec}")
    print(f"steer_max_log2_err={diag.get('steer_max_log2_err')} dose={diag.get('dose_delivered')}")
    # Rebuild E0 from the placed totals (bulk CPM vs ref) on the gated axis.
    totals = block.sum(0)[control.keep]
    s = float(block.sum())
    a0 = 1e6 * totals / s
    with np.errstate(divide="ignore", invalid="ignore"):
        e0 = np.log2(np.maximum(a0, 1e-300) / np.maximum(control.ref_lin, 1e-300))
    print(f"delivered lfc on gated axis: min {lfc.min():.3f} max {lfc.max():.3f}")
    print("(full emitter returned; compare prior_shrink_k with F1. kf<1 means protected downs moved.)")
    print(f"example gated E0 vs delivered lfc (first 8):")
    for i in range(min(8, e0.size)):
        print(f"  i={i} E0={e0[i]:+.3f} lfc={lfc[i]:+.3f}")


def prior_gate() -> None:
    section("F2. prior-gate strict < drops the tie at the quantile (numpy, same as the builder)")
    # 468 RESULTS: A/B/C threshold 2, weak 86/83/83 of 300, quantiles (0, 1, 3, 7, ~450).
    # Reconstruct a strength vector with that shape and show strength==2 is excluded,
    # and that a threshold of 0 gates nobody.
    def gate(strength, q):
        thr = float(np.quantile(strength, q))
        weak = strength < thr
        return thr, int(weak.sum()), int((strength == thr).sum())

    # Piecewise stand-in matching the reported quantiles closely enough to show the tie rule.
    # Real A: 86 weak means 86 values < 2. Use 50 zeros and 36 ones, then a mass at 2, then the tail.
    s = np.array([0] * 50 + [1] * 36 + [2] * 20 + [3] * 40 + [7] * 80 + [40] * 74, dtype=np.float64)
    # Pad/trim to 300.
    if s.size < 300:
        s = np.concatenate([s, np.full(300 - s.size, 40.0)])
    s = s[:300]
    thr, n_weak, n_tie = gate(s, 0.33)
    print(f"stand-in n=300 q=0.33 threshold={thr:g} weak={n_weak} tied_at_threshold={n_tie}")
    print(f"  weakest 33% would be {round(0.33 * 300)} targets; strict < keeps {n_weak}")
    zeros = np.zeros(300)
    thr0, n0, _ = gate(zeros, 0.33)
    print(f"all-zero strengths: threshold={thr0:g} weak={n0} (gate is a silent no-op)")
    # Exact virtual index the builder uses.
    print(f"numpy linear index for q=0.33, n=300: {0.33 * 299:.4f}")
    # 468's published cut, restated: threshold 2 and `<` means strength <= 1.
    published = np.array([0, 1, 3, 7, 466], dtype=np.float64)
    print("468 published A quantiles (0,25,50,75,100)%:", published.tolist(),
          "; code keeps strength < threshold, so a threshold of 2 keeps strength <= 1")


def zero_repair() -> None:
    section("F3. zero_repair vs the scorer's 400-cell mean CPM")
    from vcc2026.zero_repair import zero_repair

    n_cells, n_car, g = 400, 384, 30  # k=16 ballast
    # Genes 0..4 gated and zero; gene 0 is moderate CPM, gene 1 is huge (cap binds).
    ctrl_cpm = np.array([8.0, 6.0, 80.0, 400.0, 8000.0])
    gated = np.zeros(g, dtype=bool)
    gated[:5] = True
    block = np.zeros((n_cells, g), dtype=np.int64)
    # Carriers: library ~20_000 from a sub-gate gene; ballast the same.
    block[:, 10] = 20_000
    rng = np.random.default_rng(0)
    out, diag = zero_repair(block, n_car, gated, ctrl_cpm, target_col=-1, rng=rng, floor_log2=-5.0)
    print("diag", diag)
    libs = out.sum(1).astype(np.float64)
    print(f"{'gene':>6} {'ctrl':>10} {'n_counts':>8} {'target_cpm':>12} {'scorer_cpm':>12} {'lfc':>8} {'vs_-5':>8}")
    for j, cpm in enumerate(ctrl_cpm):
        n_counts = int(out[:, j].sum())
        # Scorer arithmetic: mean over ALL cells of count / (lib/1e6). Zeros contribute 0.
        mean_cpm = float(np.mean(out[:, j] / (libs / 1e6)))
        lfc = float(np.log2((mean_cpm + 1e-9) / (cpm + 1e-9)))
        target = cpm * 2.0 ** -5
        print(f"{j:6d} {cpm:10.1f} {n_counts:8d} {target:12.4f} {mean_cpm:12.4f} {lfc:8.3f} {lfc - (-5):8.3f}")
    # Under-repair when the sub-gate cannot pay: gene listed, counts stay 0.
    bare = np.zeros((20, 4), dtype=np.int64)
    bare_gated = np.array([True, True, False, False])
    # Sub-gate columns 2,3 are empty, so every payment is skipped.
    out2, diag2 = zero_repair(
        bare, n_car=20, gated=bare_gated, ctrl_cpm_gated=np.array([10.0, 10.0]),
        target_col=-1, rng=np.random.default_rng(1), floor_log2=-5.0,
    )
    print("empty sub-gate diag", diag2, "repaired col still zero", int(out2[:, 0].sum()), int(out2[:, 1].sum()))


def cpm_gates() -> None:
    section("F4. pseudobulk CPM gate vs mean-per-cell CPM gate, both at 5")
    rng = np.random.default_rng(1)
    n, g = 200, 400
    lib = rng.integers(4_000, 40_000, size=n)
    # Gene 0: expressed only in the small-library half. Gene 1: only in the large half.
    counts = np.zeros((n, g), dtype=np.float64)
    small = lib < np.median(lib)
    counts[small, 0] = 80
    counts[~small, 1] = 80
    # Background so other genes sit near the threshold with no library correlation.
    counts[:, 2:] = rng.poisson(0.05, size=(n, g - 2))
    mean_cpm = (counts * (1e6 / lib)[:, None]).mean(0)
    pseudo = 1e6 * counts.sum(0) / lib.sum()
    both = (mean_cpm > 5) & (pseudo > 5)
    only_mean = (mean_cpm > 5) & ~(pseudo > 5)
    only_pseudo = (pseudo > 5) & ~(mean_cpm > 5)
    print(f"genes mean>5 and pseudo>5: {int(both.sum())}")
    print(f"scorer-only (mean>5, pseudo<=5): {int(only_mean.sum())}")
    print(f"builder-only (pseudo>5, mean<=5): {int(only_pseudo.sum())}")
    if only_mean.any():
        i = int(np.flatnonzero(only_mean)[0])
        print(f"  scorer-only example gene {i}: mean={mean_cpm[i]:.3f} pseudo={pseudo[i]:.3f}")
    print("gene0/gene1 are the correlated extremes and land above 5 on BOTH gates in this draw;")
    print("the disagreements are sparse genes whose two CPM definitions straddle 5.")
    print("a scorer-only gene has key 0 in the emitter (emission_rank_axis) and sign forced to +1")


def ballast_cap() -> None:
    section("F5. ballast remainder can push a cell over 1_000_000 (crash, not a silent count)")
    k, g = 16, 18_533
    # Adversarial: every gene whose index % k == 0 has remainder 1; all other remainders 0.
    # Cell 0 is the only cell that receives those +1s.
    left = np.zeros(g, dtype=np.int64)
    left[np.arange(g) % k == 0] = 1
    base = np.full(g, 800, dtype=np.int64)  # irrelevant; libraries come from q = per*k + left
    q = base * 0 + left  # only the remainder mass, then shift up to the trim cap
    # Put the mass at the trim cap: q.sum() == k * 999_000, remainders as above.
    budget = k * 999_000
    per_gene = (budget - int(left.sum())) // g
    q = np.full(g, per_gene, dtype=np.int64) + left
    # Fix the sum exactly.
    q[0] += budget - int(q.sum())
    # Recompute remainder the way ballast_rows does, from this q. q[0] adjustment may change left.
    per = q // k
    left = q - per * k
    ball = np.tile(per, (k, 1))
    off = np.arange(q.size) % k
    for row in range(k):
        ball[row] += (left > ((row - off) % k)).astype(np.int64)
    libs = ball.sum(1)
    print(f"adversarial libs min/max {int(libs.min())} {int(libs.max())}  over_cap {(libs > 1_000_000).sum()}")
    # Random remainders at the same trim cap.
    rng = np.random.default_rng(0)
    worst = 0
    for _ in range(20):
        raw = rng.random(g)
        raw = raw / raw.sum() * budget
        qq = np.floor(raw).astype(np.int64)
        rem = budget - int(qq.sum())
        qq[np.argsort(-(raw - qq))[:rem]] += 1
        per = qq // k
        left = qq - per * k
        ball = np.tile(per, (k, 1))
        off = np.arange(g) % k
        for row in range(k):
            ball[row] += (left > ((row - off) % k)).astype(np.int64)
        worst = max(worst, int(ball.sum(1).max()))
    print(f"max library over 20 random compositions at the 999_000 trim: {worst}")


def target_in_callset() -> None:
    section("F6. the target's own gene can occupy a call-set head slot")
    # The emitter ranks |key| over the whole DE universe. The scorer drops the target
    # gene from every DE metric (metrics brief, Gene restrictions). One head slot spent
    # on the target is one real gene that does not enter the prescribed set.
    n = 50
    key = np.linspace(0.01, 1.0, n)
    target = 7
    key[target] = 5.0  # the knockdown is the largest |H|
    head_n, cpm_n = 5, 10
    head = np.zeros(n, dtype=bool)
    head[np.argsort(-np.abs(key), kind="stable")[:head_n]] = True
    ref = np.linspace(10, 100, n)
    cpm_rank = np.where(head, -np.inf, ref)
    called = head.copy()
    called[np.argsort(-cpm_rank, kind="stable")[:cpm_n]] = True
    real = called.copy()
    real[target] = False
    print(f"prescribed {int(called.sum())} includes target: {bool(called[target])}")
    print(f"genes the DE metrics can see: {int(real.sum())} (one short of the prescribed size)")


def main() -> None:
    section("F1. q80 does not keep delivered LFC on protected DOWN genes")
    run_prior_math()
    audit_shrink()
    try_full_emitter()
    prior_gate()
    zero_repair()
    cpm_gates()
    ballast_cap()
    target_in_callset()
    print("\nALL REQUESTED REPRODUCTIONS FINISHED")


if __name__ == "__main__":
    main()
