"""214 — arm A's exposure on the CHALLENGE call set, and the delivered-amplitude check.

`rc_callset` prescribes `top 400 by |H| U top 2,200 by control CPM` and `rc_emit.py`:920
ASSERTS the union comes out at exactly 2,600, so no change to the effect matrix can move the
call-set SIZE.  What a geometry arm can move is MEMBERSHIP, through the `|H|` head.  This
measures it on all 900 challenge blocks, the same way 210 falsifier 4 measured the CD4 dose's
exposure (head overlap min 184 / median 314 / max 399 of 400).

It also measures the two amplitude statistics 211's standing rule requires of any arm that
might move `A`: the delivered median row norm and `A = mean ||row||^2`, on the same five
nested gene sets `amp210.py` used, so the shape-only / amplitude-changing classification is
MEASURED rather than argued.

    uv run python head214.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
WORK = HERE / "work"
WORK.mkdir(exist_ok=True)
for p in ("experiments/162-retrodicting-validator", "experiments/200-production-ablation",
          "experiments/196-emitted-depth", "experiments/096-mse-floor", "src"):
    sys.path.insert(0, str(ROOT / p))
sys.path.insert(0, str(HERE))

import chalparts200                                                      # noqa: E402
import mse200                                                            # noqa: E402
import recon196                                                          # noqa: E402
from plolo import PROD_OPS, Operator                                     # noqa: E402
from pds214 import arm                                                   # noqa: E402

TS = 50_000.0
N_HEAD, N_CPM = 400, 2200
CFG = {"inc": arm(cd4_f=0.34, pm=0.50, mass=None),
       "A": arm(cd4_f=0.34, pm=0.00, mass=None),
       "f010_inc": arm(cd4_f=0.10, pm=0.50, mass=None)}


def main() -> None:
    p = chalparts200.parts()
    c4 = chalparts200.cd4()
    op = Operator(PROD_OPS / "whitener_union_pruned.npz")
    E = {k: mse200.geometry_effects(cfg, p, c4, op) for k, cfg in CFG.items()}
    del op
    targets = list(p["targets"])
    pos = {t: i for i, t in enumerate(targets)}
    print(f"challenge panel effect matrices rebuilt: {E['inc'].shape}, "
          f"{len(targets)} targets")

    out: dict = {}
    print("\n=== ARM A's CALL-SET EXPOSURE (rc_callset's prescribed set, per block) ===")
    print(f"  {'ctx':<5}{'head overlap /400':>20}{'call-set jaccard':>20}"
          f"{'prescribed size':>18}")
    rows = []
    for ctx in "ABC":
        D = recon196.load_chal(ctx)
        keep = D["mean_cpm"] > 5.0            # the SCORER's DE universe on this context
        cpm_rank = np.where(keep, D["mean_cpm"], -np.inf)
        rest_order = np.argsort(-cpm_rank, kind="stable")
        ov, jac, sizes = [], [], []
        for t in D["targets"]:
            i = pos[str(t)]
            sets = {}
            for k in ("inc", "A"):
                key = np.where(keep, np.abs(E[k][i]), -np.inf)
                head = np.zeros(key.size, bool)
                head[np.argsort(-key, kind="stable")[:N_HEAD]] = True
                # `rc_emit.CallSet`: the CPM half runs over the genes the head did not
                # take.  The head holds at most N_HEAD genes, so the first
                # N_HEAD + N_CPM by CPM always contain N_CPM non-head genes.
                called = head.copy()
                top = rest_order[:N_HEAD + N_CPM]
                called[top[~head[top]][:N_CPM]] = True
                assert int(called.sum()) == N_HEAD + N_CPM
                sets[k] = (head, called)
            ov.append(int((sets["inc"][0] & sets["A"][0]).sum()))
            u = (sets["inc"][1] | sets["A"][1]).sum()
            jac.append(float((sets["inc"][1] & sets["A"][1]).sum() / u))
            sizes.append((int(sets["inc"][1].sum()), int(sets["A"][1].sum())))
        rows.append(dict(ctx=ctx, head_min=int(min(ov)), head_med=int(np.median(ov)),
                         head_max=int(max(ov)), jac_min=float(min(jac)),
                         jac_med=float(np.median(jac)),
                         size_inc=sorted({s[0] for s in sizes}),
                         size_A=sorted({s[1] for s in sizes})))
        print(f"  {ctx:<5}{f'min {min(ov)} med {int(np.median(ov))} max {max(ov)}':>20}"
              f"{f'min {min(jac):.4f} med {np.median(jac):.4f}':>20}"
              f"{str(sorted({s[0] for s in sizes})) + '/' + str(sorted({s[1] for s in sizes})):>18}")
    out["call_set"] = rows
    print("  (`rc_emit.py`:920 asserts the prescribed union is exactly "
          f"{N_HEAD + N_CPM}; the sizes above confirm it for both arms, so arm A cannot "
          "move the challenge call-set SIZE and 210's measured n_pred 2,599 carries over.)")

    # ------------------------------------------------------ 211's amplitude requirement
    print("\n=== DELIVERED AMPLITUDE (211's standing rule: does arm A move `A`?) ===")
    print("  rendered exactly as `profile_from_effect` does, per context, then read back "
          "the way the scorer does")
    print(f"  {'gene set':<28}{'n':>7}{'median|row| inc':>18}{'arm':>10}{'ratio':>9}"
          f"{'A inc':>11}{'arm':>10}{'ratio':>9}")
    amp = {}
    for ctx in "ABC":
        D = recon196.load_chal(ctx)
        b = D["b"]
        gate = D["gate"]
        idx = np.asarray([pos[str(t)] for t in D["targets"]])
        deliv = {}
        for k in ("inc", "A"):
            H = np.where(gate[None, :], E[k][idx], 0.0)
            lin = np.clip(np.expm1(b[None, :] + H), 0.0, None)
            prof = lin / lin.sum(1, keepdims=True)
            deliv[k] = np.log1p(TS * prof) - b[None, :]
        panel_tot = D["sums"].sum(0)
        sets = {"all genes": np.ones(b.size, bool),
                "builder 10-CPM gate": gate,
                "top 2,000 by control CPM":
                    np.zeros(b.size, bool),
                "panel total >= 1,000 counts": panel_tot >= 1_000,
                "panel total >= 10,000 counts": panel_tot >= 10_000}
        sets["top 2,000 by control CPM"][np.argsort(-D["mean_cpm"])[:2000]] = True
        for name, m in sets.items():
            r = {}
            for k in ("inc", "A"):
                rn = np.linalg.norm(deliv[k][:, m], axis=1)
                r[k] = (float(np.median(rn)), float(np.mean(rn ** 2)))
            amp.setdefault(name, []).append(dict(ctx=ctx, n=int(m.sum()),
                                                 med_inc=r["inc"][0], med_A=r["A"][0],
                                                 A_inc=r["inc"][1], A_A=r["A"][1]))
    for name, lst in amp.items():
        n = "/".join(str(x["n"]) for x in lst)
        mi = float(np.mean([x["med_inc"] for x in lst]))
        ma = float(np.mean([x["med_A"] for x in lst]))
        ai = float(np.mean([x["A_inc"] for x in lst]))
        aa = float(np.mean([x["A_A"] for x in lst]))
        print(f"  {name:<28}{n:>7}{mi:>18.6f}{ma:>10.6f}{ma/mi:>9.5f}"
              f"{ai:>11.6f}{aa:>10.6f}{aa/ai:>9.5f}")
    out["amplitude"] = amp
    json.dump(out, open(WORK / "head214.json", "w"), indent=1, default=float)
    print(f"\nwrote {WORK/'head214.json'}")


if __name__ == "__main__":
    main()
