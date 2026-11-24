"""308 -- the TARGET-SET-MATCHED controls for the only cell with 3/3 positive folds (`soft:abs:p1.0:dis`).

The controls in sweep308.py are UNRESTRICTED (they re-weight all ~259 two-donor rows) while the `:dis` arm touches
only the 20 / 63 / 49 disagreeing rows, so their comparison is not matched: any perturbation of more rows costs more.
Here every control acts on EXACTLY the arm's own row set (the real disagreeing gate, held fixed) and differs only in
the weight it puts there: |kd| shuffled within the set, |kd| reversed, and entropy-matched random weights.

    PYTHONDONTWRITEBYTECODE=1 uv run python matched308.py <fold>
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
sys.path.insert(0, str(HERE))
import lib308 as Z  # noqa: E402

Y, L, X = Z.Y, Z.Y.L, Z.Y.X
POWER = 1.0
CUT = "abs"


def main(held: str) -> None:
    t0 = time.time()
    st = Y.stage(held)
    ident = X.Identity(st)
    targets, donors, avail = list(st["targets"]), list(st["sources"]), st["avail"]
    feat = Z.load_feature()
    ff = Z.fold_feature(feat, targets, donors)
    A = Z._fill_missing(ff["A"], ff["M"], avail)
    cut = ff["cut_abs"] if CUT == "abs" else ff["cut_rel"]
    two = avail.sum(0) >= 2
    works, fails = (A >= cut[:, None]) & avail, (A < cut[:, None]) & avail
    dis = two & (works.sum(0) >= 1) & (fails.sum(0) >= 1)
    rows, per_target, identity = [], {}, {}

    def run(name, W=None):
        if W is None:
            pre_sel, base = st["pre"], st["whitened"]
        else:
            pre_sel = Y.combine(st, W, norm_match=True)
            base = Y.base_from(st, pre_sel)
        e, m = Y.evaluate(st, base, ident)
        per_target[f"{name}|pds"] = np.asarray(e["_pds_per_target_a23"])
        for k in ("U", "M", "D"):
            identity[f"{name}|{k}"] = m[k]
        rows.append(dict(fold=held, arm=name, n_weighted=int(np.isfinite(W).all(1).sum()) if W is not None else 0,
                         **{k: v for k, v in e.items() if not k.startswith("_")},
                         nmae_wald89=m["nmae_wald89"], top89_precision=m["top89_precision"],
                         delivered_norm=m["delivered_norm"], audit=m["audit"]))
        print(f"[{held}/{name}] cos_agg {e['cos_agg']:.6f} pds {e['pds']:.6f} nmae89 {m['nmae_wald89']:.5f} "
              f"w {rows[-1]['n_weighted']}", flush=True)

    run("cd4svd5")
    W_arm, _ = Z.soft_weights(A, avail, cut, POWER, gate=dis)
    run(f"soft:{CUT}:p{POWER}:dis", W_arm)
    Wr, _ = Z.soft_weights(A, avail, cut, POWER, gate=dis, reverse=True)
    run(f"soft:{CUT}:p{POWER}:dis:rev", Wr)
    for s in Z.SEEDS:
        Ash = Z.shuffle_kd(A, avail, s)                    # |kd| permuted within donor; the GATE stays the real one
        Ws, _ = Z.soft_weights(Ash, avail, cut, POWER, gate=dis)
        run(f"soft:{CUT}:p{POWER}:dis:shuf:s{s}", Ws)
        run(f"soft:{CUT}:p{POWER}:dis:rand:s{s}", Z.entropy_matched_random(W_arm, avail, 1000 + s))
    (WORK / f"matched308-{held}.json").write_text(json.dumps(dict(
        meta=dict(fold=held, n_disagree=int(dis.sum()), cut=CUT, power=POWER, seconds=round(time.time() - t0, 1),
                  peak_gib=Y.peak_gib()), rows=rows), indent=2))
    np.savez_compressed(WORK / f"matchedpt308-{held}.npz", targets=np.asarray(targets), dis=dis.astype(np.int8),
                        **per_target)
    np.savez_compressed(WORK / f"matchedid308-{held}.npz", **identity)
    print(f"[{held}] {len(rows)} rows, {time.time()-t0:.0f}s, peak {Y.peak_gib():.2f} GiB", flush=True)


if __name__ == "__main__":
    main(sys.argv[1])
