"""300 -- the ROW-LEVEL collinearity of the two axes: reliability rank `u` vs top-50 leverage rank `v`.
If arm C is arm A by another name, these two ranks must agree at the row level.

    PYTHONDONTWRITEBYTECODE=1 uv run python uv300.py orion_hct116
"""
from __future__ import annotations
import gc, json, sys
import numpy as np
sys.path.insert(0, ".")
import run300 as RU

def sp(x, y):
    rx = np.argsort(np.argsort(x)).astype(float); ry = np.argsort(np.argsort(y)).astype(float)
    rx -= rx.mean(); ry -= ry.mean()
    return float((rx * ry).sum() / np.sqrt((rx * rx).sum() * (ry * ry).sum()))

def main(held: str) -> None:
    st = RU.stage300(held)
    donors = list(st["sources"]); banks = [st["banks"][nm] for nm in donors]
    st["fit_banks"] = banks; wc = st["wcols"]
    regp = RU.plolo.OPERATOR_REGISTRY[f"held_{held}"]
    with np.load(RU.ROOT / regp["path"]) as z:
        ev = np.asarray(z["evals"], dtype=np.float64)
        top = np.argsort(ev)[-RU.N_PC:][::-1]
        V50 = np.asarray(z["evecs"][:, top], dtype=np.float64)
    out = dict(fold=held, n_pc=RU.N_PC, top_eig_trace_share=float(ev[top].sum() / ev.sum()))
    rel = {nm: RU.load_rel(nm) for nm in donors}
    full = {b.name: np.arange(len(b.perts), dtype=np.int64) for b in banks}
    ess = set(RU.ESSENTIAL.read_text().split("\n")) - {""}
    for b in banks:
        u = RU.tilt_u(b, rel[b.name], full[b.name])
        v, aud = RU.spread_v(st, b, full[b.name], wc, V50)
        ise = np.asarray([p in ess for p in b.perts], bool)
        out[b.name] = dict(n=int(u.size), spearman_u_v=sp(u, v), share=aud,
                           mean_v_essential=float(v[ise].mean()), mean_v_other=float(v[~ise].mean()),
                           mean_u_essential=float(u[ise].mean()), mean_u_other=float(u[~ise].mean()))
        print(f"{b.name:14s} n {u.size:6d}  Spearman(u reliability rank, v leverage rank) {sp(u,v):+.4f}   "
              f"top-{RU.N_PC} energy share median {aud['share_median']:.4f} (p10 {aud['share_p10']:.4f} p90 {aud['share_p90']:.4f})   "
              f"essential rows: mean v {v[ise].mean():.3f} vs {v[~ise].mean():.3f}, mean u {u[ise].mean():.3f} vs {u[~ise].mean():.3f}",
              flush=True)
    del V50; gc.collect()
    (RU.WORK / f"uv300-{held}.json").write_text(json.dumps(out, indent=1, default=float))

if __name__ == "__main__":
    main(sys.argv[1])
