"""487 decision 1e: does the per-target prior (`pt`: a covered target's block prior = the same-line bank row's own sign)
beat the shipped CLEAN4 programme sign when the bank is one of our donors and the truth is another lab's screen?

450 §2's H1 evidence (Pan's row signs 0.764 vs CLEAN4 0.578 on gate genes) covers 30 strong targets. On a genome-wide
in-mix bank `pt` applies to every panel target, the weak ones included, where 467/470 found weak responses lean down
and want CLEAN4 (84 % down), not a line's own signs.
Truths: VIPerturb K562 (Satija; bank K562-GWPS, Weissman) and XuCao HEK293 (bank Orion HEK293T, Xaira; 408 cells/pert).
Per target (own gene excluded): the truth's top-100 genes by |delta| among genes all sources measure; sign accuracy of
  pt      = sign(bank row) where the bank's entry is non-zero (else CLEAN4's sign, as the builder does),
  clean4  = CLEAN4's per-gene sign (391's pinned `psign` on the challenge axis; -1 = down, +1 = up),
  down    = every gene down.
Strata: non-essential targets, by the truth row's energy (halves), and essential for reference. Paired per-target
difference pt - clean4 with a target bootstrap. Rows are read one at a time (RSS < 1 GB).
usage: uv run python experiments/487-october-decision-review/pt487.py -> work/pt487.json
"""
import csv
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
D = ROOT / "data/donors"
T225 = ROOT / "experiments/225-source-sweep/work/tables"
ESS = set((D / "k562_essential/perts.txt").read_text().split())
TOP = 100


def bank(path):
    perts = [l.strip() for l in open(path / "perts.txt")]
    return (np.load(path / "delta.npy", mmap_mode="r"), perts, [l.strip() for l in open(path / "genes.txt")],
            {p: i for i, p in enumerate(perts)})


def row(b, t, cols):
    delta, _, _, pi = b
    r = np.asarray(delta[pi[t]], dtype=np.float64)
    v = np.where(cols >= 0, r[np.maximum(cols, 0)], 0.0)
    v[~np.isfinite(v)] = 0.0
    return v


def main():
    official = [r["gene_name"] for r in csv.DictReader(open(ROOT / "data/challenge/gene_names.csv"))]
    with np.load(ROOT / "experiments/391-consensus-block-signs/work/clean4_prior_challenge.npz", allow_pickle=False) as z:
        keys = list(z.files)
        psign = np.asarray(z["psign"], dtype=np.float64)
    assert psign.shape == (len(official),), (psign.shape, keys)
    print(f"CLEAN4 psign: {len(official)} genes, fraction down {float((psign < 0).mean()):.3f} (keys {keys})")
    rng = np.random.default_rng(4877)
    out = {"clean4_frac_down": float((psign < 0).mean())}
    for name, tb, bb, cap in (("K562: truth VIPerturb, bank K562-GWPS", D / "viperturb_k562", D / "k562_gwps", 2500),
                              ("HEK293: truth XuCao, bank Orion HEK293T", T225 / "xucao_hek293", D / "orion_hek293t", 10 ** 6)):
        T, B = bank(tb), bank(bb)
        genes = [g for g in official if g in set(T[2]) and g in set(B[2])]
        oi = {g: i for i, g in enumerate(official)}
        ps = psign[[oi[g] for g in genes]]
        tpos = {g: i for i, g in enumerate(T[2])}
        bpos = {g: i for i, g in enumerate(B[2])}
        tcols = np.array([tpos[g] for g in genes])
        bcols = np.array([bpos[g] for g in genes])
        ts = sorted(set(T[1]) & set(B[1]))
        if len(ts) > cap:
            ts = sorted(np.random.default_rng(20260925).choice(ts, size=cap, replace=False).tolist())
        gi = {g: j for j, g in enumerate(genes)}
        acc = {"pt": [], "clean4": [], "down": [], "pt_cov": []}
        energy, ess = [], []
        for t in ts:
            tr, br = row(T, t, tcols), row(B, t, bcols)
            if t in gi:
                tr[gi[t]] = 0.0
                br[gi[t]] = 0.0
            top = np.argsort(-np.abs(tr))[:TOP]
            truth = np.sign(tr[top])
            ok = truth != 0
            pt = np.where(br[top] != 0, np.sign(br[top]), ps[top])
            acc["pt"].append(float(np.mean(pt[ok] == truth[ok])))
            acc["clean4"].append(float(np.mean(ps[top][ok] == truth[ok])))
            acc["down"].append(float(np.mean(truth[ok] < 0)))
            acc["pt_cov"].append(float(np.mean(br[top] != 0)))
            energy.append(float(np.linalg.norm(tr)))
            ess.append(t in ESS)
        A = {k: np.array(v) for k, v in acc.items()}
        e, es = np.array(energy), np.array(ess)
        ne = ~es
        med = np.median(e[ne]) if ne.any() else 0.0
        strata = {"non-essential, all": ne, "non-essential, weaker half": ne & (e < med),
                  "non-essential, stronger half": ne & (e >= med), "essential": es}
        res = {"targets": len(ts), "genes": len(genes)}
        print(f"\n{name}: {len(ts)} targets, {len(genes)} genes, truth top-{TOP} genes per target")
        for sname, m in strata.items():
            if m.sum() < 10:
                continue
            idx = np.flatnonzero(m)
            d = A["pt"][idx] - A["clean4"][idx]
            bs = [d[rng.integers(0, len(d), len(d))].mean() for _ in range(4000)]
            res[sname] = {"n": int(m.sum()), "pt": round(float(A["pt"][idx].mean()), 4), "clean4": round(float(A["clean4"][idx].mean()), 4),
                          "down": round(float(A["down"][idx].mean()), 4), "pt_coverage": round(float(A["pt_cov"][idx].mean()), 3),
                          "pt_minus_clean4": round(float(d.mean()), 4),
                          "ci95": [round(float(np.quantile(bs, .025)), 4), round(float(np.quantile(bs, .975)), 4)]}
            r = res[sname]
            print(f"  {sname:30s} n {r['n']:5d}  pt {r['pt']:.4f}  CLEAN4 {r['clean4']:.4f}  all-down {r['down']:.4f}  "
                  f"pt-CLEAN4 {r['pt_minus_clean4']:+.4f} [{r['ci95'][0]:+.4f}, {r['ci95'][1]:+.4f}]  (bank covers {r['pt_coverage']:.2f})")
        out[name] = res
    (HERE / "work/pt487.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
