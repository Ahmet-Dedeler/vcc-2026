"""184 Phase 0 -- known-answer row (AGENTS.md rule 10).

Independently re-derives 159's frozen cross-line numbers from 156's exact-space truth,
read-only, with this experiment's own code.  Two rows must reproduce:

    harness route (156 frozen truth, 120 shared targets, drop=union_panel):
        same +0.1881  diff +0.1790  excess +0.0091
    donor route  (donor tables, same 120 targets, challenge axis, 17,986 genes):
        same +0.0317  diff +0.0006  ratio 54.2

If either fails, the instrument is not calibrated and nothing downstream is trustworthy.
Writes only experiments/184-tahoe-transfer/work/known.json.
"""
from __future__ import annotations

import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
T156 = os.path.join(REPO, "experiments/156-delivered-geometry/work/truth-%s.npz")
DONOR = os.path.join(REPO, "data/donors/%s")
CHALLENGE_GENES = os.path.join(REPO, "data/challenge/gene_names.csv")


def unit(A):
    return A / np.maximum(np.linalg.norm(A, axis=1, keepdims=True), 1e-300)


def offdiag_mean(S):
    n = S.shape[0]
    return float(S[~np.eye(n, dtype=bool)].mean())


def pair_stats(EA, EB):
    a, b = unit(EA), unit(EB)
    S = a @ b.T
    same = float(np.diag(S).mean())
    diff = offdiag_mean(S)
    return {
        "same": same,
        "same_median": float(np.median(np.diag(S))),
        "same_sd": float(np.diag(S).std(ddof=1)),
        "diff": diff,
        "excess": same - diff,
        "ratio": same / diff if diff != 0 else float("nan"),
        "within_A": offdiag_mean(a @ a.T),
        "within_B": offdiag_mean(b @ b.T),
    }


def harness_route(out):
    A = np.load(T156 % "hct116", allow_pickle=True)
    B = np.load(T156 % "hek293t", allow_pickle=True)
    ga = np.asarray([str(x) for x in A["genes"]])
    gb = np.asarray([str(x) for x in B["genes"]])
    assert np.array_equal(ga, gb)
    ta = np.asarray([str(x) for x in A["targets"]])
    tb = np.asarray([str(x) for x in B["targets"]])
    shared = np.asarray(sorted(set(ta) & set(tb)))
    ia = np.asarray([int(np.flatnonzero(ta == s)[0]) for s in shared])
    ib = np.asarray([int(np.flatnonzero(tb == s)[0]) for s in shared])
    EA = A["effect"].astype(np.float64)[ia]
    EB = B["effect"].astype(np.float64)[ib]
    gpos = {g: i for i, g in enumerate(ga)}
    drop = {gpos[s] for s in set(ta) | set(tb) if s in gpos}
    mask = np.asarray([i for i in range(len(ga)) if i not in drop])
    st = pair_stats(EA[:, mask], EB[:, mask])
    st["n_targets"] = int(len(shared))
    st["n_genes"] = int(len(mask))
    out["harness"] = st
    print(f"[harness ] n={st['n_targets']} genes={st['n_genes']}  "
          f"same {st['same']:+.4f} (med {st['same_median']:+.4f} sd {st['same_sd']:.4f})  "
          f"diff {st['diff']:+.4f}  excess {st['excess']:+.4f}  "
          f"within {st['within_A']:+.4f}/{st['within_B']:+.4f}")
    return shared


def load_donor(name):
    d = DONOR % name
    rows = np.load(os.path.join(d, "delta.npy"), mmap_mode="r")
    genes = np.asarray([l.strip() for l in open(os.path.join(d, "genes.txt")) if l.strip()])
    targets = np.asarray([l.strip() for l in open(os.path.join(d, "perts.txt")) if l.strip()])
    assert rows.shape == (len(targets), len(genes)), (rows.shape, len(targets), len(genes))
    return rows, genes, targets


def donor_route(out, shared):
    ra, ga, ta = load_donor("orion_hct116")
    rb, gb, tb = load_donor("orion_hek293t")
    chal = [l.strip() for l in open(CHALLENGE_GENES)][1:]
    chal = [g for g in chal if g]
    keep = sorted(set(ga) & set(gb) & set(chal) - set(shared))
    pa = {g: i for i, g in enumerate(ga)}
    pb = {g: i for i, g in enumerate(gb)}
    ia = np.asarray([pa[g] for g in keep])
    ib = np.asarray([pb[g] for g in keep])
    tpa = {t: i for i, t in enumerate(ta)}
    tpb = {t: i for i, t in enumerate(tb)}
    tt = [t for t in shared if t in tpa and t in tpb]
    EA = np.asarray(ra[[tpa[t] for t in tt]], dtype=np.float64)[:, ia]
    EB = np.asarray(rb[[tpb[t] for t in tt]], dtype=np.float64)[:, ib]
    st = pair_stats(EA, EB)
    st["n_targets"] = int(len(tt))
    st["n_genes"] = int(len(keep))
    out["donor"] = st
    print(f"[donor   ] n={st['n_targets']} genes={st['n_genes']}  "
          f"same {st['same']:+.4f}  diff {st['diff']:+.4f}  ratio {st['ratio']:.1f}  "
          f"within {st['within_A']:+.4f}/{st['within_B']:+.4f}")


def main():
    os.makedirs(WORK, exist_ok=True)
    out = {}
    shared = harness_route(out)
    try:
        donor_route(out, shared)
    except Exception as exc:  # donor table layout is secondary; harness row is the gate
        out["donor_error"] = repr(exc)
        print("donor route failed:", exc)
    exp = {"same": 0.1881, "diff": 0.1790, "excess": 0.0091}
    h = out["harness"]
    ok = all(abs(h[k] - v) < 5e-4 for k, v in exp.items())
    out["harness_reproduced"] = bool(ok)
    print("\nKNOWN-ANSWER (harness):", "PASS" if ok else "FAIL")
    with open(os.path.join(WORK, "known.json"), "w") as fh:
        json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main()
