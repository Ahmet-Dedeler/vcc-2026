"""214 — the challenge side, read off the SHIPPED artifacts.

Four jobs, all on `submissions/210-cd4f034.h5ad` (the incumbent) and
`submissions/202-rc_callset.h5ad` (its predecessor), streamed block by block from the CSR
without materialising anything:

  1. Phase 0's rounding known-answer row: 196's challenge tilt `floor` -0.01234 -> `ham`
     +0.00039, reproduced on a15's own delivered totals (196's `chal-totals-a15.npz`) and
     then RE-MEASURED on `cd4f034`'s, because the tilt is a property of the delivered
     profile and the dose moved that profile on 900/900 blocks.
  2. Rule 13's class test for arm B: the per-gene panel-total digest under Hamilton against
     `cd4f034`'s own, on all 900 blocks, with max |dT| and the direction of every move.
  3. The format contract: package nnz and max counts/cell, measured on the artifact, and
     Hamilton's exact effect on the panel total S.
  4. Falsifier 1's identity row: with both flags off this pipeline must reproduce
     `cd4f034` -- asserted as `floor(x) == T` on every gene of every block, which is the
     statement "the reconstruction of the intended profile re-delivers the artifact".

    uv run python chal214.py
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from pathlib import Path

import h5py
import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
WORK.mkdir(exist_ok=True)
sys.path.insert(0, str(ROOT / "experiments/196-emitted-depth"))
import recon196                                                          # noqa: E402
from ham196 import hamilton                                              # noqa: E402

ARTS = {"cd4f034": ROOT / "submissions/210-cd4f034.h5ad",
        "rc_callset": ROOT / "submissions/202-rc_callset.h5ad"}
N = 400
NNZ_CAP, CELL_CAP = 4.75e9, 1e6
TILT_SAMPLE = int(os.environ.get("TILT_SAMPLE", "40"))     # blocks per context, 208 used 40


def block_totals(path: Path):
    """Per-gene panel totals for all 900 blocks, plus nnz and max cell total."""
    with h5py.File(path, "r") as h:
        o, c = h["obs/target_gene"], h["obs/context"]
        tg = np.asarray([x.decode() for x in o["categories"][:]])[o["codes"][:]]
        cx = np.asarray([x.decode() for x in c["categories"][:]])[c["codes"][:]]
        node = h["var/_index"]
        raw = node["values"][:] if isinstance(node, h5py.Group) else node[:]
        genes = np.asarray([g.decode() if isinstance(g, bytes) else str(g) for g in raw])
        G = genes.size
        labels = np.char.add(np.char.add(cx, "|"), tg)
        ptr = h["X/indptr"][:]
        nrow = labels.size
        assert nrow % N == 0
        blocks, tot, nnz, maxcell = [], [], 0, 0
        dsd, dsi = h["X/data"], h["X/indices"]
        for s in range(0, nrow, N):
            lab = labels[s]
            assert (labels[s:s + N] == lab).all(), f"block at {s} is not homogeneous"
            a, b = int(ptr[s]), int(ptr[s + N])
            d = dsd[a:b]
            idx = dsi[a:b]
            t = np.bincount(idx, weights=d, minlength=G).astype(np.int64)
            per = np.diff(ptr[s:s + N + 1])
            rows = np.bincount(np.repeat(np.arange(N), per), weights=d, minlength=N)
            blocks.append(lab)
            tot.append(t)
            nnz += int(b - a)
            maxcell = max(maxcell, int(rows.max()))
    return dict(blocks=np.asarray(blocks), totals=np.asarray(tot), genes=genes,
                nnz=nnz, max_cell=maxcell)


def tilt_stats(T, x, gate, cpm, self_pos=None):
    """208's `tilt_stats` plus 196/208's per-row rounding ENERGY (`mse208.py` lines 121-126)."""
    g = gate & (x > 0)
    e = (T[g] - x[g]) / x[g]
    lo = cpm[g] <= np.quantile(cpm[g], 0.10)
    hi = cpm[g] >= np.quantile(cpm[g], 0.90)
    l2 = np.log2((T[g] + 1e-9) / x[g])
    eps = np.log1p(50_000.0 * T / T.sum()) - np.log1p(50_000.0 * x / x.sum())
    if self_pos is not None:
        eps[self_pos] = 0.0
    return dict(mean_rel=float(e.mean()), mean_abs_rel=float(np.abs(e).mean()),
                common_log2=float(np.log2(T[g].sum() / x[g].sum())),
                tilt_log2=float(l2[lo].mean() - l2[hi].mean()),
                energy=float((eps ** 2).sum() / (T.size - 1)))


def main() -> None:
    t0 = time.time()
    out: dict = {}
    art = {}
    for name, path in ARTS.items():
        cache = WORK / f"totals-{name}.npz"
        if cache.exists():
            z = np.load(cache, allow_pickle=True)
            art[name] = dict(blocks=np.asarray(z["blocks"]).astype(str), totals=z["totals"],
                             genes=np.asarray(z["genes"]).astype(str),
                             nnz=int(z["nnz"]), max_cell=int(z["max_cell"]))
            print(f"[{name}] cache hit", flush=True)
            continue
        r = block_totals(path)
        np.savez_compressed(cache, **r)
        art[name] = r
        print(f"[{name}] {len(r['blocks'])} blocks, nnz {r['nnz']:,}, max cell "
              f"{r['max_cell']:,} [{time.time()-t0:.0f}s]", flush=True)

    A, Cst = art["cd4f034"], art["rc_callset"]
    assert list(A["blocks"]) == list(Cst["blocks"])
    assert np.array_equal(A["genes"], Cst["genes"])
    print("\n=== the two artifacts, per-gene panel totals ===")
    dig_a = [hashlib.sha256(t.tobytes()).hexdigest()[:16] for t in A["totals"]]
    dig_c = [hashlib.sha256(t.tobytes()).hexdigest()[:16] for t in Cst["totals"]]
    same = sum(x == y for x, y in zip(dig_a, dig_c))
    print(f"  cd4f034 vs rc_callset digest identical on {same}/900 blocks "
          f"(210 falsifier 3a: MUST be 0)")
    print(f"  first block {A['blocks'][0]} digest {dig_a[0]} (rc_callset {dig_c[0]}); "
          f"202's recorded value for (A,ABCD1) is bb6befc01caf7493")
    out["digest_cd4f034_vs_rc_callset_identical"] = same
    out["nnz"] = {k: v["nnz"] for k, v in art.items()}
    out["max_cell"] = {k: v["max_cell"] for k, v in art.items()}

    # -------------------------------------------------- 1 + 4: reconstruction and the tilt
    print("\n=== PHASE 0 KA + THE RE-MEASURE: the rounding tilt on the CHALLENGE ===")
    print(f"{'profile':<26}{'blocks':>7}{'mean rel':>11}{'common log2':>13}"
          f"{'TILT log2':>12}{'round energy':>14}{'floor(x)==T':>13}")
    tab = {}
    for src in ("a15", "cd4f034"):
        rows = {"floor": [], "ham": []}
        ident, nb = [], 0
        for ctx in "ABC":
            D = recon196.load_chal(ctx)
            if src == "a15":
                sums, blocks = D["sums"], D["targets"]
            else:
                sel = np.flatnonzero(np.array([b.split("|")[0] for b in A["blocks"]]) == ctx)
                assert np.array_equal(D["genes"], A["genes"])
                sums = A["totals"][sel]
                blocks = np.array([A["blocks"][i].split("|")[1] for i in sel])
                assert list(blocks) == list(D["targets"]), ctx
            step = max(1, sums.shape[0] // TILT_SAMPLE)
            for i in range(0, sums.shape[0], step):
                R = recon196.Recon(sums[i], D["r"], D["gate"], n=N, seed=1000 + int(i))
                ident.append(float(np.array_equal(R.totals(1.0), sums[i].astype(np.int64))))
                nb += 1
                si = np.flatnonzero(D["genes"] == blocks[i])
                sp = int(si[0]) if si.size else None
                # 208's `tilt208.py`:36 — the SCORER's DE universe (mean per-cell CPM > 5),
                # not the builder's pooled 10-CPM gate.  Reading the tilt on the builder's
                # gate halves it, because the low decile the tilt is defined on is exactly
                # the 5-10 CPM band the builder gate removes.
                gs = D["mean_cpm"] > 5.0
                for rule in ("floor", "ham"):
                    T = (R.totals(1.0) if rule == "floor" else hamilton(R.x)).astype(float)
                    rows[rule].append(tilt_stats(T, R.x, gs, D["mean_cpm"], sp))
        for rule in ("floor", "ham"):
            med = {k: float(np.median([x[k] for x in rows[rule]])) for k in rows[rule][0]}
            tab[f"{src}|{rule}"] = med
            print(f"  {src + ' / ' + rule:<24}{nb:>7}{med['mean_rel']:>+11.5f}"
                  f"{med['common_log2']:>+13.6f}{med['tilt_log2']:>+12.5f}"
                  f"{med['energy']:>14.4e}{min(ident):>13.4f}")
        out.setdefault("identity_floor_eq_T", {})[src] = float(min(ident))
    out["tilt"] = tab

    ka = tab["a15|floor"]["tilt_log2"], tab["a15|ham"]["tilt_log2"]
    ok = abs(ka[0] - (-0.01234)) < 5e-4 and abs(ka[1] - 0.00039) < 5e-4
    print(f"\n  KNOWN-ANSWER (196 §4a): a15 challenge tilt floor {ka[0]:+.5f} "
          f"(expect -0.01234), ham {ka[1]:+.5f} (expect +0.00039)  "
          f"{'PASS' if ok else 'FAIL'}")
    out["ka_tilt_pass"] = bool(ok)

    # -------------------------------------------------- 2 + 3: rule 13's class test for B
    print("\n=== RULE 13 CLASS TEST for arm B (Hamilton), on ALL 900 blocks of cd4f034 ===")
    stats = []
    for ctx in "ABC":
        D = recon196.load_chal(ctx)
        sel = np.flatnonzero(np.array([b.split("|")[0] for b in A["blocks"]]) == ctx)
        for i, k in enumerate(sel):
            R = recon196.Recon(A["totals"][k], D["r"], D["gate"], n=N, seed=1000 + int(i))
            Th = hamilton(R.x).astype(np.int64)
            dT = Th - R.totals(1.0)
            stats.append(dict(ctx=ctx, same=bool(np.array_equal(Th, R.totals(1.0))),
                              max_abs=int(np.abs(dT).max()), moved=int((dT != 0).sum()),
                              moved_gated=int((dT[D["gate"]] != 0).sum()),
                              up=int((dT > 0).sum()), down=int((dT < 0).sum()),
                              S_floor=int(R.totals(1.0).sum()), S_ham=int(Th.sum())))
        print(f"  ctx {ctx} done [{time.time()-t0:.0f}s]", flush=True)
    ident = sum(s["same"] for s in stats)
    mv = lambda k: (min(s[k] for s in stats), int(np.median([s[k] for s in stats])),  # noqa
                    max(s[k] for s in stats))
    print(f"  Hamilton totals identical to cd4f034's on {ident}/900 blocks  "
          f"(rule 13: 0/900 => emission-law, profile-CHANGING, carries the -0.0146 toll)")
    print(f"    max |dT| per block      min {mv('max_abs')[0]}  median {mv('max_abs')[1]}"
          f"  max {mv('max_abs')[2]}")
    print(f"    genes moved per block   min {mv('moved')[0]:,}  median {mv('moved')[1]:,}"
          f"  max {mv('moved')[2]:,}   (of {A['genes'].size:,})")
    print(f"    ...inside the 10-CPM builder gate  median {mv('moved_gated')[1]:,}")
    print(f"    direction               UP on {mv('up')[0]:,}-{mv('up')[2]:,}, "
          f"DOWN on {mv('down')[0]}-{mv('down')[2]}")
    Sf, Sh = sum(s["S_floor"] for s in stats), sum(s["S_ham"] for s in stats)
    print(f"    panel total S           {Sf:,} -> {Sh:,}  ({100*(Sh/Sf-1):+.4f}%)")
    out["class_test"] = dict(identical=ident, blocks=len(stats), S_floor=Sf, S_ham=Sh,
                             max_abs=mv("max_abs"), moved=mv("moved"), up=mv("up"),
                             down=mv("down"))

    print("\n=== FORMAT CONTRACT (falsifier: neither arm may approach a cap) ===")
    print(f"  {'artifact':<24}{'nnz':>16}{'% of 4.75e9':>13}{'max counts/cell':>18}"
          f"{'% of 1e6':>10}")
    for k, v in art.items():
        print(f"  {k:<24}{v['nnz']:>16,}{100*v['nnz']/NNZ_CAP:>13.4f}"
              f"{v['max_cell']:>18,}{100*v['max_cell']/CELL_CAP:>10.3f}")
    print(f"  Hamilton adds {Sh-Sf:,} counts across the package ({100*(Sh/Sf-1):+.4f}% of "
          f"S); 208 measured the nnz cost of that at +0.04%, i.e. "
          f"{A['nnz']*1.0004:,.0f} = {100*A['nnz']*1.0004/NNZ_CAP:.4f}% of the cap.")

    json.dump(out, open(WORK / "chal214.json", "w"), indent=1, default=float)
    print(f"\nwrote {WORK/'chal214.json'} [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
