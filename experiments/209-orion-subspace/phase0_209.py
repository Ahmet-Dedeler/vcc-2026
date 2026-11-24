"""209 Phase 0 — the known-answer rows.  Abort on failure.

Three independent reproductions, plus rule 11's operator table:

  KA-1  204 Phase 1c/1d's one-factor test: the Orion<->Orion structural excess at
        z = +5.1 (five-line block) to +7.6 (six-line block), and the same-cell-line
        POSITIVE CONTROL (k562_gwps <-> k562_essential) at z = +11.2.  If the positive
        control does not fire the test has no power and Phase 1 is unreadable.
  KA-2  175's addendum: the two Orion donors' top PCs at cosine 0.872, HEK293T's top-1 at
        19.7% of its own between-perturbation variance, and per-donor rank-1 deflation
        lifting consensus sign accuracy on the h1-test pool 0.5697 -> 0.5770 with k >= 2
        removing signal.
  KA-3  rule 11: a path and a sha256 for every operator and every donor input.

KA-4 (PLOLO a15 0.795784 / a23 0.789121 and the a15->a23 mass-fix price +0.004014) is the
`equal` arm of `score209.py` and is adjudicated in `analyse209.py`; it cannot be run here
without the fold machinery.

Nothing outside experiments/209-orion-subspace/ is written.  204's cached
`prepare_main_source` outputs are copied into this experiment's own work/ and read from
there (`prep204.WORK` is repointed before any call), so 204's directory is read-only.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import os
import resource
import sys
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "6")
os.environ.setdefault("OMP_NUM_THREADS", "6")

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(exist_ok=True)

sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(ROOT / "experiments/204-repression-weighted-mixing"))
sys.path.insert(0, str(ROOT / "experiments/153-de-bench"))
sys.path.insert(0, str(ROOT / "experiments/172-sign-confidence-key"))

import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"
plolo.WORK.mkdir(parents=True, exist_ok=True)

import prep204  # noqa: E402

prep204.WORK = WORK                      # every cache read lands in OUR directory
from prep204 import DONORS, prepared, shared_targets  # noqa: E402
from factor204 import BLOCKS, CELL_LINE, LINES, PLATFORM, fit_rank1  # noqa: E402
from plolo import DONOR_PATHS, PANEL_SET  # noqa: E402

SEED = 20260903
NB = 2000

BLOCKS = dict(BLOCKS)
BLOCKS["h1_block"] = ("k562_gwps", "orion_hct116", "orion_hek293t", "h1_2025")


def rss_gib() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1 << 30)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


# ------------------------------------------------------------------ KA-3, rule 11
def rule11() -> dict:
    ops = {
        "production/whitener_union_pruned": ROOT / "data/donors/whitener_union_pruned.npz",
        "held_k562_gwps": ROOT / "experiments/062-cd4-post-whiten/cloud-output/whiteners/held_k562_gwps.npz",
        "held_orion_hct116": ROOT / "experiments/062-cd4-post-whiten/cloud-output/whiteners/held_orion_hct116.npz",
        "held_orion_hek293t": ROOT / "experiments/062-cd4-post-whiten/cloud-output/whiteners/held_orion_hek293t.npz",
        "NOT-LOADED/whitener_union": ROOT / "data/donors/whitener_union.npz",
        "NOT-LOADED/whitener_union_full": ROOT / "data/donors/whitener_union_full.npz",
        "NOT-LOADED/whitener_full": ROOT / "data/donors/whitener_full.npz",
    }
    out = {"operators": {}, "donor_inputs": {}}
    for tag, p in ops.items():
        s = sha256(p)
        with np.load(p, mmap_mode="r") as z:
            shapes = {k: list(np.asarray(z[k]).shape) for k in ("evals", "evecs") if k in z}
        out["operators"][tag] = dict(path=str(p.relative_to(ROOT)), sha256_8=s[:8],
                                     sha256=s, **shapes)
        print(f"  {tag:38s} {s[:8]}  {shapes}", flush=True)
    meta = ROOT / "data/donors/whitener_union_pruned.json"
    if meta.exists():
        j = json.loads(meta.read_text())
        out["production_operator_json"] = {k: j[k] for k in j if k != "columns"}
        print(f"  whitener_union_pruned.json basis = {j.get('basis')}  "
              f"prune = {j.get('basis_prune')}", flush=True)
    for d in ("k562_gwps", "orion_hct116", "orion_hek293t"):
        out["donor_inputs"][d] = {
            f: sha256(ROOT / f"data/donors/{d}/{f}")[:16]
            for f in ("delta.npy", "b_ctrl.npy")}
        print(f"  donor {d:15s} delta {out['donor_inputs'][d]['delta.npy']}  "
              f"b_ctrl {out['donor_inputs'][d]['b_ctrl.npy']}", flush=True)
    return out


# ------------------------------------------------------------------ KA-1, the factor test
def ka1_one_factor() -> dict:
    tg = shared_targets()
    P = {}
    for name, path in LINES.items():
        DONOR_PATHS[name] = path
        tag = "shared" if name in DONORS else "allshared"
        cache = WORK / f"prep-{tag}-{name}.npz"
        if not cache.exists():
            raise SystemExit(f"missing cache {cache}; refusing to write into 204's dir")
        p = prepared(name, tg, tag)
        P[name] = (np.asarray(p["effects"], dtype=np.float32),
                   np.asarray(p["target_available"], bool))
    support = {n: (np.abs(P[n][0]).sum(0) > 0) for n in LINES}

    out: dict = {}
    for bname, block in BLOCKS.items():
        rows = np.flatnonzero(np.logical_and.reduce([P[n][1] for n in block]))
        cols = np.flatnonzero(np.logical_and.reduce([support[n] for n in block]))
        if len(rows) < 30:
            continue
        k = len(block)
        M = [P[n][0][np.ix_(rows, cols)].astype(np.float64) for n in block]
        pair = np.zeros((len(rows), k, k))
        for i in range(k):
            for j in range(i, k):
                pair[:, i, j] = pair[:, j, i] = (M[i] * M[j]).sum(1)
        del M

        def stat(G):
            d = np.sqrt(np.diag(G))
            R = G / np.outer(d, d)
            c = fit_rank1(R)
            resid = R - np.outer(c, c)
            off = ~np.eye(k, dtype=bool)
            return R, c, resid, float(np.sqrt((resid[off] ** 2).mean()))

        R0, c0, res0, rms0 = stat(pair.sum(0))
        rng = np.random.default_rng(SEED)
        bs_res = []
        for _ in range(NB):
            idx = rng.integers(0, len(rows), len(rows))
            bs_res.append(stat(pair[idx].sum(0))[2][np.triu_indices(k, 1)])
        bs_res = np.asarray(bs_res)

        pairs = []
        for n, (i, j) in enumerate(itertools.combinations(range(k), 2)):
            v = bs_res[:, n]
            pairs.append(dict(pair=[block[i], block[j]], rho=float(R0[i, j]),
                              resid=float(res0[i, j]),
                              z=float(res0[i, j] / max(v.std(), 1e-12)),
                              same_cell_line=CELL_LINE[block[i]] == CELL_LINE[block[j]],
                              same_platform=PLATFORM[block[i]] == PLATFORM[block[j]]))
        off = ~np.eye(k, dtype=bool)
        out[bname] = dict(lines=list(block), n_targets=int(len(rows)),
                          n_genes=int(len(cols)),
                          n_panel=int(sum(1 for r in rows if tg[r] in PANEL_SET)),
                          rank1_rms_resid=rms0,
                          mean_abs_rho=float(np.abs(R0[off]).mean()),
                          rel_rms=float(rms0 / np.abs(R0[off]).mean()), pairs=pairs)
        print(f"\n[{bname}] n_t={len(rows)} n_panel={out[bname]['n_panel']} "
              f"n_g={len(cols)}  rank-1 RMS resid {rms0:.5f} "
              f"({100 * out[bname]['rel_rms']:.1f}% of mean|rho|)", flush=True)
        for p in sorted(pairs, key=lambda x: -abs(x["z"]))[:4]:
            flag = ("SAME LINE" if p["same_cell_line"]
                    else ("same platform" if p["same_platform"] else ""))
            print(f"    {p['pair'][0]:15s} {p['pair'][1]:15s} resid {p['resid']:+.5f} "
                  f"z {p['z']:+6.2f}   {flag}")
    return out


# ------------------------------------------------------------------ KA-2, the rank-1 rows
def ka2_rank1() -> dict:
    from debench import Surface
    from donors172 import DonorHandle

    S = Surface(str(ROOT / "experiments/164-callset-artefact/de164.a15.npz"))
    pool = S.real_sig & S.adj
    sT = np.sign(S.real_lfc)
    genes = np.asarray([str(g) for g in S.genes])
    targets = [str(p) for p in S.perts]

    D, PCS = {}, {}
    rng = np.random.default_rng(0)                      # 175's seed, unchanged
    for n in ("k562_gwps", "orion_hct116", "orion_hek293t"):
        h = DonorHandle(n)
        Sd, carried = h.rows(targets, genes)
        D[n] = (Sd.astype(np.float64), carried)
        tcols, dcols = h.project(genes)
        sub = np.sort(rng.choice(h.n_perts, min(4000, h.n_perts), replace=False))
        M = np.zeros((sub.size, len(genes)))
        for r, i in enumerate(sub):
            M[r, tcols] = (np.asarray(h.delta[i], dtype=np.float64) - h.m)[dcols]
        _, s_, Vt = np.linalg.svd(M, full_matrices=False)
        PCS[n] = (Vt, s_ ** 2 / np.sum(s_ ** 2))
        del M
    allc = np.all([c for _, c in D.values()], axis=0)
    m3 = pool & allc[:, None] & np.all([np.sign(Sd) != 0 for Sd, _ in D.values()], axis=0)

    def acc(key, m):
        return float(np.mean(key[m] == sT[m]))

    sweep = {}
    for k in (0, 1, 2, 3, 5, 10, 20):
        defl = {}
        for n, (Sd, _) in D.items():
            Vt = PCS[n][0][:k]
            defl[n] = Sd - (Sd @ Vt.T) @ Vt if k > 0 else Sd
        cons = sum(defl.values()) / 3
        w211 = (2 * defl["k562_gwps"] + defl["orion_hct116"]
                + defl["orion_hek293t"]) / 4
        sweep[k] = dict(
            k562=acc(np.sign(defl["k562_gwps"]), m3),
            hct116=acc(np.sign(defl["orion_hct116"]), m3),
            hek293t=acc(np.sign(defl["orion_hek293t"]), m3),
            equal=acc(np.sign(cons), m3), w211=acc(np.sign(w211), m3))
        print(f"  k={k:<3d} K562 {sweep[k]['k562']:.4f}  HCT {sweep[k]['hct116']:.4f}  "
              f"HEK {sweep[k]['hek293t']:.4f}  equal {sweep[k]['equal']:.4f}  "
              f"(2,1,1) {sweep[k]['w211']:.4f}", flush=True)

    cos_top = {f"{a}x{b}": abs(float(PCS[a][0][0] @ PCS[b][0][0]))
               for a, b in itertools.combinations(PCS, 2)}
    var_top = {n: float(PCS[n][1][0]) for n in PCS}
    print(f"  top-PC cosines {cos_top}")
    print(f"  top-1 variance share {var_top}")
    return dict(sweep={str(k): v for k, v in sweep.items()},
                top_pc_cosine=cos_top, top1_variance_share=var_top,
                n_pool=int(m3.sum()))


def main() -> None:
    out: dict = {}
    print("=== KA-3 / rule 11: operators and donor inputs, by path and sha ===", flush=True)
    out["rule11"] = rule11()

    print("\n=== KA-2: 175's rank-1 deflation on the h1-test sign pool ===", flush=True)
    out["ka2_rank1"] = ka2_rank1()

    print("\n=== KA-1: 204's one-factor test and its positive control ===", flush=True)
    out["ka1_one_factor"] = ka1_one_factor()

    # -------- verdicts, hard-coded expectations from the mandate's sources ------------
    k2 = out["ka2_rank1"]
    checks = [
        ("KA-2 top-PC cos(HCT116, HEK293T) = 0.872",
         k2["top_pc_cosine"]["orion_hct116xorion_hek293t"], 0.872, 0.003),
        ("KA-2 HEK293T top-1 variance share = 0.197",
         k2["top1_variance_share"]["orion_hek293t"], 0.197, 0.002),
        ("KA-2 equal-mean sign accuracy, k=0 = 0.5697",
         k2["sweep"]["0"]["equal"], 0.5697, 0.0002),
        ("KA-2 equal-mean sign accuracy, k=1 = 0.5770",
         k2["sweep"]["1"]["equal"], 0.5770, 0.0002),
    ]
    f = out["ka1_one_factor"]
    for bname, pair, expect, tol in (
            ("five_no_sameline", ("orion_hct116", "orion_hek293t"), 5.14, 0.6),
            ("replogle_block", ("orion_hct116", "orion_hek293t"), 7.56, 0.8),
            ("replogle_block", ("k562_gwps", "k562_essential"), 11.2, 1.2)):
        z = [p["z"] for p in f[bname]["pairs"] if tuple(p["pair"]) == pair][0]
        checks.append((f"KA-1 {bname} z({pair[0]}, {pair[1]}) = {expect}", z, expect, tol))

    print("\n=== VERDICTS ===")
    ok = True
    for label, got, expect, tol in checks:
        good = abs(got - expect) <= tol
        ok &= good
        print(f"  [{'PASS' if good else 'FAIL'}] {label:58s} measured {got:+.4f} "
              f"(tol {tol})")
    out["checks"] = [dict(label=a, measured=b, expected=c, tol=d,
                          passed=bool(abs(b - c) <= d)) for a, b, c, d in checks]
    out["all_passed"] = bool(ok)
    out["peak_rss_gib"] = rss_gib()
    (WORK / "phase0_209.json").write_text(json.dumps(out, indent=2))
    print(f"\npeak RSS {rss_gib():.2f} GiB -> work/phase0_209.json")
    if not ok:
        raise SystemExit("Phase 0 known-answer row FAILED — aborting.")


if __name__ == "__main__":
    main()
