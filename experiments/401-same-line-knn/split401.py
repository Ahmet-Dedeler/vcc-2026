#!/usr/bin/env python3
"""401 -- does a SAME-LINE vote bank beat the six-bank family for the neighbourhood tilt? Nadig Jurkat split-half: targets split
into two halves (5 seeded splits); for each eval-half target, neighbours = its 25 nearest GWPS perturbations (self excluded);
votes from (a) the six-bank family (399), (b) the fit-half Jurkat rows only (same line, other targets), (c) six banks + fit-half
Jurkat rows; shrinkage toward the corresponding marginal (m = 400 for (a)/(c); m = 50 and 400 for (b)); truth = the eval-half
target's Jurkat delta sign on Jurkat's top-2,200 block (own targets excluded). Also the marginals alone. Units: block accuracy."""
import numpy as np, json, time
from pathlib import Path
REPO = Path("~/Code/virtual-cell-challenge-2026"); OUT = REPO / "experiments/401-same-line-knn/work"; OUT.mkdir(parents=True, exist_ok=True)
def load_bank(name):
    d = (REPO / "data/donors" / name) if name != "nadig_hepg2" else (REPO / "experiments/166-donor-acquisition/donors/nadig_hepg2")
    D = np.asarray(np.load(d / "delta.npy", mmap_mode="r"), dtype=np.float32)
    P = np.array([l.strip() for l in open(d / "perts.txt")]); G = np.array([l.strip() for l in open(d / "genes.txt")])
    B = np.load(d / "b_ctrl.npy") if (d / "b_ctrl.npy").exists() else None
    return D, P, G, B
t0 = time.time()
GW, GWP, GWG, _ = load_bank("k562_gwps"); gwp = {p: i for i, p in enumerate(GWP)}; gwg = {g: i for i, g in enumerate(GWG)}
GWc = GW - GW.mean(1, keepdims=True); GWn = GWc / np.maximum(np.linalg.norm(GWc, axis=1, keepdims=True), 1e-12)
JD, JP, JG, JB = load_bank("nadig_jurkat"); jp = {p: i for i, p in enumerate(JP)}; jg = {g: i for i, g in enumerate(JG)}
others = {b: load_bank(b) for b in ("k562_essential", "nadig_hepg2", "orion_hct116", "orion_hek293t")}
tset = set(JP); order = np.argsort(-JB, kind="stable")
block = [JG[j] for j in order if JG[j] in gwg and JG[j] not in tset][:2200]
bl = np.array([jg[g] for g in block]); gwi = np.array([gwg[g] for g in block])
vb = []
for b, (D, P, G, _) in others.items():
    pi = {p: i for i, p in enumerate(P)}; gi = {g: i for i, g in enumerate(G)}; vb.append((D, pi, np.array([gi.get(g, -1) for g in block])))
def votes(rows, use_six, jur_rows_allowed):
    """down/total over GWPS rows (+ other banks if use_six) (+ Jurkat rows of the same perturbations restricted to jur_rows_allowed)."""
    dn = np.zeros(len(block)); nn = np.zeros(len(block)); names = [GWP[j] for j in rows]
    if use_six:
        sub = GW[rows][:, gwi]; nz = sub != 0; dn += ((sub < 0) & nz).sum(0); nn += nz.sum(0)
        for D, pi, gm in vb:
            ridx = [pi[p] for p in names if p in pi]
            if not ridx: continue
            cv = gm >= 0; s2 = D[ridx][:, np.where(cv, gm, 0)]; nz2 = (s2 != 0) & cv[None, :]; dn += ((s2 < 0) & nz2).sum(0); nn += nz2.sum(0)
    if jur_rows_allowed is not None:
        ridx = [jp[p] for p in names if p in jp and jur_rows_allowed[jp[p]]]
        if ridx:
            s3 = JD[ridx][:, bl]; nz3 = s3 != 0; dn += ((s3 < 0) & nz3).sum(0); nn += nz3.sum(0)
    return np.where(nn > 0, dn / np.maximum(nn, 1), np.nan), nn
def sg(P, fb): return np.where(np.isfinite(P), np.where(P > 0.5, -1.0, 1.0), fb)
targets = [i for i, t in enumerate(JP) if t in gwp]
res = {}; rng = np.random.default_rng(20260916)
for rep in range(3):
    perm = rng.permutation(len(JP)); fit = np.zeros(len(JP), bool); fit[perm[: len(JP) // 2]] = True
    ev = [i for i in targets if not fit[i]]
    P6m, _ = votes(np.arange(len(GWP)), True, None)                      # six-bank marginal
    PJm, _ = votes(np.arange(len(GWP)), False, fit)                      # same-line (fit half) marginal
    P7m, _ = votes(np.arange(len(GWP)), True, fit)                       # six banks + fit-half Jurkat marginal
    truth = JD[ev][:, bl]; T_nz = truth != 0; T_s = np.sign(truth)
    arms = {"six_marginal": np.tile(sg(P6m, -1.0), (len(ev), 1)), "jurkat_marginal": np.tile(sg(PJm, -1.0), (len(ev), 1)),
            "seven_marginal": np.tile(sg(P7m, -1.0), (len(ev), 1))}
    K = 25
    for name in ("six_knn_m400", "jurkat_knn_m50", "jurkat_knn_m400", "seven_knn_m400"):
        arms[name] = np.zeros((len(ev), len(block)))
    sims = GWn @ GWn[[gwp[JP[i]] for i in ev]].T
    for j, i in enumerate(ev):
        s_ = sims[:, j].copy(); s_[gwp[JP[i]]] = -np.inf; nb = np.argsort(-s_)[:K]
        P6, n6 = votes(nb, True, None); PJ, nJ = votes(nb, False, fit); P7, n7 = votes(nb, True, fit)
        def shrink(P, nn, Pm, m): return np.where(np.isfinite(P) & np.isfinite(Pm), (nn * np.nan_to_num(P) + m * np.nan_to_num(Pm)) / (nn + m), np.where(np.isfinite(P), P, Pm))
        arms["six_knn_m400"][j] = sg(shrink(P6, n6, P6m, 400.0), -1.0)
        arms["jurkat_knn_m50"][j] = sg(shrink(PJ, nJ, PJm, 50.0), -1.0)
        arms["jurkat_knn_m400"][j] = sg(shrink(PJ, nJ, PJm, 400.0), -1.0)
        arms["seven_knn_m400"][j] = sg(shrink(P7, n7, P7m, 400.0), -1.0)
    for name, S in arms.items():
        ok = T_nz; m = (S == T_s) & ok; acc = (m.sum(1) / np.maximum(ok.sum(1), 1)).mean()
        res.setdefault(name, []).append(float(acc))
    print(f"rep {rep}: " + "  ".join(f"{k} {v[-1]:.4f}" for k, v in res.items()), flush=True)
summary = {k: {"mean": round(float(np.mean(v)), 4), "sd": round(float(np.std(v)), 4)} for k, v in res.items()}
print(json.dumps(summary, indent=1)); (OUT / "results401.json").write_text(json.dumps(summary, indent=1)); print(f"{time.time()-t0:.0f}s")
