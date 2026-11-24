"""518 — a multi-screen per-target response-STRENGTH rank as magnitude information for the bulk key, scored out of fold.

Robust strength of target p in one screen (only its panel rows are read; delta = log-pseudobulk vs that screen's NTC):
    genes = the screen's genes > 10 CPM in its control, minus the 300 panel genes;
    scale_g = 1.4826 MAD_p(delta_pg) over the screen's covered panel targets (most are weak: a noise scale);
    z = delta / scale;  noise_p = median_g |z_pg| / 0.6745 (per-target noise: cells per perturbation vary 20-1963);
    strength_p = log1p #{g : |z_pg| > 4 noise_p}.
Checked against CD4's DESeq2 nDE (a well-powered count) before use. Screens: the production sources (K562-GWPS, Orion
HCT116, Orion HEK293T; on a fold only the fold's own SOURCES), VIPerturb-seq K562 (284/300), the KOLF2.1J atlas panel rows
(273/300; dropped on the kolf2 fold, where it is the truth), CD4 T cells (DESeq2 nDE, 3 conditions, 283/300).
Features = centred percentile ranks; a set's feature = the mean of its screens' ranks over the screens that carry a target.
One multiplier exp(theta z) on the production key; theta fitted on the training folds, scored held out (lofo518's route).

    uv run python experiments/518-bulk-mse-gap/strength518.py
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import json  # noqa: E402
import warnings  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
from cd4feat518 import cd4_features  # noqa: E402
from lofo518 import crank, fit, gain  # noqa: E402

warnings.filterwarnings("ignore")
WORK = HERE / "work"
W540 = ROOT / "experiments/540-lodo-clip-scan/work"
BASE = "V1_0.5_3@R16"
FOLDS = ("k562_gwps", "orion_hct116", "orion_hek293t", "kolf2")
CANCER = FOLDS[:3]
PANEL = [l.split(",")[0] for l in (ROOT / "data/challenge/pert_counts.csv").read_text().splitlines()[1:] if l.strip()]
PSET = set(PANEL)
BANKS = {"k562_gwps": ROOT / "data/donors/k562_gwps", "orion_hct116": ROOT / "data/donors/orion_hct116",
         "orion_hek293t": ROOT / "data/donors/orion_hek293t", "viperturb_k562": ROOT / "data/donors/viperturb_k562",
         "kolf2": ROOT / "data/external/kolf2"}
CD4 = ROOT / "data/external/cd4t/cloud-compact/cd4_compact"


def robust_strength(path: Path, targets=PANEL):
    perts = [l.strip() for l in (path / "perts.txt").read_text().splitlines() if l.strip()]
    genes = [l.strip() for l in (path / "genes.txt").read_text().splitlines() if l.strip()]
    b = np.load(path / "b_ctrl.npy").astype(np.float64)
    ex = np.expm1(b)
    cpm = 1e6 * ex / ex.sum()
    keep = (cpm > 10.0) & np.asarray([g not in PSET for g in genes])
    pi = {p: i for i, p in enumerate(perts)}
    rows = [pi[t] for t in targets if t in pi]
    cov = [t for t in targets if t in pi]
    delta = np.load(path / "delta.npy", mmap_mode="r")
    R = np.asarray(delta[np.asarray(sorted(rows))], np.float64)
    order = np.argsort(np.argsort(rows))
    R = R[order][:, keep]
    med = np.median(R, 0)
    scale = 1.4826 * np.median(np.abs(R - med), 0)
    ok = scale > 0
    Z = R[:, ok] / scale[ok]
    noise = np.median(np.abs(Z), 1) / 0.6745
    cnt = (np.abs(Z) > 4.0 * noise[:, None]).sum(1)
    out = dict(zip(cov, np.log1p(cnt).tolist()))
    return np.asarray([out.get(t, np.nan) for t in targets])


def cd4_robust(targets=PANEL):
    vals = []
    for c in ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr"):
        vals.append(robust_strength(CD4 / c, targets))
    return np.nanmean(np.vstack(vals), 0)


def main():
    S = {k: robust_strength(p) for k, p in BANKS.items()}
    cd4 = cd4_features(PANEL)
    S["cd4_deseq"] = cd4["lognde_mean"]
    S["cd4_robust"] = cd4_robust()
    ok = np.isfinite(S["cd4_deseq"]) & np.isfinite(S["cd4_robust"])
    sp_cd4 = float(np.corrcoef(np.argsort(np.argsort(S["cd4_deseq"][ok])), np.argsort(np.argsort(S["cd4_robust"][ok])))[0, 1])
    print(f"validation: Spearman(robust strength on CD4 delta, CD4 DESeq2 nDE) = {sp_cd4:+.3f} (n {int(ok.sum())})")
    names = list(S)
    M = np.vstack([S[k] for k in names])
    C = np.full((len(names), len(names)), np.nan)
    for i in range(len(names)):
        for j in range(len(names)):
            m = np.isfinite(M[i]) & np.isfinite(M[j])
            C[i, j] = np.corrcoef(np.argsort(np.argsort(M[i][m])), np.argsort(np.argsort(M[j][m])))[0, 1]
    print("cross-screen Spearman of per-target strength:")
    print(" " * 16 + " ".join(f"{n[:9]:>9s}" for n in names))
    for i, n in enumerate(names):
        print(f"{n[:15]:16s}" + " ".join(f"{C[i, j]:+9.2f}" for j in range(len(names))))
    np.savez_compressed(WORK / "strength518.npz", panel=np.asarray(PANEL), **S)
    # ---------------- LOFO on the folds
    data = {}
    for f in FOLDS:
        z = np.load(W540 / f"arrays540-{f}.npz")
        x, a, d = z[f"{BASE}__x"], z[f"{BASE}__a"], z["d"]
        targets = [str(t) for t in np.load(WORK / f"dump518-{f}.npz")["targets"]]
        idx = np.asarray([PANEL.index(t) for t in targets])
        srcs = [s for s in ("k562_gwps", "orion_hct116", "orion_hek293t") if s != f]
        R = lambda k: crank(S[k][idx])

        def meanrank(keys):
            cols = []
            for k in keys:
                v = S[k][idx]
                r = crank(v)
                r[~np.isfinite(v)] = np.nan
                cols.append(r)
            m = np.nanmean(np.vstack(cols), 0)
            return np.where(np.isfinite(m), m, 0.0)
        ext = ["cd4_deseq", "viperturb_k562"] + (["kolf2"] if f != "kolf2" else [])
        feats = {"cd4_deseq": R("cd4_deseq"), "cd4_robust": R("cd4_robust"), "viperturb_k562": R("viperturb_k562"),
                 "kolf2": R("kolf2") if f != "kolf2" else np.zeros(len(targets)),
                 "sources": meanrank(srcs), "meta_ext": meanrank(ext), "meta_all": meanrank(ext + srcs),
                 "meta_ext_no_viper": meanrank([e for e in ext if e != "viperturb_k562"])}
        base = float(np.mean([x[j].sum() / np.sqrt(a[j].sum() * d.sum()) for j in range(3)]))
        data[f] = dict(x=x, a=a, d=d, base=base, feats=feats)
    out = dict(validation_cd4=sp_cd4, cross=dict(names=names, spearman=C.tolist()), lofo={})
    print(f"\n{'feature':20s} {'LOFO cancer':>12s}  per fold (k562 / hct116 / hek293t)   kolf2 (fit on cancer)   theta(cancer3)")
    for nm in data["k562_gwps"]["feats"]:
        dd = {f: (data[f]["x"], data[f]["a"], data[f]["d"], data[f]["feats"][nm][:, None], data[f]["base"]) for f in FOLDS}
        held = {}
        for h in CANCER:
            th = fit(dd, [f for f in CANCER if f != h], [0])
            held[h] = gain(th, dd[h], [0])
        th = fit(dd, list(CANCER), [0])
        kol = gain(th, dd["kolf2"], [0])
        cm = float(np.mean(list(held.values())))
        out["lofo"][nm] = dict(heldout=held, kolf2=kol, theta=float(th[0]), cancer_mean=cm)
        print(f"{nm:20s} {100 * cm:+11.2f}%  " + " / ".join(f"{100 * held[h]:+.2f}" for h in CANCER) + f"   {100 * kol:+6.2f} %   {th[0]:+.3f}")
    (WORK / "strength518.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
