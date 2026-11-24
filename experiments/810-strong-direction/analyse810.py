"""810 -- the decomposition table from decomp810's per-fold records (work/decomp/*-comp.npz, *-real.npz; ~0.2 MB each).
Light (local, < 1 s per fold). Composition level = the incumbent's composition rows (before key / sharpening / render);
realised level = engine801's incumbent rows (after key, P2g sharpening, render, ballast; frame-averaged).

Per fold and z_src tertile (0 weak, 1 mid, 2 strong; fold522.zsrc's rule over the fold's targets):
  cos            per-target cosine with the held line's truth T on its scored genes (median and mean over targets)
  contrib_k      <part_k, T> / (|fin| |T|): the additive split of cos(fin, T) into base / cd4 / h1 / part
  energy_k       |part_k|^2 / |fin|^2 (cross terms make them not sum to 1)
  X-weighted     tier aggregate X_t = sum_i x_i / sqrt(sum a sum d) with the realised x, a, d (the instrument's cos_agg
                 restricted to the tier), and the same with composition-level cos at the realised norms (direction only)
  ablations      direction-only replacements at the same per-target norm: drop one part, base only, U only, pre only,
                 the 707 blend (beta 0.25 toward U), src_j alone, and the in-sample per-target LS oracle in the span of
                 the parts (+U, pre, sources): an upper bound on how much direction those vectors carry
  donors         cos(src0, src1) (processed rows), cos(raw0, raw1), cos(src_j, T), each source's and the held line's
                 knockdown of the target, robust strength count / noise
    uv run python experiments/810-strong-direction/analyse810.py [--json OUT]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
D = HERE / "work" / "decomp"
PARTS4 = ["base", "cd4", "h1", "part"]


def load(stem):
    c = dict(np.load(D / f"{stem}-comp.npz", allow_pickle=False))
    r = dict(np.load(D / f"{stem}-real.npz", allow_pickle=False))
    assert list(c["targets"]) == list(r["targets"]), stem
    c["names"] = [str(x) for x in c["names"]]
    c["sources"] = [str(x) for x in c["sources"]]
    c.update({f"r_{k}": v for k, v in r.items() if k != "targets"})
    c["meta"] = json.loads((D / f"{stem}-comp.json").read_text())
    return c


def idx(c, k):
    return c["names"].index(k)


def comb_cos(c, alpha):
    """per-target cosine with T of v_i = sum_k alpha[i, k] V_ik on the scored genes (alpha n x K)."""
    G, b, tt = c["G"], c["bT"], c["tt"]
    num = np.einsum("nk,nk->n", alpha, b)
    nv = np.sqrt(np.maximum(np.einsum("nk,nkl,nl->n", alpha, G, alpha), 0))
    return num / np.maximum(nv * np.sqrt(tt), 1e-300), nv


def per_target(c):
    n, K = c["bT"].shape
    out = {}
    a_fin = np.zeros((n, K))
    for k in PARTS4:
        a_fin[:, idx(c, k)] = 1.0
    cos_fin, nfin = comb_cos(c, a_fin)
    out["cos_fin"], out["nfin"] = cos_fin, nfin
    nt = np.sqrt(c["tt"])
    for k in PARTS4:
        j = idx(c, k)
        out[f"contrib_{k}"] = c["bT"][:, j] / np.maximum(nfin * nt, 1e-300)
        out[f"energy_{k}"] = c["G"][:, j, j] / np.maximum(nfin ** 2, 1e-300)
        out[f"cos_{k}"] = np.where(c["G"][:, j, j] > 0, c["bT"][:, j] / np.maximum(np.sqrt(c["G"][:, j, j]) * nt, 1e-300), np.nan)
    for k in ["U", "pre"] + [f"src{j}" for j in range(len(c["sources"]))] + [f"raw{j}" for j in range(len(c["sources"]))]:
        j = idx(c, k)
        out[f"cos_{k}"] = np.where(c["G"][:, j, j] > 0, c["bT"][:, j] / np.maximum(np.sqrt(c["G"][:, j, j]) * nt, 1e-300), np.nan)
    # ablations (direction only)
    for k in PARTS4:
        a = a_fin.copy()
        a[:, idx(c, k)] = 0.0
        out[f"abl_no_{k}"] = comb_cos(c, a)[0]
    # 707's blend: |fin| normalise((1-b) fin/|fin| + b U/|U|) on the NON-PANEL axis (the norms there: Gfull)
    jU = idx(c, "U")
    Gf = c["Gfull"]
    nfin_full = np.sqrt(np.maximum(np.einsum("nk,nkl,nl->n", a_fin, Gf, a_fin), 0))
    nU_full = np.sqrt(Gf[:, jU, jU])
    for beta in (0.25, 0.5, 1.0):
        a = a_fin * ((1 - beta) / np.maximum(nfin_full, 1e-300))[:, None]
        a[:, jU] += np.where(nU_full > 0, beta / np.maximum(nU_full, 1e-300), 0.0)
        a[nU_full == 0] = a_fin[nU_full == 0]
        out[f"blend_U_b{beta:g}"] = comb_cos(c, a)[0]
    # the base alone re-directed toward the unwhitened donor mean `pre` (residual parts kept): base' = |base| normalise(
    # (1-b) base/|base| + b pre/|pre|) on the non-panel axis; then direction-only as above
    jb, jp = idx(c, "base"), idx(c, "pre")
    nb, npre = np.sqrt(Gf[:, jb, jb]), np.sqrt(Gf[:, jp, jp])
    for beta in (0.5, 1.0):
        cb = (1 - beta) / np.maximum(nb, 1e-300)
        cp = np.where(npre > 0, beta / np.maximum(npre, 1e-300), 0.0)
        nm = np.sqrt(np.maximum(cb ** 2 * Gf[:, jb, jb] + 2 * cb * cp * Gf[:, jb, jp] + cp ** 2 * Gf[:, jp, jp], 0))
        a = a_fin.copy()
        a[:, jb] = np.where(npre > 0, cb * nb / np.maximum(nm, 1e-300), 1.0)
        a[:, jp] = np.where(npre > 0, cp * nb / np.maximum(nm, 1e-300), 0.0)
        out[f"swapbase_pre_b{beta:g}"] = comb_cos(c, a)[0]
    # in-sample LS oracle in a span (least squares of T on the vectors, per target): cos = sqrt(b' G^+ b / tt)
    def oracle(keys):
        js = [idx(c, k) for k in keys]
        res = np.zeros(n)
        for i in range(n):
            Gi = c["G"][i][np.ix_(js, js)]
            bi = c["bT"][i, js]
            if not np.any(np.diag(Gi) > 0):
                res[i] = np.nan
                continue
            res[i] = np.sqrt(max(float(bi @ np.linalg.pinv(Gi, rcond=1e-10) @ bi), 0.0) / c["tt"][i])
        return res
    out["oracle_parts"] = oracle(PARTS4)
    srcs = [f"src{j}" for j in range(len(c["sources"]))]
    out["oracle_parts_U_src"] = oracle(PARTS4 + ["U", "pre"] + srcs)
    if len(srcs) >= 2:
        j0, j1 = idx(c, "src0"), idx(c, "src1")
        out["cos_src01"] = c["G"][:, j0, j1] / np.maximum(np.sqrt(c["G"][:, j0, j0] * c["G"][:, j1, j1]), 1e-300)
        r0, r1 = idx(c, "raw0"), idx(c, "raw1")
        out["cos_raw01"] = c["Gfull"][:, r0, r1] / np.maximum(np.sqrt(c["Gfull"][:, r0, r0] * c["Gfull"][:, r1, r1]), 1e-300)
        out["cos_fin_U"] = (np.einsum("nk,nk->n", a_fin, Gf[:, :, jU])) / np.maximum(nfin_full * nU_full, 1e-300)
        jb = idx(c, "base")
        out["cos_base_U"] = Gf[:, jb, jU] / np.maximum(np.sqrt(Gf[:, jb, jb]) * nU_full, 1e-300)
    # realised
    x, a, d = c["r_x"], c["r_a"], c["r_d"]
    out["cos_real"] = np.mean(x / np.maximum(np.sqrt(a * d[None, :]), 1e-300), 0)
    out["ra"] = a.mean(0)
    out["rd"] = d
    out["pds"] = c["r_pds"].mean(0)
    return out


def summarise(c, pt):
    tier = c["tier"]
    live = pt["rd"] > 0
    S = {}
    for t in (0, 1, 2):
        m = (tier == t) & live
        if m.sum() == 0:
            continue
        sa, sd = np.sqrt(pt["ra"][m]), np.sqrt(pt["rd"][m])
        den = np.sqrt((pt["ra"][m]).sum() * (pt["rd"][m]).sum())

        def X(cos):                       # direction-only aggregate at the realised norms
            v = np.nan_to_num(cos[m])
            return float((sa * sd * v).sum() / den)
        row = dict(n=int(m.sum()), X_real=float(np.mean([(c["r_x"][f, m]).sum() / np.sqrt(c["r_a"][f, m].sum() * c["r_d"][m].sum())
                                                         for f in range(c["r_x"].shape[0])])),
                   share_a=float(pt["ra"][m].sum() / pt["ra"][live].sum()), share_d=float(pt["rd"][m].sum() / pt["rd"][live].sum()),
                   med_cos_real=float(np.median(pt["cos_real"][m])), med_cos_fin=float(np.median(pt["cos_fin"][m])),
                   Xc_fin=X(pt["cos_fin"]), pds=float(np.mean(pt["pds"][m])))
        for k in PARTS4:
            row[f"contrib_{k}"] = float(np.mean(pt[f"contrib_{k}"][m]))
            row[f"energy_{k}"] = float(np.mean(pt[f"energy_{k}"][m]))
            row[f"med_cos_{k}"] = float(np.nanmedian(pt[f"cos_{k}"][m])) if np.isfinite(pt[f"cos_{k}"][m]).any() else None
            row[f"carry_{k}"] = float(np.mean(pt[f"energy_{k}"][m] > 1e-6))
            row[f"dXc_no_{k}"] = X(pt[f"abl_no_{k}"]) - row["Xc_fin"]
        for k in ["U", "pre"] + [f"src{j}" for j in range(len(c["sources"]))]:
            row[f"med_cos_{k}"] = float(np.nanmedian(pt[f"cos_{k}"][m]))
            row[f"dXc_only_{k}"] = X(pt[f"cos_{k}"]) - row["Xc_fin"]
        for beta in (0.25, 0.5, 1.0):
            row[f"dXc_blendU_b{beta:g}"] = X(pt[f"blend_U_b{beta:g}"]) - row["Xc_fin"]
        for beta in (0.5, 1.0):
            row[f"dXc_swapbase_pre_b{beta:g}"] = X(pt[f"swapbase_pre_b{beta:g}"]) - row["Xc_fin"]
        row["dXc_oracle_parts"] = X(pt["oracle_parts"]) - row["Xc_fin"]
        row["dXc_oracle_all"] = X(pt["oracle_parts_U_src"]) - row["Xc_fin"]
        if "cos_src01" in pt:
            row["med_cos_src01"] = float(np.nanmedian(pt["cos_src01"][m]))
            row["med_cos_raw01"] = float(np.nanmedian(pt["cos_raw01"][m]))
            row["med_cos_fin_U"] = float(np.nanmedian(pt["cos_fin_U"][m]))
            row["med_cos_base_U"] = float(np.nanmedian(pt["cos_base_U"][m]))
        kd = np.vstack([c[f"kd_{b}"] for b in c["sources"]])
        row["med_kd_src"] = float(np.nanmedian(kd[:, m]))
        held = c["meta"]["fold"].split("/")[1]
        if f"kd_{held}" in c:
            row["med_kd_held"] = float(np.nanmedian(c[f"kd_{held}"][m]))
        S[t] = row
    return S


def main():
    stems = sorted(p.name[:-9] for p in D.glob("*-comp.npz"))
    res, PT = {}, {}
    for st in stems:
        if not (D / f"{st}-real.npz").exists():
            continue
        c = load(st)
        pt = per_target(c)
        PT[st] = (c, pt)
        res[st] = summarise(c, pt)
    if "--json" in sys.argv:
        Path(sys.argv[sys.argv.index("--json") + 1]).write_text(json.dumps(res, indent=1))
    return res, PT


if __name__ == "__main__":
    res, _ = main()
    for st, S in res.items():
        for t, r in S.items():
            print(st, t, json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()}))
