"""480 step 2 on h1-test: pack-style and panel-donor constructions against ours, analytic nmae on the real gate.

(a)  ours: 426 base_s0_bps (known answer 0.95237) and the shipped coarse MDE (435, 0.94270)
(b)  K562-GWPS same-target response, scale 1, log2 pseudobulk ratio on the truth axis; 0 where no row / no gene
(c)  Orion HCT116, Orion HEK293T, and their consensus with K562 (mean of available), scale 1 and fitted global scale
(d)  our sign x donor magnitude
(e)  source-family-balanced consensus: families = Weissman-K562 (GWPS), Xaira-Orion (HCT116+HEK293T averaged first)
Conversions: pseudocount 0 / 1 / 5 CPM, and the naive delta/ln2 (= 20 CPM). Never uses data/donors/h1_2025*.
Writes work/h1con480.npz (per-target rows) and work/h1con480.json."""
import json

import numpy as np

from common480 import HERE, ROOT, boot_delta, load_donor_rows, load_h1, nmae, nmae_rows

h = load_h1()
T, cnt, genes, targets = h["T"], h["cnt"], h["genes"], h["targets"]
scored = np.isfinite(T).any(1)
L0 = np.nan_to_num(np.asarray(np.load(ROOT / "experiments/426-mde-magnitude/work/h1/stats-base_s0_bps.npz")["lfc"],
                              float))
Lc = np.nan_to_num(np.asarray(np.load(ROOT / "experiments/435-zs-on-coarse/work/h1/"
                                      "stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.npz")["lfc"], float))
rows = {"a_base": nmae_rows(L0, T), "a_coarse": nmae_rows(Lc, T)}
print(f"a_base {nmae(L0, T):.5f}  a_coarse {nmae(Lc, T):.5f}")

# essential-screen coverage -> panel-like subset (430/434's 38 targets)
ess = set()
for d in ["k562_essential", "rpe1", "nadig_jurkat"]:
    ess |= {l.strip() for l in open(ROOT / f"data/donors/{d}/perts.txt")}
panel_like = np.array([t not in ess for t in targets])
print("scored targets", scored.sum(), "panel-like scored", (scored & panel_like).sum())


def l1_scale(P, T):
    """Global scale k minimising the target-mean nmae of k*P (grid then refine)."""
    ks = np.concatenate([np.linspace(0, 2, 81)])
    v = [nmae(k * P, T) for k in ks]
    k0 = ks[int(np.argmin(v))]
    ks2 = np.linspace(max(0, k0 - 0.05), k0 + 0.05, 41)
    v2 = [nmae(k * P, T) for k in ks2]
    return float(ks2[int(np.argmin(v2))]), float(min(v2))


out = {}
cov = {}
DON = ["k562_gwps", "orion_hct116", "orion_hek293t"]
for conv, pc, clip in [("pc0", 0.0, 10.0), ("pc1", 1.0, None), ("pc5", 5.0, None), ("naive20", 20.0, None)]:
    R = {}
    for d in DON:
        L2, has, okc, _ = load_donor_rows(d, targets, genes, pc_cpm=pc, clip=clip)
        R[d] = L2
        cov[d] = dict(targets=int((has & scored).sum()), genes_on_axis=int(okc.sum()),
                      gate_cells_covered=int((np.isfinite(L2) & np.isfinite(T)).sum()))
    stack = np.stack([R[d] for d in DON])
    cons = np.nanmean(stack, 0)                                   # mean of available donors
    orion = np.nanmean(np.stack([R["orion_hct116"], R["orion_hek293t"]]), 0)
    fam = np.nanmean(np.stack([R["k562_gwps"], orion]), 0)       # family-balanced: Weissman vs Xaira
    cands = {"b_k562": R["k562_gwps"], "c_hct116": R["orion_hct116"], "c_hek293t": R["orion_hek293t"],
             "c_orion2": orion, "c_cons3": cons, "e_fam2": fam}
    for name, P in cands.items():
        P0 = np.nan_to_num(P)
        key = f"{conv}:{name}"
        rows[key] = nmae_rows(P0, T)
        k, v = l1_scale(P0, T)
        rows[key + "@k"] = nmae_rows(k * P0, T)
        # (d) our sign (coarse arm), donor magnitude; zero where our sign is 0
        sg = np.sign(Lc)
        Pd = sg * np.abs(P0)
        rows[key + ":oursign"] = nmae_rows(Pd, T)
        kd, vd = l1_scale(Pd, T)
        rows[key + ":oursign@k"] = nmae_rows(kd * Pd, T)
        out[key] = dict(scale1=float(np.nanmean(rows[key])), k=k, at_k=v,
                        oursign_scale1=float(np.nanmean(rows[key + ":oursign"])), oursign_k=kd, oursign_at_k=vd)
        print(f"{key:22s} scale1 {out[key]['scale1']:.4f} | k {k:.3f} -> {v:.4f} | oursign s1 "
              f"{out[key]['oursign_scale1']:.4f} k {kd:.3f} -> {vd:.4f}")

print("coverage", json.dumps(cov))
np.savez_compressed(HERE / "work/h1con480.npz", targets=targets, panel_like=panel_like, scored=scored,
                    **{k.replace(":", "__").replace("@", "_at_"): v for k, v in rows.items()})
# bootstrap: each construction vs a_coarse, all scored and panel-like
summ = {}
for key, v in rows.items():
    if key.startswith("a_"):
        continue
    d_all = boot_delta(rows["a_coarse"], v)
    d_pl = boot_delta(rows["a_coarse"], v, sub=panel_like)
    summ[key] = dict(all=d_all, panel_like=d_pl)
json.dump(dict(out=out, cov=cov, boot=summ, a_base=float(np.nanmean(rows["a_base"])),
               a_coarse=float(np.nanmean(rows["a_coarse"]))), open(HERE / "work/h1con480.json", "w"), indent=1)
print("\nDelta vs a_coarse (target bootstrap 4000): all | panel-like")
for key in ["pc0:b_k562", "pc1:b_k562", "pc5:b_k562", "naive20:b_k562", "pc1:c_cons3", "pc1:e_fam2",
            "pc1:b_k562:oursign", "pc1:c_cons3:oursign"]:
    a, p = summ[key]["all"], summ[key]["panel_like"]
    print(f"{key:24s} {a[0]:+.4f} [{a[1]:+.4f},{a[2]:+.4f}] | {p[0]:+.4f} [{p[1]:+.4f},{p[2]:+.4f}]")
