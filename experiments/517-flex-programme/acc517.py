#!/usr/bin/env python3
"""517 — accuracy rows for every candidate sign source with 462f's EXACT definitions (known answers: CLEAN4 fid-like 0.5244,
Arc full 0.5545), and the sources saved to work/sources517.npz (H1 gated axis) for later scripts.

fid-like   : all block genes (2,200 top control-CPM kept genes, panel targets excluded; 462f's block), own gene out,
             lfc != 0, sign(source) vs sign(reference lfc), pooled over h1-test's 85 targets
gate-block : the same restricted to significant genes (q < 0.05)
gate-nonblk: significant genes outside that block (the hook's territory, approximately)
usage: uv run python experiments/517-flex-programme/acc517.py      (< 2 GB)
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import rank517 as RK  # noqa: E402
P = RK.P
RC = P.RC


def main():
    t0 = time.monotonic()
    zr = np.load(RC.REF, allow_pickle=False)
    genes_full = zr["genes"]; keep = np.asarray(zr["keep"], dtype=bool); T = zr["targets"].astype(str)
    lfc = np.asarray(zr["lfc"], float); q = np.asarray(zr["q"], float); bctrl = np.asarray(zr["b_ctrl"], float)
    zr.close()
    kg = np.asarray(genes_full)[keep].astype(str)
    e = np.expm1(bctrl[keep]); cpm = 1e6 * e / e.sum()
    tset = set(T)
    order = np.argsort(-cpm)
    block = np.zeros(kg.size, bool); block[[i for i in order if kg[i] not in tset][:2200]] = True
    SRC, _ = P.load_sources(genes_full, keep)
    Pc4, Parc = SRC["c4"], SRC["arc"]
    h1 = set(T)

    def bankP(b):
        delta, perts, genes, _ = RK.S.load_bank(b)
        return RK.S.programme(delta, perts, genes, list(map(str, kg)), h1)[1]
    with np.errstate(invalid="ignore"):
        Pxa = np.nanmean(np.vstack([bankP("orion_hct116"), bankP("orion_hek293t")]), 0)
    Psa = bankP("viperturb_k562")
    ax = RK.own_axes(keep, kg, Pc4)
    zI = -ax["ident"] / np.nanstd(ax["ident"])
    zP = ax["proj"] / np.std(ax["proj"][ax["proj"] != 0])
    Z4 = RK.zc(Pc4)
    src = {"c4": Pc4, "arc": Parc, "wx": RK.zmean(Z4, RK.zc(Pxa)), "wxs": RK.zmean(Z4, RK.zc(Pxa), RK.zc(Psa)),
           "ident": np.where(np.isfinite(zI), 0.5 - 0.5 * np.tanh(zI / 2), np.nan),
           "c4id1": RK.zsum(Z4, zI), "c4id2": RK.zsum(Z4, zI, w=[1, 2]), "own": RK.zsum(Z4, zP, zI),
           "down": np.full(Pc4.shape, 1.0)}
    np.savez_compressed(HERE / "work/sources517.npz", genes_keep=kg, ident=ax["ident"], proj=ax["proj"],
                        m_cpm_h1=ax["m_cpm"], **{f"P_{k}": v for k, v in src.items()})
    res = {}
    for name, Pv in src.items():
        s = np.where(np.isfinite(Pv), np.where(Pv > 0.5, -1.0, 1.0), np.nan)
        hit = {"fidlike": [0, 0], "gate_block": [0, 0], "gate_nonblock": [0, 0]}
        for t in range(len(T)):
            own = kg == T[t]
            base = ~own & np.isfinite(lfc[t]) & (lfc[t] != 0) & np.isfinite(s)
            for key, m in (("fidlike", base & block), ("gate_block", base & block & (q[t] < 0.05)),
                           ("gate_nonblock", base & ~block & (q[t] < 0.05))):
                hit[key][0] += int((np.sign(lfc[t][m]) == s[m]).sum()); hit[key][1] += int(m.sum())
        res[name] = {k: {"acc": round(v[0] / max(v[1], 1), 4), "n": v[1]} for k, v in hit.items()}
        print(f"{name:6s} " + "  ".join(f"{k} {v['acc']:.4f} (n {v['n']})" for k, v in res[name].items()), flush=True)
    ka = {"clean4_fidlike": res["c4"]["fidlike"]["acc"], "arc_fidlike": res["arc"]["fidlike"]["acc"]}
    ok = abs(ka["clean4_fidlike"] - 0.5244) < 1e-4 and abs(ka["arc_fidlike"] - 0.5545) < 1e-4
    print(f"KA (462f): CLEAN4 fid-like {ka['clean4_fidlike']} (want 0.5244), Arc {ka['arc_fidlike']} (want 0.5545): "
          f"{'PASS' if ok else 'FAIL'}", flush=True)
    (HERE / "work/acc517.json").write_text(json.dumps({"ka": ka, "ka_pass": ok, "rows": res}, indent=1))
    print(f"done {time.monotonic() - t0:.0f}s")


if __name__ == "__main__":
    main()
