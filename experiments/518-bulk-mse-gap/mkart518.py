"""518 — contingency pert-scale artifacts for the STRENGTH KEY (written here only; nothing is registered or applied).

w = clip(key / median(key[carried]), 0.5, cap) on carried targets (1 elsewhere) x exp(theta z_p), with key = 8c82bef3's LODO
key vector and z_p = the centred percentile rank of the robust response strength (strength518) averaged over the screens
that carry target p: K562-GWPS, Orion HCT116, Orion HEK293T (the production sources), CD4 T cells (DESeq2 nDE, Marson
GWCD4i), VIPerturb-seq K562, the KOLF2.1J atlas panel rows ("meta"); 0 where none carries it.
  S1 = cap 3, theta 1.25  (the stand-alone strength key; priced +0.0016 ... +0.0024 overall vs 446bj2)
  S2 = cap 6, theta 0.5   (the rider on 495; priced +0.0001 ... +0.0003 over cap 6)
KNOWN ANSWERS (abort): theta 0 at cap 3 reproduces 8c82bef3's w bit for bit; theta 0 at cap 6 reproduces 495's dc374281's w;
the multiplier equals abc518's (the priced arm) exactly.

    uv run python experiments/518-bulk-mse-gap/mkart518.py
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import hashlib  # noqa: E402
import json  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(HERE))
from abc518 import zprod  # noqa: E402

KEY = ROOT / "experiments/420-bulk-stack/work/pert_scale_415pp_lo0.5_hi3.npz"
KEY_SHA = "8c82bef37d0a59b3ed1f4917ccad0d8fd24618a89fd9e2700c058b9aeff4339b"
CAP6 = ROOT / "experiments/495-lodo-cap6/work/pert_scale_415pp_lo0.5_hi6.npz"
CAP6_SHA = "dc374281452ed4f7ae3a25b392dbfa4ff0d5cc591bd7b5b71e2c0a2dca86a4c0"


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    assert sha(KEY) == KEY_SHA
    zk = np.load(KEY, allow_pickle=False)
    targets = [str(t) for t in zk["targets"]]
    key, car = np.asarray(zk["key"], np.float64), np.asarray(zk["carried"], bool)
    med = float(np.median(key[car]))
    w_of = lambda hi: np.where(car, np.clip(key / med, 0.5, hi), 1.0)
    assert np.array_equal(w_of(3.0), np.asarray(zk["w"], np.float64)), "KA cap 3"
    if CAP6.exists():
        assert sha(CAP6) == CAP6_SHA
        assert np.array_equal(w_of(6.0), np.asarray(np.load(CAP6)["w"], np.float64)), "KA cap 6 vs 495's artifact"
        ka6 = "PASS (bit-identical to 495's dc374281 w)"
    else:
        ka6 = "495 artifact absent"
    z = zprod(targets, "meta")
    out = {}
    for tag, hi, th in (("S1", 3.0, 1.25), ("S2", 6.0, 0.5)):
        m = np.exp(th * z)
        w = w_of(hi) * m
        path = WORK / f"pert_scale_518_{tag}_strengthmeta_th{th:g}_lo0.5_hi{hi:g}.npz"
        np.savez(path, targets=zk["targets"], w=w, key=zk["key"], carried=zk["carried"], n_donors=zk["n_donors"],
                 strength_z=z, strength_mult=m)
        s = sha(path)
        meta = dict(artifact=str(path.relative_to(ROOT)), sha256=s, base_key=str(KEY.relative_to(ROOT)), base_key_sha256=KEY_SHA,
                    rule=f"w = clip(key/median(key[carried]), 0.5, {hi:g}) on carried, 1 elsewhere; x exp({th:g} z)",
                    z="centred percentile rank of robust strength, mean over carrying screens (strength518.npz); 0 if none",
                    screens=["k562_gwps", "orion_hct116", "orion_hek293t", "cd4_deseq", "viperturb_k562", "kolf2"],
                    strength_npz_sha256=sha(WORK / "strength518.npz"),
                    w_quantiles=np.quantile(w, [0, .1, .25, .5, .75, .9, 1]).tolist(),
                    mult_quantiles=np.quantile(m, [0, .1, .5, .9, 1]).tolist(), n_z_zero=int((z == 0).sum()))
        path.with_suffix(".json").write_text(json.dumps(meta, indent=1))
        out[tag] = meta
        print(f"{tag}: {path.name} sha {s[:16]}  w q {np.round(meta['w_quantiles'], 3).tolist()}  mult {np.round(meta['mult_quantiles'], 3).tolist()}")
    print(f"KA cap 3: PASS (bit-identical to 8c82bef3's w); KA cap 6: {ka6}")
    (WORK / "mkart518.json").write_text(json.dumps(dict(ka_cap3="PASS", ka_cap6=ka6, artifacts=out), indent=1))


if __name__ == "__main__":
    main()
