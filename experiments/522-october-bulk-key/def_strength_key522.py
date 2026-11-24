"""522 — October's src-strength RIDER on the LODO key artifact. NOT APPLIED: a candidate day-of step for
scripts/def_day.sh, between `def_lodo_key.py` (which writes the LODO artifact) and the pin that registers it.

    w'_p = max(floor, w_p x exp(theta z_p))   for targets above the key floor (w_p > floor, default 0.5);
    w'_p = w_p                                  at the floor
    z_p  = mean, over the production donors (K562-GWPS, Orion HCT116, Orion HEK293T) that carry target p, of the centred
           percentile rank (over the panel's targets that donor carries) of 518's robust response strength
           (log1p #{genes > 10 CPM, off-panel, with |delta / (1.4826 MAD)| > 4 x the target's own noise}); 0 if none.
theta 0.25 (522): bounded multiplier exp(+-0.125) = 0.88 ... 1.13. The builder's apply_pert_scale reads only `targets`
and `w`, and restores the non-panel Frobenius norm, so only the relative weights change. The floor block is left intact
and nothing is pushed below it (522c): 479's P2g (`tercile:`) puts every target at the MINIMUM key ratio into LOW, so a
multiplier on the floor block, or a target pushed under it, would move the whole floor block (a third of an October
panel) to MID; S_P2 (`stercile:`) reads the unkeyed rows and is untouched either way.

Input: the npz `def_lodo_key.py` wrote (targets, w, key, carried, n_donors; any cap, e.g. --hi-pct 96.64697522769474).
Output: the same arrays with `w` replaced, plus strength_z, strength_mult, theta; a .json sidecar with the sha256.

KNOWN ANSWERS (--ka; runs on data/challenge; every failure exits 1):
  KA-z      z on the A/B/C panel == 518's abc518.zprod(targets, 'src') (strength518.npz route), bit for bit;
  KA-inert  theta 0 on 495's cap-6 artifact dc374281 -> targets, w, key, carried, n_donors np.array_equal to the input;
  KA-arm    theta 0.25 on dc374281 -> w == where(w > 0.5, max(0.5, dc374281.w x exp(0.25 z)), dc374281.w) exactly (the
            production form of 522c's priced arm `pctxS0.25fpc`: on A/B/C, P 96.64697522769474 == cap 6, so its w IS
            dc374281's); and min(w') == 0.5 (the floor block is still the minimum).

usage:
  uv run python experiments/522-october-bulk-key/def_strength_key522.py --in <lodo.npz> --bundle <dir> --theta 0.25 --out <npz>
  uv run python experiments/522-october-bulk-key/def_strength_key522.py --ka --out <scratch npz>
Peak RSS < 1.5 GB, a few seconds (reads ~300 rows of each donor bank).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))

import numpy as np  # noqa: E402

from strength522 import BANKS, panel_of, robust_strength  # noqa: E402

CAP6 = (ROOT / "experiments/495-lodo-cap6/work/pert_scale_415pp_lo0.5_hi6.npz",
        "dc374281452ed4f7ae3a25b392dbfa4ff0d5cc591bd7b5b71e2c0a2dca86a4c0")
FIELDS = ("targets", "w", "key", "carried", "n_donors")


def sha256(p) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def crank(v):
    """518's lofo518.crank: centred percentile rank in [-0.5, 0.5] over the finite entries; 0 elsewhere."""
    v = np.asarray(v, float)
    ok = np.isfinite(v)
    out = np.zeros(v.size)
    if ok.sum() > 1:
        r = np.argsort(np.argsort(v[ok]))
        out[ok] = r / (ok.sum() - 1) - 0.5
    return out


def zsrc(targets, panel_genes):
    """The production src feature: per donor, robust strength over `targets` (panel genes excluded from the gene set),
    centred rank over the targets it carries; mean over carrying donors; 0 where none carries the target."""
    cols, cover = [], {}
    for bk, path in BANKS.items():
        s, _, _ = robust_strength(path, list(targets), set(panel_genes))
        r = crank(s)
        r[~np.isfinite(s)] = np.nan
        cols.append(r)
        cover[bk] = int(np.isfinite(s).sum())
    with np.errstate(all="ignore"):
        m = np.nanmean(np.vstack(cols), 0)
    return np.where(np.isfinite(m), m, 0.0), cover


def rider(in_path: Path, bundle: Path, theta: float, out: Path, floor: float = 0.5) -> dict:
    z_in = np.load(in_path, allow_pickle=False)
    arrs = {k: np.asarray(z_in[k]) for k in FIELDS}
    targets = [str(t) for t in arrs["targets"]]
    panel = panel_of(bundle)
    assert sorted(targets) == sorted(panel) and len(targets) == 300, "artifact targets != the bundle's panel"
    z, cover = zsrc(targets, panel)
    w0 = np.asarray(arrs["w"], np.float64)
    assert float(w0.min()) >= floor - 1e-12, f"weights below the floor {floor}: min {w0.min()}"
    at_floor = w0 <= floor + 1e-12
    m = np.where(at_floor, 1.0, np.exp(theta * z))
    w = np.where(at_floor, w0, np.maximum(floor, w0 * m))
    assert np.all(np.isfinite(w)) and np.all(w > 0), "rider weights must be finite and > 0"
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez(out, **{**arrs, "w": w}, strength_z=z, strength_mult=m, theta=np.float64(theta))
    q = lambda x: np.round(np.quantile(x, [0, .1, .25, .5, .75, .9, 1]), 4).tolist()
    side = dict(sha256=sha256(out), input=str(in_path), input_sha256=sha256(in_path), bundle=str(bundle), theta=theta,
                floor=floor, n_at_floor_kept=int(at_floor.sum()),
                n_clipped_to_floor=int(np.sum(~at_floor & (w0 * m < floor))),
                rule="w' = max(floor, w x exp(theta z_src)) above the floor, w at the floor; z = mean over carrying "
                     "production donors of crank(robust strength)",
                donor_coverage=cover, n_z_zero=int((z == 0).sum()), z_quantiles=q(z), mult_quantiles=q(m),
                w_in_quantiles=q(arrs["w"]), w_out_quantiles=q(w))
    out.with_suffix(".json").write_text(json.dumps(side, indent=1))
    return side


def ka(out: Path) -> None:
    ok = True
    assert sha256(CAP6[0]) == CAP6[1], "KA reference artifact sha (495's dc374281)"
    ref = np.load(CAP6[0], allow_pickle=False)
    targets = [str(t) for t in ref["targets"]]
    panel = panel_of(ROOT / "data/challenge")
    z, _ = zsrc(targets, panel)
    sys.path.insert(0, str(ROOT / "experiments/518-bulk-mse-gap"))
    from abc518 import zprod           # 518's production src feature (strength518.npz)
    z518 = zprod(targets, "src")
    okz = np.array_equal(z, z518)
    print(f"KA-z     z == abc518.zprod(targets, 'src'): {'PASS' if okz else 'FAIL'} (max |d| {np.abs(z - z518).max():.1e})")
    ok &= okz
    s0 = rider(CAP6[0], ROOT / "data/challenge", 0.0, out)
    z0 = np.load(out, allow_pickle=False)
    oki = all(np.array_equal(np.asarray(z0[k]), np.asarray(ref[k])) for k in FIELDS)
    print(f"KA-inert theta 0 on dc374281: arrays identical {'PASS' if oki else 'FAIL'} ({s0['sha256'][:12]})")
    ok &= oki
    s1 = rider(CAP6[0], ROOT / "data/challenge", 0.25, out)
    z1 = np.load(out, allow_pickle=False)
    w_ref = np.asarray(ref["w"], np.float64)
    want = np.where(w_ref <= 0.5 + 1e-12, w_ref, np.maximum(0.5, w_ref * np.exp(0.25 * z518)))
    oka = np.array_equal(np.asarray(z1["w"]), want) and float(np.asarray(z1["w"]).min()) == 0.5 and \
        all(np.array_equal(np.asarray(z1[k]), np.asarray(ref[k])) for k in FIELDS if k != "w")
    print(f"KA-arm   theta 0.25 on dc374281: w == w_cap6 x exp(0.25 z) above the floor {'PASS' if oka else 'FAIL'} ({s1['sha256'][:12]}); "
          f"mult q {s1['mult_quantiles']}; w' q {s1['w_out_quantiles']}")
    ok &= oka
    if not ok:
        raise SystemExit(1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="inp", type=Path)
    ap.add_argument("--bundle", type=Path)
    ap.add_argument("--theta", type=float, default=0.25)
    ap.add_argument("--floor", type=float, default=0.5, help="the LODO key's floor (def_lodo_key.py --lo); kept untouched")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--ka", action="store_true")
    a = ap.parse_args()
    if a.ka:
        ka(a.out)
        return
    assert a.inp is not None and a.bundle is not None, "--in and --bundle are required"
    side = rider(a.inp.resolve(), a.bundle.resolve(), a.theta, a.out, a.floor)
    print(json.dumps(side, indent=1))


if __name__ == "__main__":
    main()
