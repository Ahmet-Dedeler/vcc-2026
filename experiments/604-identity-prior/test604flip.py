#!/usr/bin/env python3
"""604 flip KAs (625's octile flips, [7772dc]'s spec, 2026-10-02 ~14:15Z). Light: A/B/C pinned profiles only.
  F1  build604 with no mask and with --flip 00000000 rebuild the COMMITTED v3 artifacts byte for byte (lam 0 / 0.5 / 1)
  F2  lam 0 + 11111111: pdown == 1 - P_C4 on every covered gene, psign == -CLEAN4 psign there (no P == 1/2 exists),
      uncovered genes untouched (CLEAN4 psign, P NaN)
  F3  lam 0 + one octile k (k = 0..7): exactly octile k's genes change (pdown and psign), every other gene bit-identical
  F4  lam 0.5 + a mask: the v3 side rule holds on the FLIPPED prior: sign == sign(zc' + 0.5 zi), zc' = (P' - 1/2)/sd,
      on every covered gene; uncovered untouched; the provenance records the mask and the module sha
  F5  invalid masks are refused
usage: uv run python experiments/604-identity-prior/test604flip.py
"""
from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import def_identity_prior as D  # noqa: E402

RES: list[tuple[bool, str]] = []


def rec(ok, name):
    RES.append((bool(ok), name))
    print(f"{'PASS' if ok else 'FAIL'} {name}", flush=True)


def main() -> int:
    meta = json.loads((HERE / "work/flex_refs604.npz.json").read_text())
    genes, rcpm, rcalls = D.load_refs(HERE / "work/flex_refs604.npz", meta["sha256"])
    P, ps = D.load_clean4(genes)
    cov = np.isfinite(P)
    oc = D.load_octiles(genes)
    ctx = {c: rcpm[c] for c in "ABC"}
    calls = {c: rcalls[c] for c in "ABC"}
    # F1: committed artifacts rebuilt byte for byte, with no mask and with 00000000
    committed = {lam: D.file_sha256(HERE / f"work/ident604v3_flex_lam{lam}_ABC.npz") for lam in ("0", "0.5", "1")}
    with tempfile.TemporaryDirectory(prefix="t604flip-") as td:
        for flag in ([], ["--flip", "00000000"]):
            r = subprocess.run([sys.executable, str(HERE / "build604.py"), "--lam", "0", "--lam", "0.5", "--lam", "1",
                                "--out-dir", td, *flag], capture_output=True, text=True)
            got = {lam: D.file_sha256(Path(td) / f"ident604v3_flex_lam{lam}_ABC.npz") for lam in ("0", "0.5", "1")}
            rec(r.returncode == 0 and got == committed,
                f"F1 build604 {' '.join(flag) or '(no mask)'} == committed v3 artifacts {[v[:8] for v in got.values()]}")
    # F2: lam 0, all octiles
    res = D.identity_prior(genes, ctx, calls, rcpm, rcalls, 0.0, P, ps, flip_octiles="11111111")
    for c in "ABC":
        r = res[c]
        ok = (np.array_equal(r["pdown"][cov], 1.0 - P[cov]) and np.array_equal(r["psign"][cov], -ps[cov])
              and np.array_equal(r["psign"][~cov], ps[~cov]) and np.isnan(r["pdown"][~cov]).all())
        rec(ok and int((P[cov] == 0.5).sum()) == 0, f"F2 {c}: lam 0 + 11111111 -> pdown = 1 - P_C4, psign negated on all "
            f"{int(cov.sum())} covered genes; uncovered untouched; prov {r['prov'].get('flip_octiles')}/{r['prov'].get('flipped_genes')}")
    # F3: lam 0, one octile at a time
    base = D.identity_prior(genes, ctx, calls, rcpm, rcalls, 0.0, P, ps)
    for k in range(8):
        mask = "".join("1" if j == k else "0" for j in range(8))
        rk = D.identity_prior(genes, ctx, calls, rcpm, rcalls, 0.0, P, ps, flip_octiles=mask)
        oks = []
        for c in "ABC":
            ch_p = ~((rk[c]["pdown"] == base[c]["pdown"]) | (np.isnan(rk[c]["pdown"]) & np.isnan(base[c]["pdown"])))
            ch_s = rk[c]["psign"] != base[c]["psign"]
            want = cov & (oc == k)
            oks.append(np.array_equal(ch_p, want) and np.array_equal(ch_s, want))
        rec(all(oks), f"F3 lam 0 + {mask}: exactly octile {k}'s {int((cov & (oc == k)).sum())} genes change (A/B/C)")
    # F4: lam 0.5 with a mask -> the v3 side rule on the flipped prior
    mask = "10010001"
    r5 = D.identity_prior(genes, ctx, calls, rcpm, rcalls, 0.5, P, ps, flip_octiles=mask)
    fv = D.flip_vector(genes, P, mask)
    Pf = np.where(fv, 1.0 - P, P)
    for c in "ABC":
        refs = r5[c]["prov"]["references"]
        m = np.isfinite(rcpm[c]).copy()
        for rr in refs:
            m &= np.isfinite(rcpm[rr])
        I = np.full(len(genes), np.nan)
        I[m] = np.log1p(rcpm[c][m]) - np.mean([np.log1p(rcpm[rr][m]) for rr in refs], 0)
        zi = (I - np.nanmean(I)) / np.nanstd(I)
        zc = (Pf - 0.5) / np.nanstd(P - 0.5)   # the UNFLIPPED prior's scale (14:55Z fix)
        sc = zc + 0.5 * np.where(np.isfinite(zi), zi, 0.0)
        side = np.where(sc > 0, -1.0, np.where(sc < 0, 1.0, 0.0))
        mm = cov & (side != 0)
        pr = r5[c]["prov"]
        ok = (np.array_equal(r5[c]["psign"][mm], side[mm]) and np.array_equal(r5[c]["psign"][~cov], ps[~cov])
              and pr.get("flip_octiles") == mask and pr.get("modules625_sha256") == D.MODULES625_SHA256
              and pr.get("flipped_genes") == int(fv.sum()))
        rec(ok, f"F4 {c}: lam 0.5 + {mask}: sign == sign(zc' + 0.5 zi) on the flipped prior ({int(mm.sum())} genes); "
            f"uncovered untouched; prov records mask + module sha")
    # F6 (607): a rehearsal reusing ids A/B/C with OTHER data and calls. Bundle B is an A/C mix called K562, so the pinned
    # B (HeLa, call HELA) survives the call guard. Pre-fix, the pinned B overwrote bundle B in the pool dict: B's own
    # deviation was HeLa's and the pin read similarity 1.0 to "itself" and was dropped. Fixed: it is pin:B, scored
    # against the A/C mix, and kept.
    ctx_x = dict(ctx)
    ctx_x["B"] = np.sqrt(np.asarray(rcpm["A"]) * np.asarray(rcpm["C"]))
    calls_x = dict(calls)
    calls_x["B"] = "K562_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"
    rx = D.identity_prior(genes, ctx_x, calls_x, rcpm, rcalls, 0.5, P, ps)
    pb = rx["B"]["prov"]
    sim = pb["sister_similarity"].get("pin:B")
    rec("pin:B" in pb["references"] and sim is not None and sim < 0.999 and "pin:B" not in pb["sister_dropped"],
        f"F6 name collision fixed: bundle B keeps pin:B (similarity {sim}); refs {pb['references']}")
    import importlib.util  # the pre-fix module (scratch copy of 13:40Z) must show the bug on the same input
    pre = Path("/private/tmp/claude-501/-Users-ahmet-Code-virtual-cell-challenge-2026/9b142876-d1c8-4ee2-80dd-f97f8de75a75/scratchpad/def_identity_prior.preflipfix.py")
    if pre.exists():
        spec = importlib.util.spec_from_file_location("dip_pre", pre)
        M = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(M)
        ro = M.identity_prior(genes, ctx_x, calls_x, rcpm, rcalls, 0.5, P, ps)["B"]["prov"]
        rec(ro["sister_similarity"].get("B") == 1.0 and "B" in ro["sister_dropped"],
            f"F6 the pre-fix module shows the bug on this input (B sim {ro['sister_similarity'].get('B')}, dropped {ro['sister_dropped']})")
    # F7 (14:55Z): at lam 0.5 a single-octile flip may change ONLY that octile's genes (zc keeps the unflipped scale)
    base5 = D.identity_prior(genes, ctx, calls, rcpm, rcalls, 0.5, P, ps)
    for k in (0, 3, 7):
        mask = "".join("1" if j == k else "0" for j in range(8))
        rk = D.identity_prior(genes, ctx, calls, rcpm, rcalls, 0.5, P, ps, flip_octiles=mask)
        ok, detail = True, []
        for c in "ABC":
            ch = rk[c]["psign"] != base5[c]["psign"]
            outside = int((ch & ~(cov & (oc == k))).sum())
            detail.append(f"{c}: {int(ch.sum())} changed, {outside} outside")
            ok &= outside == 0
        rec(ok, f"F7 lam 0.5 + {mask}: sides change only inside octile {k} ({'; '.join(detail)})")
    # F8 (15:05Z, 7772dc's option b): at lam 0.5 a flip leaves EVERY non-octile gene's pdown and psign byte-identical to the
    # base, and the octile genes keep the base's |P - 1/2| (only their sides move), so a tau > 0 hook sees nothing new
    for k in (0, 3, 7):
        mask = "".join("1" if j == k else "0" for j in range(8))
        rk = D.identity_prior(genes, ctx, calls, rcpm, rcalls, 0.5, P, ps, flip_octiles=mask)
        ok, detail = True, []
        inside = cov & (oc == k)
        for c in "ABC":
            pa, pb = rk[c]["pdown"], base5[c]["pdown"]
            same_out = np.array_equal(pa[~inside], pb[~inside], equal_nan=True) and np.array_equal(rk[c]["psign"][~inside], base5[c]["psign"][~inside])
            mag_in = np.allclose(np.abs(pa[inside] - 0.5), np.abs(pb[inside] - 0.5), rtol=0, atol=1e-15)
            detail.append(f"{c}: outside identical {same_out}, inside |P-1/2| kept {mag_in}")
            ok &= same_out and mag_in
        rec(ok, f"F8 lam 0.5 + {mask}: non-octile pdown/psign byte-identical to base; octile magnitudes kept ({'; '.join(detail)})")
    # F5: invalid masks
    for bad in ("0000000", "000000000", "00000002", "abcdefgh", "1111 111"):
        try:
            D.identity_prior(genes, ctx, calls, rcpm, rcalls, 0.0, P, ps, flip_octiles=bad)
            rec(False, f"F5 refuses {bad!r}")
        except SystemExit:
            rec(True, f"F5 refuses {bad!r}")
    n = sum(ok for ok, _ in RES)
    print(f"\n{n}/{len(RES)} PASS" if n == len(RES) else f"\n{len(RES) - n} FAILED of {len(RES)}")
    return 0 if n == len(RES) else 1


if __name__ == "__main__":
    sys.exit(main())
