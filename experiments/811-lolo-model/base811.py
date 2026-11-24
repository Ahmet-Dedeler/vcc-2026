"""811 BASELINE on one bulk fold (EC2 job): the equal-weight consensus of all NON-HELD public lines' unit responses blended
into the composition, its sources-only control, and a same-line ceiling diagnostic. Frozen in PREREGISTRATION.md.

Arms (STANDARD FORMAT npz, priced by engine801.price against the incumbent):
  inc                     the fold's incumbent (ka804's construction)           KA: == 804's inc.npz (every array)
  707                     ka804's 707 rebuilt from ka804's own functions        KA: == 804's 707.npz (every array)
  all_l10/l25/l50         U_all (every non-held LINE, equal weight per line), all targets carrying it, lambda 0.10/0.25/0.50
  src_l10/l25/l50         U_src (the fold's composition sources only, = 707's U), all targets, lambda 0.10/0.25/0.50
  oracle_vip (K562 folds) lambda 1 swap to viperturb_k562's own rows: NOT shippable (same line as the truth), a ceiling probe
Blend (707's formula on the non-panel genes): fin'_t = |fin_t| normalise((1-l) fin_t/|fin_t| + l U_t/|U_t|).

Inputs: engine801 folds (801 registry / 518 dumps), data/donors/{k562_gwps,orion_hct116,orion_hek293t} (the composition
donors, read as ka804 reads them) and 811's compact banks (work/banks/*.npz from ext811.py, sha256 in its manifest).

    uv run python experiments/811-lolo-model/base811.py <fold-id>        # writes $BASE811_OUT/{bulk,meta}
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/801-lab-folds"))
sys.path.insert(0, str(ROOT / "experiments/522-october-bulk-key"))
sys.path.insert(0, str(ROOT / "experiments/804-known-answers"))
import numpy as np  # noqa: E402

import engine801 as E  # noqa: E402  (guard801: writes inside the repo outside 801 raise; we write under BASE811_OUT)
import ka804 as K  # noqa: E402  (804's arm constructions, reused verbatim for the KA rows)

OUT = Path(os.environ.get("BASE811_OUT", "/work/base811"))
BANKS = HERE / "work" / "banks"
REF804 = HERE / "work" / "ref804"   # byte copies (cp -c) of 804's run/bulk/<fold>/{inc,707}.npz, SHA256SUMS.txt
GENES = K.GENES
COMPOSITION = ("k562_gwps", "orion_hct116", "orion_hek293t")
GROUPS = {  # line -> banks (PREREGISTRATION.md, frozen)
    "K562": ("k562_gwps", "viperturb_k562", "jiang_K562"),
    "HCT116": ("orion_hct116",), "HEK293T": ("orion_hek293t",),
    "CD4": ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr"),
    "H1": ("h1_pan2025", "h1_vcc2025"), "KOLF2": ("kolf2",), "TeloHAEC": ("telohaec_pilot",), "DLD1": ("dld1_lowmoi",),
    "A549": ("altair_a549", "jiang_A549"), "RPE1": ("rpe1",), "Jurkat": ("nadig_jurkat",),
    "BXPC3": ("jiang_BXPC3",), "HAP1": ("jiang_HAP1",), "HT29": ("jiang_HT29",), "MCF7": ("jiang_MCF7",),
}
HELD_GROUP = {"k562_gwps": "K562", "orion_hct116": "HCT116", "orion_hek293t": "HEK293T", "kolf2": "KOLF2"}
LAMS = (0.10, 0.25, 0.50)
TOL = 1e-12


def sha(p) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def bank_rows(bank: str, targets: list[str]) -> np.ndarray:
    """(n x 18533) raw rows of `bank` for `targets` (0 where absent), official axis."""
    R = np.zeros((len(targets), len(GENES)))
    if bank in COMPOSITION:
        D = ROOT / "data/donors" / bank
        dg = (D / "genes.txt").read_text().split()
        dp = (D / "perts.txt").read_text().split()
        delta = np.load(D / "delta.npy", mmap_mode="r")
        ppos = {p: i for i, p in enumerate(dp)}
        gpos = {g: i for i, g in enumerate(GENES)}
        cols = np.array([gpos.get(g, -1) for g in dg])
        ok = cols >= 0
        for i, t in enumerate(targets):
            if t in ppos:
                R[i, cols[ok]] = np.asarray(delta[ppos[t]], np.float64)[ok]     # ka804.donor_unit_mean's mapping
        return R
    z = np.load(BANKS / f"{bank}.npz")
    pos = {str(t): i for i, t in enumerate(z["targets"])}
    rows = z["rows"]
    for i, t in enumerate(targets):
        if t in pos:
            R[i] = np.asarray(rows[pos[t]], np.float64)
    return R


def unit(R: np.ndarray, nonpanel: np.ndarray) -> np.ndarray:
    R = R.copy()
    R[:, ~nonpanel] = 0.0
    n = np.linalg.norm(R, axis=1, keepdims=True)
    return np.where(n > 0, R / np.maximum(n, 1e-300), 0.0)


def consensus(F, groups: dict) -> tuple[np.ndarray, dict]:
    """equal weight per LINE: line row = normalise(mean of its banks' unit rows carrying t); U = mean over carrying lines."""
    acc = np.zeros((F.n, len(GENES)))
    nlines = np.zeros(F.n, int)
    cov = {}
    for line, banks in groups.items():
        S = np.zeros((F.n, len(GENES)))
        for b in banks:
            S += unit(bank_rows(b, list(F.targets)), F.nonpanel)
        Lr = unit(S, F.nonpanel)
        has = np.linalg.norm(Lr, axis=1) > 0
        cov[line] = int(has.sum())
        acc += Lr
        nlines += has
    U = np.where(nlines[:, None] > 0, acc / np.maximum(nlines, 1)[:, None], 0.0)
    return U, dict(lines_per_target_q=np.quantile(nlines, [0, .25, .5, .75, 1]).tolist(), coverage=cov)


def blend(fin, U, lam, npn, rows=None):
    """707's norm-preserving blend on the non-panel genes (ka804 lines 'for i in np.flatnonzero(sel)')."""
    out = fin.copy()
    sel = np.ones(fin.shape[0], bool) if rows is None else rows
    for i in np.flatnonzero(sel):
        nb, nu = np.linalg.norm(fin[i, npn]), np.linalg.norm(U[i, npn])
        if nb > 0 and nu > 0:
            m = (1 - lam) * fin[i, npn] / nb + lam * U[i, npn] / nu
            out[i, npn] = nb * m / np.linalg.norm(m)
    return out


def compare(fold_id, name, arr, ref_dir):
    """KA: every array 804 stored for this arm must equal ours (|d| <= 1e-12; targets exactly)."""
    p = ref_dir / f"{name}.npz"
    z = dict(np.load(p, allow_pickle=True))
    worst = 0.0
    for k, v in z.items():
        if k == "targets":
            assert list(map(str, v)) == list(map(str, arr["targets"])), f"KA {fold_id} {name}: targets differ"
            continue
        assert k in arr, f"KA {fold_id} {name}: missing array {k}"
        d = float(np.max(np.abs(np.asarray(arr[k], float) - np.asarray(v, float)))) if np.size(v) else 0.0
        worst = max(worst, d)
    assert worst <= TOL, f"KA FAILED {fold_id} {name}: worst |d| {worst}"
    return dict(ref=str(p.relative_to(ROOT)), ref_sha256=sha(p), worst=worst, n_arrays=len(z) - 1)


def main(fold_id: str):
    t0 = time.time()
    log = lambda *a: print(f"[{fold_id}]", *a, f"({time.time() - t0:.0f}s)", flush=True)
    F = E.load(fold_id, check=True)
    log("engine KA:", json.dumps(F.ka, default=str)[:200])
    abc = F.panel == "ABC"
    K.OUT = OUT                                               # ka804.save writes under OUT/bulk/<fold>/
    inc0 = E.named(F, "cap6+P2g(abc)" if abc else "pct+P2g")
    inc = E.Arm(w=inc0.w, q=inc0.q, name="inc")
    fin = np.asarray(F.fin, np.float64)
    npn = F.nonpanel
    rec = dict(fold=fold_id, n=F.n, sources=F.sources, held=F.held, ka_engine=F.ka, arms={}, ka811={})
    # ---- 707, rebuilt with ka804's own functions (zsrc_fold, donor_unit_mean) and its loop
    z = K.zsrc_fold(F)
    U7, have7 = K.donor_unit_mean(F)
    tier = np.full(F.n, -1)
    r = np.argsort(np.argsort(z[have7], kind="stable"), kind="stable")
    tier[np.flatnonzero(have7)] = np.minimum((r * 3) // int(have7.sum()), 2)
    fin707 = blend(fin, U7, K.BETA, npn, rows=(tier == 2))
    rec["swap707"] = dict(n_have=int(have7.sum()), n_top=int((tier == 2).sum()),
                          change_share=float(((fin707 - fin) ** 2).sum() / (fin ** 2).sum()))
    log("707 rebuild:", rec["swap707"])
    # ---- the baseline's consensus
    groups_all = {ln: b for ln, b in GROUPS.items() if ln != HELD_GROUP.get(F.held)}
    groups_src = {f"src:{s}": (s,) for s in F.sources}
    U_all, cov_all = consensus(F, groups_all)
    U_src, cov_src = consensus(F, groups_src)
    # U_src (equal weight per source, unit rows averaged) is 707's U up to a positive per-row scale: check the direction
    nz = have7 & (np.linalg.norm(U_src, axis=1) > 0)
    cosd = np.einsum("ij,ij->i", U_src[nz], U7[nz]) / (np.linalg.norm(U_src[nz], axis=1) * np.linalg.norm(U7[nz], axis=1))
    rec["U_src_vs_707U_min_cos"] = float(cosd.min()) if nz.any() else None
    rec["U_all"], rec["U_src"] = cov_all, cov_src
    rec["groups_all"] = {k: list(v) for k, v in groups_all.items()}
    arms = {"inc": inc, "707": E.Arm(w=inc.w, q=inc.q, fin=fin707, name="707")}
    for lam in LAMS:
        tag = f"l{int(round(lam * 100)):02d}"
        arms[f"all_{tag}"] = E.Arm(w=inc.w, q=inc.q, fin=blend(fin, U_all, lam, npn), name=f"all_{tag}")
        arms[f"src_{tag}"] = E.Arm(w=inc.w, q=inc.q, fin=blend(fin, U_src, lam, npn), name=f"src_{tag}")
    if F.held == "k562_gwps":
        V = unit(bank_rows("viperturb_k562", list(F.targets)), npn)
        arms["oracle_vip"] = E.Arm(w=inc.w, q=inc.q, fin=blend(fin, V, 1.0, npn), name="oracle_vip")
        rec["oracle_vip_n"] = int((np.linalg.norm(V, axis=1) > 0).sum())
    for nm, arm in arms.items():
        if arm.fin is not None:
            ch = arm.fin - fin
            rec.setdefault("change_energy_of_fin", {})[nm] = float((ch * ch).sum() / (fin * fin).sum())
    ref_dir = REF804 / fold_id.replace("/", "__")
    for nm, arm in arms.items():
        b = None if nm == "inc" else "inc"
        s, arr = E.price(F, arm, base=arms[b] if b else None)
        arr["targets"] = np.asarray(F.targets)
        rec["arms"][nm] = K.save(fold_id, nm, s, arr, b)
        if nm in ("inc", "707"):
            keep = {k: v for k, v in arr.items() if not k.startswith("base_")}
            try:
                rec["ka811"][nm] = compare(fold_id, nm, keep, ref_dir)
                log(f"KA {nm} == 804 (worst |d| {rec['ka811'][nm]['worst']:.1e}, {rec['ka811'][nm]['n_arrays']} arrays)")
            except AssertionError as e:
                rec["ka811"][nm] = dict(FAILED=str(e))
                log(f"KA {nm} FAILED: {e}")
        log(f"{nm:10s} cos_agg {s['cos_agg']:.6f} pds {s['pds']:.5f} K {s['K_mean']:.2f}" +
            (f" | vs inc: r {100 * s['vs_base']['r']:+.3f}% dpds {s['vs_base']['dpds']:+.5f} dpdsl {s['vs_base']['dpds_lam05']:+.5f}"
             f" dK {100 * s['vs_base']['dK_rel']:+.2f}% unscored {s['vs_base']['share_unscored']:.2f}" if b else ""))
    (OUT / "meta").mkdir(parents=True, exist_ok=True)
    rec["banks_manifest_sha256"] = sha(BANKS / "manifest.json")
    rec["script_sha256"] = sha(Path(__file__))
    (OUT / "meta" / f"{fold_id.replace('/', '__')}.json").write_text(json.dumps(rec, indent=1, default=float))
    bad = [k for k, v in rec["ka811"].items() if "FAILED" in v]
    assert not bad, f"KA FAILED on {bad}: this fold's numbers are void"
    log("done")


if __name__ == "__main__":
    main(sys.argv[1])
