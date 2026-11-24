"""804 — the bulk known-answer and prospective arms on every bulk fold, through 801's engine (engine801.price), written in
800's STANDARD FORMAT for check800 (PROTOCOL.md; this file adds no arithmetic of its own beyond the arm constructions,
each copied from its production source and cited).

Per fold (L1-L6 x {k562_gwps, orion_hct116, orion_hek293t}, ABC x {k562_gwps, orion_hct116, orion_hek293t, kolf2}):
  named arms (engine801.named): cap3, cap6, pct, cap3+P2g, cap3+S_P2, pct+P2g (incumbent), cap6+P2g; on ABC folds the
    board's plain tercile `+P2g(abc)` (479p2g / 498p2g / ALLSW used the builder's plain `tercile`, arms540c's p2g_abc).
  702  : the incumbent with PIE's K562 rows blended (w 0.25, PIE-seen targets, norm-preserving on the K562 axis's non-panel
         genes; arm702.py = make_pie_dump700's arithmetic) into the composition rows. K562 folds only.
  705  : the incumbent's key weights times exp(0.25 z_src) above the floor, never below it (def_strength_key522.rider,
         the 'fpc' form), P2g re-cut on the new weights (arms522c's fpc + P2g).
  707  : the top z_src tertile's composition rows blended toward the fold SOURCES' unit mean U (swap707 / dirswap706),
         beta 0.25, norm-preserving on non-panel genes.
z_src on a fold = 522's src feature over the fold's own source donors only (strength522.robust_strength on the panel,
crank over the fold's targets, mean over the sources that carry the target): the held line never enters.

Outputs: work/run/bulk/<fold with / as __>/<arm>.npz (targets, x, a, d, pds, pdsl, K [+ dE_* vs the arm's declared base]),
work/run/meta/<fold>.json (fold KA record, sha256 of inputs, arm summaries). Run on EC2 (one process per fold):
    uv run python experiments/804-known-answers/ka804.py <fold-id>
"""
from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/801-lab-folds"))
sys.path.insert(0, str(ROOT / "experiments/522-october-bulk-key"))
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402

import engine801 as E  # noqa: E402  (imports lib801 -> guard801: writes outside 801 raise ... see OUTDIR note)
import arm702  # noqa: E402
from strength522 import robust_strength, BANKS  # noqa: E402

import os  # noqa: E402
# engine801 imports guard801, which refuses writes inside the repo outside experiments/801-lab-folds; on EC2 the job writes
# under KA804_OUT (outside the repo) and the job line copies the result into this experiment's work/ afterwards.
OUT = Path(os.environ.get("KA804_OUT", str(HERE / "work" / "run")))
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
PIE_L = HERE / "work/pie/pred_k562_L.parquet"
PIE_ABC = ROOT / "experiments/700-pie-transfer/work/cloud/20261014T225548Z/out/pred_k562.parquet"
GENES = [ln.split(",")[0] for ln in (ROOT / "data/challenge/gene_names.csv").read_text().splitlines()[1:] if ln.strip()]
FLOOR, THETA, BETA, W702 = 0.5, 0.25, 0.25, 0.25


def sha(p) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def crank(v):
    v = np.asarray(v, float)
    ok = np.isfinite(v)
    out = np.zeros(v.size)
    if ok.sum() > 1:
        r = np.argsort(np.argsort(v[ok]))
        out[ok] = r / (ok.sum() - 1) - 0.5
    return out


def panel_targets(F) -> list[str]:
    if F.panel == "ABC":
        p = ROOT / "data/challenge/pert_counts.csv"
    else:
        p = ROOT / f"experiments/801-lab-folds/work/{F.panel}/bundle/pert_counts.csv"
    lines = [ln for ln in p.read_text().splitlines() if ln.strip()]
    j = [c.strip() for c in lines[0].split(",")].index("target_gene")
    return [ln.split(",")[j].strip() for ln in lines[1:]]


def zsrc_fold(F) -> np.ndarray:
    """522's src feature on the fold: per SOURCE donor, robust strength over the panel's targets (panel genes out), indexed
    to the fold's targets, centred rank over the fold's targets; mean over the sources that carry the target; 0 if none."""
    panel = panel_targets(F)
    pset = set(panel)
    idx = np.asarray([panel.index(t) for t in F.targets])
    cols = []
    for src in F.sources:
        if src not in BANKS:
            continue
        s, _, _ = robust_strength(Path(BANKS[src]), panel, pset)
        v = np.asarray(s, float)[idx]
        r = crank(v)
        r[~np.isfinite(v)] = np.nan
        cols.append(r)
    with np.errstate(all="ignore"):
        m = np.nanmean(np.vstack(cols), 0)
    return np.where(np.isfinite(m), m, 0.0)


def donor_unit_mean(F) -> tuple[np.ndarray, np.ndarray]:
    """dirswap706.donor_rows over the fold's SOURCES: U = sum_j r_j/|r_j| (raw delta rows, official axis, panel cols 0)."""
    U = np.zeros((F.n, len(GENES)))
    gpos_off = {g: i for i, g in enumerate(GENES)}
    for src in F.sources:
        D = ROOT / "data/donors" / src
        dg = (D / "genes.txt").read_text().split()
        dp = (D / "perts.txt").read_text().split()
        delta = np.load(D / "delta.npy", mmap_mode="r")
        ppos = {p: i for i, p in enumerate(dp)}
        cols = np.array([gpos_off.get(g, -1) for g in dg])
        ok = cols >= 0
        for i, t in enumerate(F.targets):
            if t not in ppos:
                continue
            r = np.zeros(len(GENES))
            r[cols[ok]] = np.asarray(delta[ppos[t]], np.float64)[ok]
            r[~F.nonpanel] = 0.0
            nr = np.linalg.norm(r)
            if nr > 0:
                U[i] += r / nr
    return U, np.linalg.norm(U, axis=1) > 0


def p2g_for(F, w):
    return E.p2g_abc(w) if F.panel == "ABC" else E.p2g_oct(w, F.carried)


def save(fold_id, name, s, arr, base_name=None):
    d = OUT / "bulk" / fold_id.replace("/", "__")
    d.mkdir(parents=True, exist_ok=True)
    keep = {k: v for k, v in arr.items() if not k.startswith("base_")}
    keep["targets"] = np.asarray(keep.get("targets", []))
    np.savez_compressed(d / f"{name}.npz", **keep)
    return dict(name=name, base=base_name, cos_agg=s["cos_agg"], pds=s["pds"], K_mean=s["K_mean"],
                vs_base=s.get("vs_base"))


def main(fold_id: str):
    t0 = time.time()
    F = E.load(fold_id, check=True)
    log = lambda *a: print(f"[{fold_id}]", *a, f"({time.time() - t0:.0f}s)", flush=True)
    log("engine KA:", json.dumps(F.ka, default=str)[:300])
    abc = F.panel == "ABC"
    arms = {}
    for nm in ("cap3", "cap6", "pct", "cap3+S_P2"):
        arms[nm] = E.named(F, nm)
    arms["cap3+P2g"] = E.named(F, "cap3+P2g(abc)" if abc else "cap3+P2g")
    arms["cap6+P2g"] = E.named(F, "cap6+P2g(abc)" if abc else "cap6+P2g")
    arms["pct+P2g"] = E.named(F, "pct+P2g")
    inc = arms["cap6+P2g"] if abc else arms["pct+P2g"]       # the board's ALLSW bulk on A/B/C; October's incumbent elsewhere
    arms["inc"] = E.Arm(w=inc.w, q=inc.q, name="inc")
    rec = dict(fold=fold_id, n=F.n, sources=F.sources, ka=F.ka, inc=("cap6+P2g(abc)" if abc else "pct+P2g"), arms={})
    # 705: the strength rider on the incumbent's weights (fpc), P2g re-cut on the new weights
    z = zsrc_fold(F)
    w0 = np.asarray(inc.w, float)
    at_floor = w0 <= FLOOR + 1e-12
    w705 = np.where(at_floor, w0, np.maximum(FLOOR, w0 * np.exp(THETA * z)))
    arms["705"] = E.Arm(w=w705, q=p2g_for(F, w705), name="705")
    rec["z705"] = dict(n_floor=int(at_floor.sum()), mult_q=np.quantile(np.exp(THETA * z[~at_floor]), [0, .5, 1]).tolist()
                       if (~at_floor).any() else None, n_changed=int(np.sum(w705 != w0)))
    # 707: the strong-tier swap on the composition rows
    U, have = donor_unit_mean(F)
    tier = np.full(F.n, -1)
    r = np.argsort(np.argsort(z[have], kind="stable"), kind="stable")
    tier[np.flatnonzero(have)] = np.minimum((r * 3) // int(have.sum()), 2)
    sel = tier == 2
    fin = np.asarray(F.fin, np.float64)
    fin707 = fin.copy()
    npn = F.nonpanel
    for i in np.flatnonzero(sel):
        nb, nu = np.linalg.norm(fin[i, npn]), np.linalg.norm(U[i, npn])
        if nb > 0 and nu > 0:
            m = (1 - BETA) * fin[i, npn] / nb + BETA * U[i, npn] / nu
            fin707[i, npn] = nb * m / np.linalg.norm(m)
    arms["707"] = E.Arm(w=inc.w, q=inc.q, fin=fin707, name="707")
    rec["swap707"] = dict(n_have=int(have.sum()), n_top=int(sel.sum()),
                          change_share=float(((fin707 - fin) ** 2).sum() / (fin ** 2).sum()))
    # 702: PIE blend on the K562 axis (K562 folds only)
    if F.held == "k562_gwps":
        pred = PIE_ABC if abc else PIE_L
        pie = arm702.pie_rows(pred, list(F.targets), GENES)
        kg, _, kpos = arm702.k562_axis(GENES)
        cols = np.zeros(len(GENES), bool)
        cols[kpos] = True
        cols &= npn
        seen = arm702.seen_mask(list(F.targets)) & (np.abs(pie).sum(1) > 0)
        fin702 = arm702.blend_axis(fin, pie, W702, cols, seen)
        arms["702"] = E.Arm(w=inc.w, q=inc.q, fin=fin702, name="702")
        rec["pie702"] = dict(pred=str(pred.relative_to(ROOT)), pred_sha256=sha(pred), n_seen=int(seen.sum()),
                             n_rows=int((np.abs(pie).sum(1) > 0).sum()), change_share_k562axis=arm702.change_share(fin, fin702, cols),
                             change_energy_of_fin=float(((fin702 - fin) ** 2).sum() / (fin ** 2).sum()))
    # price: each arm against the base it is read against (the change energy is recorded vs that base)
    bases = {"cap3": None, "cap6": "cap3", "pct": "cap3", "cap3+S_P2": "cap3", "cap3+P2g": "cap3", "cap6+P2g": "cap6",
             "pct+P2g": "pct", "inc": None, "705": "inc", "707": "inc", "702": "inc"}
    for nm, arm in arms.items():
        b = bases[nm]
        s, arr = E.price(F, arm, base=arms[b] if b else None)
        arr["targets"] = np.asarray(F.targets)
        rec["arms"][nm] = save(fold_id, nm, s, arr, b)
        log(f"{nm:10s} cos_agg {s['cos_agg']:.6f} pds {s['pds']:.5f} K {s['K_mean']:.2f}" +
            (f" | vs {b}: r {100 * s['vs_base']['r']:+.3f}% dpds {s['vs_base']['dpds']:+.5f}" if b else ""))
    (OUT / "meta").mkdir(parents=True, exist_ok=True)
    (OUT / "meta" / f"{fold_id.replace('/', '__')}.json").write_text(json.dumps(rec, indent=1, default=float))
    log("done")


if __name__ == "__main__":
    main(sys.argv[1])
