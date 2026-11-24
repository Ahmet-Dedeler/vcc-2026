#!/usr/bin/env python3
"""636 — the council's leave-one-LINE-out instrument for DE sign / membership sources (portfolio #1; owner d4cd5f).

Scores any per-(target, gene) prediction of a held-out line's DE signs (P_down) and, optionally, DE membership (a score
whose high values mean "this gene is in the target's reference gate") on four surfaces whose gates were built with the
scorer's own rules (Wilcoxon vs the surface's own controls, BH per target over genes > 5 CPM, lfc = log2 ratio of
arithmetic mean CPM; 480's gates480.py, known answers in 480 §1):
  h1        h1-test, Arc Flex, 85 targets (323's reference; full power: 400 cells, ~20k UMI)
  hct116    Orion HCT116 thinned (132 targets, 150 cells)     -- power-limited, 3' FiCS
  hek293t   Orion HEK293T thinned (132 targets, 150 cells)    -- power-limited, 3' FiCS
  k562      K562-GWPS harness (267 targets)                   -- 3'
Why gate genes: 635 showed the shared-control-noise objection does not touch gate-gene signs, and that the reachable
"programme" prize lives there (clean oracle 0.618 vs CLEAN4 0.573 on h1).

Leave-one-line-out rule (enforced by `allowed_banks`): a source may not use the surface's own line, nor a same-study
sibling (430: same-study siblings inflate fold reads; Orion HCT116 <-> HEK293T), nor any H1 bank on h1 (oracle leakage).

API:  from inst636 import Inst;  I = Inst();  I.score(name, P_down=<dict surface -> (T,G) or (G,) array on that
      surface's axes>, member=<optional dict>)  -> rows;  I.axes(surface) gives (targets, genes).
CLI:  uv run python inst636.py            -> reference rows (all-down, CLEAN4 LOLO, per-bank and pooled unconditional
                                             programmes) + known answers. Writes work/ref636.json.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
WORK.mkdir(exist_ok=True)
G480 = ROOT / "experiments/480-nmae-instrument/work"
SURF = {"h1": "gate-h1_full_ka.npz", "hct116": "gate-orion_hct116.npz", "hek293t": "gate-orion_hek293t.npz",
        "k562": "gate-k562_gwps.npz"}
# banks on disk that a source may draw on, and which surfaces each is FORBIDDEN on (own line or same-study sibling)
BANK_FORBIDDEN = {
    "k562_gwps": {"k562"}, "k562_essential": {"k562"}, "viperturb_k562": {"k562"},
    "orion_hct116": {"hct116", "hek293t"}, "orion_hek293t": {"hct116", "hek293t"},
    "rpe1": set(), "nadig_jurkat": set(), "nadig_hepg2": set(),
    "h1_2025": {"h1"}, "h1_2025_full": {"h1"},
}
ESS = {ln.strip() for ln in open(ROOT / "data/donors/k562_essential/perts.txt") if ln.strip()}
N_BOOT = 2000


def allowed_banks(surface: str) -> list[str]:
    return [b for b, f in BANK_FORBIDDEN.items() if surface not in f]


class Inst:
    def __init__(self):
        self.S = {}
        for s, f in SURF.items():
            z = np.load(G480 / f, allow_pickle=True)
            T, G = z["targets"].astype(str), z["genes"].astype(str)
            q, lfc = z["q"].astype(float), z["lfc"].astype(float)
            gate = (q < 0.05) & np.isfinite(lfc) & (lfc != 0) & (G[None, :] != T[:, None])
            tested = np.isfinite(lfc) & (G[None, :] != T[:, None])
            self.S[s] = dict(T=T, G=G, lfc=lfc, gate=gate, tested=tested, noness=np.array([t not in ESS for t in T]),
                             gpos={g: i for i, g in enumerate(G)})

    def axes(self, s):
        return self.S[s]["T"], self.S[s]["G"]

    def gene_vec(self, s, genes, values):
        """Map a per-gene vector given on `genes` onto surface s's gene axis (NaN where absent)."""
        out = np.full(len(self.S[s]["G"]), np.nan)
        pos = self.S[s]["gpos"]
        for g, v in zip(genes, values):
            j = pos.get(str(g))
            if j is not None:
                out[j] = v
        return out

    def _sign_rows(self, s, P, label, ref=None):
        D = self.S[s]
        P2 = np.broadcast_to(P, D["lfc"].shape) if P.ndim == 1 else P
        res = []
        for scope, tm in (("all", np.ones(len(D["T"]), bool)), ("noness", D["noness"])):
            m = D["gate"] & tm[:, None] & np.isfinite(P2) & (P2 != 0.5)
            truth_dn = D["lfc"] < 0
            ok = ((P2 > 0.5) == truth_dn) & m
            w = np.abs(np.nan_to_num(D["lfc"])) * m
            ngate = int((D["gate"] & tm[:, None]).sum())
            up, dn = m & ~truth_dn, m & truth_dn
            r = dict(surface=s, source=label, targets=scope, n=int(m.sum()), coverage=float(m.sum() / max(ngate, 1)),
                     acc=float(ok.sum() / max(m.sum(), 1)), acc_w=float((w * ok).sum() / max(w.sum(), 1e-12)),
                     acc_up=float((ok & up).sum() / max(up.sum(), 1)), acc_dn=float((ok & dn).sum() / max(dn.sum(), 1)))
            if ref is not None:  # paired vs a reference source on the genes BOTH cover; target bootstrap of the difference
                R2 = np.broadcast_to(ref, D["lfc"].shape) if ref.ndim == 1 else ref
                mm = m & np.isfinite(R2) & (R2 != 0.5)
                a = (((P2 > 0.5) == truth_dn) & mm).sum(1).astype(float)
                b = (((R2 > 0.5) == truth_dn) & mm).sum(1).astype(float)
                n = mm.sum(1).astype(float)
                keep = n > 0
                a, b, n = a[keep], b[keep], n[keep]
                rng = np.random.default_rng(636)
                idx = rng.integers(0, len(n), size=(N_BOOT, len(n)))
                d = (a[idx].sum(1) - b[idx].sum(1)) / n[idx].sum(1)
                r.update(d_vs_ref=float((a.sum() - b.sum()) / max(n.sum(), 1)), d_ci=[float(np.quantile(d, .025)),
                         float(np.quantile(d, .975))], n_shared=int(n.sum()))
            res.append(r)
        return res

    def _member_rows(self, s, M, label):
        D = self.S[s]
        M2 = np.broadcast_to(M, D["lfc"].shape) if M.ndim == 1 else M
        res = []
        aucs, precs = [], []
        for i in range(len(D["T"])):
            t = D["tested"][i] & np.isfinite(M2[i])
            y = D["gate"][i][t]
            x = M2[i][t]
            if y.sum() == 0 or y.sum() == y.size:
                continue
            r = np.argsort(np.argsort(x))  # ranks (ties broken arbitrarily)
            aucs.append((r[y].sum() - y.sum() * (y.sum() - 1) / 2) / (y.sum() * (~y).sum()))
            k = int(y.sum())
            top = np.argsort(-x)[:k]
            precs.append(float(y[top].mean()))
        res.append(dict(surface=s, source=label, member_auroc=float(np.mean(aucs)), member_prec_at_nreal=float(np.mean(precs)),
                        n_targets=len(aucs)))
        return res

    def score(self, label, P_down=None, member=None, ref=None):
        rows = []
        for s in self.S:
            if P_down is not None and s in P_down:
                rows += self._sign_rows(s, np.asarray(P_down[s], float), label,
                                        None if ref is None or s not in ref else np.asarray(ref[s], float))
            if member is not None and s in member:
                rows += self._member_rows(s, np.asarray(member[s], float), label)
        return rows


def bank_uncond(bank: str):
    """Unconditional P_down per gene over ALL of a bank's perturbations (371's definition); genes."""
    d = ROOT / "data/donors" / bank
    if not (d / "delta.npy").exists():
        return None, None
    D = np.load(d / "delta.npy", mmap_mode="r")
    genes = [ln.strip() for ln in open(d / "genes.txt")]
    neg = np.zeros(D.shape[1], np.int64)
    for a in range(0, D.shape[0], 2000):
        neg += (np.asarray(D[a:a + 2000]) < 0).sum(0)
    return genes, neg / D.shape[0]


def clean4_lolo(I: Inst, s: str):
    """CLEAN4 (371's 4 banks) with any bank forbidden on surface s dropped."""
    cp = np.load(ROOT / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
    names = list(cp["bank_names"].astype(str))
    use = [b for b in ("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2") if s not in BANK_FORBIDDEN[b]]
    P = cp["P_banks"].astype(float)[[names.index(b) for b in use]].mean(0)
    return I.gene_vec(s, cp["genes"].astype(str), P), use


def main():
    I = Inst()
    out = {"surfaces": {s: dict(targets=len(I.S[s]["T"]), genes=len(I.S[s]["G"]), gate_pairs=int(I.S[s]["gate"].sum()),
                                noness=int(I.S[s]["noness"].sum()), allowed_banks=allowed_banks(s)) for s in I.S}}
    for s, v in out["surfaces"].items():
        print(f"surface {s}: {v}", flush=True)
    C4 = {}
    for s in I.S:
        C4[s], used = clean4_lolo(I, s)
        out["surfaces"][s]["clean4_banks"] = used
    rows = []
    rows += I.score("all-down", {s: np.full(len(I.S[s]["G"]), 0.9) for s in I.S})
    rows += I.score("CLEAN4 (LOLO: forbidden banks dropped)", C4)
    # KA: h1 CLEAN4 gate accuracy must equal 634/635's 0.5734 (same banks, same gate)
    ka = [r for r in rows if r["surface"] == "h1" and r["source"].startswith("CLEAN4") and r["targets"] == "all"][0]
    out["KA_h1_clean4_acc"] = ka["acc"]
    assert abs(ka["acc"] - 0.5734) < 6e-4, ka
    print(f"KA PASS: h1 CLEAN4 gate accuracy {ka['acc']:.4f} (634/635: 0.5734)", flush=True)
    # per-bank unconditional programmes, and the pooled mean of the allowed banks, each scored LOLO
    U = {}
    for b in BANK_FORBIDDEN:
        g, p = bank_uncond(b)
        if g is None:
            continue
        U[b] = (g, p)
        P = {s: I.gene_vec(s, g, p) for s in I.S if s not in BANK_FORBIDDEN[b]}
        rows += I.score(f"uncond {b}", P, ref={s: C4[s] for s in P})
    Pm = {}
    for s in I.S:
        stack = [I.gene_vec(s, *U[b]) for b in U if s not in BANK_FORBIDDEN[b]]
        Pm[s] = np.nanmean(np.stack(stack), 0)
    rows += I.score("uncond mean of allowed banks", Pm, ref=C4)
    out["rows"] = rows
    print(f"\n{'surface':8s} {'source':42s} {'tgt':6s} {'cov':>5s} {'acc':>6s} {'acc_w':>6s} {'up':>6s} {'dn':>6s} {'d_vs_C4 [95% CI]':>26s}")
    for r in rows:
        d = (f"{r['d_vs_ref']:+.4f} [{r['d_ci'][0]:+.4f},{r['d_ci'][1]:+.4f}]" if "d_vs_ref" in r else "")
        print(f"{r['surface']:8s} {r['source'][:42]:42s} {r['targets']:6s} {r['coverage']:5.2f} {r['acc']:6.4f} "
              f"{r['acc_w']:6.4f} {r['acc_up']:6.4f} {r['acc_dn']:6.4f} {d:>26s}")
    (WORK / "ref636.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    sys.exit(main())
