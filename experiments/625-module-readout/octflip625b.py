"""625 A4 evidence (10-03): per CLEAN4 octile, how often FLIPPING it would help a held-out line, i.e. CLEAN4's sign agreement
with the line's truth on that octile's gate genes < 0.5. Broad-16's 12 working lines (622's loader; KO modality, 627's
caveat) and H1 (h1-test, lib521; 625d found CLEAN4 agrees on every octile). Octiles from modules625_v8.npz.
Agreement is pooled over targets, by count and |lfc|-weighted. Analysis only, < 4 GB.
usage: uv run python experiments/625-module-readout/octflip625b.py"""
import json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments/622-identity-track"))
sys.path.insert(0, str(ROOT / "experiments/521-strength-nreal"))
import broad622 as B  # noqa: E402
import lib521 as L  # noqa: E402
V8 = np.load(ROOT / "experiments/625-module-readout/work/modules625_v8.npz", allow_pickle=True)
OCT = dict(zip(V8["genes"].astype(str), V8["octile"].astype(int)))
c4 = np.load(ROOT / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
names = list(c4["bank_names"].astype(str))
Pc = c4["P_banks"].astype(float)[[names.index(b) for b in ("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2")]].mean(0)
C4 = dict(zip(c4["genes"].astype(str), Pc))
out, res = [], {}
say = lambda s="": (print(s), out.append(s))  # noqa: E731


def agree(genes, lfc, gate):
    """per octile: (count agreement, |lfc|-weighted agreement, n cells) of CLEAN4's sign with the truth's sign on gate cells"""
    oc = np.array([OCT.get(str(g), -1) for g in genes]); p = np.array([C4.get(str(g), np.nan) for g in genes])
    c4dn = p > 0.5; cov = np.isfinite(p) & (p != 0.5)
    r = {}
    for k in range(8):
        m = gate & ((oc == k) & cov)[:, None] if gate.shape[0] == len(genes) else None
        tdn = lfc < 0
        ok = (c4dn[:, None] == tdn) & m
        w = np.abs(np.nan_to_num(lfc))
        n = int(m.sum())
        r[k] = (float(ok.sum() / max(n, 1)), float((w * ok).sum() / max((w * m).sum(), 1e-12)), n)
    return r


gl, LFC, cl = B.load_matrix(B.P / "lfc_matrix.csv", 1)
gf, FDR, cf = B.load_matrix(B.P / "fdr_matrix.csv", 1)
lines = sorted({ln for ln, _ in cl}); assert len(lines) == 12
colline = np.array([ln for ln, _ in cl]); coltgt = np.array([k for _, k in cl])
per = {k: [] for k in range(8)}
for ln in lines:
    jl = np.flatnonzero(colline == ln)
    gate = (FDR[:, jl] < 0.05) & np.isfinite(LFC[:, jl]) & (LFC[:, jl] != 0) & (np.asarray(gl)[:, None] != coltgt[jl][None, :])
    r = agree(gl, LFC[:, jl], gate)
    for k in range(8):
        per[k].append((ln,) + r[k])
S = L.load_state(0)
g1 = (S["rq"] < 0.05) & np.isfinite(S["rl"]) & (S["rl"] != 0) & ~S["own"]
rh = agree(S["genes"], S["rl"].T, g1.T)
say("CLEAN4 octile: lines (of 12) where flipping would help [count agreement < 0.5] | mean agreement (count / |lfc|-w) | H1 agreement (count / w)")
for k in range(8):
    a = np.array([x[1] for x in per[k]]); aw = np.array([x[2] for x in per[k]])
    helps = int((a < 0.5).sum()); helps_w = int((aw < 0.5).sum())
    say(f"  octile {k}: {helps}/12 (|lfc|-w {helps_w}/12) | {a.mean():.3f} / {aw.mean():.3f} (min {a.min():.3f}) | H1 {rh[k][0]:.3f} / {rh[k][1]:.3f} (n {rh[k][2]})")
    res[k] = dict(lines_flip_helps=helps, lines_flip_helps_w=helps_w, mean_agree=float(a.mean()), mean_agree_w=float(aw.mean()),
                  per_line={x[0]: [x[1], x[2], x[3]] for x in per[k]}, h1=list(rh[k]))
(Path(__file__).resolve().parent / "work/octflip625b.json").write_text(json.dumps(res, indent=1))
(Path(__file__).resolve().parent / "work/octflip625b.out").write_text("\n".join(out) + "\n")
