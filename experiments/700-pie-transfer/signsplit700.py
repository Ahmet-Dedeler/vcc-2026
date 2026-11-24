#!/usr/bin/env python3
"""700/721: sign accuracy on OUR call set split into (a) called & truth-significant and (b) called & not significant (the
a0 that carries ~95 % of fid's mass, 721), per sign source, all / seen / unseen targets, on h1-test (production emitter
435 seed 0 = lib521 state). Sources: our emitted lfc (ol), CLEAN4 prior sign, our design key sign, PIE lfc sign,
PIE (p_de>=0.5 else CLEAN4). Also a0 restricted to the CPM block only (the genes a block-prior swap touches).
usage: uv run python experiments/700-pie-transfer/signsplit700.py <pred_h1test.parquet>"""
import json, sys
from pathlib import Path
import numpy as np, pyarrow.parquet as pq
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/521-strength-nreal")); import lib521 as L
S = L.load_state(0); targets, genes = list(S["targets"]), list(S["genes"]); T, G = len(targets), len(genes)
t = pq.read_table(sys.argv[1]); pg = json.loads(t.schema.metadata[b"pie"])["genes"]["h1test"]; df = t.to_pandas()
pos = {g: i for i, g in enumerate(pg)}; have = np.array([g in pos for g in genes]); ix = np.array([pos.get(g, 0) for g in genes])
plfc = np.zeros((T, G)); pde = np.zeros((T, G))
for _, r in df.iterrows():
    i = targets.index(r["perturbation"]); plfc[i] = np.where(have, np.asarray(r["lfc_pred"])[ix], 0); pde[i] = np.where(have, np.asarray(r["p_de"])[ix], 0)
c4 = np.load(ROOT / "experiments/391-consensus-block-signs/work/clean4_prior_challenge.npz", allow_pickle=False)
cp = {g: i for i, g in enumerate(map(str, c4["genes"]))}
pdn = np.array([c4["P_down"][cp[g]] if g in cp else np.nan for g in genes]); c4s = np.where(np.isfinite(pdn), np.where(pdn > 0.5, -1.0, 1.0), 0.0)
own = S["own"]; rl = S["rl"]; adj = np.isfinite(rl) & (rl != 0) & ~own
sig = (S["rq"] < L.ALPHA) & np.isfinite(S["rq"]) & ~own
seen = set(json.load(open(HERE / "work/pie_train_perts.json"))); sm = np.array([x in seen for x in targets])
src = {"our_emitted": np.sign(np.nan_to_num(S["ol"])), "clean4": np.broadcast_to(c4s, (T, G)), "our_key": np.sign(S["key"]),
       "pie_lfc": np.sign(plfc), "pie_pde05_else_c4": np.where(pde >= 0.5, np.sign(plfc), c4s[None, :])}
out = {}
for setname, C in (("callset_head|block", S["C"]), ("block_only", S["block"]), ("head_only", S["head"])):
    for lab, tm in (("all", np.ones(T, bool)), ("seen", sm), ("unseen", ~sm)):
        for part, pm in (("a_sig", sig), ("b_nonsig", ~sig)):
            M = C & adj & pm & tm[:, None]
            row = {}
            for k, s in src.items():
                mm = M & (s != 0)
                row[k] = round(float((s[mm] == np.sign(rl[mm])).mean()), 4) if mm.sum() else None
            row["n"] = int(M.sum())
            out[f"{setname}|{lab}|{part}"] = row
for k, v in out.items():
    print(f"{k:34s} n={v.pop('n'):7d} " + " ".join(f"{a}={b}" for a, b in v.items()))
json.dump(out, open(HERE / "work/signsplit700.json", "w"), indent=1)
