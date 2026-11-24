"""269 -- the leak, stated formally and measured: for held X and action 'choose Z', the training label on fold Y (donors X, Z)
and the held label on fold X (donors Y, Z) share the term cos(x, y) -- the held line's own response against the OTHER
donor's -- which is a function of the held truth.  Partial Spearman of (training score, held label) given the held fold's
truth-side agreement with donor Y (cov263 `cos_truth_wht`, an oracle quantity used for analysis only)."""
import json
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr
ROOT = Path(__file__).resolve().parents[2]; WORK = Path(__file__).resolve().parent / "work"
E263 = ROOT / "experiments/263-donor-selection-arbiters/work"
FOLD = {"K562": "k562_gwps", "HCT116": "orion_hct116", "HEK293T": "orion_hek293t", "KOLF2": "kolf2"}
DID = {"k562_gwps_pm50": "K562", "orion_hct116": "HCT116", "orion_hek293t": "HEK293T"}
CANCER = ["K562", "HCT116", "HEK293T"]
dfc = pd.read_csv(WORK / "rows_choose269.csv")
out = {}
def rank(v): return pd.Series(v).rank().to_numpy()
def partial_spearman(x, y, z):
    X = np.column_stack([rank(x), rank(y), rank(z)]).astype(float)
    X -= X.mean(0)
    # residualise x and y on z
    bz = X[:, 2]
    rx = X[:, 0] - bz * (X[:, 0] @ bz) / (bz @ bz)
    ry = X[:, 1] - bz * (X[:, 1] @ bz) / (bz @ bz)
    return float(np.corrcoef(rx, ry)[0, 1])
for X in CANCER:
    cov = np.load(E263 / f"cov263-{FOLD[X]}.npz", allow_pickle=True)
    meta = json.load(open(E263 / f"sweep263-{FOLD[X]}.json"))["meta"]
    donors = [DID[d] for d in meta["donors"]]
    targets = cov["targets"].tolist()
    ctw = np.asarray(cov["cos_truth_wht"], float)              # held truth vs each donor's whitened row
    for Z in donors:
        Y = [d for d in donors if d != Z][0]
        held = dfc[(dfc.context == X) & (dfc.donor == Z)].set_index("target").u
        train = dfc[(dfc.context == Y) & (dfc.donor == Z)].set_index("target").u     # fold Y: donors X, Z -> 'choose Z' = drop X
        common = [t for t in targets if t in held.index and t in train.index]
        h = held.loc[common].to_numpy(); tr = train.loc[common].to_numpy()
        cy = ctw[[targets.index(t) for t in common], donors.index(Y)]    # cos(x_truth, y_row) on the held fold: the shared term
        cz = ctw[[targets.index(t) for t in common], donors.index(Z)]
        r = spearmanr(tr, h).correlation
        r_y = spearmanr(tr, cy).correlation
        pr = partial_spearman(tr, h, cy)
        pr2 = partial_spearman(tr, h, cy - cz)
        out[f"held {X}, choose {Z} (training fold {Y})"] = dict(n=len(common), spearman=float(r), spearman_train_vs_cos_x_y=float(r_y), partial_given_cos_xy=float(pr), partial_given_cos_xy_minus_xz=float(pr2))
        print(f"held {X:8s} choose {Z:8s} (training fold {Y:8s}): n {len(common)}  Spearman(train, held) {r:+.3f}; Spearman(train, cos(x_truth, y_row)) {r_y:+.3f}; "
              f"partial given cos(x,y) {pr:+.3f}; partial given cos(x,y)-cos(x,z) {pr2:+.3f}")
(WORK / "mech269.json").write_text(json.dumps(out, indent=1))
