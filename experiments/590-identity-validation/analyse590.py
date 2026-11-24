"""590 — the preregistered read: FP_reg / TP_reg / accuracy / sister rows / abstention, plain and Flex-simulated,
and the proposed gate's cost. Writes work/summary590.json and prints the tables RESULTS.md quotes."""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import def_sameline as DS  # noqa: E402
import def_identity as DI  # noqa: E402

REG = set(DS.REGISTRY)
GATE = 0.65  # proposed: min over the 8 methods of the top-1 correlation
# profiles whose gene axis could not be read (documented in RESULTS): never scored
UNREADABLE = {"sciplex3", "cui", "liang", "gehring"}


def load(tag: str) -> pd.DataFrame:
    d = pd.read_csv(HERE / f"work/match590{tag}.csv")
    return d[~d.dataset.isin(UNREADABLE) & (d["query"] != "donor:kolf2")].copy()


def flex_combined() -> pd.DataFrame:
    out = []
    for k in (0, 1):
        d = load(f"_flexbias_fold{k}")
        z = np.load(HERE / f"work/flexbias590_fold{k}.npz")
        applies = {q.split(":", 1)[1] for q in z["apply_to"].astype(str)}
        mix = d.dataset.str.startswith("mcfarland")
        lid = d["query"].str.split(":", n=1).str[1]
        keep = (mix & lid.isin(applies)) | (~mix if k == 0 else False)
        out.append(d[keep])
    return pd.concat(out, ignore_index=True)


def read(d: pd.DataFrame, gate: float | None) -> dict:
    sw = d.switch_on & (d.min_s1 >= gate) if gate is not None else d.switch_on
    t = d[d.true_in_ccle]
    _, _, _, clinical = DI.load_cache()
    lin = clinical["ONCOTREE_LINEAGE"].to_dict()
    reg_lineages = {lin.get(r) for r in REG}
    fp = d[sw & (d.call != d.true)]
    regrows = d[d.true.isin(REG)]
    sister = d[d.true_in_ccle & ~d.true.isin(REG) & d.true_lineage.isin(reg_lineages)]
    return {
        "n": int(len(d)), "n_in_ccle": int(len(t)),
        "top1_acc": round(float(t.correct.mean()), 4), "p_unanimous": round(float(d.unanimous.mean()), 4),
        "acc_given_unanimous": round(float(t[t.unanimous].correct.mean()), 4),
        "switch_on": int(sw.sum()), "false_switch_on": int(len(fp)),
        "false_switch_on_rows": [f"{q} -> {c} (min_s1 {s:.3f})" for q, c, s in zip(fp["query"], fp.call, fp.min_s1)],
        "FP_reg_rate": round(len(fp) / max(1, int((~d.true.isin(REG)).sum())), 4),
        "registry_rows": int(len(regrows)), "registry_rows_switched_on_correctly": int((sw & d.correct)[regrows.index].sum()),
        "sister_rows": int(len(sister)), "sister_false_switch_on": int((sw & (d.call != d.true))[sister.index].sum()),
        "challenge_rows": {q: {"call": c, "unanimous": bool(u), "min_s1": round(float(s), 3), "min_gap12": round(float(g), 4)}
                           for q, c, u, s, g in zip(d["query"], d.call, d.unanimous, d.min_s1, d.min_gap12)
                           if q in ("A", "B", "C")},
    }


def main() -> None:
    plain = load("")
    plain = plain[plain["query"].isin(flex_combined()["query"]) | plain.dataset.isin(["challenge", "anchor"])]
    flex = flex_combined()
    s = {"gate_min_s1": GATE}
    for name, d in (("plain", plain), ("flex_simulated", flex)):
        s[name] = {"current_rule": read(d, None), "with_gate": read(d, GATE)}
        cu = d[d.true_in_ccle & d.correct & d.unanimous]
        s[name]["correct_unanimous_min_s1_quantiles"] = {q: round(float(cu.min_s1.quantile(q)), 4) for q in (0.01, 0.05, 0.5)}
        s[name]["false_call_min_s1"] = sorted(round(float(x), 3) for x in d[d.switch_on & (d.call != d.true)].min_s1)
        s[name]["gate_rejects_correct"] = [f"{q} ({c}, {s1:.3f})" for q, c, s1 in
                                           zip(*[d[d.switch_on & d.correct & (d.min_s1 < GATE)][k] for k in ("query", "chemistry", "min_s1")])]
        s[name]["not_in_ccle"] = [f"{q} -> {c} unanimous={bool(u)} min_s1={m:.3f}" for q, c, u, m in
                                  zip(*[d[~d.true_in_ccle][k] for k in ("query", "call", "unanimous", "min_s1")])]
        s[name]["wrong_in_ccle"] = [f"{q}: true {t} -> {c} unanimous={bool(u)}" for q, t, c, u in
                                    zip(*[d[d.true_in_ccle & ~d.correct][k] for k in ("query", "true", "call", "unanimous")])]
    (HERE / "work/summary590.json").write_text(json.dumps(s, indent=1) + "\n")
    print(json.dumps(s, indent=1))


if __name__ == "__main__":
    main()
