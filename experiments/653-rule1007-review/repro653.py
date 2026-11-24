"""653: reproduce the review's findings against 537's rule1007_draft537.py (v3) itself, read-only (imports it; fake entries
built the way its own selftest builds them). Usage: uv run python experiments/653-rule1007-review/repro653.py"""
import importlib.util, json, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("r", ROOT / "experiments/537-arm-pipeline/rule1007_draft537.py")
r = importlib.util.module_from_spec(spec); spec.loader.exec_module(r)
base = {m: 0.5 for m in r.RAW.values()} | {f"score_{m}": 0.0 for m in r.RAW} | {"score_avg": 0.27}
def e(**kw):
    x = dict(base); x.update(kw); return x
fake = {r.Z: base, r.BJ2: base,
        "cpRES": e(score_fid=-0.012, score_jac=0.009, de_wilcoxon_sig_jaccard=0.503, de_wilcoxon_direction_fidelity_yield_raw=0.497),
        "hFA": e(de_wilcoxon_lfc_nmae=0.4999), "hTR": e(de_wilcoxon_lfc_nmae=0.498),
        "kPOS": e(score_fid=0.01), "kMID": e(score_fid=0.002), "cOK": e(score_avg=0.273, expr_mse_unbiased_capped_norm=0.49),
        "cphNO": e(score_fid=0.001, score_jac=0.001),
        "cphSHIPreach": e(score_fid=0.006, score_jac=0.002, score_reach=-0.04, score_avg=0.2695),
        "kLEAK": e(score_fid=0.01, pds_cosine=0.49)}
print("== R1: the expected branch (CP NULL+rescue, 496a FAILS, k16 POS). 10-07 picks CPH | K20.")
B = r.bins("cpRES", "hFA", "kPOS", "cOK", fake)
s1, s2, d8 = r.rule(B, False)
print("   10-07 lists:", s1, "|", s2)
print("   10-08 list computed on 10-07 (CPH unread, k20 arg none):", d8)
B = r.bins("cpRES", "hFA", "kPOS", "cOK", fake, cph="cphNO")
print("   10-08 list on 10-08, k20 arg = none (the docstring's 'scored on 10-06' semantics):", r.rule(B, False)[2])
print("   10-08 list on 10-08, k20 arg = K20's id (it WAS scored, on 10-07):", r.rule(B, True)[2],
      "<- 459dr, never scored, is dropped: d8 excludes the recomputed 10-07 HEADS (s2[0] became 459dr)")
print("   with 496a TRANSFERS:", r.rule(r.bins("cpRES", "hTR", "kPOS", "cOK", fake, cph="cphNO"), True)[2])
print("\n== R2: a missing capture crashes the rule instead of reading UNKNOWN")
try:
    r.bins("K1Q0y2OFYjTNRqh9jKfk", "NO_SUCH_ENTRY", "none", "none")
    print("   no crash")
except FileNotFoundError as ex:
    print("   FileNotFoundError:", Path(ex.filename).name)
print("\n== R3: no falsifier: a k16 that moved pds (bulk leak) still bins POS")
print("   ", r.bins("cpRES", "hFA", "kLEAK", "cOK", fake)["K"])
print("\n== R4: CPH 'SHIP' reads fid + jac only: an arm with fid+jac +0.0013 overall but reach -0.04 scaled (overall -0.0005) SHIPS")
B = r.bins("cpRES", "hFA", "kPOS", "cOK", fake, cph="cphSHIPreach")
print("   CPH bin:", B["CPH"], "; d8:", r.rule(B, False)[2])
print("\n== R5: the selftest's last 10-08 assertion accepts both ALLSW and 459dr, so it cannot catch R1")
