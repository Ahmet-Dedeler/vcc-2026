import numpy as np, scipy.stats as ss, sys
sys.path.insert(0, "experiments/521-strength-nreal")
import lib521 as L
S = L.load_state(0)                       # 435's production-emitter stats on H1, seed 0 (oq/op/ol) + the real reference
rl, rq, ol, oq = S["rl"], S["rq"], S["ol"], S["oq"]
own = S["own"]; ref_lin = np.asarray(S["ref_lin"], float); w = ref_lin / ref_lin.sum()
C = S["C"]; key = S["key"]
T = rl.shape[0]
out = []
for i in range(T):
    f = np.isfinite(rl[i]) & ~own[i]
    null = f & ~(rq[i] < 0.05)
    de = f & (rq[i] < 0.05)
    m_real = float(np.sum(w[de] * (2.0 ** rl[i, de] - 1)))
    fo = np.isfinite(ol[i])
    m_ours_all = float(np.nansum(w[fo] * (2.0 ** ol[i, fo] - 1)))                 # our emitted net mass change (all genes)
    m_key = float(np.sum(w * (2.0 ** np.nan_to_num(key[i]) - 1)))                  # the design key's implied mass change
    real_up = float(np.mean(rl[i, null] > 0))
    ours_up = float(np.mean(ol[i, null & fo] > 0))                                 # our emitted sign on the REAL null genes
    acc_null_all = float(np.mean(np.sign(ol[i, null & fo]) == np.sign(rl[i, null & fo])))
    callnull = null & C[i] & fo
    acc_null_calls = float(np.mean(np.sign(ol[i, callnull]) == np.sign(rl[i, callnull]))) if callnull.any() else np.nan
    up_calls_real = float(np.mean(rl[i, callnull] > 0)) if callnull.any() else np.nan
    out.append((m_real, m_ours_all, m_key, real_up, ours_up, acc_null_all, acc_null_calls, up_calls_real))
A = np.array(out)
lab = ["m_real", "m_ours", "m_key", "real_up", "ours_up", "acc_null_all", "acc_null_calls", "realup_calls"]
print("means:", {l: round(float(np.nanmean(A[:, j])), 4) for j, l in enumerate(lab)})
for a, b in (("m_ours", "m_real"), ("m_key", "m_real"), ("ours_up", "real_up"), ("m_ours", "real_up"), ("m_key", "real_up")):
    x, y = A[:, lab.index(a)], A[:, lab.index(b)]
    s = ss.spearmanr(x, y)
    print(f"spearman({a}, {b}) {s.correlation:+.3f} (p {s.pvalue:.2g})")
# rule: sign of null CALLS = -sign(predicted mass change); accuracy on real null genes among our calls
for src in ("m_real", "m_ours", "m_key"):
    m = A[:, lab.index(src)]
    acc = np.where(m < 0, A[:, lab.index("realup_calls")], 1 - A[:, lab.index("realup_calls")])
    print(f"rule '-sign({src})' on our null CALLS: p_null {np.nanmean(acc):.4f}   (ours today {np.nanmean(A[:, lab.index('acc_null_calls')]):.4f})")
