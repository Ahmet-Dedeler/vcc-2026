"""Q4 — how much of what the harness says survives power matching to the challenge regime?"""
import numpy as np
from debench import Surface, overall_delta

print(f"{'n_real regime':>16} {'n_real':>8} {'n_pred':>8} {'reach':>9} {'fid':>9} {'jac':>9} {'nmae':>9} {'signacc':>8}")
for tgt in (None, 250, 500, 1000, 2000, 3000):
    s = Surface("../095-clean-harness/de_thinned.npz", n_real_target=tgt)
    nr, npd = s.counts()
    lab = "native (BH .05)" if tgt is None else f"top-{tgt}"
    print(f"{lab:>16} {np.median(nr):8.0f} {np.median(npd):8.0f} {s.reach():9.5f} {s.fid():9.5f} "
          f"{s.jac():9.5f} {s.nmae():9.5f} {s.sign_accuracy():8.4f}")

print("\ncall-set truncation, priced in each regime (keep top-M of OUR calls by pred p_adj)")
print(f"{'regime':>10} " + "".join(f"{f'M={m}':>22}" for m in (500, 2000, 4000)))
for tgt in (None, 2000):
    s = Surface("../095-clean-harness/de_thinned.npz", n_real_target=tgt)
    row = f"{('native' if tgt is None else 'top-2000'):>10} "
    b_f, b_j = s.fid(), s.jac()
    for M in (500, 2000, 4000):
        call = np.zeros(s.pred_padj.shape, bool)
        for p in range(len(s.perts)):
            ok = np.flatnonzero(s.pred_sig[p] & s.adj[p])
            if ok.size: call[p, ok[np.argsort(s.pred_padj[p, ok])[:M]]] = True
        d = overall_delta(d_fid=s.fid(call=call) - b_f, d_jac=s.jac(call=call) - b_j)
        row += f"  Δfid{s.fid(call=call)-b_f:+.4f} Δjac{s.jac(call=call)-b_j:+.4f} ->{d:+.4f}"
    print(row)
