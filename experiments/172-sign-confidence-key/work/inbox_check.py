import numpy as np, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from harness172 import H1Test, check
T=H1Test(); check(T,verbose=False); print("known-answer rows OK")
pool=T.pool; real=T.S.real_lfc
hcorr = (np.sign(T.H)==np.sign(real))
print(f"H-oracle reach (key = design-sign correct)            {T.reach(hcorr.astype(float)):.5f}   INBOX says 0.62119")
print(f"full oracle  reach (key = emitted-sign correct)       {T.reach(T.key_oracle()):.5f}   INBOX says 0.64553")
print(f"P(emitted sign == design sign) over the pool          {float((np.sign(T.S.pred_lfc[pool])==np.sign(T.H[pool])).mean()):.5f}   INBOX says 0.945")
print(f"design-sign accuracy over the pool                    {float(hcorr[pool].mean()):.5f}   INBOX says 0.537")
Z=np.load(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),"work","allkeys.npz"),allow_pickle=True)
print()
print(f"{'key':<40}{'reach':>9}{'H-oracle gap closed':>21}{'full-oracle gap closed':>24}")
base=T.reach(None); ho=T.reach(hcorr.astype(float)); fo=T.reach(T.key_oracle())
for k in ["|cons| (un-whitened)","|H| (whitened)","M1 logistic @ C=0.1 w=top25",
          "M2 step-additive @ C=0.01 w=invrank","M4 gbt @ it=100 w=invrank",
          "M5 mlp @ h=(16,) w=none","50/50 blend with |H|"]:
    r=T.reach(np.asarray(Z[k],dtype=np.float64))
    print(f"{k:<40}{r:>9.5f}{100*(r-base)/(ho-base):>20.1f}%{100*(r-base)/(fo-base):>23.1f}%")
