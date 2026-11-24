import sys, time, numpy as np
sys.path.insert(0, "~/Code/virtual-cell-challenge-2026/experiments/099-donors")
from dkit import Ctx, BASE_DONORS
t0=time.time()
C = Ctx("orion")
print(f"loaded orion in {time.time()-t0:.1f}s  P={C.P} G={C.G} rows(nmae)={C.rows.size}")
C.set_donors(BASE_DONORS)
print("coverage", C.coverage())
E = C.effects(); dh = C.dh_exact(E)
X,A,W = C.XAW(dh)
print(f"BASE  X={X:.6g} A={A:.6g} W={W:.6g} kappa={C.kappa(dh):.6g} "
      f"pds={C.pds(dh):.6f} nmae={C.nmae(dh):.6f} signdot={C.signdot(dh):.6f}")
Xs,As,Ws = C.XAW(C.dh_ship)
print(f"SHIP  X={Xs:.6g} A={As:.6g} W={Ws:.6g} kappa={C.kappa(C.dh_ship):.6g} "
      f"pds={C.pds(C.dh_ship):.6f} nmae={C.nmae(C.dh_ship):.6f}")
print("cos(rebuilt dh, shipped dh) per row median:",
      float(np.median(np.einsum('ij,ij->i',dh,C.dh_ship)/
            (np.linalg.norm(dh,axis=1)*np.linalg.norm(C.dh_ship,axis=1)))))
