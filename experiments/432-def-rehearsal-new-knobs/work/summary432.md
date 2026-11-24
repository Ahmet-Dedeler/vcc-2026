
re-level check (b1 -> b2):
  D: scale asked 0.9938  B ratio 0.98805 vs s^2 0.98767 (ratio 1.00038)  rho 0.9916 -> 0.9797  cliff after 1.0103  EV-scale after 0.9998
  E: scale asked 1.0353  B ratio 1.07011 vs s^2 1.07185 (ratio 0.99837)  rho 0.9137 -> 0.9777  cliff after 1.0113  EV-scale after 1.0008
  F: scale asked 1.0680  B ratio 1.13720 vs s^2 1.14059 (ratio 0.99703)  rho 0.8586 -> 0.9764  cliff after 1.0120  EV-scale after 1.0015
| work dir | stage | seconds | peak RSS GiB | rc |
|---|---|---:|---:|---:|
| b0 | driver:repin-frozen | 1.2 | 0.24 | 0 |
| b0 | driver:pin-artifacts | 0.5 | 0.06 | 0 |
| b0 | driver:postbuild | 0.9 | 0.13 | 1 |
| b0 | driver:prep | 9.9 | 14.28 | 0 |
| b0 | repin-panel-support | 9.83 | 20.12 | |
| b0 | build | 1718.08 | 20.21 | |
| b0 | prep | 9.92 | 14.28 | |
| b1 | driver:repin-frozen | 1.2 | 0.24 | 0 |
| b1 | driver:pin-artifacts | 0.6 | 0.06 | 0 |
| b1 | driver:postbuild | 0.8 | 0.14 | 1 |
| b1 | repin-panel-support | 10.7 | 20.12 | |
| b1 | build | 1515.41 | 20.21 | |
| b2 | driver:repin-frozen | 1.2 | 0.24 | 0 |
| b2 | driver:pin-artifacts | 0.5 | 0.06 | 0 |
| b2 | driver:postbuild | 0.7 | 0.13 | 1 |
| b2 | driver:prep | 7.7 | 19.65 | 0 |
| b2 | repin-panel-support | 12.82 | 20.03 | |
| b2 | build | 2099.57 | 20.21 | |
| b2 | prep | 7.67 | 19.65 | |
| base | driver:preflight | 3.2 | 0.33 | 0 |
| base | driver:identity | 10.1 | 2.74 | 0 |
| base | driver:cd4-extract | 2.3 | 0.41 | 0 |
| base | driver:census | 0.3 | 0.1 | 0 |
| base | driver:whitener | 321.6 | 12.86 | 0 |
| dump_b0env | driver:repin-frozen | 1.0 | 0.24 | 0 |
| dump_b0env | driver:pin-artifacts | 0.5 | 0.06 | 0 |
| dump_b0env | repin-panel-support | 9.45 | 20.12 | |
| dump_b0env | build | 11.3 | 20.12 | |
| dump_h1 | driver:repin-frozen | 1.0 | 0.24 | 0 |
| dump_h1 | driver:pin-artifacts | 0.5 | 0.06 | 0 |
| dump_h1 | repin-panel-support | 10.15 | 20.12 | |
| dump_h1 | build | 11.27 | 20.12 | |
| dump_lit | driver:repin-frozen | 1.3 | 0.24 | 0 |
| dump_lit | driver:pin-artifacts | 0.5 | 0.06 | 0 |
| dump_lit | repin-panel-support | 10.21 | 20.12 | |
| dump_lit | build | 11.48 | 20.12 | |
| dump_stack | driver:repin-frozen | 1.5 | 0.24 | 0 |
| dump_stack | driver:pin-artifacts | 0.6 | 0.06 | 0 |
| dump_stack | repin-panel-support | 10.34 | 20.12 | |
| dump_stack | build | 12.23 | 20.12 | |
| ef_b1 | driver:repin-frozen | 1.1 | 0.24 | 0 |
| ef_b1 | driver:pin-artifacts | 0.5 | 0.06 | 0 |
| ef_b1 | driver:postbuild | 0.7 | 0.13 | 1 |
| ef_b1 | repin-panel-support | 10.03 | 20.06 | |
| ef_b1 | build | 2045.26 | 20.15 | |
| ef_dump_stack | driver:repin-frozen | 1.1 | 0.24 | 0 |
| ef_dump_stack | driver:pin-artifacts | 0.6 | 0.06 | 0 |
| ef_dump_stack | repin-panel-support | 9.88 | 20.06 | |
| ef_dump_stack | build | 11.46 | 20.06 | |
| neg_abcctx | driver:repin-frozen | 1.0 | 0.24 | 0 |
| neg_abcctx | driver:pin-artifacts | 0.5 | 0.06 | 0 |
| neg_abcctx | repin-panel-support | 9.69 | 20.12 | |
| neg_abcctx | build | 9.63 | 20.12 | |
| neg_abcdump | driver:repin-frozen | 1.0 | 0.24 | 0 |
| neg_abcdump | driver:pin-artifacts | 0.5 | 0.06 | 0 |
| neg_abcdump | repin-panel-support | 9.48 | 20.12 | |
| neg_abcdump | build | 9.64 | 20.12 | |
| neg_abckey | driver:repin-frozen | 1.0 | 0.24 | 0 |
| neg_abckey | driver:pin-artifacts | 0.5 | 0.06 | 0 |
| neg_abckey | repin-panel-support | 9.9 | 20.12 | |
| neg_abckey | build | 9.9 | 20.12 | |

b0 rc dump (900 blocks):
  D: mde_lambda med 0.9570 q10 0.9315 min 0.0342 frac0 0.000; zero_repaired med 0.0000 q10 0.0000 min 0.0000 frac0 0.983
  E: mde_lambda med 1.0000 q10 1.0000 min 0.0117 frac0 0.000; zero_repaired med 0.0000 q10 0.0000 min 0.0000 frac0 0.987
  F: mde_lambda med 1.0000 q10 1.0000 min 0.0010 frac0 0.000; zero_repaired med 0.0000 q10 0.0000 min 0.0000 frac0 0.980

b0 ntclevel: D: Cp46 14.4398 rho 0.9420 cliff 1.0303 scale 1.0196; E: Cp46 14.8057 rho 0.8716 cliff 1.0711 scale 1.0600; F: Cp46 16.4397 rho 0.8201 cliff 1.1043 scale 1.0928

b1 rc dump (900 blocks):
  D: mde_lambda med 0.9565 q10 0.9316 min 0.0342 frac0 0.000; zero_repaired med 0.0000 q10 0.0000 min 0.0000 frac0 0.983
  E: mde_lambda med 1.0000 q10 1.0000 min 0.0117 frac0 0.000; zero_repaired med 0.0000 q10 0.0000 min 0.0000 frac0 0.987
  F: mde_lambda med 1.0000 q10 1.0000 min 0.0010 frac0 0.000; zero_repaired med 0.0000 q10 0.0000 min 0.0000 frac0 0.980

b1 ntclevel: D: Cp46 14.4398 rho 0.9916 cliff 1.0042 scale 0.9938; E: Cp46 14.8057 rho 0.9137 cliff 1.0462 scale 1.0353; F: Cp46 16.4397 rho 0.8586 cliff 1.0792 scale 1.0680

b2 rc dump (900 blocks):
  D: mde_lambda med 0.9561 q10 0.9315 min 0.0342 frac0 0.000; zero_repaired med 0.0000 q10 0.0000 min 0.0000 frac0 0.983; zs_n_moved med 40.0000 q10 30.0000 min 12.0000 frac0 0.000
  E: mde_lambda med 1.0000 q10 1.0000 min 0.0117 frac0 0.000; zero_repaired med 0.0000 q10 0.0000 min 0.0000 frac0 0.987; zs_n_moved med 31.0000 q10 22.0000 min 14.0000 frac0 0.000
  F: mde_lambda med 1.0000 q10 1.0000 min 0.0010 frac0 0.000; zero_repaired med 0.0000 q10 0.0000 min 0.0000 frac0 0.980; zs_n_moved med 57.0000 q10 44.0000 min 22.0000 frac0 0.000

b2 ntclevel: D: Cp46 14.4398 rho 0.9797 cliff 1.0103 scale 0.9998; E: Cp46 14.8057 rho 0.9777 cliff 1.0113 scale 1.0008; F: Cp46 16.4397 rho 0.9764 cliff 1.0120 scale 1.0015

ef_b1 rc dump (900 blocks):
  D: mde_lambda med 0.9580 q10 0.9305 min 0.8662 frac0 0.000; zero_repaired med 0.0000 q10 0.0000 min 0.0000 frac0 0.987; zs_n_moved med 51.0000 q10 32.0000 min 12.0000 frac0 0.000
  E: mde_lambda med 1.0000 q10 0.9997 min 0.0000 frac0 0.003; zero_repaired med 0.0000 q10 0.0000 min 0.0000 frac0 0.987; zs_n_moved med 37.0000 q10 24.9000 min 13.0000 frac0 0.000
  F: mde_lambda med 1.0000 q10 1.0000 min 0.6309 frac0 0.000; zero_repaired med 0.0000 q10 0.0000 min 0.0000 frac0 0.977; zs_n_moved med 64.0000 q10 46.9000 min 21.0000 frac0 0.000

ef_b1 ntclevel: D: Cp46 14.4398 rho 1.0712 cliff 0.9662 scale 0.9561; E: Cp46 14.8057 rho 0.9925 cliff 1.0038 scale 0.9934; F: Cp46 16.4397 rho 0.9351 cliff 1.0341 scale 1.0234
