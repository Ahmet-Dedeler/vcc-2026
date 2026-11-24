
### h1-test (a15-analogue artifact, 85 perturbations)   n=400 L=20190  n*L=8,076,000  gated=10779  blocks=85
    sum(T)/(nL) = 1.002414  ->  common-mode profile shift -0.2408%  (-3.474 millilog2)

| decile | ctrl CPM | median T | share T<100 | E[rel shortfall] | E[log2 bias] | net tilt (log2) |
|---|---|---:|---:|---:|---:|---:|
| 1 | 4.5–9.4 | 56 | 99.2% | 1.007% | -0.03482 | -0.03830 |
| 2 | 9.4–15.8 | 99 | 50.5% | 0.534% | -0.01257 | -0.01604 |
| 3 | 15.8–24.1 | 159 | 0.9% | 0.324% | -0.00549 | -0.00896 |
| 4 | 24.1–34.5 | 234 | 0.0% | 0.220% | -0.00439 | -0.00786 |
| 5 | 34.5–47.7 | 328 | 0.0% | 0.156% | -0.00307 | -0.00654 |
| 6 | 47.7–65.3 | 452 | 0.0% | 0.112% | -0.00161 | -0.00509 |
| 7 | 65.3–89.0 | 617 | 0.0% | 0.083% | -0.00160 | -0.00507 |
| 8 | 89.0–127.5 | 853 | 0.0% | 0.059% | -0.00085 | -0.00433 |
| 9 | 127.7–211.4 | 1286 | 0.0% | 0.039% | -0.00056 | -0.00404 |
| 10 | 211.6–7423.0 | 2520 | 0.0% | 0.019% | -0.00028 | -0.00375 |

    gated pool: median T 385; share below 100 counts **15.1%**, below 40 1.3%, below 10 0.02%, exactly zero 0.008%
    mean relative shortfall over the gated pool **0.255%**
    tilt crossover: T = nan counts = nan control CPM  (genes below it are delivered LOW, above it HIGH)

### challenge context A (a15 as shipped, 300 perturbations)   n=400 L=20109  n*L=8,043,600  gated=9929  blocks=300
    sum(T)/(nL) = 0.999174  ->  common-mode profile shift +0.0827%  (+1.193 millilog2)

| decile | ctrl CPM | median T | share T<100 | E[rel shortfall] | E[log2 bias] | net tilt (log2) |
|---|---|---:|---:|---:|---:|---:|
| 1 | 5.0–11.6 | 63 | 98.0% | 0.828% | -0.01251 | -0.01132 |
| 2 | 11.6–20.6 | 126 | 13.8% | 0.406% | -0.00600 | -0.00481 |
| 3 | 20.6–30.9 | 202 | 0.0% | 0.251% | -0.00400 | -0.00281 |
| 4 | 31.0–43.2 | 292 | 0.0% | 0.172% | -0.00249 | -0.00130 |
| 5 | 43.2–57.4 | 399 | 0.0% | 0.126% | -0.00182 | -0.00063 |
| 6 | 57.4–75.3 | 526 | 0.0% | 0.096% | -0.00138 | -0.00019 |
| 7 | 75.4–101.6 | 697 | 0.0% | 0.072% | -0.00104 | +0.00015 |
| 8 | 101.7–141.1 | 954 | 0.0% | 0.053% | -0.00076 | +0.00043 |
| 9 | 141.1–220.4 | 1376 | 0.0% | 0.036% | -0.00053 | +0.00067 |
| 10 | 220.4–5731.3 | 2495 | 0.0% | 0.019% | -0.00027 | +0.00092 |

    gated pool: median T 459; share below 100 counts **11.2%**, below 40 0.2%, below 10 0.00%, exactly zero 0.000%
    mean relative shortfall over the gated pool **0.206%**
    tilt crossover: T = 605 counts = 75.2 control CPM  (genes below it are delivered LOW, above it HIGH)

### challenge context B (a15 as shipped, 300 perturbations)   n=400 L=19946  n*L=7,978,400  gated=9626  blocks=300
    sum(T)/(nL) = 0.999727  ->  common-mode profile shift +0.0273%  (+0.394 millilog2)

| decile | ctrl CPM | median T | share T<100 | E[rel shortfall] | E[log2 bias] | net tilt (log2) |
|---|---|---:|---:|---:|---:|---:|
| 1 | 5.0–10.2 | 57 | 99.9% | 0.918% | -0.01336 | -0.01296 |
| 2 | 10.2–17.4 | 105 | 39.6% | 0.496% | -0.00936 | -0.00896 |
| 3 | 17.4–26.3 | 170 | 0.3% | 0.302% | -0.00539 | -0.00500 |
| 4 | 26.3–37.1 | 247 | 0.0% | 0.205% | -0.00322 | -0.00282 |
| 5 | 37.1–50.2 | 343 | 0.1% | 0.148% | -0.00226 | -0.00187 |
| 6 | 50.2–68.5 | 465 | 0.0% | 0.108% | -0.00156 | -0.00117 |
| 7 | 68.6–95.9 | 641 | 0.0% | 0.078% | -0.00113 | -0.00074 |
| 8 | 95.9–138.7 | 914 | 0.0% | 0.055% | -0.00079 | -0.00040 |
| 9 | 138.8–228.7 | 1376 | 0.0% | 0.036% | -0.00053 | -0.00013 |
| 10 | 228.8–7856.7 | 2777 | 0.0% | 0.017% | -0.00025 | +0.00014 |

    gated pool: median T 399; share below 100 counts **14.0%**, below 40 0.4%, below 10 0.02%, exactly zero 0.001%
    mean relative shortfall over the gated pool **0.236%**
    tilt crossover: T = 1830 counts = 229.4 control CPM  (genes below it are delivered LOW, above it HIGH)

### challenge context C (a15 as shipped, 300 perturbations)   n=400 L=20034  n*L=8,013,600  gated=10124  blocks=300
    sum(T)/(nL) = 0.996432  ->  common-mode profile shift +0.3581%  (+5.167 millilog2)

| decile | ctrl CPM | median T | share T<100 | E[rel shortfall] | E[log2 bias] | net tilt (log2) |
|---|---|---:|---:|---:|---:|---:|
| 1 | 5.0–10.1 | 58 | 99.9% | 0.878% | -0.01287 | -0.00770 |
| 2 | 10.1–17.2 | 105 | 39.3% | 0.486% | -0.00861 | -0.00345 |
| 3 | 17.2–25.9 | 170 | 0.1% | 0.300% | -0.00482 | +0.00035 |
| 4 | 25.9–36.3 | 246 | 0.0% | 0.207% | -0.00360 | +0.00157 |
| 5 | 36.3–49.7 | 338 | 0.0% | 0.150% | -0.00253 | +0.00264 |
| 6 | 49.7–67.8 | 462 | 0.0% | 0.109% | -0.00157 | +0.00359 |
| 7 | 67.8–92.7 | 632 | 0.0% | 0.080% | -0.00115 | +0.00402 |
| 8 | 92.8–134.5 | 891 | 0.0% | 0.057% | -0.00082 | +0.00435 |
| 9 | 134.5–218.4 | 1346 | 0.0% | 0.037% | -0.00054 | +0.00463 |
| 10 | 218.7–5652.1 | 2592 | 0.0% | 0.018% | -0.00027 | +0.00490 |

    gated pool: median T 395; share below 100 counts **13.9%**, below 40 0.2%, below 10 0.00%, exactly zero 0.001%
    mean relative shortfall over the gated pool **0.232%**
    tilt crossover: T = 140 counts = 17.4 control CPM  (genes below it are delivered LOW, above it HIGH)
