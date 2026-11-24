#!/usr/bin/env bash
# 700 day-of PIE step for def_day's DEF_DAY_DUMP_HOOK (690g interface, [17adab]):
#   bash dumphook700.sh <bundle> <work> <stack_dump.npz> <identity.json> <outdir>
# Writes <outdir>/<basename of stack_dump> with bulk = the stack dump's bulk with PIE's direction blended in on SEEN
# targets (norm-preserving, w = PIE700_W, default 0.25 = 702's preregistered w*), and carriers/perts/genes/env copied
# unchanged (make_pie_dump700.py --keep-env). Provenance goes to <outdir>/<basename>.pie700.json.
# Per context: identity call (CCLE SAMPLE_ID) -> DepMap ID (cBioPortal clinical table) -> Cellosaurus accession
# (work/ccle_cellosaurus_map.tsv) -> PIE context text (data/external/pie/local_sources/context_text_ccle, precomputed, so no
# API). A context with no call, a call below the identity floor (PIE700_ID_FLOOR, default 0.65 = 590's), or no accession or
# text embedding ABSTAINS: it contributes no PIE row. If every context abstains, the dump is copied unchanged (a no-op).
# Offline: HF_HUB_OFFLINE=1; needs the warm PIE cache (checkpoint, assets, evidence 87a01957... from the cloud build).
# Logs the seen count and the abstentions per context to stdout (def_day captures it as <work>/dumphook.log).
set -euo pipefail
REPO=~/Code/virtual-cell-challenge-2026
X=$REPO/experiments/700-pie-transfer; PY=$X/.venv/bin/python
BUNDLE="$1"; WORK="$2"; DUMP="$3"; IDJ="$4"; OUT="$5"
W="${PIE700_W:-0.25}"; FLOOR="${PIE700_ID_FLOOR:-0.65}"
export PIE_DATA_ROOT=$REPO/data/external/pie/data PIE_RUNS_ROOT=$REPO/data/external/pie/runs PIE_CACHE_DIR=$REPO/data/external/pie/cache
export HF_HUB_OFFLINE=1
CT=$REPO/data/external/pie/local_sources/context_text_ccle
T="$WORK/pie700"; mkdir -p "$T" "$OUT"
ts(){ date -u +%FT%TZ; }
echo "$(ts) dumphook700: bundle $BUNDLE dump $(basename "$DUMP") w $W floor $FLOOR"
# 1. identity -> accession per context (abstentions logged)
"$PY" - "$BUNDLE" "$IDJ" "$FLOOR" "$X/work/ccle_cellosaurus_map.tsv" "$REPO/data/external/ccle_broad_2025/data_clinical_sample.txt" \
      "$CT/meta.json" "$T/ctxmap.json" <<'PY'
import csv, json, sys
bundle, idj, floor, cmap, clin, ctmeta, out = sys.argv[1:]
ident = json.load(open(idj)); floor = float(floor)
calls = ident.get("call", {}) or {}; top1 = ident.get("min_top1_score", {}) or {}
s2d = {}
with open(clin) as f:
    rows = [ln for ln in f if not ln.startswith("#")]
for r in csv.DictReader(rows, delimiter="\t"):
    s2d[r["SAMPLE_ID"]] = r["DEPMAP_ID"]
d2a = {r["depmap_id"]: r["accession"] for r in csv.DictReader(open(cmap), delimiter="\t") if r["accession"]}
keys = set(json.load(open(ctmeta))["keys"])
ctxs = json.load(open(f"{bundle}/manifest.json")).get("contexts") if False else ident.get("contexts")
res = {}
for c in ctxs:
    call = calls.get(c); sc = top1.get(c)
    why = None
    if not call: why = "no identity call"
    elif sc is not None and float(sc) < floor: why = f"identity score {float(sc):.3f} < {floor}"
    else:
        acc = d2a.get(s2d.get(call, ""), "")
        if not acc: why = f"no Cellosaurus accession for {call}"
        elif acc not in keys: why = f"no precomputed context text for {acc}"
    res[c] = {"call": call, "score": sc, "accession": None if why else d2a[s2d[call]], "abstain": why}
    print(f"  context {c}: call {call} score {sc} -> " + (f"ABSTAIN ({why})" if why else f"PIE context {res[c]['accession']}"))
json.dump(res, open(out, "w"), indent=1)
PY
USE="$("$PY" -c "import json,sys; m=json.load(open(sys.argv[1])); print(' '.join(f'{c}={v[\"accession\"]}' for c,v in m.items() if not v['abstain']))" "$T/ctxmap.json")"
if [ -z "$USE" ]; then echo "$(ts) every context abstained: copying the stack dump unchanged"; cp "$DUMP" "$OUT/$(basename "$DUMP")"; exit 0; fi
# 2. controls -> log1p(CP10k) h5ads named by accession; controls-only prep on the bundle gene axis ∩ PIE gene_text
mkdir -p "$T/ctrl"
"$PY" "$X/work/make_ctrl_h5ad700.py" "$BUNDLE" "$T/ctrl" $USE
"$PY" - "$BUNDLE" "$PIE_DATA_ROOT" "$T/genes.txt" <<'PY'
import csv, json, sys, glob
bundle, root, out = sys.argv[1:]
gt = set(json.load(open(glob.glob(f"{root}/hf/datasets/arcinstitute/pie_sources/*/gene_text/meta.json")[0]))["keys"])
genes = [r["gene_name"] for r in csv.DictReader(open(f"{bundle}/gene_names.csv"))]
open(out, "w").write("\n".join(g for g in genes if g in gt) + "\n")
PY
rm -rf "$T/prep"
"$X/.venv/bin/pie" prep dataset=arc_vcc_25 name=defday controls_only=true "h5ad=$T/ctrl/ctrl_*.h5ad" output_dir="$T/prep" genes="$T/genes.txt"
# 3. query = the bundle's panel x the non-abstaining contexts; offline inference
"$PY" - "$BUNDLE" "$T/ctxmap.json" "$T/query.json" <<'PY'
import csv, json, sys
bundle, cm, out = sys.argv[1:]
t = [r["target_gene"] for r in csv.DictReader(open(f"{bundle}/pert_counts.csv")) if r["target_gene"] != "non-targeting"]
m = json.load(open(cm))
json.dump({f"defday.{v['accession']}": t for v in m.values() if not v["abstain"]}, open(out, "w"))
PY
/usr/bin/time -l "$PY" "$X/infer700.py" --ckpt "$PIE_RUNS_ROOT/replogle_xdataset/best_auprc.ckpt" --prep "$T/prep" \
  --query "$T/query.json" --out "$T/pred.parquet" --context-text "$CT" 2>&1 | tail -4
# 4. the blend (seen rule = PIE's pinned train split; env/perts/genes/carriers copied unchanged)
CTXARGS=""; for kv in $USE; do CTXARGS="$CTXARGS --ctx $kv"; done
"$PY" "$X/make_pie_dump700.py" --base "$DUMP" --pred "$T/pred.parquet" --dataset defday $CTXARGS --bundle "$BUNDLE" \
  --w "$W" --seen "$X/work/pie_train_perts.json" --keep-env --out "$OUT/$(basename "$DUMP")"
echo "$(ts) dumphook700 done: $(cat "$OUT/$(basename "$DUMP").pie700.json" | head -c 400)"
