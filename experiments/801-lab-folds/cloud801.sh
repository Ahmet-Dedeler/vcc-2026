#!/usr/bin/env bash
# 801 -- the EC2 side of the lab folds (run by labworker.sh from /work/repo; never on the Mac). One subcommand per job line.
#
#   cloud801.sh verify            sha256 of every mirrored input vs the Mac's manifests (work/manifest/*.txt); writes
#                                 work/cloud/verify.ok or verify.fail. Every other job waits for it and aborts on a fail.
#   cloud801.sh prep <L>          one WORKING panel (L1..L6): bundle controls (copies of 199's relabelled files, sha256
#                                 checked) and CD4 extraction (cd4_801.py; every non-json file's sha256 == the Mac's).
#   cloud801.sh fold <L> <held>   one fold: the held-out operator fit (scripts/build_union_whitener.py, 534's flags,
#                                 --max-rss-gib 14), kept at s3://$BUCKET/output/lab/ops801/<L>/ the moment it exists (a later
#                                 run fetches it instead of re-fitting), then the dump (dump801.py, closure recorded).
#   cloud801.sh refit-efp <held>  re-fit 238's efp held-out operator on EC2 with --expect-sha256 (the Mac's file): does a
#                                 fit reproduce across machines bit for bit?
#   cloud801.sh ka-efp <held>     dump801 --registry efp --ka522 (dump522's known answers on the uploaded efp operators),
#                                 then engine801 ka efp/<held> (546o, 551, 522 rows) on the fresh dump; marker ka-efp-<held>.ok.
#   cloud801.sh refit-study <held>  after ka-efp <held>: re-fit efp's operator on EC2, dump with it, refit801.py's diffs.
#   cloud801.sh smoke <tag> <fold>...  smoke801.py on those folds, only once all three ka-efp markers exist (else abort).
#   cloud801.sh designs-L <L> <held>   designs801.py L (h2.5c1 beside h2c1; KA: the re-run dump == the registry dump)
#   cloud801.sh designs-ABC <held>     designs801.py ABC (h2.5c1.5 beside h2c1.5; KA: h2c1.5 == 518's fST)
#   cloud801.sh fit-abc <held>         A/B/C's held-out operator by the day recipe (--basis-prune essential, 12k cap), on S3
#   cloud801.sh abcp <held>            dump518p 062 (KA == 518's dump) -> dump518p day (ABCp/<held>) -> ka804 on ABCp/<held>
# Threads: OPENBLAS/OMP/MKL = $T801_THREADS (default 5) per job, so 6 parallel jobs fit 32 vCPU.
set -Eeuo pipefail
cd /work/repo
E=experiments/801-lab-folds
W=$E/work
C=$W/cloud
mkdir -p "$C"
T="${T801_THREADS:-5}"
export OPENBLAS_NUM_THREADS="$T" OMP_NUM_THREADS="$T" MKL_NUM_THREADS="$T"
FOLDS="k562_gwps orion_hct116 orion_hek293t"
say() { echo "[cloud801 $(date -u +%FT%TZ)] $*"; }
tag_ok() {                            # any new panel tag (the kit), never the sealed reserve
  [[ "$1" =~ ^[A-Za-z][A-Za-z0-9]{0,15}$ && "$1" != L7 && "$1" != L8 && "$1" != sealed ]]; }
# Foreign modules mkdir their own work/ at import (plolo: experiments/162-.../work). Those dirs exist on the Mac, not in the
# EC2 checkout, and guard801 refuses a NEW directory outside 801; so the empty dirs are made here, on EC2 only. Any WRITE
# into them still raises (guard801), so nothing foreign is cached.
for d in experiments/*/; do mkdir -p "${d}work"; done
export GUARD801_ALLOW_MKDIR=1         # nested ones (218's work/plolo ...): guard801 then allows NEW directories, never files
others() { case "$1" in k562_gwps) echo orion_hct116,orion_hek293t;; orion_hct116) echo k562_gwps,orion_hek293t;;
  orion_hek293t) echo k562_gwps,orion_hct116;; *) return 1;; esac; }
wait_verify() {
  local i=0
  while [[ ! -f $C/verify.ok && ! -f $C/verify.fail ]]; do sleep 10; i=$((i+1)); [[ $i -lt 360 ]] || { say "verify never finished"; exit 4; }; done
  [[ -f $C/verify.ok ]] || { say "inputs failed verification: $(cat $C/verify.fail)"; exit 5; }
}

cmd_verify() {
  local bad=0
  for m in $W/manifest/local_sha256_backup.txt $W/manifest/local_sha256_defbundle.txt $W/manifest/local_sha256_backup2.txt \
           $W/manifest/local_sha256_abc.txt; do
    [[ -f "$m" ]] || continue
    while read -r sha path; do
      [[ -n "$path" ]] || continue
      got="$(sha256sum "$path" 2>/dev/null | cut -d' ' -f1 || true)"
      if [[ "$got" == "$sha" ]]; then echo "ok $path"; else echo "MISMATCH $path (mac $sha, ec2 ${got:-missing})"; bad=$((bad+1)); fi
    done < "$m"
  done | tee $C/verify.log
  bad="$(grep -c '^MISMATCH' $C/verify.log || true)"
  if [[ "$bad" -eq 0 ]]; then say "verify: $(grep -c '^ok' $C/verify.log) files identical"; touch $C/verify.ok
  else grep '^MISMATCH' $C/verify.log > $C/verify.fail; say "verify FAILED ($bad)"; exit 6; fi
}

cmd_prep() {                          # bundle controls + CD4 banks for one WORKING panel; marker $C/<L>.prep.ok
  local L="$1"
  tag_ok "$L" || { say "bad or sealed panel tag $L"; exit 2; }
  wait_verify
  trap 'touch $C/'"$L"'.prep.fail' ERR
  local B=$W/$L/bundle
  if [[ -f $B/frames_src.txt ]]; then      # kit --frames <dir>: that bundle's own controls are the frames (no relabel, not synthetic)
    local FS; FS="$(cat $B/frames_src.txt)"
    cp "$FS/manifest.json" $B/manifest.json
    for c in $(jq -r '.contexts[]' "$FS/manifest.json"); do cp "$FS/context_$c.h5ad" $B/; done
    ( cd $B && sha256sum manifest.json $(for c in $(jq -r '.contexts[]' manifest.json); do printf 'context_%s.h5ad ' $c; done) ) > $C/$L-frames.sha256
    diff <(sort -k2 $B/frames.sha256) <(sort -k2 $C/$L-frames.sha256) > $C/$L-frames.diff || { say "$L: frames from $FS differ from the Mac's sha256"; exit 7; }
    say "$L: frames = $FS's own controls ($(jq -r '.contexts|join(",")' $B/manifest.json)), sha256 == the Mac's"
  else
  if [[ ! -f $B/manifest.json ]]; then      # the kit's bundles carry only the two csv files: A/B/C's manifest, relabelled D/E/F
    uv run python -c "import json,sys; m=json.load(open('data/challenge/manifest.json')); R={'A':'D','B':'E','C':'F'}; m.update(partition='test', panel_id='vcc2026-synthetic-lab801kit-'+sys.argv[1], contexts=[R[c] for c in m['contexts']], per_context={R[k]: v for k, v in m['per_context'].items()}); open(sys.argv[2],'w').write(json.dumps(m, indent=2)+'\\n')" "$L" "$B/manifest.json"
  fi
  for x in D E F; do cp experiments/199-def-runbook/work/def-bundle/context_$x.h5ad $B/; done
  ( cd $B && sha256sum context_D.h5ad context_E.h5ad context_F.h5ad ) | while read -r s f; do
    grep -q "^$s  experiments/199-def-runbook/work/def-bundle/$f\$" $W/manifest/local_sha256_defbundle.txt || { say "$L: $f differs"; exit 7; }
  done
  say "$L: bundle controls in place"
  fi
  uv run python $E/cd4_801.py "$L"
  ( cd $W && find $L/cd4 -type f | sort | xargs sha256sum ) > $C/$L-cd4.sha256
  # json records carry absolute paths (Mac vs EC2), so every non-json file is compared. A-panels have no Mac extraction
  # (local-job hold); the L folds proved EC2 == Mac on 6 panels, so their hashes are recorded only.
  if [[ ! "$L" =~ ^L[1-6]$ ]]; then say "$L: CD4 banks extracted ($(wc -l < $C/$L-cd4.sha256) files; no Mac reference beyond L1-L6)"
  elif diff <(grep " $L/cd4/" $W/manifest/local_sha256_cd4.txt | grep -v "\.json$" | sort -k2) <(grep -v "\.json$" $C/$L-cd4.sha256 | sort -k2) > $C/$L-cd4.diff; then
    say "$L: CD4 banks (every non-json file) byte-identical to the Mac's extraction"
  else say "$L: CD4 banks DIFFER from the Mac's extraction"; touch $C/$L.prep.fail; exit 8; fi
  flock $C/frames.lock uv run python $E/frames801.py "$L"   # serialised: the first one writes the shared A/B/C frame cache
  # the whole dump path with an identity stub operator (seconds; numbers meaningless): a path or guard error surfaces here,
  # before 18 fits are spent on it
  GUARD801_CLOSURE=$C/$L-stubclosure.txt uv run python $E/discover801.py dumpL "$L" k562_gwps > $C/$L-stub.log 2>&1
  say "$L: stub dump path OK ($(grep -c . $C/$L-stubclosure.txt) files opened)"
  echo "$L $(date -u +%FT%TZ)" | aws s3 cp - "s3://${BUCKET}/output/lab/ops801/$L/probe.txt" --only-show-errors   # the role may write there
  say "$L: S3 write probe OK"
  touch $C/$L.prep.ok
}

cmd_fold() {                          # one fold: its held-out operator (fitted, or fetched from output/lab/ops801 if a run already
  local L="$1" h="$2"                 # kept it), uploaded to S3 at once; then the dump (dump801, closure recorded)
  tag_ok "$L" || { say "bad or sealed panel tag $L"; exit 2; }
  wait_verify
  local i=0
  while [[ ! -f $C/$L.prep.ok ]]; do [[ -f $C/$L.prep.fail ]] && exit 8; sleep 10; i=$((i+1)); [[ $i -lt 180 ]] || exit 4; done
  local B=$W/$L/bundle O=$W/$L/ops S3O="s3://${BUCKET}/output/lab/ops801/$L"
  mkdir -p $O
  if aws s3 ls "$S3O/held_$h.json" >/dev/null 2>&1; then
    aws s3 cp "$S3O/held_$h.npz" $O/ --only-show-errors; aws s3 cp "$S3O/held_$h.json" $O/ --only-show-errors
    [[ "$(sha256sum $O/held_$h.npz | cut -d' ' -f1)" == "$(jq -r .sha256 $O/held_$h.json)" ]] || { say "$L/$h: kept operator corrupt"; exit 10; }
    say "$L/$h: operator fetched from $S3O (sha256 $(jq -r .sha256 $O/held_$h.json))"
  else
    local t0; t0=$(date +%s)
    local ft="${FIT_THREADS:-$T}"           # the kit fits at 8 threads (ABCp's setting), dumps at T
    OPENBLAS_NUM_THREADS=$ft OMP_NUM_THREADS=$ft MKL_NUM_THREADS=$ft uv run python scripts/build_union_whitener.py --bundle $B \
      --donors "$(others $h)" --basis-prune essential --out $O/held_$h.npz --scratch $O/.scratch_$h --max-rss-gib 14 > $C/$L-fit-$h.log 2>&1
    local kept=1
    { aws s3 cp $O/held_$h.npz "$S3O/" --only-show-errors && aws s3 cp $O/held_$h.json "$S3O/" --only-show-errors; } || kept=0
    say "$L/$h: fit done in $(( $(date +%s) - t0 )) s; sha256 $(sha256sum $O/held_$h.npz | cut -d' ' -f1); kept at $S3O/: $kept"
  fi
  GUARD801_CLOSURE=$C/$L-closure-$h.txt uv run python $E/dump801.py --tag "$L" --bundle $B --op none \
    --cd4-root $W/$L/cd4/cd4_compact --fold-ops held --held-ops $O $([[ -f $B/frames_src.txt ]] || echo --synthetic) "$h" > $C/$L-dump-$h.log 2>&1
  say "$L/$h: dump: $(tail -1 $C/$L-dump-$h.log)"
  ( cd $W && sha256sum $L/ops/held_$h.npz $L/ops/held_$h.json dumps/dump708-$L-$h.npz dumps/dump708-$L-$h.json ) > $C/$L-$h-outputs.sha256
  aws s3 ls "$S3O/held_$h.npz" >/dev/null 2>&1 && rm -f $O/held_$h.npz || say "$L/$h: operator NOT on S3; left on the instance disk"
}

cmd_refit_efp() {
  local h="$1" want
  wait_verify
  want="$(grep "efp__held_$h.npz" $W/manifest/local_sha256_efp_ops.txt | cut -d' ' -f1)"
  mkdir -p $W/efp-refit
  set +e
  uv run python scripts/build_union_whitener.py --bundle experiments/238-essential-free-panel/work/ef-bundle \
    --donors "$(others $h)" --basis-prune essential --out $W/efp-refit/efp__held_$h.npz --scratch $W/efp-refit/.scratch_$h \
    --max-rss-gib 14 --expect-sha256 "$want" > $C/refit-efp-$h.log 2>&1
  local rc=$?
  set -e
  say "refit efp held_$h: rc $rc; $(grep -E 'MATCH|sha256' $C/refit-efp-$h.log | tail -2 | tr '\n' ' ')"
  rm -f $W/efp-refit/efp__held_$h.npz
}

cmd_ka_efp() {
  local h="$1"
  wait_verify
  for x in D E F; do cp -n experiments/199-def-runbook/work/def-bundle/context_$x.h5ad experiments/238-essential-free-panel/work/ef-bundle/ 2>/dev/null || true; done
  flock $C/efp.bind.lock uv run python $E/frames801.py efp      # the three ka-efp jobs share efp's frame caches
  GUARD801_CLOSURE=$C/efp-closure-$h.txt uv run python $E/dump801.py --tag efp --registry --ka522 "$h" > $C/efp-dump-$h.log 2>&1
  say "efp dump $h: $(grep -c 'KA ' $C/efp-dump-$h.log) KA lines; $(tail -1 $C/efp-dump-$h.log)"
  uv run python $E/engine801.py ka "efp/$h" > $C/efp-engine-$h.log 2>&1
  say "efp engine $h: $(head -1 $C/efp-engine-$h.log)"
  touch $C/ka-efp-$h.ok
  touch $C/ka-efp-$h.ok
}

cmd_refit_study() {                   # an EC2 refit of efp's held-out operator, its dump, and refit801's price differences
  local h="$1"
  wait_verify
  local i=0
  while [[ ! -f $C/ka-efp-$h.ok ]]; do sleep 10; i=$((i+1)); [[ $i -lt 180 ]] || exit 9; done
  local O=$W/efprefit-ops
  mkdir -p $O
  uv run python scripts/build_union_whitener.py --bundle experiments/238-essential-free-panel/work/ef-bundle \
    --donors "$(others $h)" --basis-prune essential --out $O/held_$h.npz --scratch $O/.scratch_$h --max-rss-gib 14 > $C/efprefit-fit-$h.log 2>&1
  say "efp refit held_$h: sha256 $(sha256sum $O/held_$h.npz | cut -d' ' -f1) (Mac: $(grep "efp__held_$h.npz" $W/manifest/local_sha256_efp_ops.txt | cut -d' ' -f1))"
  uv run python $E/dump801.py --tag efprefit --bundle experiments/238-essential-free-panel/work/ef-bundle --op none \
    --cd4-root experiments/238-essential-free-panel/work/cd4/efp/cd4_compact --fold-ops held --held-ops $O --synthetic "$h" > $C/efprefit-dump-$h.log 2>&1
  uv run python $E/refit801.py "$h" > $C/refit801-$h.log 2>&1
  say "refit study $h: $(grep -c . $C/refit801-$h.log) lines"
  rm -f $O/held_$h.npz
}

cmd_designs_L() {                     # the fold re-run with --designs h2.5c1 (designs801.py L), after its panel's prep
  local L="$1" h="$2" i=0
  wait_verify
  while [[ ! -f $C/$L.prep.ok ]]; do [[ -f $C/$L.prep.fail ]] && exit 8; sleep 10; i=$((i+1)); [[ $i -lt 180 ]] || exit 4; done
  local want; want="$(jq -r ".folds[\"$L/$h\"].operator_sha256" $W/registry.json)"
  [[ "$(sha256sum $W/$L/ops/held_$h.npz | cut -d' ' -f1)" == "$want" ]] || { say "$L/$h: operator sha differs from the registry"; exit 11; }
  GUARD801_CLOSURE=$C/designs-$L-$h-closure.txt uv run python $E/designs801.py L "$L" "$h" > $C/designs-$L-$h.log 2>&1
  say "$L/$h designs: $(tail -1 $C/designs-$L-$h.log)"
}

cmd_designs_ABC() {                   # A/B/C rules through lib801.bind("abc"); serialised with every other frame-cache writer
  local h="$1"
  wait_verify
  GUARD801_CLOSURE=$C/designs-ABC-$h-closure.txt flock $C/frames.lock uv run python $E/designs801.py ABC "$h" > $C/designs-ABC-$h.log 2>&1
  say "ABC/$h designs: $(tail -1 $C/designs-ABC-$h.log)"
}

cmd_fit_abc() {                       # A/B/C's held-out operator by the DAY recipe (the L folds' flags), kept on S3
  local h="$1"
  wait_verify
  local O=$W/ABCp/ops S3O="s3://${BUCKET}/output/lab/ops801/ABCp"
  mkdir -p $O
  echo "probe $(date -u +%FT%TZ)" | aws s3 cp - "$S3O/probe.txt" --only-show-errors
  if [[ -f $O/held_$h.json && -f $O/held_$h.npz ]] && [[ "$(sha256sum $O/held_$h.npz | cut -d' ' -f1)" == "$(jq -r .sha256 $O/held_$h.json)" ]]; then
    say "ABCp/$h: operator already fitted (mirror copy of $S3O), sha256 $(jq -r .sha256 $O/held_$h.json)"
    touch $C/ABCp-fit-$h.ok; return 0
  fi
  local t0; t0=$(date +%s)
  uv run python scripts/build_union_whitener.py --bundle data/challenge --donors "$(others $h)" --basis-prune essential \
    --out $O/held_$h.npz --scratch $O/.scratch_$h --max-rss-gib 14 > $C/ABCp-fit-$h.log 2>&1
  aws s3 cp $O/held_$h.npz "$S3O/" --only-show-errors; aws s3 cp $O/held_$h.json "$S3O/" --only-show-errors
  say "ABCp/$h: fit done in $(( $(date +%s) - t0 )) s; sha256 $(sha256sum $O/held_$h.npz | cut -d' ' -f1); kept at $S3O/"
  touch $C/ABCp-fit-$h.ok
}

cmd_abcp() {                          # dump518p with 062's operators (KA vs 518's dumps), then with the day operator, then ka804
  local h="$1" i=0
  wait_verify
  for o in 062 day; do                # frames from 420's caches, as 518's dumps saw them
    mkdir -p $W/lib/abcp_plolo/$o; cp experiments/420-bulk-stack/work/plolo/bctrl-challenge-*.npy $W/lib/abcp_plolo/$o/
  done
  GUARD801_CLOSURE=$C/abc062-closure-$h.txt uv run python $E/dump518p.py "$h" 062 > $C/dump518p-062-$h.log 2>&1
  say "ABC/$h with 062's operator: $(grep 'KA every array' $C/dump518p-062-$h.log | tail -1)"
  while [[ ! -f $C/ABCp-fit-$h.ok ]]; do sleep 10; i=$((i+1)); [[ $i -lt 360 ]] || exit 4; done
  GUARD801_CLOSURE=$C/abcp-closure-$h.txt uv run python $E/dump518p.py "$h" day > $C/dump518p-day-$h.log 2>&1
  say "ABCp/$h: $(tail -1 $C/dump518p-day-$h.log)"
  rm -f $W/ABCp/ops/held_$h.npz
  mkdir -p /work/ka804 $W/ka-abcp
  KA804_OUT=/work/ka804 uv run python experiments/804-known-answers/ka804.py "ABCp/$h" > $C/ka804-ABCp-$h.log 2>&1
  cp -R /work/ka804/. $W/ka-abcp/
  say "ka804 ABCp/$h: $(tail -1 $C/ka804-ABCp-$h.log)"
}

cmd_designs_ABCp() {                  # 062 check (row sha == ABC__ file) then the day-operator designs, serialised as designs-ABC
  local h="$1"
  wait_verify
  GUARD801_CLOSURE=$C/designs-ABC062-$h-closure.txt flock $C/frames.lock uv run python $E/designs801.py ABC062check "$h" > $C/designs-ABC062-$h.log 2>&1
  say "ABC062/$h: $(tail -1 $C/designs-ABC062-$h.log)"
  GUARD801_CLOSURE=$C/designs-ABCp-$h-closure.txt flock $C/frames.lock uv run python $E/designs801.py ABCp "$h" > $C/designs-ABCp-$h.log 2>&1
  say "ABCp/$h designs: $(tail -1 $C/designs-ABCp-$h.log)"
}

cmd_prep_labc() {                     # LABC: A/B/C's panel through the L-fold recipe. CD4 banks by cd4_801 (as L1-L6 got)
  wait_verify
  trap 'touch $C/LABC.prep.fail' ERR
  uv run python $E/cd4_801.py LABC
  ( cd $W && find LABC/cd4/cd4_compact -type f ! -name '*.json' | sort | while read -r f; do
      printf '%s  %s\n' "$(sha256sum "$f" | cut -d' ' -f1)" "${f#LABC/cd4/cd4_compact/}"; done ) > $C/LABC-cd4.sha256
  ( cd data/external/cd4t/cloud-compact/cd4_compact && find . -type f ! -name '*.json' | sed 's|^\./||' | sort | while read -r f; do
      printf '%s  %s\n' "$(sha256sum "$f" | cut -d' ' -f1)" "$f"; done ) > $C/ABCprod-cd4.sha256
  if diff $C/LABC-cd4.sha256 $C/ABCprod-cd4.sha256 > $C/LABC-cd4.diff; then say "LABC: CD4 banks == A/B/C's production banks (every non-json file)"
  else say "LABC: CD4 banks differ from A/B/C's production banks in $(grep -c '^[<>]' $C/LABC-cd4.diff) lines (recorded; the L recipe's own extraction is used)"; fi
  touch $C/LABC.prep.ok
}

cmd_labc() {                          # one LABC fold: dump801 (D/E/F rules) with ABCp's day operator for this panel, then ka804
  local h="$1" i=0
  wait_verify
  while [[ ! -f $C/LABC.prep.ok ]]; do [[ -f $C/LABC.prep.fail ]] && exit 8; sleep 10; i=$((i+1)); [[ $i -lt 180 ]] || exit 4; done
  [[ "$(sha256sum $W/ABCp/ops/held_$h.npz | cut -d' ' -f1)" == "$(jq -r .sha256 $W/ABCp/ops/held_$h.json)" ]] || { say "LABC/$h: operator sha"; exit 11; }
  GUARD801_CLOSURE=$C/LABC-closure-$h.txt flock $C/frames.lock uv run python $E/dump801.py --tag LABC --bundle data/challenge \
    --op none --cd4-root $W/LABC/cd4/cd4_compact --fold-ops held --held-ops $W/ABCp/ops "$h" > $C/LABC-dump-$h.log 2>&1
  say "LABC/$h: dump: $(tail -1 $C/LABC-dump-$h.log)"
  ( cd $W && sha256sum dumps/dump708-LABC-$h.npz dumps/dump708-LABC-$h.json ABCp/ops/held_$h.npz ) > $C/LABC-$h-outputs.sha256
  mkdir -p /work/ka804labc-$h $W/ka-labc
  KA804_OUT=/work/ka804labc-$h uv run python experiments/804-known-answers/ka804.py "LABC/$h" > $C/ka804-LABC-$h.log 2>&1
  cp -R /work/ka804labc-$h/. $W/ka-labc/
  say "ka804 LABC/$h: $(tail -1 $C/ka804-LABC-$h.log)"
  touch $C/LABC-$h.ok
}

cmd_labcd() {                         # LABC rebuilt under tag LABCd with --designs h2.5c1,h3c1; KA vs the current LABC dump
  local h="$1" i=0
  wait_verify
  while [[ ! -f $C/LABC.prep.ok ]]; do [[ -f $C/LABC.prep.fail ]] && exit 8; sleep 10; i=$((i+1)); [[ $i -lt 180 ]] || exit 4; done
  [[ "$(sha256sum $W/ABCp/ops/held_$h.npz | cut -d' ' -f1)" == "$(jq -r .sha256 $W/ABCp/ops/held_$h.json)" ]] || exit 11
  GUARD801_CLOSURE=$C/LABCd-closure-$h.txt flock $C/frames.lock uv run python $E/dump801.py --tag LABCd --bundle data/challenge \
    --op none --cd4-root $W/LABC/cd4/cd4_compact --fold-ops held --held-ops $W/ABCp/ops --designs h2.5c1,h3c1 "$h" > $C/LABCd-dump-$h.log 2>&1
  say "LABCd/$h: dump: $(tail -1 $C/LABCd-dump-$h.log)"
  uv run python - "$W/dumps/dump708-LABCd-$h.npz" "$W/dumps/dump708-LABC-$h.npz" "$C/LABCd-ka-$h.json" <<'PY'
import hashlib, json, sys
import numpy as np
new, old, out = sys.argv[1:4]
a, b = np.load(new), np.load(old)
ka, bit = {}, {}
for k in b.files:                                     # every array of the current LABC dump (none are design arrays)
    x, y = a[k], b[k]
    if x.dtype.kind in "fc":
        both = np.isnan(x) & np.isnan(y)
        d = float(np.max(np.abs(np.where(both, 0, x - y)))) if x.size else 0.0
        if (np.isnan(x) != np.isnan(y)).any(): d = float("inf")
    else:
        d = 0.0 if np.array_equal(x, y) else float("inf")
    ka[k], bit[k] = d, bool(np.array_equal(x, y, equal_nan=x.dtype.kind in "fc"))
extra = sorted(set(a.files) - set(b.files))
assert extra == ["fin_h2.5c1", "fin_h3c1"], extra
rec = dict(new=new, old=old, new_sha256=hashlib.sha256(open(new, "rb").read()).hexdigest(),
           old_sha256=hashlib.sha256(open(old, "rb").read()).hexdigest(), max_abs=ka, bitwise=bit, extra=extra,
           fin_eq_old_fin=ka["fin"], all_bitwise=all(bit.values()), worst=max(ka.values()),
           design_row_sha256={k: hashlib.sha256(np.ascontiguousarray(a[k], np.float64).tobytes()).hexdigest() for k in extra},
           design_differs_from_fin={k: bool(not np.array_equal(a[k], a["fin"])) for k in extra})
json.dump(rec, open(out, "w"), indent=1)
print(f"KA worst |d| {rec['worst']:.1e}, all bitwise {rec['all_bitwise']}")
assert rec["worst"] < 1e-12, rec["worst"]
PY
  say "LABCd/$h KA: $(jq -r '"worst \(.worst) bitwise \(.all_bitwise) sha \(.new_sha256)"' $C/LABCd-ka-$h.json)"
}

cmd_regfrag() {                       # the kit's registry fragment for <TAG> (merged into work/registry.json on the Mac)
  local L="$1" i=0
  tag_ok "$L" || exit 2
  while :; do local n=0; for h in $FOLDS; do [[ -f $C/$L-$h-outputs.sha256 ]] && n=$((n+1)); done; [[ $n -eq 3 ]] && break
    i=$((i+1)); [[ $i -lt 360 ]] || exit 9; sleep 10; done
  uv run python - "$L" "$W" "$C" "${BUCKET}" <<'PY'
import json, sys
from pathlib import Path
L, W, C, bucket = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3]), sys.argv[4]
frag = {}
for h in ("k562_gwps", "orion_hct116", "orion_hek293t"):
    sh = {}
    for ln in (C / f"{L}-{h}-outputs.sha256").read_text().splitlines():
        s, p = ln.split(None, 1); sh[p.strip()] = s
    op = json.loads((W / L / "ops" / f"held_{h}.json").read_text())
    assert op["sha256"] == sh[f"{L}/ops/held_{h}.npz"]
    rec = json.loads((W / "dumps" / f"dump708-{L}-{h}.json").read_text())
    frag[f"{L}/{h}"] = dict(status="kit", panel=L, held=h, rule="panel kit (cloud801.sh; the L-fold recipe on a given panel file)",
                            panel_csv_sha256=__import__("hashlib").sha256((W / L / "bundle/pert_counts.csv").read_bytes()).hexdigest(),
                            dump=f"experiments/801-lab-folds/work/dumps/dump708-{L}-{h}.npz", dump_sha256=sh[f"dumps/dump708-{L}-{h}.npz"],
                            operator=f"s3://{bucket}/output/lab/ops801/{L}/held_{h}.npz", operator_sha256=op["sha256"],
                            operator_axis_sha256=op["axis_sha256"], operator_donors=op["donors"], n_targets=rec["n"],
                            n_carried=rec["n_carried"], coverage=f"{rec['n']}/300")
(W / f"registry-{L}.json").write_text(json.dumps(frag, indent=1))
print(f"registry fragment: {len(frag)} folds")
PY
  say "registry fragment work/registry-$L.json written"
}

cmd_kacmp() {                         # the kit's dry-run KA: <TAG>'s dumps and ka804 arrays vs a reference tag's, per array
  local L="$1" R="$2" i=0
  while :; do local n=0; for h in $FOLDS; do [[ -f $W/ka-$(echo $L | tr A-Z a-z)/meta/${L}__$h.json ]] && n=$((n+1)); done; [[ $n -eq 3 ]] && break
    i=$((i+1)); [[ $i -lt 360 ]] || exit 9; sleep 10; done
  uv run python - "$L" "$R" "$W" "$C" <<'PY'
import hashlib, json, sys
from pathlib import Path
import numpy as np
L, R, W, C = sys.argv[1], sys.argv[2], Path(sys.argv[3]), Path(sys.argv[4])
def cmp(a, b):
    ka, bit = {}, {}
    for k in b.files:
        x, y = a[k], b[k]
        if x.shape != y.shape: ka[k], bit[k] = float("inf"), False; continue
        if x.dtype.kind in "fc":
            nn = np.isnan(x) | np.isnan(y)
            ka[k] = float(np.max(np.abs(np.where(nn, 0, x - y)))) if x.size else 0.0
            if (np.isnan(x) != np.isnan(y)).any(): ka[k] = float("inf")
        else:
            ka[k] = 0.0 if np.array_equal(x, y) else float("inf")
        bit[k] = bool(np.array_equal(x, y, equal_nan=x.dtype.kind in "fc"))
    return ka, bit
refka = {"LABC": W / "ka-labc", "ABCp": W / "ka-abcp"}.get(R, W / f"ka-{R.lower()}")
out = {}
for h in ("k562_gwps", "orion_hct116", "orion_hek293t"):
    rec = {}
    a, b = np.load(W / "dumps" / f"dump708-{L}-{h}.npz"), np.load(W / "dumps" / f"dump708-{R}-{h}.npz")
    rec["dump_max_abs"], rec["dump_bitwise"] = cmp(a, b)
    opa = json.loads((W / L / "ops" / f"held_{h}.json").read_text())["sha256"]
    rjs = W / R / "ops" / f"held_{h}.json"
    rjs = rjs if rjs.exists() else W / "ABCp" / "ops" / f"held_{h}.json"
    rec["operator_sha256"], rec["ref_operator_sha256"] = opa, json.loads(rjs.read_text())["sha256"]
    rec["operator_equal"] = opa == rec["ref_operator_sha256"]
    arms = {}
    for f in sorted((W / f"ka-{L.lower()}" / "bulk" / f"{L}__{h}").glob("*.npz")):
        g = refka / "bulk" / f"{R}__{h}" / f.name
        if g.exists():
            m, bt = cmp(np.load(f), np.load(g))
            m.pop("targets", None)
            arms[f.stem] = dict(worst=max(m.values()) if m else 0.0, bitwise=all(bt.values()))
    rec["ka804_arms"] = arms
    rec["worst_dump"] = max(rec["dump_max_abs"].values())
    rec["worst_arms"] = max((v["worst"] for v in arms.values()), default=None)
    out[h] = rec
    print(f"{L}/{h} vs {R}: operator equal {rec['operator_equal']}; dump worst |d| {rec['worst_dump']:.1e} "
          f"(bitwise {all(rec['dump_bitwise'].values())}); ka804 arms worst {rec['worst_arms']:.1e} over {len(arms)} arms")
(C / f"kacmp-{L}-vs-{R}.json").write_text(json.dumps(out, indent=1))
PY
}

cmd_profile() {                       # profile801 after the three LABC folds
  local i=0
  while :; do local n=0; for h in $FOLDS; do [[ -f $C/LABC-$h.ok ]] && n=$((n+1)); done; [[ $n -eq 3 ]] && break
    i=$((i+1)); [[ $i -lt 240 ]] || exit 9; sleep 10; done
  uv run python $E/profile801.py > $C/profile801.log 2>&1
  say "profile801: $(head -1 $C/profile801.log | cut -c1-60) ..."
}

cmd_draw_a() {                        # the A-only draws (draw801a.py); bundles for A1..A6
  wait_verify
  uv run python $E/draw801a.py > $C/draw801a.log 2>&1
  say "draw801a: $(tail -1 $C/draw801a.log)"
}

cmd_ka804() {                         # ka804 (unmodified) on one fold after its dump; standard arrays -> work/ka-aonly
  local L="$1" h="$2" i=0
  while [[ ! -f $C/$L-$h-outputs.sha256 ]]; do sleep 15; i=$((i+1)); [[ $i -lt 240 ]] || exit 4; done
  local KD="$W/ka-${KA_DIR:-aonly}"            # the kit sets KA_DIR=<tag, lower case>
  mkdir -p /work/ka804-$L-$h $KD
  KA804_OUT=/work/ka804-$L-$h uv run python experiments/804-known-answers/ka804.py "$L/$h" > $C/ka804-$L-$h.log 2>&1
  cp -R /work/ka804-$L-$h/. $KD/
  say "ka804 $L/$h: $(tail -1 $C/ka804-$L-$h.log)"
}

cmd_smoke() {                         # smoke801 on the given folds, only after all three efp known answers passed
  local tag="$1"; shift
  local i=0
  while :; do
    local n=0; for h in $FOLDS; do [[ -f $C/ka-efp-$h.ok ]] && n=$((n+1)); done
    [[ $n -eq 3 ]] && break
    i=$((i+1)); [[ $i -lt 180 ]] || { say "efp known answers did not all pass: smoke not run"; exit 9; }
    sleep 10
  done
  uv run python $E/smoke801.py "$tag" "$@" > $C/smoke-$tag.log 2>&1
  say "smoke $tag: $(grep -c "^\[" $C/smoke-$tag.log) folds priced"
}

sub="${1:?subcommand}"; shift
case "$sub" in
  verify) cmd_verify;; prep) cmd_prep "$@";; fold) cmd_fold "$@";; refit-efp) cmd_refit_efp "$@";; ka-efp) cmd_ka_efp "$@";;
  smoke) cmd_smoke "$@";; refit-study) cmd_refit_study "$@";;
  designs-L) cmd_designs_L "$@";; designs-ABC) cmd_designs_ABC "$@";;
  fit-abc) cmd_fit_abc "$@";; abcp) cmd_abcp "$@";; designs-ABCp) cmd_designs_ABCp "$@";;
  draw-a) cmd_draw_a;; ka804) cmd_ka804 "$@";; prep-labc) cmd_prep_labc;; labc) cmd_labc "$@";; regfrag) cmd_regfrag "$@";; kacmp) cmd_kacmp "$@";; labcd) cmd_labcd "$@";; profile) cmd_profile;;
  *) echo "unknown $sub"; exit 2;;
esac
