"""Identity-bound local embedding generation; never load biological response labels."""

import argparse
from datetime import datetime, timezone
import hashlib
from importlib.metadata import version
import json
import os
from pathlib import Path
import re
import resource
import subprocess
import time
import zipfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE/"work"
MODEL = WORK/"model"
REV = "b79526d6ef3645e0df4530322e266f24c829f5ef"
REPO = "NeuML/pubmedbert-base-embeddings"
FILES = ("README.md", "config.json", "model.safetensors", "tokenizer.json",
         "tokenizer_config.json", "special_tokens_map.json", "vocab.txt",
         "added_tokens.json", "1_Pooling/config.json", "sentence_bert_config.json")
WEIGHT_SHA = "929ddc16369bb4ff6f8e92a6ee9f5e0748f2bef33989b5617533dd50beebe4ef"
PREVIOUS = ROOT/"experiments/334-text-conditioned-transfer"
ZIP = ROOT/"experiments/333-genept-feature-feasibility/work/GenePT_emebdding_v2.zip"
TEXT_SHAS = {"NCBI_summary_of_genes.json":
             "3db721b7cfbc35795428c6a7ab4e8a9c59fb5f5ccb8c64af50dd3c97c2d4de2e",
             "NCBI_UniProt_summary_of_genes.json":
             "319897d7e2f8a0306bccdf2d022982b3c8321dbec8273739094021d57e414c31"}


def sha(p):
    h = hashlib.sha256()
    with p.open("rb") as f:
        for block in iter(lambda:f.read(8<<20),b""):
            h.update(block)
    return h.hexdigest()


def write(p,obj):
    p.parent.mkdir(parents=True,exist_ok=True)
    q = p.with_suffix(p.suffix+".tmp")
    q.write_text(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+"\n")
    q.replace(p)


def own_inputs():
    return [HERE/n for n in ("PREREGISTRATION.md","requirements.in",
                            "requirements.lock","run335.py","test335.py")]


def paths():
    return own_inputs()+[ROOT/"uv.lock",ZIP,PREVIOUS/"selection.json",
                         PREVIOUS/"recovery-inputs.json",
                         PREVIOUS/"recovery-work/features.json",
                         ROOT/"experiments/331-tabicl-runtime/requirements.lock"]


def freeze():
    if (HERE/"inputs.json").exists():
        raise FileExistsError("already frozen")
    write(HERE/"inputs.json",{"sha256":{str(p.relative_to(ROOT)):sha(p) for p in paths()},
                             "repository":REPO,"revision":REV,"weights_sha256":WEIGHT_SHA})


def verify():
    m = json.loads((HERE/"inputs.json").read_text())
    assert set(m["sha256"]) == {str(p.relative_to(ROOT)) for p in paths()}
    for name,expected in m["sha256"].items():
        assert sha(ROOT/name) == expected,name
    own = [str(p.relative_to(ROOT)) for p in own_inputs()+[HERE/"inputs.json"]]
    for name in own:
        subprocess.run(["git","ls-files","--error-unmatch",name],cwd=ROOT,check=True,
                       stdout=subprocess.DEVNULL)
    assert not subprocess.check_output(["git","status","--porcelain","--",*own],cwd=ROOT)


def download():
    verify()
    if (WORK/"download.json").exists() or MODEL.exists():
        raise FileExistsError("do not overwrite downloads")
    url = f"https://huggingface.co/api/models/{REPO}/tree/{REV}?recursive=true&expand=true"
    raw = subprocess.check_output(["curl","--fail","--silent","--show-error","--max-time","30",url])
    tree = json.loads(raw)
    index = {r["path"]:r for r in tree}
    assert index["model.safetensors"]["lfs"]["oid"] == WEIGHT_SHA
    assert index["model.safetensors"]["size"] == 437951328
    assert sum(index[f]["size"] for f in FILES) < 500_000_000
    write(WORK/"official-tree.json",tree)
    receipts = {}
    for filename in FILES:
        p = MODEL/filename
        p.parent.mkdir(parents=True,exist_ok=True)
        subprocess.run(["curl","--fail","--silent","--show-error","--location",
                        "--max-time","900","--output",str(p),
                        f"https://huggingface.co/{REPO}/resolve/{REV}/{filename}"],check=True)
        item = index[filename]
        assert p.stat().st_size == item["size"],filename
        actual = sha(p)
        if item.get("lfs"):
            assert actual == item["lfs"]["oid"],filename
        else:
            data = p.read_bytes()
            blob = hashlib.sha1(f"blob {len(data)}\0".encode()+data).hexdigest()
            assert blob == item["oid"],filename
        receipts[filename] = {"sha256":actual,"bytes":p.stat().st_size,"git_oid":item["oid"]}
    assert "license: apache-2.0" in (MODEL/"README.md").read_text()
    config = json.loads((MODEL/"config.json").read_text())
    pool = json.loads((MODEL/"1_Pooling/config.json").read_text())
    assert config["model_type"] == "bert" and config["hidden_size"] == 768
    assert pool["pooling_mode_mean_tokens"] and not pool["pooling_mode_cls_token"]
    write(WORK/"download.json",{"repository":REPO,"revision":REV,"files":receipts,
                               "utc":datetime.now(timezone.utc).isoformat(),
                               "license":"Apache-2.0","paid_api":False})
    verify()


def body(mapping,key):
    text = mapping.get(key,"")
    if not isinstance(text,str):
        raise ValueError("non-string description")
    m = re.match(r"^Gene Symbol "+re.escape(key)+r"(?:\s+|$)",text)
    return text[m.end():].strip() if m else ""


def documents():
    s = json.loads((PREVIOUS/"selection.json").read_text())
    genes = s["targets"]+s["genes"]
    previous = json.loads((PREVIOUS/"recovery-work/features.json").read_text())["0"]
    dictionaries = []
    with zipfile.ZipFile(ZIP) as z:
        for suffix,expected in TEXT_SHAS.items():
            member, = [n for n in z.namelist() if Path(n).name == suffix]
            raw = z.read(member)
            assert hashlib.sha256(raw).hexdigest() == expected
            dictionaries.append(json.loads(raw))
    texts = [body(dictionaries[0],g) or body(dictionaries[1],g) for g in genes]
    for g,t in zip(genes,texts):
        assert hashlib.sha256(t.encode()).hexdigest() == previous["description_sha256"][g],g
    return genes,texts


def mean_pool(hidden,mask):
    weights = mask.unsqueeze(-1).to(hidden.dtype)
    return (hidden*weights).sum(1)/weights.sum(1).clamp(min=1)


def encode():
    verify()
    if (WORK/"features.npz").exists():
        raise FileExistsError("feature bank already exists")
    download = json.loads((WORK/"download.json").read_text())
    for filename,receipt in download["files"].items():
        assert sha(MODEL/filename) == receipt["sha256"]
    # Network disabled at the model-loader layer; no remote code, fallback or telemetry.
    os.environ.update(HF_HUB_OFFLINE="1",TRANSFORMERS_OFFLINE="1",
                      HF_HUB_DISABLE_TELEMETRY="1",TOKENIZERS_PARALLELISM="false")
    import numpy as np
    import torch
    from transformers import AutoModel,AutoTokenizer
    torch.set_num_threads(6)
    torch.set_num_interop_threads(1)
    torch.use_deterministic_algorithms(True)
    start = time.monotonic()
    started = datetime.now(timezone.utc).isoformat()
    tokenizer = AutoTokenizer.from_pretrained(MODEL,local_files_only=True,
                                             trust_remote_code=False)
    model = AutoModel.from_pretrained(MODEL,local_files_only=True,
                                     trust_remote_code=False,use_safetensors=True)
    model.to(device="cpu",dtype=torch.float32).eval()
    genes,texts = documents()
    available = np.array([bool(t) for t in texts])
    assert available[:1024].sum() >= .95*1024 and available[1024:].sum() >= .9*512
    values = np.zeros((len(genes),768),dtype=np.float32)
    active = np.flatnonzero(available)
    token_hashes,full_lengths = {},{}
    for i in active:
        ids = tokenizer(texts[i],truncation=False)["input_ids"]
        full_lengths[genes[i]] = len(ids)
    def batch(rows):
        inputs = tokenizer([texts[i] for i in rows],padding=True,truncation=True,
                           max_length=512,return_tensors="pt")
        for row,ids,mask in zip(rows,inputs["input_ids"],inputs["attention_mask"]):
            consumed = ids[mask.bool()].tolist()
            token_hashes[genes[row]] = hashlib.sha256(json.dumps(consumed).encode()).hexdigest()
        with torch.inference_mode():
            pooled = mean_pool(model(**inputs).last_hidden_state,inputs["attention_mask"])
            pooled = torch.nn.functional.normalize(pooled,p=2,dim=1)
        return pooled.numpy()
    for begin in range(0,len(active),8):
        if time.monotonic()-start > 2700:
            raise TimeoutError("45-minute encoding budget exceeded")
        rows = active[begin:begin+8]
        values[rows] = batch(rows)
        if begin % 256 == 0:
            print(f"encoded {min(begin+8,len(active))}/{len(active)}",flush=True)
    repeated = batch(active[:8])
    repeat_error = float(np.max(np.abs(repeated-values[active[:8]])))
    assert repeat_error <= 1e-6,repeat_error
    assert np.isfinite(values).all() and np.all(np.linalg.norm(values[available],axis=1)>0)
    assert not np.any(values[~available])
    groups = {}
    for i in active:
        groups.setdefault(hashlib.sha256(values[i].tobytes()).hexdigest(),[]).append(int(i))
    duplicates = [rows for rows in groups.values() if len(rows)>1]
    for rows in duplicates:
        assert len({token_hashes[genes[i]] for i in rows}) == 1,"different-input vector collision"
    np.savez_compressed(WORK/"features.npz",genes=np.array(genes),matrix=values,available=available)
    verify()
    write(WORK/"results.json",{
        "status":"ENGINEERING_FEASIBLE_ONLY","repository":REPO,"revision":REV,
        "rows":len(genes),"dimensions":768,"available_targets":int(available[:1024].sum()),
        "available_genes":int(available[1024:].sum()),"repeat_max_abs_difference":repeat_error,
        "exact_vector_duplicate_groups":[[genes[i] for i in rows] for rows in duplicates],
        "different_token_vector_collisions":0,
        "truncated_rows":sum(n>512 for n in full_lengths.values()),
        "full_token_lengths":full_lengths,"consumed_token_sha256":token_hashes,
        "description_sha256":{g:hashlib.sha256(t.encode()).hexdigest() for g,t in zip(genes,texts)},
        "start_utc":started,"finish_utc":datetime.now(timezone.utc).isoformat(),
        "seconds":time.monotonic()-start,
        "max_rss_bytes_macos":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        "runtime":{p:version(p) for p in ("numpy","torch","transformers","safetensors","tokenizers")},
        "features_sha256":sha(WORK/"features.npz"),"inputs_sha256":sha(HERE/"inputs.json"),
        "download_sha256":sha(WORK/"download.json"),
        "biological_response_read":False,"paid_api":False,"submission":False})
    print("ENGINEERING_FEASIBLE_ONLY",flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("action",choices=("freeze","verify","download","encode"))
    globals()[parser.parse_args().action]()
