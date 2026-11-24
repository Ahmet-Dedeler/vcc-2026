"""Narrow adapter to the already-validated current-family composer/scorer."""

import argparse
from datetime import datetime,timezone
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time

if __name__=="__main__":
    # The report imports this adapter; do not re-import after configuring P.
    sys.modules["run337"]=sys.modules[__name__]

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
WORK=HERE/"work"
sys.path.insert(0,str(ROOT/"experiments/327-coordinate-median"))
import run327 as P
from estimator337 import NAMES,estimate

DEPENDENCIES=P.input_paths()+[
    ROOT/"experiments/327-coordinate-median/run327.py",
    ROOT/"experiments/327-coordinate-median/estimator327.py"]
R,E,BASE=P.R,P.E,P.BASE


def paths():
    return sorted(set(DEPENDENCIES+[HERE/n for n in
                  ("PREREGISTRATION.md","estimator337.py","run337.py","report337.py","test337.py")]))


def freeze():
    if (HERE/"inputs.json").exists():
        raise FileExistsError("already frozen")
    R.write_json(HERE/"inputs.json",{"inputs_sha256":{str(p.relative_to(ROOT)):R.sha(p) for p in paths()},
                                    "arms":NAMES,"reused_H1_development":True})


def verify():
    m=json.loads((HERE/"inputs.json").read_text())
    assert set(m["inputs_sha256"])=={str(p.relative_to(ROOT)) for p in paths()}
    for name,expected in m["inputs_sha256"].items():
        assert R.sha(ROOT/name)==expected,name
    own=[str(p.relative_to(ROOT)) for p in paths() if p.parent==HERE]
    own+=[str((HERE/"inputs.json").relative_to(ROOT))]
    for name in own:
        subprocess.run(["git","ls-files","--error-unmatch",name],cwd=ROOT,check=True,
                       stdout=subprocess.DEVNULL)
    assert not subprocess.check_output(["git","status","--porcelain","--",*own],cwd=ROOT)
    E.verify()


def configure():
    P.HERE,P.WORK,P.NAMES=HERE,WORK,NAMES
    P.verify=verify


def compose():
    configure()
    details={}
    def construction(base,sources,observed,mask):
        out,eligible,audit=estimate(base,sources,observed,mask)
        details.update(audit)
        return out,eligible
    P.arms=construction
    P.compose()
    receipt=json.loads((WORK/"compose.json").read_text())
    receipt.pop("dose",None)
    receipt["availability"]=details
    receipt["construction_policy"]="PREREGISTRATION.md"
    receipt["borrowed_composer"]="experiments/327-coordinate-median/run327.py"
    R.write_json(WORK/"compose.json",receipt)


def all_steps():
    started=time.monotonic()
    start_utc=datetime.now(timezone.utc).isoformat()
    compose()
    for name in NAMES:
        if time.monotonic()-started>5400:
            raise TimeoutError("90-minute budget")
        P.evaluate(name)
    from report337 import report
    report()
    R.write_json(WORK/"execution.json",{"start_utc":start_utc,
          "finish_utc":datetime.now(timezone.utc).isoformat(),
          "seconds":time.monotonic()-started,"pid":os.getpid(),
          "max_rss_bytes_macos":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss})


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("action",choices=("freeze","verify","all"))
    args=p.parse_args()
    {"freeze":freeze,"verify":verify,"all":all_steps}[args.action]()
