#!/usr/bin/env python3
"""Extract per-job wall times from a completed Snakemake text log."""
import argparse, re
from datetime import datetime
from pathlib import Path
import pandas as pd

p=argparse.ArgumentParser(); p.add_argument("log"); p.add_argument("output"); a=p.parse_args()
stamp=re.compile(r"^\[([A-Z][a-z]{2} [A-Z][a-z]{2} \d+ \d\d:\d\d:\d\d \d{4})\]$")
rule=re.compile(r"^(?:local)?rule ([^:]+):$"); done=re.compile(r"Finished jobid: (\d+) \(Rule: ([^)]+)\)")
current_time=None; starts={}; pending_rule=None; rows=[]
for line in Path(a.log).read_text(errors="replace").splitlines():
    m=stamp.match(line)
    if m: current_time=datetime.strptime(m.group(1),"%a %b %d %H:%M:%S %Y"); continue
    m=rule.match(line)
    if m: pending_rule=m.group(1); continue
    if pending_rule and line.strip().startswith("jobid:"):
        jid=line.split(":",1)[1].strip(); starts[jid]=(pending_rule,current_time); pending_rule=None; continue
    m=done.search(line)
    if m and m.group(1) in starts:
        name,start=starts.pop(m.group(1)); rows.append({"jobid":m.group(1),"rule":name,"start":start.isoformat(),"finish":current_time.isoformat(),"wall_seconds":(current_time-start).total_seconds()})
out=Path(a.output); out.parent.mkdir(parents=True,exist_ok=True); pd.DataFrame(rows).to_csv(out,sep="\t",index=False)
