"""Guarded after-any collection; never launches or reruns native protocols."""
from pathlib import Path
import hashlib
import json
import re
import subprocess
import sys
from datetime import datetime,timezone

ROOT=Path(__file__).resolve().parents[2]
RUN=ROOT/'data/runs/revision_amplitude_operating_point_20261019'
ADMIN=RUN/'terminal_collection'
OUT=ROOT/'figures/length_mechanisms_20261008/dev/revision_program_20261019/figure2/amplitude_operating_point_20261019'


def write(p,d):
    p.parent.mkdir(parents=True,exist_ok=True)
    temp=p.with_suffix('.tmp.json');temp.write_text(json.dumps(d,indent=2)+'\n');temp.replace(p)


def main():
    guard=json.loads((ADMIN/'source_guard.json').read_text())
    for relative,digest in guard['hashes'].items():
        assert hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()==digest,relative
    launch=json.loads((RUN/'launch.json').read_text());job=launch['array_job_id']
    command=['sacct','-j',job,'--noheader','--parsable2','--format=JobID,State,ExitCode']
    raw=subprocess.check_output(command,text=True)
    tasks={line.split('|')[0]:line.split('|')[1:3] for line in raw.strip().splitlines()
           if re.fullmatch(re.escape(job)+r'_\d+',line.split('|')[0])}
    assert len(tasks)==20,(len(tasks),raw)
    assert not any(state in ('PENDING','RUNNING','COMPLETING','CONFIGURING') for state,code in tasks.values()),tasks
    scheduler_pass=all(state=='COMPLETED' and code=='0:0' for state,code in tasks.values())
    write(ADMIN/'scheduler_terminal.json',dict(generated_at=datetime.now(timezone.utc).isoformat(),
        array_job_id=job,command=command,raw=raw,tasks=tasks,all_20_completed_zero=scheduler_pass))
    driver=RUN/'frozen_sources/code/mechanisms/revision_amplitude_operating_point_20261019.py'
    subprocess.run([sys.executable,str(driver),'collect','--out',str(RUN)],check=True)
    status=json.loads((RUN/'results/latest_status.json').read_text())
    status.update(terminal_array_observed=True,scheduler_all_20_completed_zero=scheduler_pass,
                  producer='finished' if status['complete'] and scheduler_pass else 'terminal_array_incomplete')
    write(RUN/'results/latest_status.json',status)
    subprocess.run([sys.executable,str(ROOT/'code/mechanisms/revision_amplitude_operating_point_readouts_20261019.py')],check=True)
    if status['complete'] and scheduler_pass:
        subprocess.run([sys.executable,str(ROOT/'code/mechanisms/revision_amplitude_operating_point_review_20261019.py')],check=True)
    write(ADMIN/'execution_status.json',dict(generated_at=datetime.now(timezone.utc).isoformat(),
        source_guard_sha256=hashlib.sha256((ADMIN/'source_guard.json').read_bytes()).hexdigest(),
        all_scientific_protocols_complete=status['complete'],scheduler_all_20_completed_zero=scheduler_pass,
        derived_review_generated=bool(status['complete'] and scheduler_pass),
        independent_terminal_and_export_review='pending'))
    if not status['complete'] or not scheduler_pass:sys.exit(2)


if __name__=='__main__':main()
