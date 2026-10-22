"""Manually refreshed Spring noise status reader; backend updates only."""
from pathlib import Path
from datetime import datetime,timezone
import argparse
import html
import json
import os
import subprocess
import time


def render(run,out):
    m=json.loads((run/"manifest.json").read_text());launch=json.loads((run/"submission.json").read_text())
    queue=subprocess.run(["squeue","--noheader","--jobs",launch["job_id"],"--format=%i|%T|%M|%R"],capture_output=True,text=True)
    records=[]
    for j in m["jobs"]:
        p=run/"workers"/f"{j['slot']:04d}"/"status.json"
        r=json.loads(p.read_text()) if p.exists() else dict(status="not_started",complete=False)
        records.append(dict(j,worker_status=r.get("status"),complete=r.get("complete",False),elapsed_seconds=r.get("elapsed_seconds"),producer=r.get("producer")))
    finished=sum(r["worker_status"]=="complete" for r in records)
    allterminal=all(r.get("producer")=="terminal" for r in records)
    generated=datetime.now(timezone.utc).isoformat()
    audit=out/"supporting/metadata/spring_noise_analysis_audit.json"
    summary=dict(generated_at=generated,job_id=launch["job_id"],native_workers_finished=finished,total_workers=len(records),
        native_producer_state="complete" if finished==len(records) else ("terminal_with_errors" if allterminal else "running_or_queued"),
        analysis_audit_exists=audit.exists(),scheduler_returncode=queue.returncode,scheduler_stdout=queue.stdout,scheduler_stderr=queue.stderr,workers=records)
    target=out/"supporting/metadata/spring_live_status.json";target.parent.mkdir(parents=True,exist_ok=True)
    temp=target.with_suffix(".tmp");temp.write_text(json.dumps(summary,indent=2));temp.replace(target)
    def link(p,label):return '<a href="'+html.escape(os.path.relpath(p,out),quote=True)+'">'+html.escape(label)+'</a>'
    table=""
    for C in (1.,2.,3.):
        for kind,f in (("noise",.001),("noise",.01),("mean",.001),("mean",.01),("basal",0.)):
            rows=[r for r in records if r["resting_circumference"]==C and r["kind"]==kind and r["fraction"]==f]
            count=sum(r["worker_status"]=="complete" for r in rows)
            table+=f'<tr><td>{C:g}</td><td>{kind}</td><td>{f:g}</td><td>{count}/{len(rows)}</td></tr>'
    details=""
    for r in records:
        slot=r["slot"];p=run/"workers"/f"{slot:04d}"/"status.json"
        details+=f'<tr><td>{link(p,str(slot))}</td><td>{r["resting_circumference"]:g}</td><td>{r["kind"]}</td><td>{r["fraction"]:g}</td><td>{r["seed"]}</td><td>{html.escape(r["worker_status"])}</td></tr>'
    document=f'''<!doctype html><html lang="en"><meta charset="utf-8"><title>Spring noisy-ring generation</title>
<style>body{{font:16px/1.5 system-ui,sans-serif;max-width:1050px;margin:2em auto;padding:0 1em;color:#19313d}}table{{border-collapse:collapse;margin:1em 0;width:100%}}td,th{{text-align:left;border-bottom:1px solid #dce2e6;padding:.3em .6em}}code,pre{{background:#f3f6f8;padding:.3em}}a{{color:#175c86}}.status{{padding:.8em;background:#e8f2f5;border-radius:5px}}</style>
<h1>Spring noisy-ring generation</h1><p class="status"><b>{finished}/{len(records)} native workers finished</b> · producer state: {summary['native_producer_state']}<br>Data generated {html.escape(generated)}<br>SLURM array {launch['job_id']} · analysis audit {'available' if audit.exists() else 'pending'}</p>
<p>Manually refreshing loads the latest backend-generated snapshot. This page has no browser timer. Finished simulation generation remains distinct from the subsequent independent field/statistics audit.</p>
<p>Three exact qualified basals, both declared noise strengths and 12 independent noisy seeds per group; matched deterministic mean and no-input controls. Spring uses native time/length units, uniform stretch/dilution/diffusion and the declared mean-force global-extension ring closure.</p>
<p>Finite input law: lognormal sigma 1, mean 4f with f=.001/.01; basal physical covariance width .15 transported with material coordinates; temporal knot correlation time5, final knot spacing.125. Horizon3200; scalar cadence.5, primary profiles4; maximum analysis lag800. The field is interpolated, not exact continuous-time OU noise.</p>
<p>{link(run/'manifest.json','Frozen main manifest')} · {link(run/'protocol/spring_noise_refinement_plan.json','Prospective sensitivity gates')} · {link(run/'protocol/spring_noise_bridge_audit.json','Actual paired refinement audit')} · {link(out/'figure4_stochastic_revision.png','Preserved three-model interim figure')}</p>
<h2>Coverage</h2><table><tr><th>Resting C</th><th>Condition</th><th>f</th><th>Finished / prescribed</th></tr>{table}</table>
<h2>Latest scheduler observation</h2><pre>{html.escape(queue.stdout or ('No queued/running elements reported.' if queue.returncode==0 else queue.stderr))}</pre>
<details><summary>All original prescribed workers</summary><table><tr><th>Slot / record</th><th>Resting C</th><th>Condition</th><th>f</th><th>Seed</th><th>Status</th></tr>{details}</table></details>
</html>'''
    p=out/"spring_noise_status.html";tmp=p.with_suffix(".tmp");tmp.write_text(document);tmp.replace(p)
    return summary


def main():
    p=argparse.ArgumentParser();p.add_argument("--run",type=Path,required=True);p.add_argument("--output",type=Path,required=True);p.add_argument("--watch",action="store_true");a=p.parse_args()
    while True:
        s=render(a.run.resolve(),a.output.resolve())
        if not a.watch or s["native_producer_state"] in ("complete","terminal_with_errors"):break
        time.sleep(30)
    print(json.dumps({k:v for k,v in s.items() if k!="workers"},indent=2))


if __name__=="__main__":main()
