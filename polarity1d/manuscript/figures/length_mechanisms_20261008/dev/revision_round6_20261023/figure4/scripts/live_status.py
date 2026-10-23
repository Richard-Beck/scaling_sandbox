"""Producer-side current status, with no timed browser refresh."""
from common import *
import time,subprocess,re
from datetime import datetime,timezone

def update():
    counts=dict(complete=0,failed=0,running=0,initial_complete=0,additional_complete=0)
    for p in RUN.glob('workers/c*/seed*/status.json'):
        s=read(p)
        if s.get('complete'):
            counts['complete']+=1
            counts['initial_complete' if s.get('replicate',0)<24 else 'additional_complete']+=1
        elif s.get('error'):counts['failed']+=1
        else:counts['running']+=1
    p=subprocess.run(['squeue','-u',os.environ.get('USER','4473331'),'-h','-o','%j %T'],capture_output=True,text=True)
    jobs=[line for line in p.stdout.splitlines() if line.startswith('F4r6-') and not line.startswith('F4r6-live ')]
    now=datetime.now(timezone.utc).isoformat()
    ledger=read(RUN/'allocation_targets.json') if (RUN/'allocation_targets.json').exists() else None
    expected=sum(ledger['targets'].values()) if ledger else 4800
    status=dict(generated_utc=now,modes=counts,initial_expected=4800,currently_allocated_expected=expected,
      running_native_jobs=sum(line.endswith(' RUNNING') and ('F4r6-mode-' in line or 'F4r6-topup ' in line) for line in jobs),
      exemplar_complete_conditions=len(list((OUT/'exemplar_analysis').glob('c*/complete.json'))),
      exemplar_expected_conditions=260,scheduler_records=jobs,producer_active=bool(jobs),
      native_profiles_saved=False,explanation='Strong regime remains unexplained; J/K supplementary until independent mode evidence supports them.')
    write(META/'live_status.json',status)
    page=OUT/'index.html'
    if page.exists():
        s=page.read_text();s=re.sub(r'<!-- LIVE BEGIN -->.*?<!-- LIVE END -->','',s,flags=re.S)
        block=f'<!-- LIVE BEGIN --><p><strong>Producer status {now}</strong>: mode histories {counts["complete"]}/{expected} currently allocated complete ({counts["initial_complete"]}/4800 initial, {counts["additional_complete"]} precision additions), {counts["failed"]} failed; existing-profile conditions {status["exemplar_complete_conditions"]}/260. Producer {"active" if jobs else "stopped"}. <a href="supporting/metadata/live_status.json">Status record</a></p><!-- LIVE END -->'
        page.write_text(s.replace('</body>',block+'</body>'))
    return status

if __name__=='__main__':
    for iteration in range(480):
        s=update()
        if not s['producer_active']:break
        time.sleep(30)
