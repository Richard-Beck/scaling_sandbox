"""Update files as campaign products arrive; never refresh the browser."""
import time
from datetime import datetime,timezone
from design import *
from report import build

def main():
    previous=None
    for tick in range(120):
        paths=[p for folder in ['fits','scaling','literature','literature_validation']for p in (OUT/folder).glob('*.json')]
        retained=OUT/'retained_transfer.json'
        if retained.exists():paths.append(retained)
        state=tuple(sorted((str(p),p.stat().st_mtime_ns)for p in paths))
        if state!=previous:
            build();previous=state
        complete=(len(list((OUT/'fits').glob('*.json')))==60 and
                  all(len(list((OUT/f).glob('*.json')))==20 for f in ['scaling','literature','literature_validation']) and retained.exists())
        save_json(OUT/'producer_status.json',dict(updated=datetime.now(timezone.utc).isoformat(),
                  state='complete' if complete else 'running',poll_seconds=30))
        if complete:return
        time.sleep(30)
    save_json(OUT/'producer_status.json',dict(updated=datetime.now(timezone.utc).isoformat(),state='watcher_time_limit'))

if __name__=='__main__':main()
