"""Keep the producer-status report current without browser self-refresh."""
import argparse
from datetime import datetime, timezone
import html
import json
import os
from pathlib import Path
import subprocess
import time

ROOT=Path(__file__).resolve().parents[1]
DEST=ROOT/"manuscript_outline/revision_round2_20261022"
RUN=ROOT/"polarity1d/manuscript/data/runs/signature_holmes_retention_round2_20261022"


def update():
    observed=datetime.now(timezone.utc).isoformat()
    result=subprocess.run(["squeue","-j","21279666","-h","-o","%i %T %M"],text=True,capture_output=True,timeout=20)
    if result.returncode:
        raise RuntimeError(result.stderr)
    records=result.stdout.strip().splitlines()
    completed=len(list((RUN/"workers").glob("*/result.json")))
    phases={p.name:len(list((p/"trial").glob("phase_*.npz"))) for p in (RUN/"workers").iterdir() if p.is_dir()}
    producer="running" if records else "stopped"
    state=dict(generated_utc=observed,figure2=dict(completed_workers=completed,total_workers=32,
        producer_state=producer,live_scheduler_records=records,saved_phase_archives=phases),
        figure3="complete; all nine approved extensions independently audited",
        figure4="complete; all native scalars and bootstrap draws independently checked",
        final_integration="waits for complete Figure2 ensemble and final audit")
    temporary=DEST/"completion_progress.json.tmp"
    temporary.write_text(json.dumps(state,indent=2)+"\n");temporary.replace(DEST/"completion_progress.json")
    status=(f"Figure 2: {completed}/32 native histories complete; {len(records)} producer tasks {producer}. "
        "Figures 3 and 4 are complete. Final integration follows the complete Figure 2 ensemble audit.")
    body=["<!doctype html><meta charset='utf-8'><title>Round-2 completion progress</title>",
        "<style>body{max-width:1400px;margin:2rem auto;font:17px/1.5 system-ui}img{width:100%}</style>",
        "<h1>Figures 2–4 completion progress</h1>",
        "<p>Data generated "+observed+". Native producer: "+producer+".</p>",
        "<p>"+html.escape(status)+"</p>",
        "<p>Figure 2 is a preview of independently checked deterministic measurements until its full noise audit is complete. The page reads the latest saved report when opened or manually refreshed.</p>"]
    for number,stem in [(2,"figure2_round2_preview"),(3,"figure3_main"),(4,"figure4_main")]:
        path=ROOT/f"polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round2_20261022/figure{number}/{stem}.png"
        body.extend([f"<h2>Figure {number}</h2>",'<img alt="Figure '+str(number)+'" src="'+os.path.relpath(path,DEST)+'">'])
    temporary=DEST/"completion_progress.html.tmp"
    temporary.write_text("\n".join(body));temporary.replace(DEST/"completion_progress.html")
    print(json.dumps(dict(generated_utc=observed,completed_workers=completed,live_tasks=len(records),producer_state=producer)),flush=True)
    return bool(records)


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--watch",action="store_true");args=parser.parse_args()
    while True:
        active=update()
        if not args.watch or not active:
            break
        time.sleep(30)
