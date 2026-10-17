"""Measured field timings and storage; no inference or downloads."""
import json,statistics
from common import BASE,profile_fields,output_dir,read_json,write_json,csv_write
rows=[]
for key in profile_fields('smoke'):
    out=output_dir('smoke',key);c=read_json(out/'completed.json');s=read_json(out/'segmentation.json');a=read_json(out/'analysis.json')
    pairs={f:0 for f in s['frames']}
    for t in s['timings']:pairs[t['frame']]+=t['seconds']
    rows.append(dict(key=key,inference_job=s['job_id'],latest_pipeline_job=c['job_id'],latest_pipeline_reused_masks=s['job_id']!=c['job_id'],frames=len(pairs),detections=a['detections'],assigned_cells=a['assigned_cells'],assigned_nuclei=a['assigned_nuclei'],median_pair_seconds=statistics.median(pairs.values()),mean_pair_seconds=statistics.mean(pairs.values()),prepare_seconds=c['prepare_seconds'],analysis_seconds=a['analysis_qc_seconds'],wall_seconds=c['total_seconds'],peak_rss_MiB=c['max_rss_KiB']/1024,masks_bytes=sum((out/f).stat().st_size for f in ['cells.tif','nuclei.tif']),isolated_runs=a['isolated_runs'],longest_run=a['longest_eligible_run']))
csv_write(BASE/'results/smoke/benchmark.csv',rows)
print(json.dumps(rows,indent=2))
