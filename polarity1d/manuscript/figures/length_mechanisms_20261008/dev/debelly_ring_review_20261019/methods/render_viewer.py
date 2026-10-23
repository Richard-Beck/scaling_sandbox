#!/usr/bin/env python3
"""Regenerate the retained viewer using only files in this review package."""
from datetime import datetime,timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
sys.dont_write_bytecode=True
import numpy as np

PACKAGE=Path(__file__).resolve().parents[1]
OUT=PACKAGE/'viewer'
FIELDS=[dict(key='R',label='Active Rac',unit='concentration',split=False),
        dict(key='Q',label='Active Rho',unit='concentration',split=False),
        dict(key='b',label='Bound MCA',unit='concentration',split=True),
        dict(key='c',label='Unbound MCA',unit='concentration',split=False),
        dict(key='u',label='Membrane displacement',unit='µm',split=True),
        dict(key='v',label='Cortical velocity',unit='µm/s',split=True),
        dict(key='T',label='Membrane tension',unit='model units',split=True)]

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def compact(values,digits=10):return [float(f'{float(v):.{digits}g}') for v in values]
def cuts_for(cuts):return [dict(c,face=int(face),position=float(c['x'])) for face,c in sorted(cuts.items(),key=lambda x:int(x[0]))]
def sample(archive,row,metadata):
    start,end=map(int,archive['offsets'][row:row+2]);h=archive['widths'][start:end]
    return dict(time=float(archive['time'][row]),stage=metadata['stage'],next_id=metadata['next_id'],
        edges=compact(np.r_[0.,np.cumsum(h)],12),
        fields={q['key']:compact(archive[q['key']][start:end]) for q in FIELDS},cuts=cuts_for(metadata['cuts']))

def payload(key,segment,samples,**kwargs):
    ranges={q['key']:[min(min(s['fields'][q['key']]) for s in samples),max(max(s['fields'][q['key']]) for s in samples)] for q in FIELDS}
    value=dict(case_key=key,segment_key=segment,samples=samples,ranges=ranges,reference_index=0,
        pre_index=None,post_index=None,after_step_index=None,event_time=None);value.update(kwargs)
    folder=OUT/'data';folder.mkdir(exist_ok=True)
    path=folder/f'{key}__{segment}.js'
    text=json.dumps(value,separators=(',',':'),ensure_ascii=False,allow_nan=False).replace('<','\\u003c')
    path.write_text('window.DEBELLY_EVENT_DATA=window.DEBELLY_EVENT_DATA||{};\n'
        +'window.DEBELLY_EVENT_DATA['+json.dumps(key+'__'+segment)+']='+text+';\n')
    return dict(key=segment,src='data/'+path.name),path

def render():
    raw_manifest=read(OUT/'raw/manifest.json')
    for name,digest in raw_manifest['files'].items():
        if sha(OUT/name)!=digest:raise RuntimeError('Retained local input changed: '+name)
    overview=read(OUT/'raw/overview.json')
    if not overview['validated_replays']:raise RuntimeError('Retained replay validation is missing')
    spec=importlib.util.spec_from_file_location('retained_viewer_ui',PACKAGE/'methods/viewer_ui.py')
    ui=importlib.util.module_from_spec(spec);spec.loader.exec_module(ui)
    cases=[];outputs={};checks=[];default_segment=None
    for branch in overview['branches']:
        key=branch['key'];folder=OUT/branch['folder'];meta=read(folder/'frames.json');events=read(folder/'events.json')
        with np.load(folder/'high_resolution.npz',allow_pickle=False) as z:archive={q:z[q] for q in z.files}
        frames=meta['frames']
        with np.load(folder/'context.npz',allow_pickle=False) as z:context={q:z[q] for q in z.files}
        context_meta=read(folder/'context_frames.json')['frames']
        context_samples=[sample(context,i,row) for i,row in enumerate(context_meta)]
        item,path=payload(key,'overview',context_samples);item.update(label='Whole trajectory · coarse 250-s context',event=None)
        segments=[item];outputs[str(path.relative_to(OUT))]=sha(path)
        for number,t in enumerate(sorted(set(float(e['time']) for e in events))):
            selected=[i for i,f in enumerate(frames) if max(0.,t-30.)<=f['time']<=t+30.]
            samples=[sample(archive,i,frames[i]) for i in selected]
            def where(time,stage):return [j for j,i in enumerate(selected) if frames[i]['time']==time and frames[i]['stage_code']==stage]
            pre,post,after=where(t,0),where(t,1),where(t+.5,0)
            if not(len(pre)==len(post)==len(after)==1):raise RuntimeError('Incomplete retained event stages')
            ordinary=[f['time'] for f in (frames[i] for i in selected) if f['stage_code']==0]
            if not np.all(np.diff(ordinary)==.5):raise RuntimeError('Incomplete native half-second sequence')
            grouped=[e for e in events if e['time']==t]
            description='; '.join(('heal' if e['event']=='close' else 'open')+f" gap {e['id']} at {float(e['x']):.4f} µm" for e in grouped)
            segment_key=f'event_{number:02d}'
            item,path=payload(key,segment_key,samples,reference_index=pre[0],pre_index=pre[0],post_index=post[0],
                after_step_index=after[0],event_time=t,event_description=description)
            item.update(label=f't = {t:g} s · {description}',event=dict(time=t,description=description))
            segments.append(item);outputs[str(path.relative_to(OUT))]=sha(path)
            checks.append(dict(branch=key,event=t,frames=len(samples),native_dt=.5,pre_index=pre[0],post_index=post[0],after_step_index=after[0],all_seven_fields=True))
            if key==overview['default_branch'] and t>1000. and default_segment is None:default_segment=segment_key
        cases.append(dict(key=key,label=branch['label'],cue_location=branch['cue_location'],circumference=50.,segments=segments,
            raw_data=str((folder/'high_resolution.npz').relative_to(OUT)),raw_metadata=str((folder/'frames.json').relative_to(OUT))))
    stamp=datetime.now(timezone.utc).isoformat(timespec='seconds')
    config=dict(generated_utc=stamp,status='retained review; original replay and independent validation passed',
        fields=FIELDS,cases=cases,default_case=overview['default_branch'],default_segment=default_segment,
        display_significant_digits=10,native_dt=.5,event_window_half_width=30.,
        validation_href='../methods/provenance/original_replay_validation.json',report_manifest_href='report_manifest.json')
    html=ui.HTML_TEMPLATE.replace('@@CONFIG@@',json.dumps(config,separators=(',',':'),ensure_ascii=False).replace('<','\\u003c'))
    if '@@CONFIG@@' in html:raise RuntimeError('Viewer configuration was not inserted')
    (OUT/'index.html').write_text(html);outputs['index.html']=sha(OUT/'index.html')
    manifest=dict(generated_utc=stamp,producer_state='retained_review_complete',
        local_inputs={**raw_manifest['files'],'raw/manifest.json':sha(OUT/'raw/manifest.json')},
        sources={str(p.relative_to(PACKAGE)):sha(p) for p in [Path(__file__),PACKAGE/'methods/viewer_ui.py']},
        outputs=outputs,event_export_checks=checks,display_significant_digits=10,
        full_precision_event_and_context_archives_retained=True,no_interpolation=True,no_timed_browser_refresh=True,
        relocatable=True,external_runtime_dependencies=False)
    (OUT/'report_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(dict(viewer='viewer/index.html',branches=len(cases),events=len(checks),payloads=len(outputs)-1,
        logical_payload_bytes=sum((OUT/p).stat().st_size for p in outputs if p.endswith('.js')))))

if __name__=='__main__':render()
