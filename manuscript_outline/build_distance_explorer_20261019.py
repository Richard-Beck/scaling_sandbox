#!/usr/bin/env python3
"""Build an offline exploratory reader for the distance-response source tables.

The publication plots remain Matplotlib exports. This reader exposes the same
point estimates, sampled crossing brackets and missingness for discussion.
"""
from __future__ import annotations

import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / 'polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_program_20261019'
DATA = PACKAGE / 'figure3/distance_response/supporting/data'
VIEWER = PACKAGE / 'supporting/viewers/distance_response.html'

METRICS = {
    'angle_final_deg': 'Retained activity angle after withdrawal (degrees; R ≥ 0.1)',
    'angle_withdrawal_deg': 'Activity angle at cue withdrawal (degrees; R ≥ 0.1)',
    'fraction_requested_turn_final': 'Retained angle / requested cue-centre angle',
    'target_error_final_deg': 'Final signed angular error relative to cue (degrees)',
    'R_withdrawal': 'Relative first-moment magnitude R at withdrawal',
    'fixed_region_excess_withdrawal': 'Original material-region excess at withdrawal / initial excess',
    'fixed_region_excess_final': 'Final original material-region excess / initial excess',
    'fixed_region_excess_minimum': 'Minimum original material-region excess / initial excess',
    'fixed_region_excess_minimum_time': 'Saved time of minimum original material-region excess',
    'final_resolved_peaks': 'Final resolved peak count',
    'overlap_fraction_cue_support': 'Initial cue-support fraction overlapping the original core',
    'M0_ratio_withdrawal': 'Total physical active amount at withdrawal / initial amount',
    'M0_ratio_final': 'Final total physical active amount / initial amount',
    'M0_ratio_sampled_minimum': 'Minimum sampled total active amount / initial amount',
    'M0_ratio_sampled_maximum': 'Maximum sampled total active amount / initial amount',
    'R_sampled_minimum': 'Minimum sampled relative first-moment magnitude R',
    'release_signed_turn_deg': 'Signed angular change after withdrawal (degrees)',
    'fixed_region_excess_cue_on_sampled_minimum': 'Minimum original-region excess during cue / initial excess',
    'fixed_region_excess_release_sampled_minimum': 'Minimum original-region excess during release / initial excess',
    'fixed_region_release_net_change': 'Original-region excess: final minus withdrawal',
    'fixed_region_rebound_from_sampled_minimum': 'Original-region excess: final minus sampled minimum',
    'M0_initial': 'Initial total physical active amount (native units)',
    'M1_initial_magnitude': 'Initial physical first-moment magnitude (native units)',
    'q1_initial': 'Initial first shape moment / total activity',
    'M0_ratio_at_R_sampled_minimum': 'Total active amount / initial at minimum sampled R',
    'q2_at_R_sampled_minimum': 'Second shape moment / total activity at minimum sampled R',
    'old_core_excess_at_withdrawal': 'Original-region excess at withdrawal / initial excess',
    'old_core_excess_final': 'Final original-region excess / initial excess',
    'old_core_excess_minimum': 'Minimum original-region excess / initial excess',
    'old_core_excess_minimum_cue_on': 'Minimum original-region excess during cue',
    'old_core_excess_rebound_after_withdrawal': 'Original-region rebound after withdrawal',
    'old_core_excess_rebound': 'Original-region rebound after withdrawal',
    'old_core_raw_final': 'Final raw original-region amount / initial amount',
    'cue_patch_excess_final': 'Final cue-region excess / initial incumbent excess',
    'cue_fixed_material_excess_final': 'Final cue material-region excess / initial incumbent excess',
    'final_R': 'Final relative first-moment magnitude R',
    'final_angle_deg': 'Final activity angle relative to incumbent (degrees)',
    'final_cue_angle_error_deg': 'Final signed angular error relative to cue (degrees)',
    'half_loss_upper': 'First half-loss during cue: upper bound (fixed-duration references)',
    'old_core_excess_fraction': 'Original-region excess / initial excess',
    'old_core_raw_fraction': 'Raw original-region amount / initial amount',
    'cue_patch_excess_over_basal_incumbent': 'Cue-region excess / initial incumbent excess',
    'R': 'Relative first-moment magnitude R',
    'q1': 'First shape moment / total activity',
    'q2': 'Second shape moment / total activity',
    'q1_final': 'Final first shape moment / total activity',
    'q2_final': 'Final second shape moment / total activity',
    'angle_deg': 'Activity angle relative to incumbent (degrees)',
    'cue_angle_error_deg': 'Signed angular error relative to cue (degrees)',
    'upper': 'Event time: upper saved-time bound',
    'lower': 'Event time: lower saved-time bound',
    'original_core_amount_fraction': 'Raw original-region amount / initial amount',
    'resolved_peak_count': 'Resolved peak count',
    'R_at_withdrawal': 'Relative first-moment magnitude R at withdrawal',
    'R_final': 'Final relative first-moment magnitude R',
    'angle_at_withdrawal': 'Activity angle at withdrawal (degrees)',
    'angle_final': 'Final activity angle (degrees)',
    'cue_angle_error_final': 'Final signed angular error relative to cue (degrees)',
    'resolved_peaks_final': 'Final resolved peak count',
    'physical_C_min': 'Minimum physical circumference',
    'physical_C_max': 'Maximum physical circumference',
    'physical_C_final': 'Final physical circumference',
    'cue_on_loss_0p5_upper': 'First half-loss during cue: upper bound (duration branches)',
    'loss_0p5_upper': 'First half-loss time in complete observation: upper bound',
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def clean(value):
    if value is None or value == '':
        return None
    try:
        number = float(value)
        return number if math.isfinite(number) else None
    except (ValueError, TypeError):
        return value


def dataset(path, label, note, *, preferred_x=None, preferred_metric=None):
    with path.open(newline='') as stream:
        rows = [{key: clean(value) for key, value in row.items()}
                for row in csv.DictReader(stream)]
    for row in rows:
        if row.get('resting_circumference') is not None:
            row['model'] = 'spring_ring'
            row['circumference'] = row['resting_circumference']
        if 'physical_distance' not in row:
            row['physical_distance'] = row.get('physical_distance_at_onset')
            if row['physical_distance'] is None:
                row['physical_distance'] = row['circumference'] * row['fraction']
        row.setdefault('angular_separation_deg', 360 * row['fraction'])
        # Compact embedded data retain observables and their interpretation
        # fields. Full native paths/hashes stay in the linked original table.
        for key in list(row):
            if key.endswith(('_path', '_sha256')) or key in ('checkpoint_sha256', 'history_path'):
                del row[key]
    return dict(label=label, note=note, preferred_x=preferred_x, preferred_metric=preferred_metric,
                path=str(path.relative_to(ROOT)),
                href=os.path.relpath(path, VIEWER.parent), sha256=sha(path), rows=rows)


PAGE = r'''<!doctype html><html lang="en"><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Distance-response explorer</title>
<style>
body{font:16px/1.5 system-ui;max-width:1180px;margin:2rem auto;padding:0 1rem;color:#172d40}
h1{line-height:1.15}label{display:inline-flex;flex-direction:column;margin:.4rem .8rem .4rem 0;font-size:14px}
select{font:inherit;max-width:400px;padding:.4rem}#plot{width:100%;height:auto;border:1px solid #d9e3e8}
.muted{color:#536575}table{border-collapse:collapse;font-size:13px;white-space:nowrap}
td,th{text-align:left;border-bottom:1px solid #d9e3e8;padding:.35rem .6rem}.scroll{overflow:auto;max-height:380px}
#filters{background:#f3f6f8;padding:.5rem}a{color:#146d83}.point{cursor:crosshair}
</style>
<h1>How does cue distance alter polarity response?</h1>
<p class="muted">Generated __TIME__. Static snapshot of the current source tables; manually reload after regeneration.</p>
<p>Select a model, size and observable to inspect retained deterministic conditions. Chemical-ring sweeps use a 1000-unit cue;
older duration anchors and Spring sweeps have different coverage. The standardized stationary-cue experiment uses
the same four cue-centre angles in every model and retains initial-overlap cases. Time and length units are native to each model.</p>
<p>Horizontal distance and normal angle describe the imposed cue placement relative to the initial incumbent at cue onset;
they do not measure subsequent activity-site displacement. Spring uses its prepared physical circumference for onset distance,
while its size selector reports resting circumference C₀. Activity angle summarizes the whole field rather than a tracked peak position.</p>
<label>Evidence set<select id="dataset"></select></label>
<label>Model<select id="model"></select></label>
<label>Size<select id="size"></select></label>
<label>Horizontal axis<select id="x"></select></label>
<label>Observable<select id="metric"></select></label>
<div id="filters"></div><p id="note"></p><p id="counts" class="muted"></p>
<svg id="plot" viewBox="0 0 1100 540" role="img" aria-label="Selected distance response points"></svg><p id="legend"></p>
<p id="hover" aria-live="polite">Hover over a point for its condition and status.</p>
<p>Vertical intervals show saved-time crossing brackets where applicable. Open triangles mark right-censored times at their
observation limit; they are not measured event times. Missing or low-magnitude angles are omitted. Connecting lines are omitted:
the points are tested conditions. Cue duration uses a logarithmic horizontal axis. Open circles flag a failed separate full-field
tail gate; unavailable gates remain explicit in the table. A fixed original region is not a tracked site identity;
angle alone does not establish replacement.</p>
<details><summary>Selected source rows</summary><div class="scroll"><table id="table"></table></div></details>
<p><a href="../../RESULTS.md">Revision status and figure evidence</a></p>
<script id="evidence" type="application/json">__DATA__</script>
<script>
'use strict';
const D=JSON.parse(document.getElementById('evidence').textContent),M=__METRICS__;
const $=id=>document.getElementById(id),NS='http://www.w3.org/2000/svg';
const names={wave_pinning:'Wave-pinning',goryachev:'Goryachev',holmes_model3:'Holmes 3',spring_ring:'Spring ring'};
const filters=['cue_duration','time','cue_free_time','observable','threshold','window','event','magnitude_gate'];
const labels={cue_duration:'Cue duration',time:'Time from onset',cue_free_time:'Time after withdrawal',observable:'Event observable',threshold:'Threshold',window:'Observation window',event:'Event rule',magnitude_gate:'Magnitude gate'};
let E,chosen={};
const unique=a=>[...new Set(a.filter(x=>x!==null&&x!==undefined))].sort((a,b)=>typeof a==='number'&&typeof b==='number'?a-b:String(a).localeCompare(String(b)));
function options(select,values,display=x=>x,preferred=null){select.replaceChildren();for(const value of values){const o=document.createElement('option');o.value=String(value);o.textContent=display(value);select.append(o);}if(preferred!==null&&values.some(x=>String(x)===String(preferred)))select.value=String(preferred);}
function metricKeys(rows){const keys=Object.keys(M).filter(k=>rows.some(r=>typeof r[k]==='number'));for(const group of [['final_R','R_final'],['final_angle_deg','angle_final'],['final_cue_angle_error_deg','cue_angle_error_final'],['old_core_excess_rebound_after_withdrawal','old_core_excess_rebound'],['resolved_peak_count','resolved_peaks_final']]){const present=group.filter(k=>keys.includes(k));if(present.length>1){const best=present.reduce((a,b)=>rows.filter(r=>typeof r[b]==='number').length>rows.filter(r=>typeof r[a]==='number').length?b:a);for(const k of present)if(k!==best)keys.splice(keys.indexOf(k),1);}}return keys;}
function element(tag,attrs={},text=null){const e=document.createElementNS(NS,tag);for(const[k,v]of Object.entries(attrs))e.setAttribute(k,v);if(text!==null)e.textContent=text;return e;}
function init(){options($('dataset'),D.map((_,i)=>i),i=>D[i].label);$('dataset').onchange=changeDataset;$('model').onchange=changeModel;$('size').onchange=changeAxes;$('x').onchange=changeAxes;$('metric').onchange=draw;changeDataset();}
function changeDataset(){E=D[Number($('dataset').value)];options($('model'),unique(E.rows.map(r=>r.model)),x=>names[x]||x);changeModel();}
function changeModel(){let rows=E.rows.filter(r=>r.model===$('model').value);options($('size'),unique(rows.map(r=>r.circumference)),x=>`${x} ${$('model').value==='spring_ring'?'(resting C₀)':'(C)'}`,40);const axes=['physical_distance','angular_separation_deg','fraction'];if(unique(rows.map(r=>r.cue_duration)).length>1)axes.push('cue_duration');options($('x'),axes,x=>({physical_distance:'Initial physical cue distance',angular_separation_deg:'Initial cue normal angle (degrees)',fraction:'Initial cue offset / circumference',cue_duration:'Cue duration (log scale)'})[x],E.preferred_x);options($('metric'),metricKeys(rows),x=>M[x],E.preferred_metric);changeAxes();}
function changeAxes(){chosen={};buildFilters();}
function buildFilters(){let rows=E.rows.filter(r=>r.model===$('model').value&&String(r.circumference)===$('size').value);$('filters').replaceChildren();for(const key of filters){if(key===$('x').value)continue;const vals=unique(rows.map(r=>r[key]));if(vals.length<=1)continue;const lab=document.createElement('label');lab.append(labels[key]);const sel=document.createElement('select');sel.id='filter_'+key;let pref=vals[0];if(key==='observable'&&vals.includes('frozen_original_core_excess_fraction'))pref='frozen_original_core_excess_fraction';if(key==='threshold'&&vals.includes(.5))pref=.5;if(key==='window'&&vals.includes('cue_on'))pref='cue_on';if(key==='event')pref=vals.find(x=>String(x).startsWith('first'))||pref;if(key==='time'&&vals.includes(14000))pref=14000;if(key==='cue_free_time'&&vals.includes(13000))pref=13000;if(vals.some(x=>String(x)===String(chosen[key])))pref=chosen[key];options(sel,vals,x=>x,pref);chosen[key]=sel.value;sel.onchange=()=>{chosen[key]=sel.value;buildFilters();};lab.append(sel);$('filters').append(lab);rows=rows.filter(r=>r[key]==null||String(r[key])===sel.value);}draw();}
function draw(){let rows=E.rows.filter(r=>r.model===$('model').value&&String(r.circumference)===$('size').value);for(const key of filters){const sel=$('filter_'+key);if(sel)rows=rows.filter(r=>r[key]==null||String(r[key])===sel.value);}
 const metric=$('metric').value,xkey=$('x').value;const pointRows=[];let omitted=0,censored=0;for(const row of rows){let y=row[metric],lo=null,hi=null,status=row.status||'';let stem=null;if(metric==='upper'||metric==='lower')stem='';else if(metric.endsWith('_upper'))stem=metric.slice(0,-5);if(stem!==null){status=row[stem+'status']||row[stem+'_status']||status;lo=row[stem+'lower']??row[stem+'_lower'];hi=row[stem+'upper']??row[stem+'_upper'];if(stem==='half_loss_'){lo=row.half_loss_lower;hi=row.half_loss_upper;status=row.half_loss_status;}}
 const isCensored=typeof status==='string'&&status.includes('censor');if(typeof y!=='number'&&isCensored&&typeof lo==='number')y=lo;if(typeof y!=='number'||typeof row[xkey]!=='number'){omitted++;continue;}if(isCensored)censored++;pointRows.push({row,x:row[xkey],y,lo,hi,censored:isCensored});}
 $('note').textContent=E.note;const link=document.createElement('a');link.href=E.href;link.textContent='Complete source CSV';$('counts').replaceChildren(`${rows.length} selected rows; ${pointRows.length} displayed, ${omitted} unavailable, ${censored} censored. `,link);
 const svg=$('plot');svg.replaceChildren();const left=105,top=35,right=1060,bottom=450;const extent=a=>{let lo=Math.min(...a),hi=Math.max(...a);if(!a.length)return[0,1];if(lo===hi){lo-=.5;hi+=.5;}const pad=(hi-lo)*.06;return[lo-pad,hi+pad];};const logX=xkey==='cue_duration'&&pointRows.every(p=>p.x>0),transformX=x=>logX?Math.log10(x):x;const[xmin,xmax]=extent(pointRows.map(p=>transformX(p.x))),[ymin,ymax]=extent(pointRows.flatMap(p=>[p.y,...[p.lo,p.hi].filter(v=>typeof v==='number')]));const X=x=>left+(transformX(x)-xmin)/(xmax-xmin)*(right-left),Y=y=>bottom-(y-ymin)/(ymax-ymin)*(bottom-top);const line=(x1,y1,x2,y2,color='#d9e3e8')=>svg.append(element('line',{x1,y1,x2,y2,stroke:color}));
 for(let i=0;i<=5;i++){let transformed=xmin+(xmax-xmin)*i/5,x=logX?10**transformed:transformed,y=ymin+(ymax-ymin)*i/5;line(X(x),top,X(x),bottom);line(left,Y(y),right,Y(y));svg.append(element('text',{x:X(x),y:bottom+25,'text-anchor':'middle','font-size':14},x.toPrecision(3)));svg.append(element('text',{x:left-12,y:Y(y)+4,'text-anchor':'end','font-size':14},y.toPrecision(3)));}line(left,bottom,right,bottom,'#172d40');line(left,top,left,bottom,'#172d40');
 const seriesKeys=unique(pointRows.map(p=>xkey==='cue_duration'?p.row.fraction:0));const colors=['#287da8','#d47732','#238a82','#9766b4','#b04f73'];$('legend').replaceChildren();if(xkey==='cue_duration')for(let i=0;i<seriesKeys.length;i++){const span=document.createElement('span');span.style.color=colors[i%colors.length];span.textContent=`● Offset fraction ${seriesKeys[i]}   `;$('legend').append(span);}for(const p of pointRows){const color=colors[seriesKeys.indexOf(xkey==='cue_duration'?p.row.fraction:0)%colors.length];if(typeof p.lo==='number'&&typeof p.hi==='number')line(X(p.x),Y(p.lo),X(p.x),Y(p.hi),color);const marker=p.censored?element('path',{d:`M ${X(p.x)} ${Y(p.y)-7} l -6 12 h 12 Z`,fill:'white',stroke:color,'stroke-width':2}):element('circle',{cx:X(p.x),cy:Y(p.y),r:5,fill:String(p.row.full_field_tail_gate??p.row.tail_stationary).toLowerCase()==='false'?'white':color,stroke:color,'stroke-width':2});marker.setAttribute('class','point');marker.addEventListener('mouseenter',()=>{$('hover').textContent=Object.entries(p.row).filter(([k,v])=>v!==null&&!k.startsWith('trace_')).map(([k,v])=>`${k}: ${v}`).join(' · ');});svg.append(marker);}
 svg.append(element('text',{x:(left+right)/2,y:505,'text-anchor':'middle','font-size':17},$('x').selectedOptions[0]?.textContent||xkey));svg.append(element('text',{x:24,y:(top+bottom)/2,transform:`rotate(-90 24 ${(top+bottom)/2})`,'text-anchor':'middle','font-size':16},M[metric]||metric));
 const eventStem=metric.endsWith('_upper')?metric.slice(0,-5):null;const keys=unique([xkey,'cue_duration','fraction',metric,'status','half_loss_status','lower','upper','observable','threshold','window','event','magnitude_gate','source_kind','tail_stationary','field_tail_status','full_field_tail_status','initial_overlap','angle_label',...(eventStem?[eventStem+'status',eventStem+'lower',eventStem+'upper']:[]),'outcome','location_association'].filter(k=>rows.some(r=>r[k]!==undefined)));const table=$('table');table.replaceChildren();let tr=document.createElement('tr');for(const k of keys){const th=document.createElement('th');th.textContent=k;tr.append(th);}table.append(tr);for(const row of rows.slice(0,200)){tr=document.createElement('tr');for(const k of keys){const td=document.createElement('td');td.textContent=row[k]??'unavailable';tr.append(td);}table.append(tr);}}
init();
</script></html>'''


def main():
    datasets = []
    size_summary = PACKAGE / 'figure3/standard_angles_sizes/supporting/data/stationary_sizes_summary.csv'
    if size_summary.exists():
        with size_summary.open(newline='') as stream:
            size_rows = list(csv.DictReader(stream))
        complete = sum(row.get('status') == 'complete' for row in size_rows)
        datasets.append(dataset(size_summary, 'Stationary cue: three sizes × shared angles × duration',
            f'{complete}/{len(size_rows)} requested conditions analyzed in this snapshot. Chemical C=20/30/40; Spring resting C₀=1/2/3. All four cue-centre angles and five native durations are retained, including all 125 initial-overlap cases. Fixed physical support changes angular footprint with size. Compare withdrawal with the response after 13,000 cue-free units; inspect R and total activity alongside direction. Missing readouts remain unavailable; failed or unavailable full-field tail gates stay explicit. Initial magnitude, local architecture, material totals and native clocks differ across sizes/models; this is not an isolated causal geometry contrast or a refined minimum-duration search.',
            preferred_x='cue_duration', preferred_metric='angle_final_deg'))
    standard = PACKAGE / 'figure3/standard_angles/supporting/data/standard_angles_summary.csv'
    rich_standard = PACKAGE / 'figure3/stationary_characterization/supporting/data/stationary_whole_grid_summary.csv'
    if rich_standard.exists():
        standard = rich_standard
    if standard.exists():
        datasets.append(dataset(standard, 'Stationary cue: shared angles × duration (all four models)',
            'All 80 completed conditions at π, π/2, π/4 and π/8. Chemical C=40; Spring resting C₀=2. The 35 initial-overlap conditions remain included. Compare angles and R at withdrawal with those after 13,000 cue-free native units. Shared centre angles retain native cue widths, peaks and clocks; they are not matched-dose comparisons. Original material-region loss can reflect translation or redistribution. Final direction is unavailable below R=0.1; separate full-field tail flags remain visible. No refined minimum cue duration is inferred.',
            preferred_x='cue_duration', preferred_metric='angle_final_deg'))
    for name, label, note in [
        ('distance_response_summary.csv', 'Complete fixed-duration distance sweep',
         'All 652 completed chemical-ring conditions. Cue duration 1000; final readout 13,000 units after withdrawal. Half-loss times use original-region excess relative to its frozen basal background.'),
        ('distance_response_readouts.csv', 'Fixed-duration responses at selected times',
         'All completed chemical-ring conditions at the selected native time. Negative time after withdrawal denotes cue-on sampling. Angles require relative first-moment magnitude R ≥ 0.1.'),
        ('distance_crossing_brackets.csv', 'Threshold and orientation-entry times',
         'Saved-time first crossings and entries retained at every later saved sample are distinct events. Select threshold, window and event rule. No crossing at the horizon is right-censored.'),
        ('existing_duration_anchor_readouts.csv', 'Existing sparse duration anchors',
         'Nine selected chemical-ring geometries, each at durations 500/1000/2000. This is not a distance–duration surface or a minimum-duration search.'),
    ]:
        datasets.append(dataset(DATA / name, label, note))
    spring = ROOT / 'polarity1d/manuscript/data/runs/revision_spring_distance_duration_20261019/results/distance_duration_readouts.csv'
    if spring.exists():
        datasets.append(dataset(spring, 'Spring distance × duration sweep',
            'Spring native resting sizes C₀=1/2/3, qualified basal actual circumferences recorded in the table. Cue width 0.2 and peak 4; duration grid is 0.1/1/10/100/1000. Initial core overlap is excluded. Physical quadrature follows C(t); the original material region is fixed. A shortest passing tested duration is not a refined minimum.'))
    chemical_duration = DATA / 'chemical_distance_duration_with_reference.csv'
    if chemical_duration.exists():
        datasets.append(dataset(chemical_duration, 'Chemical distance × duration extension',
            'Seventy-five C40 chemical-ring conditions: 60 new branches at durations 1/10/100/300 plus 15 exact retained 1000-unit references, all read 13,000 native units after withdrawal. Source kind distinguishes new from reused evidence. Three Holmes 300-unit branches fail the separate full-field tail gate; finite values remain displayed. That gate is unavailable for the compact reused references. This tested duration grid does not identify a refined minimum.'))
    generated = datetime.now(timezone.utc).isoformat()
    embedded = json.dumps(datasets, separators=(',', ':'), allow_nan=False).replace('<', '\\u003c')
    page = PAGE.replace('__TIME__', generated).replace('__DATA__', embedded).replace('__METRICS__', json.dumps(METRICS))
    VIEWER.parent.mkdir(parents=True, exist_ok=True)
    VIEWER.write_text(page)
    record = dict(generated_utc=generated, viewer=str(VIEWER.relative_to(ROOT)),
        viewer_sha256=sha(VIEWER), source_tables=[{k:v for k,v in d.items() if k!='rows'} | {'rows':len(d['rows'])} for d in datasets],
        numerical_content='Source values copied without fitting or interpolation. Censoring and unavailable values retained.',
        scope='Exploratory local reader; Matplotlib figures remain the manuscript artwork.')
    output = PACKAGE / 'supporting/metadata/distance_explorer.json'
    output.write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
