"""Local HTML views of existing audited F4 aggregate exports, without reanalysis."""
from pathlib import Path
import hashlib
import json
import math
import pandas as pd


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


JS=r'''(() => {
const DATA=JSON.parse(document.getElementById('f4-data').textContent), $=id=>document.getElementById(id);
const colors=['#2576a5','#23856b','#c47825'];
const metricLabels={qualified_multiregion_fraction:'Contrast-qualified multiple regions',multiregion_fraction:'RAW multiple regions',low_contrast_fraction:'Low contrast',mean_R:'Mean activity-vector magnitude R',mean_A:'Mean active-amount ratio A',mean_Q:'Mean coherence ratio Q',memory:'Directional memory and eligible-pair coverage'};
const descriptions={qualified_multiregion_fraction:'Fraction of all saved times with >1 current half-contrast region AND contrast ≥10% of own basal contrast.',multiregion_fraction:'No amplitude floor: counts in nearly flat fields remain included. Read qualified occupancy and low contrast alongside this value.',low_contrast_fraction:'Fraction of all saved times with contrast <10% of the case’s initial basal contrast.',mean_R:'R=|P| relative to the initial physical first moment; R can exceed1 and is distinct from heading retention.',mean_A:'Physical active amount M₀ relative to its own basal value; this excludes inactive material.',mean_Q:'q₁/q₁(0), where q₁=|M₁|/M₀. Coherence alone does not establish one peak or unchanged shape.',memory:'Direction uses both endpoint R≥0.1. Coverage retains all available lag pairs; an unavailable direction estimate is not zero correlation.'};
const models=Object.keys(DATA.labels), grids=[.0001,.0003,.001,.003,.01,.03,.1,.3,1];
const numeric=v=>v!==null&&v!==undefined&&Number.isFinite(+v), ftext=f=>Number(f).toString(), value=(v,pct=false)=>numeric(v)?(pct?(100*v).toFixed(1)+'%':Number(v).toFixed(3)):'unavailable';
$('f4-model').innerHTML=models.map(k=>`<option value="${k}">${DATA.labels[k]}</option>`).join('');
$('f4-metric').innerHTML=Object.entries(metricLabels).map(([k,v])=>`<option value="${k}">${v}</option>`).join('');
$('f4-factor').innerHTML=grids.map(f=>`<option value="${f}">${ftext(f)}</option>`).join('');$('f4-factor').value='.03';
function updateLags(){const spring=$('f4-model').value==='spring_ring',lags=spring?[.5,4,10,50,200,800]:[5,40,100,500,2000,8000];$('f4-lag').innerHTML=lags.map(f=>`<option value="${f}">${f}</option>`).join('');$('f4-lag').value=String(lags[lags.length-1]);}
function frame(){const model=$('f4-model').value,window=$('f4-window').value;return {model,window,sizes:model==='spring_ring'?[1,2,3]:[20,40,60],lag:+$('f4-lag').value,metric:$('f4-metric').value,f:+$('f4-factor').value};}
function cond(s,size,f){return DATA.conditions.find(r=>r.model===s.model&&r.size===size&&r.fraction===f&&r.window===s.window);}
function mem(s,size,f){return DATA.memory.find(r=>r.model===s.model&&r.size===size&&r.fraction===f&&r.window===s.window&&r.lag===s.lag);}
const sx=f=>75+(Math.log10(f)+4)*230;
function section(rows,key,lo,hi,top,height,label,range=false){const sy=y=>top+height-(y-lo)/(hi-lo)*height;let out=`<text x="75" y="${top-12}" font-size="16" fill="#334956">${label}</text>`;
 for(let i=0;i<5;i++){const y=lo+(hi-lo)*i/4,py=sy(y);out+=`<line x1="75" x2="995" y1="${py}" y2="${py}" stroke="#e1e8ec"/><text x="64" y="${py+5}" text-anchor="end" font-size="14">${y.toFixed(hi<=1?2:1)}</text>`;}
 rows.forEach((rr,si)=>{let path='',connected=false;rr.forEach(r=>{const y=r? r[key]:null;if(!numeric(y)){connected=false;return;}const x=sx(r.fraction),py=sy(+y);path+=(connected?'L':'M')+x+','+py;connected=true;});out+=`<path d="${path}" fill="none" stroke="${colors[si]}" stroke-width="2.4"/>`;
 rr.forEach(r=>{if(!r||!numeric(r[key]))return;const x=sx(r.fraction),y=sy(+r[key]),n=r.available_cases??r.cases,square=n===12;
  if(range&&numeric(r[key+'_seed_min'])&&numeric(r[key+'_seed_max']))out+=`<line x1="${x}" x2="${x}" y1="${sy(r[key+'_seed_min'])}" y2="${sy(r[key+'_seed_max'])}" stroke="${colors[si]}" opacity=".5"/>`;
  const title=`f=${r.fraction}; value=${Number(r[key]).toFixed(6)}; seeds=${n}`;
  out+=square?`<rect x="${x-4}" y="${y-4}" width="8" height="8" fill="${colors[si]}" onclick="f4Pick(${r.fraction})"><title>${title}</title></rect>`:`<circle cx="${x}" cy="${y}" r="4" fill="${colors[si]}" onclick="f4Pick(${r.fraction})"><title>${title}</title></circle>`;
 });});return out;}
function render(){const s=frame(),isMemory=s.metric==='memory';$('f4-description').textContent=descriptions[s.metric];$('f4-lag-label').style.display=isMemory?'block':'none';
 const rows=s.sizes.map(size=>grids.map(f=>isMemory?mem(s,size,f):cond(s,size,f))),limit=isMemory||s.metric.endsWith('fraction')?1:Math.max(1,...rows.flat().map(r=>numeric(r[s.metric])?+r[s.metric]:0),...rows.flat().map(r=>numeric(r[s.metric+'_seed_max'])?+r[s.metric+'_seed_max']:0))*1.08;
 let body=`<svg viewBox="0 0 1060 ${isMemory?530:390}" role="img" aria-label="${DATA.labels[s.model]}: ${metricLabels[s.metric]}"><rect width="1060" height="${isMemory?530:390}" fill="white"/>`;
 body+=`<line x1="${sx(s.f)}" x2="${sx(s.f)}" y1="38" y2="${isMemory?465:312}" stroke="#abbac4" stroke-dasharray="4 4"/>`;
 if(isMemory){body+=section(rows,'direction_cosine_reported',-1,1,60,225,`Conditional direction cosine, native lag ${s.lag}`);body+=section(rows,'direction_coverage',0,1,350,115,'Eligible fraction of all lag pairs');}
 else body+=section(rows,s.metric,0,limit,55,255,metricLabels[s.metric],true);
 grids.forEach((f,i)=>{if([0,2,4,6,8].includes(i))body+=`<text x="${sx(f)}" y="${isMemory?491:337}" text-anchor="middle" font-size="15">${ftext(f)}</text>`;});
 body+=`<text x="535" y="${isMemory?520:375}" text-anchor="middle" font-size="16">Whole-signal strength factor f (logarithmic axis)</text></svg>`;
 $('f4-chart').innerHTML=body;$('f4-size-key').innerHTML=s.sizes.map((size,i)=>`<span><b style="color:${colors[i]}">━━</b> ${s.model==='spring_ring'?'C₀':'C'}=${size}</span>`).join('');
 $('f4-detail-heading').textContent=`Selected f=${ftext(s.f)} · ${DATA.labels[s.model]} · ${s.window.replace('_',' ')} window · native lag ${s.lag}`;
 $('f4-details').innerHTML=s.sizes.map((size,i)=>{const r=cond(s,size,s.f),m=mem(s,size,s.f);return `<tr><th>${s.model==='spring_ring'?'C₀':'C'}=${size}</th><td>${r.available_cases}</td><td>${value(r.multiregion_fraction,true)}</td><td>${value(r.qualified_multiregion_fraction,true)}</td><td>${value(r.low_contrast_fraction,true)}</td><td>${value(r.mean_R)}</td><td>${value(r.mean_A)}</td><td>${value(r.mean_Q)}</td><td>${value(m.direction_cosine_reported)}</td><td>${value(m.direction_coverage,true)}</td><td>${m.seeds_with_pairs}/${m.cases}</td></tr>`;}).join('');
 $('f4-window-note').textContent=s.model==='spring_ring'?'Spring: full0–3200, first0–1600, last1600–3200; scalar cadence0.5.':'Chemicals: full0–32000, first0–16000, last16000–32000; scalar cadence5.';
}
window.f4Pick=f=>{$('f4-factor').value=String(f);render();};
$('f4-model').addEventListener('change',()=>{updateLags();render();});['f4-metric','f4-factor','f4-window','f4-lag'].forEach(k=>$(k).addEventListener('change',render));
updateLags();$('f4-model').value='wave_pinning';$('f4-metric').value='qualified_multiregion_fraction';$('f4-factor').value='0.03';render();
})();'''


def render(out):
    out=Path(out);data=out/'supporting/data'
    condition=pd.read_csv(data/'condition_window_snapshot.csv');condition=condition[condition.kind=='noise']
    memory=pd.read_csv(data/'memory_all_conditions_windows_lags.csv');memory=memory[memory.kind=='noise']
    def clean(df):
        return [{k:None if isinstance(v,float) and not math.isfinite(v) else v for k,v in r.items()} for r in df.to_dict('records')]
    fields=['model','size','fraction','window','available_cases','target_cases','multiregion_fraction','qualified_multiregion_fraction','low_contrast_fraction','mean_R','mean_A','mean_Q']
    for key in list(fields[6:]):fields.extend([key+'_seed_min',key+'_seed_max'])
    mfields=['model','size','fraction','window','lag','cases','direction_cosine_reported','direction_coverage','seeds_with_pairs','seeds_with_declared_support','group_support']
    payload=dict(labels=dict(wave_pinning='Wave-pinning',goryachev='Goryachev',holmes_model3='Holmes 3',spring_ring='Spring'),conditions=clean(condition[fields]),memory=clean(memory[mfields]))
    assert len(payload['conditions'])==324 and len(payload['memory'])==1944
    viewers=out/'supporting/viewers';viewers.mkdir(exist_ok=True)
    (viewers/'figure4_explorer.js').write_text(JS)
    (data/'interactive_metric_payload.json').write_text(json.dumps(payload,allow_nan=False,separators=(',',':')))
    jsontext=json.dumps(payload,allow_nan=False,separators=(',',':')).replace('</','<\\/')
    markup='''<section id="explorer"><h2>Compare a readout across all three sizes</h2>
<p>Choose a model, readout and observation window. Each point is an existing audited condition; lines guide the eye between the nine tested factors. New levels have two seeds, while 0.001 and 0.01 have twelve. Whiskers show observed seed ranges, not confidence intervals.</p>
<div class="f4-controls"><label for="f4-model">Model<select id="f4-model"></select></label><label for="f4-metric">Readout<select id="f4-metric"></select></label><label for="f4-window">Window<select id="f4-window"><option value="full">Full</option><option value="first_half">First half</option><option value="last_half">Last half</option></select></label><label id="f4-lag-label" for="f4-lag">Memory lag<select id="f4-lag"></select></label><label for="f4-factor">Selected factor<select id="f4-factor"></select></label></div>
<p id="f4-description" class="caption"></p><p id="f4-size-key" class="f4-key"></p><div id="f4-chart" class="f4-plot"></div>
<h3 id="f4-detail-heading">Selected condition</h3><div class="scroll"><table><thead><tr><th>Size</th><th>Seeds</th><th>RAW regions</th><th>Qualified regions</th><th>Low contrast</th><th>R</th><th>A</th><th>Q</th><th>Direction cosine</th><th>Pair coverage</th><th>Seeds with pairs</th></tr></thead><tbody id="f4-details"></tbody></table></div>
<p class="caption" id="f4-window-note"></p><p class="caption">Both halves retain the midpoint and the same lag checkpoints; they are not independent replicates. Direction estimates with fewer than100 eligible pooled pairs or less than1% coverage are unavailable, not zero. Marker shapes: squares for twelve seeds, circles for two. Select a factor or click a point for exact context. The three sizes share this readout axis; native time/forcing are not calibrated across models.</p>
<noscript><p>The complete static overview and the memory/coverage table below remain available without JavaScript.</p></noscript></section>'''
    markup+=f'<script type="application/json" id="f4-data">{jsontext}</script><script src="supporting/viewers/figure4_explorer.js"></script>'
    return markup,dict(condition_records=324,memory_records=1944,source_condition_sha256=sha(data/'condition_window_snapshot.csv'),source_memory_sha256=sha(data/'memory_all_conditions_windows_lags.csv'),payload_sha256=sha(data/'interactive_metric_payload.json'),javascript_sha256=sha(viewers/'figure4_explorer.js'),helper_source_sha256=sha(__file__))


CSS='''.f4-controls{display:flex;flex-wrap:wrap;gap:1rem;background:#f1f6f7;padding:1rem}.f4-controls label{font-size:.85rem;color:#4c626d}.f4-controls select{display:block;font:inherit;color:#233c48;background:white;border:1px solid #b5c8cf;border-radius:4px;padding:.45rem;max-width:330px}.f4-key{display:flex;gap:1.5rem;flex-wrap:wrap}.f4-plot{overflow-x:auto}.f4-plot svg{width:100%;min-width:720px}.f4-plot circle,.f4-plot rect[onclick]{cursor:pointer}.equation{position:relative;margin:1.2rem 0;padding-right:2rem;text-align:center}.equation math{font-size:1.18rem}.eqno{position:absolute;right:0;top:35%}'''


MATH='''<p>The same physical-activity vector used in Figure3 weights every cell by its physical arc length and outward-normal direction. For the circular models:</p>
<div class="equation"><math display="block"><mrow><msub><mi>M</mi><mn>0</mn></msub><mo>=</mo><munder><mo>∑</mo><mi>j</mi></munder><msub><mi>u</mi><mi>j</mi></msub><mi>Δ</mi><msub><mi>s</mi><mi>j</mi></msub><mo>,</mo><mspace width="1em"/><msub><mi>M</mi><mn>1</mn></msub><mo>=</mo><munder><mo>∑</mo><mi>j</mi></munder><msub><mi>u</mi><mi>j</mi></msub><msup><mi>e</mi><mrow><mi>i</mi><msub><mi>θ</mi><mi>j</mi></msub></mrow></msup><mi>Δ</mi><msub><mi>s</mi><mi>j</mi></msub></mrow></math><span class="eqno">(1)</span></div>
<div class="equation"><math display="block"><mrow><mi>P</mi><mo>(</mo><mi>t</mi><mo>)</mo><mo>=</mo><mfrac><mrow><msub><mi>M</mi><mn>1</mn></msub><mo>(</mo><mi>t</mi><mo>)</mo></mrow><mrow><msub><mi>M</mi><mn>1</mn></msub><mo>(</mo><mn>0</mn><mo>)</mo></mrow></mfrac><mo>,</mo><mspace width="1em"/><msub><mi>q</mi><mn>1</mn></msub><mo>=</mo><mfrac><mrow><mo>|</mo><msub><mi>M</mi><mn>1</mn></msub><mo>|</mo></mrow><msub><mi>M</mi><mn>0</mn></msub></mfrac><mo>,</mo><mspace width="1em"/><mi>R</mi><mo>=</mo><mo>|</mo><mi>P</mi><mo>|</mo><mo>=</mo><mfrac><mrow><msub><mi>M</mi><mn>0</mn></msub><mo>(</mo><mi>t</mi><mo>)</mo></mrow><mrow><msub><mi>M</mi><mn>0</mn></msub><mo>(</mo><mn>0</mn><mo>)</mo></mrow></mfrac><mfrac><mrow><msub><mi>q</mi><mn>1</mn></msub><mo>(</mo><mi>t</mi><mo>)</mo></mrow><mrow><msub><mi>q</mi><mn>1</mn></msub><mo>(</mo><mn>0</mn><mo>)</mo></mrow></mfrac><mo>=</mo><mi>A</mi><mi>Q</mi></mrow></math><span class="eqno">(2)</span></div>
<p class="caption">P divides by the initial complex M₁: this rotates the initial axis and scales its magnitude. It does not divide by the current amount. Mean R is the mean of paired framewise A×Q, not the product of separate population means. Actual Spring C(t) enters Δs; angles use its material positions. A large R does not imply one unchanged peak, and a weak R does not establish low total activity or uniformity.</p>'''
