#!/usr/bin/env python3
"""Build an offline SVG contour viewer; Python standard library only."""
from pathlib import Path
from collections import defaultdict, Counter
import json
import math
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'reports/lfct_20x_contour_viewer'
PALETTE = ['#56B4E9', '#E69F00', '#009E73', '#CC79A7', '#F0E442', '#D55E00', '#8172B3', '#80CDC1']
CLIP_LENGTH = 50
LABELS = {'hek293': 'HEK293', 'mcf10a': 'MCF10A', 'mdamb': 'MDA-MB-231', 'u87': 'U87'}


def build_set(key, label):
    path = min((ROOT / 'data/lfct_trackmate_20x' / key).glob('*.xml'), key=lambda p: int(p.stem.split('_')[-1]))
    doc = ET.parse(path).getroot()
    meta = doc.find('Settings/ImageData').attrib
    width, height, nframes = (int(meta[k]) for k in ('width', 'height', 'nframes'))
    assert nframes >= CLIP_LENGTH
    start = (nframes - CLIP_LENGTH) // 2
    crop = [0, 0, width, height]
    left, top, cw, ch = crop
    mapping = {}
    for track in doc.findall('Model/AllTracks/Track'):
        tid = track.get('TRACK_ID')
        for edge in track.findall('Edge'):
            for attr in ('SPOT_SOURCE_ID', 'SPOT_TARGET_ID'):
                sid = edge.get(attr)
                assert sid not in mapping or mapping[sid] == tid
                mapping[sid] = tid
    frames = [[] for _ in range(CLIP_LENGTH)]
    missing = Counter()
    for block in doc.findall('Model/AllSpots/SpotsInFrame'):
        f = int(block.get('frame'))
        if not start <= f < start + CLIP_LENGTH:
            continue
        for spot in block.findall('Spot'):
            a = spot.attrib
            x = float(a['POSITION_X']) / float(meta['pixelwidth'])
            y = float(a['POSITION_Y']) / float(meta['pixelheight'])
            vals = [float(v) for v in (spot.text or '').split()]
            if len(vals) < 6 or len(vals) % 2:
                if left <= x <= left + cw and top <= y <= top + ch:
                    missing[f] += 1
                points = []
                box = (x, y, x, y)
            else:
                assert all(math.isfinite(v) for v in vals)
                assert len(vals) // 2 == int(a.get('ROI_N_POINTS', len(vals) // 2))
                points = [(x + vals[i] / float(meta['pixelwidth']), y + vals[i+1] / float(meta['pixelheight'])) for i in range(0, len(vals), 2)]
                xs, ys = zip(*points)
                box = (min(xs), min(ys), max(xs), max(ys))
            if box[2] < left or box[0] > left+cw or box[3] < top or box[1] > top+ch:
                continue
            frames[f-start].append(dict(id=a['ID'], track=mapping.get(a['ID']), points=points, box=box, x=x, y=y))
    # Weighted adjacency across the entire clip: close bounding boxes plus three
    # nearest centers per detection. Same TrackMate ID (including branches) stays one color.
    graph = defaultdict(Counter)
    for cells in frames:
        for c in cells:
            if c['track'] is not None:
                graph[c['track']]
        nearest_sets = [set(sorted((j for j in range(len(cells)) if j != i), key=lambda j: (a['x']-cells[j]['x'])**2+(a['y']-cells[j]['y'])**2)[:3]) for i, a in enumerate(cells)]
        for i, a in enumerate(cells):
            for j in range(i+1, len(cells)):
                b = cells[j]
                if a['track'] is None or b['track'] is None or a['track'] == b['track']:
                    continue
                aa, bb = a['box'], b['box']
                gap = math.hypot(max(0, aa[0]-bb[2], bb[0]-aa[2]), max(0, aa[1]-bb[3], bb[1]-aa[3]))
                if gap <= 12 or j in nearest_sets[i] or i in nearest_sets[j]:
                    weight = 10 if gap <= 12 else 1
                    graph[a['track']][b['track']] += weight
                    graph[b['track']][a['track']] += weight
    colors = {}
    usage = Counter()
    while len(colors) < len(graph):
        node = max((n for n in graph if n not in colors), key=lambda n: (len({colors[v] for v in graph[n] if v in colors}), sum(graph[n].values()), -int(n)))
        color = min(range(len(PALETTE)), key=lambda c: (sum(w for v, w in graph[node].items() if colors.get(v) == c), usage[c], c))
        colors[node] = color
        usage[color] += 1
    conflicts = sum(1 for a in graph for b in graph[a] if int(a)<int(b) and colors[a]==colors[b])
    folder = OUT / key
    folder.mkdir(parents=True, exist_ok=True)
    svgs = []
    for offset, cells in enumerate(frames):
        polygons = []
        markers = []
        for c in cells:
            fill = PALETTE[colors[c['track']]] if c['track'] is not None else '#B8BEC6'
            points = ' '.join(f'{x:.2f},{y:.2f}' for x,y in c['points'])
            title = f"Track {c['track']} · Spot {c['id']}" if c['track'] is not None else f"Untracked · Spot {c['id']}"
            if c['points']:
                polygons.append(f'<polygon points="{points}" fill="{fill}" stroke="#263344" stroke-width="0.85" stroke-linejoin="round"><title>{title}</title></polygon>')
            markers.append(f'<circle cx="{c["x"]:.3f}" cy="{c["y"]:.3f}" r="2.3" fill="#152335" stroke="#FFFFFF" stroke-width="0.8" vector-effect="non-scaling-stroke"><title>{title} · Tracked centroid' + (' · No stored outline' if not c['points'] else '') + '</title></circle>')
        svg = f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{left} {top} {cw} {ch}" width="{cw:g}" height="{ch:g}"><rect x="{left}" y="{top}" width="{cw}" height="{ch}" fill="#FAFBFD"/>'+''.join(polygons)+''.join(markers)+'</svg>'
        (folder / f'frame_{start+offset:04d}.svg').write_text(svg)
        svgs.append(svg)
    result = dict(key=key, label=label, source=str(path.relative_to(ROOT)), field=path.stem, width=width, height=height, totalFrames=nframes, start=start, crop=crop, tracks=len(colors), counts=[sum(bool(c['points']) for c in f) for f in frames], centroidCounts=[len(f) for f in frames], clipLength=CLIP_LENGTH, nuclearBoundariesAvailable=False, missingContours=sum(missing.values()), untracked=sum(c['track'] is None for f in frames for c in f), neighborConflicts=conflicts, trackColors={t: PALETTE[c] for t,c in colors.items()}, svgs=svgs)
    assert len(svgs) == CLIP_LENGTH and all(result['counts']), f'Empty frames in {key}'
    print(json.dumps({k:v for k,v in result.items() if k not in ('svgs','trackColors')}), flush=True)
    return result


HTML = '''<!doctype html>
<html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>LFCT · Cell contours at 20×</title>
<style>
*{box-sizing:border-box}body{margin:0;background:#edf1f5;color:#19283b;font:15px system-ui,sans-serif}main{max-width:1500px;margin:auto;padding:32px}h1{font-size:30px;margin:6px 0 12px}h2{font-size:20px;margin:0}.eyebrow{letter-spacing:.13em;font-size:12px;color:#53667d}.intro{max-width:1000px;line-height:1.65}.grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:22px;margin-top:24px}.card{background:white;border:1px solid #d7dfe8;border-radius:14px;overflow:hidden;box-shadow:0 3px 12px #22334408}.head{padding:18px 20px}.meta{font-size:13px;color:#58697c;margin-top:6px}.viewport{aspect-ratio:612/452;background:#fafbfd;overflow:hidden;touch-action:none;cursor:grab;border-block:1px solid #e0e6ed}.viewport:focus-visible{outline:3px solid #3277a8;outline-offset:-3px}.viewport:active{cursor:grabbing}.viewport svg{display:block;width:100%;height:100%;user-select:none}.controls{padding:14px 18px}.row{display:flex;align-items:center;gap:9px;flex-wrap:wrap}.slider{width:100%;margin:12px 0;accent-color:#27678b}button{font:inherit;border:1px solid #c5d0dc;border-radius:6px;background:#fff;padding:5px 10px;color:#22384d;cursor:pointer}button:hover{background:#eaf2f8}button:disabled{opacity:.4;cursor:default}.readout{font-size:13px;font-variant-numeric:tabular-nums;margin-left:auto}.zoom{font-size:13px;min-width:42px;text-align:center}details{background:#fff;border:1px solid #d7dfe8;border-radius:10px;padding:16px 20px;margin-top:24px;line-height:1.65}summary{cursor:pointer;font-weight:600}.swatch{display:inline-block;width:16px;height:16px;border:1px solid #263344;border-radius:3px;vertical-align:middle;margin-right:5px}.legend{margin-top:12px;color:#526276;font-size:13px}table{border-collapse:collapse;width:100%;font-size:13px}th,td{text-align:left;padding:8px;border-bottom:1px solid #e0e6ed}code{font-size:12px;overflow-wrap:anywhere}@media(max-width:850px){main{padding:18px}.grid{grid-template-columns:1fr}}a{color:#246687}
</style>
<main><div class="eyebrow">LFCT / TRACKMATE / 20×</div><h1>Cell contours through time</h1>
<p class="intro">Four cell lines, one field each. Each view shows the <b>full 1224 × 904 pixel field</b> over <b>50 consecutive frames</b> around the movie midpoint. Track colors stay fixed through each clip. Dark dots with white rims mark the stored tracking centroids.</p>
<p class="intro"><b>Explore:</b> move a slider or scroll over a view to change frames. Use − / + to zoom, then drag to pan. Ctrl/⌘ + wheel also zooms. Focus a view and use ← / → for frames. Hover a cell to see its Track ID.</p>
<div class="legend" id="legend"></div><div class="grid" id="grid"></div>
<details><summary>Sampling, colors, and source data</summary><p>Fields were selected deterministically by lowest numeric field number, without screening for density or appearance. Each view includes the entire 1224 × 904 pixel image extent. The same fields as the earlier report are retained.</p>
<p>The first selected zero-based frame is floor((total frames − 50) / 2). For odd-length movies, the clip midpoint is half a frame earlier than the movie midpoint. Frame counters show both one-based movie frame and original zero-based XML frame; timing is shown in frames (the XML metadata specifies a 1-second interval).</p>
<p>Original ROI polygon offsets are added to spot positions. Track membership is recovered from all TrackMate track edges, including tracks outside the saved filtered-track selection. Eight fill colors are assigned using a weighted neighbor graph accumulated over the clip: bounding-box gaps ≤12 pixels receive priority, supplemented by nearby centers (up to three nearest neighbors per detection). A greedy saturation-based coloring favors different colors for neighbors; colors are reused across distant tracks and independently across fields. A shared TrackMate ID, including branches, retains one color. Hover labels identify individual tracks. All outlines use the same dark color; gray indicates detections without a track ID. Dots mark POSITION_X / POSITION_Y from the tracking file, without recomputing polygon centroids. Detections lacking valid ROI polygons remain visible as centroid dots only. Separate nuclear boundaries are not available in these four XML files: each detection stores at most one ROI polygon, and PERIMETER is a scalar measurement of that detection. Second-channel intensity features do not contain nuclear boundary coordinates.</p><div id="methods"></div><p>This report embeds all 200 SVG images and requires no network access. Individual SVG files and a sampling/color manifest accompany the report in <code>lfct_20x_contour_viewer/</code>.</p></details></main>
<script id="data" type="application/json">__DATA__</script>
<script>
const data=JSON.parse(document.getElementById('data').textContent);
document.getElementById('legend').innerHTML=data.palette.map(c=>`<span class="swatch" style="background:${c}"></span>`).join('')+' Fixed track palette · dark outlines · white-rimmed dots = tracked centroids · gray = untracked';
for(const d of data.sets){
 const card=document.createElement('section');card.className='card';
 card.innerHTML=`<div class="head"><h2>${d.label}</h2><div class="meta">${d.field} · ${d.width} × ${d.height} px · ${d.tracks} tracks in clip</div></div><div class="viewport" tabindex="0" aria-label="${d.label} interactive contour view"></div><div class="controls"><div class="row"><button class="prev" aria-label="Previous frame">←</button><button class="play" aria-label="Play frames">Play</button><button class="next" aria-label="Next frame">→</button><span class="readout"></span></div><input class="slider" aria-label="${d.label} frame" type="range" min="0" max="${d.svgs.length-1}" value="${Math.floor((d.svgs.length-1)/2)}"><div class="row"><button class="minus" aria-label="Zoom out">−</button><span class="zoom"></span><button class="plus" aria-label="Zoom in">+</button><button class="reset">Reset view</button><button class="save">Save SVG</button><span class="meta count"></span></div></div>`;
 document.getElementById('grid').append(card);
 const q=s=>card.querySelector(s),viewport=q('.viewport'),slider=q('.slider');
 let frame=Math.floor((d.svgs.length-1)/2), zoom=1, cx=d.crop[0]+d.crop[2]/2, cy=d.crop[1]+d.crop[3]/2, timer=null, drag=null, wheel=0;
 const stop=()=>{clearInterval(timer);timer=null;q('.play').textContent='Play';q('.play').setAttribute('aria-label','Play frames');};
 function transform(){const w=d.crop[2]/zoom,h=d.crop[3]/zoom;cx=Math.max(d.crop[0]+w/2,Math.min(d.crop[0]+d.crop[2]-w/2,cx));cy=Math.max(d.crop[1]+h/2,Math.min(d.crop[1]+d.crop[3]-h/2,cy));viewport.firstElementChild.setAttribute('viewBox',`${cx-w/2} ${cy-h/2} ${w} ${h}`);q('.zoom').textContent=`${zoom.toFixed(1)}×`;q('.minus').disabled=zoom<=1;q('.plus').disabled=zoom>=8;}
 function show(i){frame=Math.max(0,Math.min(d.svgs.length-1,i));slider.value=frame;viewport.innerHTML=d.svgs[frame];q('.readout').textContent=`${frame+1}/${d.svgs.length} · Frame ${d.start+frame+1}/${d.totalFrames} · XML ${d.start+frame}`;q('.count').textContent=`${d.counts[frame]} contours · ${d.centroidCounts[frame]} centroids`;q('.prev').disabled=frame===0;q('.next').disabled=frame===d.svgs.length-1;transform();}
 function step(delta){stop();show(frame+delta);}
 function scale(factor){zoom=Math.max(1,Math.min(8,zoom*factor));transform();}
 slider.oninput=()=>{stop();show(Number(slider.value));};q('.prev').onclick=()=>step(-1);q('.next').onclick=()=>step(1);
 q('.play').onclick=()=>{if(timer){stop();return;}q('.play').textContent='Pause';q('.play').setAttribute('aria-label','Pause frames');timer=setInterval(()=>show((frame+1)%d.svgs.length),160);};
 q('.minus').onclick=()=>scale(1/1.25);q('.plus').onclick=()=>scale(1.25);q('.reset').onclick=()=>{zoom=1;transform();};
 q('.save').onclick=()=>{const url=URL.createObjectURL(new Blob([d.svgs[frame]],{type:'image/svg+xml'}));const a=document.createElement('a');a.href=url;a.download=`${d.key}_${d.field}_frame_${d.start+frame}.svg`;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);};
 viewport.addEventListener('wheel',e=>{e.preventDefault();if(e.ctrlKey||e.metaKey){scale(Math.exp(-e.deltaY*.003));return;}wheel+=e.deltaY*(e.deltaMode===1?16:e.deltaMode===2?400:1);if(Math.abs(wheel)>=25){const delta=Math.trunc(wheel/25);wheel-=delta*25;step(Math.max(-5,Math.min(5,delta)));}},{passive:false});
 viewport.addEventListener('keydown',e=>{if(['ArrowRight','ArrowDown','ArrowLeft','ArrowUp','+','=','-','Home'].includes(e.key)){e.preventDefault();if(e.key.startsWith('Arrow'))step(['ArrowRight','ArrowDown'].includes(e.key)?1:-1);else if(e.key==='Home'){zoom=1;transform();}else scale(e.key==='-'?1/1.25:1.25);}});
 viewport.addEventListener('pointerdown',e=>{if(e.button!==0)return;viewport.focus();drag={x:e.clientX,y:e.clientY,cx,cy};viewport.setPointerCapture(e.pointerId);});
 viewport.addEventListener('pointermove',e=>{if(!drag)return;const r=viewport.getBoundingClientRect();cx=drag.cx-(e.clientX-drag.x)*d.crop[2]/zoom/r.width;cy=drag.cy-(e.clientY-drag.y)*d.crop[3]/zoom/r.height;transform();});
 for(const event of ['pointerup','pointercancel','lostpointercapture'])viewport.addEventListener(event,()=>drag=null);
 show(frame);
}
document.getElementById('methods').innerHTML='<table><thead><tr><th>Cell line / source</th><th>XML frames</th><th>Missing outlines</th><th>Untracked detections</th><th>Same-color neighbor pairs</th></tr></thead><tbody>'+data.sets.map(d=>`<tr><td>${d.label}<br><code>${d.source}</code></td><td>${d.start}–${d.start+d.svgs.length-1}</td><td>${d.missingContours}</td><td>${d.untracked}</td><td>${d.neighborConflicts}</td></tr>`).join('')+'</tbody></table><p>Missing and untracked counts are detection counts summed across the 50 frames. Neighbor-pair counts are unique track pairs across the clip; zero means no same-color pair in the constructed neighbor graph.</p>';
</script></html>'''

if __name__ == '__main__':
    sets = [build_set(k,v) for k,v in LABELS.items()]
    payload = dict(palette=PALETTE, sets=sets)
    report = ROOT / 'reports/lfct_20x_contour_viewer.html'
    report.write_text(HTML.replace('__DATA__', json.dumps(payload, separators=(',', ':')).replace('</', '<\\/')))
    (OUT / 'manifest.json').write_text(json.dumps(dict(palette=PALETTE, sets=[{k:v for k,v in d.items() if k!='svgs'} for d in sets]), indent=2))
    print(f'Wrote {report} ({report.stat().st_size / 1024**2:.2f} MiB)')
