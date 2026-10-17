#!/usr/bin/env python3
"""Verify and extract paired FLD_102 TIFFs, retaining exact source/frame provenance."""
from pathlib import Path
import hashlib, json, re, zipfile, xml.etree.ElementTree as ET
import numpy as np
import tifffile
from PIL import Image, ImageDraw
ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data/lfct_mcf10a_fld102'
SOURCES={'phase':'https://clemson.box.com/s/bgrt5z8g56vrdev9np7ruu2zigb3mff3','rfp':'https://clemson.box.com/s/07cvk67q6x8ymim344tz538qf54romn0'}
def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
    return h.hexdigest()
def display(a):
    lo,hi=np.percentile(a,[1,99.5]);return (np.clip((a.astype(float)-lo)/max(hi-lo,1),0,1)*255).astype('uint8')
def main():
    manifest={'dataset_page':'https://davislab.github.io/LFCT-dataset/mcf10a.html','sources':{},'frames':[]}
    channels={}
    for channel,url in SOURCES.items():
        archive=DATA/'source'/f'{channel}.zip'
        out=DATA/'images'/channel;out.mkdir(parents=True,exist_ok=True)
        files={};excluded=[]
        with zipfile.ZipFile(archive) as z:
            for info in z.infolist():
                name=info.filename
                if name.startswith('__MACOSX/') or info.is_dir() or not name.lower().endswith(('.tif','.tiff')):
                    excluded.append(name);continue
                index=int(re.search(r'_(\d+)\.tif',name)[1])
                assert index not in files
                dst=out/Path(name).name
                if not dst.exists():dst.write_bytes(z.read(info))
                else:assert dst.stat().st_size==info.file_size
                files[index]=dst
        channels[channel]=files
        manifest['sources'][channel]={'url':url,'archive':str(archive.relative_to(ROOT)),'bytes':archive.stat().st_size,'sha256':sha(archive),'excluded_archive_entries':len(excluded),'exclusion_reason':'macOS resource forks and directory entries; no image frames excluded'}
    xml=ROOT/'data/lfct_trackmate_20x/mcf10a/FLD_102.xml'
    doc=ET.parse(xml);meta=doc.find('Settings/ImageData').attrib;n=int(meta['nframes'])
    assert set(channels['phase'])==set(channels['rfp'])==set(range(1,n+1))
    manifest['xml']={'path':str(xml.relative_to(ROOT)),'sha256':sha(xml),'image_data':meta}
    manifest['frame_mapping']='TIFF numeric suffix 1..385 maps to XML FRAME 0..384 by subtracting 1; channels matched by numeric suffix, never ZIP order.'
    for i in range(1,n+1):
        row={'frame':i-1,'image_number':i}
        for channel in SOURCES:
            p=channels[channel][i]
            with tifffile.TiffFile(p) as t:
                assert len(t.pages)==1
                assert t.pages[0].shape==(int(meta['height']),int(meta['width']))
                row[channel]={'path':str(p.relative_to(ROOT)),'dtype':str(t.pages[0].dtype)}
        manifest['frames'].append(row)
    (DATA/'source_manifest.json').write_text(json.dumps(manifest,indent=2))
    frames=[0,192,384];canvas=Image.new('RGB',(1224*2,904*3),'white');draw=ImageDraw.Draw(canvas)
    for j,f in enumerate(frames):
        for k,c in enumerate(SOURCES):
            a=tifffile.imread(channels[c][f+1]);img=Image.fromarray(display(a)).convert('RGB')
            canvas.paste(img,(1224*k,904*j));draw.text((1224*k+12,904*j+12),f'{c} XML frame {f}',fill='red')
    canvas.resize((1224,1356)).save(DATA/'source_preview.jpg',quality=92)
    print(f'Verified {n} matched frame pairs and wrote source_manifest.json',flush=True)
if __name__=='__main__':main()
