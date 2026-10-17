#!/usr/bin/env python3
"""Catalog public LFCT 20x sources and stage only requested TIFF frames in scratch."""
import argparse,csv,gzip,hashlib,io,json,re,shutil,subprocess,time,urllib.request,zipfile
from pathlib import Path
import xml.etree.ElementTree as ET
from common import BASE,REPO,sha,write_json,read_json,csv_write,field_row,frame_selection,config


def fetch_text(url):
    return subprocess.check_output(['curl','--silent','--show-error','--location','--fail','--retry','3','--max-time','90',url],text=True)

def discover():
    rows=[]
    for cell in ['hek293','mcf10a','mdamb','u87']:
        page=f'https://davislab.github.io/LFCT-dataset/{cell}.html'
        text=fetch_text(page)
        images=text.split('<summary>Images</summary>',1)[1].split('<summary>Annotation',1)[0]
        section=images.split('<summary>20',1)[1]
        cards=re.findall(r'<h3>(FLD_\d+)</h3>\s*<p>(.*?)</p>',section,re.S)
        links={name:dict((label,url) for url,label in re.findall(r'<a href="([^"]+)">(pc|rfp)</a>',body)) for name,body in cards}
        for xml in sorted((REPO/'data/lfct_trackmate_20x'/cell).glob('*.xml'),key=lambda p:int(p.stem.split('_')[-1])):
            data=xml.read_bytes();meta=ET.fromstring(re.search(rb'<ImageData\b[^>]*/>',data)[0]).attrib
            assert set(links[xml.stem])=={'pc','rfp'}
            rows.append(dict(key=f'{cell}/{xml.stem}',cell_line=cell,field=xml.stem,nframes=meta['nframes'],width=meta['width'],height=meta['height'],phase_url=links[xml.stem]['pc'],rfp_url=links[xml.stem]['rfp'],dataset_page=page,original_xml=str(xml.relative_to(REPO)),xml_sha256=hashlib.sha256(data).hexdigest()))
    csv_write(BASE/'data/sources.csv',rows)
    print(f'Cataloged {len(rows)} fields / {sum(int(r["nframes"]) for r in rows)} frames',flush=True)


def annotation(row):
    dst=BASE/'data/annotations'/row['cell_line']/(row['field']+'.xml.gz')
    if dst.exists():
        with gzip.open(dst,'rb') as f:data=f.read()
        assert hashlib.sha256(data).hexdigest()==row['xml_sha256'];return dst
    src=REPO/row['original_xml'];assert sha(src)==row['xml_sha256']
    dst.parent.mkdir(parents=True,exist_ok=True);tmp=dst.with_suffix('.gz.part')
    with src.open('rb') as source,gzip.open(tmp,'wb',compresslevel=6) as target:shutil.copyfileobj(source,target)
    tmp.replace(dst);return dst


def box_download(url,dst):
    before=time.monotonic();text=fetch_text(url)
    m=re.search(r'"itemID"\s*:\s*"?(\d+)',text)
    if not m:raise ValueError('No public Box item ID: '+url)
    shared=url.rstrip('/').split('/')[-1]
    direct=f'https://clemson.box.com/index.php?rm=box_download_shared_file&shared_name={shared}&file_id=f_{m[1]}'
    tmp=dst.with_suffix('.zip.part')
    subprocess.run(['curl','--silent','--show-error','--location','--fail','--retry','3','--max-time','1200','--output',str(tmp),direct],check=True)
    assert zipfile.is_zipfile(tmp),'Download is not a ZIP archive'
    tmp.replace(dst)
    return dict(url=url,sha256=sha(dst),bytes=dst.stat().st_size,download_seconds=time.monotonic()-before)


def source_indices(key,n,available):
    override=config().get('source_mapping_overrides',{}).get(key)
    if override:
        assert available==set(range(override['expected_suffix_count'])), 'Override archive changed'
        skipped=override['skip_suffixes'];indices=sorted(available-set(skipped))
        assert len(indices)==n
        return indices,skipped
    if available==set(range(1,n+1)):return list(range(1,n+1)),[]
    if available==set(range(n)):return list(range(n)),[]
    raise ValueError(f'{key}: image indices differ from XML frame range: expected {n}; actual {sorted(available)}')


def prepare(row,frames,scratch):
    scratch=Path(scratch);scratch.mkdir(parents=True,exist_ok=True)
    ann=annotation(row);before=time.monotonic();manifest=dict(key=row['key'],frames=frames,width=int(row['width']),height=int(row['height']),annotation=str(ann),xml_sha256=row['xml_sha256'],sources={},images={})
    # Reuse only the explicitly known MCF10A local source, verifying its source record.
    legacy=REPO/'data/lfct_mcf10a_fld102/source_manifest.json'
    if row['key']=='mcf10a/FLD_102' and legacy.exists():
        old=read_json(legacy)
        assert old['xml']['sha256']==row['xml_sha256']
        for c in ['phase','rfp']:
            assert old['sources'][c]['url']==row[c+'_url']
            manifest['sources'][c]={**old['sources'][c],'reuse':'existing source TIFFs; no duplicate download'}
            manifest['images'][c]={str(f):str(REPO/old['frames'][f][c]['path']) for f in frames}
    else:
        for c in ['phase','rfp']:
            archive=scratch/(c+'.zip');meta_path=scratch/(c+'_download.json')
            if archive.exists() and meta_path.exists():
                info=read_json(meta_path);assert info['url']==row[c+'_url'] and info['sha256']==sha(archive)
            else:info=box_download(row[c+'_url'],archive);write_json(meta_path,info)
            image_dir=scratch/c;image_dir.mkdir(exist_ok=True);mapping={};excluded=0
            with zipfile.ZipFile(archive) as z:
                for obj in z.infolist():
                    name=obj.filename
                    if obj.is_dir() or name.startswith('__MACOSX/') or Path(name).name.startswith('._') or not name.lower().endswith(('.tif','.tiff')):excluded+=1;continue
                    m=re.search(r'_(\d+)\.tiff?$',name,re.I)
                    if not m:raise ValueError('Unrecognized time suffix: '+name)
                    number=int(m[1]);assert number not in mapping,'Duplicate image index'
                    mapping[number]=name
                n=int(row['nframes'])
                indices,skipped=source_indices(row['key'],n,set(mapping))
                manifest['images'][c]={}
                for f in frames:
                    dst=image_dir/f'{f:04d}.tif';dst.write_bytes(z.read(mapping[indices[f]]));manifest['images'][c][str(f)]=str(dst)
            manifest['sources'][c]={**info,'excluded_entries':excluded,'xml_to_numeric_suffix':indices,'skipped_image_suffixes':skipped,'available_frames':len(mapping)}
    import tifffile
    for c in ['phase','rfp']:
        for f,path in manifest['images'][c].items():
            with tifffile.TiffFile(path) as t:
                assert len(t.pages)==1 and t.pages[0].shape==(manifest['height'],manifest['width'])
    manifest['seconds']=time.monotonic()-before
    write_json(scratch/'prepared.json',manifest)
    return manifest

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--discover',action='store_true');a=p.parse_args()
    if a.discover:discover()
