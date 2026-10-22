"""Portable S11 PDF volumes with independently decoded lossless image proof."""
from pathlib import Path
from datetime import datetime,timezone
import argparse
import csv
import hashlib
import html
import json
import os
import shutil
import zlib
import numpy as np
from PIL import Image
from pypdf import PdfReader,PdfWriter
from pypdf.generic import DictionaryObject,NameObject,NumberObject,EncodedStreamObject,DecodedStreamObject


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def read(path):return json.loads(Path(path).read_text())


def add_exact_png(writer,path):
    with Image.open(path) as image:
        pixels=np.asarray(image.convert('RGBA'))
        if not np.all(pixels[:,:,3]==255):raise ValueError('PDF exporter requires opaque native PNG; alpha may not be discarded')
        width,height=image.size;raw=pixels[:,:,:3].tobytes()
    encoded=EncodedStreamObject();encoded._data=zlib.compress(raw,9)
    encoded.update({NameObject('/Type'):NameObject('/XObject'),NameObject('/Subtype'):NameObject('/Image'),
        NameObject('/Filter'):NameObject('/FlateDecode'),NameObject('/Width'):NumberObject(width),NameObject('/Height'):NumberObject(height),
        NameObject('/ColorSpace'):NameObject('/DeviceRGB'),NameObject('/BitsPerComponent'):NumberObject(8)})
    physical_width,physical_height=width*72/170,height*72/170
    page=writer.add_blank_page(width=physical_width,height=physical_height)
    page[NameObject('/Resources')]=DictionaryObject({NameObject('/XObject'):DictionaryObject({NameObject('/NativePNG'):writer._add_object(encoded)})})
    content=DecodedStreamObject();content.set_data(f'q\n{physical_width:.9f} 0 0 {physical_height:.9f} 0 0 cm\n/NativePNG Do\nQ\n'.encode())
    page[NameObject('/Contents')]=writer._add_object(content)
    return dict(width=width,height=height,decoded_RGB_sha256=hashlib.sha256(raw).hexdigest(),opaque_alpha_255=True,encoded_bytes=len(encoded._data))


def main():
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);args=p.parse_args();run=args.run.resolve()
    if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('Complete atlas PDF export uses SLURM resources')
    out=run/'results/export';complete=read(out/'complete.json');index=read(out/'native_atlas_index.json');assert complete['passed']
    target=out/'pdf';target.mkdir(exist_ok=False)
    shutil.copy2(__file__,target/'frozen_pdf_export.py')
    context=complete['stage']
    volumes=[];page_rows=[];writer=PdfWriter();writer.add_metadata({'/Title':f'Figure4 S11 {context} native atlas; predeclared medoid and distinct random realizations'})
    current=[];estimated_bytes=0
    def flush():
        nonlocal writer,current,estimated_bytes
        if not current:return
        number=len(volumes)+1;name=f'S11_native_atlas_{context}_volume{number:02d}.pdf';path=target/name
        writer.write(path);assert path.stat().st_size<100*1024**2
        reader=PdfReader(path);assert len(reader.pages)==len(current)
        for page_number,(page,record) in enumerate(zip(reader.pages,current),1):
            obj=page['/Resources']['/XObject']['/NativePNG'].get_object()
            assert obj['/Filter']=='/FlateDecode' and obj['/ColorSpace']=='/DeviceRGB' and obj['/BitsPerComponent']==8
            pixels=obj.get_data();assert hashlib.sha256(pixels).hexdigest()==record['decoded_RGB_sha256']
            assert obj['/Width']==record['width'] and obj['/Height']==record['height']
            record.update(volume=name,page=page_number,PDF_decoded_pixel_equal=True)
            page_rows.append(record)
        volumes.append(dict(file=name,sha256=sha(path),bytes=path.stat().st_size,pages=len(current),lossless_flate=True,decoded_image_pixels_equal=True))
        writer=PdfWriter();writer.add_metadata({'/Title':f'Figure4 S11 {context} native atlas; predeclared medoid and distinct random realizations'})
        current=[];estimated_bytes=0
    for group in index['groups']:
        artifact=group['artifact'];condition=group['condition'];png=Path(artifact['png'])
        assert sha(png)==artifact['png_sha256']
        # Compressed stream-size estimate with explicit per-page overhead.
        # Each volume also caps its page count at50 for a usable download.
        with Image.open(png) as image:
            estimate=len(zlib.compress(image.convert('RGB').tobytes(),9))+10000
        if current and (len(current)>=50 or estimated_bytes+estimate>90*1024**2):flush()
        data=add_exact_png(writer,png)
        roles={d['role']:d['case_id'] for d in artifact['displayed']}
        data.update(condition,png_local_path=str(png),png_sha256=artifact['png_sha256'],webp='../webp/'+Path(artifact['webp']).name,
            webp_sha256=artifact['webp_sha256'],medoid_case_id=roles['Summary-feature medoid'],random_case_id=roles.get('Deterministic distinct random'),
            candidate_physical_runs=len(group['cases']),selected_native_histories=len(artifact['displayed']),
            independent_stochastic_realizations_in_condition=len(group['cases']) if condition['kind']=='noise' else 0,
            deterministic_controls=len(group['cases']) if condition['kind']!='noise' else 0,
            native_origin_stages=json.dumps(sorted({c.get('origin',c.get('experiment')) for c in group['cases']})),
            no_distinct_random=len(roles)==1,main_example_selection=False,native_integrations=0)
        current.append(data);estimated_bytes+=estimate
    flush()
    with (target/'page_selection_inventory.csv').open('w',newline='') as stream:
        csv_writer=csv.DictWriter(stream,fieldnames=list(page_rows[0]));csv_writer.writeheader();csv_writer.writerows(page_rows)
    receipt=dict(passed=True,generated_utc=datetime.now(timezone.utc).isoformat(),analysis_arm=context,groups=len(page_rows),volumes=volumes,
        noisy_candidate_physical_runs=sum(r['candidate_physical_runs'] for r in page_rows if r['kind']=='noise'),
        deterministic_control_runs=sum(r['deterministic_controls'] for r in page_rows),
        independence_scope='Noisy counts are per fixed condition/forcing arm. Paired sigma branches, observation companions and retained long prefixes have explicit canonical lineage accounting in the final three-context master.',
        opaque_PNG_alpha_verified=True,all_decoded_PDF_RGB_pixels_equal_to_source_PNG=True,source_field_transformations=0,lossy_image_compression=False,
        maximum_artifact_bytes=max(v['bytes'] for v in volumes),native_NPZ_and_original_PNG='Retained local/HPC archives; no deletion',
        main_example_selection=False,native_integrations=0,source_atlas_complete_sha256=sha(out/'complete.json'),
        page_inventory_sha256=sha(target/'page_selection_inventory.csv'),export_source_sha256=sha(__file__))
    (target/'lossless_pdf_audit.json').write_text(json.dumps(receipt,indent=2)+'\n')
    links=''.join(f'<li><a href="{html.escape(v["file"])}">Volume{i+1}: {v["pages"]} condition pages</a> ({v["bytes"]/1024**2:.1f}MiB)</li>' for i,v in enumerate(volumes))
    (target/'index.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><title>S11 '+html.escape(context)+' lossless PDF atlas</title><body style="font:18px/1.5 system-ui;max-width:1000px;margin:40px auto"><h1>S11 · '+html.escape(context)+' native atlas PDF volumes</h1><p>'+html.escape(receipt['generated_utc'])+'. Predeclared summary-feature medoid and distinct random case for every available condition. No main example is selected. All native fields use physical arc; Spring uses actual circumference. Each decoded PDF image is pixel-identical to its opaque PNG source, with lossless Flate compression.</p><ul>'+links+'</ul><p><a href="page_selection_inventory.csv">Page and selection inventory</a> · <a href="lossless_pdf_audit.json">Decoded-pixel proof</a> · <a href="frozen_pdf_export.py">Frozen export source</a> · <a href="../index.html">Searchable WebP atlas and native-source links</a></p></body></html>')
    print(json.dumps(dict(passed=True,pages=len(page_rows),volumes=volumes)))


if __name__=='__main__':main()
