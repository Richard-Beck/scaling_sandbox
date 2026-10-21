"""Exact17-digit per-page S11 source CSVs and binary roundtrip proof.

Full native field matrices remain HPC source archives. Portable companions
carry q grids, every plotted scalar/snapshot and explicit source identities.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse
import csv
import gzip
import hashlib
import json
import os
import shutil
import numpy as np


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def read(path):return json.loads(Path(path).read_text())
def write(path,value):Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
def binary(value):return hashlib.sha256(np.asarray(value,dtype='<f8').tobytes(order='C')).hexdigest()


class WideShards:
    def __init__(self,directory,stem,header,limit=64*1024**2):
        self.directory=Path(directory);self.directory.mkdir(parents=True,exist_ok=True)
        self.stem=stem;self.header=header+'\n';self.limit=limit;self.files=[];self.stream=None;self.rows=0;self.byte_count=0

    def write(self,text):
        data=text.encode('ascii')
        if self.stream is None or self.byte_count+len(data)>self.limit:
            self.close();path=self.directory/f'{self.stem}_part{len(self.files)+1:02d}.csv.gz'
            handle=path.open('wb');self.handle=handle;self.stream=gzip.GzipFile(filename='',mode='wb',fileobj=handle,mtime=0,compresslevel=6)
            self.stream.write(self.header.encode('ascii'));self.byte_count=len(self.header);self.rows=0
            self.files.append(dict(path=str(path),rows=0))
        self.stream.write(data);self.byte_count+=len(data);self.rows+=text.count('\n');self.files[-1]['rows']=self.rows

    def close(self):
        if self.stream is not None:
            self.stream.close();self.handle.close();self.stream=None


def save_numeric(path,header,data,fmt):
    with Path(path).open('wb') as raw:
        with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0,compresslevel=6) as stream:
            np.savetxt(stream,data,fmt=fmt,delimiter=',',header=header,comments='')


def worker(run,slot):
    index=read(run/'results/export/native_atlas_index.json');group=index['groups'][slot]
    assert group['condition']==read(run/'manifest.json')['jobs'][slot]['condition']
    context=group['condition']['experiment']
    stem=f'S11_{context}_page_{slot:04d}';archive=run/'results/native_field_CSV_archives';portable=run/'results/export/page_data'
    portable.mkdir(exist_ok=True);(run/'results/csv_workers').mkdir(exist_ok=True)
    selected=[];native_cells=None;global_high=0.
    for number,display in enumerate(group['artifact']['displayed'],1):
        case=next(c for c in group['cases'] if c['id']==display['case_id'])
        assert sha(case['profile_path'])==case['profile_sha256'] and sha(case['scalar_path'])==case['scalar_sha256']
        with np.load(case['profile_path']) as z:
            time=z['time'].copy();field=z['primary'].copy();C=np.broadcast_to(z['circumference'],time.shape).copy()
            q=z['material_positions'].copy() if 'material_positions' in z else z['x'].copy()/C[0]
        assert np.isfinite(field).all() and np.isfinite(C).all()
        if native_cells is None:native_cells=field.shape[1]
        assert native_cells==field.shape[1]
        global_high=max(global_high,float(field.max()))
        with np.load(case['scalar_path']) as z:
            st=z['time'].copy();P=abs(z['P_real']+1j*z['P_imag']);forcing=z['local_input_mean'].copy()
        features=group['features'][group['cases'].index(case)]
        selected.append(dict(number=number,display=display,case=case,time=time,field=field,C=C,q=q,st=st,P=P,forcing=forcing,features=features))
    fields=WideShards(archive,stem,'case_index,native_frame,time,actualC,'+','.join(f'active_cell_{j:04d}' for j in range(native_cells)))
    expected={};selection=[];snapshots=[];scalar=[];grids=[]
    try:
        for s in selected:
            number=s['number'];time,field,C,q=[s[k] for k in ['time','field','C','q']];case=s['case'];display=s['display']
            np.savetxt(fields,np.column_stack([np.full(len(time),number),np.arange(len(time)),time,C,field]),
                fmt=['%d','%d','%.17g','%.17g']+['%.17g']*native_cells,delimiter=',')
            expected[number]=dict(field_sha256=binary(field),time_sha256=binary(time),actualC_sha256=binary(C),rows=len(time),cells=native_cells)
            scalar.append(np.column_stack([np.full(len(s['st']),number),np.arange(len(s['st'])),s['st'],s['P'],s['forcing']]))
            grids.append(np.column_stack([np.full(len(q),number),np.arange(len(q)),q]))
            for snapshot_number,snapshot in enumerate(display['snapshots']):
                i=snapshot['native_frame'];snapshots.append(np.r_[number,snapshot_number,i,time[i],C[i],field[i]])
            row=dict(case_index=number,case_id=case['id'],role=display['role'],seed=case['seed'],independent_realization_id=case['independent_realization_id'],
                source_observation_id=case.get('source_observation_id',case['id']),origin=case.get('origin',case.get('experiment')),
                physical_run_id=case.get('physical_run_id',case['independent_realization_id']),latent_Gaussian_cohort_id=case.get('latent_Gaussian_cohort_id'),
                independent_stochastic_seed_count=int(case['kind']=='noise'),deterministic_control=case['kind']!='noise',
                prefix_and_topup_relation=json.dumps(case.get('prefix_and_topup_relation'),sort_keys=True),
                condition_forced=case['condition_forced'],native_frames=len(time),native_cells=native_cells,scalar_frames=len(s['st']),
                scalar_cadence=case['scalar_cadence'],minimum_actualC=float(C.min()),maximum_actualC=float(C.max()),
                colormap='viridis',concentration_vmin=0.,concentration_vmax=global_high,active_field='active Rac' if case['model']=='holmes_model3' else 'native active',
                field_archive_distribution='HPC source archive; full numeric values exist and are verified locally, outside compact Git exports',
                native_profile_NPZ_local_path=case['profile_path'],native_profile_NPZ_sha256=case['profile_sha256'],
                native_scalar_NPZ_local_path=case['scalar_path'],native_scalar_NPZ_sha256=case['scalar_sha256'],
                basal_profile_NPZ_local_path=case['basal_profile_path'],basal_profile_NPZ_sha256=case['basal_profile_sha256'],
                source_status_local_path=case['source_status_path'],source_status_sha256=case['source_status_sha256'],
                q_coordinate='Native q; physical plotted cell centres are q times actualC at each saved frame',
                snapshot_fractions='0;0.5;0.75;1',snapshot_actual_times=';'.join(format(x['native_saved_time'],'.17g') for x in display['snapshots']),
                first_half_shaded=True,selection_window_start_fraction=.5,main_example_selection=False,
                **group['condition'],**s['features'])
            selection.append(row)
    finally:fields.close()
    checks={n:dict(field=hashlib.sha256(),time=hashlib.sha256(),C=hashlib.sha256(),rows=0) for n in expected}
    for part in fields.files:
        path=Path(part['path']);assert path.stat().st_size<90*1024**2
        with gzip.open(path,'rt') as stream:data=np.loadtxt(stream,delimiter=',',skiprows=1,ndmin=2)
        assert data.shape[1]==native_cells+4 and len(data)==part['rows']
        for number in np.unique(data[:,0]).astype(int):
            rows=data[data[:,0]==number];check=checks[number]
            np.testing.assert_array_equal(rows[:,1],np.arange(check['rows'],check['rows']+len(rows)))
            check['rows']+=len(rows);check['field'].update(np.asarray(rows[:,4:],dtype='<f8').tobytes(order='C'))
            check['time'].update(np.asarray(rows[:,2],dtype='<f8').tobytes());check['C'].update(np.asarray(rows[:,3],dtype='<f8').tobytes())
        part.update(sha256=sha(path),bytes=path.stat().st_size,allocated_bytes=path.stat().st_blocks*512,roundtrip_float64_equal=True)
    for number,expected_case in expected.items():
        actual=checks[number];assert actual['rows']==expected_case['rows']
        assert actual['field'].hexdigest()==expected_case['field_sha256'] and actual['time'].hexdigest()==expected_case['time_sha256'] and actual['C'].hexdigest()==expected_case['actualC_sha256']
    companions={
        'snapshots':(portable/(stem+'_snapshots.csv.gz'),'case_index,snapshot_index,native_frame,time,actualC,'+','.join(f'active_cell_{j:04d}' for j in range(native_cells)),np.stack(snapshots),['%d']*3+['%.17g']*(native_cells+2)),
        'scalar_traces':(portable/(stem+'_scalar_traces.csv.gz'),'case_index,scalar_frame,time,basal_M1_normalized_P_magnitude,native_spatial_mean_input',np.vstack(scalar),['%d','%d','%.17g','%.17g','%.17g']),
        'q_grid':(portable/(stem+'_q_grid.csv.gz'),'case_index,native_cell_index,q',np.vstack(grids),['%d','%d','%.17g'])}
    companion_records={}
    for name,(path,header,data,fmt) in companions.items():
        save_numeric(path,header,data,fmt)
        with gzip.open(path,'rt') as stream:actual=np.loadtxt(stream,delimiter=',',skiprows=1,ndmin=2)
        np.testing.assert_array_equal(actual,data)
        companion_records[name]=dict(path=str(path),sha256=sha(path),rows=len(data),columns=data.shape[1],bytes=path.stat().st_size,float64_roundtrip_equal=True)
    selection_path=portable/(stem+'_selection.csv')
    with selection_path.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(selection[0]));writer.writeheader();writer.writerows(selection)
    receipt=dict(passed=True,slot=slot,condition=group['condition'],selected_cases=[s['case']['id'] for s in selected],
        native_cell_time_values=sum(s['field'].size for s in selected),full_field_CSV_distribution='HPC-only verified source archives; compact Git carries exact inventory, hashes and selection CSV',
        native_field_CSV_shards=fields.files,expected_native_binary_hashes=expected,companions=companion_records,
        selection_CSV=dict(path=str(selection_path),sha256=sha(selection_path)),source_atlas_PNG_sha256=group['artifact']['png_sha256'],
        decimal_significant_digits=17,float64_roundtrip_equal=True,sourcefield_resampling=False,sourcefield_quantization=False,export_source_sha256=sha(__file__))
    write(run/'results/csv_workers'/f'{slot:04d}.json',receipt)
    return dict(slot=slot,values=receipt['native_cell_time_values'],archive_bytes=sum(p['bytes'] for p in fields.files))


def collect(run):
    index=read(run/'results/export/native_atlas_index.json');out=run/'results/export/page_data';out.mkdir(exist_ok=True)
    pdf_rows=list(csv.DictReader((run/'results/export/pdf/page_selection_inventory.csv').open()))
    assert len(pdf_rows)==len(index['groups']);records=[];flat=[]
    for slot,(group,pdf_page) in enumerate(zip(index['groups'],pdf_rows)):
        receipt=read(run/'results/csv_workers'/f'{slot:04d}.json');assert receipt['passed'] and receipt['float64_roundtrip_equal']
        assert receipt['source_atlas_PNG_sha256']==pdf_page['png_sha256']==group['artifact']['png_sha256']
        for path in receipt['native_field_CSV_shards']:
            assert sha(path['path'])==path['sha256']
            flat.append(dict(slot=slot,pdf_volume=pdf_page['volume'],pdf_page=int(pdf_page['page']),**group['condition'],
                local_HPC_source_CSV_path=path['path'],source_CSV_sha256=path['sha256'],source_CSV_bytes=path['bytes'],allocated_bytes=path['allocated_bytes'],
                native_rows=path['rows'],native_cells=group['artifact']['displayed'][0]['native_cells'],decimal_significant_digits=17,float64_roundtrip_equal=True,
                distribution='HPC-only full-field source archive; use this inventory for retrieval, not an absent GitHub link',
                selection_CSV=Path(receipt['selection_CSV']['path']).name,selection_CSV_sha256=receipt['selection_CSV']['sha256']))
        for item in list(receipt['companions'].values())+[receipt['selection_CSV']]:assert sha(item['path'])==item['sha256']
        records.append(dict(receipt,pdf_volume=pdf_page['volume'],pdf_page=int(pdf_page['page'])))
    with (out/'native_field_archive_inventory.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(flat[0]));writer.writeheader();writer.writerows(flat)
    shutil.copy2(__file__,out/'frozen_csv_export.py')
    summary=dict(passed=True,generated_utc=datetime.now(timezone.utc).isoformat(),analysis_arm=index['summary']['stage'],pages=len(records),
        selected_histories=sum(len(r['selected_cases']) for r in records),native_cell_time_values=sum(r['native_cell_time_values'] for r in records),
        full_field_archive_bytes=sum(r['source_CSV_bytes'] for r in flat),full_field_archive_allocated_bytes=sum(r['allocated_bytes'] for r in flat),
        maximum_archive_file_bytes=max(r['source_CSV_bytes'] for r in flat),native_rows=sum(r['native_rows'] for r in flat),
        sourcefield_resampling=False,sourcefield_quantization=False,float64_roundtrip_equal_every_field_time_actualC=True,
        full_native_field_CSVs='HPC-only source archives retained outside compact export; paths and exact hashes are indexed',
        source_atlas_complete_sha256=sha(run/'results/export/complete.json'),source_PDF_audit_sha256=sha(run/'results/export/pdf/lossless_pdf_audit.json'),
        archive_inventory_sha256=sha(out/'native_field_archive_inventory.csv'),export_source_sha256=sha(__file__),native_integrations=0,main_example_selection=False)
    write(out/'page_data_index.json',dict(summary=summary,pages=records));write(out/'complete.json',dict(summary,page_data_index_sha256=sha(out/'page_data_index.json')))
    return summary


def main():
    p=argparse.ArgumentParser();p.add_argument('action',choices=['worker','collect']);p.add_argument('--run',type=Path,required=True);p.add_argument('--slot',type=int);args=p.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('Full exact native CSV export uses SLURM resources')
    print(json.dumps(worker(args.run.resolve(),args.slot) if args.action=='worker' else collect(args.run.resolve())))


if __name__=='__main__':main()
