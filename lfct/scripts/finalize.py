"""Audit every field after the array settles; publish only a complete validated report."""
import sys
from common import BASE,profile_fields,field_row,output_dir,read_json,sha,write_json
from report import build

def main():
    issues=[];frames=0
    for key in profile_fields('full'):
        out=output_dir('full',key)
        try:
            c=read_json(out/'completed.json')
            assert c['complete_field'] and c['frames']==list(range(int(field_row(key)['nframes'])))
            assert c['source']==field_row(key),'Source catalog changed'
            for name,h in c['output_hashes'].items():assert sha(out/name)==h,f'Checksum mismatch: {name}'
            for name,h in c['qc_asset_hashes'].items():assert sha(BASE/'reports/full'/name)==h,f'PNG mismatch: {name}'
            assert read_json(out/'validation.json')['status']=='passed'
            frames+=len(c['frames'])
        except Exception as e:
            issues.append(dict(key=key,error=f'{type(e).__name__}: {e}',failure=read_json(out/'failure.json') if (out/'failure.json').exists() else None))
    result=dict(status='failed' if issues else 'passed',validated_frames=frames,issues=issues)
    write_json(BASE/'results/full/final_audit.json',result)
    if issues:
        print(result,flush=True);return 1
    build('full',bundle=True)
    print(f'FINAL AUDIT PASSED: {frames} pairs; reports/full/index.html and reports/full.zip',flush=True)
    return 0
if __name__=='__main__':sys.exit(main())
