"""Inspect native Holmes profiles behind a finite half-share calibration limit."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
from scipy.signal import find_peaks
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def digest(path):
    result=hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda:handle.read(1024*1024),b""):result.update(block)
    return result.hexdigest()


def profile(path,index):
    with np.load(path) as data:
        return {"x":data["x"].copy(),"rac":data["field_rac_active"][index].copy(),
                "rho":data["field_rho_active"][index].copy()}


def describe(x,y,length):
    floor=float(y.min());excess=y-floor;amplitude=float(excess.max());dx=length/len(y)
    peaks,_=find_peaks(np.r_[floor,y,floor],prominence=.05*amplitude)
    peaks=peaks[(peaks>0)&(peaks<=len(y))]-1
    regions={}
    for fraction in (.1,.5,.75):
        mask=excess>=fraction*amplitude
        changes=np.diff(np.r_[False,mask,False].astype(int))
        starts=np.flatnonzero(changes==1);stops=np.flatnonzero(changes==-1)
        regions[str(fraction)]=[[float(a*dx),float(b*dx)] for a,b in zip(starts,stops)]
    # A 201-cell domain has a midpoint-straddling cell; split its amount equally.
    edges=np.linspace(0.,length,len(y)+1)
    left_weights=np.maximum(0.,np.minimum(edges[1:],length/2)-edges[:-1])
    left=float(np.sum(excess*left_weights));right=float(np.sum(excess*(dx-left_weights)))
    return dict(floor=floor,peak=float(y.max()),peak_positions=x[peaks].tolist(),
        left_max=float(y[x<length/2].max()),right_max=float(y[x>length/2].max()),
        midpoint_value=float(y[len(y)//2]),regions=regions,
        excess_centroid=float(np.sum(x*excess)/np.sum(excess)),
        left_excess=left,right_excess=right,D_half=(right-left)/(right+left),
        interpretation_scope="Threshold components and native profile geometry; region identity is not inferred from a half-domain share.")


def audit(root,output):
    root=root.resolve();output=output.resolve();output.mkdir(parents=True,exist_ok=True)
    manifest=json.loads((root/"manifest.json").read_text());rows=[];sources={}
    fig,axes=plt.subplots(4,3,figsize=(12,10),constrained_layout=True)
    for row,case in enumerate(manifest["cases"]):
        length=case["native_length"];worker=root/"workers"/f"{case['task']:04d}"
        folder=worker/"evaluations/100000"
        paths=[root/case["checkpoint"],folder/"phase_000.npz",folder/"phase_001.npz"]
        snapshots=[]
        for col,(path,label) in enumerate(zip(paths,("Qualified basal","Cue ON: 100,000","Cue-free recovery: 10,000"))):
            native=profile(path,-1);rac=describe(native["x"],native["rac"],length)
            rho=describe(native["x"],native["rho"],length)
            sources[str(path)]=digest(path)
            snapshots.append(dict(label=label,source=str(path),rac=rac,rho=rho,
                profiles={key:value.tolist() for key,value in native.items()}))
            ax=axes[row,col];ax.plot(native["x"],native["rac"],label="Active Rac",color="#1764a0")
            ax.plot(native["x"],native["rho"],label="Active Rho",color="#d77913",ls="--")
            ax.axvline(length/2,color=".6",lw=.7)
            ax.set_title(f"L{length:g}: {label}\nRac D_half={rac['D_half']:.3f}; 50% regions={len(rac['regions']['0.5'])}",fontsize=9)
            if col==0:ax.set_ylabel("Native concentration")
            if row==3:ax.set_xlabel("Native coordinate")
            if row==0 and col==0:ax.legend(fontsize=8)
        tail=json.loads((folder/"native_tail.json").read_text())
        rows.append(dict(model=case["model"],length=length,snapshots=snapshots,full_field_tail=tail))
    fig.savefig(output/"native_profiles.png",dpi=170);fig.savefig(output/"native_profiles.pdf");plt.close(fig)
    report=dict(manifest_sha256=digest(root/"manifest.json"),source_sha256=sources,
        audit_source_sha256=digest(Path(__file__)),rows=rows,
        scope="Native profile/topology diagnostic of the longest tested continuous pulse and exact recovery; original D_half and all frozen data remain unchanged.")
    (output/"native_profiles.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(json.dumps([dict(L=r["length"],cue_end=r["snapshots"][1]["rac"],
                         recovery_end=r["snapshots"][2]["rac"]) for r in rows],indent=2))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root",type=Path,required=True);parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args();audit(args.run_root,args.output)
