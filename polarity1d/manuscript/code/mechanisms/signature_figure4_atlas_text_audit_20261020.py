"""Audit exact frozen title/footer text extents at native PNG export DPI."""
from pathlib import Path
from datetime import datetime,timezone
import argparse
import ast
import hashlib
import importlib.util
import json
import matplotlib
matplotlib.use('Agg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);args=p.parse_args();run=args.run.resolve()
    source=run/'inputs/frozen_atlas.py';plan=json.loads((run/'manifest.json').read_text());assert sha(source)==plan['code_sha256']
    module_spec=importlib.util.spec_from_file_location('frozen_atlas',source);module=importlib.util.module_from_spec(module_spec);module_spec.loader.exec_module(module)
    plot=next(n for n in ast.parse(source.read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='plot')
    title=next(n for n in ast.walk(plot) if isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and n.func.attr=='suptitle')
    caption=next(n.value for n in ast.walk(plot) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='caption' for t in n.targets))
    footer=next(n for n in ast.walk(plot) if isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and isinstance(n.func.value,ast.Name) and n.func.value.id=='fig' and n.func.attr=='text')
    evaluate=lambda node,variables:eval(compile(ast.Expression(node),str(source),'eval'),{'__builtins__':{'len':len}},variables)
    records=[]
    for job in plan['jobs']:
        cols=min(2,len(job['cases']));env=dict(job=job,condition=job['condition'],model=job['condition']['model'],C=job['condition']['C_label'],cols=cols,NAMES=module.NAMES)
        fig=Figure(figsize=(14 if cols==2 else 8,12),dpi=170);canvas=FigureCanvasAgg(fig)
        t=fig.suptitle(evaluate(title.args[0],env),fontsize=13)
        f=fig.text(.5,.015,evaluate(caption,env),ha='center',fontsize=evaluate(next(k.value for k in footer.keywords if k.arg=='fontsize'),env))
        canvas.draw();renderer=canvas.get_renderer();width,height=canvas.get_width_height()
        for role,artist in [('title',t),('footer',f)]:
            b=artist.get_window_extent(renderer);assert b.x0>=0 and b.y0>=0 and b.x1<=width and b.y1<=height,(job['slot'],role,b.bounds)
            if role=='footer':assert b.y1<.075*height,'Footer must fit below the declared subplot layout rectangle'
        records.append(dict(slot=job['slot'],columns=cols,title_and_footer_inside_native_PNG=True,footer_below_subplot_layout=True))
    result=dict(passed=True,created_utc=datetime.now(timezone.utc).isoformat(),pages=len(records),scope='Exact title/footer AST from the frozen plot source, rendered at170DPI for every condition; source figure text bounds and footer reservation checked. Native data and axes are covered by separate source/pixel audits.',
        source_renderer_sha256=sha(source),audit_source_sha256=sha(__file__),pages_checked=records)
    (run/'results/text_layout_audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='pages_checked'}))


if __name__=='__main__':main()
