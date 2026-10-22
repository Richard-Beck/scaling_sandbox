"""Read-only display/source checks for the complete-condition supplement."""
from pathlib import Path
from datetime import datetime, timezone
from html.parser import HTMLParser
import hashlib
import importlib.util
import json
import re
import subprocess
import os
for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"
import numpy as np
import pandas as pd
from pypdf import PdfReader
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.quiver import Quiver

ROOT = Path(__file__).resolve().parents[2]
TOPIC = ROOT/"figures/length_mechanisms_20261008/dev/revision_program_20261019/figure3/narrow_arc_experiment"
OUT = TOPIC/"supplement"


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


class Links(HTMLParser):
    def __init__(self):
        super().__init__(); self.paths = []

    def handle_starttag(self, tag, attrs):
        self.paths.extend(v for k, v in attrs if k in ("href", "src"))


def main():
    manifest = json.loads((OUT/"supporting/metadata/manifest.json").read_text())
    for p, digest in manifest["source_sha256"].items():
        assert sha(p) == digest, p
    assert sha(manifest["renderer"]) == manifest["renderer_sha256"]
    for rel, digest in manifest["artifact_sha256"].items():
        assert sha(OUT/rel) == digest, rel
    table = pd.read_csv(OUT/"supporting/data/displayed_conditions.csv")
    source = pd.read_csv(TOPIC/"supporting/data/narrow_arc_summary.csv")
    assert len(table) == len(source) == 645
    np.testing.assert_array_equal(table.condition_id, source.condition_id)
    scalar_checks = 0
    max_error = 0.
    for col in source.columns:
        if pd.api.types.is_numeric_dtype(source[col]):
            a, b = source[col].to_numpy(float), table[col].to_numpy(float)
            np.testing.assert_allclose(a, b, rtol=1e-13, atol=1e-13, equal_nan=True)
            scalar_checks += int(np.isfinite(a).sum())
            if np.isfinite(a).any(): max_error = max(max_error, float(np.nanmax(abs(a-b))))
        else:
            assert source[col].fillna("").tolist() == table[col].fillna("").tolist(), col
    for endpoint, col in (("withdrawal", "withdrawal_error_closed"), ("final", "retained_error_closed")):
        P = source[f"P_{endpoint}_real"].to_numpy()+1j*source[f"P_{endpoint}_imag"].to_numpy()
        angle = np.degrees(np.angle(P)); target = source.target_angle_in_P_coordinates_deg.to_numpy()
        value = 1-abs((angle-target+180)%360-180)/abs((-target+180)%360-180)
        value[(abs(P)<.1)|(source.status.to_numpy()!="complete")] = np.nan
        np.testing.assert_allclose(source[col], value, rtol=1e-12, atol=1e-12, equal_nan=True)
    text = (OUT/"index.html").read_text()
    assert not re.search(r"http-equiv\s*=\s*['\"]refresh|setInterval|setTimeout", text, re.I)
    assert not any(ord(c)<32 and c not in "\n\r\t" for c in text)
    links = Links(); links.feed(text)
    local = [p for p in links.paths if not p.startswith(("https:", "http:", "#"))]
    assert all((OUT/p.split("#")[0]).exists() for p in local)
    for model in ("wave_pinning", "goryachev", "holmes_model3", "spring"):
        for label, count in (("complete_supplement", 20), ("core_supplement", 5), ("native_profile_atlas", 15)):
            assert len(PdfReader(OUT/"supporting/vectors"/(model+"_"+label+".pdf")).pages) == count
    spec = importlib.util.spec_from_file_location("supplement_renderer", manifest["renderer"])
    renderer = importlib.util.module_from_spec(spec); spec.loader.exec_module(renderer)
    status = json.loads((TOPIC/"supporting/metadata/latest_status.json").read_text())
    runtime = {}
    for model in renderer.NAMES:
        fig, checks = renderer.vector_page(model, table, status)
        fig.canvas.draw()
        arrows = [a for ax in fig.axes for a in ax.collections if isinstance(a, Quiver)]
        expected = table[(table.model==model)&(table.status=="complete")]
        actual = []
        for q in arrows:
            actual.append([float(q.U[0]), float(q.V[0])])
            assert q.angles == "xy" and q.scale_units == "xy" and q.scale == 1
        want = []
        for size in sorted(expected.size_parameter.unique()):
            for r in expected[expected.size_parameter==size].to_dict("records"):
                for endpoint in ("withdrawal", "final"):
                    if r[f"R_{endpoint}"] >= .1:
                        want.append([.30*r[f"P_{endpoint}_real"], -.30*r[f"P_{endpoint}_imag"]])
        np.testing.assert_allclose(actual, want, rtol=0, atol=2e-14)
        for ax in fig.axes:
            origin = ax.transData.transform([0, 0])
            sx = np.linalg.norm(ax.transData.transform([1, 0])-origin)
            sy = np.linalg.norm(ax.transData.transform([0, 1])-origin)
            assert abs(sx-sy) <= 1e-10
        assert checks == manifest["vector_artist_checks"][model]
        runtime[model] = checks
        plt.close(fig)
    # Exercise actual browser code using a minimal Canvas/DOM harness; no network.
    data = re.search(r'<script id="data" type="application/json">(.*?)</script>', text, re.S).group(1)
    script = re.findall(r"<script>(.*?)</script>", text, re.S)[0]
    harness = r'''
const elementStore={};
const ctx={clearRect(){},fillText(){},fillRect(){}};
function element(id){if(!elementStore[id])elementStore[id]={value:id==='regime'?'all':id==='edge'?'10':'',textContent:'',innerHTML:'',options:[],appendChild(o){this.options.push(o);if(!this.value)this.value=o.value},addEventListener(){},getContext(){return ctx}};return elementStore[id]}
const document={getElementById:element,createElement(){return {value:'',textContent:''}}};
element('data').textContent=JSONDATA;
'''.replace("JSONDATA", json.dumps(data))
    trailer = r'''
let states=0;for(const m of Object.keys(names))for(const metric of metrics)for(const edge of ['10','50'])for(const regime of ['all','inside','straddling','outside']){el('model').value=m;el('metric').value=metric;el('edge').value=edge;el('regime').value=regime;draw();if(!el('matrixnote').textContent.includes('complete'))throw Error('no coverage');states++}console.log(JSON.stringify({selector_states:states,embedded_rows:rows.length,models:Object.keys(names).length}));
'''
    node = "/home/4473331/.conda/envs/codex/bin/node"
    js = subprocess.check_output([node], input=harness+script+trailer, text=True)
    browser = json.loads(js)
    assert browser["embedded_rows"] == 645
    result = dict(generated_at=datetime.now(timezone.utc).isoformat(), passed=True, complete=int((source.status=="complete").sum()),
                  source_hashes=len(manifest["source_sha256"]), artifact_hashes=len(manifest["artifact_sha256"]),
                  numeric_source_checks=scalar_checks, maximum_numeric_error=max_error, local_links=len(local),
                  vector_artists=runtime, browser_harness=browser, pdf_pages=160,
                  scope="Local readout/display consistency and source hashes; independent native validation remains separate.",
                  audit_source_sha256=sha(__file__), supplement_manifest_sha256=sha(OUT/"supporting/metadata/manifest.json"))
    (OUT/"supporting/metadata/delivery_check.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
