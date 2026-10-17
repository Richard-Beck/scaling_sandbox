"""Figure metadata consumed by the single aggregate HTML report."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path


@dataclass(frozen=True)
class Figure:
    path: str
    title: str
    caption: str


STYLE = """
:root { color-scheme: light; --ink:#202124; --muted:#5f6368; --line:#dadce0;
  --paper:#fff; --wash:#f7f8fa; --accent:#175a8a; }
* { box-sizing:border-box; }
body { margin:0; color:var(--ink); background:var(--wash);
  font:16px/1.55 system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif; }
main { max-width:1320px; margin:0 auto; padding:2.5rem 2rem 5rem; }
h1,h2,h3 { line-height:1.2; } h1 { margin-bottom:.4rem; }
h2 { margin-top:0; } h3 { margin-top:1.5rem; }
.lead { max-width:80ch; color:var(--muted); }
.contents,.guide,.section-guide { background:var(--paper); border:1px solid var(--line);
  border-radius:8px; padding:1.2rem 1.4rem; margin:1.5rem 0; }
.contents ol { columns:2; margin-bottom:0; } .contents a { color:var(--accent); }
.equations { padding:1rem 1.2rem; background:#f1f3f4; border-left:4px solid var(--accent);
  line-height:1.9; margin:1rem 0; }
.note { display:block; margin-top:1.2rem; padding:1rem 1.2rem; background:#e8f1f8;
  border-left:4px solid var(--accent); }
.result { display:block; margin:1.2rem 0; padding:1rem 1.2rem; background:#edf7f1;
  border-left:4px solid #26845b; }
.table-wrap { overflow-x:auto; }
table { width:100%; border-collapse:collapse; background:var(--paper); }
th,td { border:1px solid var(--line); padding:.65rem .75rem; text-align:left; vertical-align:top; }
thead th { background:#eef2f5; } table.summary th { width:34%; font-family:ui-monospace,monospace;
  font-size:.85rem; font-weight:600; } table.summary td { font-family:ui-monospace,monospace;
  font-size:.85rem; }
.metrics { display:grid; grid-template-columns:minmax(180px,260px) 1fr; gap:.55rem 1rem; }
.metrics dt { font-weight:700; } .metrics dd { margin:0; }
.swatch { display:inline-block; width:.8rem; height:.8rem; margin-right:.4rem; border-radius:2px; }
.swatch.low { background:#9E9E9E; } .swatch.high { background:#F0E442; }
.swatch.one { background:#0072B2; } .swatch.two { background:#009E73; }
.swatch.multi { background:#56B4E9; } .swatch.transient { background:#D55E00; }
.swatch.dynamic { background:#CC79A7; }
.profile-example { margin:1.5rem 0; padding:1rem; background:#f8fafb; border:1px solid var(--line); }
.profile-example svg { display:block; width:100%; height:auto; min-width:760px; }
.profile-example { overflow-x:auto; }
.profile-example .axis { stroke:#555; stroke-width:2; }
.profile-example .profile { fill:none; stroke:#175a8a; stroke-width:5; }
.profile-example .threshold { stroke:#D55E00; stroke-width:2; stroke-dasharray:7 5; }
.profile-example .kline { stroke:#555; stroke-width:2; stroke-dasharray:4 4; }
.profile-example .domain { fill:#56B4E9; opacity:.18; }
.profile-example text { fill:var(--ink); font:14px system-ui,sans-serif; }
.profile-example .svg-title { font-weight:700; font-size:15px; }
.figures { display:grid; grid-template-columns:repeat(auto-fit,minmax(440px,1fr)); gap:1.3rem; }
figure { margin:0; padding:1rem; background:var(--paper); border:1px solid var(--line);
  border-radius:8px; box-shadow:0 1px 3px #00000012; }
figure img { display:block; width:100%; height:auto; }
figcaption { margin-top:.8rem; color:var(--muted); }
figcaption strong,figcaption span { display:block; }
figcaption strong { color:var(--ink); margin-bottom:.2rem; }
.interpret { margin-top:.55rem; color:var(--ink); }
pre { overflow:auto; padding:1rem; background:#f1f3f4; border-radius:6px; font-size:.82rem; }
.files { display:flex; flex-wrap:wrap; gap:.6rem; padding:0; list-style:none; }
.files a,.report-link { display:inline-block; padding:.45rem .7rem; color:var(--accent);
  background:#e8f1f8; border-radius:5px; text-decoration:none; }
.analysis { margin:3.5rem 0; padding-top:1rem; border-top:3px solid var(--accent); }
.section-intro { max-width:90ch; font-size:1.05rem; }
details { margin:1.5rem 0; } details summary { cursor:pointer; font-weight:700; color:var(--accent); }
details .table-wrap { margin-top:.8rem; max-height:32rem; }
@media (max-width:620px) { main { padding:1.2rem .8rem 3rem; }
  .figures { grid-template-columns:1fr; } .contents ol { columns:1; }
  .metrics { grid-template-columns:1fr; } .metrics dd { margin-bottom:.5rem; } }
"""


def write_figure_manifest(
    output: Path,
    *,
    title: str,
    introduction: str,
    figures: list[Figure],
    data_files: tuple[str, ...] = ("summary.json",),
) -> Path:
    output.mkdir(parents=True, exist_ok=True)
    manifest = {
        "title": title,
        "introduction": introduction,
        "figures": [asdict(figure) for figure in figures],
        "data_files": list(data_files),
    }
    path = output / "figure_manifest.json"
    path.write_text(json.dumps(manifest, indent=2) + "\n")
    return path
