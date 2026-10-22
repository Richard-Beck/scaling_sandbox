"""Render figure 2–7 editorial revisions from audited, unchanged native evidence.

Only presentation changes: grouped panel labels, concise legends, no finite-pulse
markers on uninterrupted-cue mass-median histories. No integrations or changes to
retained scientific inputs. The original renderer is a read-only data adapter.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, importlib.util, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "figures/length_mechanisms_20261008/revision_20261010/editorial"
ASSETS = OUT / "assets"
ADAPTER = ROOT / "code/mechanisms/length_revision_20261009_figures.py"
spec = importlib.util.spec_from_file_location("length_revision_adapter", ADAPTER)
a = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a)
a.FIG = ASSETS
a.DATA = OUT
current = None
chapters = {}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def panel(ax, letter, title=None):
    mapping = {
        "material": {v:v for v in "ABCD"},
        "initiation": {"A":"A", "C":"B", "E":"C", "G":"D"},
        "restoration": {"A":"A", "C":"B", "E":"C", "F":"D"},
        "sink": {"A":"A", "B":"B", "C":"C", "D":"D", "E":"E", "G":"F"},
        "interface": {"A":"A", "C":"B"},
        "growth": {"A":"A", "C":"B"},
    }
    display = mapping[current].get(letter)
    if display:
        ax.text(-.12, 1.04, display, transform=ax.transAxes,
                fontsize=13, fontweight="bold", fontfamily="DejaVu Sans")
    if title:
        ax.set_title(title, pad=10)
    ax.grid(alpha=.12, lw=.6)


def save(fig, stem, ignored_caption, methods):
    chapter = chapters[current]
    # Pair each model/dose as a facet of the same scientific panel.
    # Finite orientation-pulse markers add an unrelated outcome to Fig. 7.
    if current == "growth":
        for ax in fig.axes[:2]:
            for line in list(ax.lines):
                if len(line.get_xdata()) == 1 and line.get_marker() == "o":
                    line.remove()
    for legend in list(fig.legends):
        legend.remove()
    lengths = [Line2D([], [], color=c, lw=2, label=f"L = {L}") for L,c in a.LC.items()]
    if current == "material":
        handles = [Line2D([],[],color="#0072B2",label="Total amount / mean concentration"),
                   Line2D([],[],color="#D55E00",ls="--",label="Excess amount / peak concentration"),
                   Line2D([],[],color="#888888",marker="x",ls="none",label="Ineligible")]
    elif current == "initiation":
        handles = lengths
    elif current == "restoration":
        handles = lengths + [Line2D([],[],color=mc,lw=2,label=label) for label,model,mc in a.MODELS.values()]
    elif current == "sink":
        handles = [Line2D([],[],color=c,lw=2,label=label) for label,model,c in a.MODELS.values()]
        handles += [Line2D([],[],color="#333333",ls="-",label="Fixed gap"),
                    Line2D([],[],color="#333333",ls="--",label="Opposite end"),
                    Line2D([],[],color="#777777",label="Initial"),
                    Line2D([],[],color="#0072B2",label="10% core loss"),
                    Line2D([],[],color="#D55E00",label="30% core loss")]
    elif current == "interface":
        handles = [Line2D([],[],color="#777777",label="Initial"),
                   Line2D([],[],color="#0072B2",label="Mid-pulse / 250 after withdrawal"),
                   Line2D([],[],color="#D55E00",label="Withdrawal"),
                   Line2D([],[],color="#CC79A7",label="25 after withdrawal"),
                   Line2D([],[],color="#222222",label="10000 after withdrawal")]
    elif current == "growth":
        handles = lengths
    fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False)
    # Readable ticks; use model/dose column headings and panel letters only.
    for ax in fig.axes:
        ax.tick_params(labelsize=9)
    caption = chapter["caption"]
    outputs = {}
    for ext in ("png", "pdf", "svg"):
        path = ASSETS / f"{stem}.{ext}"
        fig.savefig(path, dpi=240, facecolor="white", bbox_inches="tight", pad_inches=.16)
        outputs[str(path.relative_to(ROOT))] = sha(path)
    plt.close(fig)
    caption_path = ASSETS / f"{stem}_caption.txt"
    caption_path.write_text(caption+"\n")
    outputs[str(caption_path.relative_to(ROOT))] = sha(caption_path)
    artifacts = {Path(k).name:v for k,v in outputs.items()}
    record = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "producer": "finished", "passed": True,
        "description": caption, "legend": caption,
        "artifact_root": ".", "artifacts": artifacts,
        "sha256": {Path(k).name:v for k,v in outputs.items()},
        "sources": dict(a.sources), "source_sha256": dict(a.sources),
        "renderer": str(Path(__file__).relative_to(ROOT)),
        "renderer_sha256": sha(__file__), "adapter_sha256": sha(ADAPTER),
        "methods": methods + " Presentation revision groups matched facets under one panel label. Retained native inputs are unchanged.",
        "audit": {"retained_native_audit": "data/runs/length_revision_20261009/independent_audit.json", "passed": True},
        "reproduction": ["OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python " + str(Path(__file__))],
    }
    if current == "growth":
        record["methods"] += " Finite-pulse gate markers are omitted because all plotted histories have uninterrupted cues. Scalar outcomes and native-field flags remain in the retained source endpoint table, which is hashed."
    (ASSETS/f"{stem}.manifest.json").write_text(json.dumps(record,indent=2)+"\n")
    chapter.update(image=str((ASSETS/f"{stem}.png").relative_to(ROOT/"figures/length_mechanisms_20261008")),
                   shape={"material":[11.7,12], "initiation":[11.7,20.5], "restoration":[11.7,14], "sink":[11.7,20.5], "interface":[11.7,12.8], "growth":[11.7,13]}[current],
                   sources=[(str((ASSETS/f"{stem}.pdf").relative_to(ROOT/"figures/length_mechanisms_20261008")), "Vector figure"),
                            (str((ASSETS/f"{stem}.manifest.json").relative_to(ROOT/"figures/length_mechanisms_20261008")), "Legend, hashes and reproduction")])


def main():
    global current, chapters
    OUT.mkdir(parents=True,exist_ok=True)
    ASSETS.mkdir(parents=True,exist_ok=True)
    chapters = json.loads((OUT/"chapter_overrides.json").read_text())
    a.panel = panel
    a.save = save
    assert a.read(a.RUN/"independent_audit.json")["passed"]
    manifest = a.read(a.RUN/"manifest.json")
    for path,digest in manifest["sources"].items():
        assert a.sha(ROOT/path) == digest, path
    trials = {}
    for slot, job in enumerate(manifest["jobs"]):
        path=a.RUN/f"workers/{slot:02d}.npz"
        meta=a.read(path.with_suffix(".json"))
        assert meta["job"] == job and a.sha(path) == meta["archive_sha256"]
        z=a.load(path)
        checkpoint=a.load(ROOT/job["case"]["checkpoint"])
        np.testing.assert_array_equal(z["originalstate"],checkpoint["state"])
        trials[slot]=z
    initial_rows=a.rows(ROOT/"data/runs/initial_polarity_length_corrected_20261008/steady_state_metrics.csv")
    for row in initial_rows:
        assert sha(ROOT/row["checkpoint"]) == row["checkpoint_sha256"]
    assert len(initial_rows) == 97
    for chapter, render in [("material",a.initial_figure),("initiation",a.chemical_figure),
                            ("restoration",a.restoration_figure),("sink",a.sink_figure),
                            ("interface",lambda:a.withdrawal_figure(trials)),
                            ("growth",lambda:a.half_mass_figure(trials))]:
        current=chapter
        a.sources={}
        a.register(ADAPTER)
        a.register(a.RUN/"independent_audit.json")
        if chapter in ("interface", "growth"):
            a.register(a.RUN/"manifest.json")
            slots = (6,7) if chapter == "interface" else range(6)
            for slot in slots:
                a.register(a.RUN/f"workers/{slot:02d}.npz")
                a.register(a.RUN/f"workers/{slot:02d}.json")
                a.register(ROOT/manifest["jobs"][slot]["case"]["checkpoint"])
        if chapter == "material":
            for row in initial_rows:
                a.register(ROOT/row["checkpoint"])
        render()
    np.savez_compressed(OUT/"plot_evidence.npz",**a.evidence)
    (OUT/"chapter_overrides.json").write_text(json.dumps(chapters,indent=2)+"\n")
    audit = {"passed": True, "generated_at": datetime.now(timezone.utc).isoformat(),
             "producer": "finished", "new_integrations": 0,
             "checks": "Retained audit, eight replay archives, source hashes and exact checkpoint identity; all 97 initial checkpoint hashes; native reaction closed forms and finite-volume mass median rerun by the audited data adapter.",
             "renderer_sha256":sha(__file__), "adapter_sha256":sha(ADAPTER),
             "artifacts": {str(p.relative_to(OUT)):sha(p) for p in ASSETS.iterdir() if p.suffix in (".png",".pdf",".svg")}}
    (OUT/"audit.json").write_text(json.dumps(audit,indent=2)+"\n")
    print(OUT)


if __name__ == "__main__":
    main()
