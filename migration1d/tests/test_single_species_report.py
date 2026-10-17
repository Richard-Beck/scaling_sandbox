import json

from migration1d.single_species.report import build_report
from migration1d.single_species.reporting import Figure, write_figure_manifest


def test_single_aggregate_html_embeds_all_figure_assets(tmp_path):
    analysis = tmp_path / "calibration" / "example"
    (analysis / "figure.png").parent.mkdir(parents=True)
    (analysis / "figure.png").write_bytes(b"not-an-image-needed-for-html-test")
    write_figure_manifest(
        analysis,
        title="Example <analysis>",
        introduction="A & B",
        figures=[Figure("figure.png", "Outcome", "A useful caption.")],
    )
    aggregate = build_report(tmp_path, tmp_path / "report.html")
    aggregate_html = aggregate.read_text()
    assert not (analysis / "report.html").exists()
    assert "Example &lt;analysis&gt;" in aggregate_html
    assert "calibration/example/figure.png" in aggregate_html
    manifest = json.loads((analysis / "figure_manifest.json").read_text())
    assert manifest["figures"][0]["caption"] == "A useful caption."
