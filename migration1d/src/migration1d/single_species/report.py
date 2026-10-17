"""Assemble completed single-species analyses into one self-contained HTML report."""

from __future__ import annotations

import argparse
import csv
from html import escape
import json
import os
from pathlib import Path

from .reporting import STYLE


MODEL_GUIDE = """
<section id="guide" class="guide">
<h2>What is being modeled?</h2>
<p>This is a one-dimensional wave-pinning model of one conserved GTPase. “Single
species” means one molecular polarity system, not one state variable: the GTPase
switches between slowly diffusing active membrane-bound material <code>u(x,t)</code>
and rapidly diffusing inactive cytoplasmic material <code>v(x,t)</code>. Position
<code>x=0</code> is the left end of the cell strip and <code>x=L</code> is the right
end. A polarized state has active GTPase concentrated in one or more spatial
regions rather than distributed nearly uniformly.</p>
<div class="equations"><code>∂u/∂t = D<sub>u</sub> ∂²u/∂x² + J(u,v)</code><br>
<code>∂v/∂t = D<sub>v</sub> ∂²v/∂x² − J(u,v)/h</code><br>
<code>J(u,v) = δv[κ + Γu<sup>n</sup>/(K<sup>n</sup>+u<sup>n</sup>) + stimulus] − δu</code></div>
<p>The first two terms describe diffusion and exchange between compartments. The
reaction <code>J</code> activates cytoplasmic GTPase onto the membrane through a
basal rate and positive feedback, while <code>δu</code> returns active material to
the inactive pool. With no external stimulus, the conserved quantity is
<code>∫(u+h v)dx</code>.</p>

<h3>Symbols and control parameters</h3>
<div class="table-wrap"><table><thead><tr><th>Symbol</th><th>Plain-language meaning</th><th>Definition and interpretation</th></tr></thead><tbody>
<tr><td><code>u</code></td><td>Active membrane GTPase</td><td>The spatial profile used to decide whether the model is polarized.</td></tr>
<tr><td><code>v</code></td><td>Inactive cytoplasmic GTPase</td><td>The shared, faster-diffusing reservoir that supplies active regions.</td></tr>
<tr><td><code>μ</code> (mu)</td><td>Total normalized GTPase abundance</td><td><code>μ = mean(u+h v)/K</code>. Moving upward in a phase map supplies more total GTPase.</td></tr>
<tr><td><code>Λ</code> (capital lambda)</td><td>Dimensionless domain length</td><td><code>Λ = L√(δ/D<sub>u</sub>)</code>. Moving right means a longer domain relative to the distance active GTPase diffuses during one inactivation time. It can represent a larger physical cell, slower membrane diffusion, or faster inactivation—not physical length alone.</td></tr>
<tr><td><code>h</code></td><td>Cytoplasmic reservoir per membrane area</td><td><code>h = A<sub>cytoplasm</sub>/s<sub>membrane</sub></code> in the strip reduction. Larger <code>h</code> gives each unit of membrane more adjacent cytoplasmic volume. Each map panel is one <code>h</code> value.</td></tr>
<tr><td><code>ρ</code> (rho)</td><td>Diffusion ratio</td><td><code>ρ = D<sub>v</sub>/D<sub>u</sub></code>. Here it is 100, so inactive cytoplasmic material diffuses 100 times faster.</td></tr>
<tr><td><code>κ</code> (kappa)</td><td>Normalized basal activation</td><td><code>κ = k<sub>0</sub>/δ</code>. Activation that occurs without positive feedback.</td></tr>
<tr><td><code>Γ</code> (capital gamma)</td><td>Normalized feedback strength</td><td><code>Γ = γ/δ</code>. Maximum positive-feedback activation relative to inactivation.</td></tr>
<tr><td><code>n</code></td><td>Hill coefficient</td><td>Controls how switch-like the positive feedback is.</td></tr>
<tr><td><code>K</code></td><td>Half-feedback concentration</td><td>The active-GTPase level at which feedback is half maximal; used as the concentration scale. <code>K=1</code> here.</td></tr>
<tr><td><code>δ</code> (delta)</td><td>Inactivation rate</td><td>Sets the time scale; <code>δ=1 s⁻¹</code> here.</td></tr>
<tr><td><code>D<sub>u</sub>, D<sub>v</sub></code></td><td>Membrane and cytoplasmic diffusivities</td><td><code>D<sub>u</sub></code> is implied by <code>Λ</code>; <code>D<sub>v</sub>=ρD<sub>u</sub></code>.</td></tr>
<tr><td><code>ε</code> (epsilon)</td><td>Initial two-pole imbalance</td><td>In the dual-pole test, left and right pulse amplitudes are <code>A(1+ε)</code> and <code>A(1−ε)</code>. Zero is perfectly symmetric; larger values deliberately seed competition.</td></tr>
<tr><td>stimulus</td><td>Externally imposed activation cue</td><td>A dimensionless addition to the activation rate in units of <code>δ</code>. Amplitude 1 therefore adds an activation rate equal to the baseline inactivation-rate scale at the cue center.</td></tr>
</tbody></table></div>
<aside class="note"><strong>How to read every phase map.</strong> The horizontal
axis is <code>Λ</code>, the vertical axis is <code>μ</code>, and the side-by-side
panels are <code>h=1, 2, 4</code>. A pixel is one 300-second simulation. Colors in
categorical maps are outcome labels; colors in continuous maps are values read
from that figure’s color bar. Comparisons across fixed-<code>μ</code> panels change
how total material is divided between membrane and cytoplasmic volume.</aside>
<p><strong>Common numerical setup.</strong> The strip is fixed at model length
<code>L=10</code> with no-flux boundaries, so material cannot leave through either
end. It is discretized into 81 finite-volume nodes. Results are saved every 2 s,
the adaptive integrator takes steps no larger than 0.5 s, and no random noise or
mechanical motion is applied. Changing <code>Λ</code> in these maps changes
<code>D<sub>u</sub></code> while the represented strip length remains 10.</p>
</section>
"""


CLASSIFICATION_GUIDE = """
<section id="classification" class="guide">
<h2>How simulation outcomes are classified</h2>
<p>Labels are deterministic analysis categories, not informal visual judgments.
At each saved time, polarization amplitude is
<code>P=(max(u)−min(u))/mean(u)</code>. A “high domain” is a connected region where
<code>u &gt; min(u)+0.5[max(u)−min(u)]</code>. The final label uses the last profile
and the trajectory’s final 50 seconds.</p>
<div class="table-wrap"><table><thead><tr><th>Outcome shown in legends</th><th>Exact operational meaning</th><th>Biological reading</th></tr></thead><tbody>
<tr><td><span class="swatch low"></span><strong>homogeneous low</strong></td><td>Late normalized change ≤0.01, <code>P&lt;0.5</code>, and final mean <code>u&lt;K</code>.</td><td>No substantial spatial polarity; most GTPase remains on the low-active side of the feedback scale. “Low” refers to mean active membrane GTPase relative to <code>K</code>, not low total abundance.</td></tr>
<tr><td><span class="swatch high"></span><strong>homogeneous high</strong></td><td>Late normalized change ≤0.01, <code>P&lt;0.5</code>, and final mean <code>u≥K</code>.</td><td>The membrane is broadly active but not polarized. This is whole-domain activation, not a strong pole.</td></tr>
<tr><td><span class="swatch one"></span><strong>single polarized</strong></td><td><code>P≥0.5</code>, one connected high domain, and late dynamics below the change threshold.</td><td>One persistent active region: the desired monopolar wave-pinned-like outcome.</td></tr>
<tr><td><span class="swatch two"></span><strong>dual polarized</strong></td><td><code>P≥0.5</code>, two connected high domains, and late dynamics below the change threshold.</td><td>Two active regions remain at 300 s. This demonstrates persistence over the observation window, not mathematical stability.</td></tr>
<tr><td><span class="swatch multi"></span><strong>multi peak</strong></td><td><code>P≥0.5</code>, more than two connected high domains, and late dynamics below the change threshold.</td><td>Several spatial activity sites persist.</td></tr>
<tr><td><span class="swatch transient"></span><strong>transient coarsening</strong></td><td>Late normalized change exceeds 0.01 and the number of peaks falls during the final 50 s.</td><td>The pattern is still merging or selecting winners; the endpoint is not settled.</td></tr>
<tr><td><span class="swatch dynamic"></span><strong>oscillatory/dynamic</strong></td><td>Late normalized change exceeds 0.01 without a late loss of peaks.</td><td>The endpoint remains time-dependent and should not be interpreted as stationary.</td></tr>
</tbody></table></div>
<div class="profile-example"><h3>What does a connected spatial region mean?</h3>
<p>The sketches below apply the rule to one-dimensional profiles. A region is a
continuous interval along the strip—not one mesh point and not merely a local
maximum. The dashed orange line is halfway between that profile’s minimum and
maximum; each separate blue interval above it is one region. The contrast rule
<code>P≥0.5</code> must also be met, so tiny numerical ripples count as zero.</p>
<svg viewBox="0 0 1120 330" role="img" aria-label="Example active-GTPase profiles with zero, one, and two counted polarized domains">
<g transform="translate(20 20)"><text x="120" y="18" text-anchor="middle" class="svg-title">homogeneous low: 0</text><line x1="25" y1="250" x2="250" y2="250" class="axis"/><line x1="25" y1="40" x2="25" y2="250" class="axis"/><line x1="25" y1="135" x2="250" y2="135" class="kline"/><text x="245" y="128" text-anchor="end">K</text><path d="M25 210 L250 210" class="profile"/><text x="138" y="292" text-anchor="middle">position x</text></g>
<g transform="translate(300 20)"><text x="120" y="18" text-anchor="middle" class="svg-title">homogeneous high: 0</text><line x1="25" y1="250" x2="250" y2="250" class="axis"/><line x1="25" y1="40" x2="25" y2="250" class="axis"/><line x1="25" y1="135" x2="250" y2="135" class="kline"/><text x="245" y="128" text-anchor="end">K</text><path d="M25 90 L250 90" class="profile"/><text x="138" y="292" text-anchor="middle">position x</text></g>
<g transform="translate(580 20)"><text x="120" y="18" text-anchor="middle" class="svg-title">one connected region: 1</text><line x1="25" y1="250" x2="250" y2="250" class="axis"/><line x1="25" y1="40" x2="25" y2="250" class="axis"/><rect x="52" y="55" width="70" height="195" class="domain"/><line x1="25" y1="150" x2="250" y2="150" class="threshold"/><path d="M25 225 C42 225 45 200 55 145 C66 72 82 52 92 52 C105 52 116 83 125 150 C132 205 145 225 250 225" class="profile"/><text x="87" y="277" text-anchor="middle">region 1</text></g>
<g transform="translate(860 20)"><text x="120" y="18" text-anchor="middle" class="svg-title">two disconnected regions: 2</text><line x1="25" y1="250" x2="250" y2="250" class="axis"/><line x1="25" y1="40" x2="25" y2="250" class="axis"/><rect x="42" y="55" width="55" height="195" class="domain"/><rect x="178" y="55" width="55" height="195" class="domain"/><line x1="25" y1="150" x2="250" y2="150" class="threshold"/><path d="M25 225 C35 220 40 170 47 130 C55 70 66 52 76 52 C88 52 95 93 101 155 C106 215 128 225 148 225 C168 225 173 175 180 130 C188 70 198 52 208 52 C220 52 228 100 234 160 C239 210 244 225 250 225" class="profile"/><text x="70" y="277" text-anchor="middle">region 1</text><text x="205" y="277" text-anchor="middle">region 2</text></g>
</svg></div>
<p><strong>Important boundary:</strong> <code>P=0.5</code> is the operational
polarization cutoff. It means the final max–min range is half the spatial mean.
“Homogeneous” therefore includes weak gradients below this cutoff; it does not
require every mesh point to be identical.</p>

<h3>Quantities plotted in the figures</h3>
<dl class="metrics">
<dt>Polarized-domain count</dt><dd>Zero when <code>P&lt;0.5</code>; otherwise the number of separate blue intervals illustrated above. The time-dependent tests report it directly at named times.</dd>
<dt>Endpoint polarity contrast</dt><dd><code>(u<sub>right</sub>−u<sub>left</sub>)/(|u<sub>right</sub>|+|u<sub>left</sub>|)</code>. Values near −1 favor the left, near +1 favor the right, and near 0 have similar endpoints.</dd>
<dt>Late-time normalized change</dt><dd>Largest root-mean-square change between saved profiles during the last 50 s, divided by the late mean activity. Values above 0.01 trigger a dynamic label.</dd>
<dt>Conservation error</dt><dd>Largest relative numerical drift in <code>∫(u+h v)dx</code>. This is a solver-quality check; smaller is better.</dd>
<dt>Reversal time</dt><dd>Seconds after the opposing right-hand cue begins until a single right-biased domain persists for 10 s. It is only defined for conditions that had one clearly left-biased domain before that cue.</dd>
</dl>
</section>
"""


SECTION_GUIDES = {
    "local_bistability": """
<div class="section-guide"><h3>Question and protocol</h3>
<p><strong>Question:</strong> can the local reaction, even before diffusion is
considered, support both low- and high-active states? For each fixed inactive
concentration <code>v</code>, this analysis solves <code>J(u,v)=0</code>. Three
roots—two locally stable roots separated by one unstable root—define a
bistability window. This is necessary for classical wave-pinning but is not
sufficient: conservation, diffusion, domain size, and stimulation still decide
whether a PDE simulation pins a front.</p></div>""",
    "single_pulse_regime_map": """
<div class="section-guide"><h3>Question and protocol</h3>
<p><strong>Question:</strong> after a localized cue, which combinations of total
abundance, effective domain length, and geometry create sustained polarity? Each
run starts from the lowest stable homogeneous equilibrium. A Gaussian activation
cue centered at 10% of the strip, with width 5% of the strip and amplitude 1,
acts from 10 to 30 s. The cue is then removed and the simulation continues to
300 s. The grid contains 17 <code>μ</code> values × 17 <code>Λ</code> values × 3
<code>h</code> values = 867 runs. Fixed kinetics are <code>ρ=100</code>,
<code>κ=0.067</code>, <code>Γ=1</code>, and <code>n=2</code>.</p></div>""",
    "dual_pole_domain_counts": """
<div class="section-guide"><h3>Question and protocol</h3>
<p><strong>Question:</strong> if two activity sites are created, do both persist
or does one consume the shared cytoplasmic pool and win? Each run begins with
Gaussian enrichments centered at 20% and 80% of the strip. Material added to
<code>u</code> is removed from <code>v</code>, preserving mass. The initial amplitude
is 2 and width is 4% of the strip. Maps are repeated for
<code>ε=0, 0.001, 0.01, 0.05</code>; the latter three break exact numerical
symmetry. The grid contains 13 <code>μ</code> × 13 <code>Λ</code> × 3 <code>h</code>
× 4 <code>ε</code> values = 2,028 runs, each observed for 300 s.</p>
<p>The output is deliberately literal: 0, 1, 2, or 3+ qualifying connected
domains at each observation time. A transition <code>2→2→1</code> means two seeded
domains at 0 and 100 s followed by one surviving domain at 300 s. A
<code>2→2→2</code> transition means both persist over 300 s, but does not prove
infinite-time stability.</p></div>""",
    "polarity_reversal_map": """
<div class="section-guide"><h3>Question and protocol</h3>
<p><strong>Question:</strong> can an established left-hand pole be reoriented by
an equally strong cue on the right? A left cue acts from 10–30 s. A right cue
acts from 150–170 s, and the state is observed to 300 s. This is now a full
17×17×3 sweep over <code>μ</code>, <code>Λ</code>, and <code>h</code> (867 runs),
matched to the single-pulse calibration grid.</p>
<p>A condition is eligible for reversal only if it has exactly one polarized
domain and endpoint contrast ≤−0.5 immediately before the right cue. It is
“reversed” only if it ends with one domain and contrast ≥+0.5; otherwise the map
separates retention, polarity loss, multiple domains, and ambiguous direction.
Gray “not initially left polarized” cells are controls where reversal could not
be tested.</p></div>""",
}


FIGURE_GUIDES = {
    "single_pulse_phase_diagram.png": "Read the legend using the exact outcome definitions above. This is the primary map of where the model is unpolarized, monopolar, or multipolar.",
    "single_pulse_outcome_counts.png": "Each bar contains the same 289 sampled (μ, Λ) combinations. It makes the strong loss of polarization at h=4 explicit without requiring pixel-by-pixel inspection.",
    "local_bistability.png": "At a fixed horizontal-axis v, three vertical roots mean a stable low-u state and stable high-u state separated by an unstable state.",
    "local_root_count.png": "The green three-root interval is the local bistability window; outside it only one local active state exists.",
    "polarity_reversal_outcome.png": "Green is successful reversal. Blue retains the original left pole. Gray never established the required initial pole, so it is not evidence for or against reversibility.",
    "polarity_reversal_domain_counts.png": "Compare rows to see whether the opposing cue changes the number of qualifying domains, independently of their left-right direction.",
    "polarity_reversal_time.png": "Only successful reversals are shown; color gives the delay after the right-hand cue begins.",
}


def _figure_guide(filename: str) -> str:
    if filename in FIGURE_GUIDES:
        return FIGURE_GUIDES[filename]
    if filename.startswith("polarized_domain_counts_epsilon_"):
        return "Read downward within a column to follow the same (μ, Λ, h, ε) condition from initialization through 100 s to 300 s."
    return ""


def _flatten(value: object, prefix: str = "") -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = []
    if isinstance(value, dict):
        for key, child in value.items():
            label = f"{prefix}.{key}" if prefix else str(key)
            rows.extend(_flatten(child, label))
    elif isinstance(value, list):
        rows.append((prefix, ", ".join(str(item) for item in value)))
    else:
        rows.append((prefix, "not reached" if value is None else str(value)))
    return rows


def _summary_table(summary: dict[str, object]) -> str:
    rows = "".join(
        f"<tr><th>{escape(key)}</th><td>{escape(value)}</td></tr>"
        for key, value in _flatten(summary)
    )
    return f"<div class='table-wrap'><table class='summary'><tbody>{rows}</tbody></table></div>"


def _observed_result(
    analysis_name: str, summary: dict[str, object], analysis_dir: Path
) -> str:
    """Build a short, data-derived interpretation without changing the analyses."""
    if analysis_name == "local_bistability":
        interval = summary.get("three_root_inactive_interval")
        if isinstance(interval, list) and len(interval) == 2:
            text = (
                f"A three-root local bistability window was found for inactive "
                f"concentration v from {interval[0]:g} to {interval[1]:g}."
            )
        else:
            text = "No three-root local bistability window was found."
    elif analysis_name == "single_pulse_regime_map":
        counts = summary.get("state_counts", {})
        if not isinstance(counts, dict):
            return ""
        total = sum(int(count) for count in counts.values())
        pieces = [
            f"{int(count)} {str(state).replace('_', ' ')}"
            for state, count in counts.items()
        ]
        text = f"Across the {total} sampled conditions: " + ", ".join(pieces) + "."
        path = analysis_dir / "regime_map.csv"
        if path.exists():
            with path.open(newline="") as handle:
                rows = list(csv.DictReader(handle))
            h4 = [row for row in rows if float(row["h"]) == 4.0]
            if h4 and all(row["terminal_polarity_status"] == "homogeneous_low" for row in h4):
                text += " Every sampled h=4 condition was homogeneous low at fixed μ."
        text += " Counts describe this deliberately uniform parameter grid; they are not probabilities of biological behavior."
    elif analysis_name == "dual_pole_domain_counts":
        transitions = summary.get("transition_counts", {})
        if not isinstance(transitions, dict):
            return ""
        combined: dict[str, int] = {}
        for epsilon_counts in transitions.values():
            if not isinstance(epsilon_counts, dict):
                continue
            for transition, count in epsilon_counts.items():
                combined[str(transition)] = combined.get(str(transition), 0) + int(count)
        total = sum(combined.values())
        ordered = sorted(combined.items(), key=lambda item: item[1], reverse=True)
        leading = ", ".join(f"{name}: {count}" for name, count in ordered[:5])
        text = (
            f"Across {total} runs, the most frequent 0→100→300 s domain-count "
            f"transitions were {leading}. These counts state the result directly; "
            f"they do not imply infinite-time stability."
        )
    elif analysis_name == "polarity_reversal_map":
        counts = summary.get("outcome_counts", {})
        if not isinstance(counts, dict):
            return ""
        total = sum(int(count) for count in counts.values())
        ineligible = int(counts.get("not_initially_left_polarized", 0))
        eligible = total - ineligible
        reversed_count = int(counts.get("reversed", 0))
        retained = int(counts.get("retained_original_pole", 0))
        lost = int(counts.get("lost_polarity", 0))
        multiple = int(counts.get("multiple_polarized_domains", 0))
        ambiguous = int(counts.get("ambiguous", 0))
        text = (
            f"Of {total} sampled conditions, {eligible} established the required "
            f"initial left pole. Among those eligible conditions: {reversed_count} "
            f"reversed, {retained} retained the left pole, {lost} lost polarity, "
            f"{multiple} ended with multiple domains, and {ambiguous} had ambiguous "
            f"direction. The remaining {ineligible} conditions did not test reversal."
        )
    else:
        return ""
    return f"<aside class='result'><strong>Observed result.</strong> {text}</aside>"


def build_report(results_root: Path, output: Path) -> Path:
    manifests = sorted(results_root.glob("**/figure_manifest.json"))
    sections = []
    navigation = []
    for index, manifest_path in enumerate(manifests, start=1):
        analysis_dir = manifest_path.parent
        manifest = json.loads(manifest_path.read_text())
        summary_path = analysis_dir / "summary.json"
        summary = json.loads(summary_path.read_text()) if summary_path.exists() else {}
        analysis_name = str(summary.get("analysis", ""))
        section_id = f"analysis-{index}"
        navigation.append(f"<li><a href='#{section_id}'>{escape(manifest['title'])}</a></li>")
        figures = []
        for figure in manifest.get("figures", []):
            source = analysis_dir / figure["path"]
            relative = os.path.relpath(source, output.parent)
            guide = _figure_guide(figure["path"])
            guide_html = (
                f'<span class="interpret"><b>How to read it:</b> {escape(guide)}</span>'
                if guide else ""
            )
            figures.append(
                f"<figure><a href='{escape(relative)}'><img src='{escape(relative)}' "
                f"alt='{escape(figure['title'])}' loading='lazy'></a><figcaption>"
                f"<strong>{escape(figure['title'])}</strong>"
                f"<span>{escape(figure['caption'])}</span>{guide_html}</figcaption></figure>"
            )
        files = []
        for data_file in manifest.get("data_files", []):
            relative = os.path.relpath(analysis_dir / data_file, output.parent)
            files.append(f"<li><a href='{escape(relative)}'>{escape(data_file)}</a></li>")
        sections.append(
            f"<section class='analysis' id='{section_id}'><h2>{index}. {escape(manifest['title'])}</h2>"
            f"<p class='section-intro'>{escape(manifest.get('introduction', ''))}</p>"
            f"{SECTION_GUIDES.get(analysis_name, '')}"
            f"{_observed_result(analysis_name, summary, analysis_dir)}"
            f"<div class='figures'>{''.join(figures)}</div>"
            f"<details><summary>Exact run configuration and numerical results</summary>"
            f"{_summary_table(summary)}</details>"
            f"<h3>Downloadable data</h3><ul class='files'>{''.join(files)}</ul></section>"
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    document = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width"><title>Single-species wave-pinning report</title>
<style>{STYLE}</style></head><body><main><header><h1>Single-species wave-pinning report</h1>
<p class="lead">A self-contained calibration and testing report for the
membrane–cytoplasm strip model. It asks where one localized cue produces a
persistent pole, whether two poles coexist or compete, and whether an established
pole can reverse.</p></header>
<nav class="contents"><strong>Contents</strong><ol><li><a href="#guide">Model and parameter definitions</a></li>
<li><a href="#classification">Outcome and metric definitions</a></li>{''.join(navigation)}</ol></nav>
{MODEL_GUIDE}{CLASSIFICATION_GUIDE}
<section class="guide"><h2>Scope and limitations</h2><p>These results are numerical
regime maps, not a full bifurcation analysis. “Persistent” means present at the
300-second endpoint. In particular, a symmetric or long-lived dual-pole state may
still be metastable and coarsen later. The parameter-mapping work needed to
translate concentrations and geometry from a particular published model remains
a separate deferred analysis.</p></section>
{''.join(sections) if sections else '<p>No completed analyses were found.</p>'}
</main></body></html>"""
    output.write_text(document)
    return output


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    default_root = Path(__file__).resolve().parents[3] / "results" / "single_species"
    parser.add_argument("--results-root", type=Path, default=default_root)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    destination = args.output or args.results_root / "report.html"
    print(build_report(args.results_root.resolve(), destination.resolve()))


if __name__ == "__main__":
    main()
