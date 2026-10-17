"""Build a short, figure-led report from the completed, audited mechanism run.

This renderer changes presentation only. Native trials, scalar thresholds,
spatial measurements and the independent audit remain the source of evidence.
"""
from pathlib import Path
from datetime import datetime, timezone
from contextlib import contextmanager
import argparse
import csv
import html
import json
import shutil
import sys
import textwrap

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.transforms import Bbox

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'code'))
import figure_builder as f
import three_test_analysis as analysis

OUT = ROOT / 'data/runs/three_test_mechanisms_20261006'
FIG = ROOT / 'figures/three_test_mechanisms'
COLORS = ['#0072B2', '#D55E00', '#009E73', '#CC79A7']
GEOMETRIES = ['fixed_physical', 'fixed_relative', 'fixed_slope']
NAMES = {'fixed_physical': 'Fixed patch', 'fixed_relative': '20% patch',
         'fixed_slope': 'Fixed-slope gradient'}
CONTROL_NAMES = {'native': 'Original model', 'mobility_4x': 'Inactive diffusion ×4',
                 'matched_gradient_peak': 'Hold gradient peak fixed',
                 'mca_frozen_during_pulse': 'Hold attachments fixed during cue',
                 'feedback_off_during_pulse': 'Remove extension feedback during cue',
                 'mesh_401': '401-point grid'}
CONTROL_MARKERS = {'native': 'o', 'mobility_4x': 's', 'matched_gradient_peak': '^',
                   'mca_frozen_during_pulse': 'D', 'feedback_off_during_pulse': 'v',
                   'mesh_401': 'X'}


def read(path):
    return json.loads(path.read_text())


class Evidence:
    def __init__(self, out):
        self.out = out
        self.manifest = read(out / 'manifest.json')
        self.campaign = ROOT / self.manifest['campaign']
        self.status = read(out / 'results/status.json')
        audit = read(out / 'results/spatial_audit.json')
        assert self.status['complete'] and audit['passed'] and audit['complete']
        assert read(self.campaign / 'results/status.json')['complete']
        self.ends = f.read_endpoint_rows(self.campaign / 'results/endpoints.csv')
        self.records = {}
        for job in read(out / 'jobs.json'):
            path = out / 'trajectories' / f'trajectory_{job["slot"]:04d}.json'
            record = read(path)
            assert record['job'] == job
            self.records[job['slot']] = record
        self.cells = read(out / 'results/cell_evidence.json')
        with (out / 'results/trajectory_metrics.csv').open() as stream:
            self.metrics = {int(row['slot']): row for row in csv.DictReader(stream)}
        self.selected = set()
        self.verified = set()

    def anchors(self, model, geometry):
        return sorted([r for r in self.records.values()
                       if r['job']['case']['model'] == model
                       and r['job']['case']['condition'] == geometry
                       and r['job']['kind'] == 'native'
                       and any(label in r['job']['labels']
                               for label in ['orientation_pass', 'orientation_censored'])],
                      key=lambda r: r['job']['case']['length'])

    def anchor(self, model, geometry, length):
        return next(r for r in self.anchors(model, geometry)
                    if np.isclose(r['job']['case']['length'], length))

    def controls(self, model, geometry, kind):
        return sorted([r for r in self.records.values()
                       if r['job']['case']['model'] == model
                       and r['job']['case']['condition'] == geometry
                       and r['job']['kind'] == kind],
                      key=lambda r: r['job']['case']['length'])

    @contextmanager
    def history(self, record):
        slot = record['job']['slot']
        path = self.out / 'trajectories' / f'trajectory_{slot:04d}.npz'
        if slot not in self.verified:
            assert f.digest(path) == record['archive_sha256'], path
            self.verified.add(slot)
        self.selected.add(slot)
        with np.load(path) as z:
            yield z

    def number(self, record, key):
        value = self.metrics[record['job']['slot']][key]
        return float(value) if value else None


class Report:
    def __init__(self, path, stamp, data_stamp):
        self.path = path
        self.pdf = PdfPages(path)
        self.pdf.infodict()['Title'] = 'Length and polarity reversal: a reader\'s guide'
        self.pdf.infodict()['Subject'] = 'Five models, three temporary cues: starting-state burden, cue geometry and retained memory'
        self.pages = []
        self.stamp = stamp
        self.data_stamp = data_stamp
        self.text_boxes = []

    def new(self, title, kicker, panels=0):
        fig = plt.figure(figsize=(12, 9))
        fig.text(.065, .968, kicker.upper(), fontsize=8, color='#526674')
        self.text_boxes.append(fig.text(.065, .92, title, fontsize=19, weight='bold'))
        axes = []
        if panels:
            grid = fig.add_gridspec(1, panels, left=.08, right=.955,
                                    top=.81, bottom=.49, wspace=.38)
            axes = [fig.add_subplot(grid[0, i]) for i in range(panels)]
        return fig, axes

    def prose(self, fig, paragraphs, top=.40, width=145, size=10):
        y = top
        for heading, body in paragraphs:
            if heading:
                self.text_boxes.append(fig.text(.065, y, heading, fontsize=size,
                                               weight='bold', va='top'))
                y -= .026
            lines = textwrap.wrap(body, width=width)
            self.text_boxes.append(fig.text(.065, y, '\n'.join(lines),
                                           fontsize=size, va='top', linespacing=1.4))
            y -= .022 * len(lines) + .018
        assert y > .035, (paragraphs, y)

    def save(self, fig, title, paragraphs, slug, has_plot=True):
        number = len(self.pages) + 1
        for ax in fig.axes:
            ax.spines[['top', 'right']].set_visible(False)
            ax.tick_params(labelsize=9)
        fig.text(.065, .019, 'Length and polarity reversal · five models, three temporary cues',
                 fontsize=8, color='#526674')
        fig.text(.95, .019, str(number), fontsize=9, ha='right', color='#526674')
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        for artist in self.text_boxes:
            bbox = artist.get_window_extent(renderer=renderer)
            assert bbox.y0 > .035 * fig.bbox.height, ('Text overflow', title)
            assert bbox.x1 < .978 * fig.bbox.width, ('Text too wide', title)
        self.text_boxes = []
        image_name = f'reader_{number:02d}_{slug}.png' if has_plot else None
        if has_plot:
            fig.savefig(FIG / image_name, dpi=170,
                        bbox_inches=Bbox.from_bounds(0, 3.87, 12, 3.85))
        self.pdf.savefig(fig)
        plt.close(fig)
        self.pages.append(dict(number=number, title=title, paragraphs=paragraphs,
                               image=image_name))

    def overview(self):
        title = 'How does cell length change the cue needed to reverse polarity?'
        fig = plt.figure(figsize=(14.5, 17))
        fig.text(.055, .967, title, fontsize=19, weight='bold')
        fig.text(.055, .94, 'Start with the measured curves; the following pages explain their shapes.',
                 fontsize=12)
        ax = fig.add_axes([.025, .085, .95, .845])
        ax.imshow(plt.imread(ROOT / 'figures/supplement_width_response.png'))
        ax.axis('off')
        caption = ('Read down the rows to compare models and across the columns to compare cues. '
                   'The curves report cue duration, not the time needed to finish reversing. '
                   'Blue measures persistent right orientation; orange adds recovery of the original polarity strength.')
        fig.text(.055, .059, '\n'.join(textwrap.wrap(caption, width=149)),
                 fontsize=10.5, va='top', linespacing=1.4)
        fig.text(.95, .014, '1', fontsize=10, ha='right')
        self.pdf.savefig(fig)
        plt.close(fig)
        self.pages.append(dict(number=1, title=title, image='../supplement_width_response.png',
                               paragraphs=[('', caption)]))

    def finish(self):
        assert len(self.pages) <= 25
        self.pdf.close()
        links = ''.join(f'<li><a href="#page-{p["number"]}">{p["title"]}</a></li>'
                        for p in self.pages)
        sections = []
        for p in self.pages:
            picture = f'<img src="{p["image"]}" alt="Selected evidence for {html.escape(p["title"])}">' if p['image'] else ''
            paragraphs = ''.join((f'<h3>{html.escape(h)}</h3>' if h else '') +
                                 f'<p>{html.escape(t)}</p>' for h, t in p['paragraphs'])
            sections.append(f'<section id="page-{p["number"]}"><h2>{p["number"]}. {html.escape(p["title"])}</h2>{picture}{paragraphs}</section>')
        diagnostics = []
        for model in f.MODELS:
            for geometry in GEOMETRIES:
                stem = f'{model}_{geometry}'
                label = f'{f.LABELS[model]} · {NAMES[geometry]}'
                if model == 'spring' and geometry == 'fixed_physical':
                    label += ' (width 0.2)'
                images = ''.join(f'<img loading="lazy" src="{stem}_{suffix}.png" alt="Detailed {suffix} for {html.escape(label)}">'
                                 for suffix in ['evidence', 'fields', 'mechanics']
                                 if (FIG / f'{stem}_{suffix}.png').exists())
                diagnostics.append(f'<details><summary>{label}</summary>{images}</details>')
        (FIG / 'index.html').write_text(f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Length and polarity reversal: a reader's guide</title>
<style>body{{max-width:1140px;margin:36px auto;padding:0 28px;font:17px/1.65 system-ui;color:#202c35}}
h1,h2,h3{{line-height:1.25}}h1{{font-size:32px}}h2{{font-size:25px}}h3{{font-size:19px}}
img{{display:block;max-width:100%;height:auto;margin:20px auto}}section{{border-top:1px solid #dbe2e6;margin-top:40px;padding-top:20px}}
a{{color:#176487}}.status{{font-size:14px;color:#526674}}details{{margin:14px 0}}summary{{cursor:pointer;font-weight:600}}</style></head>
<body><h1>Length and polarity reversal: a reader's guide</h1>
<p>Three questions connect the five models: what a longer cell must overcome, how the cue helps, and what retains the old direction.</p>
<p><a href="three_test_mechanisms.pdf">Read the {len(self.pages)}-page PDF</a></p>
<p class="status">Evidence generated {self.data_stamp}; simulation producer finished. Report revised {self.stamp}.</p>
<details><summary>Contents</summary><ol>{links}</ol></details>
{''.join(sections)}
<section><h2>Detailed evidence</h2><p>The selected figures above draw on all fifteen model/test combinations. Full diagnostic figures remain available below, alongside the
<a href="../../data/runs/three_test_mechanisms_20261006/results/trajectory_metrics.csv">measurements</a>,
<a href="../../data/runs/three_test_mechanisms_20261006/results/cell_evidence.json">per-test evidence</a>, and
<a href="../../data/runs/three_test_mechanisms_20261006/results/spatial_audit.json">independent spatial audit</a>.</p>
{''.join(diagnostics)}</section></body></html>''')


def profile(ax, evidence, record, title, normalized=True):
    with evidence.history(record) as z:
        length = record['job']['case']['length']
        x = z['pulse_x'] / length if normalized else z['pulse_x']
        for index, name, style, color in [(0, 'Before cue', ':', '#667985'),
                                         (-1, 'At cue withdrawal', '--', '#D55E00')]:
            ax.plot(x, z['pulse_active'][index], ls=style, color=color, lw=2, label=name)
        ax.plot(x, z['recovery_active'][-1], color='#0072B2', lw=2,
                label='After recovery')
    outcome = 'reverses' if record['metrics']['orientation_success'] else 'does not reverse'
    applied_duration = record['job']['duration']
    duration_label = f'{applied_duration:,.0f}' if applied_duration >= 100 else f'{applied_duration:.3g}'
    ax.set(title=f'{title}\nCue duration {duration_label}; {outcome}',
           xlabel='Position along cell (fraction)', ylabel='Active concentration')
    ax.legend(frameon=False, fontsize=8, loc='best')


def controls(ax, evidence, model, geometry, kinds, title):
    for i, kind in enumerate(kinds):
        records = evidence.anchors(model, geometry) if kind == 'native' else evidence.controls(model, geometry, kind)
        color = COLORS[i]
        xs = [r['job']['case']['length'] for r in records]
        ys = [r['metrics']['final_orientation'] for r in records]
        ax.plot(xs, ys, color=color, lw=1, alpha=.65)
        for r, x, y in zip(records, xs, ys):
            ax.plot(x, y, marker=CONTROL_MARKERS[kind], color=color,
                    mfc=color if r['metrics']['orientation_success'] else 'white', ms=7, mew=1.5)
        ax.plot([], [], marker=CONTROL_MARKERS[kind], color=color, ls='-',
                label=CONTROL_NAMES[kind])
    ax.axhline(.05, ls=':', lw=1, color='#667985')
    ax.axhline(0, lw=.6, color='#c4cdd2')
    ax.set(title=title, xlabel='Resting length L0' if model == 'spring' else 'Cell length',
           ylabel='Polarity after recovery')
    ax.legend(frameon=False, fontsize=8, loc='best')


def duration(ax, evidence, model, geometry, title):
    for gate, color, marker, _ in f.GATES:
        rows = sorted([r for r in evidence.ends if r['model'] == model and r['condition'] == geometry
                       and r['definition'] == gate], key=lambda r: r['length'])
        label = 'Right orientation' if gate == 'orientation' else 'Original strength'
        ax.plot([r['length'] for r in rows], [r['endpoint'] if r['status'] == 'success' else np.nan for r in rows],
                color=color, ls='-' if gate == 'orientation' else '--', label=label)
        for r in rows:
            if r['status'] == 'success':
                ax.plot(r['length'], r['endpoint'], marker=marker, color=color, ms=5,
                        mfc=color if r['selected_field_tail_pass'] else 'white')
            elif r['status'] == 'censored':
                ax.plot(r['length'], r['maximum_tested_duration'], '^', color=color,
                        mfc='white', ms=6)
    ax.set(title=title, xlabel='Cell length', ylabel='Shortest passing cue duration', yscale='log')
    ax.legend(frameon=False, fontsize=8)


def render(out=OUT):
    """Retained entry point now rebuilds the current corrected mechanism reader."""
    from length_mechanisms_report_20261008 import main
    main()


def render_historical(out=OUT):
    """Archived 6 October design; contains deprecated centered deBelly evidence."""
    FIG.mkdir(exist_ok=True)
    evidence = Evidence(out)
    pdf_path = FIG / 'three_test_mechanisms.pdf'
    old_provenance = read(out / 'results/report_provenance.json')
    if pdf_path.exists() and old_provenance.get('report_style') != 'figure_led_reader_guide':
        archive = FIG / 'three_test_mechanisms_detailed_20261006.pdf'
        if not archive.exists():
            shutil.copy2(pdf_path, archive)
    stamp = datetime.now(timezone.utc).isoformat()
    report = Report(pdf_path, stamp, evidence.status['generated_at'])
    report.overview()

    title = 'One starting cell, three temporary challenges'
    fig, axes = report.new(title, 'Protocol', 3)
    q = np.linspace(0, 1, 401)
    for i, ax in enumerate(axes):
        for length, color, ls in [(15, COLORS[0], '-'), (30, COLORS[1], '--')]:
            signal = f.cosine(q, 3 / length if i == 0 else .2, end='right') if i < 2 else (length / 15) * q
            ax.plot(q, signal, color=color, ls=ls, lw=2, label=f'Length {length}')
        ax.set(title=['A  Fixed width', 'B  Fixed fraction (20%)', 'C  Fixed gradient slope'][i],
               xlabel='Position along cell (fraction)', ylabel='Added cue / reference peak')
        ax.legend(frameon=False, fontsize=9)
    paragraphs = [('', 'At every model and length, all three tests start from the same left-polarized cell, settled without added cue. The preparation protocol uses a left patch covering 20% of the cell; Spring can also polarize without a preparation pulse. Each independent trial applies a right-favouring cue for a chosen duration, removes it, and follows recovery for 10,000 model time units.'),
                  ('What changes between tests', 'A fixed patch covers 3 model units (0.2 in Spring). The relative patch covers 20% of cell length. Both keep the same peak as length changes. A gradient keeps its physical slope, so its peak grows in a longer cell (C). Goryachev uses a patch peak of 6.4, chosen once and held fixed across lengths; the other models retain their nominal peaks.'),
                  ('What the curves measure', 'Polarity weights where active material lies: negative points left and positive points right. Blue asks for persistent right orientation; orange adds a magnitude within 1% of the starting magnitude. Duration excludes preparation and recovery. On duration curves, open circles/squares fail the native-field stability check; triangles mark no success in the tested range. On control plots, filled symbols reverse and hollow symbols fail. Time and concentration use each model\'s native units.')]
    report.prose(fig, paragraphs)
    report.save(fig, title, paragraphs, 'protocol')

    title = 'Meet the five models'
    fig, _ = report.new(title, 'The parts needed to read the figures')
    paragraphs = [
        ('Wave-pinning: a broad chemical domain', 'A slowly spreading active molecule reinforces its own activation and exchanges material with a faster inactive pool. Their conserved total allows one broad active region to persist after the cue is removed. A new active region competes with the old one for that shared material.'),
        ('Goryachev: a concentrated chemical pole', 'Active Cdc42 promotes further activation and recruits material from a shared inactive pool. In the regime studied here, this feedback builds a narrow pole whose height rises markedly with cell length. The cue can activate only the inactive material available at the opposite end.'),
        ('Holmes 3: coupled Rac and Rho domains', 'Rac and Rho inhibit one another, and each exchanges material with its own inactive pool. Polarization is a spatial contest between coupled chemical domains. A rise in right-end Rac activity does not necessarily replace the old direction.'),
        ('deBelly: chemicals coupled to local attachments', 'Rac and Rho interact with membrane–cortex attachments (MCA), which link the membrane to the underlying cortex, and with mechanical tension. The attachment pattern can retain the old direction even while Rac responds to the new cue. Bound and unbound attachments can exchange and redistribute.'),
        ('Spring: chemicals coupled to whole-cell extension', 'A polarity circuit changes the cell\'s moving boundaries. Extension increases deactivation throughout the cell, while motion also dilutes concentrations. The mechanical feedback is therefore global, unlike deBelly\'s spatial attachment pattern.'),
        ('Three questions organize the evidence', 'What must a longer cell overcome? Compare the starting profiles on page 4. How does cue geometry help? Follow coverage, available material and gradient strength on pages 5–10. What retains the old direction? Compare attachment memory and extension feedback on pages 11–14, then separate cue duration from completed reorientation on page 15.')]
    report.prose(fig, paragraphs, top=.84, width=130, size=11)
    report.save(fig, title, paragraphs, 'model_introductions', has_plot=False)

    title = 'Longer cells add different burdens to the starting pole'
    fig, axes = report.new(title, '1 · What must a longer cell overcome?', 3)
    for ax, model, lengths, letter in zip(axes, ['wave_pinning', 'holmes_model3', 'goryachev'],
                                        [[5, 21, 39], [10, 28, 50], [3, 19, 39]], 'ABC'):
        for length, color in zip(lengths, COLORS):
            record = evidence.anchor(model, 'fixed_relative', length)
            with evidence.history(record) as z:
                ax.plot(z['pulse_x'], z['pulse_active'][0], color=color, lw=2,
                        label=f'Length {length}')
        ax.set(title=f'{letter}  {f.LABELS[model]}', xlabel='Position (model units)',
               ylabel='Initial active concentration')
        ax.legend(frameon=False, fontsize=9)
    paragraphs = [
        ('Read width and height separately', 'In Wave-pinning (A), the left active region widens from about 1.5 units at length 5 to 11.9 at length 39. In Holmes (B), it widens from about 4 units at length 10 to 21 at length 50. The height changes much less. A fixed-width cue therefore stimulates a decreasing fraction of a cell whose existing active region grows.'),
        ('Goryachev changes in another way', 'Its pole stays about 0.7–0.8 units wide between lengths 3 and 39, but its peak rises from about 6 to 117 (C). Length adds concentrated material rather than a proportionally wider active region. The following pages show why that concentration also weakens the cue\'s initial action at the other end.'),
        ('The same preparation underlies all three tests', 'These are actual starting profiles from one shared preparation per length. Width is measured at half the initial pole contrast above the mid-cell background. Each starting profile is used unchanged for the fixed patch, the 20% patch and the gradient.')]
    report.prose(fig, paragraphs)
    report.save(fig, title, paragraphs, 'starting_profiles')

    for model, length, label in [('wave_pinning', 39, 'Wave-pinning'), ('holmes_model3', 50, 'Holmes 3')]:
        title = f'{label}: a wider patch changes the long-cell outcome'
        fig, axes = report.new(title, '2 · How does cue geometry help?', 3)
        profile(axes[0], evidence, evidence.anchor(model, 'fixed_physical', length), 'A  Fixed-width patch')
        profile(axes[1], evidence, evidence.anchor(model, 'fixed_relative', length), 'B  Patch covering 20%')
        controls(axes[2], evidence, model, 'fixed_physical', ['native', 'mobility_4x'],
                 'C  Fixed patch: change transport')
        if model == 'wave_pinning':
            paragraphs = [
                ('The geometry comparison', 'At length 39, the fixed patch still fails after a 5,120-unit cue (A), whereas the 20% patch leaves persistent right orientation after a 204-unit cue (B). The original peak is the same in both tests. The broader patch stimulates more of the cell and helps the new active region compete for shared material.'),
                ('Transport provides a second route', 'Fourfold inactive diffusion rescues the fixed-patch case at the unchanged 5,120-unit duration (C). This intervention changes movement of the shared material, while the cue itself stays fixed. Coverage and redistribution therefore make distinct contributions to the long-cell response.'),
                ('Read the recovered shapes', 'Both blue profiles retain activity at both ends. A small asymmetry separates the fixed-patch polarity of about 0.041 from the 20% value of about 0.051, just above the 0.05 orientation threshold. This explains why similar-looking profiles receive different classifications. Filled control symbols meet the complete right-orientation criterion; hollow symbols fail it.')]
        else:
            paragraphs = [
                ('Coverage is especially important here', 'The fixed patch reverses only 10 of 23 eligible lengths. At length 50 it fails even after a 5,120-unit cue (A); the 20% patch reverses this cell after a 1,020-unit cue and reverses all 23 eligible lengths (B). This broader stimulation helps the right Rac domain compete with the left Rac/Rho organization.'),
                ('Faster transport alone does not rescue these patch failures', 'At sampled lengths 28 and 50, fourfold inactive diffusion still fails at the original 5,120-unit limit (C). The result differs from Wave-pinning. Increasing spatial coverage is more effective than this particular transport change in the long Holmes cells.'),
                ('Direction and strength remain visibly different', 'The recovered length-50 profile retains two poles with a small rightward asymmetry, just above the 0.05 orientation threshold. It does not recover the original polarity magnitude within the tested range. That explains the separated blue and orange curves in the opening overview. Filled control symbols reverse; hollow symbols fail.')]
        report.prose(fig, paragraphs)
        report.save(fig, title, paragraphs, f'{model}_patches')

    title = 'Goryachev: less inactive material is available to the cue'
    fig, axes = report.new(title, '2 · How does cue geometry help?', 3)
    anchors = evidence.anchors('goryachev', 'fixed_relative')
    xs = [r['job']['case']['length'] for r in anchors]
    axes[0].plot(xs, [evidence.number(r, 'initial_inactive_right') for r in anchors], 'o-', color=COLORS[0], lw=2)
    axes[0].set(title='A  Available material at the right end', xlabel='Cell length',
                ylabel='Initial inactive concentration')
    for geometry, color in zip(GEOMETRIES, COLORS):
        records = evidence.anchors('goryachev', geometry)
        values = [evidence.number(r, 'initial_right_cue_activation_source') for r in records]
        axes[1].plot([r['job']['case']['length'] for r in records], np.array(values) / values[0],
                     'o-', color=color, lw=2, label=NAMES[geometry])
    axes[1].set(title='B  Initial action of the cue', xlabel='Cell length',
                ylabel='Right-end activation / length-3 value', yscale='log')
    axes[1].legend(frameon=False, fontsize=9)
    for geometry, color in zip(GEOMETRIES, COLORS):
        rows = sorted([r for r in evidence.ends if r['model'] == 'goryachev'
                       and r['condition'] == geometry and r['definition'] == 'orientation'
                       and r['status'] == 'success'], key=lambda r: r['length'])
        axes[2].plot([r['length'] for r in rows], [r['endpoint'] for r in rows],
                     color=color, lw=1.5, label=NAMES[geometry])
        for r in rows:
            axes[2].plot(r['length'], r['endpoint'], 'o', color=color, ms=4,
                         mfc=color if r['selected_field_tail_pass'] else 'white')
    axes[2].set(title='C  Persistent right orientation', xlabel='Cell length',
                ylabel='Shortest passing cue duration', yscale='log')
    axes[2].legend(frameon=False, fontsize=8)
    paragraphs = [
        ('A tall pole depletes the opposite end', 'As the initial pole concentrates material, inactive abundance at the right end falls from 0.88 at length 3 to 0.088 at length 39 (A). The native activation term multiplies cue strength by this inactive concentration. Thus the same patch peak produces about ten times less initial right-end activation in the long cell (B).'),
        ('A larger footprint does not remove that local depletion', 'Both patches reverse all 21 eligible lengths at the chosen, length-independent peak of 6.4 (C). Nevertheless, the fixed-patch duration rises from about 34 to 3,540, and the 20% duration from about 11 to 2,680. Broadening stimulation helps, but it does not restore the substrate available at the right endpoint.'),
        ('The gradient partly compensates', 'Because its peak grows with length, the gradient offsets much of the decrease in inactive abundance: the initial right-end source is about 1.3 times its length-3 value at length 39. B compares local activation rates, not total stimulation over the cell. The next pages test whether this input-strength contribution explains the whole response.')]
    report.prose(fig, paragraphs)
    report.save(fig, title, paragraphs, 'goryachev_substrate')

    title = 'The gradient curves show two different kinds of length response'
    fig, axes = report.new(title, '2 · Put the duration curves beside the mechanism tests', 3)
    for ax, model, letter in zip(axes, ['wave_pinning', 'holmes_model3', 'goryachev'], 'ABC'):
        duration(ax, evidence, model, 'fixed_slope', f'{letter}  {f.LABELS[model]}')
    paragraphs = [
        ('Wave-pinning and Holmes have an intermediate minimum', 'Wave-pinning requires about 31 units at length 5, reaches a minimum near length 13, and rises to 39 units at length 39 (A). Holmes likewise reaches a minimum near length 28 before increasing again (B). Both start with broad active regions that grow with length, while the gradient becomes stronger at the opposite end.'),
        ('Goryachev continues to slow with length', 'Its orientation duration rises from 104 units at length 3 to 2,020 at length 39 (C), despite the growing gradient peak. The narrow, taller incumbent and the depleted inactive substrate change the balance differently. At the longest length, recovering original strength also requires a longer pulse than obtaining a small persistent rightward bias.'),
        ('Now test the competing explanations', 'On the next two pages, a fixed-peak gradient tests the contribution from stronger input. Faster inactive diffusion tests the contribution from redistribution. The blue and orange curves are separate readouts; open circles or squares retain scalar successes whose native fields fail the stability check.')]
    report.prose(fig, paragraphs)
    report.save(fig, title, paragraphs, 'gradient_duration_curves')

    title = 'Wave-pinning and Holmes: gradient strength and transport both help'
    fig, axes = report.new(title, '2 · How does cue geometry help?', 2)
    for ax, model, letter in zip(axes, ['wave_pinning', 'holmes_model3'], 'AB'):
        controls(ax, evidence, model, 'fixed_slope',
                 ['native', 'matched_gradient_peak', 'mobility_4x'], f'{letter}  {f.LABELS[model]}')
    paragraphs = [
        ('Hold the gradient peak fixed', 'The original gradient becomes taller as length increases. The orange control adjusts its slope so that its peak stays at the length-15 value (Spring uses reference length 1). At the original selected durations, the longer Wave-pinning cells (lengths 21 and 39) and Holmes cells (28 and 50) then retain the old direction. The growing input helps the original-model long-cell reversals.'),
        ('Let inactive material move faster', 'The green control increases inactive diffusion fourfold, during both the cue and recovery. Reversal is retained at every sampled gradient duration, with substantially stronger recovered polarity in the longer cells. In Holmes this changes both inactive Rac and inactive Rho transport.'),
        ('Why a minimum can appear in the duration curve', 'A longer cell receives stronger gradient initiation but also contains a larger original domain and a larger pool to redistribute. The opposite effects help explain the minimum near length 13 in Wave-pinning and length 28 in Holmes. These panels separate the contributions: the cue can initiate a competitor quickly while competition and recovery still determine its final strength. Each control uses the original selected cue duration.')]
    report.prose(fig, paragraphs)
    report.save(fig, title, paragraphs, 'gradient_strength_transport')

    title = 'Goryachev: a taller gradient still leaves slow competition'
    fig, axes = report.new(title, '2 · How does cue geometry help?', 2)
    for record, color in zip(evidence.anchors('goryachev', 'fixed_slope'), COLORS):
        with evidence.history(record) as z:
            _, traces = analysis.measurements(z, record['job'])
            length = record['job']['case']['length']
            axes[0].plot(traces['time'], traces['old_fraction'], color=color, lw=2,
                         label=f'Length {length:g}')
            cut = traces['cut']
            axes[0].plot(traces['time'][cut], traces['old_fraction'][cut], 's', color=color, ms=6)
    axes[0].set(title='A  Removal of old-pole activity', xlabel='Time since cue onset',
                ylabel='Old-pole activity / initial amount', xscale='symlog')
    axes[0].axhline(.5, ls=':', lw=1, color='#667985')
    axes[0].legend(frameon=False, fontsize=9)
    controls(axes[1], evidence, 'goryachev', 'fixed_slope',
             ['native', 'matched_gradient_peak'], 'B  Change gradient strength')
    paragraphs = [
        ('The redistribution clock lengthens', 'In A, old-pole activity reaches half its initial amount after about 53 units at length 3, 223 at length 19 and 1,108 at length 39. Squares mark cue withdrawal. The gradient duration likewise rises from 104 to 2,020 units across these lengths, even though the growing gradient partly compensates for inactive-substrate depletion.'),
        ('The growing peak still contributes', 'Holding the gradient peak fixed prevents the middle and long sampled cells from reversing at their original durations (B). Stronger initial activation is therefore helpful, but it does not eliminate the subsequent competition for concentrated material.'),
        ('A narrow pole is not a constant-speed erosion problem', 'The starting width barely changes with length, and measured retreat speeds change markedly. A width divided by one universal edge speed would miss both the growing pole height and the depleted inactive pool. At length 39 the selected orientation pulse leaves only a small positive polarity; reaching the original reversed magnitude requires a different, longer selected pulse.')]
    report.prose(fig, paragraphs)
    report.save(fig, title, paragraphs, 'goryachev_competition')

    title = 'deBelly: changing Rac is not enough to rewrite the attachments'
    fig, axes = report.new(title, '3 · What retains the old direction?', 3)
    record = evidence.anchor('debelly', 'fixed_slope', 39)
    with evidence.history(record) as z:
        q = z['pulse_x'] / 39
        bound = z['pulse_field_mca_bound']
        moment = (bound @ (2 * q - 1)) / bound.sum(axis=1)
        axes[0].plot(z['pulse_time'] / record['job']['duration'], z['pulse_orientation'],
                     color=COLORS[0], lw=2, label='Rac polarity')
        axes[0].plot(z['pulse_time'] / record['job']['duration'], moment,
                     color=COLORS[1], lw=2, ls='--', label='Bound-attachment polarity')
        axes[1].plot(q, bound[0], ':', color='#667985', lw=2, label='Before cue')
        axes[1].plot(q, bound[-1], color=COLORS[1], lw=2, label='At withdrawal')
    axes[0].set(title='A  Chemical and attachment responses', xlabel='Fraction of cue duration',
                ylabel='Spatial polarity')
    axes[1].set(title='B  Bound attachments at length 39', xlabel='Position along cell (fraction)',
                ylabel='Bound MCA concentration')
    axes[0].legend(frameon=False, fontsize=8)
    axes[1].legend(frameon=False, fontsize=8)
    for i, geometry in enumerate(GEOMETRIES):
        native = evidence.anchors('debelly', geometry)
        frozen = evidence.controls('debelly', geometry, 'mca_frozen_during_pulse')
        offset = [-.12, 0, .12][i]
        for records, marker, filled in [(native, 'o', True), (frozen, 'D', False)]:
            for r in records:
                axes[2].plot(r['job']['case']['length'] + offset, r['metrics']['final_orientation'],
                             marker=marker, color=COLORS[i], mfc=COLORS[i] if filled else 'white', ms=6)
        axes[2].plot([], [], color=COLORS[i], lw=2, label=NAMES[geometry])
    axes[2].axhline(.05, ls=':', lw=1, color='#667985')
    axes[2].axhline(0, lw=.6, color='#c4cdd2')
    axes[2].set(title='C  Hold attachments fixed during cue', xlabel='Cell length', ylabel='Polarity after recovery')
    axes[2].legend(frameon=False, fontsize=8)
    paragraphs = [
        ('A second state must change', 'At length 39 under the gradient, Rac rapidly shifts rightward while the bound-attachment pattern changes on a different trajectory (A–B). MCA denotes membrane–cortex attachments, the links to the underlying cortex. Their distribution records an additional state; attachment polarity describes where linkers lie, so its sign need not match Rac polarity.'),
        ('Test that record directly', 'The control in C holds both bound and unbound attachment fields fixed during the cue, then restores their original dynamics for recovery. All twelve sampled cases that reverse in the original model (filled circles) fail to retain right orientation under this intervention (hollow diamonds), across all three cue geometries.'),
        ('What to measure', 'A rightward Rac response during stimulation does not by itself identify a lasting reversal. Rac and attachment profiles at withdrawal must be read together, then followed through recovery. The next page compares recovery after changing attachment mobility and gradient strength.')]
    report.prose(fig, paragraphs)
    report.save(fig, title, paragraphs, 'debelly_attachments')

    title = 'deBelly: attachment transport can erase the retained polarity'
    fig, axes = report.new(title, '3 · What retains the old direction?', 2)
    kinds = ['native', 'matched_gradient_peak', 'mobility_4x', 'mca_frozen_during_pulse']
    for kind, color in zip(kinds, COLORS):
        record = evidence.anchor('debelly', 'fixed_slope', 39) if kind == 'native' else next(
            r for r in evidence.controls('debelly', 'fixed_slope', kind)
            if r['job']['case']['length'] == 39)
        with evidence.history(record) as z:
            name = 'Unbound-attachment diffusion ×4' if kind == 'mobility_4x' else CONTROL_NAMES[kind]
            axes[0].plot(z['recovery_time'], z['recovery_orientation'], lw=2, color=color, label=name)
    axes[0].set(title='A  Length 39: same gradient duration', xlabel='Time after cue withdrawal',
                ylabel='Rac polarity', xscale='symlog')
    axes[0].axhline(.05, ls=':', lw=1, color='#667985')
    axes[0].axhline(0, lw=.6, color='#c4cdd2')
    axes[0].legend(frameon=False, fontsize=8)
    controls(axes[1], evidence, 'debelly', 'fixed_slope',
             ['native', 'matched_gradient_peak'], 'B  Hold gradient peak fixed')
    paragraphs = [
        ('Different interventions lead to different retained states', 'At the unchanged length-39 gradient pulse, the original model recovers right polarity of about 0.069. Holding the gradient peak fixed also leaves right polarity (about 0.072). Holding attachments fixed during the cue instead returns the cell left (about −0.073), while faster unbound-attachment diffusion leaves it nearly unpolarized (about 0.0035). A separates those outcomes in time.'),
        ('The gradient peak is not the sampled bottleneck here', 'All four sampled fixed-peak controls still reverse at their original durations (B), unlike Wave-pinning, Holmes and Goryachev. Attachment redistribution offers a more direct explanation of why a transient chemical response can be retained or lost.'),
        ('Transport sets memory as well as speed', 'Increasing unbound-attachment diffusion often removes persistent polarization rather than simply producing a faster reversed state. This is why the mobility intervention must be interpreted through the recovered profile and polarity, not as an automatic acceleration of the same process.')]
    report.prose(fig, paragraphs)
    report.save(fig, title, paragraphs, 'debelly_recovery')

    title = 'deBelly: the duration minimum also depends on the readout'
    fig, axes = report.new(title, '3 · Read the curve alongside the recovered state', 3)
    for ax, geometry, letter in zip(axes, GEOMETRIES, 'ABC'):
        duration(ax, evidence, 'debelly', geometry, f'{letter}  {NAMES[geometry]}')
    paragraphs = [
        ('The same model has three differently placed minima', 'All three cues reverse all twelve eligible lengths. The orientation-duration minimum lies near length 19 for the fixed patch, 21 for the 20% patch, and 23 for the gradient. Broader stimulation lowers the long-cell durations, but the attachment state still has to be changed and the recovered polarity has to satisfy the final observation criterion.'),
        ('The short-cell arm is close to the polarity boundary', 'At length 18 the initial polarity magnitude is only just above 0.05. The response must therefore relax into a narrow region that both passes the right-orientation requirement and, for the orange curve, matches the initial magnitude. The long pulse at this boundary cannot be read simply as a long time to create right-end Rac activity.'),
        ('Keep the relevant numerical issue local', 'Open circles or squares mark successful scalar readouts whose full profiles are still changing. At length 39, the selected fixed-patch pulse reverses on the original grid but fails on a 401-point grid. Its exact duration is therefore resolution sensitive. That issue belongs to A, rather than to every explanation in the report.')]
    report.prose(fig, paragraphs)
    report.save(fig, title, paragraphs, 'debelly_duration_curves')

    title = 'Spring: extension feedback depends on cue geometry'
    fig, axes = report.new(title, '3 · What retains the old direction?', 3)
    for ax, geometry, letter in zip(axes, GEOMETRIES, 'ABC'):
        kinds = ['native', 'feedback_off_during_pulse']
        if geometry == 'fixed_slope':
            kinds.append('matched_gradient_peak')
        controls(ax, evidence, 'spring', geometry, kinds,
                 f'{letter}  ' + ('Fixed width 0.2' if geometry == 'fixed_physical' else NAMES[geometry]))
    paragraphs = [
        ('A whole-cell feedback state', 'Spring couples its polarity circuit to moving boundaries. Cell extension increases the deactivation rate throughout the domain, and boundary motion dilutes concentrations. The orange intervention removes only the extension-dependent deactivation term during the cue; the original feedback resumes during recovery.'),
        ('The gradient depends on both contributions at these durations', 'Feedback removal prevents all four sampled gradient cases from reversing (C). Holding the gradient peak fixed also prevents reversal at the three larger sampled lengths. The original gradient therefore benefits from both its growing peak and extension feedback. Its duration is nearly flat: about 0.19 at resting length 0.6 and 0.125 near 2.2.'),
        ('Longer sampled patch pulses can reverse without cue-time feedback', 'For both patches, feedback removal prevents reversal in the shortest sampled cell but permits the two longer cells to reverse (A–B). Original fixed-width durations rise from about 0.15 to 5.82; 20% durations rise from 0.22 to 1.22. The role of extension feedback therefore depends on geometry and the selected pulse duration. Filled symbols pass the complete right-orientation criterion; hollow symbols fail it.')]
    report.prose(fig, paragraphs)
    report.save(fig, title, paragraphs, 'spring_feedback')

    title = 'A successful cue can end while much of the old pole remains'
    fig, axes = report.new(title, 'Cue duration and completed reorientation are different clocks', 2)
    profile(axes[0], evidence, evidence.anchor('wave_pinning', 'fixed_physical', 5),
            'A  Wave-pinning · length 5')
    profile(axes[1], evidence, evidence.anchor('spring', 'fixed_relative', 2.6),
            'B  Spring · resting length 2.6')
    anchors = [r for r in evidence.records.values() if r['job']['kind'] == 'native'
               and 'orientation_pass' in r['job']['labels']]
    before = sum(evidence.metrics[r['job']['slot']]['pulse_ends_before_old_half_removed'] == 'True'
                 for r in anchors)
    paragraphs = [
        ('Follow the profile beyond withdrawal', 'In A, the Wave-pinning cue ends after about 5.9 units with 59% of the old-pole activity still present. In B, the Spring 20% cue ends after about 1.22 units with 66% still present; half-removal is first recorded about 3.72 units after cue onset. Both cells subsequently retain the reversed direction.'),
        ('The measured duration asks how long the cue must act', f'Across the selected successful orientation replays, {before} of {len(anchors)} cues end before half the old-pole activity has disappeared. Withdrawal can therefore occur before the cell finishes reorienting. Subsequent redistribution and mechanical recovery are part of the trajectory even though their time is excluded from the duration curve.'),
        ('What would support a constant-speed removal explanation?', 'The old region would need to widen with cell length, its physical edge would need to retreat at similar speeds, and the removal time would need to explain the plotted response time. The examples here show why a positive duration slope alone is insufficient: the cue can commit the cell while the old region remains, and some models add pole height or mechanical memory instead of width.')]
    report.prose(fig, paragraphs)
    report.save(fig, title, paragraphs, 'commitment_and_recovery')

    title = 'What the comparison explains'
    fig, _ = report.new(title, 'Bring the three questions together')
    paragraphs = [
        ('What must a longer cell overcome?', 'Wave-pinning and Holmes add a broader original active region. A fixed patch covers less of that region as length grows. Goryachev instead concentrates material into a narrow, taller pole and depletes inactive material at the right end. Width and height identify different burdens.'),
        ('How does the cue help?', 'In longer cells, the 20% patch stimulates a larger region than the fixed patch, without changing peak amplitude. The fixed-slope gradient also grows stronger with length, helping initiation in Wave-pinning, Holmes, Goryachev and Spring. Faster inactive transport changes the subsequent competition: it rescues the longest fixed-patch Wave-pinning case but does not rescue the sampled long Holmes patch failures.'),
        ('What retains the old direction?', 'deBelly stores an additional spatial state in membrane–cortex attachments. Holding that state fixed during the cue defeats all twelve sampled native reversals. Spring instead couples the chemical contest to global extension-dependent deactivation; removing that feedback defeats the sampled gradient reversals but leaves the longer patch reversals possible.'),
        ('Which measurements distinguish these explanations?', 'For a growing domain, measure its width and the cue footprint. For a concentrated pole, measure peak height and inactive abundance at the opposite end. For a gradient, compare initiation with recovery and test a fixed peak. For deBelly, follow attachments alongside Rac. For Spring, follow physical extension alongside polarity. Each measurement addresses a process that a single duration-versus-length curve cannot reveal.'),
        ('Predictions worth testing next', 'Increasing cue coverage should preferentially help long cells with broad original domains. Restoring inactive availability should change Goryachev\'s initial activation. Perturbing attachment transport should change deBelly\'s retained state, whereas perturbing extension feedback should affect Spring differently under patches and gradients. These are distinct experiments, even when their unperturbed duration curves look similar.')]
    report.prose(fig, paragraphs, top=.84, width=130, size=11)
    report.save(fig, title, paragraphs, 'synthesis', has_plot=False)

    title = 'Measurements and evidence behind the selected figures'
    fig, _ = report.new(title, 'One reference page')
    paragraphs = [
        ('The two outcome criteria', 'Right orientation requires every recorded sample in the final 500 recovery units to have polarity above 0.05, active contrast at least 0.05, and a polarity range below 0.005. Original-strength replacement additionally requires polarity within 1% of the absolute initial left magnitude. A separate profile check requires pointwise variation and projected drift below 0.1% of fixed field scales. A scalar pass alone does not establish a new full-field steady state.'),
        ('Profiles, amounts and timing', 'Polarity is the normalized spatial first moment of active material. Old-pole activity integrates positive excess above the current mid-cell background in the left half; the new pole uses the right half. The background is the median in the central 5% of the cell, preventing broad gradient activation from counting as a localized pole. Domain width uses half the initial pole contrast above that background. First crossings use recorded samples.'),
        ('How interventions are compared', 'The controls reuse the original selected orientation pulse, or the tested maximum when the original cell did not reverse. Faster inactive diffusion acts during the cue and recovery; in deBelly it changes unbound-attachment diffusion, and in Holmes it changes both inactive chemical pools. Attachment freezing and Spring feedback removal act only during the cue. These comparisons identify effects at the same duration; newly optimized control thresholds were not searched.'),
        ('Coverage and numerical checks', 'The completed mechanism run contains 128 native spatial replays and 107 diagnostic controls across all fifteen model/test combinations. Every archive, starting checkpoint, scalar outcome and conserved pool passed independent audit. The longest-length 401-point checks preserve 14 of 15 orientation classifications and all 15 strength classifications at the sampled doses; deBelly\'s fixed patch is the orientation exception. These are selected parameter regimes in native units, rather than calibrated cross-model biological times.'),
        ('Full evidence remains available', 'The accompanying HTML report links all fifteen detailed diagnostic figure sets, spatial-history figures, trajectory measurements and the independent audit. Full native arrays remain in the recorded mechanism run. The short PDF selects the evidence needed for the argument; the underlying trials and measurements have not been changed.')]
    report.prose(fig, paragraphs, top=.84, width=130, size=10.5)
    report.save(fig, title, paragraphs, 'evidence_reference', has_plot=False)
    report.finish()

    selection = []
    for slot in sorted(evidence.selected):
        r = evidence.records[slot]
        selection.append(dict(slot=slot, archive_sha256=r['archive_sha256'],
                              model=r['job']['case']['model'], geometry=r['job']['case']['condition'],
                              length=r['job']['case']['length'], duration=r['job']['duration'],
                              kind=r['job']['kind']))
    source_paths = [out / 'results/trajectory_metrics.csv', out / 'results/cell_evidence.json',
                    out / 'results/spatial_audit.json', out / 'results/mesh_diagnostic.json',
                    out / 'jobs.json', evidence.campaign / 'results/endpoints.csv',
                    ROOT / 'figures/supplement_width_response.png']
    provenance = dict(generated_at=stamp, evidence_generated_at=evidence.status['generated_at'],
                      producer='finished', report_style='figure_led_reader_guide',
                      page_count=len(report.pages), maximum_pages=25,
                      scientific_data_changed=False,
                      input_hashes={str(p.relative_to(ROOT)): f.digest(p) for p in source_paths},
                      renderer_sha256=f.digest(Path(__file__)),
                      measurement_renderer_sha256=f.digest(HERE / 'three_test_analysis.py'),
                      pdf_sha256=f.digest(pdf_path), html_sha256=f.digest(FIG / 'index.html'),
                      selected_histories=selection,
                      pages=[dict(number=p['number'], title=p['title'], image=p['image'])
                             for p in report.pages])
    (out / 'results/report_provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    print(json.dumps(dict(pdf=str(pdf_path), pages=len(report.pages), selected_histories=len(selection)), indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=Path, default=OUT)
    render(parser.parse_args().out.resolve())
