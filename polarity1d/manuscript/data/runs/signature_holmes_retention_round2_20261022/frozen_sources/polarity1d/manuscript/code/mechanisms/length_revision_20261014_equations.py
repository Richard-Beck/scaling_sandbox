"""Publication-style display mathematics for the retained length report.

Only the two pages containing equations are re-typeset. Every other PDF page
is reused from the preserved pre-typesetting report. Mathtext embeds vector
STIX glyphs, so neither a TeX installation nor raster equation images is needed.
"""
from functools import lru_cache
from pathlib import Path
import hashlib
import html
import io
import json
import re

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.mathtext import MathTextParser
from matplotlib.textpath import TextPath, TextToPath
from pypdf import PdfReader, PdfWriter, Transformation

ROOT = Path(__file__).resolve().parents[2]
FIG = ROOT / 'figures/length_mechanisms_20261008'
OUT = FIG / 'revision_20261014_draft_alignment/equation_typesetting'
PRIOR = FIG / 'archive_20261014_before_equation_typesetting'
WIDTH, MARGIN = 842.4, 38.
USABLE = WIDTH - 2 * MARGIN
PARSER = MathTextParser('path')
METRICS = TextToPath()

# One numbered block for each original equation. Secondary lines define the
# same activation functions or auxiliary quantities included in the original.
EQUATIONS = {
    1: [
        r'\frac{\partial R}{\partial t}=k_R\left[\frac{\alpha}{1+Q^2}-R+\ell_R^2\frac{\partial^2 R}{\partial x^2}\right]',
        r'\alpha=1+s+0.2H_M(b),\qquad H_M(b)=\frac{1}{2}\left[1-\tanh\left(\frac{b-0.068}{0.01}\right)\right]',
    ],
    2: [
        r'V_{\mathrm{left}}=-\frac{V_0 H_R(R_{\mathrm{left}})}{1+T_{\mathrm{left}}^2/T_p^2},\qquad V_{\mathrm{right}}=+\frac{V_0 H_R(R_{\mathrm{right}})}{1+T_{\mathrm{right}}^2/T_p^2}',
        r'H_R(R)=\frac{1}{2}\left[1+\tanh\left(\frac{R-R_{\mathrm{th}}}{w_R}\right)\right]',
    ],
    3: [
        r'\frac{\partial Q}{\partial t}=k_Q\left[\frac{\beta}{1+R^2}-Q+\ell_Q^2\frac{\partial^2 Q}{\partial x^2}\right]',
        r'\beta=1+0.2H_T(T),\qquad H_T(T)=\frac{1}{2}\left[1+\tanh\left(\frac{T-5}{5}\right)\right]',
    ],
    4: [r'\eta\frac{\partial^2 v}{\partial x^2}+\sigma_A\frac{\partial Q}{\partial x}-\chi b\left(v-v_{\mathrm{m}}\right)=0'],
    5: [
        r'K\frac{\partial^2 u}{\partial x^2}-\chi b\left(v_{\mathrm{m}}-v\right)=0',
        r'v_{\mathrm{m}}=\frac{\partial u}{\partial t},\qquad T=K\frac{\partial u}{\partial x}',
    ],
    6: [r'\frac{\partial b}{\partial t}=-v\frac{\partial b}{\partial x}-b\frac{\partial v}{\partial x}+D_b\frac{\partial^2 b}{\partial x^2}+k_{\mathrm{on}}c-k_{\mathrm{off}}b+S_b'],
    7: [r'\frac{\partial c}{\partial t}=-\frac{\partial\left(v_{\mathrm{m}}c\right)}{\partial x}+D_c\frac{\partial^2 c}{\partial x^2}-k_{\mathrm{on}}c+k_{\mathrm{off}}b+S_c'],
    8: [
        r'S_b=\frac{f\,\Delta M}{L\,\Delta t},\qquad S_c=\frac{(1-f)\,\Delta M}{L\,\Delta t}',
        r'f=\frac{k_{\mathrm{on}}}{k_{\mathrm{on}}+k_{\mathrm{off}}}=0.78125',
    ],
}

INLINE = {
    'R': 'R', 'Q': 'Q', 'b': 'b', 'c': 'c', 's': 's', 'α': r'\alpha',
    'β': r'\beta', 'k_R': 'k_R', 'k_Q': 'k_Q', 'ℓ_R²': r'\ell_R^2',
    'ℓ_Q²': r'\ell_Q^2', 'V₀': 'V_0', 'R_th': r'R_{\mathrm{th}}',
    'w_R': 'w_R', 'T': 'T', 'T_p²': 'T_p^2',
    'V_left': r'V_{\mathrm{left}}', 'V_right': r'V_{\mathrm{right}}',
    'v': 'v', 'v_m': r'v_{\mathrm{m}}', 'η': r'\eta', 'σ_A': r'\sigma_A',
    'χ': r'\chi', 'u': 'u', 'K': 'K', 'D_b': 'D_b', 'D_c': 'D_c',
    'k_on': r'k_{\mathrm{on}}', 'k_off': r'k_{\mathrm{off}}',
    'S_b': 'S_b', 'S_c': 'S_c', 'ΔM': r'\Delta{}M', 'L': 'L',
    'Δt': r'\Delta{}t', 'f': 'f',
}
TOKEN = re.compile(r'(?<!\w)(' + '|'.join(re.escape(k) for k in sorted(INLINE, key=len, reverse=True)) + r')(?!\w)')


def configure():
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'pdf.fonttype': 42,
                         'mathtext.fontset': 'stix', 'svg.fonttype': 'path'})


@lru_cache(maxsize=1)
def originals():
    chapters = json.loads((PRIOR / 'chapters.json').read_text())
    return {int(re.match(r'Eq\. (\d+)\.', p).group(1)): p
            for c in chapters for p in c['paragraphs'] if p.startswith('Eq. ')}


def block(paragraph):
    number = int(re.match(r'Eq\. (\d+)\.', paragraph).group(1))
    assert paragraph == originals()[number], 'Equation source changed; review display math before rendering'
    markers = {6: 'The terms are ', 7: 'Unbound MCA follows '}
    marker = markers.get(number, 'Here ')
    description = paragraph[paragraph.index(marker):]
    return number, EQUATIONS[number], description


def inline_pdf(text):
    return TOKEN.sub(lambda m: '$' + INLINE[m.group()] + '$', text)


def inline_html(text):
    # Accessible inline MathML avoids showing underscores or TeX source.
    def tag(token):
        if token in {'ℓ_R²', 'ℓ_Q²'}:
            field = token[2]
            return f'<msubsup><mi>ℓ</mi><mi>{field}</mi><mn>2</mn></msubsup>'
        if token in {'T_p²'}:
            return '<msubsup><mi>T</mi><mi>p</mi><mn>2</mn></msubsup>'
        if token == 'V₀':
            return '<msub><mi>V</mi><mn>0</mn></msub>'
        if '_' in token:
            base, sub = token.split('_', 1)
            subtag = 'mtext' if sub in {'left', 'right', 'th', 'm', 'on', 'off'} else 'mi'
            return f'<msub><mi>{base}</mi><{subtag}>{sub}</{subtag}></msub>'
        if token in {'ΔM', 'Δt'}:
            return f'<mrow><mi mathvariant="normal">Δ</mi><mi>{token[1]}</mi></mrow>'
        return '<mi>' + html.escape(token) + '</mi>'
    pieces, previous = [], 0
    for match in TOKEN.finditer(text):
        pieces.append(html.escape(text[previous:match.start()]))
        pieces.append('<math xmlns="http://www.w3.org/1998/Math/MathML">' + tag(match.group()) + '</math>')
        previous = match.end()
    pieces.append(html.escape(text[previous:]))
    return ''.join(pieces)


@lru_cache(maxsize=24000)
def measure(text, size, weight='normal'):
    font = FontProperties(family='DejaVu Sans', weight=weight, size=size)
    return METRICS.get_text_width_height_descent(text, font, ismath=('$' in text))[0]


def wrap(text, width, size, weight='normal'):
    lines, current = [], ''
    for word in text.split():
        candidate = (current + ' ' + word).strip()
        if current and measure(candidate, size, weight) > width:
            lines.append(current)
            current = word
        else:
            current = candidate
    if current:
        lines.append(current)
    return lines


def equation_rows(number):
    rows = []
    for i, source in enumerate(EQUATIONS[number]):
        size = 14.5 if i == 0 else 13.
        parsed = PARSER.parse('$' + source + '$', dpi=72, prop=FontProperties(size=size))
        assert parsed.width < USABLE - 85, (number, parsed.width)
        rows.append(dict(tex=source, size=size, width=float(parsed.width), height=float(parsed.height),
                         row_height=max(24., float(parsed.height) + 12.)))
    return rows


def export_equation_svgs():
    OUT.mkdir(parents=True, exist_ok=True)
    for number in EQUATIONS:
        rows = equation_rows(number)
        height = sum(r['row_height'] for r in rows) + 8
        fig = plt.figure(figsize=(USABLE / 72, height / 72))
        pos = 4.
        for row in rows:
            fig.text(.5, 1 - pos / height, '$' + row['tex'] + '$', ha='center', va='top',
                     fontsize=row['size'], color='#233d48')
            pos += row['row_height']
        target = OUT / f'equation_{number}.svg'
        fig.savefig(target, format='svg', transparent=True,
                    metadata={'Title': f'Equation ({number})', 'Description': originals()[number]})
        plt.close(fig)


def html_equation(paragraph):
    number, _, description = block(paragraph)
    source = f'revision_20261014_draft_alignment/equation_typesetting/equation_{number}.svg'
    return (f'<div class="equation-display" id="equation-{number}" role="group" aria-label="Equation {number}">'
            f'<img class="equation-math" src="{source}" alt="{html.escape(paragraph, quote=True)}">'
            f'<span class="equation-number">({number})</span></div><p>{inline_html(description)}</p>')


def pdf_report(content, stamp, base):
    configure()
    previous = PdfReader(PRIOR / 'length_mechanisms.pdf')
    writer, layouts, equations = PdfWriter(), [], []
    for n, chapter in enumerate(content, 1):
        if n not in {8, 9}:
            writer.add_page(previous.pages[n - 1])
            writer.add_outline_item(f'Figure {n}: {chapter["title"]}', n - 1)
            layouts.append(dict(figure=n, preserved_page=True))
            continue
        vector = (FIG / chapter['image']).with_suffix('.pdf')
        source = PdfReader(vector).pages[0]
        fw, fh = float(source.mediabox.width), float(source.mediabox.height)
        titles = wrap(chapter['title'], USABLE, 17., 'bold')
        blocks = []
        for paragraph in chapter['paragraphs']:
            if paragraph.startswith('Eq. '):
                number, _, description = block(paragraph)
                rows = equation_rows(number)
                blocks.append(dict(type='equation', number=number, rows=rows,
                                   height=sum(r['row_height'] for r in rows) + 22.))
                lines = wrap(inline_pdf(description), USABLE, 11.)
            else:
                lines = wrap(paragraph, USABLE, 11.)
            blocks.append(dict(type='prose', lines=lines, height=len(lines) * 15.6 + 11.))
        prefix = f'Figure {n}'
        legends = wrap(prefix + '  ' + chapter['caption'][len(prefix):].lstrip(' .'), USABLE, 9.5)
        figure_height = USABLE * fh / fw
        prose_top = 43 + len(titles) * 23 + 14
        plot_top = prose_top + sum(b['height'] for b in blocks) + 9
        plot_bottom = plot_top + figure_height
        height = plot_bottom + 22 + len(legends) * 13.5 + 48
        fig = plt.figure(figsize=(WIDTH / 72, height / 72))
        def put(text, top, size=11., color='#233d48', weight='normal', style='normal', x=MARGIN):
            fig.text(x / WIDTH, 1 - top / height, text, va='top', fontsize=size,
                     color=color, weight=weight, style=style, linespacing=1.42)
        put(f'POLARITY REVERSAL  ·  FIGURE {n}', 22, 8.5, '#63747c', 'bold')
        put('\n'.join(titles), 43, 17., '#173c49', 'bold')
        pos = prose_top
        for item in blocks:
            if item['type'] == 'prose':
                put('\n'.join(item['lines']), pos)
            else:
                top = pos + 8
                row_records = []
                for row in item['rows']:
                    fig.text(.5, 1 - top / height, '$' + row['tex'] + '$', ha='center', va='top',
                             fontsize=row['size'], color='#233d48')
                    row_records.append(dict(center_x_points=WIDTH / 2, top_points=top,
                                            width_points=row['width'], height_points=row['height'],
                                            fontsize_points=row['size']))
                    top += row['row_height']
                number_y = pos + (item['height'] - 6) / 2
                fig.text((WIDTH - MARGIN) / WIDTH, 1 - number_y / height,
                         f'({item["number"]})', ha='right', va='center', fontsize=11., color='#233d48')
                equations.append(dict(number=item['number'], page=n, top_points=pos,
                                      bottom_points=pos + item['height'], rows=row_records,
                                      number_right_points=WIDTH - MARGIN, number_center_y_points=number_y))
            pos += item['height']
        legend_y = plot_bottom + 22
        put(prefix, legend_y, 9.5, '#42565f', style='italic')
        prefix_w = TextPath((0, 0), prefix + '  ', size=9.5,
                            prop=FontProperties(family='DejaVu Sans', style='italic')).get_extents().width + 3
        put(legends[0][len(prefix):].lstrip(), legend_y, 9.5, '#42565f', x=MARGIN + prefix_w)
        if len(legends) > 1:
            put('\n'.join(legends[1:]), legend_y + 13.5, 9.5, '#42565f')
        fig.text(MARGIN / WIDTH, 17 / height, 'Generated ' + stamp[:19] + ' UTC · native-data producers finished',
                 fontsize=7.5, color='#63747c')
        fig.text((WIDTH - MARGIN) / WIDTH, 17 / height, str(n), ha='right', fontsize=8, color='#63747c')
        buffer = io.BytesIO()
        fig.savefig(buffer, format='pdf')
        plt.close(fig)
        buffer.seek(0)
        page = PdfReader(buffer).pages[0]
        scale = USABLE / fw
        page.merge_transformed_page(source, Transformation().scale(scale).translate(
            MARGIN - float(source.mediabox.left) * scale,
            height - plot_bottom - float(source.mediabox.bottom) * scale))
        writer.add_page(page)
        writer.add_outline_item(f'Figure {n}: {chapter["title"]}', n - 1)
        layouts.append(dict(figure=n, page_width_points=WIDTH, page_height_points=height,
                            plot_bounds=[MARGIN, height - plot_bottom, WIDTH - MARGIN, height - plot_top],
                            preserved_page=False))
    writer.add_metadata({'/Title': base.TITLE, '/Author': 'Polarity manuscript working report',
                         '/Subject': 'Chemical reorganization and membrane–cortex attachment memory'})
    with (FIG / 'length_mechanisms.pdf').open('wb') as stream:
        writer.write(stream)
    export_equation_svgs()
    record = dict(font='STIX', equations=equations, layout=layouts, producer='finished',
                  equation_tex=EQUATIONS, original_equations_unchanged=True,
                  renderer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (OUT / 'equation_layout.json').write_text(json.dumps(record, indent=2) + '\n')
    (base.DATA / 'report_layout.json').write_text(json.dumps(layouts, indent=2) + '\n')
