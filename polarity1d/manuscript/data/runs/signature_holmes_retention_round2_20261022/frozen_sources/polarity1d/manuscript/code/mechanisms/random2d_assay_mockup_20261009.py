"""Measurement-first 2D random-migration assay mockup; no model integrations.

Real HT1080 contours provide only a drawing template. Size groups, added
protrusions, event evolution and every numerical comparison are synthetic.
"""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Polygon
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[2]
REPO = ROOT.parents[1]
FIG = ROOT / 'figures/random2d_assay_mockup_20261009'
DATA = ROOT / 'data/random2d_assay_mockup_20261009'
FIG.mkdir(parents=True, exist_ok=True)
DATA.mkdir(parents=True, exist_ok=True)

INK = '#203342'
MUTED = '#61717c'
SMALL = '#008d91'
LARGE = '#d66c32'
EVENT = '#8054a1'
PALE = '#f6f8fa'
RULE = '#dbe3e8'
OLD = '#526d7a'
plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 11,
    'text.color': INK, 'axes.labelcolor': INK, 'xtick.color': MUTED,
    'ytick.color': MUTED, 'axes.edgecolor': '#a8b5be',
    'axes.spines.top': False, 'axes.spines.right': False,
    'axes.titlesize': 13, 'axes.labelsize': 10.5,
    'xtick.labelsize': 9, 'ytick.labelsize': 9,
    'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
    'lines.linewidth': 2.5, 'savefig.facecolor': 'white',
})

sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
outlines_path = DATA / 'ht1080_outline_examples.json'
outlines = json.loads(outlines_path.read_text())
raw = np.array(outlines['records'][0]['outline_xy_um'][2], dtype=float)
raw[:, 1] *= -1
raw -= raw.mean(axis=0)
_, eig = np.linalg.eigh(np.cov(raw.T))
base = raw @ eig[:, ::-1]
if np.std(base[base[:, 0] > 0, 1]) < np.std(base[base[:, 0] < 0, 1]):
    base[:, 0] *= -1
base /= np.ptp(base[:, 0])
base *= 2.0
angles = np.arctan2(base[:, 1], base[:, 0])

fig = plt.figure(figsize=(17.8, 14.8), facecolor='white')
tracked_text = []


def text(x, y, s, **kw):
    artist = fig.text(x, y, s, **kw)
    tracked_text.append(artist)
    return artist


def box(x, y, w, h, face=PALE):
    p = FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.008,rounding_size=0.009',
                      transform=fig.transFigure, facecolor=face,
                      edgecolor=RULE, linewidth=.75, zorder=-5)
    fig.add_artist(p)


def header(x, y, letter, title, subtitle):
    text(x, y, letter, fontsize=18, fontweight='bold', va='top', color=INK)
    text(x + .024, y, title, fontsize=14, fontweight='bold', va='top')
    text(x + .024, y - .024, subtitle, fontsize=10.2, color=MUTED, va='top')


def arrow(ax, a, b, color=INK, lw=1.4, style='-|>', **kw):
    ax.annotate('', xy=b, xytext=a,
                arrowprops=dict(arrowstyle=style, lw=lw, color=color,
                                shrinkA=0, shrinkB=0, **kw))


def cell(ax, center, scale=1.0, rotation=0, color=OLD,
         bump=0, old_strength=1, alpha=.10, previous=False):
    p = base.copy()
    # Construct a lateral competitor on a real contour template. These
    # deformations are illustrative and are not native HT1080 time samples.
    delta = np.angle(np.exp(1j*(angles - 1.25)))
    p[:, 1] += bump * np.exp(-.5*(delta/.31)**2)
    front = np.exp(-.5*(np.angle(np.exp(1j*angles))/.40)**2)
    p[:, 0] -= (1-old_strength)*.32*front
    r = np.deg2rad(rotation)
    rot = np.array([[np.cos(r), -np.sin(r)], [np.sin(r), np.cos(r)]])
    p = (p @ rot.T) * scale + center
    ax.add_patch(Polygon(p, closed=True, facecolor=color if not previous else 'none',
                         edgecolor=color, alpha=alpha if not previous else 1,
                         linewidth=0, zorder=2))
    ax.plot(*np.vstack([p, p[0]]).T, color=color, lw=1.6,
            ls=(0, (3, 2)) if previous else '-', alpha=.55 if previous else 1,
            zorder=3)
    if not previous:
        idx = np.where(abs(delta) < .32)[0]
        ax.plot(*p[idx].T, color=EVENT, lw=4.0, solid_capstyle='round', zorder=5)
        idx_old = np.where(abs(np.angle(np.exp(1j*angles))) < .37)[0]
        ax.plot(*p[idx_old].T, color=OLD, lw=1+2*old_strength,
                alpha=max(.15, old_strength), solid_capstyle='round', zorder=5)
    return p


def axes_clean(ax):
    ax.spines[['top', 'right']].set_visible(False)
    ax.spines[['left', 'bottom']].set_linewidth(.7)
    ax.tick_params(length=3, width=.65, pad=4)
    ax.grid(axis='y', color=RULE, linewidth=.65, alpha=.7, zorder=0)
    ax.set_axisbelow(True)


def lines(ax, x, small, large):
    ax.plot(x, small, color=SMALL, label='Small', zorder=3)
    ax.plot(x, large, color=LARGE, label='Large', zorder=3)


def card(col, row, letter, title, measure, expectation, mechanism, refs):
    x = .038 + col * .321
    y = .359 if row == 0 else .071
    w, h = .300, .266
    box(x, y, w, h, 'white')
    header(x + .010, y + h - .011, letter, title, measure)
    ax = fig.add_axes([x + .045, y + .092, w - .065, .112])
    axes_clean(ax)
    text(x + .012, y + .058, expectation, fontsize=10.6, fontweight='bold', va='top')
    text(x + .012, y + .039, mechanism, fontsize=9.8, color=EVENT, va='top')
    text(x + .012, y + .022, refs, fontsize=9.0, color=MUTED, va='top')
    return ax


text(.038, .970, 'What makes a new protrusion become a lasting turn?',
     fontsize=25, fontweight='bold', va='top')
text(.038, .936, 'A measurement-first design for random migration of small and large cells on a 2D substrate',
     fontsize=13.5, color=MUTED, va='top')
text(.038, .910, 'HYPOTHETICAL ASSAY  /  SYNTHETIC COMPARISONS', fontsize=10.5,
     color=EVENT, fontweight='bold', va='top')
for xpos, col, label in [(.765, SMALL, 'Small'), (.844, LARGE, 'Large')]:
    fig.add_artist(Line2D([xpos, xpos+.026], [.912, .912],
                         transform=fig.transFigure, color=col, lw=3))
    text(xpos+.032, .912, label, fontsize=11, va='center')

# A: source-derived geometry and observable annotation.
box(.038, .659, .325, .221)
header(.050, .869, 'A', 'Track the contour and motion',
       'HT1080-derived geometry; constructed size contrast')
ax = fig.add_axes([.047, .698, .309, .126])
ax.set(xlim=(0, 10), ylim=(-.2, 4.4)); ax.set_aspect('equal'); ax.axis('off')
p_small = cell(ax, (2.1, 1.60), scale=1.10, color=SMALL, bump=.23)
p_large_old = cell(ax, (6.9, 1.75), scale=1.72, color='#a2afb8',
                   previous=True)
p_large = cell(ax, (6.9, 1.75), scale=1.72, color=LARGE, bump=.25)
for c in [(2.1, 1.60), (6.9, 1.75)]:
    ax.scatter(*c, color=INK, s=12, zorder=8)
arrow(ax, (1.65, .15), (2.65, .15), color=SMALL)
arrow(ax, (6.0, .15), (7.5, .15), color=LARGE)
ax.text(2.1, 3.9, 'Small', color=SMALL, ha='center', fontsize=11, fontweight='bold')
ax.text(6.9, 3.9, 'Large', color=LARGE, ha='center', fontsize=11, fontweight='bold')
# Choose representative event/old-front positions from drawn contours.
new_i = np.argmin(abs(np.angle(np.exp(1j*(angles-1.25)))))
old_i = np.argmax(base[:, 0])
newpoint, oldpoint = p_large[new_i], p_large[old_i]
arc_idx = np.where(abs(np.angle(np.exp(1j*(angles-1.25)))) < .48)[0]
ax.add_patch(Polygon(np.vstack([p_large[arc_idx], p_large_old[arc_idx][::-1]]),
                     closed=True, facecolor=EVENT, alpha=.25, edgecolor='none', zorder=4))
arrow(ax, newpoint, oldpoint, color=MUTED, lw=1, style='<->')
ax.text(*(0.5*(newpoint+oldpoint)+[.03, -.02]), 'd', fontsize=11,
        color=MUTED, fontstyle='italic', zorder=10,
        bbox=dict(facecolor=PALE, edgecolor='none', pad=.4))
# A physical footprint is the purple contour arc; end ticks identify its extent.
width_idx = np.where(abs(np.angle(np.exp(1j*(angles-1.25)))) < .32)[0]
for ii in [width_idx[0], width_idx[-1]]:
    q = p_large[ii]
    ax.plot([q[0], q[0]], [q[1]-.09, q[1]+.09], color=EVENT, lw=1.3, zorder=9)
ax.annotate('width w', xy=newpoint, xytext=(4.0, 3.15), fontsize=10,
            color=EVENT, arrowprops=dict(arrowstyle='-', color=EVENT, lw=1), va='center')
ax.annotate('area gain', xy=.5*(newpoint+p_large_old[new_i]), xytext=(4.0, 2.55), fontsize=9.3,
            color=EVENT, arrowprops=dict(arrowstyle='-', color=EVENT, lw=.8), va='center')
ax.annotate('old front', xy=oldpoint, xytext=(8.2, 3.1), fontsize=10,
            color=OLD, arrowprops=dict(arrowstyle='-', color=OLD, lw=1), va='center')
text(.051, .682, 'Contour: width w, area gain, front separation d', fontsize=10.4, va='center')
text(.051, .667, 'Motion: heading, speed, rear retraction; retain shape history', fontsize=9.6, color=MUTED, va='center')

# B: a constructed event timeline based on the same source contour.
box(.387, .659, .576, .221)
header(.399, .869, 'B', 'Follow every attempt through commitment and recovery',
       'Align events at competitor onset; retain attempts that retract or turn back')
ax = fig.add_axes([.400, .726, .550, .101])
ax.set(xlim=(-.45, 10.1), ylim=(-1.1, 1.65)); ax.set_aspect('equal'); ax.axis('off')
stages = [(0.8, 0, .08, 1), (3.2, 0, .42, .90),
          (5.7, 35, .68, .45), (8.4, 78, .85, .15)]
for i, (x, rotation, bump, old_strength) in enumerate(stages):
    p = cell(ax, (x, -.05), scale=.77, rotation=0, color=OLD,
             bump=bump, old_strength=old_strength, alpha=.07)
    direction = np.array([np.cos(np.deg2rad(rotation)), np.sin(np.deg2rad(rotation))])
    arrow(ax, np.array([x, -.79]), np.array([x, -.79]) + .48*direction, color=INK)
    ax.text(x, 1.46, ['Onset', 'Competition', 'First deflection', 'Sustained turn'][i],
            ha='center', fontsize=10.3, fontweight='bold')
    if i < 3:
        arrow(ax, (x+1.05, .15), (stages[i+1][0]-1.1, .15), color='#adb9c1', lw=1.1)

tx = [.432, .563, .699, .847]
ty = .706
fig.add_artist(Line2D([tx[0], tx[-1]], [ty, ty], transform=fig.transFigure,
                     color=RULE, lw=2))
for x, label in zip(tx, ['onset', 'old front responds', 'heading changes', 'new heading holds']):
    fig.add_artist(Line2D([x], [ty], marker='o', markersize=5,
                         transform=fig.transFigure, color=EVENT))
    text(x, ty-.012, label, fontsize=9.0, color=MUTED, ha='center', va='top')
text(.405, .670, 'Measure separate clocks: early growth  →  competition  →  turn  →  recovery',
     fontsize=10.3, fontweight='bold', va='center')

text(.038, .643, 'Six observable signatures', fontsize=13, fontweight='bold', va='center')
text(.258, .643, 'Each curve pair illustrates a separate mechanism-compatible scenario to test.',
     fontsize=10.6, color=MUTED, va='center')

synthetic = []


def record(panel, group, x, y, x_name, y_name):
    for xx, yy in zip(x, y):
        synthetic.append(dict(panel=panel, size_group=group, x_name=x_name,
                              x=float(xx), y_name=y_name, y=float(yy), synthetic=True))


# C: direct early-event proxy, absolute area gain at matched onset width.
ax = card(0, 0, 'C', 'Early competitor growth',
          'Area gain over time; compare a fixed early window',
          'Large: weaker early expansion',
          'Available-pool depletion / harder competitor growth',
          'Goryachev-like · Fig. 3 · match onset width and angle')
t = np.linspace(0, 2, 121)
ys = 15*(1-np.exp(-t/1.3))
yl = 7.5*(1-np.exp(-t/1.7))
lines(ax, t, ys, yl)
ax.axvspan(0, .5, color=EVENT, alpha=.055, linewidth=0, zorder=-1)
ax.set(xlim=(0, 2), ylim=(0, 14), xticks=[0, 1, 2], yticks=[0, 5, 10],
       xlabel='Time after onset (min)', ylabel='New protruded area (µm²)')
ax.text(.95, .94, 'include aborted attempts', transform=ax.transAxes,
        ha='right', va='top', fontsize=9, color=MUTED)
record('C', 'small', t, ys, 'time_min', 'new_area_um2')
record('C', 'large', t, yl, 'time_min', 'new_area_um2')

# D: initial width avoids conditioning the predictor on a later successful state.
ax = card(1, 0, 'D', 'Width needed for commitment',
          'Lasting-turn probability versus initial event width',
          'Large: wider seeds needed for the same success',
          'A broader polarized domain must establish and persist',
          'Wave-pinning-like · Figs. 2, 6 · also test width / perimeter')
w = np.linspace(0, 14, 141)
ys = .04+.88/(1+np.exp(-(w-4.2)/1.2))
yl = .04+.88/(1+np.exp(-(w-8.2)/1.6))
lines(ax, w, ys, yl)
ax.set(xlim=(0, 14), ylim=(0, 1), xticks=[0, 5, 10], yticks=[0, .5, 1],
       xlabel='Initial competitor width (µm)', ylabel='P(lasting turn)')
ax.axhline(.5, color=RULE, lw=.9, ls=':')
arrow(ax, (4.3, .52), (8.3, .52), color=EVENT, lw=1.2, style='<->')
record('D', 'small', w, ys, 'initial_width_um', 'probability_lasting_turn')
record('D', 'large', w, yl, 'initial_width_um', 'probability_lasting_turn')

# E: geometric extent is defined even in a steadily translating cell.
ax = card(2, 0, 'E', 'Old-front resistance and rebound',
          'Old-front extent during and after a transient dip',
          'Large: less disruption; greater rebound per loss',
          'Stronger local restoration of the incumbent pole',
          'Goryachev-like · Fig. 4 · early rebound ≠ full recovery')
t = np.linspace(-1.5, 6, 241)
ys = np.where(t < 0, 1-.50*np.exp(t/.23), 1-.50*np.exp(-t/1.05))
yl = np.where(t < 0, 1-.22*np.exp(t/.23),
              1-.16*np.exp(-t/.20)-.06*np.exp(-t/4))
lines(ax, t, ys, yl)
ax.axvline(0, color=RULE, lw=1, ls=':')
ax.set(xlim=(-1.5, 6), ylim=(.45, 1.07), xticks=[0, 2, 4, 6], yticks=[.5, .75, 1],
       xlabel='Time from old-front minimum (min)', ylabel='Old-front extent / baseline')
ax.text(.97, .09, 'dip → rebound', transform=ax.transAxes, ha='right',
        fontsize=9, color=MUTED)
record('E', 'small', t, ys, 'time_from_minimum_min', 'old_front_extent_normalized')
record('E', 'large', t, yl, 'time_from_minimum_min', 'old_front_extent_normalized')

# F: no asserted diffusion exponent; arbitrary smooth increasing candidate curves.
ax = card(0, 1, 'F', 'Response delay versus separation',
          'Competitor onset → old-front recession',
          'Large: longer delay, even at a matched gap',
          'Separation penalty + residual local / material burden',
          'Goryachev / Wave-pinning · Figs. 1, 2, 5 · compare shared gaps')
ds = np.linspace(4, 21, 81)
dl = np.linspace(4, 36, 121)
ys = 1.2+.11*ds+.002*ds**2
yl = 3.5+.11*dl+.002*dl**2
ax.plot(ds, ys, color=SMALL)
ax.plot(dl, yl, color=LARGE)
gap = 14
y0 = 1.2+.11*gap+.002*gap**2
arrow(ax, (gap, y0+.10), (gap, y0+2.20), color=EVENT, style='<->', lw=1.2)
ax.text(gap+1, y0+1.35, 'size effect at fixed d', fontsize=8.6, color=EVENT,
        va='center', bbox=dict(facecolor='white', edgecolor='none', pad=1))
ax.set(xlim=(0, 38), ylim=(0, 11), xticks=[0, 10, 20, 30], yticks=[0, 5, 10],
       xlabel='Old–new front separation d (µm)', ylabel='Old-front response delay (min)')
record('F', 'small', ds, ys, 'separation_um', 'old_front_response_delay_min')
record('F', 'large', dl, yl, 'separation_um', 'old_front_response_delay_min')

# G: cumulative return, not persistence or a return hazard.
ax = card(1, 1, 'G', 'Return after an initial turn',
          'Return to the old heading after a qualifying deflection',
          'Large, slow-memory regime: more delayed returns',
          'Visible redirection can precede mechanical commitment',
          'deBelly-like · Figs. 8–11 · follow shape / speed recovery too')
t = np.linspace(0, 12, 121)
ys = .20*(1-np.exp(-t/1.3))
yl = .58*(1-np.exp(-t/4.3))
lines(ax, t, ys, yl)
ax.set(xlim=(0, 12), ylim=(0, .65), xticks=[0, 4, 8, 12], yticks=[0, .3, .6],
       xlabel='Time after initial turn (min)', ylabel='Cumulative P(return to old heading)')
record('G', 'small', t, ys, 'time_after_turn_min', 'cumulative_return_probability')
record('G', 'large', t, yl, 'time_after_turn_min', 'cumulative_return_probability')

# H: a deliberately exploratory size × history example, no established sign.
ax = card(2, 1, 'H', 'Shape history predicts event fate?',
          'Match current shape; vary preceding extension history',
          'Test whether history changes the size effect',
          'Extension feedback + coupled polarity–mechanics',
          'Spring / deBelly · Figs. 1, 9–11 · sign is to be measured')
h = np.linspace(-.25, .30, 121)
ys = .40+.10*np.tanh(h/.12)
yl = .40+.29*np.tanh(h/.12)
lines(ax, h, ys, yl)
ax.axvline(0, color=RULE, lw=1, ls=':')
ax.set(xlim=(-.25, .3), ylim=(0, .8), xticks=[-.2, 0, .2], yticks=[0, .4, .8],
       xlabel='Prior length change / current length', ylabel='P(lasting turn)')
ax.text(.02, .94, 'one possible interaction', transform=ax.transAxes,
        fontsize=9, color=MUTED, va='top')
record('H', 'small', h, ys, 'preceding_length_change_normalized', 'probability_lasting_turn')
record('H', 'large', h, yl, 'preceding_length_change_normalized', 'probability_lasting_turn')

text(.038, .044, 'Shared design', fontsize=11, fontweight='bold', va='center')
text(.127, .044, 'Count attempts per cell and per perimeter; retain failures; analyze size continuously; match event geometry and current shape.',
     fontsize=10.1, va='center')
text(.038, .025, 'All curves and numbers are illustrative, not fitted data. Mechanisms motivate behavioral hypotheses; contour outputs do not identify molecular fields.',
     fontsize=9.8, color=MUTED, va='center')
text(.038, .011, 'Outline template: HT1080 source 0005, track 0, frame 146. Event concept: Bosgraaf & Van Haastert, PLOS ONE (2009), doi:10.1371/journal.pone.0005253.',
     fontsize=8.5, color=MUTED, va='center')

CAPTION = '''Hypothetical reporter-free measurements in a 2D random-migration assay. All size contrasts, event evolution, numerical values and curves are constructed illustrations, not HT1080 measurements, simulations of the characterized models, or fitted predictions. The contour template alone comes from previously segmented HT1080 data (source 0005, track 0, frame 146; 2-µm spatial smoothing); the same template is rescaled and deformed for the drawings.

A. Track contours, a body reference/centroid, heading and speed. Record protrusion footprint width in physical units and as a perimeter fraction, new protruded area, physical separation from the incumbent front, front–rear length, area and aspect ratio. Purple denotes the competing protrusive region, slate the incumbent front, and dashed contours a constructed earlier state. Rigid translation/rotation and deformation should be stored separately; contour correspondence does not imply permanent membrane material tracking.

B. Align all detected competitors at onset. Follow expansion, the old-front response, first qualifying deflection, a sustained new heading and subsequent recovery. Retain aborted attempts and returns. The stages are schematic and are not a time sequence from the source HT1080 movie. Keep waiting, turn-completion and recovery clocks separate. A workable pilot definition is a deflection of at least 60 degrees from the pre-event heading, a sustained new heading within a 30-degree band for 60 s, and a speed floor calibrated to tracking noise; freeze the actual thresholds before size comparisons. These are illustrative analysis choices, not biologically validated constants.

C. New protruded area relative to the registered onset contour over a fixed early window, including aborted events. A weaker large-cell curve is motivated by lower available inactive concentration and harder competitor growth in the characterized Goryachev regime (length-mechanisms Fig. 3). Match initial footprint width, event angle and incumbent state. A visible extension response is a downstream proxy, not a measurement of inactive concentration or molecular activation.

D. Lasting-turn probability versus initial competitor footprint width. The larger-cell right shift illustrates the possibility that a larger polarized domain needs a wider seed to establish a retained replacement (Wave-pinning, Figs. 2 and 6). Measure subsequent footprint expansion and lifetime too. Test whether using width/perimeter reduces the size difference; the report does not establish a universal fractional threshold or collapse. Initial rather than final width avoids defining the predictor by a later successful state.

E. Old-front projected extent relative to a consistent body reference, normalized to its own pre-event baseline. A smaller disruption and stronger early rebound illustrate the downstream behavioral signature of stronger incumbent restoration (Goryachev, Fig. 4). Match competitor geometry and early growth. Measure rebound relative to the lost extent above a noise floor. The large-cell illustration retains a slower tail: initial restoration and complete recovery need not have the same length dependence. Visible extent is not chemical pole concentration.

F. Time from competitor onset to sustained old-front recession versus physical old–new sector separation. Larger separation and a residual size difference at equal gap are motivated by the separate distance and fixed-gap effects in inactive-sink controls (Goryachev and Wave-pinning, Fig. 5), together with greater material and reflected-target distance (Figs. 1–2). The sink is a diagnostic perturbation, not an observed competing protrusion; the displayed minutes and curve shapes are arbitrary extrapolations to an assay readout. Record onset-to-stable-heading and its component clocks as well. In real data, unresolved old-front responses remain censored rather than being discarded.

G. Cumulative probability of returning to a predeclared old-heading band after an initial qualifying turn. One slow-memory scenario yields more delayed returns in large cells. This is motivated by rapid Rac response and slower membrane–cortex attachment/mechanical redistribution and recovery in corrected deBelly trials (Figs. 8–11), but is not a universal monotonic prediction. Compare returns toward the old heading with equal-angle control heading bands to distinguish directional memory from broadly increased turning. Measure stable-heading, rear-retraction, shape-recovery and speed-recovery clocks separately. Membrane–cortex attachment is distinct from cell–substrate adhesion. No prescribed model endpoint-relaxation time is transferred to this assay.

H. Conditional lasting-turn probability versus the preceding fractional change in front–rear length over a fixed pre-onset window (for example 5 min). Compare history at matched current dimensionless shape descriptors, speed and event geometry while retaining current size as a continuous predictor and history interaction. Current size is not matched away between small/large display groups. The illustrated interaction is only one possible association; neither its sign nor its strength is established. Spring extension-feedback controls and deBelly coupled mechanical redistribution motivate testing whether shape history adds predictive information beyond current shape, and whether it changes the apparent size effect (Figs. 1 and 9–11). Recent extension/contraction, cumulative positive extension, orientation history and front–rear coordination are useful complementary predictors.

Small/large colors are display bins; the intended analysis uses continuous area/size, physical length and perimeter separately. Count initiation per cell and per perimeter, and separate event frequency from event success. Use identical noise-aware detection thresholds, retain overlapping events as a declared class, and use cell-level resampling or hierarchical analysis so events are not treated as independent biological replicates. Synthetic curves carry no error bars, significance tests or experimental uncertainty. Different cards depict alternative mechanism-compatible regimes and do not describe a single jointly simulated cell population. A larger cell is not assumed to be uniformly slower: faster mechanical redistribution and larger burdens can compete, and size trends may be nonmonotonic. Molecular active half-mass motion in Fig. 7 is not equated with cell migration speed.

Reference: Bosgraaf L, Van Haastert PJM (2009), The Ordered Extension of Pseudopodia by Amoeboid Cells in the Absence of External Cues. PLOS ONE 4(4): e5253. https://doi.org/10.1371/journal.pone.0005253
'''

(FIG / 'CAPTION.txt').write_text(CAPTION)
with (DATA / 'synthetic_curves.csv').open('w') as f:
    writer = csv.DictWriter(f, fieldnames=list(synthetic[0]))
    writer.writeheader()
    writer.writerows(synthetic)

specifications = [
    ('C', 'Early competitor growth', 'New protruded area over a fixed early window',
     'area_um2; width_um; normal_speed_um_per_min', 'matched onset width; angle; old-front state',
     'weaker large-cell growth', 'Goryachev inactive depletion', '3'),
    ('D', 'Commitment versus initial width', 'Probability of a sustained new heading',
     'initial_width_um; width_perimeter_fraction; outcome', 'event angle; early amplitude; current state',
     'larger physical seed needed', 'Wave-pinning retained domain expansion', '2;6'),
    ('E', 'Old-front rebound', 'Projected incumbent extent and initial rebound per extent loss',
     'extent_um; extent_ratio; rebound_per_min', 'competitor width; angle; early growth; pre-event extent',
     'smaller dip; stronger initial rebound', 'Goryachev local restoration', '4'),
    ('F', 'Delay versus separation', 'Competitor-onset to sustained old-front recession',
     'separation_um; old_front_delay_min; censored', 'physical separation; current shape; competitor geometry',
     'longer delay; residual size difference at matched gap', 'distance and local/material burden', '1;2;5'),
    ('G', 'Return after an initial turn', 'Cumulative return to old heading after qualifying deflection',
     'initial_turn_time; return_time; censored; recovery_times', 'initial deflection; speed; shape; event geometry',
     'more delayed returns in one slow-memory regime', 'deBelly mechanical memory', '8;9;10;11'),
    ('H', 'Shape-history dependence', 'Outcome dependence on preceding extension at matched current state',
     'length_history; current_shape; event_geometry; outcome',
     'current dimensionless shape; speed; event geometry; continuous size covariate',
     'test whether history changes size effect; no established sign', 'Spring/deBelly coupled mechanics', '1;9;10;11'),
]
with (DATA / 'measurement_specifications.csv').open('w') as f:
    writer = csv.writer(f)
    writer.writerow(['panel', 'measurement', 'definition', 'retain', 'condition_on',
                     'illustrated_signature', 'compatible_mechanism', 'length_report_figures'])
    writer.writerows(specifications)

fig.canvas.draw()
renderer = fig.canvas.get_renderer()
overflow = []
for artist in tracked_text:
    bb = artist.get_window_extent(renderer).transformed(fig.transFigure.inverted())
    if bb.x0 < 0 or bb.y0 < 0 or bb.x1 > 1 or bb.y1 > 1:
        overflow.append(artist.get_text())
assert not overflow, overflow

output_hashes = {}
for extension in ['pdf', 'svg', 'png']:
    p = FIG / ('random2d_measurements.' + extension)
    fig.savefig(p, dpi=220, metadata={'Creator': 'scaling_sandbox; hypothetical assay mockup'}
                if extension == 'pdf' else None)
    output_hashes[p.name] = sha(p)
plt.close(fig)

manifest = dict(
    generated_at=datetime.now(timezone.utc).isoformat(), producer='finished',
    hypothetical=True, all_numeric_comparisons_synthetic=True,
    curve_generation='Deterministic illustrative functions; no stochastic sampling, fitted data, or model integrations.',
    source_template=dict(source='HT1080 0005', track=0, frame=146,
                         smoothing_um=2, size_groups_observed=False,
                         transformations='centered, principal-axis aligned, rescaled; constructed local protrusion and front recession; rotation'),
    sources_sha256={
        str(Path(__file__).relative_to(REPO)): sha(__file__),
        str(outlines_path.relative_to(REPO)): sha(outlines_path),
        'polarity1d/manuscript/figures/length_mechanisms_20261008/length_mechanisms.pdf':
            sha(ROOT/'figures/length_mechanisms_20261008/length_mechanisms.pdf'),
        'literature_review/bosgraaf_2009.md': sha(REPO/'literature_review/bosgraaf_2009.md'),
    },
    outputs_sha256=output_hashes,
    data_sha256={p.name: sha(p) for p in [DATA/'synthetic_curves.csv', DATA/'measurement_specifications.csv']},
    caption_sha256=sha(FIG/'CAPTION.txt'), figure_text_within_page=True,
    reproduction='OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python polarity1d/manuscript/code/mechanisms/random2d_assay_mockup_20261009.py',
)
(FIG/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
(FIG/'README.txt').write_text('''Hypothetical 2D random-migration assay figure

Open random2d_measurements.pdf for the vector figure, random2d_measurements.svg
for editable vector artwork, or random2d_measurements.png for a quick preview.
CAPTION.txt defines the observables, mechanism links, example event rules and
the status of every synthetic element. This figure does not modify the current
length-mechanism report or any model results.

Synthetic curve coordinates and measurement specifications:
../../data/random2d_assay_mockup_20261009/

Reproduce from the repository root:
'''+manifest['reproduction']+'\n')
print(json.dumps({'figure_directory': str(FIG), 'outputs': list(output_hashes),
                  'synthetic_curve_rows': len(synthetic), 'text_bounds_pass': True}))
