"""Draw a plainly unmeasured main-figure placeholder; no simulated observations."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

OUT = Path(__file__).resolve().parent / 'assets'
OUT.mkdir(exist_ok=True)
fig = plt.figure(figsize=(11.7, 8.3), facecolor='white')
fig.text(.06, .93, 'RESERVED MAJOR EXPERIMENTAL RESULTS', fontsize=20,
         weight='bold', color='#173e54')
fig.text(.06, .885, 'New WGD motility assays — all data pending', fontsize=16)
panels = [
    ('A', 'Paired biological system',
     'SUM-159 chemically induced WGD+/− pilot\nConfirm group identity, size and morphology'),
    ('B', 'Imaging and observed tracks',
     'Optimize acquisition and tracking\nDeclare cadence, duration, eligibility and replicates'),
    ('C', 'Measured migration behaviour',
     'Speed, path length, displacement and directness\nLag-dependent MSD under declared observation rules'),
    ('D', 'Evaluate the persistence hypothesis',
     'Relate the acquired observations to WGD and size\nRetain null or opposite outcomes; specify controls'),
]
for (letter, title, text), (x, y) in zip(panels, [(0.06, .51), (.53, .51), (.06, .15), (.53, .15)]):
    fig.add_artist(FancyBboxPatch((x, y), .41, .28, transform=fig.transFigure,
                   boxstyle='round,pad=0.012', facecolor='#f3f6f7',
                   edgecolor='#b4c4ca', linestyle='--', linewidth=1.5))
    fig.text(x+.012, y+.225, letter, fontsize=20, weight='bold', color='#173e54')
    fig.text(x+.05, y+.23, title, fontsize=12.5, weight='bold')
    fig.text(x+.012, y+.13, text, fontsize=10.5, linespacing=1.8)
    fig.text(x+.012, y+.035, 'PENDING — no observations available',
             fontsize=11, color='#906012', weight='bold')
fig.text(.06, .055, 'Proposed panel plan only. Design and analysis remain to be agreed before interpreting data.', fontsize=11)
for ext in ['png', 'pdf']:
    fig.savefig(OUT / ('planned_WGD_results.'+ext), dpi=200)
plt.close(fig)
