"""Validated fit tables and shared manuscript panel drawing."""
from pathlib import Path
import sys, json, csv, hashlib
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.ticker import LogLocator, FuncFormatter
ROOT=Path(__file__).resolve().parents[3]
CAMPAIGN=ROOT/'wgd_langevin/refined_ou_20261003'
RESULTS_DIR=CAMPAIGN/'results'
FIGURE_DIR=CAMPAIGN/'figures'
RELATIVE_SURFACE_CSV=RESULTS_DIR/'hosny_sensitivity/hosny_surface.csv'
NAVY='#19394d'
STUDY_LABELS={'Xuan2018':'Xuan 2018','Mallin2023':'Mallin 2023','Jemaa2017':'Jemaa 2017','Jemaa2023':'Jemaa 2023','Hosny_baseline':'Hosny'}
STUDY_ORDER=list(STUDY_LABELS)
METRIC_ORDER={'Xuan2018':['speed','msd'],'Mallin2023':['path_length','directness'],'Jemaa2017':['speed'],'Jemaa2023':['speed'],'Hosny_baseline':['path_length','directness','parallel']}
METRIC_LABELS={'speed':'speed','msd':'MSD at last shared population lag','path_length':'path length','directness':'directness','parallel':'signed migration index'}
def load_fit_scenarios():
    out={}
    for path in sorted(RESULTS_DIR.glob('fit_*.json')):
        record=json.loads(path.read_text())
        if record.get('success') and (record['scenario'] not in out or record['total_score']<out[record['scenario']]['record']['total_score']):
            out[record['scenario']]={'path':path,'record':record,'z':record['z']}
    return out

def build_observation_rows(primary_fit):
    rows=[];record=primary_fit['record']
    with (CAMPAIGN/'data/hosny_observations.csv').open() as source_handle:
        hosny=list(csv.DictReader(source_handle))
    for study in STUDY_ORDER:
        for metric in METRIC_ORDER[study]:
            vals=[]
            for group in (0,1):
                candidates=[r for r in record['rows'] if r['study']==study and r['metric']==('msd_curve' if metric=='msd' else metric) and int(r['group'])==group]
                if not candidates: raise ValueError((study,metric,group))
                vals.append(max(candidates,key=lambda r:float(r.get('lag_h') or r.get('lag') or 0)))
            if metric=='msd':
                common_lag=min(float(r['lag_h']) for r in vals)
                for group in (0,1):
                    native=sorted([r for r in record['rows'] if r['study']==study and r['metric']=='msd_curve' and int(r['group'])==group],key=lambda r:r['lag_h'])
                    vals[group]=dict(vals[group],lag_h=common_lag,obs=float(np.interp(common_lag,[r['lag_h'] for r in native],[r['obs'] for r in native])),predicted=float(np.interp(common_lag,[r['lag_h'] for r in native],[r['predicted'] for r in native])))
            a,b=vals
            rows.append({'study':study,'metric':metric,'observed_WGDminus':a['obs'],'observed_WGDplus':b['obs'],'predicted_WGDminus':a['predicted'],'predicted_WGDplus':b['predicted'],'observed_fold_WGDplus_over_WGDminus':b['obs']/a['obs'],'predicted_fold_WGDplus_over_WGDminus':b['predicted']/a['predicted'],'summary_lag_WGDminus_h':a.get('lag_h') or '', 'summary_lag_WGDplus_h':b.get('lag_h') or '','fit_convention':'median' if study=='Hosny_baseline' else a.get('statistic','reviewed source central estimate'),'ratio_error_bars':False})
            for group,source_row in zip(('WGDminus','WGDplus'),vals):
                for field in ('statistic','low','high','raster_low','raster_high','sem_lower','sem_upper'):
                    rows[-1][group+'_'+field]=source_row.get(field,'')
                for q in ('q1','q3'):
                    hs=[r for r in hosny if r['metric']==metric and r['group']==('reference' if group=='WGDminus' else 'polyploid') and r['statistic'].lower()==q]
                    rows[-1][group+'_'+q]=hs[0]['value'] if study=='Hosny_baseline' and hs else ''
    return rows,{'note':'The common final overlapping-lag MSD fold is a visual summary interpolated from both native curve traces; the full native curves enter the balanced fit.'}

def write_observation_table():
    rows,meta=build_observation_rows(load_fit_scenarios()['primary'])
    path=FIGURE_DIR/'main_03_observation_folds.csv'
    with path.open('w',newline='') as handle:
        w=csv.DictWriter(handle,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    mp=path.with_suffix('.json');mp.write_text(json.dumps(meta,indent=2))
    return path,mp

def draw_observed_predicted(ax, rows: list[dict], title: str | None = None,
                            observed_color: str = "#5b7380",
                            predicted_color: str = "#d36c32"):
    """Draw all nine study/metric fold pairs on one shared logarithmic axis."""
    lookup = {(r["study"], r["metric"]): r for r in rows}
    order = [(study, metric) for study in STUDY_ORDER for metric in METRIC_ORDER[study]]
    y_positions = np.arange(len(order), dtype=float)
    for idx, (study, metric) in enumerate(order):
        r = lookup[(study, metric)]
        yobs, ypred = idx + .12, idx - .12
        observed = float(r["observed_fold_WGDplus_over_WGDminus"])
        predicted = float(r["predicted_fold_WGDplus_over_WGDminus"])
        ax.plot([observed, predicted], [yobs, ypred], color="#aebbc1", lw=.8, zorder=1)
        ax.scatter([observed], [yobs], marker="o", s=28, facecolor=observed_color,
                   edgecolor="white", linewidth=.45, zorder=3)
        ax.scatter([predicted], [ypred], marker="D", s=25, facecolor=predicted_color,
                   edgecolor="white", linewidth=.45, zorder=4)
    labels = [f"{STUDY_LABELS[s]} · {METRIC_LABELS[m]}" for s, m in order]
    ax.set_yticks(y_positions, labels)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.axvline(1., color="#33454e", linestyle="--", linewidth=.85, zorder=0)
    ax.set_xlabel("WGD+ / WGD− fold")
    if title:
        ax.set_title(title, loc="left", color=NAVY)
    ax.grid(axis="x", which="major", color="#dce4e8", linewidth=.6)
    ax.xaxis.set_major_locator(LogLocator(base=10, subs=(1., 2., 5.)))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{x:g}"))
    ax.tick_params(axis="y", labelsize=7.6)
    ax.tick_params(axis="x", labelsize=7.2)
    ax.set_ylim(len(order) - .4, -.6)
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([0], [0], marker="o", color="none", markerfacecolor=observed_color,
                              markeredgecolor="white", label="Observed fold", markersize=6),
                       Line2D([0], [0], marker="D", color="none", markerfacecolor=predicted_color,
                              markeredgecolor="white", label="Predicted fold", markersize=5.5),
                       Line2D([0], [0], color="#33454e", linestyle="--", label="Equal groups")],
              frameon=False, ncol=3, fontsize=7, loc="lower right")
    return {f"{s}:{m}": float(y) for y, (s, m) in zip(y_positions, order)}


def draw_relative_surface(ax, csv_path: Path | None = None,
                          score_column: str = "fixed_drift_delta_data",
                          norm: Normalize | None = None, add_colorbar: bool = False,
                          title: str | None = None, cmap: str = "RdBu_r"):
    """Draw the saved WGD−-relative score surface; performs no refit."""
    path = csv_path or RELATIVE_SURFACE_CSV
    with Path(path).open(newline="") as f:
        rows = list(csv.DictReader(f))
    rows = [r for r in rows if r.get("status") == "valid" and r.get(score_column, "") != ""]
    tau = np.asarray(sorted({float(r["tau_ratio"]) for r in rows}), dtype=float)
    sigma = np.asarray(sorted({float(r["sigma_ratio"]) for r in rows}), dtype=float)
    if len(rows) != len(tau) * len(sigma):
        raise ValueError(f"Saved relative surface is not a complete rectangular grid: {path}")
    grid = np.full((len(tau), len(sigma)), np.nan)
    ix = {v: i for i, v in enumerate(tau)}; iy = {v: i for i, v in enumerate(sigma)}
    for r in rows:
        grid[ix[float(r["tau_ratio"])],
             iy[float(r["sigma_ratio"])]] = float(r[score_column])
    if not np.isfinite(grid).all():
        raise ValueError(f"Relative surface contains non-finite displayed scores: {score_column}")
    if norm is None:
        lim = float(np.max(np.abs(grid)))
        norm = Normalize(-lim, lim)
    mesh = ax.pcolormesh(tau, sigma, grid.T, shading="auto", cmap=cmap, norm=norm, rasterized=True)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.scatter([1.], [1.], marker="o", s=37, facecolor="white", edgecolor="#242f35",
               linewidth=.8, zorder=5)
    fitpars=np.asarray(load_fit_scenarios()['primary']['z'][:20]).reshape(5,2,2)
    fitfold=np.exp((fitpars[:,1,:]-fitpars[:,0,:]).mean(axis=0))
    ax.scatter([fitfold[0]], [fitfold[1]], marker="*", s=105,
               facecolor="#111111", edgecolor="white", linewidth=.65, zorder=6)
    ax.set_xlabel("WGD+ / WGD− persistence fold")
    ax.set_ylabel("WGD+ / WGD− noise fold")
    if title:
        ax.set_title(title, loc="left", color=NAVY)
    ax.grid(which="major", color="white", linewidth=.55, alpha=.52)
    if add_colorbar:
        cb = ax.figure.colorbar(mesh, ax=ax, fraction=.045, pad=.035)
        cb.set_label(score_column.replace("_", " "))
    return mesh


