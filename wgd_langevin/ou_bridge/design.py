"""Frozen step 1–4 pilot design; no prospective experiment selection."""
from pathlib import Path
import hashlib
import json
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = HERE / 'results'
MODELS = ['wave_pinning', 'otsuji', 'goryachev', 'legi', 'debelly']
FAMILIES = ['gaussian', 'memory', 'coupled', 'self_activation']
TRAIN = {'default': [5., 15., 39.], 'debelly': [21., 29., 39.]}
TEST = {'default': [10., 25.], 'debelly': [25., 33.]}
ASSAYS = ['pulse', 'weak_plus', 'weak_minus', 'tracking_train',
          'tracking_test', 'reversal', 'competition']
CALIBRATION = ['pulse', 'weak_plus', 'weak_minus', 'tracking_train']
DEBELLY = dict(assay_mode='reflection_symmetric_extension',
               initialization_mode='algebraic_fixed_point', input_sampling_mode='cell_centred',
               rac_basal_activation=1., rac_mechanical_feedback=.2,
               rho_basal_activation=1., rho_mechanical_feedback=.2)

def lengths(model, train=False):
    key = 'debelly' if model == 'debelly' else 'default'
    return TRAIN[key] if train else sorted(TRAIN[key] + TEST[key])

def jobs():
    return [(m, L) for m in MODELS for L in lengths(m)]

def fit_jobs():
    return [(m, L, f) for m in MODELS for L in lengths(m, True) for f in FAMILIES]

def save_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temp.replace(path)

def digest(paths):
    return {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}

def clocks(model):
    # Seconds; dimensionless fitting clock prevents poorly scaled parameters.
    if model == 'otsuji': return (100., 20000., 40., .2, .1)
    return (1000., 10000., 160., 3.6, .2) if model == 'debelly' else (100., 4000., 40., .2, .1)

class Protocol:
    """Fixed peak, relative geometry, independent of L. No hidden dose/L effect."""
    def __init__(self, model, assay):
        self.model, self.assay = model, assay
        self.unit, self.condition, self.pulse, self.peak, self.width = clocks(model)
        self.period = self.unit * (4 if assay == 'tracking_test' else 1)

    @property
    def end(self):
        if self.assay in ('pulse', 'competition'):
            return self.condition + 10 * self.unit
        if self.assay.startswith('tracking'):
            return self.condition + 6 * self.period
        return self.condition + 10 * self.unit

    @property
    def event_times(self):
        return (0., self.pulse/2, self.pulse, self.condition,
                self.condition + self.unit)

    def __call__(self, t, x, length):
        from polarity1d.protocols import LocalizedPulse
        if self.assay in ('pulse', 'competition'):
            return LocalizedPulse(amplitude=self.peak, start=0., stop=self.pulse,
                                  taper=self.pulse/2, width_fraction=self.width,
                                  sides=('left', 'right') if self.assay == 'competition' else ('left',),
                                  asymmetry=.02 if self.assay == 'competition' else 0.)(t, x, length)
        # Nonnegative tonic gradient. Weak probes modulate its slope by +/-5%,
        # holding its spatial mean fixed; they probe a cue-maintained state.
        orientation = -1.
        if t >= self.condition:
            if self.assay.startswith('tracking'):
                orientation = -np.cos(2*np.pi*(t-self.condition)/self.period)
            elif self.assay == 'reversal':
                orientation = 1.
            elif t < self.condition + self.unit:
                orientation += .05 if self.assay == 'weak_plus' else -.05
        # 0.8 leaves room for weak modulation without negative stimulus.
        return self.peak/2 * (1. + .8*orientation*(2*x/length-1))

    def projected(self, t):
        # Exact continuous right-half minus left-half means.
        if self.assay in ('pulse', 'competition'):
            if t < 0 or t >= self.pulse:
                return 0.
            envelope = 1. if t <= self.pulse/2 else .5*(1+np.cos(np.pi*(t-self.pulse/2)/(self.pulse/2)))
            return -self.peak*self.width*envelope*(.02 if self.assay == 'competition' else 1.)
        direction = -1.
        if t >= self.condition:
            if self.assay.startswith('tracking'):
                direction = -np.cos(2*np.pi*(t-self.condition)/self.period)
            elif self.assay == 'reversal':
                direction = 1.
            elif t < self.condition+self.unit:
                direction += .05 if self.assay == 'weak_plus' else -.05
        return .4*self.peak*direction

def manifest():
    return dict(version=1, models=MODELS, families=FAMILIES, train=TRAIN, test=TEST,
                calibration=CALIBRATION, assays=ASSAYS, starts=8, seed=20260925,
                source_hashes=digest([HERE/'design.py']),
                observable='interpolated half-mass/L - 0.5; zero-mass profile maps to zero, flagged',
                cue='continuous right-half minus left-half mean; input enters p only',
                max_dx_um=.2, stationary_noise='not inferred from deterministic spatial residuals')
