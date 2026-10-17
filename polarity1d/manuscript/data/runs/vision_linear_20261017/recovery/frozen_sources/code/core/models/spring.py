from __future__ import annotations
import numpy as np
import sys
from typing import Callable, Any
from scipy.integrate import solve_ivp
AMP=4.0
WIDTH=0.2
PREP=1500.0

def f(u, v, b, gamma, n, RT, delta):
    return (b+gamma*u**n/(1+u**n))*v - delta*u


def F(R):
    sharp = 10
    switch = 1
    magnitude = 0.001
    return magnitude/(1+np.exp(-2*sharp*(R-switch)))


auth={"f":f,"F":F}
def rhs(t,y,L0=1.,d1=80.,b=4.,cue=0.,side=1,author_dx=False):
    u=y[:-2:2];v=y[1:-2:2];n=len(u);x=np.linspace(0,1,n)
    length=y[-1]-y[-2];extension=length-L0
    left=.01*extension-auth['F'](u[0]);right=-.01*extension+auth['F'](u[-1])
    rate=(right-left)/length
    activation=b+cue*np.exp(-((x-(1 if side==1 else 0))/.2)**2)
    delta=3+d1*extension
    reaction=auth['f'](u,v,activation,5,6,2,delta)
    dx=1/n if author_dx else 1/(n-1)
    out=np.empty_like(y)
    for a,D,j,sign in [(u,.01,0,1),(v,10,1,-1)]:
        lap=np.r_[2*(a[1]-a[0]),np.diff(a,2),2*(a[-2]-a[-1])]/dx**2
        out[j:-2:2]=sign*reaction+D/length**2*lap-a*rate
    out[-2:]=left,right
    return out


def initial(n,L0,perturb=1e-3):
    x=np.linspace(0,1,n);u=.05645+perturb*np.cos(np.pi*x)
    y=np.empty(2*n+2);y[:-2:2]=u;y[1:-2:2]=2-u;y[-2:]=0,L0
    return y


def metrics(y):
    u=y[:,:-2:2];v=y[:,1:-2:2];length=y[:,-1]-y[:,-2]
    x=np.linspace(0,1,u.shape[1]);mass=np.trapezoid(u+v,x,axis=1)*length
    polarity=np.trapezoid(u*(2*x-1),x,axis=1)/np.trapezoid(u,x,axis=1)
    return polarity,mass,length

class _CallableStimulus:
    def __init__(self, function: Callable):
        self.function = function
        self.event_times = tuple(getattr(function, "event_times", ()))

    def __call__(self, time: float, x: np.ndarray, length: float) -> np.ndarray:
        return self.function(time, x, length)


def _as_stimulus(stimulus: Callable | None) -> Callable:
    if stimulus is None:
        return _CallableStimulus(lambda t, x, length: np.zeros_like(x, dtype=float))
    if not callable(stimulus):
        raise TypeError("stimulus must be a callable stimulus(t, x, length)")
    return _CallableStimulus(stimulus)


def run_spring(
    *,
    length: float,
    stimulus: Callable | None,
    end_time: float,
    output_interval: float = 5.0,
    cells: int = 201,
    seed: int = 1,
    initial_state: np.ndarray | None = None,
    rtol: float = 2e-7,
    atol: float = 2e-9,
) -> dict[str, Any]:
    """Run the Zmurchok et al. author model at native resting length ``L0``.

    Length is the author's resting-length parameter (their Figure 4 default is
    L0=1), not microns. Time is uncalibrated author model time. The source's
    stiffness/drag/reaction coefficients, basal activation b=4, d1=80, and
    diffusivities are unchanged. Only the assay's cue profile varies in time.
    """
    if length <= 0 or cells < 5 or end_time <= 0 or output_interval <= 0:
        raise ValueError("length, cells, end_time, and output_interval must be positive")
    spring = sys.modules[__name__]
    cue = _as_stimulus(stimulus)
    n = int(cells)
    xfrac = np.linspace(0.0, 1.0, n)
    state0 = spring.initial(n, float(length), perturb=1e-3) if initial_state is None else np.array(initial_state, dtype=float, copy=True)
    expected = 2 * n + 2
    if state0.shape != (expected,):
        raise ValueError(f"spring initial_state must have shape ({expected},)")
    times = np.arange(0.0, end_time, output_interval, dtype=float)
    if len(times) == 0 or not np.isclose(times[-1], end_time):
        times = np.r_[times, float(end_time)]
    else:
        times[-1] = float(end_time)

    def rhs(t, y):
        # The author equations use normalized material coordinates. Map those
        # positions onto the instantaneous physical domain for the shared cue.
        domain_length = float(y[-1] - y[-2])
        x = xfrac * domain_length
        profile = np.asarray(cue(float(t), x, domain_length), dtype=float)
        if profile.ndim == 0:
            profile = np.full(n, float(profile))
        if profile.shape != (n,):
            raise ValueError(f"spring stimulus returned {profile.shape}, expected ({n},)")
        out = spring.rhs(t, y, L0=float(length), d1=80.0, b=4.0, cue=0.0)
        if np.any(profile):
            u = y[:-2:2]
            v = y[1:-2:2]
            extension = domain_length - float(length)
            delta = 3.0 + 80.0 * extension
            transfer = spring.auth["f"](u, v, 4.0 + profile, 5, 6, 2, delta) - spring.auth["f"](u, v, 4.0, 5, 6, 2, delta)
            out[:-2:2] += transfer
            out[1:-2:2] -= transfer
        return out

    from scipy.integrate import solve_ivp
    # Stop the adaptive solver exactly at every cue discontinuity/cross-fade
    # boundary. This matters for square-wave periodic following and staged tests.
    boundaries = [0.0]
    boundaries.extend(sorted({float(event) for event in getattr(cue, "event_times", ())
                              if 0.0 < float(event) < float(end_time)}))
    boundaries.append(float(end_time))
    out_t = [0.0]
    out_y = [state0]
    state = state0
    for start, stop in zip(boundaries[:-1], boundaries[1:]):
        sampled = times[(times > start + 1e-12) & (times < stop - 1e-12)]
        segment_times = np.r_[start, sampled, stop]
        sol = solve_ivp(rhs, (start, stop), state, method="LSODA", t_eval=segment_times,
                        rtol=rtol, atol=atol, max_step=min(5.0, output_interval))
        if not sol.success:
            raise RuntimeError(f"spring integration failed: {sol.message}")
        out_t.extend(sol.t[1:].tolist())
        out_y.extend(sol.y.T[1:])
        state = sol.y[:, -1].copy()
    out_t = np.asarray(out_t, dtype=float)
    y = np.asarray(out_y, dtype=float)
    active, inactive = y[:, :-2:2], y[:, 1:-2:2]
    domain_length = y[:, -1] - y[:, -2]
    # Author polarity is the normalized first moment on material coordinate.
    polarity, mass, _ = spring.metrics(y)
    xphysical = xfrac[None, :] * domain_length[:, None]
    left_peak = np.max(active[:, : max(1, n // 3)], axis=1)
    right_peak = np.max(active[:, -(max(1, n // 3)):], axis=1)
    baseline = np.min(active, axis=1)
    return {
        "model": "spring_zmurchok2020",
        "time": out_t,
        "x_fraction": xfrac,
        "x_physical": xphysical,
        "active": active,
        "inactive": inactive,
        "orientation": np.asarray(polarity),
        "length": domain_length,
        "left_peak_excess": left_peak - baseline,
        "right_peak_excess": right_peak - baseline,
        "mass": np.asarray(mass),
        "mass_relative_error": float(np.max(np.abs(np.asarray(mass) / float(mass[0]) - 1.0))),
        "native_length": float(length),
        "native_length_units": "author model length (L0=1 is source default)",
        "native_time_units": "author model time (uncalibrated)",
        "parameter_defaults": {"L0": float(length), "reference_L0": 1.0, "basal_activation_b": 4.0,
                               "tension_feedback_d1": 80.0, "diffusion_active": 0.01,
                               "diffusion_inactive": 10.0, "rate_parameters": [5, 6, 2, 3],
                               "perturbation": 1e-3, "source": "Zmurchok et al. (2020), author Fig4c.py"},
    }

def gaussian(side: str, amplitude: float = AMP, width: float = WIDTH):
    center = 1.0 if side == "right" else 0.0
    return lambda q: float(amplitude) * np.exp(-((q - center) / float(width)) ** 2)


def _run(length, state, cue, duration, *, sample=5.0, cells=201):
    return run_spring(length=float(length), initial_state=state, stimulus=cue,
                      end_time=float(duration), output_interval=float(sample), cells=int(cells))


def _end(result):
    active, inactive = result["active"][-1], result["inactive"][-1]
    state = np.empty(2 * len(active) + 2, dtype=float)
    state[:-2:2] = active
    state[1:-2:2] = inactive
    state[-2:] = 0.0, result["length"][-1]
    return state


def _initial_left(length, cells):
    """Use the model's default perturbation, then orient the spontaneous front left."""
    initial = run_spring(length=length, stimulus=None, end_time=PREP,
                         output_interval=5.0, cells=cells)
    state = _end(initial)
    if float(initial["orientation"][-1]) > 0:
        n = int(cells)
        active = state[:-2:2].copy()
        inactive = state[1:-2:2].copy()
        state[:-2:2] = active[::-1]
        state[1:-2:2] = inactive[::-1]
        # Reflect the two moving boundaries around the resting length.
        edges = state[-2:].copy()
        state[-2:] = float(length) - edges[1], float(length) - edges[0]
    return state, initial

