#!/usr/bin/env python3
"""Generate the small, versioned dataset for the 1D De Belly mechanics report."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from debelly_1d.model import (  # noqa: E402
    active_stress_response,
    membrane_pull,
    variable_friction_response,
)


OUT = ROOT / "reports" / "assets" / "debelly_1d" / "v1"
BLUE = "#2166ac"
RED = "#b2182b"
GOLD = "#d89400"
GREEN = "#1b7837"
PURPLE = "#762a83"


def style():
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.titleweight": "bold",
            "figure.dpi": 150,
            "savefig.dpi": 180,
            "savefig.bbox": "tight",
        }
    )


def save(fig, name):
    fig.savefig(OUT / name, facecolor="white")
    plt.close(fig)


def figure_diffusion(rows):
    times = np.array([1.0, 5.0, 20.0, 100.0])
    result = membrane_pull(mu=2.0, times=times, length=10.0, n=321)
    fig, axes = plt.subplots(1, 3, figsize=(12.0, 3.25))
    colors = plt.cm.viridis(np.linspace(0.12, 0.9, len(times)))
    for i, (time, color) in enumerate(zip(times, colors)):
        axes[0].plot(result["s"], result["displacement"][i], color=color, label=f"{time:g} s")
        axes[1].plot(result["s"], -result["tension"][i], color=color)
        normalized = result["tension"][i] / result["tension"][i, 0]
        scaled_s = result["s"] / np.sqrt(result["diffusivity"] * time)
        use = scaled_s <= 5
        # The similarity solution applies before the front senses the distal
        # boundary. Later curves are retained in A/B to show finite-size
        # equilibration, but deliberately excluded from this collapse panel.
        if time <= 5:
            axes[2].plot(scaled_s[use], normalized[use], color=color, label=f"{time:g} s")
        for s, disp, tension, zeta, norm_t in zip(
            result["s"], result["displacement"][i], result["tension"][i], scaled_s, normalized
        ):
            rows.append(
                {
                    "figure": "diffusion",
                    "case": f"t={time:g}",
                    "time": time,
                    "position": s,
                    "scaled_position": zeta,
                    "displacement": disp,
                    "tension": tension,
                    "normalized_tension": norm_t,
                }
            )
    axes[0].set(title="A  Step pull spreads", xlabel="position s (µm)", ylabel="displacement x")
    axes[1].set(title="B  Tension front", xlabel="position s (µm)", ylabel="tension magnitude −T")
    similarity_x = np.linspace(0, 5, 200)
    axes[2].plot(similarity_x, np.exp(-similarity_x**2 / 4), "--", color="0.25", label="semi-infinite limit")
    axes[2].set(title="C  Early-time collapse", xlabel=r"$s/\sqrt{(k/\mu)t}$", ylabel=r"$T(s,t)/T(0,t)$", xlim=(0, 5))
    axes[2].legend(frameon=False)
    axes[0].legend(frameon=False, title="time")
    fig.suptitle(r"Membrane-only limit: $x_t=(k/\mu)x_{ss}$ ($k/\mu=0.5\,\mu m^2/s$)", y=1.03)
    fig.tight_layout()
    save(fig, "01_tension_diffusion.png")


def figure_friction(rows, summary):
    mus = np.logspace(-1, 2, 25)
    exogenous = []
    endogenous = []
    for mu in mus:
        pull = membrane_pull(mu=mu, times=[0.5, 20.0], length=10.0, n=181)
        tension = np.abs(pull["tension"][-1])
        proximal = np.mean(tension[1:4])
        distal = np.mean(tension[-4:-1])
        exogenous.append(distal / proximal)
        active = active_stress_response(mu=mu, times=[0.5, 20.0], length=10.0, n=256)
        endogenous.append(np.sqrt(np.mean(active["tension"][-1] ** 2)))
    exogenous = np.asarray(exogenous)
    endogenous = np.asarray(endogenous)
    endogenous /= endogenous.max()

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.45))
    axes[0].semilogx(mus, exogenous, "o-", color=BLUE, ms=4, label="membrane step pull")
    axes[0].set(xlabel="friction µ", ylabel="distal / proximal tension", title="A  Exogenous membrane loading")
    axes[0].annotate("rapid redistribution", xy=(mus[1], exogenous[1]), xytext=(0.12, 0.6),
                     arrowprops={"arrowstyle": "->", "color": BLUE}, color=BLUE)
    axes[1].semilogx(mus, endogenous, "o-", color=RED, ms=4, label="cortical active stress")
    axes[1].set(xlabel="friction µ", ylabel="normalized membrane tension RMS", title="B  Endogenous cortical loading")
    axes[1].annotate("stronger loading", xy=(mus[-3], endogenous[-3]), xytext=(1.0, 0.45),
                     arrowprops={"arrowstyle": "->", "color": RED}, color=RED)
    for ax in axes:
        ax.set_ylim(bottom=0)
        ax.grid(alpha=0.2)
    fig.suptitle("Friction has opposite roles because the force is applied to a different layer", y=1.04)
    fig.tight_layout()
    save(fig, "02_opposite_friction.png")

    for mu, exo, endo in zip(mus, exogenous, endogenous):
        rows.extend(
            [
                {"figure": "friction", "case": "membrane_pull", "mu": mu, "response": exo},
                {"figure": "friction", "case": "cortex_stress", "mu": mu, "response": endo},
            ]
        )
    summary.extend(
        [
            {"quantity": "exogenous distal/proximal at mu=0.1", "value": exogenous[0]},
            {"quantity": "exogenous distal/proximal at mu=100", "value": exogenous[-1]},
            {"quantity": "endogenous normalized RMS at mu=0.1", "value": endogenous[0]},
            {"quantity": "endogenous normalized RMS at mu=100", "value": endogenous[-1]},
        ]
    )


def figure_force_chain(rows, summary):
    result = active_stress_response(
        mu=50.0,
        times=[0.5, 50.0, 100.0, 150.0, 250.0],
        length=10.0,
        pulse_stop=100.0,
        n=512,
    )
    index = 1
    fields = [
        ("cortex_velocity", r"cortex $v$", RED),
        ("membrane_velocity", r"membrane $x_t$", BLUE),
        ("slip", r"slip $v-x_t$", PURPLE),
        ("friction", r"friction $\mu(v-x_t)$", GOLD),
        ("tension", r"tension $T$", GREEN),
    ]
    fig, axes = plt.subplots(2, 1, figsize=(8.8, 6.1), sharex=True)
    for key, label, color in fields:
        target = axes[0] if key in {"cortex_velocity", "membrane_velocity", "slip"} else axes[1]
        target.plot(result["s"], result[key][index], color=color, lw=2, label=label)
        for s, value in zip(result["s"], result[key][index]):
            rows.append(
                {"figure": "force_chain", "case": key, "time": result["time"][index], "position": s, "value": value}
            )
    axes[0].set(ylabel="velocity", title="A  Motion and relative slip")
    axes[1].set(xlabel="position s (µm)", ylabel="force or tension", title="B  Friction loads the elastic membrane")
    for ax in axes:
        ax.axvspan(0.95, 2.05, color="0.9", zorder=-2, label="active-stress zone" if ax is axes[0] else None)
        ax.axhline(0, color="0.6", lw=0.7)
        ax.legend(frameon=False, ncol=3)
    fig.suptitle("A localized cortical contraction at t = 50 s", y=1.01)
    fig.tight_layout()
    save(fig, "03_force_chain.png")

    ds = result["s"][1] - result["s"][0]
    q = 2 * np.pi * np.fft.fftfreq(len(result["s"]), d=ds)
    tension_gradient = np.fft.ifft(1j * q * np.fft.fft(result["tension"][index])).real
    residual = tension_gradient + result["friction"][index]
    summary.append({"quantity": "max membrane force-balance residual", "value": np.max(np.abs(residual))})

    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.3), sharex=True)
    for time_index, color in zip([1, 2, 3, 4], [RED, GOLD, BLUE, PURPLE]):
        label = f"{result['time'][time_index]:g} s"
        axes[0].plot(result["s"], result["cortex_velocity"][time_index], color=color, label=label)
        axes[1].plot(result["s"], result["membrane_velocity"][time_index], color=color)
        axes[2].plot(result["s"], result["tension"][time_index], color=color)
    axes[0].set(title="A  Cortical velocity", ylabel="v")
    axes[1].set(title="B  Membrane velocity", ylabel=r"$x_t$")
    axes[2].set(title="C  Membrane tension", ylabel="T")
    for ax in axes:
        ax.set_xlabel("position s (µm)")
        ax.axhline(0, color="0.7", lw=0.7)
    axes[0].legend(frameon=False, title="pulse ends at 100 s")
    fig.suptitle("Pulse and recovery: motion stops first, elastic stress relaxes afterward", y=1.04)
    fig.tight_layout()
    save(fig, "04_pulse_recovery.png")


def figure_modes(rows, summary):
    length, stiffness, viscosity, mu = 10.0, 1.0, 100.0, 4.0
    tau = viscosity / stiffness
    lam = np.sqrt(viscosity / mu)
    times = np.linspace(0, 180, 181)
    modes = np.arange(1, 6)
    q = 2 * np.pi * modes / length
    rates = (1 + lam**2 * q**2) / tau
    amplitudes = np.exp(-times[:, None] * rates[None, :])
    fit = np.polyfit(q**2, rates, 1)
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.5))
    for mode, rate, amp in zip(modes, rates, amplitudes.T):
        axes[0].semilogy(times, amp, lw=2, label=f"mode {mode}")
        for time, value in zip(times, amp):
            rows.append(
                {"figure": "modes", "case": f"mode_{mode}", "time": time, "q2": q[mode - 1] ** 2, "amplitude": value, "decay_rate": rate}
            )
    axes[0].set(xlabel="time (s)", ylabel="amplitude A(t)", title="A  Fine modes disappear first", ylim=(1e-7, 1.1))
    axes[0].legend(frameon=False, ncol=2)
    axes[1].plot(q**2, rates, "o", color=PURPLE, ms=7)
    q2line = np.linspace(0, 1.05 * q[-1] ** 2, 100)
    axes[1].plot(q2line, fit[0] * q2line + fit[1], color="0.25")
    axes[1].set(xlabel=r"wave number squared $q^2$", ylabel=r"decay rate $-s$ (s$^{-1}$)", title="B  Mechanical low-pass filter")
    axes[1].text(0.04, 0.94, r"intercept $=1/\tau$" + "\n" + r"slope $=\lambda^2/\tau$", transform=axes[1].transAxes, va="top")
    fig.suptitle(r"Fourier prediction: $A(t)=A_0\exp[-(1+\lambda^2q^2)t/\tau]$", y=1.04)
    fig.tight_layout()
    save(fig, "05_wavelength_filter.png")
    summary.extend(
        [
            {"quantity": "mode-fit intercept", "value": fit[1]},
            {"quantity": "predicted 1/tau", "value": 1 / tau},
            {"quantity": "mode-fit slope", "value": fit[0]},
            {"quantity": "predicted lambda^2/tau", "value": lam**2 / tau},
        ]
    )


def endpoint_transmission(length_over_lambda, time_over_tau, modes=700):
    """Exact series for a step pull on a finite interval."""
    a = time_over_tau / length_over_lambda**2
    n = np.arange(1, modes + 1, dtype=float)
    decay = np.exp(-(np.pi * n) ** 2 * a)
    proximal = 1.0 + 2.0 * decay.sum()
    distal = 1.0 + 2.0 * np.sum((-1.0) ** n * decay)
    return np.clip(distal / proximal, 0.0, 1.0)


def figure_regime(rows):
    xvals = np.logspace(-1, 1.5, 120)
    yvals = np.logspace(-3, 1, 120)
    transmission = np.array([[endpoint_transmission(x, y) for x in xvals] for y in yvals])
    fig, ax = plt.subplots(figsize=(7.1, 4.7))
    mesh = ax.pcolormesh(xvals, yvals, transmission, shading="auto", cmap="viridis", vmin=0, vmax=1)
    contours = ax.contour(xvals, yvals, transmission, levels=[0.1, 0.5, 0.9], colors="white", linewidths=0.8)
    ax.clabel(contours, fmt="%.1f", fontsize=8)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set(xlabel=r"cell size / hydrodynamic length $L/\lambda$", ylabel=r"observation time $t_{obs}/\tau$", title="Distal / proximal tension after a membrane step pull")
    colorbar = fig.colorbar(mesh, ax=ax, label="transmission")
    colorbar.outline.set_visible(False)
    ax.text(0.14, 3.0, "global / equilibrated", color="white", weight="bold")
    ax.text(8.0, 0.002, "apparently local", color="white", weight="bold", ha="center")
    fig.tight_layout()
    save(fig, "06_regime_map.png")
    for yi in range(0, len(yvals), 3):
        for xi in range(0, len(xvals), 3):
            rows.append(
                {"figure": "regime", "case": "step_pull", "L_over_lambda": xvals[xi], "t_over_tau": yvals[yi], "transmission": transmission[yi, xi]}
            )


def figure_mca(rows):
    n = 241
    s = np.linspace(0, 10, n)
    base_mu = 5.0
    profiles = {
        "uniform MCA": np.full(n, base_mu),
        "front depletion": base_mu * (1 - 0.8 * np.exp(-0.5 * ((s - 1.2) / 0.55) ** 2)),
        "middle depletion": base_mu * (1 - 0.8 * np.exp(-0.5 * ((s - 5.0) / 0.55) ** 2)),
    }
    colors = {"uniform MCA": "0.25", "front depletion": RED, "middle depletion": BLUE}
    times = [0.5, 20.0]
    cortex_velocity = np.exp(-0.5 * ((s - 1.2) / 0.55) ** 2)
    cortex_velocity[[0, -1]] = 0
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.35))
    for label, profile in profiles.items():
        pull = variable_friction_response(mu_profile=profile, times=times, pull=1.0)
        drive = variable_friction_response(mu_profile=profile, times=times, cortex_velocity=cortex_velocity)
        axes[0].plot(s, profile / base_mu, color=colors[label], lw=2, label=label)
        axes[1].plot(s, -pull["tension"][-1], color=colors[label], lw=2)
        axes[2].plot(s, drive["tension"][-1], color=colors[label], lw=2)
        for position, relative_mu, pull_t, drive_t in zip(s, profile / base_mu, pull["tension"][-1], drive["tension"][-1]):
            rows.append(
                {"figure": "mca", "case": label, "position": position, "relative_mu": relative_mu, "membrane_pull_tension": pull_t, "cortex_drive_tension": drive_t}
            )
    axes[0].set(title=r"A  $\mu(s)=\chi\rho_b(s)$", ylabel=r"relative friction $\mu/\mu_0$")
    axes[1].set(title="B  Stress already in membrane", ylabel="tension magnitude −T")
    axes[2].set(title="C  Prescribed cortical motion", ylabel="membrane tension T")
    for ax in axes:
        ax.set_xlabel("position s (µm)")
        ax.axhline(0, color="0.75", lw=0.6)
    axes[0].legend(frameon=False)
    fig.suptitle("An MCA-depleted patch is a local slip zone; its effect depends on where loading begins", y=1.04)
    fig.tight_layout()
    save(fig, "07_mca_friction.png")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    style()
    rows = []
    summary = []
    figure_diffusion(rows)
    figure_friction(rows, summary)
    figure_force_chain(rows, summary)
    figure_modes(rows, summary)
    figure_regime(rows)
    figure_mca(rows)
    pd.DataFrame(rows).to_csv(OUT / "plot_data.csv", index=False)
    pd.DataFrame(summary).to_csv(OUT / "validation_summary.csv", index=False)
    metadata = {
        "schema_version": 1,
        "generator": "scripts/generate_debelly_1d_report.py",
        "model_scope": "linear 1D mechanics only; no Rac/Rho kinetics",
        "coordinate_units": "micrometres",
        "time_units": "seconds",
        "figures": sorted(path.name for path in OUT.glob("*.png")),
        "animations": sorted(path.name for path in OUT.glob("*.gif")),
    }
    (OUT / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"Wrote {len(rows):,} plot rows and {len(summary)} validation values to {OUT}")


if __name__ == "__main__":
    main()
