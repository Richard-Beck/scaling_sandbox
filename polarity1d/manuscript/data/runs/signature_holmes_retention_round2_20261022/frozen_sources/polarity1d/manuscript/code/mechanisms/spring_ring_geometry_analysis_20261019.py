"""Read-only geometry/dilution accounting of the verified fresh Spring ring.

No trajectories or ablations are generated. Exact geometric integrals use
circumference endpoint ratios; sampled quadrature is retained separately and
must not replace that identity. Phase archives and frozen executable sources
are hash checked before analysis.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import sys

sys.dont_write_bytecode = True
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DEFAULT = ROOT / "data/runs/spring_ring_validation_20261019"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main(run):
    manifest = json.loads((run / "manifest.json").read_text())
    for relative, digest in manifest["source_sha256"].items():
        assert sha(run / "frozen_sources" / relative) == digest
    sys.path.insert(0, str(run / "frozen_sources/code/mechanisms"))
    from spring_ring_closure_20261017 import Closure, native

    record_path = run / "workers/00/status.json"
    record = json.loads(record_path.read_text())
    assert record["complete"] and record["status"] == "complete"
    assert record["job"]["name"] == "fresh_200"
    closure = Closure(manifest["resting_circumference"], record["job"]["cells"])
    n = closure.cells
    reference_C = closure.resting_circumference
    rows, summaries, inputs = [], [], {"manifest.json": sha(run / "manifest.json"),
        "workers/00/status.json": sha(record_path)}
    elapsed = 0.0
    reference_amount = None
    previous = None

    for phase in record["phases"]:
        archive = run / phase["path"]
        assert sha(archive) == phase["sha256"]
        inputs[phase["path"]] = phase["sha256"]
        with np.load(archive) as z:
            times, states = z["time"].copy(), z["states"].copy()
        assert np.isfinite(states).all() and (np.diff(times) > 0).all()
        if previous is not None:
            np.testing.assert_array_equal(states[0], previous)
        previous = states[-1].copy()
        u, v, circumference = states[:, :n], states[:, n:-1], states[:, -1]
        mean_force = np.mean(native.auth["F"](u), axis=1)
        speed = 2.0 * (mean_force - 0.01 * (circumference - reference_C))
        alpha = speed / circumference
        active_dilution, inactive_dilution = -u*alpha[:, None], -v*alpha[:, None]
        # Independently validate the frozen closure's circumference derivative.
        for index in np.unique(np.r_[0, len(states)//2, len(states)-1, np.argmax(alpha), np.argmin(alpha)]):
            np.testing.assert_allclose(closure.rhs(states[index], 0.0)[-1], speed[index], atol=1e-18, rtol=1e-13)
        mean_active, mean_inactive = np.mean(u, axis=1), np.mean(v, axis=1)
        mean_total = mean_active + mean_inactive
        amount = circumference * mean_total
        if reference_amount is None:
            reference_amount = float(amount[0])
        assert np.max(np.abs(amount/reference_amount - 1)) < manifest["amount_relative_tolerance"]
        integrated = np.log(circumference/circumference[0])
        relative_C = circumference/circumference[0]
        rest_relative_C = circumference/reference_C
        delta = 3.0 + 80.0*(circumference - reference_C)
        for index, time in enumerate(times):
            rows.append(dict(phase=phase["kind"], phase_time=float(time),
                protocol_time=float(elapsed + time), circumference=float(circumference[index]),
                circumference_rate=float(speed[index]), fractional_dilation_rate=float(alpha[index]),
                signed_dilation_integral_since_phase_start=float(integrated[index]),
                geometric_concentration_factor_since_phase_start=float(1/relative_C[index]),
                geometric_diffusion_factor_since_phase_start=float(relative_C[index]**-2),
                geometric_concentration_factor_from_C0=float(1/rest_relative_C[index]),
                geometric_diffusion_factor_from_C0=float(rest_relative_C[index]**-2),
                deactivation_delta=float(delta[index]), mean_active=float(mean_active[index]),
                mean_inactive=float(mean_inactive[index]), mean_total_concentration=float(mean_total[index]),
                physical_total_amount=float(amount[index]),
                mean_active_dilution_rate=float(np.mean(active_dilution[index])),
                mean_inactive_dilution_rate=float(np.mean(inactive_dilution[index])),
                local_active_dilution_rate_min=float(np.min(active_dilution[index])),
                local_active_dilution_rate_max=float(np.max(active_dilution[index])),
                local_inactive_dilution_rate_min=float(np.min(inactive_dilution[index])),
                local_inactive_dilution_rate_max=float(np.max(inactive_dilution[index]))))

        summaries.append(dict(phase=phase["kind"], duration=float(times[-1]),
            initial_circumference=float(circumference[0]), final_circumference=float(circumference[-1]),
            circumference_range=[float(circumference.min()), float(circumference.max())],
            maximum_absolute_C_excursion_from_phase_start=float(np.max(np.abs(circumference-circumference[0]))),
            maximum_fractional_C_excursion_from_phase_start=float(np.max(np.abs(relative_C-1))),
            fractional_dilation_rate_range=[float(alpha.min()), float(alpha.max())],
            fractional_dilation_rate_min_phase_time=float(times[np.argmin(alpha)]),
            fractional_dilation_rate_max_phase_time=float(times[np.argmax(alpha)]),
            exact_signed_dilation_integral=float(integrated[-1]),
            sampled_trapezoid_dilation_integral=float(np.trapezoid(alpha, times)),
            sampled_quadrature_minus_exact=float(np.trapezoid(alpha, times)-integrated[-1]),
            geometric_concentration_factor_end=float(1/relative_C[-1]),
            geometric_concentration_percent_change_end=float(100*(1/relative_C[-1]-1)),
            geometric_diffusion_factor_end=float(relative_C[-1]**-2),
            geometric_diffusion_percent_change_end=float(100*(relative_C[-1]**-2-1)),
            deactivation_delta_range=[float(delta.min()), float(delta.max())],
            local_active_dilution_rate_range=[float(active_dilution.min()), float(active_dilution.max())],
            local_inactive_dilution_rate_range=[float(inactive_dilution.min()), float(inactive_dilution.max())],
            mean_total_concentration_start=float(mean_total[0]),
            mean_total_concentration_end=float(mean_total[-1]),
            physical_total_amount_range=[float(amount.min()), float(amount.max())],
            max_physical_total_amount_relative_error=float(np.max(np.abs(amount/reference_amount-1)))))
        elapsed += float(times[-1])

    basal = next(item for item in summaries if item["phase"] == "basal")
    prepared_C = basal["final_circumference"]
    prepared = dict(initial_circumference=reference_C, prepared_circumference=prepared_C,
        fractional_expansion=prepared_C/reference_C-1,
        exact_signed_dilation_integral=float(np.log(prepared_C/reference_C)),
        geometric_concentration_factor=reference_C/prepared_C,
        geometric_concentration_percent_change=100*(reference_C/prepared_C-1),
        geometric_diffusion_factor=(reference_C/prepared_C)**2,
        geometric_diffusion_percent_change=100*((reference_C/prepared_C)**2-1),
        deactivation_delta_initial=3.0, deactivation_delta_prepared=3+80*(prepared_C-reference_C),
        mean_total_concentration_initial=reference_amount/reference_C,
        mean_total_concentration_prepared=basal["mean_total_concentration_end"],
        conserved_physical_total_amount=reference_amount)

    authors = ROOT.parents[1] / "reports/trial_report_assets/zmurchok2020/authors/Fig4/Fig4c.py"
    active = ROOT / "code/core/models/spring.py"
    result = dict(generated_at=datetime.now(timezone.utc).isoformat(), passed=True,
        trajectory_integrations=0, scientific_scope="Read-only geometry/dilution accounting of one audited C0=2 fresh N200 protocol",
        input_sha256=inputs,
        analysis_source_sha256=sha(Path(__file__)),
        source_evidence=[dict(path=str(authors), sha256=sha(authors),
                lines="25–33 endpoint velocities; 48 total length; 53 affine physical coordinate; 55 global deactivation; 65–71 global diffusion and dilution"),
            dict(path=str(active), sha256=sha(active),
                lines="22–35 native scalar length rate and global diffusion/dilution; 92–107 normalized-to-physical coordinate")],
        user_decision_status="User expressed conditional comfort if native interval already applies endpoint expansion/dilution globally; that condition is verified by source. No blanket approval of all ring mechanics is inferred.",
        prepared_geometry=prepared, phases=summaries,
        derivation=dict(
            affine_map_latex=r"x(q,t)=x_-(t)+qL(t),\quad q\in[0,1]",
            velocity_latex=r"w(x(q,t),t)=\dot x_-+q\dot L,\qquad \partial_x w=\dot L/L",
            concentration_latex=r"\left.\partial_t c\right|_q=\frac{D}{L^2}\partial_q^2c+R_c-c\frac{\dot L}{L}",
            geometric_integral_latex=r"\int_{t_a}^{t_b}\frac{\dot C}{C}\,dt=\log\frac{C(t_b)}{C(t_a)},\qquad c_{\mathrm{geometry}}(t_b)/c(t_a)=C(t_a)/C(t_b)",
            amount_latex=r"M=C\int_0^1(u+v)\,dq,\qquad \dot M=0",
            meaning="The affine material-coordinate treatment is implied by the source chemical equations. It is not a solved local elasticity field. Endpoint asymmetry changes translation as well as total length; translation does not remove the uniform dilation term.",
            ring_carryover="Periodic chemistry carries the global length-scale diffusion, dilution and extension-dependent deactivation into C. Replacing two endpoint forces by a mean of local nonlinear Rac forces remains the distinct new circumference closure.",
            chemical_amount_distinction="The geometric concentration factor is a mathematical contribution, not an isolated evolved trajectory or observed change of either chemical species. Reaction exchanges active/inactive material and diffusion redistributes it. Their combined physical amount is conserved; their mean combined concentration changes inversely with C.",
            no_ablation="Removing dilution while retaining C(t) changes the physical-amount equation; it would not isolate the same mass-conserving model."),
        sampling_limit="Rates and extrema are evaluated at saved samples; exact integrated signed dilation is obtained from endpoint log ratios. Trapezoidal sampled quadrature is a diagnostic approximation, not the exact integral.")
    output = run / "results/geometry_dilution_20261019"
    output.mkdir(parents=True, exist_ok=False)
    analysis_sources = output / "analysis_sources"
    analysis_sources.mkdir()
    shutil.copy2(Path(__file__), analysis_sources / Path(__file__).name)
    with (output / "timecourses.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    result["timecourse_sha256"] = sha(output / "timecourses.csv")
    (output / "summary.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    print(json.dumps(dict(output=str(output), prepared_geometry=prepared, phases=summaries), indent=2, allow_nan=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, default=DEFAULT)
    args = parser.parse_args()
    main(args.run.resolve())
