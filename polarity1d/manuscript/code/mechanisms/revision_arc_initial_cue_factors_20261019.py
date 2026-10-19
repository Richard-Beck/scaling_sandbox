"""Exact initial finite-width cue projection for the six selected rate examples."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"
import numpy as np
import pandas as pd
import ring_distance_prototype_20261017 as ring

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "figures/length_mechanisms_20261008/dev/revision_program_20261019/figure3/arc_illustrative_histories"
GROUP = OUT.parent


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    old_sources = json.loads((OUT / "supporting/metadata/source_sha256.json").read_text())
    for path, digest in old_sources.items():
        assert sha(path) == digest
    summary_path = GROUP / "standard_arclength_sizes/supporting/data/stationary_arclength_summary.csv"
    summary = pd.read_csv(summary_path)
    selected_path = OUT / "supporting/data/selected_conditions.csv"
    selected = pd.read_csv(selected_path)
    records = []
    for row in selected.itertuples():
        native = summary[summary.condition_id == row.condition_id].iloc[0]
        path = Path(native.actual_cue_path)
        assert sha(path) == native.actual_cue_sha256
        with np.load(path) as z:
            state, signal, x = z["states"][0], z["stimulus"], z["x"]
        model = ring.native_model(native.model)
        grid = ring.RingGrid(native.circumference, len(x))
        C, q = native.circumference, x / native.circumference
        rhs0 = ring.native_rhs(model, grid, state, np.zeros(len(x)))
        rhsS = ring.native_rhs(model, grid, state, signal)
        activation = rhsS[0] - rhs0[0]
        norm = getattr(model, "rac_total", 1.)
        expected = signal * state[1] / norm
        channel_error = float(abs(activation - expected).max())
        assert channel_error < 2e-15
        M1 = C * np.mean(state[0] * np.exp(2j * np.pi * q))
        moment = lambda du: C * np.mean(du * np.exp(2j * np.pi * q))
        omega = lambda dm: float(np.degrees(np.imag(np.conj(M1) * dm) / abs(M1)**2))
        J, A = moment(activation), float(C * activation.mean())
        I = float(C * signal.mean())
        assert A > 0 and I > 0
        phase = float(np.angle(J / M1))
        kappa = float(abs(J) / A)
        actual_target = float(native.target_angle_in_P_coordinates_deg)
        factor_product = float(np.degrees((A / abs(M1)) * kappa * np.sin(phase)))
        direct, basal, total = omega(J), omega(moment(rhs0[0])), omega(moment(rhsS[0]))
        assert abs(factor_product - direct) < 1e-12
        assert abs(total - direct - basal) < 1e-12
        assert abs(total - row.maximum_saved_native_RHS_angular_rate) < 1e-10  # selected initial peak rates
        records.append(dict(condition_id=int(row.condition_id), model=native.model, size_parameter=float(native.size_parameter),
            initial_arc=float(native.physical_distance_at_onset), actual_target_deg=actual_target,
            M1_initial_magnitude=float(abs(M1)), M1_initial_absolute_phase_deg=float(np.degrees(np.angle(M1))),
            prescribed_input_coefficient_integral=I, physical_primary_activation_amount_rate=A,
            input_weighted_native_inactive=float(signal @ state[1] / signal.sum()), native_channel_normalizer=float(norm),
            activation_per_coefficient_integral=float(A / I), exact_finite_width_coherence=kappa,
            activation_moment_angle_relative_to_initial_M1_deg=float(np.degrees(phase)),
            sin_actual_target=float(np.sin(np.radians(actual_target))), sin_activation_moment=float(np.sin(phase)),
            cue_moment_real=float(J.real), cue_moment_imag=float(J.imag),
            exact_cue_projection_rate_deg_per_native_time=direct, basal_rate_deg_per_native_time=basal,
            total_initial_rate_deg_per_native_time=total, exact_factor_product_deg_per_native_time=factor_product,
            point_centre_approximation_deg_per_native_time=float(np.degrees(A / abs(M1) * np.sin(np.radians(actual_target)))),
            native_channel_increment_identity_max_error=channel_error, actual_cue_path=str(path), actual_cue_sha256=sha(path)))
    frame = pd.DataFrame(records)
    frame.to_csv(OUT / "supporting/data/initial_cue_geometry_factors.csv", index=False)
    comparisons = []
    for label, left, right in (("Holmes same arc, increasing size", 132, 172),
                              ("Wave-pinning same size, increasing arc", 43, 58),
                              ("Goryachev same size, increasing arc", 113, 118)):
        a, b = [frame[frame.condition_id == cid].iloc[0] for cid in (left, right)]
        observed = b.exact_cue_projection_rate_deg_per_native_time / a.exact_cue_projection_rate_deg_per_native_time
        factors = dict(activation_amount_rate_ratio=float(b.physical_primary_activation_amount_rate / a.physical_primary_activation_amount_rate),
            inverse_initial_moment_ratio=float(a.M1_initial_magnitude / b.M1_initial_magnitude),
            finite_width_coherence_ratio=float(b.exact_finite_width_coherence / a.exact_finite_width_coherence),
            exact_sine_lever_ratio=float(b.sin_activation_moment / a.sin_activation_moment),
            centre_sine_lever_ratio=float(b.sin_actual_target / a.sin_actual_target))
        product = factors["activation_amount_rate_ratio"] * factors["inverse_initial_moment_ratio"] * factors["finite_width_coherence_ratio"] * factors["exact_sine_lever_ratio"]
        assert abs(product - observed) < 1e-12
        comparisons.append(dict(label=label, numerator_condition_id=right, denominator_condition_id=left,
            exact_initial_rate_ratio=float(observed), exact_factor_product=float(product), **factors))
    result = dict(generated_at=datetime.now(timezone.utc).isoformat(), status="passed", selected_initial_states=6,
        formula="omega_cue = (A / abs(M1_initial)) * kappa * sin(psi); A = integral physical native primary activation increment, kappa = abs(cue moment)/A, psi = cue moment phase relative to initialM1",
        angular_rate_units="degrees per native time unit", comparisons=comparisons,
        maximum_basal_rate_deg_per_native_time=float(abs(frame.basal_rate_deg_per_native_time).max()),
        interpretation="Exact instantaneous decomposition at each frozen initial state. It explains geometric lever arms and actual inactive-weighted activation; it is not an evolved pathway intervention or a causal allocation of the later reorientation trajectory.",
        source_sha256={str(summary_path):sha(summary_path),str(selected_path):sha(selected_path),str(Path(__file__)):sha(__file__)},
        unchanged_rate_card_sources_verified=len(old_sources), rate_card_images_and_PDFs_unchanged=True)
    (OUT / "supporting/metadata/initial_cue_geometry_factors.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__": main()
