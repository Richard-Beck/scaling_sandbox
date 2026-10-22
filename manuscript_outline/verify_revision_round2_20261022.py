"""Independently verify the round-2 plotted algebra and protected outline sections."""
import csv
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round2_20261022/figure4"


def rows(name):
    with (FIG / "supporting/data" / name).open() as stream:
        return list(csv.DictReader(stream))


def num(row, key):
    value = row.get(key)
    return float(value) if value not in (None, "") else math.nan


def key(row):
    return row["arm"], row["model"], row["view"], num(row, "f")


def close(a, b):
    assert math.isclose(a, b, rel_tol=2e-10, abs_tol=1e-12), (a, b)


def main():
    measurements = rows("condition_lag_metrics.csv")
    index = {(key(r), num(r, "C_label"), num(r, "lag")): r for r in measurements}
    assert len(index) == len(measurements)
    ratios = rows("size_ratios.csv")
    lag_ratios = rows("all_lag_size_ratios.csv")
    by_condition = {}
    for row in lag_ratios:
        by_condition.setdefault(key(row), []).append(row)
    checked, unavailable = 0, 0
    for row in ratios:
        if row["status"] == "unavailable_unmatched_designed_input":
            unavailable += 1
            continue
        small, large = num(row, "small_C_label"), num(row, "large_C_label")
        candidates = []
        for lag_row in sorted(by_condition[key(row)], key=lambda r: num(r, "lag")):
            lag = num(lag_row, "lag")
            s, l = index[key(row), small, lag], index[key(row), large, lag]
            usable = (num(s, "coverage") >= .5 and num(l, "coverage") >= .5
                      and (num(s, "one_minus_C_ci_lower") > 0 or num(l, "one_minus_C_ci_lower") > 0))
            assert usable == (lag_row["usable"] == "True"), (key(row), lag)
            if usable and 0 < num(s, "C") <= .98 and num(l, "C") > 0:
                candidates.append((lag, s, l))
        if row["evaluable"] == "True":
            assert candidates
            lag, s, l = candidates[0]
            close(lag, num(row, "primary_lag"))
            if row["status"] == "resolved":
                expected = math.log(num(s, "C")) / math.log(num(l, "C"))
                close(expected, num(row, "rho_k"))
                close(expected / num(row, "geometry_expectation"), num(row, "rho_geo"))
            else:
                assert row["status"] == "one_sided_lower_bound"
                assert not math.isfinite(num(row, "rho_k"))
                assert math.isinf(num(row, "rho_k_upper"))
            checked += 1
        else:
            assert not candidates
            unavailable += 1
    tolerances = rows("noise_tolerance.csv")
    crossing_checks = 0
    for row in tolerances:
        if row["status"] != "resolved":
            continue
        selected = sorted((r for r in measurements if r["arm"] == row["arm"] and r["model"] == row["model"]
                           and r["view"] == row["view"] and num(r, "C_label") == num(row, "C_label")
                           and num(r, "lag") == num(row, "reference_lag")), key=lambda r: num(r, "f"))
        for i, upper in enumerate(selected):
            eligible = (num(upper, "coverage") >= .5 and (row["forced_inclusive"] == "True" or upper["forced_combined"] == "False"))
            if eligible and num(upper, "C") <= .5:
                assert i > 0
                lower = selected[i - 1]
                assert num(lower, "C") > .5 and num(lower, "coverage") >= .5
                assert row["forced_inclusive"] == "True" or lower["forced_combined"] == "False"
                fraction = (num(lower, "C") - .5) / (num(lower, "C") - num(upper, "C"))
                expected = math.exp(math.log(num(lower, "f")) + fraction * math.log(num(upper, "f") / num(lower, "f")))
                close(expected, num(row, "fhalf"))
                crossing_checks += 1
                break
        else:
            raise AssertionError((key(row), "resolved crossing missing"))
    receipt = json.loads((FIG / "supporting/metadata/analysis_receipt.json").read_text())
    for name, expected in receipt["outputs"].items():
        assert hashlib.sha256((FIG / "supporting/data" / name).read_bytes()).hexdigest() == expected, name
    decision = json.loads((FIG / "supporting/metadata/S0_decision_replay_audit.json").read_text())
    assert decision["passed"] and decision["rows"] == decision["old_estimand_exact_replay_rows"] == 308
    old, new = [json.loads((ROOT / f"manuscript_outline/outline_v{v}.json").read_text()) for v in (4, 5)]
    for ident in ("figure1", "figure5", "figure6"):
        assert next(s for s in old["sections"] if s["id"] == ident) == next(s for s in new["sections"] if s["id"] == ident)
    result = dict(passed=True, independently_checked_primary_ratios=checked,
                  explicitly_unavailable_ratios=unavailable, independently_checked_log_input_crossings=crossing_checks,
                  output_hashes_verified=len(receipt["outputs"]), original_S0_replay_rows=308,
                  protected_sections=["figure1", "figure5", "figure6"],
                  scope="Independent CSV algebra, primary-lag selection, missing-data status, output identities and protected sections; scalar-history and bootstrap checks are in the producer/test receipts.")
    destination = ROOT / "manuscript_outline/revision_round2_20261022/independent_validation.json"
    destination.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
