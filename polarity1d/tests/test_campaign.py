import json

import numpy as np
import pytest

from polarity1d.campaign import (
    _single_cue_kinetics_summary,
    collect_campaign,
    make_manifest,
    run_manifest_task,
    validate_manifest,
)
from polarity1d.contract import Trajectory


def small_task():
    return {
        "tier": "test",
        "model": "wave_pinning",
        "assay": "single_site",
        "protocol": {
            "name": "localized_pulse",
            "amplitude": 0.1,
            "stop": 5.0,
            "taper": 2.0,
            "sides": ["left"],
        },
        "simulation": {
            "length": 10.0,
            "cells": 41,
            "end_time": 20.0,
            "output_interval": 2.0,
            "backend": "bdf",
        },
        "metrics": {
            "analysis_start": 5.0,
            "sustained_duration": 2.0,
        },
        "replicate": 0,
        "save_trajectory": False,
    }


def test_manifest_task_execution_is_resumable_and_collectable(tmp_path):
    manifest = make_manifest("test", [small_task()], description="unit test")
    summary_path = run_manifest_task(manifest, 0, tmp_path)
    first = json.loads(summary_path.read_text())
    assert first["status"] == "success"
    assert first["task_digest"] == manifest["tasks"][0]["task_digest"]
    assert not (summary_path.parent / "trajectory.npz").exists()
    assert run_manifest_task(manifest, 0, tmp_path) == summary_path

    report = collect_campaign(manifest, tmp_path)
    assert report["complete"]
    assert report["success_count"] == 1
    assert (tmp_path / "task_summary.csv").is_file()
    assert (tmp_path / "group_summary.csv").is_file()


def test_manifest_validation_rejects_configuration_tampering():
    manifest = make_manifest("test", [small_task()], description="unit test")
    manifest["tasks"][0]["simulation"]["end_time"] = 21.0
    with pytest.raises(ValueError, match="digest"):
        validate_manifest(manifest)


def test_single_cue_kinetics_reports_sustained_gain_and_loss():
    time = np.arange(0.0, 101.0)
    x = np.linspace(0.5, 9.5, 10)
    rac = np.ones((len(time), len(x)))
    rac[(time >= 10) & (time < 60), :2] = 1.2
    trajectory = Trajectory(
        model="synthetic", time=time, x=x, fields={"rac": rac},
        observable_name="rac",
    )
    summary = _single_cue_kinetics_summary(
        trajectory, cue_start=0.0, cue_stop=40.0,
        cue_width_fraction=0.1, cue_amplitude=1.0,
        orientation_threshold=0.05, polarization_sustain_s=10.0,
        depolarization_sustain_s=20.0, cue_side="left",
    )
    assert summary["cue_kinetics_status"] == "depolarized"
    assert summary["polarization_time_s"] == 10.0
    assert summary["depolarization_time_s"] == 60.0
    assert summary["depolarization_latency_s"] == 20.0
