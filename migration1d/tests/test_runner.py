from pathlib import Path
from types import SimpleNamespace

from migration1d.runner import discover, load, submit_slurm


ROOT = Path(__file__).resolve().parents[1]


def test_nested_campaign_discovery_excludes_shared_standards():
    paths = discover(ROOT / "configs")
    assert paths
    assert all("_standards" not in path.parts for path in paths)


def test_campaign_paths_match_declared_taxonomy_and_have_unique_targets():
    records = [(path, load(path)) for path in discover(ROOT / "configs")]
    names = [record["name"] for _, record in records]
    outputs = [record["output"] for _, record in records]
    assert len(names) == len(set(names))
    assert len(outputs) == len(set(outputs))

    for path, record in records:
        relative = path.relative_to(ROOT / "configs")
        assert relative.parts[0] == record["chemistry_class"]
        assert relative.parts[1] == record["stage"]
        assert record["output"].startswith(
            f"results/{record['chemistry_class']}/{record['stage']}/"
        )


def test_slurm_submission_redirects_array_and_merge_logs(tmp_path, monkeypatch):
    calls = []

    def fake_run(command, **kwargs):
        calls.append((command, kwargs))
        return SimpleNamespace(returncode=0, stdout=f"{100 + len(calls)}\n", stderr="")

    monkeypatch.setattr("migration1d.runner.subprocess.run", fake_run)
    config = {
        "name": "spatial-test",
        "batch": {"array_tasks": 3, "merge_arguments": ["--mode", "merge"]},
    }
    submit_slurm(tmp_path / "config.json", config, tmp_path)

    task_command, merge_command = (call[0] for call in calls)
    log_directory = tmp_path / "logs" / "slurm"
    assert log_directory.is_dir()
    assert f"--output={log_directory}/%x-%A_%a.out" in task_command
    assert f"--error={log_directory}/%x-%A_%a.err" in task_command
    assert f"--output={log_directory}/%x-%j.out" in merge_command
    assert f"--error={log_directory}/%x-%j.err" in merge_command
