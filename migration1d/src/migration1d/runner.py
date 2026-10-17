"""Run migration campaigns described by JSON configuration files."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
from typing import Any


def discover(path: Path) -> list[Path]:
    if path.is_file():
        return [path]
    if not path.is_dir():
        raise FileNotFoundError(path)
    return sorted(
        item
        for item in path.rglob("*.json")
        if item.is_file()
        and not any(part.startswith("_") for part in item.relative_to(path).parts)
    )


def load(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text())
    required = {
        "schema_version", "name", "chemistry_class", "stage", "protocol",
        "parameterization", "module", "output",
    }
    missing = required - data.keys()
    if missing:
        raise ValueError(f"{path}: missing keys {sorted(missing)}")
    if data["schema_version"] != 1:
        raise ValueError(f"{path}: unsupported schema_version {data['schema_version']!r}")
    if data["chemistry_class"] not in {"single_species", "three_species"}:
        raise ValueError(f"{path}: invalid chemistry_class {data['chemistry_class']!r}")
    if data["stage"] not in {"parameter_mapping", "calibration", "testing"}:
        raise ValueError(f"{path}: invalid stage {data['stage']!r}")
    if not isinstance(data.get("arguments", []), list):
        raise ValueError(f"{path}: arguments must be a list")
    return data


def command(config: dict[str, Any]) -> list[str]:
    return [sys.executable, "-m", str(config["module"]), *map(str, config.get("arguments", []))]


def submit_slurm(config_path: Path, config: dict[str, Any], project: Path) -> None:
    resources = config.get("resources", {})
    child = [sys.executable, str(project / "run.py"), str(config_path), "--local"]
    array_tasks = int(config.get("batch", {}).get("array_tasks", 0))
    log_directory = project / "logs" / "slurm"
    log_directory.mkdir(parents=True, exist_ok=True)
    task_token = "%A_%a" if array_tasks else "%j"
    sbatch = [
        "sbatch", "--parsable", f"--job-name=m1d-{config['name']}",
        f"--cpus-per-task={resources.get('cpus', 1)}",
        f"--mem={resources.get('memory', '4G')}",
        f"--time={resources.get('time', '02:00:00')}",
        f"--output={log_directory}/%x-{task_token}.out",
        f"--error={log_directory}/%x-{task_token}.err",
    ]
    if resources.get("qos"):
        sbatch.append(f"--qos={resources['qos']}")
    if array_tasks:
        sbatch.append(f"--array=0-{array_tasks - 1}")
    sbatch.append(f"--wrap={shlex.join(child)}")
    completed = subprocess.run(sbatch, cwd=project, text=True, capture_output=True)
    if completed.returncode:
        raise RuntimeError(
            f"sbatch rejected {config['name']}: "
            f"{(completed.stderr or completed.stdout).strip()}"
        )
    job_id = completed.stdout.strip()
    print(f"{config['name']}: submitted SLURM job {job_id}")
    merge_arguments = config.get("batch", {}).get("merge_arguments")
    if array_tasks and merge_arguments:
        merge_resources = config.get("batch", {}).get("merge_resources", {})
        merge_child = [
            sys.executable, str(project / "run.py"), str(config_path),
            "--local", "--batch-merge",
        ]
        merge_sbatch = [
            "sbatch", "--parsable", f"--job-name=m1d-{config['name']}-merge",
            f"--dependency=afterok:{job_id}",
            f"--cpus-per-task={merge_resources.get('cpus', 1)}",
            f"--mem={merge_resources.get('memory', '4G')}",
            f"--time={merge_resources.get('time', '00:30:00')}",
            f"--output={log_directory}/%x-%j.out",
            f"--error={log_directory}/%x-%j.err",
        ]
        if merge_resources.get("qos"):
            merge_sbatch.append(f"--qos={merge_resources['qos']}")
        merge_sbatch.append(f"--wrap={shlex.join(merge_child)}")
        merged = subprocess.run(
            merge_sbatch, cwd=project, text=True, capture_output=True
        )
        if merged.returncode:
            raise RuntimeError(
                f"sbatch rejected {config['name']} merge: "
                f"{(merged.stderr or merged.stdout).strip()}"
            )
        print(
            f"{config['name']}: submitted dependent merge job "
            f"{merged.stdout.strip()} afterok:{job_id}"
        )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("configs", type=Path, help="a JSON file or directory of JSON files")
    parser.add_argument("--select", action="append", default=[], help="run only this config name; repeatable")
    parser.add_argument("--list", action="store_true", help="show selected configs without running")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--local", action="store_true", help="run in this process environment")
    mode.add_argument("--submit", action="store_true", help="submit each config through SLURM")
    parser.add_argument("--continue-on-error", action="store_true")
    parser.add_argument(
        "--batch-merge", action="store_true", help=argparse.SUPPRESS
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    project = Path(__file__).resolve().parents[2]
    selected = set(args.select)
    records = [(path.resolve(), load(path)) for path in discover(args.configs)]
    records = [item for item in records if item[1].get("enabled", True)]
    if selected:
        records = [item for item in records if item[1]["name"] in selected]
    if not records:
        raise SystemExit("no enabled configs selected")
    if args.batch_merge:
        for _, config in records:
            merge_arguments = config.get("batch", {}).get("merge_arguments")
            if not merge_arguments:
                raise ValueError(f"{config['name']}: no batch.merge_arguments configured")
            config["arguments"] = list(merge_arguments)
    for path, config in records:
        print(f"{config['name']}: {config['module']} -> {config['output']}")
    if args.list:
        return 0
    failures = 0
    for path, config in records:
        if args.submit:
            submit_slurm(path, config, project)
            continue
        output = (project / str(config["output"])).resolve()
        env = os.environ.copy()
        env["MIGRATION1D_OUTPUT_DIR"] = str(output)
        env["PYTHONPATH"] = os.pathsep.join(
            filter(None, (str(project / "src"), env.get("PYTHONPATH", "")))
        )
        try:
            subprocess.run(command(config), cwd=project, env=env, check=True)
        except subprocess.CalledProcessError:
            failures += 1
            if not args.continue_on_error:
                raise
    return failures
