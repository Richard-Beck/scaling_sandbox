# Agent Notes

- This repository is used on an HPC system. Check existing Conda environments before creating or installing a new one.
- Use `agentRrunner.R` for R work when available; it provides the expected execution environment for LLM-run R scripts.
- Use SLURM for compute-heavy work (roughly >4 cores, >10 minutes, GPU work, or >8 GB RAM). Select a QoS suited to the job shape rather than defaulting to `normal`.
- Do not add a SLURM array `%` concurrency throttle unless resource, filesystem, scheduler, or user constraints justify it.
