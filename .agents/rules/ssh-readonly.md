# SSH Execution & Read-Only Policy (CSCS Alps Eiger)

This project uses an Antigravity `PreToolUse` hook (`.agents/scripts/ssh_gate.py`) to automatically permit read-only SSH and Slurm status queries on `eiger`, while blocking direct login-node compute and gating state-changing or Slurm submission commands.

## Rules:
1. **Autonomous Read-Only Inspection**:
   - You are authorized to run read-only inspection commands over SSH without manual confirmation:
     - Directory and file inspection: `ssh eiger "ls -la /ritom/scratch/cscs/lhaeuser/data-test"`, `ssh eiger "df -h"`
     - File viewing: `ssh eiger "cat ..."`, `ssh eiger "head -n 50 ..."`, `ssh eiger "tail -n 50 ..."`
     - Slurm queue and job monitoring: `ssh eiger "squeue -u \$USER"`, `ssh eiger "sacct -j <jobid>"`, `ssh eiger "sinfo -s"`
     - Snakemake dry-runs: `ssh eiger "cd /ritom/scratch/cscs/lhaeuser/code && /users/lhaeuser/.local/bin/uv run snakemake --dry-run ..."`
     - Connection checks: `ssh -q -o BatchMode=yes -o ConnectTimeout=5 eiger exit`
2. **Strict Ban on Direct Login-Node Compute & Internet Downloads**:
   - Never run `python`, `uv run`, `pytest`, or bare `snakemake` directly in bash on `eiger` login nodes; every compute job must be submitted via Slurm (`srun`, `sbatch`, or Snakemake `--profile workflow/profiles/eiger`).
   - Never download datasets from the internet on CSCS Eiger.
3. **Strict Guard on Write / State-Changing / Slurm Submission Commands**:
   - Remote write, delete, staging, or Slurm job submission commands (`sbatch`, `srun`, `scancel`, `rm`, `mv`, `mkdir`, `make sync`, `make stage-in`, `make stage-out`, `make run-rule`) require user approval.
