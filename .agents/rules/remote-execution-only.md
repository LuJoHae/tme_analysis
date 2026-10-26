# Remote Execution on CSCS Alps Eiger (Mandatory Rule)

This project strictly requires that **all data processing, pipeline execution, script execution, and Snakemake runs be performed on CSCS Alps Eiger (`eiger`, Slurm account `em25`) using the provided `Makefile` and Slurm**.

## Strict Rules:
1. **100% Mandatory Slurm Execution (Never Run Compute in Bash)**:
   - You must **NEVER** run compute-intensive commands (`python`, `uv run`, `pytest`, or bare `snakemake --cores ...`) directly in bash on an Eiger login node (`eiger-ln*`).
   - All jobs MUST be scheduled by **Slurm** (`sbatch`, `srun -A em25`, or Snakemake with `--profile workflow/profiles/eiger`).
2. **Strict `$SCRATCH`-Only Runtime Reads & Writes (`/ritom/scratch/cscs/lhaeuser`)**:
   - Persistent master storage on CSCS is `$STORE` (`/capstor/store/cscs/empa/em25/halu/data-test` and `/capstor/store/cscs/empa/em25/halu/code`).
   - Compute jobs must **NEVER** read from or write to `$STORE` (`/capstor/store/...`) directly.
   - All runtime file reads and writes MUST occur exclusively on `$SCRATCH` (`/ritom/scratch/cscs/lhaeuser/data-test` and `/ritom/scratch/cscs/lhaeuser/code`).
   - Data that needs to be read must first be copied from `$STORE` to `$SCRATCH` (`make stage-in`), and files written to `$SCRATCH` that should be saved permanently must be copied back to `/capstor/store/cscs/empa/em25/halu/data-test` (`make stage-out`).
3. **Zero Internet Downloads on CSCS**:
   - Data must **NEVER** be downloaded from the internet on CSCS Eiger (`raw/` is intentionally not stored on CSCS).
   - New raw datasets must be downloaded and preprocessed on `olm` first and transferred manually via Globus to `/capstor/store/cscs/empa/em25/halu/data-test`.
4. **Triggering Execution via `Makefile`**:
   - Verify SSH connectivity to `eiger` (`ssh -q -o BatchMode=yes -o ConnectTimeout=5 eiger exit`).
   - Run `make sync` to push code to `$STORE/halu/code`, mirror it to `$SCRATCH/code`, and sync `.venv` via Slurm.
   - Run `make stage-in` to stage required datasets from `$STORE/halu/data-test` to `$SCRATCH/data-test` via Slurm.
   - Run `make run-rule RULE=<rule_name>` to execute via Snakemake's Slurm executor (`--profile workflow/profiles/eiger`) and stage permanent results back to `$STORE`.
   - Run `make pull-results` to download results back to `output/results/`.
