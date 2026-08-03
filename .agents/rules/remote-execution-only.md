# Remote Execution Only (Mandatory Rule)

This project strictly requires that **all data processing, pipeline execution, script execution, and Snakemake runs be performed on the remote server (`olm`) using the provided `Makefile`**. You must **never** run data processing scripts or Snakemake directly on the local machine.

## Strict Rules:
1. **Mandatory Remote Execution via Makefile**: All pipeline execution MUST be invoked via the `Makefile` (`make sync`, `make run-rule RULE=<rule_name>`, `make pull-results`), which synchronizes code and triggers Snakemake on the remote server (`olm`).
2. **No Local Pipeline Execution**: Do not run `python script.py`, `uv run ...`, or `snakemake` locally to process project data or run pipelines. Local execution is strictly restricted to unit tests (`pytest`).
3. **Triggering Execution**: When executing a pipeline step:
   - Run `make sync` to push code to `olm`.
   - Run `make run-rule RULE=<rule_name>` (e.g. `make run-rule RULE=plot_unified_embeddings` or `make run-unified`) to invoke Snakemake on `olm`.
   - Run `make pull-results` to download generated results back to `output/results/`.
4. **Connection Check**: Before running remote commands, verify SSH connectivity to `olm` (`ssh -q -o BatchMode=yes -o ConnectTimeout=5 olm exit`). If connection fails, stop immediately and ask the user to establish SSH access to `olm`.
