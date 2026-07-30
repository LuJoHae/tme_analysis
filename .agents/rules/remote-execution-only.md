# Remote Execution Only

This project strictly requires that all data processing, script execution, and pipeline building be done on the remote server (`olm`) using the provided `Makefile`. You must **never** run scripts directly on the local machine to process project data.

## Rules:
1. **No Local Execution**: Do not run `python script.py`, `uv run ...`, or `snakemake` directly on the local machine to process data or run the main pipeline. All heavy lifting must be done on the remote server.
2. **Use the Makefile**: The only permitted way to execute the pipeline or run project code is via the `Makefile` located in the root of the project.
3. **Triggering Execution**: If asked to run the code, trigger the `/run-remote-pipeline` workflow or execute the specific Makefile target `make run-rule RULE=<rule_name>` for the exact rule you are working on (along with `sync` and `pull-results` if needed). You must **never** use `make run-all` or execute the entire pipeline without an explicit rule target.
4. **Configuration Check**: Before attempting to run the pipeline, you must ensure that `REMOTE_DATA_DIR` in the `Makefile` and `DATA_DIR` in the `Snakefile` are correctly configured. If they are still set to TODO placeholders, you must ask the user for the correct paths before proceeding.
5. **Connection Check**: If there is a problem or failure when interacting with the remote server, your very first troubleshooting step must be to check if there is an active SSH connection to `olm` (e.g., by running `ssh -q -o BatchMode=yes -o ConnectTimeout=5 olm exit`). If the connection fails, **abort immediately**. Do not attempt to fix or establish the connection yourself. Tell the user that they must establish the connection to `olm` manually, as this is a common prerequisite that is often missed.
