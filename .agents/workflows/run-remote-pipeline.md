---
name: run-remote-pipeline
description: Safely synchronizes code, runs the Snakemake pipeline on the remote server, and pulls the results back.
---

# Remote Pipeline Execution Workflow

When the user asks you to "run the pipeline", "execute the code", or explicitly calls `/run-remote-pipeline`, you MUST follow these steps exactly to ensure the code is executed safely on the remote server.

## Step 1: Pre-flight Verification
1. **Check Connection**: Verify that a connection to the remote server `olm` can be established (e.g., by running `ssh -q -o BatchMode=yes -o ConnectTimeout=5 olm exit`). If this fails, **STOP** immediately and instruct the user to establish the connection to `olm`. You cannot do this for them.
2. Inspect the `Makefile` at the root of the project. Verify that the `REMOTE_DATA_DIR` variable has been updated from the placeholder (`/TODO/UPDATE/THIS/PATH/TO/DATA`) to a valid remote path.
3. Inspect the `Snakefile`. Verify that the `DATA_DIR` variable has been updated to match the path in the `Makefile`.
4. If either path is still a placeholder, **STOP** and ask the user to provide the correct absolute path on the remote server for the data directory before proceeding. Do NOT attempt to run the pipeline.

## Step 2: Synchronize and Execute
Execute the pipeline using the `Makefile`. You can either run the steps individually or all at once.

**Run All (Recommended)**:
Execute the command: `make run-all`
This will automatically sync the codebase, run the remote snakemake command, and pull the results back.

**Step-by-step**:
If the user only wants to perform specific actions:
1. Run `make sync` to push local changes to the remote server.
2. Run `make run-remote` to trigger Snakemake on the remote server.
3. Run `make pull-results` to download the generated files back to the local `output/results/` directory.

## Step 3: Verify and Report
1. After the `make` command finishes successfully, verify that the expected files have been populated in the local `output/results/` directory using your file-listing tools.
2. Present a summary of the execution to the user, highlighting any logs from Snakemake or files that were successfully pulled.
