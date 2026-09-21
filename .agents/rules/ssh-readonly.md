# SSH Execution & Read-Only Policy

This project uses an Antigravity `PreToolUse` hook to automatically permit and accelerate remote SSH queries on the remote server (`olm`), provided they do not write or mutate state.

## Rules:
1. **Autonomous Read-Only Inspection**:
   - You are authorized to run read-only inspection commands over SSH without asking for manual confirmation.
   - Examples of permitted commands:
     - Directory and file inspection: `ssh olm "ls -la /storage/halu/data"`, `ssh olm "find ..."`
     - File viewing: `ssh olm "cat ..."` , `ssh olm "head -n 50 ..."` , `ssh olm "grep ..."`
     - Process monitoring: `ssh olm "ps aux | grep snakemake"`
     - Snakemake dry-runs: `ssh olm "cd ~/python-venv/tme_analysis && uv run snakemake --dry-run ..."`
     - Connection checks: `ssh -q -o BatchMode=yes -o ConnectTimeout=5 olm exit`
2. **Strict Guard on Write / State-Changing Commands**:
   - You must **never** execute remote write, delete, or modifying commands without explicit user consent.
   - Disallowed operations include:
     - File deletions or overwrites: `rm`, `mv`, `touch`, `mkdir`, `chmod`, `chown`
     - Output redirection: `>`, `>>`
     - In-place modifications: `sed -i`
     - Pipeline execution or data generation (e.g. `snakemake` without `--dry-run`, `make run-rule`, `make sync`)
   - Any attempt to run mutating commands will be caught by the `PreToolUse` gate and require manual user approval.
