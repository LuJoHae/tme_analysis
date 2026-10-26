#!/usr/bin/env python3
"""SSH Command Gate for Antigravity PreToolUse Hook.

Inspects command line invocations before tool execution:
1. Auto-approves read-only SSH inspection and Slurm queue queries (squeue, sacct, sinfo, ls, etc.).
2. Strictly denies direct login-node compute execution (python, uv run, pytest, bare snakemake) on Eiger
   when not wrapped in Slurm (srun, sbatch, or snakemake --profile / --executor slurm / --dry-run).
3. Strictly denies internet dataset download commands on Eiger.
4. Requires confirmation ('ask') for commands that submit Slurm jobs or modify/write remote resources.
"""

from __future__ import annotations

import json
import re
import sys
from dataclasses import dataclass
from typing import Any, Mapping

# Direct login-node compute or download patterns strictly forbidden on CSCS Eiger
DIRECT_LOGIN_COMPUTE_PATTERNS: tuple[str, ...] = (
    r"\b(download_cohorts\.py|download_and_preprocess_all\.py|01_download_single_cell_immuno\.py|01_download_datasets\.py)\b",
)

# Disallowed write/mutation/job-submission patterns on remote or local execution
WRITE_PATTERNS: tuple[str, ...] = (
    # Output redirection operators (ignoring 2>/dev/null, 2>&1)
    r"(?<!2)(?<!&)>",
    r">>",
    # File system modifications & deletions
    r"\b(rm|rmdir|mv|cp|touch|mkdir|chmod|chown|chgrp|dd|mkfs|truncate|ln)\b",
    # Process, system & Slurm job alterations / submissions
    r"\b(kill|pkill|killall|reboot|shutdown|systemctl|service|sbatch|srun|salloc|scancel)\b",
    # In-place file edits & interactive editors
    r"\bsed\s+-[a-zA-Z]*i\b",
    r"\b(vim?|nano|emacs)\b",
    r"\btee\b",
    # Version control mutations
    r"\bgit\s+(push|commit|checkout\s+-b|merge|rebase|reset|clean|tag)\b",
    # Remote sync & execution targets that modify remote host or submit jobs
    r"\bmake\s+(sync|stage-in|stage-out|pull-results|run-remote|run-rule|run-all)\b",
)

# Known safe read-only inspection commands & Slurm status tools
SAFE_PATTERNS: tuple[str, ...] = (
    r"\b(cat|head|tail|less|more|grep|egrep|fgrep|awk|jq|cut|sort|uniq|wc)\b",
    r"\b(ls|find|stat|file|du|df)\b",
    r"\b(ps|top|uptime|whoami|id|uname|pwd|hostname|env|which)\b",
    r"\b(squeue|sacct|sinfo)\b",
    r"\b(scontrol\s+show|sacctmgr\s+show)\b",
    r"\b(git\s+(status|log|diff|branch|show))\b",
    r"\b(snakemake\s+.*--dry-run)\b",
    r"\bexit\b",
)


@dataclass(frozen=True)
class GateDecision:
    decision: str
    reason: str
    overwrite: Mapping[str, Any] | None = None

    def to_dict(self) -> dict[str, Any]:
        base: dict[str, Any] = {
            "decision": self.decision,
            "reason": self.reason,
        }
        if self.overwrite is not None:
            base["overwrite"] = dict(self.overwrite)
        return base


def is_ssh_command(cmd: str) -> bool:
    """Checks if the command invokes an SSH / remote connection binary."""
    return bool(re.search(r"(^|\s)(ssh|sftp|scp)\b", cmd))


def has_direct_login_compute_or_download(cmd: str) -> bool:
    """Detects forbidden direct compute execution in bash on Eiger without Slurm, or forbidden downloads."""
    if any(re.search(p, cmd) for p in DIRECT_LOGIN_COMPUTE_PATTERNS):
        return True

    # Check if targeting eiger and running python/uv/pytest/snakemake directly without srun/sbatch
    if re.search(r"\bssh\b.*\beiger\b", cmd):
        has_slurm_wrapper = bool(re.search(r"\b(srun|sbatch|salloc)\b", cmd))
        if not has_slurm_wrapper:
            # Direct python / pytest / uv run python execution on login node is forbidden
            if re.search(r"\b(python3?|pytest)\b", cmd) or re.search(r"\buv\s+(run|sync)\b", cmd):
                # Allow only snakemake --dry-run or snakemake --unlock or snakemake --profile / --executor slurm
                if re.search(r"\bsnakemake\b", cmd) and re.search(
                    r"(--dry-run|--unlock|--profile\b|--executor\s+slurm)", cmd
                ):
                    return False
                return True
            # Bare snakemake without --dry-run, --unlock, --profile, or --executor slurm
            if re.search(r"\bsnakemake\b", cmd) and not re.search(
                r"(--dry-run|--unlock|--profile\b|--executor\s+slurm)", cmd
            ):
                return True
    return False


def has_write_operation(cmd: str) -> bool:
    """Detects whether the command contains modifying, mutating, or Slurm submission tokens."""
    return any(re.search(pattern, cmd) for pattern in WRITE_PATTERNS)


def has_safe_read(cmd: str) -> bool:
    """Checks if the command matches known safe inspection utilities."""
    return any(re.search(pattern, cmd) for pattern in SAFE_PATTERNS)


def evaluate_command(cmd: str) -> GateDecision:
    """Pure evaluation of a command line string to determine the gate decision."""
    if not is_ssh_command(cmd):
        return GateDecision(
            decision="ask",
            reason="Non-SSH command: deferred to standard tool permissions.",
        )

    if has_direct_login_compute_or_download(cmd):
        return GateDecision(
            decision="deny",
            reason=(
                "Direct compute execution in bash on Eiger login nodes or internet dataset downloading "
                "is strictly forbidden. All compute jobs must run via Slurm (srun, sbatch, or Snakemake "
                "--profile workflow/profiles/eiger)."
            ),
        )

    if has_write_operation(cmd):
        return GateDecision(
            decision="ask",
            reason="SSH command contains potential write, delete, mutating, or Slurm job submission operations.",
        )

    if has_safe_read(cmd):
        return GateDecision(
            decision="allow",
            reason="Read-only SSH inspection command auto-approved.",
            overwrite={"BypassSandbox": True},
        )

    return GateDecision(
        decision="ask",
        reason="SSH command could not be verified as strictly read-only.",
    )


def process_hook_payload(payload_str: str) -> dict[str, Any]:
    """Processes incoming hook JSON string and returns response dictionary."""
    try:
        data = json.loads(payload_str)
        cmd = str(data.get("toolCall", {}).get("args", {}).get("CommandLine", ""))
        decision = evaluate_command(cmd)
        return decision.to_dict()
    except Exception as err:
        return {
            "decision": "ask",
            "reason": f"SSH gate error parsing payload: {err}",
        }


def main() -> None:
    input_text = sys.stdin.read()
    response = process_hook_payload(input_text)
    sys.stdout.write(json.dumps(response) + "\n")


if __name__ == "__main__":
    main()
