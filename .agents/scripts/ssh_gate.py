#!/usr/bin/env python3
"""SSH Command Gate for Antigravity PreToolUse Hook.

Inspects command line invocations before tool execution, auto-approving
read-only SSH inspection commands while requiring confirmation for commands
that may modify, write, or delete remote resources.
"""

from __future__ import annotations

import json
import re
import sys
from dataclasses import dataclass
from typing import Any, Mapping

# Disallowed write/mutation patterns on remote or local execution
WRITE_PATTERNS: tuple[str, ...] = (
    # Output redirection operators (ignoring 2>/dev/null, 2>&1)
    r"(?<!2)(?<!&)>",
    r">>",
    # File system modifications & deletions
    r"\b(rm|rmdir|mv|cp|touch|mkdir|chmod|chown|chgrp|dd|mkfs|truncate)\b",
    # Process & system alterations
    r"\b(kill|pkill|killall|reboot|shutdown|systemctl|service)\b",
    # In-place file edits & interactive editors
    r"\bsed\s+-[a-zA-Z]*i\b",
    r"\b(vim?|nano|emacs)\b",
    r"\btee\b",
    # Version control mutations
    r"\bgit\s+(push|commit|checkout\s+-b|merge|rebase|reset|clean|tag)\b",
    # Remote sync targets that modify remote host
    r"\bmake\s+(sync|pull-results|run-remote|run-rule|run-all)\b",
)

# Known safe read-only inspection commands & tools
SAFE_PATTERNS: tuple[str, ...] = (
    r"\b(cat|head|tail|less|more|grep|egrep|fgrep|awk|jq|cut|sort|uniq|wc)\b",
    r"\b(ls|find|stat|file|du|df)\b",
    r"\b(ps|top|uptime|whoami|id|uname|pwd|hostname|env)\b",
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


def has_write_operation(cmd: str) -> bool:
    """Detects whether the command contains modifying or mutating tokens."""
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

    if has_write_operation(cmd):
        return GateDecision(
            decision="ask",
            reason="SSH command contains potential write, delete, or mutating operations.",
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
