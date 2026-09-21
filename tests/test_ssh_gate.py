"""Unit tests for Antigravity SSH PreToolUse gate."""

from __future__ import annotations

import json
import os
import sys
import unittest

# Add .agents/scripts to sys.path to import ssh_gate
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".agents", "scripts")))

from ssh_gate import (  # noqa: E402
    evaluate_command,
    has_safe_read,
    has_write_operation,
    is_ssh_command,
    process_hook_payload,
)


class TestSSHGate(unittest.TestCase):
    def test_is_ssh_command(self) -> None:
        positive_cases = [
            "ssh olm ls -la",
            "ssh -q -o BatchMode=yes olm exit",
            "scp file.txt olm:/tmp/",
            "sftp olm",
        ]
        negative_cases = [
            "pytest tests/",
            "python script.py",
            "ls -la",
        ]
        for cmd in positive_cases:
            with self.subTest(cmd=cmd):
                self.assertTrue(is_ssh_command(cmd))
        for cmd in negative_cases:
            with self.subTest(cmd=cmd):
                self.assertFalse(is_ssh_command(cmd))

    def test_has_write_operation(self) -> None:
        write_cases = [
            'ssh olm "rm -rf /tmp/data"',
            'ssh olm "touch new_file.txt"',
            'ssh olm "mkdir -p output/dir"',
            'ssh olm "mv a.txt b.txt"',
            'ssh olm "cp a.txt b.txt"',
            'ssh olm "chmod 755 script.sh"',
            'ssh olm "chown root:root file"',
            'ssh olm "echo data > results.txt"',
            'ssh olm "echo data >> results.txt"',
            'ssh olm "sed -i \'s/a/b/\' file.txt"',
            'ssh olm "kill -9 1234"',
            'ssh olm "pkill snakemake"',
            'ssh olm "systemctl restart service"',
            'ssh olm "git push origin main"',
            "make sync",
            "make run-rule RULE=preprocess",
        ]
        safe_cases = [
            'ssh olm "ls -la"',
            'ssh olm "cat output.log"',
            'ssh olm "head -n 20 output.log"',
            'ssh olm "grep -i error output.log"',
            'ssh olm "snakemake --dry-run"',
        ]
        for cmd in write_cases:
            with self.subTest(cmd=cmd):
                self.assertTrue(has_write_operation(cmd))
        for cmd in safe_cases:
            with self.subTest(cmd=cmd):
                self.assertFalse(has_write_operation(cmd))

    def test_has_safe_read(self) -> None:
        safe_cases = [
            'ssh olm "ls -la"',
            'ssh olm "cat output.log"',
            'ssh olm "head -n 10 file.csv"',
            'ssh olm "tail -f file.log"',
            'ssh olm "grep pattern file.txt"',
            'ssh olm "df -h"',
            'ssh olm "ps aux"',
            'ssh olm "snakemake download_dataset --dry-run"',
            "ssh -q -o BatchMode=yes olm exit",
        ]
        for cmd in safe_cases:
            with self.subTest(cmd=cmd):
                self.assertTrue(has_safe_read(cmd))

        self.assertFalse(has_safe_read('ssh olm "unknown_custom_binary"'))

    def test_evaluate_command_read_only_allowed(self) -> None:
        decision = evaluate_command('ssh olm "ls -la /storage/halu/data"')
        self.assertEqual(decision.decision, "allow")
        self.assertEqual(decision.overwrite, {"BypassSandbox": True})

    def test_evaluate_command_write_gated(self) -> None:
        decision = evaluate_command('ssh olm "rm -rf /storage/halu/data/bad"')
        self.assertEqual(decision.decision, "ask")
        self.assertIn("potential write", decision.reason)

    def test_evaluate_command_redirection_gated(self) -> None:
        decision = evaluate_command('ssh olm "cat input.txt > output.txt"')
        self.assertEqual(decision.decision, "ask")

    def test_evaluate_command_non_ssh_deferred(self) -> None:
        decision = evaluate_command("pytest tests/")
        self.assertEqual(decision.decision, "ask")
        self.assertIn("Non-SSH command", decision.reason)

    def test_evaluate_command_ambiguous_ssh_gated(self) -> None:
        decision = evaluate_command('ssh olm "execute_custom_pipeline"')
        self.assertEqual(decision.decision, "ask")
        self.assertIn("could not be verified", decision.reason)

    def test_process_hook_payload_valid(self) -> None:
        payload = json.dumps({
            "toolCall": {
                "name": "run_command",
                "args": {"CommandLine": 'ssh olm "head -n 50 log.txt"'},
            }
        })
        result = process_hook_payload(payload)
        self.assertEqual(result["decision"], "allow")
        self.assertEqual(result["overwrite"], {"BypassSandbox": True})

    def test_process_hook_payload_malformed(self) -> None:
        result = process_hook_payload("invalid json {")
        self.assertEqual(result["decision"], "ask")
        self.assertIn("error parsing payload", result["reason"])


if __name__ == "__main__":
    unittest.main()
