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
    has_direct_login_compute_or_download,
    has_safe_read,
    has_write_operation,
    is_ssh_command,
    process_hook_payload,
)


class TestSSHGate(unittest.TestCase):
    def test_is_ssh_command(self) -> None:
        positive_cases = [
            "ssh eiger ls -la",
            "ssh -q -o BatchMode=yes eiger exit",
            "scp file.txt eiger:/ritom/scratch/cscs/lhaeuser/",
            "sftp eiger",
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
            'ssh eiger "rm -rf /ritom/scratch/cscs/lhaeuser/tmp"',
            'ssh eiger "touch new_file.txt"',
            'ssh eiger "mkdir -p output/dir"',
            'ssh eiger "mv a.txt b.txt"',
            'ssh eiger "cp a.txt b.txt"',
            'ssh eiger "chmod 755 script.sh"',
            'ssh eiger "chown root:root file"',
            'ssh eiger "echo data > results.txt"',
            'ssh eiger "echo data >> results.txt"',
            'ssh eiger "sed -i \'s/a/b/\' file.txt"',
            'ssh eiger "kill -9 1234"',
            'ssh eiger "pkill snakemake"',
            'ssh eiger "scancel 12345"',
            'ssh eiger "sbatch job.sh"',
            'ssh eiger "srun -A em25 hostname"',
            "make sync",
            "make stage-in",
            "make stage-out",
            "make run-rule RULE=preprocess",
        ]
        safe_cases = [
            'ssh eiger "ls -la"',
            'ssh eiger "cat output.log"',
            'ssh eiger "head -n 20 output.log"',
            'ssh eiger "grep -i error output.log"',
            'ssh eiger "squeue -u $USER"',
            'ssh eiger "snakemake --dry-run"',
        ]
        for cmd in write_cases:
            with self.subTest(cmd=cmd):
                self.assertTrue(has_write_operation(cmd))
        for cmd in safe_cases:
            with self.subTest(cmd=cmd):
                self.assertFalse(has_write_operation(cmd))

    def test_has_safe_read(self) -> None:
        safe_cases = [
            'ssh eiger "ls -la"',
            'ssh eiger "cat output.log"',
            'ssh eiger "head -n 10 file.csv"',
            'ssh eiger "tail -f file.log"',
            'ssh eiger "grep pattern file.txt"',
            'ssh eiger "df -h"',
            'ssh eiger "ps aux"',
            'ssh eiger "squeue -u lhaeuser"',
            'ssh eiger "sacct -j 12345"',
            'ssh eiger "sinfo -s"',
            'ssh eiger "snakemake report_single_cell_immuno_datasets --dry-run"',
            "ssh -q -o BatchMode=yes eiger exit",
        ]
        for cmd in safe_cases:
            with self.subTest(cmd=cmd):
                self.assertTrue(has_safe_read(cmd))

        self.assertFalse(has_safe_read('ssh eiger "unknown_custom_binary"'))

    def test_direct_login_compute_denied(self) -> None:
        denied_cases = [
            'ssh eiger "python scripts/05_run_milopy_da.py"',
            'ssh eiger "uv run python scripts/05_run_milopy_da.py"',
            'ssh eiger "uv sync"',
            'ssh eiger "snakemake --cores all"',
            'ssh eiger "srun -A em25 python scripts/download_cohorts.py"',
        ]
        for cmd in denied_cases:
            with self.subTest(cmd=cmd):
                self.assertTrue(has_direct_login_compute_or_download(cmd))
                decision = evaluate_command(cmd)
                self.assertEqual(decision.decision, "deny")

    def test_slurm_wrapped_compute_gated_not_denied(self) -> None:
        gated_cases = [
            'ssh eiger "srun -A em25 -p normal uv run python scripts/05_run_milopy_da.py"',
            'ssh eiger "cd /ritom/scratch/cscs/lhaeuser/code && /users/lhaeuser/.local/bin/uv run snakemake --profile workflow/profiles/eiger run_milopy_da"',
        ]
        for cmd in gated_cases:
            with self.subTest(cmd=cmd):
                self.assertFalse(has_direct_login_compute_or_download(cmd))

    def test_evaluate_command_read_only_allowed(self) -> None:
        decision = evaluate_command('ssh eiger "ls -la /ritom/scratch/cscs/lhaeuser/data-test"')
        self.assertEqual(decision.decision, "allow")
        self.assertEqual(decision.overwrite, {"BypassSandbox": True})

    def test_evaluate_command_write_gated(self) -> None:
        decision = evaluate_command('ssh eiger "rm -rf /ritom/scratch/cscs/lhaeuser/data-test/bad"')
        self.assertEqual(decision.decision, "ask")
        self.assertIn("potential write", decision.reason)

    def test_evaluate_command_redirection_gated(self) -> None:
        decision = evaluate_command('ssh eiger "cat input.txt > output.txt"')
        self.assertEqual(decision.decision, "ask")

    def test_evaluate_command_non_ssh_deferred(self) -> None:
        decision = evaluate_command("pytest tests/")
        self.assertEqual(decision.decision, "ask")
        self.assertIn("Non-SSH command", decision.reason)

    def test_evaluate_command_ambiguous_ssh_gated(self) -> None:
        decision = evaluate_command('ssh eiger "execute_custom_pipeline"')
        self.assertEqual(decision.decision, "ask")
        self.assertIn("could not be verified", decision.reason)

    def test_process_hook_payload_valid(self) -> None:
        payload = json.dumps({
            "toolCall": {
                "name": "run_command",
                "args": {"CommandLine": 'ssh eiger "head -n 50 log.txt"'},
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
