#!/usr/bin/env python3
"""
Diagnostic & Verification Script for CIBERSORTx Docker Container.

Verifies:
1. Docker CLI installation and daemon status
2. Stanford CIBERSORTx credential availability (.env, environment variables, or config)
3. End-to-end containerized execution (smoke test) if credentials and daemon are present

Adheres to strict functional style, Pydantic, and Returns.
"""

from __future__ import annotations

import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import numpy as np
from returns.result import Success, Failure

from bayesprism.adapters.cibersortx_docker import (
    CibersortXDockerConfig,
    is_docker_available,
    resolve_cibersortx_credentials,
    run_cibersortx_docker,
)


def verify_cibersortx_environment() -> int:
    """Run comprehensive verification of CIBERSORTx Docker requirements."""
    print("=" * 80)
    print(" CIBERSORTx Docker Integration & Credential Verification")
    print("=" * 80)

    # 0. Host & Execution Environment Details
    host_node = platform.node()
    os_sys = platform.system()
    arch = platform.machine()
    py_exec = sys.executable
    print(f"[INFO] Host: {host_node} ({os_sys} {arch})")
    print(f"[INFO] Python Interpreter: {py_exec}")

    # 1. Check Docker Binary
    docker_bin = shutil.which("docker")
    if docker_bin:
        print(f"[OK] Docker CLI found: {docker_bin}")
        try:
            ver_proc = subprocess.run([docker_bin, "--version"], capture_output=True, text=True, check=True)
            print(f"     Version: {ver_proc.stdout.strip()}")
        except Exception:
            pass
    else:
        print("[FAIL] Docker CLI not found in PATH.")

    # 2. Check Docker Daemon Reachability
    docker_reachable = is_docker_available(docker_bin or "docker")
    if docker_reachable:
        print("[OK] Docker daemon is running and accessible.")
    else:
        print("[WARN] Docker daemon is not accessible from current environment.")
        print("       (Note: When running inside restricted sandboxes, socket access is denied.")
        print("        Execute directly in your terminal: `docker info` to verify).")

    # 3. Check Credentials
    config = CibersortXDockerConfig()
    cred_res = resolve_cibersortx_credentials(config)

    has_credentials = False
    match cred_res:
        case Success((username, token)):
            has_credentials = True
            masked_token = token[:3] + "..." + token[-3:] if len(token) > 6 else "***"
            print(f"[OK] CIBERSORTx credentials configured:")
            print(f"     Username: {username}")
            print(f"     Token:    {masked_token}")
        case Failure(err):
            print(f"[INFO] CIBERSORTx credentials NOT configured yet:")
            print(f"       {err}")
            print("\n" + "-" * 80)
            print(" HOW TO PROVIDE CIBERSORTx CREDENTIALS:")
            print("-" * 80)
            print("1. Option A (Recommended): Create or edit `.env` in the repository root:")
            print("   CIBERSORTX_USERNAME=your_registered_email@domain.com")
            print("   CIBERSORTX_TOKEN=your_stanford_api_token")
            print("\n2. Option B: Set environment variables in your terminal shell:")
            print("   export CIBERSORTX_USERNAME=\"your_registered_email@domain.com\"")
            print("   export CIBERSORTX_TOKEN=\"your_stanford_api_token\"")
            print("\n3. Option C: Pass CLI arguments directly to the benchmark:")
            print("   python scripts/regularized_deconv/benchmark_collinearity_regularized_deconv.py \\")
            print("     --cibersortx-username \"your_email\" --cibersortx-token \"your_token\"")
            print("\n4. Option D: Run via Makefile:")
            print("   make benchmark-cibersortx CIBERSORTX_USERNAME=\"your_email\" CIBERSORTX_TOKEN=\"your_token\"")
            print("-" * 80)
            print("To obtain a token: Register for a free academic account at https://cibersortx.stanford.edu")
            print("-" * 80 + "\n")

    # 4. Optional Smoke Test (if Docker daemon is reachable and credentials exist)
    if docker_reachable and has_credentials:
        print("[INFO] Running 2-sample mini smoke test with `cibersortx/fractions` container...")
        genes = tuple(f"Gene_{i:02d}" for i in range(10))
        cell_types = ("TypeA", "TypeB")
        sample_names = ("Sample_01", "Sample_02")
        sig = np.array([
            [10.0, 1.0], [8.0, 2.0], [5.0, 5.0], [2.0, 8.0], [1.0, 10.0],
            [12.0, 0.5], [0.5, 12.0], [4.0, 4.0], [9.0, 2.0], [1.0, 9.0],
        ])
        mix = np.array([
            [6.0, 7.0], [5.0, 5.0], [5.0, 5.0], [5.0, 5.0], [5.0, 5.0],
            [6.5, 6.0], [6.0, 6.5], [4.0, 4.0], [5.5, 5.5], [5.0, 5.0],
        ])
        res = run_cibersortx_docker(
            mixture=mix,
            signature=sig,
            gene_names=genes,
            cell_types=cell_types,
            sample_names=sample_names,
            config=config,
        )
        match res:
            case Success(ciber_out):
                print("[OK] Container smoke test PASSED successfully!")
                print(ciber_out.proportions)
            case Failure(err):
                print(f"[FAIL] Container smoke test failed: {err}")
                return 1

    print("=" * 80)
    print(" Verification complete.")
    print("=" * 80)
    return 0


if __name__ == "__main__":
    sys.exit(verify_cibersortx_environment())
