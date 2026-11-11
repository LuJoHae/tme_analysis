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
    default_cibersortx_tar_path,
    ensure_cibersortx_image_loaded,
    is_cscs_environment,
    is_docker_available,
    resolve_cibersortx_credentials,
    resolve_container_binary,
    run_cibersortx_docker,
)


def verify_cibersortx_environment() -> int:
    """Run comprehensive verification of CIBERSORTx Docker / Podman requirements."""
    print("=" * 80)
    print(" CIBERSORTx Container (Docker / Rootless Podman) & Credential Verification")
    print("=" * 80)

    # 0. Host & Execution Environment Details
    host_node = platform.node()
    os_sys = platform.system()
    arch = platform.machine()
    py_exec = sys.executable
    cscs_mode = is_cscs_environment()
    print(f"[INFO] Host: {host_node} ({os_sys} {arch}) | CSCS offline mode: {cscs_mode}")
    print(f"[INFO] Python Interpreter: {py_exec}")

    # 1. Check Container Runtime Binary (Docker or Podman)
    resolved_bin = resolve_container_binary("auto")
    bin_path = shutil.which(resolved_bin)
    if bin_path:
        print(f"[OK] Container CLI found ({resolved_bin}): {bin_path}")
        try:
            ver_proc = subprocess.run([bin_path, "--version"], capture_output=True, text=True, check=True)
            print(f"     Version: {ver_proc.stdout.strip()}")
        except Exception:
            pass
    else:
        print("[FAIL] Neither Docker nor Podman CLI found in PATH.")

    # 2. Check Container Runtime Reachability & Offline Image Archive
    runtime_reachable = is_docker_available(resolved_bin)
    if runtime_reachable:
        print(f"[OK] Container runtime '{resolved_bin}' is running and accessible.")
    else:
        print(f"[WARN] Container runtime '{resolved_bin}' is not accessible from current environment.")
        print("       (On CSCS Eiger, run inside a compute allocation via `make verify-cibersortx`.)")

    image_ready = True
    if runtime_reachable and cscs_mode:
        tar_path = default_cibersortx_tar_path()
        print(f"[INFO] Checking CSCS offline image / archive ({tar_path})...")
        match ensure_cibersortx_image_loaded(resolved_bin, "cibersortx/fractions"):
            case Success(_):
                print(f"[OK] 'cibersortx/fractions' image is loaded in {resolved_bin} storage.")
            case Failure(img_err):
                image_ready = False
                print(f"[WARN] {img_err}")

    # 3. Check Credentials
    config = CibersortXDockerConfig(docker_binary=resolved_bin, pull_never=cscs_mode)
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

    # 4. Optional Smoke Test (if Docker/Podman runtime is reachable and credentials exist)
    if runtime_reachable and image_ready and has_credentials:
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
