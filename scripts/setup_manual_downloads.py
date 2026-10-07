#!/usr/bin/env python3
"""Setup and populate data/manual_download/ and data/preprocessed/ from existing storage.

This script physically copies (not symlinks):
1. All 14 published ICI response paper H5AD files from cluster storage into data/manual_download/.
2. Maynard NSCLC H5AD into data/manual_download/Maynard_NSCLC.h5ad.
3. Genentech IMvigor210 alignment directory into data/manual_download/EGAD00001006631-align/.
4. Existing cluster preprocessed single-cell H5AD files (GSE179994, GSE123813, GSE125449)
   into data/preprocessed/ for instant sub-second loading.
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path

from tme_datasets.paths import (
    find_repo_root,
    get_manual_download_dir,
    get_preprocessed_h5ad_path,
)

PAPER_H5AD_FILES = (
    "Auslander.h5ad",
    "Chen-CTLA4.h5ad",
    "Chen-PD1.h5ad",
    "Freeman.h5ad",
    "Gide.h5ad",
    "Hugo.h5ad",
    "Lauss.h5ad",
    "Liu.h5ad",
    "Prat.h5ad",
    "Ravi.h5ad",
    "Riaz.h5ad",
    "Rose.h5ad",
    "Snyder.h5ad",
    "VanAllen.h5ad",
)

CLUSTER_LAIR_DIR = Path(
    "/storage/halu/lair/ImmuneCheckpointTherapyResponseProcessedGeneNormalizedClinicalDataNormalized"
)
CLUSTER_EGAD_DIR = Path("/storage/halu/manual-download/EGAD00001006631-align")
CLUSTER_PREPROCESSED_DIR = Path("/storage/halu/data-test/preprocessed")



def setup_manual_downloads(verbose: bool = True) -> None:
    """Populate data/manual_download/ and data/preprocessed/ via physical file copies."""
    repo_root = find_repo_root()
    manual_dir = get_manual_download_dir(repo_root=repo_root)
    manual_dir.mkdir(parents=True, exist_ok=True)

    if verbose:
        print(f"Target manual download directory: {manual_dir}")

    # 1. Copy Paper H5AD cohorts
    if CLUSTER_LAIR_DIR.exists():
        if verbose:
            print(f"\n[1/4] Copying paper H5AD cohorts from {CLUSTER_LAIR_DIR}...")
        for fname in PAPER_H5AD_FILES:
            src = CLUSTER_LAIR_DIR / fname
            dst = manual_dir / fname
            if src.exists():
                if not dst.exists() or dst.stat().st_size != src.stat().st_size:
                    if verbose:
                        print(f"  Copying {fname} ({src.stat().st_size // 1024} KB)...")
                    shutil.copy2(src, dst)
                else:
                    if verbose:
                        print(f"  Already exists: {fname}")
            else:
                if verbose:
                    print(f"  Warning: Source not found for {fname}")
    else:
        if verbose:
            print(f"\n[1/4] Cluster lair directory not present at {CLUSTER_LAIR_DIR} (skipping remote copy)")

    # 2. Copy Maynard NSCLC H5AD
    maynard_src_candidates = [
        repo_root / "jupyter/data/maynard2020_3k.h5ad",
        Path("/storage/halu/data-test/preprocessed/maynard2020_3k.h5ad"),
        Path("/storage/halu/data/preprocessed/maynard2020_3k.h5ad"),
    ]

    maynard_dst = manual_dir / "Maynard_NSCLC.h5ad"
    maynard_src = next((p for p in maynard_src_candidates if p.exists()), None)
    if maynard_src:
        if not maynard_dst.exists() or maynard_dst.stat().st_size != maynard_src.stat().st_size:
            if verbose:
                print(f"\n[2/4] Copying Maynard NSCLC H5AD from {maynard_src} -> {maynard_dst}...")
            shutil.copy2(maynard_src, maynard_dst)
        else:
            if verbose:
                print(f"\n[2/4] Maynard NSCLC H5AD already present at {maynard_dst}")
    else:
        if verbose:
            print("\n[2/4] Maynard NSCLC source file not found")

    # 3. Copy EGAD00001006631 alignment directory
    egad_dst = manual_dir / "EGAD00001006631-align"
    if CLUSTER_EGAD_DIR.exists():
        if not egad_dst.exists():
            if verbose:
                print(f"\n[3/4] Copying EGAD alignment directory from {CLUSTER_EGAD_DIR} -> {egad_dst}...")
            shutil.copytree(CLUSTER_EGAD_DIR, egad_dst)
        else:
            if verbose:
                print(f"\n[3/4] EGAD alignment directory already present at {egad_dst}")
    else:
        if verbose:
            print(f"\n[3/4] EGAD source directory not found at {CLUSTER_EGAD_DIR}")

    # 4. Copy existing cluster preprocessed single-cell H5ADs
    sc_mapping = {
        "GSE179994_processed.h5ad": "GSE179994",
        "GSE123813_processed.h5ad": "GSE123813",
        "GSE125449_processed.h5ad": "GSE125449",
        "GSE115978_processed.h5ad": "GSE115978",
        "GSE120575_processed.h5ad": "GSE120575",
    }
    if CLUSTER_PREPROCESSED_DIR.exists():
        if verbose:
            print(f"\n[4/4] Copying preprocessed single-cell H5AD caches from {CLUSTER_PREPROCESSED_DIR}...")
        for src_name, dataset_id in sc_mapping.items():
            src_file = CLUSTER_PREPROCESSED_DIR / src_name
            dst_file = get_preprocessed_h5ad_path(dataset_id, repo_root=repo_root)
            dst_file.parent.mkdir(parents=True, exist_ok=True)
            if src_file.exists():
                if not dst_file.exists() or dst_file.stat().st_size != src_file.stat().st_size:
                    if verbose:
                        print(f"  Copying {src_name} -> {dst_file.name} ({src_file.stat().st_size // (1024 * 1024)} MB)...")
                    shutil.copy2(src_file, dst_file)
                else:
                    if verbose:
                        print(f"  Already exists in preprocessed: {dst_file.name}")
    else:
        if verbose:
            print("\n[4/4] Cluster preprocessed directory not found (skipping sc cache copy)")

    if verbose:
        print("\nSetup complete! manual_download directory populated successfully.")


if __name__ == "__main__":
    setup_manual_downloads(verbose=True)
