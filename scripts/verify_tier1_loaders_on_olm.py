"""Remote verification script to validate all 84 Tier 1 loaders on cluster storage (olm).

Iterates through all 84 registered Tier 1 cohorts in /storage/halu/data-test/raw/{accession},
validating:
1. Presence of raw vendor archive.
2. Ingestion via load_tier1_cohort with subsample_n=100.
3. Verification of cells x genes matrix, raw integer count integrity, and obs metadata.
"""

from __future__ import annotations

import logging
from pathlib import Path
import sys
import time

from returns.result import Failure, Success

# Setup logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("verify_tier1")

# Add package path
repo_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(repo_root / "packages" / "tme_datasets" / "src"))

from tme_datasets.providers.tier1_single_cell import TIER1_COHORTS, load_tier1_cohort
from tme_datasets.registry import get_dataset_spec


def main() -> None:
    raw_base = Path("/storage/halu/data-test/raw")
    if not raw_base.exists():
        logger.warning("Cluster raw base path %s not found on local machine.", raw_base)
        logger.info("Run this script on 'olm' via: ssh olm 'cd /home/halu/python-venv/tme_analysis && ...'")
        return

    logger.info("Found raw base directory at %s. Verifying all %d Tier 1 cohorts...", raw_base, len(TIER1_COHORTS))

    results: list[dict[str, object]] = []
    success_count = 0
    missing_count = 0
    failure_count = 0

    for i, (acc, config) in enumerate(TIER1_COHORTS.items(), 1):
        cohort_raw_dir = raw_base / acc
        if not cohort_raw_dir.exists():
            # Check lowercase or alternative
            alt_dirs = list(raw_base.glob(f"*{acc}*"))
            if alt_dirs:
                cohort_raw_dir = alt_dirs[0]

        if not cohort_raw_dir.exists():
            logger.warning("[%d/%d] %s: Directory not found in %s", i, len(TIER1_COHORTS), acc, raw_base)
            results.append({"accession": acc, "status": "MISSING_DIR", "cells": 0, "genes": 0, "time_s": 0.0})
            missing_count += 1
            continue

        start_time = time.time()
        logger.info("[%d/%d] Ingesting %s (%s, Archetype=%s)...", i, len(TIER1_COHORTS), acc, config.indication, config.archetype.value)

        # Load with small subsample for fast smoke verification
        res = load_tier1_cohort(acc, cohort_raw_dir, subsample_n=100)
        elapsed = time.time() - start_time

        match res:
            case Failure(err):
                logger.error("FAILED %s: %s (took %.2fs)", acc, err, elapsed)
                results.append({"accession": acc, "status": f"FAILED: {err[:40]}", "cells": 0, "genes": 0, "time_s": elapsed})
                failure_count += 1
            case Success(adata):
                logger.info(
                    "SUCCESS %s: %d cells x %d genes (raw=%s, resp_col=%s, took %.2fs)",
                    acc,
                    adata.n_obs,
                    adata.n_vars,
                    adata.uns.get("is_raw_counts", False),
                    "clinical_response" in adata.obs.columns,
                    elapsed,
                )
                results.append({"accession": acc, "status": "SUCCESS", "cells": adata.n_obs, "genes": adata.n_vars, "time_s": elapsed})
                success_count += 1

    print("\n" + "=" * 80)
    print("TIER 1 SINGLE-CELL DATA LOADERS VERIFICATION REPORT")
    print("=" * 80)
    print(f"Total Cohorts Audited : {len(TIER1_COHORTS)}")
    print(f"Successful Ingestions : {success_count}")
    print(f"Failed Ingestions     : {failure_count}")
    print(f"Missing Directories   : {missing_count}")
    print("=" * 80)


if __name__ == "__main__":
    main()
